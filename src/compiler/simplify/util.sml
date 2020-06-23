(* util.sml
 *
 * Utility code for Simplification.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure Util : sig

  (* the standard reductions extended with pseudo reductions *)
    datatype reduction = MEAN | VARIANCE | RED of Reductions.t

  (* identify a basis variable that specifies a reduction *)
    val identifyReduction : Var.t -> reduction

  (* return information about how to compute a reduction operator *)
    val reductionInfo : Reductions.t -> {
            rator : Var.t,                      (* primitive operator *)
            init : Literal.t,                   (* identity element to use for initialization *)
            mvs : SimpleTypes.meta_arg list     (* meta-variable arguments for primitive application *)
          }

  (* convert a block, which is the map part of a map-reduce, into a function by closing over
   * its free variables.  The simple variable argument is the strand object, which should always
   * be the first parameter (even if it is not referenced in the block).
   *)
    val makeFunction : string * Simple.var * Simple.block * SimpleTypes.ty
          -> Simple.func_def * Simple.var list

  (* return true if an AST constant expression is "small" *)
    val isSmallExp : AST.expr -> bool

    val deFemInput : AST.expr * Types.ty -> AST.expr list * AST.expr
    val reFem: Var.t * AST.expr list * Types.ty * Types.ty -> AST.expr

    val getNewPosVar : Var.t -> SimpleVar.t option
    val makePosAssign : Var.t * SimpleVar.t -> Simple.stmt list option
    val cvtWorldPos : SimpleVar.t * SimpleVar.t * bool -> Simple.stmt list
  end = struct

    structure S = Simple
    structure BV = BasisVars
    structure L = Literal
    structure R = RealLit
    structure VMap = SimpleVar.Map
    structure TU = TypeUtil
    structure Ty = Types

    datatype reduction = MEAN | VARIANCE | RED of Reductions.t

    fun identifyReduction rator =
          if Var.same(BV.red_all, rator)
            then RED(Reductions.ALL)
          else if Var.same(BV.red_exists, rator)
            then RED(Reductions.EXISTS)
          else if Var.same(BV.red_max_i, rator)
            then RED(Reductions.IMAX)
          else if Var.same(BV.red_max_r, rator)
            then RED(Reductions.RMAX)
          else if Var.same(BV.red_mean, rator)
            then MEAN
          else if Var.same(BV.red_min_i, rator)
            then RED(Reductions.IMIN)
          else if Var.same(BV.red_min_r, rator)
            then RED(Reductions.RMIN)
          else if Var.same(BV.red_product_i, rator)
            then RED(Reductions.IPRODUCT)
          else if Var.same(BV.red_product_r, rator)
            then RED(Reductions.RPRODUCT)
          else if Var.same(BV.red_sum_i, rator)
            then RED(Reductions.ISUM)
          else if Var.same(BV.red_sum_r, rator)
            then RED(Reductions.RSUM)
(* FIXME: variance not supported yet
          else if Var.same(BV.red_variance, rator)
            then VARIANCE
*)
            else raise Fail(Var.uniqueNameOf rator ^ " is not a reduction operator")

    val mvsReal = [SimpleTypes.SHAPE[]]

    fun reductionInfo redOp = let
          val id = Reductions.identity redOp
          in
            case redOp
             of Reductions.ALL => {rator = BV.op_and, init = id, mvs = []}
              | Reductions.EXISTS => {rator = BV.op_or, init = id, mvs = []}
              | Reductions.IMAX => {rator = BV.fn_max_i, init = id, mvs = []}
              | Reductions.RMAX => {rator = BV.fn_max_r, init = id, mvs = []}
              | Reductions.IMIN => {rator = BV.fn_min_i, init = id, mvs = []}
              | Reductions.RMIN => {rator = BV.fn_min_r, init = id, mvs = []}
              | Reductions.IPRODUCT => {rator = BV.mul_ii, init = id, mvs = []}
              | Reductions.RPRODUCT => {rator = BV.mul_rr, init = id, mvs = mvsReal}
              | Reductions.ISUM => {rator = BV.add_ii, init = id, mvs = []}
              | Reductions.RSUM => {rator = BV.add_tt, init = id, mvs = mvsReal}
            (* end case *)
          end

    local
      val n = ref 0
      fun mkFuncId (name, ty) = let val id = !n
            in
              n := id + 1;
              SimpleVar.new(name ^ Int.toString id, SimpleVar.FunVar, ty)
            end
    in
    fun makeFunction (name, strand, blk, resTy) = let
          val freeVars = ref []
          fun cvtVar (env, x) = (case VMap.find(env, x)
                 of SOME x' => (env, x')
                  | NONE => let
                      val x' = SimpleVar.copy(x, SimpleVar.FunParam)
                      in
                        freeVars := (x, x') :: !freeVars;
                        (VMap.insert(env, x, x'), x')
                      end
                (* end case *))
          fun cvtVars (env, xs) = let
                fun cvt (x, (env, xs')) = let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, x'::xs')
                      end
                in
                  List.foldr cvt (env, []) xs
                end
          fun newVar (env, x) = let
                val x' = SimpleVar.copy(x, SimpleVar.LocalVar)
                in
                  (VMap.insert(env, x, x'), x')
                end
          fun cvtBlock (env, S.Block{props, code}) = let
                fun cvtStms (env, [], stms') = (env, S.Block{props = props, code = List.rev stms'})
                  | cvtStms (env, stm::stms, stms') = let
                      val (env, stm') = cvtStm (env, stm)
                      in
                        cvtStms (env, stms, stm'::stms')
                      end
                in
                  cvtStms (env, code, [])
                end
          and cvtStm (env, stm) = (case stm
                 of S.S_Var(x, NONE) => let
                      val (env, x') = newVar (env, x)
                      in
                        (env, S.S_Var(x', NONE))
                      end
                  | S.S_Var(x, SOME e) => let
                      val (env, e') = cvtExp (env, e)
                      val (env, x') = newVar (env, x)
                      in
                        (env, S.S_Var(x', SOME e'))
                      end
                  | S.S_Assign(x, e) => let
                      val (env, e') = cvtExp (env, e)
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.S_Assign(x', e'))
                      end
                  | S.S_IfThenElse(x, b1, b2) => let
                      val (env, x') = cvtVar (env, x)
                      val (env, b1') = cvtBlock (env, b1)
                      val (env, b2') = cvtBlock (env, b2)
                      in
                        (env, S.S_IfThenElse(x', b1', b2'))
                      end
                  | S.S_Foreach(x, xs, b) => let
                      val (env, x') = cvtVar (env, x)
                      val (env, xs') = cvtVar (env, xs)
                      val (env, b') = cvtBlock (env, b)
                      in
                        (env, S.S_Foreach(x', xs', b'))
                      end
                  | S.S_New(name, args) => let
                      val (env, args') = cvtVars (env, args)
                      in
                        (env, S.S_New(name, args'))
                      end
                  | S.S_KillAll => (env, stm)
                  | S.S_StabilizeAll => (env, stm)
                  | S.S_Continue => (env, stm)
                  | S.S_Die => (env, stm)
                  | S.S_Stabilize => (env, stm)
                  | S.S_Return x => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.S_Return x')
                      end
                  | S.S_Print xs => let
                      val (env, xs') = cvtVars (env, xs)
                      in
                        (env, S.S_Print xs')
                      end
                  | S.S_MapReduce _ => raise Fail "unexpected nested MapReduce"
                (* end case *))
          and cvtExp (env, exp) = (case exp
                 of S.E_Var x => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_Var x')
                       end
                  | S.E_Lit _ => (env, exp)
                  | S.E_Kernel _ => (env, exp)
                  | S.E_Select(x, fld) => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_Select(x', fld))
                       end
                  | S.E_Apply(f, args) => let
                      val (env, args') = cvtVars (env, args)
                      in
                        (env, S.E_Apply(f, args'))
                      end
                  | S.E_Prim(f, mvs, args, ty) => let
                      val (env, args') = cvtVars (env, args)
                      in
                        (env, S.E_Prim(f, mvs, args', ty))
                      end
                  | S.E_Tensor(args, ty) => let
                      val (env, args') = cvtVars (env, args)
                      in
                        (env, S.E_Tensor(args', ty))
                      end
                  | S.E_Seq(args, ty) => let
                      val (env, args') = cvtVars (env, args)
                      in
                        (env, S.E_Seq(args', ty))
                      end
                  | S.E_Tuple xs => let
                      val (env, xs') = cvtVars (env, xs)
                      in
                        (env, S.E_Tuple xs')
                      end
                  | S.E_Project(x, i) => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_Project(x', i))
                      end
                  | S.E_Slice(x, indices, ty) => let
                      fun cvt (NONE, (env, idxs)) = (env, NONE::idxs)
                        | cvt (SOME x, (env, idxs)) = let
                            val (env, x') = cvtVar (env, x)
                            in
                              (env, SOME x' :: idxs)
                            end
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_Slice(x', indices, ty))
                      end
                  | S.E_Coerce{srcTy, dstTy, x} => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_Coerce{srcTy=srcTy, dstTy=dstTy, x=x'})
                       end
                  | S.E_BorderCtl(BorderCtl.Default x, y) => let
                      val (env, x') = cvtVar (env, x)
                      val (env, y') = cvtVar (env, y)
                      in
                        (env, S.E_BorderCtl(BorderCtl.Default x', y'))
                      end
                  | S.E_BorderCtl(ctl, x) => let
                      val (env, x') = cvtVar (env, x)
                      in
                        (env, S.E_BorderCtl(ctl, x'))
                      end
                  | S.E_LoadSeq _ => (env, exp)
                  | S.E_LoadImage _ => (env, exp)
                  | S.E_InsideImage(pos, img, s) => let
                      val (env, pos') = cvtVar (env, pos)
                      val (env, img') = cvtVar (env, img)
                      in
                        (env, S.E_InsideImage(pos', img', s))
                      end
                  | S.E_FieldFn _ => (env, exp)
		  | S.E_ExtractFemItem(v, ty, opt) =>
		    let
		     val (env, v') = cvtVar(env, v)
		    in
		     (env, S.E_ExtractFemItem(v', ty, opt))
		    end
		  | _ => raise Fail ("While simplifying \"" ^ name ^ "\" ran into an error in cvtStm while processing something?\n")
                (* end case *))
        (* the initial environment always includes the strand variable *)
          val (env, _) = cvtVar (VMap.empty, strand)
          val (env, blk) = cvtBlock (env, blk)
          val (args, params) = ListPair.unzip (List.rev (! freeVars))
          val f = SimpleFunc.new (name, resTy, List.map SimpleVar.typeOf params)
          in
            (S.Func{f=f, params=params, body=blk}, args)
          end
    end (* local *)

    fun isSmallExp (AST.E_Lit _) = true
      | isSmallExp (AST.E_Tensor(exps, _)) = (List.length exps <= 4)
      | isSmallExp _ = false

			 

    fun deFemInput(e, originTy) =
	let
	 val datas = ref []
	 fun doit (e, Ty.T_Bool, false) = e
	   | doit (e, Ty.T_Int, false) = e
	   | doit (e, Ty.T_String, false) = e
	   | doit (e, Ty.T_Sequence(ty, r), false) = e
	   | doit (e, Ty.T_Tensor _, false) = e
	   | doit (AST.E_Seq(vs, t), Ty.T_Sequence(ty, r), true) =
	     (AST.E_Seq(List.map (fn x => doit(x, ty, true)) vs, #1 (TypeUtil.normalilzeFemInputTy(t))))
	   | doit (e, Ty.T_Tuple(tys), false) = e
	   | doit (AST.E_Tuple(vals, tys'), t as Ty.T_Tuple(tys), true) =
	     let val Ty.T_Tuple(newtys) = #1 (TypeUtil.normalilzeFemInputTy(t)) in 
	      AST.E_Tuple(ListPair.map (fn (x,y) => doit(x, y, TU.hasFem y)) (vals, tys), newtys)
	     end
	   | doit (e, Ty.T_Named(_, ty), b) = doit(e, ty, b)
	   | doit (AST.E_LoadFem(data, SOME(femDep), NONE), Ty.T_Fem(data', _), true) = raise Fail "impossible: non-value fem in value input!"
	   | doit (AST.E_LoadFem(data, SOME(femDep), SOME(v)), Ty.T_Fem(data', _), true) =
	     (datas := (femDep :: !datas); v)
	   | doit (AST.E_ExtractFemItemN([mesh], _, _, (FemOpt.InvalidBuild, _), _), Ty.T_Fem(data', _), true) =

	     let val dim = FemData.underlyingDim data'
		 val nan = AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(RealLit.nan))), Ty.vecTy dim)
		 val int = AST.E_Lit(Literal.intLit (~1))
		 val ret = AST.E_Tuple([nan, int, int], [Ty.vecTy dim, Ty.T_Int, Ty.T_Int])
	     in (datas := mesh :: !datas; ret)
	     end
	   | doit (AST.E_Coerce {srcTy, dstTy, e}, _, b) =
	     let
	      val srcTy' = #1 (TypeUtil.normalilzeFemInputTy(srcTy))
	      val dstTy' = #1 (TypeUtil.normalilzeFemInputTy(dstTy))
	     in
	      AST.E_Coerce {srcTy=srcTy', dstTy=dstTy', e=doit(e, dstTy, b)}
	     end
	 val ret = doit(e, originTy, TU.hasFem originTy)
	in
	 (!datas, ret)
	end

    fun reFem(srcVar, datas, inputTy, targetTy) =
	let
	 val fems = ref (List.rev datas) (*TODO: Unclear this is correct choice.*)
	 val span = Var.locationOf srcVar
	 fun getVar() = (case !fems
			  of [] => raise Fail "impossible"
			   | f::fms => (fems := fms; f))
	 fun doit (e, _, _, false) = e
	   | doit (e, Ty.T_Sequence(ty, _), t as Ty.T_Sequence(ty', NONE), true) =
	     let
	      (*create itter var, *)
	      (*TODO: correct this.*)
	      val itVar = Var.new(Atom.atom "0_itter", span, Var.IterVar, ty)
	      val iter = (itVar, e)
	      val cleaned = doit(AST.E_Var((itVar, span)), ty, ty', true)
	     in
	      AST.E_Comprehension(cleaned, iter, t)
	     end
	       (*tuple, invidual selection rules, named*)
	   | doit (e, Ty.T_Sequence(ty, SOME _), t as Ty.T_Sequence(ty', SOME(k)), true) =
	     let
	      val Ty.DimConst(r) = k
	      val slicef = fn idx => doit(AST.E_Slice(e, [SOME(AST.E_Lit(Literal.intLit idx))], ty), ty, ty', true)
	     in
	      AST.E_Seq(List.tabulate(r, slicef), t)
	     end
	   | doit (e, Ty.T_Tuple(tys), t as Ty.T_Tuple(tys'), true) =
	     let
	      val k = List.length tys
	      val slicef = fn idx => doit(AST.E_Slice(e, [SOME(AST.E_Lit(Literal.intLit idx))], List.nth(tys, idx)),
					  List.nth(tys, idx),
					  List.nth(tys', idx),
					  TU.hasFem (List.nth(tys', idx)))
	     in
	      AST.E_Tuple(List.tabulate(k, slicef), tys')
	     end
	   | doit(e, Ty.T_Int, Ty.T_Fem(data as FemData.MeshCell(m), name), true) =
	     let
	      val femDep = getVar()
	     in
	      AST.E_LoadFem(data, SOME(femDep), SOME(e))
	     end
	   | doit (e, Ty.T_Int, Ty.T_Fem(data as FemData.FuncCell(f), name), true) =
	     let
	      val femDep = getVar()
	     in
	      AST.E_LoadFem(data, SOME(femDep), SOME(e))
	     end
	   | doit (e, Ty.T_Tuple _, t as Ty.T_Fem(data as FemData.MeshPos(m), name), true) =
	     (*TODO: assumption tuple is only possible is dangerous*)
	     let
	      val femDep = getVar()
	     in
	      AST.E_ExtractFemItemN([femDep], [TypeOf.expr femDep], t, (FemOpt.InvalidBuild, data), NONE)
	     end
	in
	 doit(AST.E_Var(srcVar, span), inputTy, targetTy, TU.hasFem targetTy)
	end

    local
     structure STy = SimpleTypes

     fun cvtSimple(x : SimpleVar.t) : SimpleVar.t = raise Fail "no worldpos set for this one"

     val {setFn, peekFn, ...} = SimpleVar.newProp cvtSimple
     val {setFn = markStatePosVar, getFn=getStatePosVar, ...} = SimpleVar.newProp cvtSimple
						  
     fun newTemp (ty as STy.T_Image _) = SimpleVar.new ("img", SimpleVar.LocalVar, ty)
       | newTemp ty = SimpleVar.new ("_t", SimpleVar.LocalVar, ty)

     fun createV(v, e) = S.S_Var(v, SOME(e))
						  

     fun cvtWorldPos(posVar : SimpleVar.t, worldPosAssign : SimpleVar.t, assign : bool) =
	 (case (peekFn posVar)
	   of SOME(unpos) =>
	      if assign
	      then [S.S_Assign(worldPosAssign, S.E_Var(unpos))]
	      else [S.S_Var(worldPosAssign, SOME(S.E_Var(unpos)))]
	    | NONE =>
	      let
	       val STy.T_Fem(ms as FemData.MeshPos(meshData)) = SimpleVar.typeOf posVar
	       val x' = posVar
	       val meshHolder = FemData.Mesh(meshData)
	       val meshCellHolder = FemData.MeshCell(meshData)
	       val meshTy = STy.T_Fem(meshHolder)
	       val dim = FemData.meshDim (meshData)
	       val newTensor = STy.T_Tensor([dim])
					   
	       val refPosTemp = newTemp newTensor;
	       val meshTemp = newTemp meshTy;
	       val cellTemp = newTemp STy.T_Int
	       val negOneTemp = newTemp STy.T_Int
	       val boolTemp = newTemp STy.T_Bool
	       val assignTemp = newTemp newTensor

	       val fieldTy = STy.T_Field({diff=NONE, dim=dim, shape=[dim]});
	       val fieldTemp = newTemp fieldTy;

	       val getCell = createV(cellTemp, S.E_ExtractFemItem(x', STy.T_Int, (FemOpt.CellIndex, ms)))
	       val getMesh = createV(meshTemp, S.E_ExtractFem(x', meshHolder))
	       val fieldAssign = createV(fieldTemp, S.E_FemField(meshTemp, meshTemp, SOME(cellTemp), fieldTy, FemOpt.Transform, NONE))
	       val getRefCell = createV(refPosTemp, S.E_ExtractFemItem(x', newTensor, (FemOpt.RefPos, ms)))
	       val metaArgs = [STy.DIFF(NONE), STy.DIM(dim), STy.SHAPE([dim])]
	       val worldPos = createV(assignTemp, S.E_Prim(BasisVars.op_probe, metaArgs, [fieldTemp, refPosTemp], newTensor))
	       val worldPosReplace = if assign
				     then S.S_Assign(worldPosAssign, S.E_Var(assignTemp))
				     else S.S_Var(worldPosAssign, SOME(S.E_Var(assignTemp)))

	       val stms = [getCell, getRefCell, getMesh, fieldAssign, worldPos, worldPosReplace]

	      in
	       stms
	      end
	 )
			          
     fun cvt x =
	 let
	  val kind = Var.kindOf x
	  val name = Var.nameOf x
	  val ty = Var.monoTypeOf x
	  val dim = (case ty
		      of Ty.T_Tensor(Ty.Shape([Ty.DimConst(d)])) => SOME(d)
		       | Ty.T_Fem(FemData.MeshPos(meshData), _) => SOME(FemData.meshDim meshData)
		       | _ => NONE
		    (*end case*))
	 in if (kind = Var.StrandStateVar
		orelse kind = Var.StrandOutputVar)
	       andalso name = FemName.pos
	       andalso Option.isSome dim
	    then SOME(SimpleVar.new ("_pos", kind, STy.T_Tensor([Option.valOf dim])))
	    else NONE
	 end
     val {getFn=getNewPosVar, ...} = Var.newProp cvt

     fun makePosAssign(x, x') =
	 (case getNewPosVar x
	   of NONE => NONE
	    | SOME(x'') => (case SimpleVar.typeOf x'
			     of STy.T_Tensor([d])
				=> (SOME([S.S_Assign(x'',S.E_Var(x'))]))
			      | STy.T_Fem(ms as FemData.MeshPos(meshData)) =>
				let
				 val meshHolder = FemData.Mesh(meshData)
				 val meshCellHolder = FemData.MeshCell(meshData)
				 val meshTy = STy.T_Fem(meshHolder)
				 val dim = FemData.meshDim (meshData)
				 val newTensor = STy.T_Tensor([dim])
							     
				 val refPosTemp = newTemp newTensor;
				 val meshTemp = newTemp meshTy;
				 val cellTemp = newTemp STy.T_Int
				 val negOneTemp = newTemp STy.T_Int
				 val boolTemp = newTemp STy.T_Bool
				 val assignTemp = newTemp newTensor
				 val startAssign = S.S_Var(assignTemp, NONE)

				 val fieldTy = STy.T_Field({diff=NONE, dim=dim, shape=[dim]});
				 val fieldTemp = newTemp fieldTy;

				 val getCell = createV(cellTemp, S.E_ExtractFemItem(x', STy.T_Int, (FemOpt.CellIndex, ms)))
				 val getMesh = createV(meshTemp, S.E_ExtractFem(x', meshHolder))
				 val fieldAssign = createV(fieldTemp, S.E_FemField(meshTemp, meshTemp, SOME(cellTemp), fieldTy, FemOpt.Transform, NONE))
				 val getRefCell = createV(refPosTemp, S.E_ExtractFemItem(x', newTensor, (FemOpt.RefPos, ms)))
				 val negOneAsssign = createV(negOneTemp, S.E_Lit(Literal.intLit (~1)))
				 val getValid = createV(boolTemp, S.E_Prim(BV.neq_ii, [], [cellTemp, negOneTemp], STy.T_Bool))

				 (* S.E_ExtractFemItem(x', STy.T_Bool, (FemOpt.Valid, ms))) *)
							  
				 val inf = S.E_Lit(Literal.Real(RealLit.posInf))
				 val infinityInits = List.tabulate(dim, fn x =>
									   let val v = newTemp STy.realTy
									   in (v, createV(v, inf)) end)
				 val inits = List.map (fn (x,y) => y) infinityInits
						      
				 val infity = S.E_Tensor(List.map (fn (x, y) => x) infinityInits, newTensor)
				 val badAssign = S.S_Assign(assignTemp, infity)
				 val metaArgs = [STy.DIFF(NONE), STy.DIM(dim), STy.SHAPE([dim])]
				 val worldPos = S.S_Assign(assignTemp, S.E_Prim(BasisVars.op_probe, metaArgs, [fieldTemp, refPosTemp], newTensor))


				 val ifStm = S.S_IfThenElse(boolTemp,
							    S.Block{props = PropList.newHolder(), code =  [getRefCell, getMesh, fieldAssign, worldPos]},
							    S.Block{props = PropList.newHolder(), code =  [badAssign]})
				 val actualFin = S.S_Assign(x'', S.E_Var(assignTemp))


							   (*test valid*)
							   (*var ret;*)
							   (*build if condition*)
							   (**)

				in
				 (* SOME([fin, fieldAssign, getRefCell, getMesh, getCell]) *)
				 SOME([actualFin, ifStm, getValid]@inits@[getCell, startAssign, negOneAsssign])
				end
			   (*end case*))
			     
	 (*end case*))
	   
    in
    val getNewPosVar = getNewPosVar
    val makePosAssign = makePosAssign
    val cvtWorldPos = cvtWorldPos

    val markStatePosVar = markStatePosVar
    val getStatePosVar = getStatePosVar
    end
	  


  end
