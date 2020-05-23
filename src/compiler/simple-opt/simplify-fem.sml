(* simplify-fem.sml
 *
 *
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)


structure SimplifyFem : sig
	   val transform : Simple.program -> Simple.program
	  end = struct



 (*
Structure of simple

Inputs
GlobalInits/Globals
Create
GlobalStart

(Strands:
stateInit
startM
updateM
stabilizeM
Do News
)
Global Update 

return to Strands until Fixed

 *)

  (*Step 0: Build a datastructure for tracking FEM Types and Inputs
-Find all base FEM inputs and categorize by type (hashmap)
-Via loadFems in globalinit, annotate current inputs with new fems to make the graph -> props on inputs
   *)

  (*step 1: in globalInit/globals, figure out mapping to inputs -> this should be via props*)

  (*step 2: Run throug the create area and register to strand params props*)

  (*step 3: run initial strand Init to popular strand vars*)

  (*step 4: run global start*)

  (*step 5: Strand
5.1: run any inits if we need to (if params have changed via new or whatever)
5.2: run startM any strats if we need too (if params changed via new or whatever)
5.3: run updateM
5.4: run stabilizeM
5.5. Record news
   *)

  (*step 6: run global update*)
  (*step 7: if changed, go back to strand*)


  (*Questions we need to ask: what are they, what is the dependent val and how do we convert it to an array for the global*)
  (*Map from type -> vals, map from vals -> dependent vals*)

  (*details: reductions (boring as just going inside) - functions - lots of inling.*)

  structure S = Simple
  structure Ty = SimpleTypes
  structure V = SimpleVar
  structure F = SimpleFunc
  structure VTbl = SimpleVar.Tbl
  structure FD = FemData
  structure FO = FemOpt
  structure HT = HashTable
  (*plan: as above*)
  (*1. Proc globals and global inits*)
  (*1.5 make tool for updating femPress (merge two togeather)*)
  (*2. make prop for assigning femPres and make sure props make sense - they do because see inline*)
  (*3. Make expression rules for pushing things - skip apply*)
  (*4. Go through stages igoning apply and reduction*)
  (*5. Fix functions *)
  (*6. Do rewrite*)

  (*add property for function*)

  local
   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn v => NONE : S.func_def option)
  in
  val registerFuncDef = setFn
  val getfuncDef = getFn
  fun registerFunctions((F as S.Func{f,...})::funcs) = (registerFuncDef(f, SOME(F));registerFunctions(funcs))
    | registerFunctions([]) = ()
  end
		   
  fun collectVarsExp (e) =
      (case e
	of S.E_Var x => [x]
	 | S.E_Lit _ => []
	 | S.E_Kernel _ => []
	 | S.E_Select(y, z) => (
	  [y])
	 | S.E_Apply(f, xs) => xs
	 | S.E_Prim(_, _, xs, _) => xs
	 | S.E_Tensor(xs, _) => xs
	 | S.E_Field(xs, _) => xs
	 | S.E_Seq(xs, _) => xs
	 | S.E_Tuple xs => xs
	 | S.E_Project(x, i) => [x]
	 | S.E_Slice(x, _, _) => [x]
	 | S.E_Coerce{x, ...} => [x]
	 | S.E_BorderCtl(BorderCtl.Default x, y) => [x, y]
	 | S.E_BorderCtl(_, x) => [x]
	 | S.E_LoadSeq _ => []
	 | S.E_LoadImage _ => []
	 | S.E_LoadFem(_,x,y) => [x,y]
	 | S.E_ExtractFem(v,_) => [v]
	 | S.E_ExtractFemItem(v,_,_) => [v]
	 | S.E_ExtractFemItem2(v1, v2,_,_, _) => [v1,v2]
	 | S.E_ExtractFemItemN(vars,_,_, _, _) => vars
	 | S.E_FemField(v1,v1',v2,_,_,_) => v1::v1'::(Option.getOpt(Option.map (fn x => [x]) v2, []))
	 | S.E_InsideImage(pos, img, _) => [pos, img]
	 | S.E_FieldFn f => []
      (* end case *))
  and collectStm(s) =
      (case s
	of S.S_Var(v, SOME(x)) => collectVarsExp x
	 | S.S_Var(v, NONE) => []
	 | S.S_Assign(v, x) => collectVarsExp x
	 | S.S_IfThenElse(v, b1, b2) => v::((collectBlock b1)@((collectBlock b2)))
	 | S.S_Foreach(x, xs, b) => (x::xs:: (collectBlock b))
	 | S.S_New(_, vs) => vs
	 | S.S_Return(v) => [v]
	 | S.S_Print(vs) => vs
	 | S.S_MapReduce(mrs) => raise Fail "fixme"
	 | S.S_KillAll => []
	 | S.S_StabilizeAll => []
	 | S.S_Continue => []
	 | S.S_Die => []
	 | S.S_Stabilize => []
      (*end case*))
  and collectBlock(S.Block{code,...}) = List.foldr (fn (x,y) => (collectStm x) @ y) [] code



  (*Make data structures for FEM possiblie assocations: list, seq, Tuple, ... - yay for SSA in conversions*)

  datatype femPres = Base of FD.femType * V.t option | ALL of FD.femType | NOTHING (* all, nothing, or one thing *)
		     | Tuple of femPres list | Array of femPres (* []  -> [seqTy]*)
		     | Seq of femPres (* must be homogenous i.e no partitioning on an array*)


  fun merge(a1 : femPres, a2 : femPres, change : bool ref) =
      let
       fun changed v = (change := true; v)
       fun merge'(t1,t2) = (case (t1, t2)
		   of (Base(f1, SOME(v1)), Base(f2, SOME(v2))) => if FemData.same(f1, f2)
								  then if V.same(v1, v2)
								       then Base(f1, SOME(v1))
								       else changed (ALL(f1))
								  else raise Fail "impossible: fem types different"
		    | (Base(f1, SOME(v1)), Base(f2, NONE)) => if FemData.same(f1, f2)
							      then changed(Base(f1, SOME(v1)))
							      else raise Fail "impossible: fem types different"
		    | (Base(f1, NONE), Base(f2, SOME(v1))) => if FemData.same(f1, f2)
							      then changed(Base(f1, SOME(v1)))
							      else raise Fail "impossible: fem types different"
		    | (Base(_, NONE), Base(_, NONE)) => raise Fail "can't merge unit femPres; should be replaced"
		    | (NOTHING, NOTHING) => NOTHING
		    | (Base _, ALL(f)) => changed(ALL(f))
		    | (ALL(f), Base _) => changed(ALL(f))
		    (*end possible base cases*)
		    | (Tuple(tys1), Tuple(tys2)) => Tuple(ListPair.map merge' (tys1, tys2))
		    | (Array(ty1), Array(ty2)) => Array(merge' (ty1, ty2))
		    | (Seq(t1), Seq(t2)) => Seq(merge'(t1, t2))
		    | _ => raise Fail "impossible: merging incompatible femPress - must of come from different types"
		 (*end case*))
      in
       (merge'(a1, a2))
      end

  fun merge'(a1, a2) = (merge(a1, a2, ref false))

  fun tyToFemPres (Ty.T_Bool) = NOTHING
    | tyToFemPres (Ty.T_Int) = NOTHING
    | tyToFemPres (Ty.T_String) = NOTHING
    | tyToFemPres (Ty.T_Tensor _) = NOTHING
    | tyToFemPres (Ty.T_Sequence(t, SOME k)) = Array(tyToFemPres t)
    | tyToFemPres (Ty.T_Sequence(t, NONE)) = Seq(tyToFemPres t)
    | tyToFemPres (Ty.T_Tuple(tys)) = Tuple(List.map tyToFemPres tys)
    | tyToFemPres (Ty.T_Strand _) = NOTHING
    | tyToFemPres (Ty.T_Field _) = NOTHING
    | tyToFemPres (Ty.T_Fem(f)) = if FD.baseFem f
				  then Base(f, NONE)
				  else Base(Option.valOf(FD.dependencyOf(f)), NONE)

  fun tyToFemPresSeq (Ty.T_Bool) = NOTHING
    | tyToFemPresSeq (Ty.T_Int) = NOTHING
    | tyToFemPresSeq (Ty.T_String) = NOTHING
    | tyToFemPresSeq (Ty.T_Tensor _) = NOTHING
    | tyToFemPresSeq (Ty.T_Sequence(t, SOME k)) = Array(tyToFemPresSeq t)
    | tyToFemPresSeq (Ty.T_Sequence(t, NONE)) = Seq(tyToFemPresSeq t)
    | tyToFemPresSeq (Ty.T_Tuple(tys)) = Tuple(List.map tyToFemPresSeq tys)
    | tyToFemPresSeq (Ty.T_Strand _) = NOTHING
    | tyToFemPresSeq (Ty.T_Field _) = NOTHING
    | tyToFemPresSeq (Ty.T_Fem(f)) = if FD.baseFem f
				     then ALL(f) (* handle empty seq case *)
				     else Base(Option.valOf(FD.dependencyOf(f)), NONE)					   
	
  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => tyToFemPres(V.typeOf(v)))
  in
  fun updateFemPres(v : V.t, p : femPres) = (case peekFn(v)
					      of SOME(p') => (let val change = ref false
								  val p'' = merge(p, p', change)
							      in setFn(v, p''); !change end)
					       | NONE => (setFn(v, p); true)
					    (*end case*))
  fun updateFemPresRef(v : V.t, p : femPres, changed) = (case peekFn(v)
					      of SOME(p') => (let val p'' = merge(p, p', changed)
							      in setFn(v, p'') end)
					       | NONE => (setFn(v, p); changed := true)
					    (*end case*))
					      
  fun setInput(v : V.t) =
      (case V.kindOf v
	of Var.InputVar => (case tyToFemPres(V.typeOf(v))
			     of Base(f, NONE) => setFn(v, Base(f, SOME(v)))
			      | alt => setFn(v, alt) (*only Base(f, NONE) can appear in base fem inputs and only base fem inputs exist for fem inputs*)
			   (*end case*))
	 | _ => raise Fail "invalid use of setInput"
      (*end case*))

  fun checkExistence msg v = (case peekFn v
			       of NONE => raise Fail (msg ^ " var not defined")
				| SOME _ => ()
			     (*end case*))

  val getFemPres = getFn
  end


  datatype input_init = datatype Inputs.input_init

  datatype input = datatype Inputs.input

  (*goes through inputs and adds base feminputs to a registery*)
  fun procBaseInputs(inputs : V.t input list, modifyInputTable : FD.femType * V.t -> unit) =
      let
       fun doit(S.INP{var, name, ty, desc, init}) =
	   ((case ty
	     of APITypes.FemData(f) => if FemData.baseFem f
				       then modifyInputTable(f, var)
				       else raise Fail "impossible: non base fems use converted inputs!"
	      | t => (case List.filter (FemData.baseFem) (APITypes.allFems t)
		       of [] => ()
			| _ => raise Fail "impossible: compound types with base fem not allowed in inputs!")
	    (*end case*)); setInput(var))
      in
       List.app doit inputs
      end


  (*Write expressions/statements with - trickies are function calls, for each - recall if/thenelse handled by property nature -> for each, a block needs to consider uniformly -> not really because itter comes from somewhere - either a sequence or a thing
  In an assign - basic outline is peet at expressions -> reason based on type of expression what the fempres is -> continue on.
  Step:
  make merge better ergonomics
  
  
   *)


  type callSite =  V.t list  (* structure for idying a call site *)
  type functionCall  = femPres list * femPres list (* arguments and globals *)
					      
  (* think about functions...*)
  (*preProc globals in each func*)
  (*call site to args+globals femPres and args+globals femPres to new function defs*)
  (*
  1. Callsite can be found via traversing htings and map into a femPres for args + femPres for global
  2. From the femPres for args (non-globs) + femPress for globals
  3. From the args + fem for globals, we can find a functoin def/copy one over.
  *)
  (*SOL: 
  1. For each function, we isolate globals and relize that as a prop (see analyzeSimple)
  2. For each function, we have a callSite -> femPres for args + femPres for globs
  3. For each femPess for args + femPress for globals -> new function definition (copied) * ret 
  --if only baseFem in args, skip -> use ty 
   *)

  fun expToFemPres(exp : S.exp,
		   retTy : Ty.ty,
		   getDep : V.t -> V.t option,
		   addDep : V.t * V.t -> unit,
		   changeOuter : bool ref) : femPres  =
      let
       val change = ref false
       val retHasFem = Ty.hasFem retTy 
       val mergeFemPreses = List.foldr (fn (x,y) => merge(x, y, change))
       fun doit (S.E_Var(v)) = (checkExistence ("_var_" ^ (V.nameOf v)) v; getFemPres v)
	 | doit (S.E_Lit(l)) = NOTHING
	 | doit (S.E_Kernel(k)) = NOTHING
	 | doit (S.E_Select(v1, v2)) = (
	  checkExistence ("_strand_var_" ^ (V.nameOf v2)) v2;
	  getFemPres v2 (* read strand state *)
	 )
	 | doit (S.E_Apply(f, vs)) = (
	  List.app (checkExistence ((F.nameOf f) ^ "_arg")) vs;
	  raise Fail "fixme"; NOTHING
	 )
	 | doit (S.E_Prim(v, [], vs, t)) =
	   (List.app (checkExistence ("_arg")) vs;
	    if Var.same(v, BasisVars.subscript) andalso retHasFem
	    then let val [sy, index] = vs (*this merge might be pointless*)
		     val Array(p') = getFemPres sy
		     (* val retP = tyToFemPres retTy *)
		     (* val merged =  merge(retP, ps, change) *)
		 in p' end
	    else if Var.same(v, BasisVars.dynSubscript) andalso retHasFem
	    then let val [sy, index] = vs
		     val Seq(p) = getFemPres sy
		 in
		  p
		 end
	    else if Var.same(v, BasisVars.at_dT)
	    then let val [sy, a] = vs
		     val Seq(p) = getFemPres sy
		     val p' = getFemPres a
		     val p'' = merge(p, p', change)
		 in Seq(p'') end
	    else if Var.same(v, BasisVars.at_Td)
	    then let val [a, sy] = vs
		     val Seq(p) = getFemPres sy
		     val p' = getFemPres a
		     val p'' = merge(p, p', change)
		 in Seq(p'') end
	    else if Var.same(v, BasisVars.at_dd)
	    then let val [sy1, sy2] = vs
		     val Seq(p1) = getFemPres sy1
		     val Seq(p2) = getFemPres sy2
		     val p = merge(p1, p2, change)
		 in Seq(p) end
	    else if retHasFem
	    then raise Fail "impossible fem basis ret"
	    else tyToFemPres retTy 
	   )
	 | doit (S.E_Tensor(vs, _)) = (List.app (checkExistence ( "_ten")) vs; NOTHING)
	 | doit (S.E_Field(vs, _)) = (List.app (checkExistence ("_fld")) vs; NOTHING)
	 | doit (S.E_Seq(vs, t)) =
	   (List.app (checkExistence ("_seq")) vs;
	   (case t
	     of Ty.T_Sequence(t', NONE) => if List.length vs <> 0
					   then Seq(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs))
					   else Seq(tyToFemPresSeq t') (*{} -> ALL*)
	      | Ty.T_Sequence(t', SOME(k)) => if k <> 0
					      then Array(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs))
					      else Array(tyToFemPresSeq t')
	   (*end case*)))
	 | doit (S.E_Tuple(vs)) =
	   (List.app (checkExistence ("_tpl")) vs;
	    Tuple(List.map getFemPres vs)
	   )
	 | doit (S.E_Project(v, i)) =
	   let val _ =  checkExistence ("_project_" ^ (Int.toString i) ^ "_") v
	    val Tuple(vs) = getFemPres v
	   in List.nth(vs, i)
	   end
	 | doit (S.E_Slice(v, _, _)) = (((checkExistence ("_slice")) v; NOTHING))
	 | doit (S.E_Coerce {srcTy, dstTy, x}) =
	   (checkExistence "_coerce" x;
	   (case (srcTy, dstTy)
	     of (Ty.T_Int, Ty.T_Tensor _) => NOTHING
	      | (Ty.T_Sequence(ty, SOME n), Ty.T_Sequence(ty', NONE)) =>
		let 
		 val Array(p') = getFemPres x
		 val p = if n = 0
			 then tyToFemPresSeq ty  (*QUESTION: this is a bet weird as we have the p' - for {}, p' would already be all so this is redudant*)
			 else merge((tyToFemPres ty), p', change) (*if there is a FEM, there is an init so Base(x, SOME k) gets used*)
		in
		 Seq(p)
		end
	      | _ => NOTHING (*kernel to kernel, sequence to sequence, int to tensor, same; tuple handled earlier*)
	   (*end case*)))
	 | doit (S.E_BorderCtl(b, v)) = ((checkExistence ( "_b")) v; NOTHING)
	 | doit (S.E_LoadSeq _) = NOTHING
	 | doit (S.E_LoadImage _) = NOTHING
	 | doit (S.E_LoadFem(f, v1, v2)) =
	   (List.app (checkExistence ("_loadFem")) [v1, v2]; 
	    if FD.baseFem f (*v1 is the dest data; v2 is the dep; *)
	    then ((case getFemPres v2
		    of Base(f', SOME(v2')) => (addDep(v1, v2'); Base(f, SOME(v1)))
		     | Base(f', NONE) => raise Fail "base fem load not set; impossible!" 
		     | ALL(f') => raise Fail "base fem load unable to trace!" (*Should be impossible*)
		     | _ => raise Fail "impossible!"
		 (* end case*)))
	    else Base(f, SOME(v1))) (*v1 is base data (func or mesh); v2 is an int for a cell?*)
	 | doit (S.E_ExtractFem(v, f)) = ( (*if the return has fem, this is getting a func or space or mesh*)
	  let
	   val _ = checkExistence "_extractFem" v
	   val retFem = (case retTy
			  of Ty.T_Fem(f') => f'
			   | _ => raise Fail "invalid ExtractFem")
	   val _ = if (FD.baseFem retFem) andalso (FD.baseFem f) (*extractFem's return type is f*)
		   then ()
		   else raise Fail "invalid extractFem"
	   val vty = V.typeOf v
	   val vtyf = (case vty
			of Ty.T_Fem(f') => f'
			 | _ => raise Fail "invalid ExtractFem Arg")
			
	  in
	   if FD.baseFem vtyf
	   then (case getFemPres v (*This is either taking something out of a func/space so we get a dep if possible*)
		  of Base(vtyf', SOME(v')) => (case (FD.dependencyOf vtyf', getDep v')
						of (SOME(f'), SOME(v'')) => Base(f',SOME(v''))
						 | (SOME(f'), NONE) => ALL(f') 
						 | (NONE, _) => raise Fail "impossible mesh extract in simple"
					      (* end case*))
		   | Base(vtyf', NONE) => raise Fail "extracting from non init"
		   | ALL(vtyf') => ALL(vtyf')
		   | _ => raise Fail "impossible extract FEM"
		(* end case*))
	   else getFemPres v (*take base out of cell or pos; just propogate same info*)
	  end)
	 | doit (S.E_ExtractFemItem(v, t, fo)) = (checkExistence "_EFI1" v;
						  if retHasFem
						  then raise Fail "impossible"
						  else tyToFemPres retTy (*should not have Base/ALL*))
	 (*above could be:c cells, refCell -> but these should be eliminated in prior phases*)
	 | doit (S.E_ExtractFemItem2(v1, v2, t1, t2, fo)) = (checkExistence "_EFI2_1" v1;
							     checkExistence "_EFI2_2" v2;
							     if retHasFem
							     then raise Fail "impossible"
							     else tyToFemPres retTy)
	 | doit (S.E_FemField(v1, v2, v3o, t, fof, func)) = (checkExistence "_FField_1" v1;
							     checkExistence "_FField_2" v2;
							     Option.app (checkExistence "_FField_3") v3o;
							    (*FIXME: handle function*)
							     tyToFemPres retTy)
	 | doit (S.E_ExtractFemItemN(vs, tys, t, fo, func)) = ((List.app (checkExistence "_FFIN") vs);
	   if Bool.not retHasFem
	   then tyToFemPres retTy
	   else
	    (*refBuild, InvalidBuild, InvalidbuildBoundary, AllBuild  *)
            (* -> scan vars for baseFem -> kind it -> huh.*)
	    let
	     (*from fo, we find key fem -> find relevant ty -> get var*)
	     (*FIXME: handle function*)
	     val (opt, data) = fo
	     val n = List.length tys
	     val tabed = List.tabulate(n, fn x => (x, List.nth(tys, x)))
	     val pos  = List.filter
			  (fn (_, Ty.T_Fem(f')) => FD.same(data, f') | _ => false) tabed


	    in
	     (case pos
	       of [(idx, _)] => getFemPres (List.nth(vs, idx)) (*propogate Base/ALL*)
		| _ => raise Fail "not planned in extractFemItemN")
	    end)
							       
	 | doit (S.E_InsideImage(v1, v2, _)) = (List.app (checkExistence ("_in")) [v1, v2]; NOTHING)
	 | doit (S.E_FieldFn _) = NOTHING (*FIXME/QUESTION, should the func be procedded? Probably, but ignore/FIXME:*)
      in
       doit(exp)
      end
  and procStatement(stm : S.stmt, changed : bool ref,
		    ret : callSite option,
		    news : (Atom.atom * V.t list) list ref, expToFemPres) =
      let
       fun doit (S.S_Var(v, NONE)) = () (*ignore but we should check for problems here*)
	 | doit (S.S_Var(v, SOME(e))) = updateFemPresRef(v, expToFemPres(e), changed)
	 | doit (S.S_Assign(v, e)) = updateFemPresRef(v, expToFemPres(e), changed) (*QUESTION: should we check existence?*)
	 | doit (S.S_IfThenElse(v, b1, b2)) = (checkExistence "condition" v; doBlock(b1); doBlock(b2))
	 | doit (S.S_New(a,vs)) = news := (a, vs) :: !news
	 | doit (S.S_Foreach(itter, src, blk)) = (checkExistence "itter" itter;
						  checkExistence "itterSrc" src;
						  doBlock(blk))
	 | doit (S.S_KillAll) = ()
	 | doit (S.S_Continue) = ()
	 | doit (S.S_Die) = ()
	 | doit (S.S_Stabilize) = ()
	 | doit (S.S_Return(v)) = (case ret
				 of NONE => ()
				  | SOME(site) => raise Fail "use call site")
	 | doit (S.S_Print(vs)) = List.app (checkExistence "print") vs
	 | doit (S.S_MapReduce(maps)) = raise Fail "mapreduce"
       and doBlock(S.Block{code, ...}) = List.app doit code
      in
       doit(stm)
      end
		   

	

  fun transform prog = let
   val S.Program{props, consts, inputs, constInit, globals,
		 globInit, funcs, strand, create, start, update} = prog

   val _ = registerFunctions(funcs)
   val S.Strand{name, params, spatialDim,  state, stateInit, startM, updateM, stabilizeM} = strand
   (*sml docs said: explain the use and semantics of mkTable HERE.*)
   val femTable : (FD.femType, V.t list) HT.hash_table = HT.mkTable(FD.hash, FD.same) (32, Fail "femtype not added yet")

   fun addType (ty, v) = (case HT.find femTable ty
			   of SOME(vs) => HT.insert femTable (ty, v::vs)
			    | NONE => HT.insert femTable (ty, [v])
			 (*end case*))

   (*find index (ty, v)*)


   val femDep = VTbl.mkTable(32, Fail "fem dep not found")
   fun addDep (v1 : V.t, v2 : V.t) = VTbl.insert femDep (v1, v2)
   fun findDep ( v1 : V.t) : V.t option = VTbl.find femDep v1


   (*Figure out the base fem inputs.*)
   val _ = procBaseInputs(inputs, addType)

   (*Go through global init to find dependecies between fems   fun procGlobalInput - base expression, but we register dependencies*)
   (*Idea: go through this as normal, but if you have a global var, we check the assignment for its type.
   If clean, we add a dependency
   handleGlobalAssign callback - 
    *)
   val _ = ()

  in
   prog
  end
  end
