(* check-simple.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure CheckSimple : sig

    val check : string * Simple.program -> bool

  end = struct
structure TU = TypeUtil
structure II = ImageInfo
    structure S = Simple
    structure V = SimpleVar

    structure SF = SimpleFunc
    structure F = SimpleFunc
    structure Ty = SimpleTypes

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

    fun cvtTy ty = (case ty
		     of Types.T_Var(Types.TV{bind, ...}) => (case !bind
							of NONE => raise Fail "unresolved type variable"
							 | SOME ty => cvtTy ty
						      (* end case *))
		      | Types.T_Bool => Ty.T_Bool
		      | Types.T_Int => Ty.T_Int
		      | Types.T_String => Ty.T_String
		      | Types.T_Sequence(ty, NONE) => Ty.T_Sequence(cvtTy ty, NONE)
		      | Types.T_Sequence(ty, SOME dim) => Ty.T_Sequence(cvtTy ty, SOME(TU.monoDim dim))
		      | Types.T_Tuple(tys) => Ty.T_Tuple(List.map cvtTy tys)
		      | Types.T_Strand id => Ty.T_Strand id
		      | Types.T_Kernel _ => Ty.T_Kernel
		      | Types.T_Tensor shape => Ty.T_Tensor(TU.monoShape shape)
		      | Types.T_Image{dim, shape} =>
			Ty.T_Image(II.mkInfo(TU.monoDim dim, TU.monoShape shape))
		      | Types.T_Field{diff, dim, shape} => Ty.T_Field{
							diff = TU.monoDiff diff,
							dim = TU.monoDim dim,
							shape = TU.monoShape shape
						       }
		      | Types.T_Named(name, def) => cvtTy def
		      | Types.T_Fem(femData, _) => Ty.T_Fem(femData)
		      | Types.T_Fun(tys1, ty2) => raise Fail "unexpected T_Fun in Simplify"

		      | Types.T_Error => raise Fail "unexpected T_Error in Simplify"
		   (* end case *))
	  

    fun checkUnusedVar(vs, defined, msgs, place) = let
     val usedVars =  V.Set.addList(V.Set.empty, vs)
     val difs = V.Set.listItems (V.Set.difference(usedVars, defined))
     val msg = "[" ^ (String.concatWith ", " (List.map V.uniqueNameOf difs)) ^ "]"
     val () = if V.Set.isSubset(usedVars, defined)
	      then ()
	      else msgs := (place ^ " contains undefined vars:" ^ msg) :: !msgs

    in
     ()
    end
    fun checkExp(ass, e, msgs, defined) = let
     val vname = V.uniqueNameOf ass
     fun femFirst v1 = (case V.typeOf v1
			 of Ty.T_Fem _ => ()
			  | _ => msgs := (vname ^ " allows load fem with no fem data") :: !msgs
		       (* end case*))
     fun baseFemFirst v = (case V.typeOf v
			    of Ty.T_Fem(f) => if FemData.baseFem f
					      then ()
					      else msgs := (vname ^ " allows femfield with no base field") :: !msgs
			     | _ => ())(*cause earlier*)
     val _ = checkUnusedVar(collectVarsExp e, defined, msgs, vname)

     val _ = if Ty.same(V.typeOf ass, S.typeOf e)
	     then ()
	     else (msgs := (vname ^ " allows invalid type to be assigned!") :: !msgs)

    in
     (case e
       of S.E_Var _ => ()
	| S.E_Lit _ => ()
	| S.E_Select(strand, v) =>
	  (case (V.typeOf(strand), V.kindOf(v))
	    of (Ty.T_Strand _, Var.StrandStateVar) => ()
	     | (Ty.T_Strand _, Var.StrandOutputVar) => ()
	     | _ => msgs := (vname ^ " allows invalid strand select: " ^ (V.uniqueNameOf v) ^ "\n") :: ( !msgs)
	  (* end case*))
	| S.E_Apply(f, vs) => let val (_, argTys) = SF.typeOf f
			      in if ListPair.all (Ty.same) (argTys, List.map V.typeOf vs)
				 then ()
				 else msgs := (vname ^ " allows poorly typed function call via " ^ (F.nameOf f)) :: !msgs
			      end
	| S.E_Prim(v, mal, vs, ty) => ()
	(* let val Types.T_Fun(args, ret) = Var.monoTypeOf v *)
	(* 	  val args' = List.map cvtTy args *)
	(* in if ListPair.all (Ty.same) (args', List.map V.typeOf vs) *)
	(* 	 then () *)
	(* 	 else msgs := (vname ^ " allows poorly typed prim call via " ^ (Var.nameOf v)) :: !msgs *)
	(* end *)
	(*FIXME: what should be checked for tensor, field, and tuple?*)
	| S.E_Tensor(vs, ty) => ()
	| S.E_Field(vs, ty) => ()
	| S.E_Tuple(vs) => ()
	| S.E_Seq(vs, ty) => (case ty
			       of Ty.T_Sequence(ty', _) => if List.all (fn x => Ty.same(V.typeOf x, ty')) vs
							   then ()
							   else msgs := (vname ^ " allows poorly typed sequence") :: !msgs
				| _ => msgs := (vname ^ " allows poorly typed sequence call" :: !msgs)
			     )
	| S.E_Project(v, i) => (case V.typeOf v
				 of Ty.T_Tuple _ => ()
				  | _ => msgs := (vname ^ " allows project of non-tuple " ^ (V.nameOf v)) :: !msgs
			       (* end case*))
	(*FIXME: slice incorrect*)
	| S.E_Slice _ => ()
	| S.E_Coerce {srcTy, dstTy, x} => if Ty.same(srcTy, V.typeOf x)
					  then ()
					  else msgs := (vname ^ " allows coerce with incorrect inner type") :: !msgs
	(*FIXME: border control and co*)
	| S.E_BorderCtl _ => ()
	| S.E_LoadSeq _ => ()
	| S.E_LoadImage _ => ()
	| S.E_LoadFem(f, v1, v2) => femFirst v1
	(*Could more precisely check fem constraints:*)
	| S.E_ExtractFem(v, _) => femFirst v
	| S.E_ExtractFemItem(v, _, _) => femFirst v
	| S.E_ExtractFemItem2(v1, v2, otherTy, outTy, _) =>
	  (femFirst v1; if Ty.same(V.typeOf v2, otherTy)
			then ()
			else msgs := (vname ^ " allows invalid typed extractFemItem2 via v2="^(V.uniqueNameOf v2)) :: !msgs
	  (*end case*))
	| S.E_ExtractFemItemN(v::vs, ty::tys, outTy, _, NONE) => (
	 femFirst v; if ListPair.all (fn (x,y) => Ty.same(V.typeOf x, y)) (v::vs, ty::tys)
		     then ()
		     else msgs := (vname ^ " allows poorly typed extractFemItemN") :: !msgs
	)
	| S.E_ExtractFemItemN(_, _, _, _, SOME _) => msgs := (vname ^ " allows invalid extractFemItemN with function call!") :: !msgs
	| S.E_FemField(v1, v2, pi, ty, _, _) =>
	  (femFirst v1; femFirst v2; baseFemFirst v1; baseFemFirst v2 (*better checking of this*)
	  )
	| S.E_InsideImage _ => ()
	| S.E_FieldFn _ => ()
     (*end case*))
    end

    fun checkStm(s, msgs, params, defined,
		 allowGlobal : bool, retTy : Ty.ty option) = let
     fun checkExp'(dst, e) = checkExp(dst, e, msgs, defined)
     fun addVar(v : V.t) : V.Set.set = V.Set.add(defined, v)
     fun doBlock'(b) = doBlock(b, msgs, params, defined, allowGlobal, retTy)
     fun checkReduce(S.MapReduce{result, reduction, mapf, args, source, domain}, defined) = (*TODO:check reduction types*)
	 let
	  val resultTy = V.typeOf result
	  val resultName  = V.uniqueNameOf result
	  val S.Func{f, params, body} = mapf
	  val (fTy, fArgsTy) = F.typeOf f
	  val _ = if Ty.same(resultTy, fTy)
		  then ()
		  else msgs := ("reduction " ^ resultName ^ " has wrong result type") :: !msgs
	  val _ = if ListPair.all Ty.same (List.map V.typeOf args, fArgsTy)
		  then ()
		  else msgs := ("reduction " ^ resultName ^ " has wrong argument types" :: !msgs)

	  val defined' = V.Set.addList(V.Set.add(defined, source), params)
	  val _ = doBlock(body, msgs, NONE, defined', false, SOME(fTy))
	 in
	  V.Set.add(defined, result)
	 end
    in
     (case s
       of S.S_Var(v, NONE) => (addVar(v))
	| S.S_Var(v, SOME e) => (checkExp'(v, e); addVar(v))
	| S.S_Assign(v, e) => (checkExp'(v,e); addVar(v))
	| S.S_IfThenElse(v, b1, b2) => (checkUnusedVar([v], defined, msgs, "if"); doBlock' b1; doBlock' b2; defined)
	| S.S_Foreach(v, vs, b) =>
	  (checkUnusedVar([vs], defined, msgs, "itter"); doBlock(b, msgs, params, V.Set.add(defined,v), allowGlobal, retTy); defined) (*remove itter*)
	| S.S_New(_, vs) => (case params
			      of NONE => (msgs := " invalid new statement outside strand" :: !msgs;defined)
			       | SOME(vs') => if ListPair.all (fn (x,y) => Ty.same(V.typeOf x, V.typeOf y)) (vs, vs')
					      then (defined)
					      else (msgs := " invalid new statement types" :: !msgs ;defined)
			    (* end case*))
	| S.S_KillAll => if allowGlobal
		       then defined
		       else (msgs := " invalid kill all " :: !msgs; defined)
	| S.S_StabilizeAll => if allowGlobal
			      then defined
			      else (msgs := " invalid stabilize all" :: !msgs; defined)
	| S.S_Continue => if Option.isSome params
			  then defined
			  else (msgs := "invalid continue" :: !msgs; defined)
	| S.S_Die => if Option.isSome params
			  then defined
			  else (msgs := "invalid die" :: !msgs; defined)
	| S.S_Stabilize => if Option.isSome params
			  then defined
			   else (msgs := "invalid stabilize" :: !msgs; defined)
	| S.S_Return ret => (case retTy
			      of NONE => (msgs := "invalid return outside function" :: !msgs; defined)
			       | SOME(t) => if Ty.same(V.typeOf ret, t)
					    then defined
					    else (msgs := ((V.uniqueNameOf ret) ^ " has invalid return type.") :: !msgs; defined)
			    (* end case*))
	| S.S_Print(vs) => (checkUnusedVar(vs, defined, msgs, "print"); defined)
	| S.S_MapReduce(maps) => List.foldr checkReduce defined maps
     (*end case*))
    end
    and doBlock(S.Block{code,...}, msgs, params, defined, ag, retTy) =
	let fun pd(S.S_Var(x, _): S.stmt, s) = (s)
	      | pd (_, s) = s
	in List.foldl (fn (x,y) => pd(x, checkStm(x, msgs, params, y, ag, retTy))) defined (code) end


    fun checkProgramTypes(prog) =
	let
	 val S.Program{
          props, consts, inputs, constInit, globals, funcs,
          globInit, strand, create, start, update
         } = prog

	 val errMsgs = ref []
	 val defined = V.Set.empty
	 val defined = V.Set.addList(defined, List.map Inputs.varOf inputs)
	 val defined = V.Set.addList(defined, consts)
	 val defined = doBlock(constInit, errMsgs, NONE, defined, false, NONE)
	 val defined = V.Set.addList(defined, globals)
	 val defined = doBlock(globInit, errMsgs, NONE, defined, false, NONE)

	 val _ = List.map (fn S.Func{f, params, body} => doBlock(body, errMsgs, NONE, V.Set.addList(defined,params), false, SOME(F.resultTypeOf f))) funcs

	 val S.Strand{name, params, spatialDim, state, stateInit, startM, updateM, stabilizeM} = strand

	 val defined = doBlock(Create.createCode create, errMsgs, SOME(params), defined, false, NONE)
	 val defined = Option.getOpt(Option.map (fn x => doBlock(x, errMsgs, NONE, defined, true, NONE)) start, defined)
	 val defined = Option.getOpt(Option.map (fn x => doBlock(x, errMsgs, NONE, defined, true, NONE)) update, defined)

	 val defined = V.Set.addList(defined, params)
	 val defined = doBlock(stateInit, errMsgs, SOME(params), defined, false, NONE)
	 val defined = V.Set.addList(defined, state)
	 val defined = doBlock(updateM, errMsgs, SOME(params), defined, false, NONE)
	 val defined = Option.getOpt(Option.map (fn x => doBlock(x, errMsgs, SOME(params), defined, false, NONE)) startM, defined)
	 val defined = Option.getOpt(Option.map (fn x => doBlock(x, errMsgs, SOME(params), defined, false, NONE)) stabilizeM, defined)

	in
	 List.rev(!errMsgs)
	end
										  
    fun check (phase, prog) = let
          val S.Program{
                  props, consts, inputs, constInit, globals, funcs,
                  globInit, strand, create, start, update
          } = prog


	  val errMsgs = checkProgramTypes(prog)
	  val errors = Bool.not(List.length errMsgs = 0)
	  val _ = if List.length errMsgs = 0
		  then ()
		  else print(String.concatWith "\n" errMsgs)
          in
(* FIXME *)errors
          end

  end
