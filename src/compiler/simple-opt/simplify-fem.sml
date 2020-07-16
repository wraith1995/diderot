(* simplify-fem.sml
 *
 *
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)


structure SimplifyFem : sig
	   val transform : Simple.program -> Simple.program

	   (*debug only:*)
	   val analysisOnly : Simple.program -> Simple.program

	   val getPropString : SimpleVar.t -> string
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

  (*plan: as above*)
  (*1. Proc globals and global inits*)
  (*1.5 make tool for updating femPress (merge two togeather)*)
  (*2. make prop for assigning femPres and make sure props make sense - they do because see inline*)
  (*3. Make expression rules for pushing things - skip apply*)
  (*4. Go through stages igoning apply and reduction*)
  (*5. Fix functions *)
  (*6. Do rewrite*)

  (*Now we encounter a call site:
   We have a call site-> we add globs and add params
   We copy, register new props, and hold the func_def there*)
  (*call site to args+globals femPres and args+globals femPres to new function defs*)
  (*
  1. Callsite can be found via traversing htings and map into a femPres for args + femPres for global
  2. From the femPres for args (non-globs) + femPress for globals we get a function def (which is copied if it doesn't exist)
  3. From the args + fem for globals, we can find a functoin def/copy one over.

  *)
  (*SOL: 
  1. For each function, we isolate globals and relize that as a prop (see analyzeSimple)
  2. For each function, we have a callSite -> femPres for args + femPres for globs
  3. For each femPess for args + femPress for globals -> new function definition (copied) * ret 
  --if only baseFem in args, skip -> use ty 
   *)
  (*Function should take:
   -depManage functions 
   -inner_change flag
   -block

   =>
   news
   *)
  (*^changes, getFn, Returns, applies, no fem function def
  Is getFn working correctly?
  Are the changes operating correctly? Being recordered correctly?
  Are callsites and news being managed correctly?
  Are stored functions working correctly? (check fem gen)
   *)


  structure S = Simple
  structure Ty = SimpleTypes
  structure V = SimpleVar
  structure F = SimpleFunc
  structure VTbl = SimpleVar.Tbl
  structure FD = FemData
  structure FO = FemOpt
  structure HT = HashTable

  (*Functions for managing function defs and globals used by a function:*)
  local
   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn v => NONE : S.func_def option)
  in
  val registerFuncDef = setFn
  val getfuncDef = getFn
  fun registerFunctions((F as S.Func{f,...})::funcs) = (registerFuncDef(f, SOME(F));registerFunctions(funcs))
    | registerFunctions([]) = ()
  end
  local	   
   fun collectVarsExp (e) =
       (case e
	 of S.E_Var x => [x]
	  | S.E_Lit _ => []
	  | S.E_Kernel _ => []
	  | S.E_Select(y, z) => (
	   [y])
	  | S.E_Apply(f, xs) =>
	    let val SOME(S.Func{f, params, body}) = getfuncDef f
		val more = collectBlock(body)
	    in more@xs end
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

   fun globsInFunc(f) =
       let val SOME(S.Func{f, params, body}) = getfuncDef f
	   val more = collectBlock(body)
	   val more' = List.filter (fn x => (case (V.kindOf x)
					      of Var.GlobalVar => true
					       | Var.InputVar => true
					       | Var.ConstVar => false
					       | _ => false
				   (*end case*))) more
       in more' end
   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn f => [] : V.t list)
  in
  fun registerGlobVarsInFunc(F as S.Func{f, params, body}) = setFn(f,globsInFunc(f))
  val getFuncGlobs = getFn
  end

  (*function to copy function defs for symbolic inlining effectively*)
  fun functionCopy(F as S.Func{f, params,body}) =
      let
       val env = V.Tbl.mkTable(10 * List.length params, Fail "copy table fail")
       val params' = List.map (fn x => V.copy(x, V.kindOf x)) params
       val _ = ListPair.app  (fn (x, x') => V.Tbl.insert env (x, x')) (params, params')

       fun rename x = (case V.Tbl.find env x
			of SOME x' => x'
			 | NONE => if V.hasGlobalScope x
				   then x
				   else raise Fail ("unknown var " ^ V.uniqueNameOf x)
		      (*end case*))
       fun add x x' = V.Tbl.insert env (x, x')
       fun copy x = let val x' = (V.copy(x, V.kindOf x)) in (add x x'; x') end
       fun remove x = V.Tbl.remove env x

       fun doStm(stm) =
	   (case stm
	     of S.S_Var(x, optE) => let val x' = copy x
				       val optE' = Option.map (doExp) optE
				   in S.S_Var(x', optE') end
	      | S.S_Assign(x, e) => (S.S_Assign(rename x, doExp e))
	      | S.S_IfThenElse(x, b1, b2) => S.S_IfThenElse(rename x, doBlock b1, doBlock b2)
	      | S.S_Foreach(x, xs, blk) => let val x' = copy x (* add itter *)
					       val xs' = rename xs
					       val blk' = doBlock blk
					   in (S.S_Foreach(x', xs', blk')) end
	      | S.S_New _ => raise Fail "new impossible in function!"
	      | S.S_KillAll => raise Fail "kill impossible in function!"
	      | S.S_StabilizeAll => raise Fail "stab impossible in function!"
	      | S.S_Continue  => raise Fail "continue impossible in function!"
	      | S.S_Die => raise Fail "die impossible in function"
	      | S.S_Return x => S.S_Return(rename x)
	      | S.S_Print(vs) => S.S_Print(List.map rename vs)
	      | S.S_MapReduce _ => raise Fail "map reduce impossible in function!"
	   (*end case*))
       and doBlock(S.Block{props, code}) = S.Block{props = props, code = List.map doStm code}
       and doExp e =
	   (case e
	     of S.E_Var x => S.E_Var(rename x)
	      | S.E_Select(x, fld) => S.E_Select(rename x, fld)
	      | S.E_Lit _ => e
	      | S.E_Kernel _ => e
	      | S.E_Apply(f, xs) => S.E_Apply(f, List.map rename xs) (*FIXME:?*)
	      | S.E_Prim(f, tys, xs, ty) => S.E_Prim(f, tys, List.map rename xs, ty)
	      | S.E_Tensor(xs, ty) => S.E_Tensor(List.map rename xs, ty)
	      | S.E_Seq(xs, ty) => S.E_Seq(List.map rename xs, ty)
	      | S.E_Tuple(xs) => S.E_Tuple(List.map rename xs)
	      | S.E_Project(x, i) => S.E_Project(rename x, i)
	      | S.E_Slice(x, idxs, ty) => S.E_Slice(rename x, idxs, ty)
	      | S.E_Coerce{srcTy, dstTy, x} => S.E_Coerce{srcTy=srcTy, dstTy=dstTy, x=rename x}
	      | S.E_BorderCtl(ctl, x) =>
		S.E_BorderCtl(BorderCtl.map rename ctl, rename x)
	      | S.E_LoadSeq _ => raise Fail "unexpected load seq during inlining"
              | S.E_LoadImage _ => raise Fail "unexpected load image during inlining"
              | S.E_InsideImage _ => raise Fail "unexpected InsideImage during inlining"
              | S.E_FieldFn _ => e
	      | S.E_ExtractFem(v, data) => S.E_ExtractFem(rename v, data)
	      | S.E_ExtractFemItem(v, ty, data) => S.E_ExtractFemItem(rename v, ty, data)
	      | S.E_ExtractFemItem2(v1, v2, ty, outTy, data) => S.E_ExtractFemItem2(rename v1, rename v2, ty, outTy, data)
	      | S.E_LoadFem(data, v1, v2) => S.E_LoadFem(data, rename v1, rename v2)
	      | S.E_FemField(v1, v1', v2, ty, field, func) => S.E_FemField(rename v1,rename v1', Option.map (rename) v2, ty, field, func)
	      | S.E_ExtractFemItemN(vars, tys, outTy, opt, NONE) => S.E_ExtractFemItemN(List.map rename vars, tys, outTy, opt, NONE) 
	   (* end case *))
      in
       S.Func{f=f, params=params',body=doBlock(body)}
      end

  (*Property and function to record order of functions*)
  local
   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn v => [] : F.t list)
   fun isolate [] = []
     | isolate (l as x::xs) =
       let fun remove (x,[]) = []
             | remove (x,l as y::ys) = if F.same(x,y)
                                       then remove(x,ys)
                                       else y::remove(x,ys)
       in
        x::isolate(remove(x,xs))
       end
   fun doStm(f, S.S_Var(_, SOME(S.E_Apply(f', _)))) = let
    val fs = getFn(f')
   in
    setFn(f, isolate(f'::fs))
   end
     | doStm (f, S.S_Assign(_, S.E_Apply(f', _))) = let
      val fs = getFn(f')
     in
      setFn(f, isolate(f'::fs))
     end
     | doStm (f, S.S_IfThenElse(_, b1, b2)) = (doBlock(f, b1); doBlock(f, b2))
     | doStm (f, S.S_Foreach(_,_, b)) = doBlock(f, b)
     | doStm _ = ()
   and doBlock(f, S.Block{code,...}) = List.app (fn x => doStm(f,x)) code
   fun defToName(S.Func{f,...}) = f

   fun getAllDeps(f : F.t) =
       let
	val involved =  getFn f
	val rest = isolate(List.concatMap getAllDeps involved)
       in
	involved@rest
       end

   fun orderFunction(def1, def2) =
       let
	val f1 = defToName(def1)
	val f2 = defToName(def2)
	val deps1 = getAllDeps(f1)
	val deps2 = getAllDeps(f2)
	val f1ltf2 = List.exists (fn x => F.same(x, f1)) deps2
	val f2ltf1 = List.exists (fn x => F.same(x, f2)) deps1

       in
	Bool.not(if f1ltf2
	then true
	else false)
       end

  in
  val defToName = defToName
  fun procFuncDef (S.Func{f, body,...}) = doBlock(f, body) (* apply this to defs *)
  fun sortFuncDef defs =  ListMergeSort.sort orderFunction defs
  end


  (*Data structure and helper functions for symbolically reasoning about fem data*)
  datatype femPres = Base of FD.femType * V.t option | ALL of FD.femType | NOTHING (* all, nothing, or one thing *)
		     | Tuple of femPres list | Array of femPres 
		     | Seq of femPres (* must be homogenous i.e no partitioning on an array*)

  fun presToString(Base(f, NONE)) = "Base(" ^ (FD.toString f) ^ ")"
    | presToString(Base(f, SOME(v))) =  "Base(" ^ (FD.toString f) ^ ", " ^ (V.uniqueNameOf v) ^ ")"
    | presToString (ALL f) = "All( " ^ (FD.toString f) ^ ")"
    | presToString (NOTHING) = "nothing"
    | presToString (Tuple ts) = "(" ^ (String.concatWith "," (List.map presToString ts)) ^ ")"
    | presToString (Seq(t)) = "Seq(" ^ (presToString t) ^ ")"
    | presToString (Array(t)) = "Array(" ^ (presToString t) ^ ")"

  fun valid' (ALL(f)) = FD.baseFem f
    | valid' (Base(f, NONE)) = FD.baseFem f
    | valid' (Base(f, SOME(v))) = FD.baseFem f andalso
				 (case V.typeOf v
				   of Ty.T_Fem(f') => FD.same(f, f')
				    | _ => raise Fail "invalid fempress: wrong var type for a base")
				 andalso (case V.kindOf v
					   of Var.InputVar => true
					    | _ => raise Fail "invalid fempress: var is not input")
    | valid' (Tuple(ts)) = List.all valid' ts
    | valid' (Seq(t)) = valid' t
    | valid' (Array(t)) = valid' t
    | valid' (NOTHING) = true

  fun valid d = if valid' d
		then d
		else raise Fail "invalid femPres"

  fun anyFem (Base _) = true
    | anyFem (ALL _) = true
    | anyFem (NOTHING) = false
    | anyFem (Tuple(ps)) = List.exists anyFem ps
    | anyFem (Array p) = anyFem p
    | anyFem (Seq p) = anyFem p


  fun hash (NOTHING) = 0w1
    | hash (ALL(f)) = 0w2 * FD.hash f
    | hash (Base(f, NONE)) = 0w3
    | hash (Base(f, SOME(v))) = 0w5 * (FD.hash f) + 0w7 * (V.hash v)
    | hash (Array(g)) = 0w11 * (hash g)
    | hash (Seq(g)) = 0w13 * (hash g)
    | hash (Tuple(gs)) = List.foldl (fn (ty, s) => hash ty + s) 0w17 gs

  fun same (NOTHING, NOTHING) = true
    | same (ALL(f1), ALL(f2)) = FD.same(f1, f2)
    | same (Base(f1, NONE), Base(f2, NONE)) = FD.same(f1, f2)
    | same (Base(f1, SOME(v1)), Base(f2, SOME(v2))) = FD.same(f1, f2) andalso V.same(v1, v2)
    | same (Array(p1), Array(p2)) = same(p1, p2)
    | same (Seq(p1), Seq(p2)) = same(p1, p2)
    | same (Tuple(ps1), Tuple(ps2)) = ListPair.all same (ps1, ps2)
    | same _ = false

  fun merge(a1 : femPres, a2 : femPres, change : bool ref) =
      let
       fun changed v = (change := true; v)
       fun failDiff(f1, f2) = raise Fail ("impossible Fem types different: " ^ (FD.toString f1) ^ (FD.toString f2))
       fun merge'(t1,t2) = (case (t1, t2)
		   of (Base(f1, SOME(v1)), Base(f2, SOME(v2))) => if FemData.same(f1, f2)
								  then if V.same(v1, v2)
								       then Base(f1, SOME(v1))
								       else changed (ALL(f1))
								  else raise failDiff(f1, f2)
		    | (Base(f1, SOME(v1)), Base(f2, NONE)) => if FemData.same(f1, f2)
							      then changed(Base(f1, SOME(v1)))
							      else failDiff(f1, f2)
		    | (Base(f1, NONE), Base(f2, SOME(v1))) => if FemData.same(f1, f2)
							      then changed(Base(f1, SOME(v1)))
							      else failDiff(f1, f2)
		    | (Base(_, NONE), Base(_, NONE)) => raise Fail "can't merge unit femPres; should be replaced"
		    | (NOTHING, NOTHING) => NOTHING
		    | (Base _, ALL(f)) => changed(ALL(f))
		    | (ALL(f), Base _) => changed(ALL(f))
		    | (ALL(f1), ALL(f2)) => if FD.same(f1, f2)
					    then ALL(f1)
					    else failDiff(f1, f2)
		    (*end possible base cases*)
		    | (Tuple(tys1), Tuple(tys2)) => Tuple(ListPair.map merge' (tys1, tys2))
		    | (Array(ty1), Array(ty2)) => Array(merge' (ty1, ty2))
		    | (Seq(t1), Seq(t2)) => Seq(merge'(t1, t2))
		    | _ => raise Fail ("impossible: merging incompatible femPress - must of come from different types: (" ^ (presToString t1) ^ ")" ^ " vs. (" ^ (presToString t2) ^ ")\n")
		 (*end case*))
      in
       valid (merge'(a1, a2)) handle exn => raise exn
      end

  fun merge'(a1, a2) = (merge(a1, a2, ref false)) handle exn => raise exn

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
    | tyToFemPres (Ty.T_Image _) = NOTHING
    | tyToFemPres (Ty.T_Kernel) = NOTHING

  fun tyToFemPresSeq onlyOne (t) = let
  fun tyToFemPresSeq' (Ty.T_Bool) = NOTHING
    | tyToFemPresSeq' (Ty.T_Int) = NOTHING
    | tyToFemPresSeq' (Ty.T_Image _) = NOTHING
    | tyToFemPresSeq' (Ty.T_Kernel) = NOTHING
    | tyToFemPresSeq' (Ty.T_String) = NOTHING
    | tyToFemPresSeq' (Ty.T_Tensor _) = NOTHING
    | tyToFemPresSeq' (Ty.T_Sequence(t, SOME k)) = Array(tyToFemPresSeq' t)
    | tyToFemPresSeq' (Ty.T_Sequence(t, NONE)) = Seq(tyToFemPresSeq' t)
    | tyToFemPresSeq' (Ty.T_Tuple(tys)) = Tuple(List.map tyToFemPresSeq' tys)
    | tyToFemPresSeq' (Ty.T_Strand _) = NOTHING
    | tyToFemPresSeq' (Ty.T_Field _) = NOTHING
    | tyToFemPresSeq' (Ty.T_Fem(f)) =
      if FD.baseFem f
      then (case onlyOne f
	     of SOME(v) => Base(f, SOME(v))
	      | _ => ALL(f)
	   (*end case*)) (* handle empty seq case *)
      else
       (case onlyOne (Option.valOf(FD.dependencyOf(f)))
	     of SOME(v) => Base(Option.valOf(FD.dependencyOf(f)), SOME(v)) 
	      | _ =>  Base(Option.valOf(FD.dependencyOf(f)), NONE) 
       (* end case*))
       handle ex => raise ex
      
  in
   tyToFemPresSeq' t handle ex => raise ex
  end

  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => false)
  in
  val isUsedFemData = getFn
  val useFemData = fn x => (setFn(x, true))
  end
  (*Propery for mapping vars to their symbolic fem info; inputs handled specially*)
  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => tyToFemPres(V.typeOf(v)))
  in
  fun updateFemPres(v : V.t, p : femPres) = (case peekFn(v)
					      of SOME(p') => (let val change = ref false
								  val p'' = merge(p, p', change) handle ex => raise ex
							      in setFn(v, p''); !change end)
					       | NONE => (setFn(v, p); true)
					    (*end case*))
  fun updateFemPresRef(v : V.t, p : femPres, changed) = (case peekFn(v)
							  of SOME(p') => (let val p'' = merge(p, p', changed) handle ex => (print("unable to merge:"^(V.uniqueNameOf v)^"\n"); raise ex)
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
			       of NONE => raise Fail (msg ^ " var " ^ (V.uniqueNameOf v) ^ " not defined")
				| SOME _ => ()
			     (*end case*))

  fun getFemPres x = let
   val r =  peekFn x
   val rr =(Option.valOf r)
  in
   rr
  end handle ex => raise ex
  val defaultFemPres = (fn x => (getFn x; ()))
  fun testParam x = (case getFemPres x
		      of Base(_, SOME _) => (case V.typeOf x
					      of Ty.T_Fem(d) => Bool.not(FD.baseFem d)
					       | _ => true
					    (* end case *))
		       | _ => true
		    (* end case *))
  fun getPropString v =
      (case peekFn v
	of NONE => "error"
	 | SOME(t) => presToString(t)
      (* end case*))
  end

  (* function to handle femField case where fictional vars are needed to handle the accompany function*)
  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => (v,v))
  in
  fun defineInOut(v, (ficIn, ficOut)) = setFn(v, (ficIn, ficOut))
  fun getInOut v = (case peekFn v
		     of SOME((a,b)) => (a,b)
		      | NONE => raise Fail "no In Out set"
		   (* end case*))
  val testInOut = peekFn
		     
  end

  (*Function to manage inputs*)
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

  (*data structures and functions to manage symbolic manipulation of functions at various callsites:*)
  type callSite =  V.t list  (* structure for idying a call site *)
  type functionCall  = femPres list * femPres list (* arguments and globals: structure for noting fem info at callsite *)
  type callResult = S.func_def * femPres option (* management for result of a callsite*)
  local
   fun sameSite(vs1, vs2) = ListPair.all (V.same) (vs1, vs2)
   fun hashSite(vs) = List.foldl (fn (v, s) => V.hash v + s) 0w3 vs
   fun sameCall((args1, globs1), (args2, globs2)) = ListPair.all same (args1, args2)
						    andalso ListPair.all same (globs1, globs2)
   fun hashCall((args1, globs1)) = (List.foldl (fn (v, s) => 0w5 * hash v + s) 0w3 args1) +
				   (List.foldl (fn (v, s) => 0w7 * hash v + s) 0w3 globs1)

  fun mkCallSiteTable(j : int) : (callSite, functionCall) HT.hash_table =
      HT.mkTable (hashSite, sameSite) (j,  Fail "missing call site")
  fun mkFuncDefTable(j : int) : (functionCall, callResult) HT.hash_table =
      HT.mkTable (hashCall, sameCall) (j,  Fail "missing function call data")

   fun mkFuncVarTable(j : int) : (functionCall, F.t) HT.hash_table =
       HT.mkTable (hashCall, sameCall) (j,  Fail "missing function call cvt")

   val cvtTbl = (F.newProp(fn f => mkFuncVarTable 16))
   val getCvtTbl = #getFn (cvtTbl)


   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn f => mkCallSiteTable 16)
   val defTblProp = F.newProp(fn f => mkFuncDefTable 16)
   val getFnTbl = #getFn defTblProp

   val needsCopyFlag = F.newFlag()

												
  in
  (*operation:we need to from a f, 
    and info about the calls (args +glob femPress) -> SOME(function_def)	
    (SOME if we need to do something and NONE if it is already done -> how to get ret var then?)  -> function_def can be processed *)
  (*from the procced thing, we get the accepted one, which we get back here.*)
  fun testFunction(args, globs, retTy) =
      let
       val femsInvolvedArgs = List.concatMap (Ty.allFems o V.typeOf) (args)
       val femsInvolvedGlobs = List.concatMap (Ty.allFems o V.typeOf) (globs)
       val anyFem = 0 <> (List.length (Ty.allFems retTy))
       val anyNonBase = (List.length femsInvolvedArgs <> 0)
			orelse (List.length femsInvolvedGlobs <> 0)
			orelse anyFem
      in
       anyNonBase
      end

  fun markForCopy(f : F.t) = (#setFn needsCopyFlag)(f, true)
  fun markForNonCopy(f : F.t) = (#setFn needsCopyFlag)(f, false)
  fun cleanBody(f : F.t) = Bool.not((#getFn needsCopyFlag)(f))

  
  fun addCall(f : F.t, call : callSite, callArgs : functionCall) : callResult =
      let
       val callTable = getFn f handle exn => raise exn
       val defsTable = getFnTbl f

       val callLookup : functionCall option = HT.find callTable call handle ex => raise ex
       val callExists = (case callLookup
		of SOME(callArgs') => sameCall(callArgs', callArgs)
		 | NONE => (HT.insert callTable (call, callArgs); false) handle exn => raise exn
			(* end case*))
      in
       if callExists
       then HT.lookup defsTable callArgs
	    handle exn => raise exn
       else let
	val _ = HT.insert callTable (call, callArgs) handle exn => raise exn
	val SOME(oldDef) = getfuncDef f
	val copyF = functionCopy(oldDef)
	val ret = (copyF, NONE)
			    (*check if the function with these callArgs exists anyway*)
       in (case HT.find defsTable callArgs handle exn => raise exn
	    of SOME(r) => r
	     | NONE => (HT.insert defsTable (callArgs, ret); ret) handle exn => raise exn
	  (* end case *))
       end
      end handle exn => raise exn

  fun updateCall(f : F.t, call : callSite, callArgs : functionCall,
		 newDef : S.func_def, retF : femPres) : unit =
      let
       val defsTable = getFnTbl f
      in
       HT.insert defsTable (callArgs, (newDef, SOME(retF)))
      end
  fun cvtFun(fvar : F.t, args : callSite, retVar : V.t,
	     newDefs : S.func_def list ref,
	     cvtBody : S.block * S.func_def list ref-> S.block,
	     cvtVar : V.t -> V.t, reductionOverride) =
      let
       val argsp = List.map getFemPres args
       val globs = getFuncGlobs fvar
       val globsp = List.map getFemPres globs
       val callArgs = (argsp, globsp)
       val anyNonBase = testFunction(args, globs, F.resultTypeOf fvar) orelse reductionOverride
      in
       if Bool.not anyNonBase
       then (markForNonCopy(fvar);NONE )
       else
	let
	 val _ = markForCopy(fvar)
	 val defsTable = getFnTbl fvar
	 val someDef = HT.find defsTable callArgs
	 val (def, retV) = (case someDef
			     of SOME(def, SOME(v)) => (def, v)
			      | NONE => raise Fail "impossible non recorded function"
			   (* end case*))
	 val retPres = retV
	in
	 (case HT.find (getCvtTbl fvar) callArgs
	   of SOME(f') => (markForCopy f'; F.use f'; F.decCnt fvar; SOME(f', retPres))
	    | NONE =>
	      let
	       val retTy' = V.typeOf(cvtVar retVar)
	       val S.Func{f, params, body} = def
	       val filterParams = List.filter testParam params
	       val params' = List.map cvtVar filterParams
	       val paramsTys = List.map V.typeOf params'
	       val body' = cvtBody(body, newDefs)
	       val f' = F.new(F.uniqueNameOf fvar, retTy', paramsTys)
	       val def' = S.Func{f=f', params=params', body=body'}
	       val _ = newDefs := (def' :: !newDefs)
	       val _ = HT.insert (getCvtTbl fvar) (callArgs, f')
	      in
	       (markForCopy f'; F.use f'; F.decCnt fvar; SOME(f', retPres))
	      end
	 (* end case*))
	end
      end
  end (* end local *)
  
  (*function to actually analyze a block of code*)
  fun analyzeBlock(b as S.Block{code, ...},
		   changeOuter : bool ref,
		   newsOuter : (Atom.atom * V.t list) list ref,
		   getDep : V.t -> V.t option,
		   addDep : V.t * V.t -> unit,
		   onlyOne : FemData.femType -> V.t option) =
      let
       fun expToFemPres(exp : S.exp, retTy : Ty.ty, call : callSite) : femPres  =
	   let
	    val change = ref false (* changes in here don't matter?*)
	    val merge = (fn (x,y) => merge(x,y,change)) handle ex => raise ex
	    (*These three are the only place where retTy is used:*)
	    val retHasFem = Ty.hasFem retTy
	    val nonFemRet = tyToFemPres retTy
	    val possibleFemData = (case retTy
				    of Ty.T_Fem(f') => SOME(f')
				     | _ => NONE)
	    val mergeFemPreses = List.foldr (fn (x,y) => merge(x, y)) handle ex => raise ex
	    fun doit (S.E_Var(v)) = (checkExistence ("_var_" ^ (V.nameOf v)) v; getFemPres v)
	      | doit (S.E_Lit(l)) = NOTHING
	      | doit (S.E_Kernel(k)) = NOTHING
	      | doit (S.E_Select(v1, v2)) = (
	       checkExistence ("_strand_" ^ (V.nameOf v1)) v1;
	       checkExistence ("_strand_var_" ^ (V.nameOf v2)) v2;
	       getFemPres v2 (* read strand state *) handle ex => raise ex
	      )
	      | doit (S.E_Apply(f, vs)) = (
	       List.app (checkExistence ((F.nameOf f) ^ "_arg")) vs;
	       procApply(f, vs, call, false) handle ex => raise ex; (*these are the only places where call is used.*)
	       getFemPres (List.hd call) handle ex => raise ex
	      )
	      | doit (S.E_Prim(v, _, vs, t)) =
		(List.app (checkExistence ("_arg")) vs;
		 if Var.same(v, BasisVars.subscript) andalso retHasFem
		 then let val [sy, index] = vs (*this merge might be pointless*)
			  val Array(p') = getFemPres sy handle ex => raise ex
		      in p' end
		 else if Var.same(v, BasisVars.dynSubscript) andalso retHasFem
		 then let val [sy, index] = vs
			  val Seq(p) = getFemPres sy handle ex => raise ex
		      in
		       p
		      end
		 else if Var.same(v, BasisVars.at_dT) andalso retHasFem
		 then let val [sy, a] = vs
			  val Seq(p) = getFemPres sy
			  val p' = getFemPres a handle ex => raise ex
			  val p'' = merge(p, p') handle ex => raise ex
		      in Seq(p'') end
		 else if Var.same(v, BasisVars.at_Td) andalso retHasFem
		 then let val [a, sy] = vs
			  val Seq(p) = getFemPres sy handle ex => raise ex
			  val p' = getFemPres a
			  val p'' = merge(p, p') handle ex => raise ex
		      in Seq(p'') end
		 else if Var.same(v, BasisVars.at_dd) andalso retHasFem
		 then let val [sy1, sy2] = vs
			  val Seq(p1) = getFemPres sy1 handle ex => raise ex
			  val Seq(p2) = getFemPres sy2 handle ex => raise ex
			  val p = merge(p1, p2) handle ex => raise ex
		      in Seq(p) end
		 else if retHasFem
		 then raise Fail "impossible fem basis ret"
		 else nonFemRet 
		)
	      | doit (S.E_Tensor(vs, _)) = (List.app (checkExistence ( "_ten")) vs; NOTHING)
	      | doit (S.E_Field(vs, _)) = (List.app (checkExistence ("_fld")) vs; NOTHING)
	      | doit (S.E_Seq(vs, t)) =
		(List.app (checkExistence ("_seq")) vs;
		 (case t
		   of Ty.T_Sequence(t', NONE) => if List.length vs <> 0
						 then Seq(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs)) handle ex => raise ex
						 else Seq(tyToFemPresSeq onlyOne t') (*{} -> ALL?*)
		    | Ty.T_Sequence(t', SOME(k)) => if k <> 0
						    then Array(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs)) handle ex => raise ex
						    else Array(tyToFemPresSeq onlyOne t') (* {} -> ALL?*)
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
			       then tyToFemPresSeq onlyOne ty  (*QUESTION: this is a bit weird as we have the p' - for {}, p' would already be all so this is redudant*)
			       else merge((tyToFemPres ty), p') handle ex => raise ex (*if there is a FEM, there is an init so Base(x, SOME k) gets used*)
		      in
		       Seq(p)
		      end
		    | _ => tyToFemPres dstTy (*kernel to kernel, int to tensor, same; tuple handled earlier*)
		(*end case*)))
	      | doit (S.E_BorderCtl(b, v)) = ((checkExistence ( "_b")) v; NOTHING)
	      | doit (S.E_LoadSeq _) = NOTHING
	      | doit (S.E_LoadImage _) = NOTHING
	      | doit (S.E_LoadFem(f, v1, v2)) =
		(List.app (checkExistence ("_loadFem")) [v1, v2];
		 if FD.baseFem f (*v1 is the dest data; v2 is the dep; *)
		 then ((case getFemPres v2
			 of Base(f', SOME(v2')) => (useFemData(v1);addDep(v1, v2'); valid (Base(f, SOME(v1))) handle exn => raise exn) 
			  | Base(f', NONE) => raise Fail "base fem load not set; impossible!" 
			  | ALL(f') => raise Fail "base fem load unable to trace!" (*Should be impossible*)
			  | _ => raise Fail ("impossible: " ^ (V.uniqueNameOf v2))
		      (* end case*)))
		 else let val g = (Option.valOf o FD.dependencyOf) f handle ex => raise ex in
		       (valid(getFemPres v1))
		       (* (case getDep v1 (*v1 is the dep, v2 is an int*) *)
		       (* 	 of SOME(v1')=> valid (Base(g, SOME(v1')))  *)
		       (* 	  | NONE => valid((ALL(g))) *)
		       (* (*end case*)) *)
		      end)
	      (*v1 is base data (func or mesh); v2 is an int for a cell?*)
	      | doit (S.E_ExtractFem(v, f)) = ( (*if the return has fem, this is getting a func or space or mesh*)
	       let
		val _ = checkExistence "_extractFem" v
		val retFem = (case possibleFemData
			       of SOME(f') => f'
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
						     of (SOME(f'), SOME(v'')) => valid (Base(f',SOME(v'')))
						      | (SOME(f'), NONE) => valid (ALL(f') )
						      | (NONE, _) => raise Fail "impossible mesh extract in simple"
						   (* end case*))
			| Base(vtyf', NONE) => raise Fail "extracting from non init"
			| ALL(vtyf') => valid (ALL(vtyf'))
			| _ => raise Fail "impossible extract FEM"
		     (* end case*))
		else getFemPres v (*take base out of cell or pos; just propogate same info*)
	       end)
	      | doit (S.E_ExtractFemItem(v, t, fo)) = (checkExistence "_EFI1" v;
						       if retHasFem
						       then raise Fail "impossible"
						       else nonFemRet (*should not have Base/ALL*))
	      (*above could be:c cells, refCell -> but these should be eliminated in prior phases*)
	      | doit (S.E_ExtractFemItem2(v1, v2, t1, t2, fo)) = (checkExistence "_EFI2_1" v1;
								  checkExistence "_EFI2_2" v2;
								  if retHasFem
								  then raise Fail "impossible"
								  else nonFemRet)
	      | doit (S.E_FemField(v1, v2, v3o, t, fof, func)) =
		(checkExistence "_FField_1" v1;
		 checkExistence "_FField_2" v2;
		 Option.app (checkExistence "_FField_3") v3o;
		 let
		  val Ty.T_Fem(femData) = V.typeOf v1
		  val dim  = FD.underlyingDim femData
		  val fictionalEvalTy = Ty.vecTy' dim
		  val fictionalEvalVar = V.new("fiction", Var.LocalVar, fictionalEvalTy)
		  val _ = defaultFemPres fictionalEvalVar

		  val outVarTy = Option.map F.resultTypeOf func
		  val outVarFiction = Option.map
					(fn x => [V.new("fiction", Var.LocalVar, x)]) outVarTy

		  val (fictionalEvalVar, outVarFiction) = (case testInOut (List.hd call)
							    of NONE => (fictionalEvalVar, outVarFiction)
							     | SOME(ficEvalVar, ficOutVar) => (ficEvalVar, SOME([ficOutVar]))
							  (*end case*))
					
		 in
		  ((case (func, fof, v3o, Option.map (V.typeOf) v3o)
		    of (NONE, _, _, _) => ()
		     | (SOME(f), FemOpt.RefField, SOME(v3), SOME(Ty.T_Fem(mesh))) =>
		       (defineInOut(List.hd call, (fictionalEvalVar, List.hd (Option.valOf outVarFiction)));
			procApply(f, [v3, fictionalEvalVar],
				  Option.valOf outVarFiction,
				  false))
		     | (SOME(f), FemOpt.InvTransform, SOME(v3),SOME(Ty.T_Fem(mesh))) =>
		       (defineInOut(List.hd call, (fictionalEvalVar, List.hd (Option.valOf outVarFiction)));
			procApply(f, [v3, fictionalEvalVar],
				  Option.valOf outVarFiction,
				  false))
		     | (SOME(f), FemOpt.InvTransform, SOME(v3), SOME(Ty.T_Int)) =>
		       (defineInOut(List.hd call, (fictionalEvalVar, List.hd (Option.valOf outVarFiction)));
			procApply(f, [fictionalEvalVar, v3, v1],
				  Option.valOf outVarFiction,
				  false)
			)
		     | _ => raise Fail "impossible"
		   (*end case*)))
		 end;
		 (*NOTE: function here is FEM free modulo a mesh arg!*)
		 nonFemRet)
	      | doit (S.E_ExtractFemItemN(vs, tys, t, fo, NONE)) = ((List.app (checkExistence "_FFIN") vs);
								    if Bool.not retHasFem
								    then nonFemRet
								    else
								     (*refBuild, InvalidBuild, InvalidbuildBoundary, AllBuild  *)
								     let
								      val (opt, data) = fo
								      val data' = (case FD.dependencyOf data
										   of SOME(d) => d
										    | _ => data
										  (* end case*))
								      val n = List.length tys
								      val tabed = List.tabulate(n, fn x => (x, List.nth(tys, x)))
								      val pos  = List.filter
										   (fn (_, Ty.T_Fem(f')) => FD.same(data', f')
										   | _ => false) tabed
								      val poses = "[" ^ (String.concatWithMap "," Ty.toString (tys)) ^ "]"
								     in
								      (case pos
									of [(idx, _)] => getFemPres (List.nth(vs, idx)) (*propogate Base/ALL*)
									 | _ => raise Fail ("not planned in extractFemItemN:"^poses)
								      (*end case*))
								     end)
								     
	      | doit (S.E_InsideImage(v1, v2, _)) = (List.app (checkExistence ("_in")) [v1, v2]; NOTHING)
	      | doit (S.E_FieldFn _) = NOTHING (* FIXME:handle function *)

	   in
	    doit(exp) handle ex => raise ex
	   end
       and procStatement(stm : S.stmt, changed : bool ref,
			 call : callSite,
			 news : (Atom.atom * V.t list) list ref) =
	   let
	    fun doit (S.S_Var(v, NONE)) = () (*Question: should we do anything: Answer: NO.*)
	      | doit (S.S_Var(v, SOME(e))) = ((* print("var:"^(V.uniqueNameOf v) ^"\n"); *) updateFemPresRef(v, expToFemPres(e, V.typeOf v, v::call), changed))
	      | doit (S.S_Assign(v, e)) = ((* print("assign:"^(V.uniqueNameOf v) ^"\n"); *)updateFemPresRef(v, expToFemPres(e, V.typeOf v, v::call), changed) )
	      | doit (S.S_IfThenElse(v, b1, b2)) = (checkExistence "condition" v; doBlock'(b1); doBlock'(b2))
	      | doit (S.S_New(a,vs)) = (news := (a, vs) :: !news)
	      | doit (S.S_Foreach(itter, src, blk)) =
		let
		 val _ = checkExistence "itterSrc" src
		 val itterPres = (case getFemPres src
				   of Seq(a) => a
				    | Array(a) => a
				    | _ => raise Fail "impossible: itter with wrong seq type"
				 (* end case*))
		in
		(updateFemPresRef(itter, itterPres, changed);
		checkExistence "itter" itter;
		 doBlock'(blk))
		end
	      | doit (S.S_KillAll) = ()
	      | doit (S.S_StabilizeAll) = ()
	      | doit (S.S_Continue) = ()
	      | doit (S.S_Die) = ()
	      | doit (S.S_Stabilize) = ()
	      | doit (S.S_Return(site)) = let
	       val sitePres = getFemPres site
	       val retVar = (case call
			      of v::_ => v
			       | [] => raise Fail "return outside of function!"
			    (* end case*))
	      in
	       updateFemPresRef(retVar, sitePres, changed) (*NOTE: preempts the last step of processing apply*)
	      end
	      | doit (S.S_Print(vs)) = List.app (checkExistence "print") vs
	      | doit (S.S_MapReduce(maps)) = let
	       fun doReduce(S.MapReduce{result, reduction, mapf, args, source, domain}) =
		   let
		    fun markStrand x = (case V.typeOf x
					 of Ty.T_Strand _ => (updateFemPres(x, NOTHING); ())
					  | _ => () (*end case*))
		    val S.Func{f, ...} = mapf
		    val _ = registerFunctions([mapf])
		    val _ = registerGlobVarsInFunc mapf
		    val _ = List.app markStrand args
		    val _ = List.app (checkExistence "_map_args_") args
		    (*QUESTION: val _ = (checkExistence "_strand") source: we generaly ingore the strand type *)
		    (*NOTE: args are just args to the function*)
		    (*NOTE: reduction and domain can be ignored *)
		    (*NOTE: This is a case of an apply where args are the params and result is the callsite *)
		    val _ = procApply(f, args, [result], true)
		   in
		    ()
		   end
	      in
	       List.app doReduce maps
	      end
	    and doBlock'(S.Block{code, ...}) = List.app doit ( code)
	   in
	    doit(stm)
	   end
       and doBlock(S.Block{code, ...}, changed : bool ref,
		   ret : callSite,
		   news : (Atom.atom * V.t list) list ref) =
	   let
	    val g =  (fn s => procStatement(s, changed, ret, news))
	   in
	    List.app g (code)
	   end
       and procApply(f : F.t, args : V.t list,
		     call : callSite, reductionOverride) = let

	val argsp = List.map getFemPres args
	val globs = getFuncGlobs f
	val globsp = List.map getFemPres globs
	val callArgs : functionCall = (argsp, globsp)
	val anyNonBase = testFunction(args, globs, F.resultTypeOf f)
	val v : V.t = List.hd call
       in
	if anyNonBase orelse reductionOverride
	then
	 let
	  val (fdef, possibleRet) = addCall(f, call, callArgs) 
	 in
	(case possibleRet
	  of SOME(p) => (case call
			  of v ::_ => (updateFemPresRef(v, p, changeOuter)  handle ex => raise ex)
			   | _ => raise Fail "call not in call"
			(*end case*))
	   | NONE =>
	     let
	      (* in this path, we can't find an already inlined version of this function so there was a change.*)
	      val _ = changeOuter := true 
	      val S.Func{params, body, ...} = fdef
	      val _ = ListPair.map updateFemPres (params, argsp)
	      val _ = doBlock(body, ref false, call, ref []) handle ex => raise ex
	      (*we can ignore changes/new here as we found them/they are impossible.*)

	      val femPres = getFemPres v
	     in
	      updateCall(f, call, callArgs, fdef, femPres) : unit
	     end
	(* end case*))
	 end
	else (updateFemPresRef(v, tyToFemPresSeq onlyOne (V.typeOf v), changeOuter)) (*no fem involved, but we mark this.*)
       end
      in
       doBlock(b, changeOuter, [], newsOuter)
      end
										      

  fun analysis prog = let
   val S.Program{props, consts, inputs, constInit, globals,
		 globInit, funcs, strand, create, start, update} = prog

   val _ = registerFunctions(funcs)
   val _ = List.app registerGlobVarsInFunc funcs
   (*sml docs said: explain the use and semantics of mkTable HERE.*)
   val femTable : (FD.femType, V.t list) HT.hash_table = HT.mkTable(FD.hash, FD.same) (32, Fail "femtype not added yet")
   fun addType (ty, v) = (case HT.find femTable ty
			   of SOME(vs) => HT.insert femTable (ty, v::vs)
			    | NONE => HT.insert femTable (ty, [v])
			 (*end case*))
   fun onlyOne (ty) = (case HT.find femTable ty
			of SOME([v]) => SOME(v)
			 | _ => NONE)
   val femDep = VTbl.mkTable(32, Fail "fem dep not found")
   fun addDep (v1 : V.t, v2 : V.t) = (VTbl.insert femDep (v1, v2))
   fun findDep ( v1 : V.t) : V.t option = VTbl.find femDep v1

   (*Figure out the base fem inputs.*)
   val _ = procBaseInputs(inputs, addType)

   (* register const inputs - no fem*)
   val _ = List.app defaultFemPres consts

   fun analyze(b, change, news) = analyzeBlock(b, change, news, findDep, addDep, onlyOne)

   (*Analyze the global init - add deps here.*)
   val _ = analyze(constInit, ref false, ref [])
   val _ = analyze(globInit, ref false, ref [])

   (*Analyze the creation block, grabing the news especially*)
   val Create.Create {dim, code} = create
   val news = ref []
   val _ = analyze(code, ref false, news)		  

   (* do the strand create: params,*)

   fun doStrandStart(strand : S.strand, news, change) =
       let
	val S.Strand{name, params, spatialDim,  state, stateInit, startM, updateM, stabilizeM} = strand
	val news' = !news

	val _ = (news := [])
	val paramArgs = List.map (fn (x,y) => y) news'

	val _ = List.map (fn p => ListPair.app (fn (pam, pres) => updateFemPresRef(pam, getFemPres pres, change)) (params, p)) paramArgs
	val _ = analyze(stateInit, change, news)
	val _ = Option.app (fn x => analyze(x, change, news)) startM
       in
	()
       end

   val _ = doStrandStart(strand, news, ref false)

   (*Do global start:*)
   val _ = Option.app (fn x => analyze(x, ref false, ref [])) start

   fun doStrands(strand, news, change) =
       let
	val _ = doStrandStart(strand, news, change)
	val S.Strand{name, params, spatialDim,  state, stateInit, startM, updateM, stabilizeM} = strand
	val _ = analyze(updateM, change, news)
	val _ = Option.app (fn x => analyze(x, change, news)) stabilizeM
       in
	()
       end

   fun doGlobal(change) = Option.app (fn x => analyze(x, change, ref [])) update

   (*first run:*)
   val  _ = doStrands(strand, news, ref false)
   val _ = doGlobal(ref false)

		    
   fun loop(change : bool ref, d) = if d > 500 then print( "way too much looping") (*debug*)
			 else
			  if Bool.not (!change)
			  then ()
			  else (change := false;
				doStrands(strand, news, change);
				doGlobal(change);
				loop(change, d+1))

   val change = ref true
   val _ = loop(change, 0)


(*
NOTE: peekFn doesn't create the thing -> getFn or setFn will though ->all vars need to be accounted for then.
NOTE: updateCall/those tables won't be around (peekFn => NONE) for functions with no fem/no new defs
NOTE: add single function/modify input/globals transformation (avoid {})
NOTE: think about {}s and def of rep (not SSA): comprehensions, {}s
*)
   val _ =()
	    (*Now do transform:
	    0. Two erros: map reduce has error and functions with no calls are filled with errors...
	    ---Reason is that we don't replace the function def.......could cause problems later.
	    ---Note: the no calls thing we might want to change: if we have a mesh1 and mesh2, we might comopare equality...
	    ---empy seq? Correct order?
	    ??
	    1. insert array of fems - mesh, space, func
	    2. start doing conversions of vars via properties and fem via extract fem - both depend on thing:
	    2.1. Base something -> just get the global and convert to the base one
	    2.2 ALL -> extra int and access into a thing
	    2.3 Conversion from base to all -> access into global (Is this possible?)
	    3. inverse, position, cancel, inside.
	     *)

  in
   (prog, femTable, femDep)
  end

  (*infustructure to do translation: TOMORROW: (do all up to ID field, then cancel/pos, then affine)
  1.we need five arrays to hold info: meshes, spaces, funcs, func -> space, space -> mesh
  2. We need to do a cvtVar (just convert the type) if needed
  3. cvt block -> cvt stm -> cvt exp + applies 
  (most things are just immediate - just convert
   next thing are load/extract -> pretty clear
   Watch out for conversions: seq merge - so we need a convert between preses -> where could this occur -> seqs merge
   ..
   temp
   for x in origin
   ret = covert{}, ..,{}
   end
   use ret for the rest...
  )
  4. Functions
  5. tomorrow we do cancel/identity field and pos and other stuff (maybe inside)
   *)

			
  fun buildGlobalFemTables(femTable : (FD.femType, V.t list) HT.hash_table, femDep : SimpleVar.t VTbl.hash_table) =
      let
       fun findDep ( v1 : V.t) : V.t option = VTbl.find femDep v1 handle ex => raise ex
       val initsAndVarTbl = HT.mapi (fn (d, vs) =>
				     let
				      val n = List.length vs
				      val ty = Ty.T_Sequence(Ty.T_Fem(d), SOME(n))
				      val glob = V.new("femGlobArray", Var.GlobalVar, ty)
				      val init = S.S_Assign(glob, S.E_Seq(vs, ty))
				     in
				      (glob, init)
				     end) femTable

       fun femItoGlob(femType, i) = (case HT.find femTable (femType)
				      of SOME(ls) => if (0 <= i) andalso (i <= ((List.length ls) - 1))
						     then SOME(List.nth(ls, i))
						     else NONE
				       | NONE => NONE
				    (* end case*))

       fun find' test vs=
	   let
	    val n = List.length vs
	    val tabed = List.tabulate(n, fn x => (x, List.nth(vs, x)))
	   in
	    Option.map (fn (x,y) => x) (List.find (fn (x,y) => test y) tabed)
	   end
       fun femGlobToI (femType, v) : int option = Option.mapPartial (find' (fn x => V.same(x,v))) (HT.find femTable femType)


				      
       (*each global is an array acting as map from Int to GlobalVar of this type.*)
       val (newGlobals : V.t list, newInit: S.stmt list) = ListPair.unzip (HT.listItems initsAndVarTbl)

       val counts = List.map List.length (HT.listItems femTable)
       val filteredTable = HT.copy(femTable)
       val _ =  (HT.filteri (Option.isSome o FD.dependencyOf o (fn (x,y) => x))) filteredTable
       val globDepTable = HT.mapi (fn (d, vs) =>
				      let
				       val d' = Option.valOf (FD.dependencyOf d)
				       val n = List.length vs
				       val _ = print(String.concatWith "," (List.map V.uniqueNameOf vs) )
				       val _ = print("\n")
				       val deps = List.map (fn v => VTbl.find femDep v) vs handle ex => raise ex
				       val depsFilt = List.map (Option.valOf) (List.filter Option.isSome deps)
				       val ty = Ty.T_Sequence(Ty.T_Fem(d'), SOME(List.length depsFilt))
				       val depNums = List.map (fn x =>
								  (case x
								    of SOME(x') => Option.valOf (femGlobToI(d', x'))
								     | NONE => (~1) (*NOTE: should be impossible to hit b/c no uses*)
							      (* end case*))) deps
								    
				       val ty' = Ty.T_Sequence(Ty.T_Int, SOME(n))


				       val globD = V.new("femGloblDepArray", Var.GlobalVar, ty)
				       val globI = V.new("femIntoToDepIntArray", Var.GlobalVar, ty')
				       val ises = List.tabulate(n, fn x => V.new("globFemI_"^(Int.toString x), Var.LocalVar, Ty.T_Int))
				       val isesAssign = ListPair.map (fn (x,y) => S.S_Var(x, SOME(S.E_Lit(Literal.intLit y)))) (ises, depNums)
				      in
				       ((globD, S.S_Assign(globD, S.E_Seq(depsFilt, ty))),
					(globI, isesAssign@[S.S_Var(globI, SOME(S.E_Seq(ises, ty')))]))
				      end
				  ) filteredTable
       val (depVars, depIdx) = ListPair.unzip (HT.listItems globDepTable)
       val ((newGlobals', newInit'), (newGlobals'', newInit'')) = (ListPair.unzip depVars, ListPair.unzip depIdx)
						      

       val newGlobals = newGlobals@newGlobals'@newGlobals''
       val newInits = newInit@newInit'@(List.concat newInit'')

       fun femItoDepIImpl(fem, dFem, v) = 
	   let
	    val depVar = HT.lookup globDepTable fem handle ex => raise ex
	    val ((_, _), (idxi, _)) = depVar
	    val accExp = S.E_Prim(BasisVars.subscript, [], [idxi, v], Ty.T_Int)
	    val accVar = V.new("depIndex", Var.LocalVar, Ty.T_Int)
	    val getI = S.S_Var(accVar, SOME(accExp))

	    val (depArray, _) = HT.lookup initsAndVarTbl dFem handle ex => raise ex
	    val accExp' = S.E_Prim(BasisVars.subscript, [], [depArray, accVar], Ty.T_Fem(dFem))
	    val accVar' = V.new("depFem", Var.LocalVar, Ty.T_Fem(dFem))
	    val getD = S.S_Var(accVar', SOME(accExp'))
	   in
	    (accVar', accVar, [getI, getD])
	   end
       fun globDepToDep(fem, var) =
	   let
	    val dep = Option.valOf (FD.dependencyOf fem)
	    val (_, depInt, stms) = femItoDepIImpl(fem, dep, var) handle ex => raise ex
	   in
	    (stms, depInt)
	   end

       fun intGlobToVar(fem, var) = (*given an idx, produce the glob assocaited to it*)
	   let
	    val (femTable, _) = HT.lookup initsAndVarTbl fem handle ex => raise ex
	    val accExp = S.E_Prim(BasisVars.subscript, [], [femTable, var], Ty.T_Fem(fem))
	    val accVar = V.new("depFem", Var.LocalVar, Ty.T_Fem(fem))
	    val getD = S.S_Var(accVar, SOME(accExp))
	   in
	    ([getD], accVar)
	   end
       fun femGlobToIImpl (fem, v) = (*given a global, find the idx associated to it.*)
	   let
	    val i = Option.valOf (femGlobToI(fem, v)) handle ex => raise ex
	    val lit = Literal.intLit i
	    val litExp = S.E_Lit(lit)
	    val litVar = V.new("femGlobI", Var.LocalVar, Ty.T_Int)
	    val litDec = S.S_Var(litVar, SOME(litExp))
	   in
	    (litVar, litDec)
	   end
      in
       (*femGlobToI, femItoGlob, femItoDepIImpl*)
       (newGlobals, newInits,femGlobToIImpl, findDep, globDepToDep, intGlobToVar)
      end

  fun femTy(d, isBase) =
      if FD.baseFem d
      then if isBase
	   then Ty.T_Fem(d)
	   else Ty.T_Int
      else (case d
	     of FD.MeshPos(m) =>
		if isBase
		then Ty.T_Tuple([Ty.vecTy' (FD.underlyingDim d), Ty.T_Int, Ty.T_Int])
		else Ty.T_Tuple([Ty.vecTy' (FD.underlyingDim d), Ty.T_Int, Ty.T_Int, Ty.T_Int])
	      | FD.FuncCell(_) => if isBase
				  then Ty.T_Int
				  else Ty.T_Tuple([Ty.T_Int, Ty.T_Int])
	      | FD.MeshCell(_) => if isBase
				  then Ty.T_Int
				  else Ty.T_Tuple([Ty.T_Int, Ty.T_Int])
	      | FD.RefCell _ => raise Fail "impossible"
	      | FD.Mesh _ => raise Fail "impossible"
	      | FD.Space _ => raise Fail "impossible"
	      | FD.Func _ => raise Fail "impossible"
	   (* end case*))

  fun cvtTy(ty, pres) = if Bool.not(Ty.hasFem ty orelse anyFem pres)
			then ty
			else
			 (case (ty, pres)
			   of (Ty.T_Bool, NOTHING) => ty
			    | (Ty.T_Int, NOTHING) => ty
			    | (Ty.T_String, NOTHING) => ty
			    | (Ty.T_Tensor _, NOTHING) => ty
			    | (Ty.T_Sequence(t', NONE), Seq(p)) => Ty.T_Sequence(cvtTy(t', p), NONE)
			    | (Ty.T_Sequence(t', SOME(k)), Array(p)) => Ty.T_Sequence(cvtTy(t', p), SOME(k))
			    | (Ty.T_Tuple(ts), Tuple(ps)) => Ty.T_Tuple(ListPair.map cvtTy (ts, ps))
			    | (Ty.T_Strand _, NOTHING) => ty
			    | (Ty.T_Kernel, NOTHING) => ty
			    | (Ty.T_Field _, NOTHING) => ty
			    | (Ty.T_Fem(d), Base(d', NONE)) => raise Fail "cvtTy: analysis marked base none"
			    | (Ty.T_Fem(d), Base(d', SOME(v))) => femTy(d, true)
			    | (Ty.T_Fem(d), ALL(d')) => femTy(d, false)
			    | _ => raise Fail ("impossible: " ^ (Ty.toString ty) ^ " vs. " ^ (presToString pres))
			 (* end case*))
  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => (case V.kindOf v
							   of Var.ConstVar => v
							    | Var.InputVar => v
							    | _ => V.new(V.nameOf v, V.kindOf v, cvtTy(V.typeOf v, getFemPres v))
							 (* end case*)))
  in
  val cvtVar = getFn
  end

  fun checkMerge(pres1, pres2) =
      let
       val ret = ref false
       val _ = merge(pres1, pres2, ret)
      in
       if same(pres1, pres2)
       then false
       else !ret
      end
  fun checkMerges(preses, retPres) = let
   val ret = ref false
   val m = fn (x,y) => merge(x, y, ret)
   val _ = List.foldr m retPres preses
  in
   !ret 
  end

  fun makeConversions(globCvt, globDep, globDepToDep, depToV) =
      let
       fun projectCvt v i = ([], S.E_Project(cvtVar v, i))
       fun toAll'(v : V.t, ty, pres, targetPres) = 
	   (case (ty, pres, targetPres)
	     of (_, NOTHING, NOTHING) => ([], v)
	      | (Ty.T_Sequence(ty', SOME(k)), Array(p), Array(p')) =>
		let
		 val accTy = cvtTy(ty', p) (*new type for these arrays*)
		 (*build index into array:*)
		 val accAddr = List.tabulate(k, fn x => V.new("arrayAccIdx_"^(Int.toString x), Var.LocalVar, Ty.T_Int))
		 val accAddrInit = List.tabulate(k, fn x => S.S_Var(List.nth(accAddr, k), SOME(S.E_Lit(Literal.intLit k))))

		 (*build acc of array:*)
		 fun acc k = S.E_Prim(BasisVars.subscript, [], [v, List.nth(accAddr, k)], ty')
		 val accVars = List.tabulate(k, fn x => V.new("arrayAcc_"^(Int.toString x), Var.LocalVar, ty'))
		 val accese = List.tabulate(k, fn x => S.S_Var(List.nth(accVars, k), SOME(acc k)))
		 (*convert each acc:*)
		 val rets = List.map (fn vv => toAll'(vv, ty', p, p')) accVars (*recursive call*)
		 (*resulting stms:*)
		 val newStms = List.concatMap (fn (x, y) => x) rets
		 val newVars = List.map (fn (x, y) => y) rets
		 (*seq replacement:*)
		 val newVar = V.new(V.nameOf v, V.kindOf v, Ty.T_Sequence(accTy, SOME(k)))
		 val ret = S.S_Var(newVar, SOME(S.E_Seq(newVars, Ty.T_Sequence(accTy, SOME(k)))))
		in
		 if checkMerge(p, p')
		 then (accAddrInit@accese@newStms@[ret], newVar)
		 else ([], v)
		end
	      | (Ty.T_Sequence(ty', NONE), Seq(p), Seq(p')) =>
		if checkMerge(p, p') 
		then ([], v)
		else let
		 val newTy = cvtTy(ty', p') (* new inner ty*)
		 val retSeqTy = Ty.T_Sequence(newTy, NONE)
					     
		 val destVar = V.new(V.nameOf v, V.kindOf v, newTy) (*seq we store in*)
		 val varStart = S.S_Var(destVar, SOME(S.E_Seq([], retSeqTy))) (*init for this seq*)
		 val itter = V.new("itter", V.IterVar, ty') 
		 val (innerStms, innerVar) = toAll'(itter, ty', p, p')
		 (*^create body of loop that converts inner, 
			store innerVar in seq:*)
		 val appendExp = S.E_Prim(BasisVars.at_dT, [], [destVar, innerVar], retSeqTy)
		 val appendStm = S.S_Assign(destVar, appendExp)
		in
		 (varStart :: S.S_Foreach(itter, v, S.Block{code=innerStms@[appendStm], props=PropList.newHolder()}) :: [],
		  destVar)
		end
	      | (Ty.T_Tuple(tys), Tuple(ps), Tuple(ps')) =>
		let
		 val k = List.length tys
		 val tys' = ListPair.map cvtTy (tys, ps')
		 val accVars = List.tabulate(k, fn x => V.new("tpl_"^(Int.toString x), Var.LocalVar, List.nth(tys, x)))
		 val accStms = List.tabulate(k, fn x => S.S_Var(List.nth(accVars, x), SOME(S.E_Project(v, x))))
		 val converts = List.tabulate(k, fn x => let val nth = fn l => List.nth(l, x) in 
							  toAll'(nth accVars, nth tys, nth ps, nth ps') end)
		 val convertStms = List.concatMap (fn (x,y) => x) converts
		 val convertVars = List.map (fn (x,y) => y) converts
		 val finVar = V.new(V.nameOf v, V.kindOf v, Ty.T_Tuple(tys'))
		 val finVarStm = S.S_Var(finVar, SOME(S.E_Tuple(convertVars)))
					
		in
		 if ListPair.exists checkMerge (ps, ps')
		 then (accStms@convertStms@[finVarStm], finVar)
		 else ([], v)
		end
	      | (Ty.T_Fem(d), Base(d', SOME v1), Base(d'', SOME v2)) => if FD.same(d', d'')
									then if V.same(v1, v2)
									     then ([], v)
									     else raise Fail "impossible convert: wrong base var"
									else raise Fail "impossible convert: wrong fems"
	      | (Ty.T_Fem(d), ALL(d'), ALL(d'')) => if FD.same(d', d'')
						    then ([], v)
						    else raise Fail "impossible convert:wrong fems"
	      | (Ty.T_Fem(d), Base(d', SOME v1), ALL(d'')) =>
		let
		 val dim = FD.underlyingDim d
		 val (newVar, newStm)  = globCvt(d', v1)
		 val newTy = cvtTy(ty, ALL(d'))
		 val newVar' = V.new(V.nameOf v, V.kindOf v, newTy)
		 val posVecVar = V.new("refPos", Var.LocalVar, Ty.vecTy' dim)
		 val posCellVar = V.new("cell", Var.LocalVar, Ty.T_Int)
		 val faceCellVar = V.new("cell", Var.LocalVar, Ty.T_Int)
					(*get a few accs*)
		in
		 if FD.same(d', d'')
		 then if FD.baseFem d
		      then ([newStm], newVar)
		      else (case d
			     of FD.MeshCell _ => ([newStm, S.S_Var(newVar', SOME(S.E_Tuple([v, newVar])))],
						  newVar')
			      | FD.FuncCell _ => ([newStm, S.S_Var(newVar', SOME(S.E_Tuple([v, newVar])))],
						  newVar')
			      | FD.MeshPos _ => ([newStm, S.S_Var(posVecVar, SOME(S.E_Project(v, 0))),
						  S.S_Var(posCellVar, SOME(S.E_Project(v, 1))),
						  S.S_Var(faceCellVar, SOME(S.E_Project(v, 2))),
						  S.S_Var(newVar', SOME(S.E_Tuple([posVecVar, posCellVar,faceCellVar, newVar])))],
						 newVar')
			   (* end case*))
		 else raise Fail "impossible convert: wrong fems"
		end
	      | (Ty.T_Fem(d), Base(_, NONE), _) => raise Fail "not allowed NONE"
	      | (Ty.T_Fem(d), _, Base(_, NONE)) => raise Fail "not allowed NONE"
	      | (Ty.T_Fem(d), ALL(d'), Base(d'', SOME v1)) => raise Fail "impossible conversion: ALL -> base; could be implemented; could be useful for apply?; if base type, then find it via indexing; if regular type, discard...."
	      | _ => raise Fail ("impossible convert: " ^ (Ty.toString ty) ^ " and " ^ (presToString pres) ^ " and " ^ (presToString targetPres))
	   (* end case*))
       fun toAll(v : V.t, retPres) =
	   if checkMerge(getFemPres v, retPres)
	   then toAll'(cvtVar v, V.typeOf v, getFemPres v,retPres)
	   else ([], cvtVar v)

       fun toAlls(vs : V.t list, target : femPres) =
	   let val (stmss, vs) = ListPair.unzip (List.map (fn x => toAll(x, target)) vs)
	   in (List.concat stmss, vs) end
      in
       (toAll', toAll, toAlls)
      end
	
  (*introduce check pres function -> use in toAll -> introduce in cvtExp to see what is needed.
   use inductive assumption: if a var is called somwhere, it has already conformed to its femPres*)
  fun cvtBody(block : S.block, newDefs : S.func_def list ref, strandParams : femPres list option,
	      globCvt : FD.femType * V.t -> (V.t * S.stmt), (*given a glob, convert it to an intVar and a stm setting up the int*)
	      globDep : V.t -> V.t option, (* get a dependency of a global*)
	      globDepToDep : FD.femType * V.t -> (S.stmt list * V.t), (*given a global,*)
	      depToV : FD.femType  * V.t -> (S.stmt list * V.t)) = (**)
      let
       val S.Block{code, props} = block
       val code' = List.concatMap (fn x => cvtStm(x, newDefs, strandParams, globCvt, globDep, globDepToDep, depToV)) code
      in
       S.Block{code=code', props=PropList.newHolder()}
      end
  and cvtStm(stm : S.stmt, newDefs : S.func_def list ref, strandParams, globCvt, globDep, globDepToDep, depToV) =
      let
       val cvtExp = fn (y,z) => fn x => cvtExp(x, newDefs, strandParams, y, z, globCvt, globDep, globDepToDep, depToV)
       val cvtBody' = fn x => cvtBody(x, newDefs, strandParams, globCvt, globDep, globDepToDep, depToV)
       val lst = fn x => [x]
       fun doit(S.S_Var(v, optE)) = (case Option.map (cvtExp (getFemPres v, v)) optE
				      of SOME((stmts, vr)) => stmts@[S.S_Var(cvtVar v, SOME(vr))]
				       | NONE => [S.S_Var(cvtVar v, NONE)]
				    (* end case*))
	 | doit (S.S_Assign(v, exp)) = let val (stms,exp') = cvtExp (getFemPres v, v) exp in stms@[S.S_Assign(cvtVar v, exp')] end
	 | doit (S.S_Foreach(v1, v2, block)) = [S.S_Foreach(cvtVar v1, cvtVar v2, cvtBody' block)]
	 | doit (S.S_IfThenElse(v, b1, b2)) = [S.S_IfThenElse(cvtVar v, cvtBody' b1, cvtBody' b2)]
	 | doit (S.S_New(a,vs)) =
	   let
	    val (toAll', toAll, toAlls) = makeConversions(globCvt, globDep, globDepToDep, depToV)
	    val stmsXvs = (case strandParams
			of SOME(strandParams') => (ListPair.map toAll (vs, strandParams'))
			 | NONE => raise Fail "impossible"
			   (* end case*))
	    val (stmss, vs') = ListPair.unzip stmsXvs
	    val stms = List.concat stmss
	   in
	    stms@[S.S_New(a, vs')]
	   end
	 | doit (S.S_KillAll) = lst (S.S_KillAll)
	 | doit (S.S_StabilizeAll) = lst(S.S_StabilizeAll)
	 | doit (S.S_Continue) = lst (S.S_Continue)
	 | doit (S.S_Die) = lst (S.S_Die)
	 | doit (S.S_Stabilize) = lst (S.S_Stabilize)
	 | doit (S.S_Return(v)) = lst (S.S_Return(cvtVar v))
	 | doit (S.S_Print(vs)) = lst (S.S_Print(List.map cvtVar vs))
	 | doit (S.S_MapReduce(mrlst)) =
	   let
	    fun translateReduce(S.MapReduce{result, reduction, mapf,args, source, domain}) =
		let
		 val args' = List.map cvtVar args
		 val S.Func{f, params, body} = mapf
		 val defs = ref []
		 val test = Option.isSome(cvtFun(f, args, result, defs, (fn (x,y) => cvtBody(x,y, NONE, globCvt, globDep, globDepToDep, depToV)), cvtVar, true))
		 val mapf' = if test
			     then let val def::vs = !defs (* our def is on top*)
				      val _ = newDefs := vs@ (!newDefs)
				  in
				   def
				  end
			     else mapf
		in
		 S.MapReduce{result=cvtVar result,
			     reduction=reduction,
			     mapf = mapf',
			     source=cvtVar source,
			     args=args',
			     domain=domain}
		end
	   in
	    lst (S.S_MapReduce(List.map translateReduce mrlst))
	   end
      in
       doit(stm)
      end
  and cvtExp(e : S.exp, newDefs : S.func_def list ref, strandParams,
	     retPres : femPres, retVar : V.t,
	     globCvt : FD.femType * V.t -> (V.t * S.stmt),
	     globDep : V.t -> V.t option,
	     globDepToDep : FD.femType * V.t -> (S.stmt list * V.t),
	     depToV : FD.femType  * V.t -> (S.stmt list * V.t)) : S.stmt list * S.exp =
      let
       fun projectCvt v i = ([], S.E_Project(cvtVar v, i))
       val cvtBody = fn (x,y) => cvtBody(x, y, strandParams, globCvt, globDep, globDepToDep, depToV)
       val (toAll', toAll, toAlls) = makeConversions(globCvt, globDep, globDepToDep, depToV)

       fun acquireGlobal(v) = (case getFemPres v
				of Base(_, SOME(v')) => ([],v')
				 | ALL(f) => depToV(f, cvtVar v)
				 | _ => raise Fail "invalid femPres for base fem global"
			      (*end case*))				
						    
       val none = fn s => ([],s)
       fun doit (S.E_Var(v)) = if checkMerge(getFemPres v, retPres)
			       then let val (stms, v') = toAll(v, retPres)
				    in (stms, S.E_Var(v')) end
			       else none (S.E_Var(cvtVar v))
	 | doit (l as S.E_Lit(_)) = none l
	 | doit (k as S.E_Kernel _) = none k
	 | doit (S.E_Select(v1, v2)) =
	   let
	    val actualRetPres = getFemPres v2
	   in
	    if checkMerge(retPres, actualRetPres)
	    then let
	     val oldTy = cvtTy(V.typeOf v2, actualRetPres)
	     val newTy = cvtTy(V.typeOf v2, retPres)
	     val temp = V.new("cvttemp", Var.LocalVar, oldTy)
	     val tempAssign = S.S_Var(temp, SOME(S.E_Select(cvtVar v1, cvtVar v2)))
	     val (cvtStms, temp') = toAll'(temp, oldTy, actualRetPres, retPres) handle ex => raise ex
	    in (tempAssign::cvtStms, S.E_Var(temp')) end
	    else none (S.E_Select(cvtVar v1, cvtVar v2))
	   end
	 | doit (S.E_Apply(f, vs)) = (* possible to all convert*)
	   let
	    val filteredVs = List.filter testParam vs
	    val vs' = List.map cvtVar filteredVs
	    val callSite = vs
	   in
	    (case cvtFun(f, callSite, retVar, newDefs, cvtBody, cvtVar, false)
			(*maybe cvt doesn't make sense here...*)
	      of SOME(f', actualRetPres) => if checkMerge(actualRetPres, retPres)
					    then let
					     val oldTy = cvtTy(F.resultTypeOf f, actualRetPres)
					     val newTy = cvtTy(F.resultTypeOf f, retPres)
					     val temp = V.new("cvtTemp", Var.LocalVar, oldTy)
					     val tempAssign = S.S_Var(temp, SOME(S.E_Apply(f', vs')))
					     val (cvtStms, temp') = toAll'(temp, oldTy, actualRetPres, retPres) handle ex => raise ex
					    in (tempAssign::cvtStms, S.E_Var(temp')) end
					    else ([], S.E_Apply(f', vs'))
	       | NONE => none (S.E_Apply(f, vs'))
	    (* end case*))
	   end
	 | doit (S.E_Prim(f, marg, vs, t)) =
	   (* here we need to check for any Bases that are appended to an all seq or an all being appaned to a base or so forth 
	      - we conver the base to an all via toAll and then we are done. *)
	   (*TODO: double check*)
	   let
	    val vs' = List.map cvtVar vs
	    val t' = cvtTy(t, retPres) handle ex => raise ex
	   in
	    if Var.same(BasisVars.at_dT, f)
	    then let
	     val [syOld, aOld] = vs
	     val [sy, a] = vs'

	     val test1 = checkMerge(getFemPres syOld, retPres)
	     val test2 = checkMerge(Seq(getFemPres aOld), retPres)
	     val Seq(aPesDest) = retPres
	     val (stms1, sy') = if test1
				then toAll(syOld, retPres)
				else ([], sy)
	     val (stms2, a') = if test2
				then toAll(aOld, aPesDest)
			       else ([], aOld)
	     val Ty.T_Sequence(t'', _) = t'
	     val [Ty.TY ty] = marg
	     val marg' = [Ty.TY t'']
	    in
	     (stms1@stms2, S.E_Prim(f, marg', [sy', a], t'))
	    end
	    else if Var.same(BasisVars.at_Td, f)
	    then let
	     val [aOld, syOld] = vs
	     val [a, sy] = vs'

	     val test1 = checkMerge(getFemPres syOld, retPres)
	     val test2 = checkMerge(Seq(getFemPres aOld), retPres)
	     val Seq(aPesDest) = retPres
	     val (stms1, sy') = if test1
				then toAll(syOld, retPres)
				else ([], sy)
	     val (stms2, a') = if test2
				then toAll(aOld, aPesDest)
				else ([], aOld)
	     val Ty.T_Sequence(t'', _) = t'
	     val [Ty.TY ty] = marg
	     val marg' = [Ty.TY t'']			    
	    in
	     (stms1@stms2, S.E_Prim(f, marg', [a', sy'], t'))
	    end
	    else if Var.same(BasisVars.at_dd, f)
	    then let
	     val [sy1Old, sy2Old] = vs
	     val [sy1, sy2] = vs'
	     val test1 = checkMerge(getFemPres sy1Old, retPres)
	     val test2 = checkMerge(getFemPres sy2Old, retPres)
	     val (stms1, sy1') = if test1
				 then toAll(sy1Old, retPres)
				 else ([], sy1)
	     val (stms2, sy2') = if test2
				 then toAll(sy2Old, retPres)
				 else ([], sy2)
	     val Ty.T_Sequence(t'', _) = t'
	     val [Ty.TY ty] = marg
	     val marg' = [Ty.TY t'']
	    in
	     (stms1@stms2, S.E_Prim(f, marg', [sy1', sy2'], t'))
	    end
	    else none (S.E_Prim(f, marg, vs', t'))
	   end
	 | doit (S.E_Field(vs, ty)) = none (S.E_Field(List.map cvtVar vs, ty))
	 | doit (S.E_Tensor(vs,ty)) = none (S.E_Tensor(List.map cvtVar vs, ty))
	 (*Seq though coerce need maybeRunAll too*)
	 | doit (S.E_Seq(vs, ty)) = let
	  val innerPres = (case retPres
			    of Seq(a) => a
			     | Array(a) => a
			     | _ => raise Fail "invalid seq pres" (*end case*))
	  val outerTy = cvtTy(ty, retPres)
	  val anyMerge = checkMerges(List.map getFemPres vs, innerPres)
	 in if anyMerge
	    then let val (stms, newvars) = toAlls(vs, innerPres)
		 in (stms, S.E_Seq(newvars, outerTy)) end
	    else none(S.E_Seq(List.map cvtVar vs, outerTy)) end
	 | doit (S.E_Tuple(vs)) =
	   let
	    val innerPreses = (case retPres
				of Tuple(ts) => ts
				 | _ => raise Fail "impossible" (* end case*))
	    val innerPreses' = List.map getFemPres vs
	    val (stmss, vs') = ListPair.unzip (ListPair.map toAll (vs, innerPreses))
	    val stms = List.concat stmss
	   in (stms, S.E_Tuple(vs')) end
	 | doit (S.E_Project(v,i)) =
	   let
	    val retPres' = (case retPres
			     of Tuple(ts) => List.nth(ts, i)
			      | _ => raise Fail "impossible femPres" (* end case*))
	    val retTyTemp = (case V.typeOf v
			      of Ty.T_Tuple(ts) => List.nth(ts, i)
			       | _ => raise Fail "impossible tuple ty" (* end case*))
	    val retTyTempCvt = (case V.typeOf (cvtVar v)
			      of Ty.T_Tuple(ts) => List.nth(ts, i)
			       | _ => raise Fail "impossible tuple ty" (* end case*))

	    
	   in
	    if checkMerge(retPres', retPres)
	    then let val v' = V.new("tempProject", V.LocalVar, retTyTempCvt)
		     val stm = S.S_Var(v', SOME((S.E_Project(cvtVar v, i))))
		     val (stms, v'') = toAll'(v', retTyTemp, retPres', retPres)
		 in (stm::stms, S.E_Var(v'')) end
	    else none (S.E_Project(cvtVar v, i))
	   end
	 | doit (S.E_Slice(v, s, t)) = none (S.E_Slice(cvtVar v, s, t)) 
	 | doit (S.E_Coerce{srcTy, dstTy, x}) = 
	   let
	    val srcPres = getFemPres x
	   in
	    none (S.E_Coerce {srcTy = cvtTy(srcTy, srcPres),
			      dstTy = cvtTy(dstTy, retPres),
			      x = cvtVar x})
	   end
	 | doit (S.E_BorderCtl(v1, v2)) = none(S.E_BorderCtl(v1, v2))
	 | doit (S.E_LoadSeq(v1, s)) = none(S.E_LoadSeq(v1, s))
	 | doit (S.E_LoadImage(t, s, i)) = none(S.E_LoadImage(t,s,i))
	 | doit (S.E_InsideImage(v1,v2,i)) = none(S.E_InsideImage(cvtVar v1, cvtVar v2, i))
	 | doit (S.E_FieldFn _) = raise Fail "FIXME: field function"
	 | doit (S.E_LoadFem(data, v1, v2)) =
	   if FD.baseFem data
	   then (case (getFemPres v1, retPres)
		  of (Base(d', SOME(v)), Base(d'', SOME(v'))) => none (S.E_Var(v))
		   | (Base(d', SOME(v)), ALL(d'')) => let val (newVar, newStm)  = globCvt(d', v1)
						      in ([newStm], S.E_Var(newVar)) end
		   | (ALL(d'), ALL(d'')) => none (S.E_Var(cvtVar v1))
		   | _ => raise Fail "conversion not allowed for LoadFem"
		(* end case*))
	   else (case (getFemPres v1, retPres, data)
		  of (Base(d', SOME(v)), Base(d'', SOME(v')), FD.MeshCell _) => none (S.E_Var(cvtVar v2))
		   | (Base(d', SOME(v)), Base(d'', SOME(v')), FD.FuncCell _) => none (S.E_Var(cvtVar v2))
		   | (Base(d', SOME(v)), ALL(d''), FD.MeshCell _) =>
		     let val (newVar, newStm)  = globCvt(d', v1)
			 val tupleExp = S.E_Tuple([cvtVar v2, newVar])
		     in ([newStm], tupleExp) end
		   | (Base(d', SOME(v)), ALL(d''), FD.FuncCell _) =>
		     let val (newVar, newStm)  = globCvt(d', v1)
			 val tupleExp = S.E_Tuple([cvtVar v2, newVar])
		     in ([newStm], tupleExp) end
		   | (_, _, FD.MeshPos _) => raise Fail "impossible"
		   | (_, _, FD.RefCell _) => raise Fail "impossible"
		   | (ALL(d'), ALL(d''), _) => none (S.E_Tuple([cvtVar v2, cvtVar v1]))
		(* end case*))
	 | doit (S.E_ExtractFem(v, data)) =
	   let val Ty.T_Fem(data') = V.typeOf v
	   in
	    if FD.baseFem data'
	    then (case (getFemPres v, retPres)
		   of (Base(d, SOME(v')), Base(d', SOME(v''))) =>
		      (case globDep v'
			of SOME v''' => ([], S.E_Var(v'''))
			 | _ => raise Fail ("globDep: " ^ (V.uniqueNameOf v')))
		    | (Base(d, SOME(v')), ALL(d')) =>
		      let val SOME(depV) = globDep v'
			  val (newVar, newStm) = globCvt(data, depV)
		      in ([newStm], S.E_Var(newVar)) end
		    | (ALL(d), ALL(d')) =>
		      let val (stms, newVar) = globDepToDep(d, cvtVar v)
		      in (stms, S.E_Var(newVar)) end
		 (* end case*))
	    else (case (getFemPres v, retPres, data)
		   of (_, _, FD.RefCell _) => raise Fail "impossible"
		    | (Base(d', SOME(v1)), Base(d'', SOME(v1')),  _) => none (S.E_Var(cvtVar v1))
		    | (Base(d', SOME(v1)), ALL(d''),  _) =>
		      let val (newVar, newStm)  = globCvt(d', v1)
		      in ([newStm], S.E_Var(newVar)) end
		    | (ALL(d'), ALL(d''), FD.MeshCell _) =>
		      ([], S.E_Project(cvtVar v, 1))
		    | (ALL(d'), ALL(d''), FD.FuncCell _) =>
		      ([], S.E_Project(cvtVar v, 1))
		    | (ALL(d'), ALL(d''), FD.MeshPos  _) =>
		      ([], S.E_Project(cvtVar v, 4))
		 (* end case*))
	   end
	 | doit (S.E_ExtractFemItem(v, ty, (opt, data))) =
	   if FO.arity opt <> 1
	   then raise Fail "validation error"
	   else if FD.baseFem data
	   then (case getFemPres v
		  of Base(_, SOME(v')) =>
		     none (S.E_ExtractFemItem(v', ty, (opt, data)))
		   | ALL(_) =>
		     let val (stms, newVar) = depToV(data, cvtVar v)
		     in (stms, S.E_ExtractFemItem(newVar, ty, (opt, data)))
		     end
		(* end case*))
	   else
	    (case (opt, getFemPres v, data)
	      of (FO.CellIndex, Base _, FD.MeshCell _) => none (S.E_Var(cvtVar v))
	      |  (FO.CellIndex, Base _, FD.FuncCell _) => none (S.E_Var(cvtVar v))
	      |  (FO.CellIndex, _, FD.MeshPos _) => projectCvt v 1
	      |  (FO.CellIndex, ALL _, FD.MeshCell _) =>
		 projectCvt v 0
	      |  (FO.CellIndex, ALL _, FD.FuncCell _ ) =>
		 projectCvt v 0
	      |  (FO.RefPos, _, FD.MeshPos _) =>
		 projectCvt v 0
	      |  (FO.InvalidBuild, _, _) => raise Fail "use with N"
	      |  (FO.RefCell, _, _) => raise Fail "invalid"
	      |  (FO.PosEntryFacet, _, _) =>
		 projectCvt v 2
	    (* end case*))
	 | doit (S.E_ExtractFemItem2(v1, v2, ty, outTy, (opt, data))) =
	   if FO.arity opt <> 2
	   then raise Fail "validation error"
	   else if FD.baseFem data
	   then (case getFemPres v1
		  of Base(_, SOME(v1')) =>
		     none (S.E_ExtractFemItem2(cvtVar v1', cvtVar v2, ty, outTy, (opt, data)))
		   | ALL(_) =>
		     let val (stms, newVar) = depToV(data, cvtVar v1)
		     in (stms, S.E_ExtractFemItem2(newVar, v2, ty, outTy, (opt, data)))
		     end
		(* end case*))
	   else raise Fail "no extractFemItem2s access non base fem data"
	 (*RefBuild, InvalidBuild, InvalidBuildBoundary, AllBuild, *CellFaceCell*)
	 | doit (S.E_ExtractFemItemN(vs, tys, ty, (opt, data), NONE)) = (*we never use the function here - for now*)
	   if FO.arity opt <> (List.length vs)
	   then raise Fail ("invalid arity " ^ (Int.toString (FO.arity opt)) ^ " with args =" ^ (Int.toString (List.length vs)) ^ " for fem opt: "  ^ (FO.toString (opt, data)))
		      (*function to get the int of the fem involved...*)
	   else
	    (case opt
	      of FO.ExtractDofs =>
		 let
		  val [data1, data2, int1] = vs
		  val acquireGlobals = List.map acquireGlobal [data1, data2]
		  val acquireGlobalsStms = List.concatMap (fn (x,y) => x) acquireGlobals
		  val acquireGlobalsVars = List.map (fn (x,y) => y) acquireGlobals
		  val [Ty.T_Fem(dataty1), Ty.T_Fem(dataTy2), Ty.T_Int] = tys
									   
		 in
		  (acquireGlobalsStms, S.E_ExtractFemItemN(acquireGlobalsVars@[cvtVar int1], tys, ty, (opt, data), NONE))
		 end
	       | _ => 
	    let
	     val femArg :: args = vs
	     val femArg' = cvtVar femArg
	     val args' = List.map cvtVar args
	     val nanVar = V.new("nan", Var.LocalVar, Ty.realTy')
	     val nanDec = S.S_Var(nanVar, SOME(S.E_Lit(Literal.Real(RealLit.nan))))
	     val nanVecVar = V.new("nanVec", Var.LocalVar, Ty.T_Tensor([3], 0))
	     val nanVecDec = S.S_Var(nanVecVar, SOME(S.E_Tensor([nanVar, nanVar, nanVar], Ty.T_Tensor([3], 0))))
	     val neg1 = V.new("neg1", Var.LocalVar, Ty.T_Int)
	     val neg1Dec = S.S_Var(neg1, SOME(S.E_Lit(Literal.intLit (~1))))

	     val presRecover = (case (getFemPres femArg, retPres)
				 of (Base(_, NONE), _) => raise Fail "impossible"
				  | (_, Base(_, NONE)) => raise Fail "impossible"
				  | (Base(_, SOME _), Base(_, SOME _)) => NONE
				  | (Base(_, SOME v), ALL(data')) =>
				    let val (newVar, newStm)  = globCvt(data', cvtVar v)
				    in SOME([newStm], newVar) 
				    end
				  | (ALL(_), ALL(_)) => SOME([], femArg')
				  | (ALL _, Array(NOTHING)) =>
				    let val (stms, newVar) = depToV(data, cvtVar femArg) (*TODO:is this right*)
				    in SOME(stms, newVar)
				    end
				  | (Base(_, SOME v), Array(NOTHING)) => SOME([], cvtVar v)
				  | (a,b) => raise Fail ("impossibe: " ^ (presToString a) ^ " vs. " ^ (presToString b))
			       (* end case*))
	     val startStms = nanDec :: nanVecDec :: neg1Dec :: []

	    in
	     (case opt
	       of FO.RefBuild =>
		  let
		   val [_, cell, refPos, facet] = (case List.map cvtVar vs
						    of [mesh, cell, refTensor] => [mesh, cell, refTensor, neg1]
						     | [mesh, cell, refTensor, facet] => [mesh, cell, refTensor, facet]
						     | _ => raise Fail "impossible"
						  (*end case*))
		  in
		   (case presRecover
		     of NONE =>  (startStms, S.E_Tuple([refPos, cell, facet]))
		      | SOME(stms, dataVar) => (startStms@stms, S.E_Tuple([refPos, cell, facet, dataVar]))
		   (* end case*))
		  end
		| FO.InvalidBuild =>
		  let
		   val [] = List.map cvtVar args
		  in
		   (case presRecover
		     of NONE =>  (startStms, S.E_Tuple([nanVecVar, neg1, neg1]))
		      | SOME(stms, dataVar) => (startStms@stms, S.E_Tuple([nanVecVar, neg1, neg1, dataVar]))
		   (* end case*))
		  end
		| FO.InvalidBuildBoundary =>
		  let
		   val [refPos, facet] = List.map cvtVar args
		  in
		   (case presRecover
		     of NONE => (startStms, S.E_Tuple([refPos, neg1, facet]))
		     | SOME(stms, dataVar) => (startStms@stms, S.E_Tuple([refPos, neg1, facet, dataVar]))
		   (* end case*))
		  end
		| FO.AllBuild =>
		  let
		   val [refPos, cell, facet] = (case List.map cvtVar args
						 of [cellInt, refPos, worldPos, facetIdExp] => [refPos, cellInt, facetIdExp]
						  | [cellIntExp, refPos, worldPos] => [refPos, cellIntExp, neg1]
						  | _ => raise Fail "invalid AllBuild config"
					       (* end case*))
		  in
		   (case presRecover
		     of NONE => (startStms, S.E_Tuple([refPos, cell, facet]))
		      | SOME(stms, dataVar) => (startStms@stms, S.E_Tuple([refPos, cell, facet, dataVar]))
		   (* end case*))
		  end
		| FO.CellFaceCell =>
		  (case presRecover
		    of NONE => raise Fail "impossible: CellFaceCell "
		     | SOME(stms, dataVar) => (startStms@stms, S.E_ExtractFemItemN(dataVar::args', tys, ty, (opt, data), NONE))
		  (* end case*))
		| _ => raise Fail "impossible FO"
	     (* end case*))
	    end
	    (* end case *))
	 | doit (S.E_FemField(v1, v2, v1opt, ty, fo, NONE)) =
	   let
	    val (stms1, v1') = acquireGlobal v1
	    val (stms2, v2') = acquireGlobal v2
	    val (stms3, v1opt') = (case Option.map getFemPres v1opt
				    of NONE => ([], NONE)
				     | SOME(NOTHING) => ([], Option.map cvtVar v1opt)
				     | SOME(_) => let val SOME(base) = v1opt
						      val (stms, base') = acquireGlobal base
						  in (stms, SOME(base')) end
				  )
	   in
	    (stms1@stms2@stms3, S.E_FemField(v1', v2', v1opt', ty, fo, NONE))
	   end
	 | doit (S.E_FemField(v1, v2, SOME(v3), ty, fo, SOME(func))) =
	   let
	    val (stms1, v1') = acquireGlobal v1
	    val (stms2, v2') = acquireGlobal v2
	    val (stms3, v1opt') = (case Option.map getFemPres (SOME(v3))
				    of NONE => ([], NONE)
				     | SOME(NOTHING) => ([], Option.map cvtVar (SOME(v3)))
				     | SOME(_) => let val SOME(base) = SOME(v3)
						      val (stms, base') = acquireGlobal base
						  in (stms, SOME(base')) end
				  )

	    local
	     val (fakeIn, fakeOut) = getInOut retVar
	     val args = (case V.typeOf v3
			  of Ty.T_Int => [fakeIn, v3, v1]
			   | _ => [v3, fakeIn])
	     val result = fakeOut
	     val SOME(f', _) = cvtFun(func, args, result, newDefs, cvtBody, cvtVar, false)
	    in
	    val funcopt = SOME(f')
	    end
	   in
	    (stms1@stms2@stms3, S.E_FemField(v1', v2', v1opt', ty, fo, funcopt))
	   end

	     
	 | doit _ = raise Fail "impossible STM"

	      
      in
       doit(e)
      end
			
  fun translate (tbs as (femTable : (FD.femType, V.t list) HT.hash_table, femDep : SimpleVar.t VTbl.hash_table)) prog =
      let
       val S.Program{props, consts, inputs, constInit, globals,
		     globInit, funcs, strand, create, start, update} = prog

       (*setup new globals to manage fems*)
       val (newGlobals, newInits,
	    globCvt, findDep, globDepToDep, depToV) = buildGlobalFemTables tbs handle ex => raise ex
       val _ = print("building....\n")
       fun doBlock(block : S.block, newDefs : S.func_def list ref, params)
	   = cvtBody(block, newDefs, params, globCvt, findDep, globDepToDep, depToV) handle ex => raise ex
       val newDefs = ref []
       fun doBlock'(b) = doBlock(b, newDefs, NONE) handle ex => raise ex
       val globals' = List.map cvtVar globals
       val globals'' = newGlobals@globals'
       local
	val S.Block({code,props}) = doBlock'(globInit) handle ex => raise ex
       in
       val globalInit' = S.Block{code=newInits@code,props=PropList.newHolder()}
       end
       val constInit' = doBlock'(constInit) handle ex => raise ex

       val start' = Option.map doBlock' start handle ex => raise ex
       val update' = Option.map doBlock' update handle ex => raise ex

       (*We get:
	 newGlobals: list of global arrays of fems and ints - these allow us to translate globals to ints and get dependencies and fems from inits; we have essentialy maps from fem -> int and int -> int
	 newInits: initiatlization of these globals
	 function: globCvt: that takes a known global and produces the int for it - essetnialy inverts the fem-> int array
	 function: globDep: takes a fem var and produces the dep var -> goes from fem -> int to int-> int and int -> fem without hassel
	 function: globDepToDep: takes a fem and an int representing it to an int representing the dep
	 function: depToV: takes an int representing a FEM and translates it from an int into a the actual FEM
	 Note: the last two basically just build access to an array.
	*)
       local
	val S.Strand{name, params, spatialDim, state, stateInit, startM, updateM, stabilizeM} = strand

	fun doBlock''(b) = doBlock(b, newDefs, SOME(List.map getFemPres params)) handle ex => raise ex
	val codeStms : S.block = Create.createCode create 
	val codeStms' = doBlock'' codeStms handle ex => raise ex

	val params' = List.map cvtVar params
	val state' = List.map cvtVar state
	val stateInit' = doBlock'' stateInit handle ex => raise ex
	val updateM' = doBlock'' updateM handle ex => raise ex
	val startM' = Option.map doBlock'' startM handle ex => raise ex
	val stabilizeM' = Option.map doBlock'' stabilizeM handle ex => raise ex

       in
       val create' = Create.Create {dim = Create.arrayDim create, code=codeStms'}
       val strand' = S.Strand{name=name, spatialDim=spatialDim, params=params', state=state',
			      stateInit=stateInit', startM=startM',
			      updateM=updateM', stabilizeM=stabilizeM'}
       end

       local
	val allDefs = funcs@ (!newDefs)
	val _ = List.app procFuncDef allDefs
	val allDefs' = sortFuncDef allDefs
	val funcsFilt = List.filter (fn f => F.useCount (defToName f) > 0) allDefs'
	fun cvtVar'(v) = (case V.kindOf v
			   of Var.GlobalVar => cvtVar v
			    | Var.InputVar => cvtVar v
			    | Var.ConstVar => cvtVar v
			    | _ => v)
	val funcs' = List.map (fn (func as S.Func{f, params, body}) => if cleanBody(f)
								     then S.Func{f=f, params=params, body=S.mapBlock(body, cvtVar')}
								     else func) funcsFilt
			       
       in
       val funcs' = funcs'
       end
      in
       S.Program{props=props, consts=List.map cvtVar consts, inputs=inputs, constInit=constInit', globals=globals'',
		 globInit=globalInit',
		 funcs=funcs', strand=strand', create=create', start=start', update=update'}
      end

  fun analysisOnly prog =
      let
       val (prog', femTable, femDep) = analysis prog handle ex => raise ex
									(*modify program for new defs...*)
      in
       prog'
      end

	
  fun transform prog =
      let
       val (prog', femTable, femDep) = analysis prog handle ex => raise ex
       val prog'' = translate (femTable, femDep) prog' handle ex => raise ex
       val _ = print("woopdydo\n\n")
      in
       prog''
      end

  end
