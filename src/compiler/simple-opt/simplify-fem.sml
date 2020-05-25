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
       val env = V.Tbl.mkTable(10 * List.length params, raise Fail "copy table fail")
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
					   in (remove x; S.S_Foreach(x', xs', blk')) end
	      | S.S_New _ => raise Fail "new impossible in function!"
	      | S.S_KillAll => raise Fail "kill impossible in function!"
	      | S.S_StabilizeAll => raise Fail "stab impossible in function!"
	      | S.S_Continue  => raise Fail "continue impossible in function!"
	      | S.S_Die => raise Fail "die impossible in function"
	      | S.S_Return x => S.S_Return(rename x)
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
	      (*^NOTE: no globals are in these functions*)
	      | S.E_ExtractFemItemN(vars, tys, outTy, opt, NONE) => S.E_ExtractFemItemN(List.map rename vars, tys, outTy, opt, NONE) 
	   (* end case *))
      in
       S.Func{f=f, params=params',body=doBlock(body)}
      end


  (*Data structure and helper functions for symbolically reasoning about fem data*)
  datatype femPres = Base of FD.femType * V.t option | ALL of FD.femType | NOTHING (* all, nothing, or one thing *)
		     | Tuple of femPres list | Array of femPres 
		     | Seq of femPres (* must be homogenous i.e no partitioning on an array*)

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

  (*Propery for mapping vars to their symbolic fem info; inputs handled specially*)
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

  val getFemPres = Option.valOf o peekFn
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
  type functionCall  = femPres list * femPres list (* arguments and globals *)
  type callResult = S.func_def * femPres option
  local
   fun sameSite(vs1, vs2) = ListPair.all (V.same) (vs1, vs2)
   fun hashSite(vs) = List.foldl (fn (v, s) => V.hash v + s) 0w3 vs
   fun sameCall((args1, globs1), (args2, globs2)) = ListPair.all same (args1, args2)
						    andalso ListPair.all same (globs1, globs2)
   fun hashCall((args1, globs1)) = (List.foldl (fn (v, s) => 0w5 * hash v + s) 0w3 args1) +
				   (List.foldl (fn (v, s) => 0w7 * hash v + s) 0w3 globs1)

  fun mkCallSiteTable(j : int) : (callSite, functionCall) HT.hash_table =
      HT.mkTable (hashSite, sameSite) (j, raise Fail "missing call site")
  fun mkFuncDefTable(j : int) : (functionCall, callResult) HT.hash_table =
      HT.mkTable (hashCall, sameCall) (j, raise Fail "missing call site")

   val {clrFn, getFn, peekFn, setFn} = F.newProp(fn f => mkCallSiteTable 8)
   val defTblProp = F.newProp(fn f => mkFuncDefTable 8)
   val getFnTbl = #getFn defTblProp
												
  in
  (*operation:we need to from a f, 
    and info about the calls (args +glob femPress) -> SOME(function_def)	
    (SOME if we need to do something and NONE if it is already done -> how to get ret var then?)  -> function_def can be processed *)
  (*from the procced thing, we get the accepted one, which we get back here.*)
  fun addCall(f : F.t, call : callSite, callArgs : functionCall) : callResult =
      let
       val callTable = getFn f
       val defsTable = getFnTbl f
       val callLookup : functionCall option = HT.find callTable call
       val callExists = (case callLookup
		of SOME(callArgs') => sameCall(callArgs', callArgs)
		 | NONE => (HT.insert callTable (call, callArgs); false)
			(* end case*))
			  
      in
       if callExists
       then HT.lookup defsTable callArgs
       else let
	val _ = HT.insert callTable (call, callArgs)
	val SOME(oldDef) = getfuncDef f
		val copyF = functionCopy(oldDef)
		val ret = (copyF, NONE)
			    (*check if the function with these callArgs exists anyway*)
       in (case HT.find defsTable callArgs
	    of SOME(r) => r
	     | NONE => (HT.insert defsTable (callArgs, ret); ret)
	  (* end case *))
       end
      end

  fun updateCall(f : F.t, call : callSite, callArgs : functionCall,
		 newDef : S.func_def, retF : femPres) : unit =
      let
       val defsTable = getFnTbl f
      in
       HT.insert defsTable (callArgs, (newDef, SOME(retF)))
      end
  end
  
  (*function to actually analyze a block of code*)
  fun analyzeBlock(b as S.Block{code, ...},
		   changeOuter : bool ref,
		   newsOuter : (Atom.atom * V.t list) list ref,
		   getDep : V.t -> V.t option,
		   addDep : V.t * V.t -> unit) =
      let
       fun expToFemPres(exp : S.exp, retTy : Ty.ty, call : callSite) : femPres  =
	   let
	    val change = ref false (* changes in here don't matter?*)
	    val merge = (fn (x,y) => merge(x,y,change))
	    (*These three are the only place where retTy is used:*)
	    val retHasFem = Ty.hasFem retTy
	    val nonFemRet = tyToFemPres retTy
	    val possibleFemData = (case retTy
				    of Ty.T_Fem(f') => SOME(f')
				     | _ => NONE)
	    val mergeFemPreses = List.foldr (fn (x,y) => merge(x, y))
	    fun doit (S.E_Var(v)) = (checkExistence ("_var_" ^ (V.nameOf v)) v; getFemPres v)
	      | doit (S.E_Lit(l)) = NOTHING
	      | doit (S.E_Kernel(k)) = NOTHING
	      | doit (S.E_Select(v1, v2)) = (
	       checkExistence ("_strand_" ^ (V.nameOf v1)) v1;
	       checkExistence ("_strand_var_" ^ (V.nameOf v2)) v2;
	       getFemPres v2 (* read strand state *)
	      )
	      | doit (S.E_Apply(f, vs)) = (
	       List.app (checkExistence ((F.nameOf f) ^ "_arg")) vs;
	       procApply(f, vs, call); (*these are the only places where call is used.*)
	       getFemPres (List.hd call)
	      )
	      | doit (S.E_Prim(v, [], vs, t)) =
		(List.app (checkExistence ("_arg")) vs;
		 if Var.same(v, BasisVars.subscript) andalso retHasFem
		 then let val [sy, index] = vs (*this merge might be pointless*)
			  val Array(p') = getFemPres sy
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
			  val p'' = merge(p, p')
		      in Seq(p'') end
		 else if Var.same(v, BasisVars.at_Td)
		 then let val [a, sy] = vs
			  val Seq(p) = getFemPres sy
			  val p' = getFemPres a
			  val p'' = merge(p, p')
		      in Seq(p'') end
		 else if Var.same(v, BasisVars.at_dd)
		 then let val [sy1, sy2] = vs
			  val Seq(p1) = getFemPres sy1
			  val Seq(p2) = getFemPres sy2
			  val p = merge(p1, p2)
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
						 then Seq(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs))
						 else Seq(tyToFemPresSeq t') (*{} -> ALL*)
		    | Ty.T_Sequence(t', SOME(k)) => if k <> 0
						    then Array(mergeFemPreses (tyToFemPres t') (List.map getFemPres vs))
						    else Array(tyToFemPresSeq t') (* {} -> ALL*)
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
			       then tyToFemPresSeq ty  (*QUESTION: this is a bit weird as we have the p' - for {}, p' would already be all so this is redudant*)
			       else merge((tyToFemPres ty), p') (*if there is a FEM, there is an init so Base(x, SOME k) gets used*)
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
						       else nonFemRet (*should not have Base/ALL*))
	      (*above could be:c cells, refCell -> but these should be eliminated in prior phases*)
	      | doit (S.E_ExtractFemItem2(v1, v2, t1, t2, fo)) = (checkExistence "_EFI2_1" v1;
								  checkExistence "_EFI2_2" v2;
								  if retHasFem
								  then raise Fail "impossible"
								  else nonFemRet)
	      | doit (S.E_FemField(v1, v2, v3o, t, fof, func)) = (checkExistence "_FField_1" v1;
								  checkExistence "_FField_2" v2;
								  Option.app (checkExistence "_FField_3") v3o;
								  (*NOTE: function here is FEM free modulo a mesh arg!*)
								  nonFemRet)
	      | doit (S.E_ExtractFemItemN(vs, tys, t, fo, NONE)) = ((List.app (checkExistence "_FFIN") vs);
								    if Bool.not retHasFem
								    then nonFemRet
								    else
								     (*refBuild, InvalidBuild, InvalidbuildBoundary, AllBuild  *)
								     let
								      (*NOTE: function here is currently impossible.*)
								      val (opt, data) = fo
								      val n = List.length tys
								      val tabed = List.tabulate(n, fn x => (x, List.nth(tys, x)))
								      val pos  = List.filter
										   (fn (_, Ty.T_Fem(f')) => FD.same(data, f')
										   | _ => false) tabed
								     in
								      (case pos
									of [(idx, _)] => getFemPres (List.nth(vs, idx)) (*propogate Base/ALL*)
									 | _ => raise Fail "not planned in extractFemItemN")
								     end)
								     
	      | doit (S.E_InsideImage(v1, v2, _)) = (List.app (checkExistence ("_in")) [v1, v2]; NOTHING)
	      | doit (S.E_FieldFn _) = NOTHING (* FIXME:handle function *)
	   in
	    doit(exp)
	   end
       and procStatement(stm : S.stmt, changed : bool ref,
			 call : callSite,
			 news : (Atom.atom * V.t list) list ref) =
	   let
	    fun doit (S.S_Var(v, NONE)) = () (*Question: should we do anythin?*)
	      | doit (S.S_Var(v, SOME(e))) = updateFemPresRef(v, expToFemPres(e, V.typeOf v, v::call), changed)
	      (*QUESTION: should we check existence @ an assign?*)
	      | doit (S.S_Assign(v, e)) = updateFemPresRef(v, expToFemPres(e, V.typeOf v, v::call), changed) 
	      | doit (S.S_IfThenElse(v, b1, b2)) = (checkExistence "condition" v; doBlock(b1); doBlock(b2))
	      | doit (S.S_New(a,vs)) = news := (a, vs) :: !news
	      | doit (S.S_Foreach(itter, src, blk)) = (checkExistence "itter" itter;
						       checkExistence "itterSrc" src;
						       doBlock(blk))
	      | doit (S.S_KillAll) = ()
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
		    val S.Func{f, ...} = mapf
		    val _ = registerFunctions([mapf])
		    val _ = registerGlobVarsInFunc mapf
		    val _ = List.app (checkExistence "_map_args_") args
		    (*QUESTION: val _ = (checkExistence "_strand") source: we generaly ingore the strand type *)
		    (*NOTE: args are just args to the function*)
		    (*NOTE: reduction and domain can be ignored *)
		    (*NOTE: This is a case of an apply where args are the params and result is the callsite *)
		   in
		    procApply(f, args, [result])
		   end
	      in
	       List.app doReduce maps
	      end
	    and doBlock(S.Block{code, ...}) = List.app doit code
	   in
	    doit(stm)
	   end
       and doBlock(S.Block{code, ...}, changed : bool ref,
		   ret : callSite,
		   news : (Atom.atom * V.t list) list ref) =
	   let
	    val g =  (fn s => procStatement(s, changed, ret, news))
	   in
	    List.app g code
	   end
       and procApply(f : F.t, args : V.t list,
		     call : callSite) = let
	val argsp = List.map getFemPres args
	val globs = getFuncGlobs f
	val globsp = List.map getFemPres globs
	val callArgs : functionCall = (argsp, globsp)
	val femsInvolved = List.concatMap (Ty.allFems o V.typeOf) (args@globs)
	val anyNonBase = List.exists (Bool.not o FD.baseFem) femsInvolved
	val v : V.t = List.hd call
       in
	if anyNonBase
	then
	 let
	  val (fdef, possibleRet) = addCall(f, call, callArgs) 
	 in
	(case possibleRet
	  of SOME(_) => () 
	   | NONE =>
	     let
	      (* in this path, we can't find an already inlined version of this function so there was a change.*)
	      val _ = changeOuter := true 
	      val S.Func{params, body, ...} = fdef
	      val _ = ListPair.map updateFemPres (params, argsp)
	      val _ = doBlock(body, ref false, call, ref []) (*we can ignore changes/new here as we found them/they are impossible.*)

	      val femPres = getFemPres v
	     in
	      updateCall(f, call, callArgs, fdef, femPres) : unit
	     end
	(* end case*))
	 end
	else updateFemPresRef(v, tyToFemPresSeq (V.typeOf v), changeOuter) (*no fem involved, but we mark this.*)
       end
      in
       doBlock(b, changeOuter, [], newsOuter)
      end
										      

  fun transform prog = let
   val S.Program{props, consts, inputs, constInit, globals,
		 globInit, funcs, strand, create, start, update} = prog

   val _ = registerFunctions(funcs)
   val _ = List.app registerGlobVarsInFunc funcs
   val S.Strand{name, params, spatialDim,  state, stateInit, startM, updateM, stabilizeM} = strand
   (*sml docs said: explain the use and semantics of mkTable HERE.*)
   val femTable : (FD.femType, V.t list) HT.hash_table = HT.mkTable(FD.hash, FD.same) (32, Fail "femtype not added yet")
   fun addType (ty, v) = (case HT.find femTable ty
			   of SOME(vs) => HT.insert femTable (ty, v::vs)
			    | NONE => HT.insert femTable (ty, [v])
			 (*end case*))
   val femDep = VTbl.mkTable(32, Fail "fem dep not found")
   fun addDep (v1 : V.t, v2 : V.t) = VTbl.insert femDep (v1, v2)
   fun findDep ( v1 : V.t) : V.t option = VTbl.find femDep v1

   (*Figure out the base fem inputs.*)
   val _ = procBaseInputs(inputs, addType)

   (* register const inputs - no fem*)
   val _ = List.map getFemPres consts


   (*

NOTE: peekFn doesn't create the thing -> getFn or setFn will though ->all vars need to be accounted for then.
NOTE: need to prevent analysis of non-fem functions.
NOTE: updateCall/those tables won't be around (peekFn => NONE) for functions with no fem/no new defs

Global question: what about globalInit/things that are not inited?
     What about function use later on?
     GetFn/setFn question
     globInit -> create -> start -> doStrands -> update -> loop from doStrands as above
    setupFn()
    strands()
    globs()
    loop (strands(); globs())*)
   val _ = ()

  in
   prog
  end
  end
