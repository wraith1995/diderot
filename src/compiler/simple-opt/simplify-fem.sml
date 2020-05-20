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



  (*Make data structures for FEM possiblie assocations: list, seq, Tuple, ... - yay for SSA in conversions*)

  datatype femPres = Base of FD.femType * V.t option | ALL of FD.femType | NOTHING (* all, nothing, or one thing *)
		     | Tuple of femPres list | Array of femPres list (* variety *)
		     | Seq of femPres (* must be homogenous i.e no partitioning on an array*)


  fun merge(a1 : femPres, a2 : femPres) =
      let
       val change : bool ref  = ref false
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
		    | (Base(_, NONE), Base(_, NONE)) => raise Fail "can't merge unit femPres."
		    | (NOTHING, NOTHING) => NOTHING
		    | (Base _, ALL(f)) => changed(ALL(f))
		    | (ALL(f), Base _) => changed(ALL(f))
		    (*end possible base cases*)
		    | (Tuple(tys1), Tuple(tys2)) => Tuple(ListPair.map merge' (tys1, tys2))
		    | (Array(tys1), Array(tys2)) => Array(ListPair.map merge' (tys1, tys2))
		    | (Seq(t1), Seq(t2)) => Seq(merge'(t1, t2))
		    | _ => raise Fail "impossible: merging incompatible femPress - must of come from different types"
		 (*end case*))
      in
       (merge'(a1, a2), !change)
      end


  fun tyToFemPres (Ty.T_Bool) = NOTHING
    | tyToFemPres (Ty.T_Int) = NOTHING
    | tyToFemPres (Ty.T_String) = NOTHING
    | tyToFemPres (Ty.T_Tensor _) = NOTHING
    | tyToFemPres (Ty.T_Sequence(t, SOME k)) = Array(List.tabulate(k, fn x => tyToFemPres t))
    | tyToFemPres (Ty.T_Sequence(t, NONE)) = Seq(tyToFemPres t)
    | tyToFemPres (Ty.T_Tuple(tys)) = Tuple(List.map tyToFemPres tys)
    | tyToFemPres (Ty.T_Strand _) = NOTHING
    | tyToFemPres (Ty.T_Field _) = NOTHING
    | tyToFemPres (Ty.T_Fem(f)) = if FD.baseFem f
				then Base(f, NONE)
				else Base(Option.valOf(FD.dependencyOf(f)), NONE)
	
  local
   val {clrFn, getFn, peekFn, setFn} = V.newProp(fn v => tyToFemPres(V.typeOf(v)))
  in
  fun updateFemPres(v : V.t, p : femPres) = (case peekFn(v)
					      of SOME(p') => (let val (p'', change) = merge(p, p')
							      in setFn(v, p''); change end)
					       | NONE => (setFn(v, p); true)
					    (*end case*))
  fun setInput(v : V.t) =
      (case V.kindOf v
	of Var.InputVar => (case tyToFemPres(V.typeOf(v))
			     of Base(f, NONE) => setFn(v, Base(f, SOME(v)))
			      | alt => setFn(v, alt) (*only Base(f, NONE) can appear in inputs*)
			   (*end case*))
	 | _ => raise Fail "invalid use of setInput"
      (*end case*))
  end

  type callSite =  V.t list  (* structure for idying a call site *)
  (* think about functions...*)



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



	

  fun transform prog = let
   val S.Program{props, consts, inputs, constInit, globals,
		 globInit, funcs, strand, create, start, update} = prog

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
