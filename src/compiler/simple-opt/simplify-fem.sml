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

  (*Proc globals and global inits*)
  (*Make data structures for FEM possiblie assocations: list, seq, Tuple, ... - yay for SSA in conversions*)
  (*make prop for assinging these*)

  datatype femPres = Base of FD.femType * V.t option | ALL of FD.femType | NOTHING (* all, nothing, or one thing *)
		     | Tuple of femPres list | Array of femPres list (* variety *)
		     | Seq of femPres (* must be homogenous i.e no partitioning on an array*)

  fun toFemPres (Ty.T_Bool) = NOTHING
    | toFemPres (Ty.T_Int) = NOTHING
    | toFemPres (Ty.T_String) = NOTHING
    | toFemPres (Ty.T_Tensor _) = NOTHING
    | toFemPres (Ty.T_Sequence(t, SOME k)) = Array(List.tabulate(k, fn x => toFemPres t))
    | toFemPres (Ty.T_Sequence(t, NONE)) = Seq(toFemPres t)
    | toFemPres (Ty.T_Tuple(tys)) = Tuple(List.map toFemPres tys)
    | toFemPres (Ty.T_Strand _) = NOTHING
    | toFemPres (Ty.T_Field _) = NOTHING
    | toFemPres (Ty.T_Fem(f)) = if FD.baseFem f
				then Base(f, NONE)
				else Base(Option.valOf(FD.dependencyOf(f)), NONE)
	

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

  in
   prog
  end
  end
