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

structure S = Simple
structure TY = SimpleTypes
structure V = SimpleVar
structure VSet = SimpleVar.Set
structure VMap = SimpleVar.Map
structure FT = FemData
structure FO = FemOpt

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



		 


    fun transform prog = prog
end
