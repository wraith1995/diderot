(* fem-name.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * The central control for various names of fem stuff.
 *)

structure FemName  = struct
(*mesh:*)

val dim = "dim"
val maxDegree = "degree"
val meshMapDim = "meshMapDim"
val tds = "transformDofShape"
val transform = "transform"
(*space: *)
val spaceMapDim = "spaceMapDim"
val dim = "dim"
val basisFunctionShape = "basisFunctionShape" (*not used*)
val rangeShape = "shape"
val sds = "spaceDofShape"

(*func:*)
val rangeShape = "shape"
val fds = "funcDofShape" (*not really used now as this is under current spec same as above*)
	    

end
