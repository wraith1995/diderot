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
val tds = "transformDofShape" (*computed*)
val transform = "transform" (*field*)
val transform' = "T" (*field*)
val invTransform = "inverseTransform"
val invTransform' = "invT"
val meshPos = "findPos"
val meshPosHidden = "$findPos"
val hiddenNewtonInverse = "newtonInverse"
val cellEnter = "enter"
val cellEnterPos = "enterPos"
val invalidPos = "invalidPos"
val invalidCell = "invalidCell"
(*ref cell:*)
val refCellIsInside = "isInside"
val isInsideMesh = "isInside"
val isInsideMeshCell = "isInside"
val isValidCell = "isValid"
val refVerts = "vertices"
val refExit = "exit"
val refEnter = "enter"
val refExitPos = "exitPos"
val refMeshPos = "meshPos"
(*space: *)
val spaceMapDim = "spaceMapDim"
val dim = "dim"
val rangeShape = "shape"
val basisFunctionShape = "basisFunctionShape" (*not used*)
val sds = "spaceDofShape" (*computed*)

(*func:*)
val rangeShape = "shape"
val fds = "funcDofShape" (*not really used now as this is under current spec same as above*)
val funcDofs = "dofs"
val refField = "refField"
val refField' = "F"
val trf = "transformedRefField";
val trf' = "TRF";
val finvt = "FinvT";
val field = "F";
val funcCell = "funcCell"
val space = "space"

(*meshPos:*)
val cell = "mc"
val refPos = "refPos"
val worldPos = "worldPos"
val worldPosHidden = "$worldPos"
val valid = "isValid"

val posExit = "exit"
val posExitPos = "exitPos"

val pos = "pos"
end
