(* basis-data.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * Information about a basis for a function space.
*)


structure BasisData : sig
	   type t

	   val empty : t
	   val same : t * t -> bool

	   val hash : t -> word
			 
	  end = struct

datatype t = BasisFunc of {
	  def : Atom.atom,
	  mono : real list (*TODO: Figure out polynomial conventions*)
	 }
val empty = BasisFunc {def = Atom.atom "", mono = []}

fun polySame (p :: ps) (q :: qs) = Real.==(p,q) andalso (polySame ps qs)
  | polySame [] [] = true
  | polySame (p :: ps) [] = false
  | polySame [] (q :: qs) = false


fun getCoeffs(BasisFunc{def , mono} : t) : real list = mono

fun same(b1,b2) = polySame (getCoeffs b1) (getCoeffs b2)
fun hash( BasisFunc{def, mono}) = Atom.hash def


    

end
