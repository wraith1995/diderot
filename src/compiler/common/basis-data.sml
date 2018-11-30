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
			 
	  end = struct

datatype t = BasisFunc of {
	  def : string,
	  mono : real list
	 }
val empty = BasisFunc {def = "", mono = []}

end
