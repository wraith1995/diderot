(* api-types.sml
 *
 * A representation of the types of values that can be communicated to and from a
 * Diderot program.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *)

structure APITypes =
  struct

    datatype t
      = IntTy
      | BoolTy
      | TensorTy of int list
      | StringTy
      | ImageTy of int * int list
      | FemData of FemData.femType
      | SeqTy of t * int option
      | TupleTy of t list

    val realTy = TensorTy[]

    fun hasFem (SeqTy(ty,_)) = hasFem ty
      | hasFem (FemData(_)) = true
      | hasFem (TupleTy(tys)) = not (List.all (not o hasFem) tys)
      | hasFem _ = false
    fun depth ty =
	      let
	       fun depth'(SeqTy(ty, SOME(s)), ds) = depth'(ty, s::ds)
		 | depth'(_, ds) = List.rev ds
	      in
	       depth'(ty, [])
	      end
    (*NOTE: IN A JUST WORLD, THESE WOULD NOT BE HERE.*)
    fun toSOA (SeqTy(s, NONE)) =
	let
	 fun toSOA' ty =
	     (case ty
	       of SeqTy(ty', NONE) => raise Fail "[] [] detected"
		| SeqTy(ty', SOME(n)) => SeqTy(toSOA' ty', SOME(n))
		| TupleTy(tys) => TupleTy(List.map toSOA' tys)
		| _ => SeqTy(ty, NONE)
	     (*End Case*))
	in
	 toSOA' s
	end
      | toSOA _ = raise Fail "invalid input to SOA"
			
    fun toAOS ty =
	let
	 fun cleanTy (TupleTy(tys)) = TupleTy(List.map cleanTy tys)
	   | cleanTy (SeqTy(ty', NONE)) = ty'
	   | cleanTy (SeqTy(ty', SOME(n))) = SeqTy(cleanTy ty', SOME(n)) 
	in
	 SeqTy(cleanTy ty, NONE)
	end		

    fun toString IntTy = "int"
      | toString BoolTy = "bool"
      | toString (TensorTy[]) = "real"
      | toString (TensorTy[2]) = "vec2"
      | toString (TensorTy[3]) = "vec3"
      | toString (TensorTy[4]) = "vec4"
      | toString (TensorTy dd) = concat["tensor[", String.concatWithMap "," Int.toString dd, "]"]
      | toString (TupleTy(tys)) = concat ["(", concat (List.map toString tys), ")"]
      | toString StringTy = "string"
      | toString (ImageTy(d, dd)) = concat[
            "image(", Int.toString d, ")[", String.concatWithMap "," Int.toString dd, "]"
          ]
      | toString (SeqTy(ty, NONE)) = toString ty ^ "[]"
      | toString (SeqTy(ty, SOME d)) = concat[toString ty, "[", Int.toString d, "]"]
      | toString (FemData data) = "FemData:" ^ (FemData.toString data)

  (* does a type have a non-static size? *)
    fun hasDynamicSize StringTy = true
      | hasDynamicSize (ImageTy _) = true
      | hasDynamicSize (FemData _) = true (*TODO: Clarify the use of this function for this entry.*)
      | hasDynamicSize (SeqTy(_, NONE)) = true
      | hasDynamicSize _ = false

  end
