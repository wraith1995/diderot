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

    fun isOutputAble(ty) =
	(case ty
	  of IntTy => true
	   | BoolTy => true
	   | TensorTy _ => true
	   | StringTy => true
	   | ImageTy(_, _) => raise Fail "ImageTy should not be considered in output"
	   | FemData _ => raise Fail "FemData should not be considered in output"
	   | SeqTy(ty', SOME(n)) => isOutputAble ty'
	   | SeqTy(ty', NONE) => isOutputAble ty'
	   | TupleTy _ => false 
	(*end case*))

    fun buildAccessPattern(bTy) =
	(*Given a type, figure out a way to access all elements of the form [n][n][n](real or tensor or int or bool or stringTy)*)
	(*Returns a [(path, ty)] where the nth element contains a path to that element through the tree and the result type; ~1 indicates going through a sequence*)
	
	let
	(*basically a sort of pre-order traversal*)
	 fun bap(ty, path) =
	      if isOutputAble ty
	      then [(List.rev path, ty)]
	      else
	       (case ty
		 of SeqTy(ty', SOME _) => bap(ty', ~1 :: path)
		  | SeqTy(ty', NONE) => bap(ty', ~2 :: path)
		  | TupleTy(tys) => List.concatMap bap (List.tabulate(List.length tys, fn x => (List.nth(tys, x), x :: path)))
		  | _ => raise Fail "Impossible"
	       (*end case*))
	in
	 bap(bTy, [])
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

    fun toOutputAbleType(ty) =
	(*Converts types to the form Tuple(base[]) so that base = base Ty ([n][n]...) ad so that base is in pre-order traversal order *)
	let
	 fun isTuple (n, TupleTy _) = true
	   | isTuple _ = false
	 (*inefficient but simple:
	  *First, flatten a nested tupple i.e Tuple(..., Tuple, ...)
	  *Second, traverse tuple to analyze next level
	  *If a seq type is ever found with a tuple directly inside it, switch them.
	  *Otherwise, continue down
	  *preverse other types, except fem types wher we fail
	 **By applying thes rules, we put the type in the suitable form and we maintain pre-order rules on base types as nodes (i.e if you wrote out the pre-order, took out all non-base types, the orders between the original type and toOutpuTableType would be the same.)
	 *)
	 fun toat(ty) =
	      (case ty
		of TupleTy(tys) => (case List.find isTuple (List.tabulate(List.length tys, fn x => (x, List.nth(tys, x))))
				     of SOME((n, TupleTy(tys'))) => let
				     in
				      TupleTy(List.take(tys, n) @ tys' @ List.drop(tys, n + 1))
				     end
				      | NONE => TupleTy(List.map toat tys)
				   (*end case*))
		 | SeqTy(ty', r) => (case ty'
				     of TupleTy(tys) => TupleTy(List.map (fn t => SeqTy(t, r)) tys)
				      | _ => SeqTy(toat ty', r)
				   (*end case*))
		| FemData _ => raise Fail "unexpected Fem Data"
		| _ => ty 
	      (*end case*))
	 fun loop ty = let
	  val ty' = toat ty
	 in
	  if (toString ty) = (toString ty') (*TODO: Fix stupid lazy hack -> will fail on FEM data, but this will crash earlier when making access pattern*)
	  then ty'
	  else loop ty'
	 end
	in
	 loop ty
	end



  (* does a type have a non-static size? *)
    fun hasDynamicSize StringTy = true
      | hasDynamicSize (ImageTy _) = true
      | hasDynamicSize (FemData _) = false (*TODO: Clarify the use of this function for this entry.*)
      | hasDynamicSize (SeqTy(_, NONE)) = true
      | hasDynamicSize (TupleTy(tys)) = (List.exists hasDynamicSize tys)
      | hasDynamicSize _ = false

  end
