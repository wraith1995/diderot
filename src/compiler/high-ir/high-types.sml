(* high-types.sml
 *
 * Types for the HighIR.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure HighTypes =
  struct

    datatype ty
      = BoolTy | StringTy | IntTy
      | TensorTy of int list * int option            (* tensor types, which include reals, vectors, etc. *)
      | TupleTy of ty list              (* tuples; used for multiple return values *)
      | SeqTy of ty * int option
      | ImageTy of ImageInfo.t
      | FemData of FemData.femType
      | StrandTy of Atom.atom
      | KernelTy
      | FieldTy

    val intTy = IntTy
    val realTy' = TensorTy([], NONE)
    val realTy = TensorTy([], NONE)
    fun vecTy' 1 = realTy'
      | vecTy' n = TensorTy([n],NONE)
    val vec2Ty' = vecTy' 2
    val vec3Ty' = vecTy' 3
    val vec4Ty' = vecTy' 4

  (* smart constructor for tensor type that prunes out dimensions with size 1 *)
    fun tensorTy' dd = TensorTy(List.mapPartial (fn 1 => NONE | d => SOME d) dd, NONE)
    fun tensorTy dd k = TensorTy(List.mapPartial (fn 1 => NONE | d => SOME d) dd, k)

    fun same (BoolTy, BoolTy) = true
      | same (StringTy, StringTy) = true
      | same (IntTy, IntTy) = true
      | same (TensorTy (dd1, NONE), TensorTy (dd2, NONE)) = (dd1 = dd2)
      | same (TensorTy (dd1, SOME k1), TensorTy (dd2, SOME k2)) = (dd1 = dd2) andalso k1 = k2
      | same (TupleTy tys1, TupleTy tys2) = ListPair.allEq same (tys1, tys2)
      | same (SeqTy(ty1, NONE), SeqTy(ty2, NONE)) = same(ty1, ty2)
      | same (SeqTy(ty1, SOME n1), SeqTy(ty2, SOME n2)) = (n1 = n2) andalso same(ty1, ty2)
      | same (ImageTy info1, ImageTy info2) = ImageInfo.sameShape(info1, info2)
      | same (StrandTy n1, StrandTy n2) = Atom.same(n1, n2)
      | same (KernelTy, KernelTy) = true
      | same (FieldTy, FieldTy) = true
      | same (FemData f1, FemData f2) = FemData.same(f1, f2)
      | same _ = false

    fun hash BoolTy = 0w1
      | hash StringTy = 0w2
      | hash IntTy = 0w3
      | hash (TensorTy (dd, NONE)) = List.foldl (fn (d, s) => 0w11 * Word.fromInt d + s) 0w5 dd
      | hash (TupleTy tys) = List.foldl (fn (ty, s) => hash ty + s) 0w7 tys
      | hash (SeqTy(ty, NONE)) = hash ty + 0w11
      | hash (SeqTy(ty, SOME n)) = Word.fromInt n * hash ty + 0w13
      | hash (ImageTy info) = 0w17 + 0w3 * ImageInfo.hash info
      | hash KernelTy = 0w19
      | hash FieldTy = 0w23
      | hash (StrandTy n) = Atom.hash n
      | hash (FemData(data)) = 0w29 + 0w31 * (FemData.hash data)
      | hash (TensorTy (dd, SOME j)) = 0w37 * (Word.fromInt j) + List.foldl (fn (d, s) => 0w11 * Word.fromInt d + s) 0w5 dd

    fun pre(NONE) =""
      | pre (SOME 1) = "interval "
      | pre (SOME k) = "aff[" ^ (Int.toString k) ^ "] "
    fun toString BoolTy = "bool"
      | toString StringTy = "string"
      | toString IntTy = "int"
      | toString (TensorTy([], a)) = (pre a) ^ "real"
      | toString (TensorTy (dd, a)) = String.concat[ (pre a),
            "tensor[", String.concatWithMap "," Int.toString dd, "]"
          ]
      | toString (TupleTy tys) = String.concat[
            "(", String.concatWithMap " * " toString tys, ")"
          ]
      | toString (SeqTy(ty, NONE)) = toString ty ^ "[]"
      | toString (SeqTy(ty, SOME n)) = concat[toString ty, "[", Int.toString n, "]"]
      | toString (ImageTy info) = concat["image(", ImageInfo.toString info, ")"]
      | toString (StrandTy n) = Atom.toString n
      | toString KernelTy = "kernel"
      | toString FieldTy = "field"
      | toString (FemData(data)) = "femData:" ^ (FemData.toString data)

  end
