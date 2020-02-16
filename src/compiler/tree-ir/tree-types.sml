(* tree-types.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *
 * Types for the TreeIR.  For now, these are identical to the LowIR types.
 *)

structure TreeTypes =
  struct

  (* the layout of a vector onto target-supported vectors. *)
    type vec_layout = VectorLayout.t

    datatype t
      = BoolTy | StringTy | IntTy
      | VecTy of int * int              (* a vector value that fits into a target-supported
                                         * vector. The first int is the width of the vector
                                         * and the second is the target width, which will be
                                         * >= the vector width.
                                         *)
      | TensorTy of int list            (* in-memory tensor type.  For 0th order and 1st
                                         * order tensors, the local-variable types will
                                         * be RealTy and VecTy (resp.).
                                         *)
      | TensorRefTy of int list         (* a shared reference to tensor *)
      | TupleTy of t list               (* tuples; used for multiple return values *)
      | SeqTy of t * int option
      | ImageTy of ImageInfo.t
      | StrandIdTy of Atom.atom         (* index into strand-state array *)
      | FemData of FemData.femType

    fun containsFem (FemData(_)) = true
      | containsFem (SeqTy(ty, _)) = containsFem ty
      | containsFem (TupleTy(tys)) = List.foldr (fn (x,y) => (containsFem x) orelse y) false tys
      | containsFem _ = false

    fun findFem (ty as FemData(_)) = SOME(ty)
      | findFem (SeqTy(ty', _)) = findFem ty'
      | findFem _ = NONE
					  

    fun replaceFem (SeqTy(ty, s)) = SeqTy(replaceFem ty, s)
      | replaceFem (FemData(FemData.MeshCell(_))) = IntTy
      | replaceFem (FemData(FemData.FuncCell(_))) = IntTy
      | replaceFem (FemData(FemData.MeshPos(data))) =
	let val dim = FemData.meshDim data
	in TensorTy([dim * 2 + 1]) end
      | replaceFem _ = raise Fail "impossible"

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
			     
    fun fromAPI (APITypes.IntTy) = IntTy
      | fromAPI (APITypes.BoolTy) = BoolTy
      | fromAPI (APITypes.TensorTy(s)) = TensorTy(s) (*ignore possible issue*)
      | fromAPI (APITypes.StringTy) = StringTy
      | fromAPI (APITypes.FemData(data)) = FemData(data)
      | fromAPI (APITypes.SeqTy(ty, opt)) = SeqTy(fromAPI ty, opt)
      | fromAPI (APITypes.ImageTy _) = raise Fail "impossible"

    (* (*FIXME: VecTys and TensorRefTy*) *)
    (* fun toAPI (IntTy) = APITypes.IntTy *)
    (*   | toAPI (BoolTy) = APITypes.BoolTy *)
    (*   | toAPI (TensorTy(s)) = APITypes.TensorTy(s) *)
    (*   | toAPI (TensorRefTy(s)) = APITypes.TensorTy(s) (*WARNING*) *)
    (*   | toAPI (FemData(data)) = APITypes.FemData(data) *)
    (*   | toAPI (SeqTy(ty, opt)) = (APITypes.SeqTy(ty, opt)) *)
    (*   | toAPI _ = raise Fail "impossible" *)

    val intTy = IntTy

    val realTy = VecTy(1, 1)

    fun iVecTy 1 = IntTy
      | iVecTy n = SeqTy(IntTy, SOME n)

  (* return the vector types that comprise the pieces of a vec_layout *)
    fun piecesOf (layout : vec_layout) = let
          fun toVecs (_, []) = raise Fail(concat["piecesOf(", VectorLayout.toString layout, ")"])
            | toVecs (wid, [w]) = [VecTy(wid, w)]
            | toVecs (wid, w::ws) = VecTy(w, w) :: toVecs(wid-w, ws)
          in
            toVecs (#wid layout, #pieces layout)
          end

  (* return the unpadded width of the n'th component of a vector layout *)
    fun nthWidth (layout : vec_layout, n) = let
     fun get (_, _, []) = raise Fail(concat["nthWidth(", VectorLayout.toString layout, ", ", Int.toString n, ")"])
            | get (0, wid, _) = wid
            | get (n, wid, w::ws) = get (n-1, wid-w, ws)
	  val ret = get (n, #wid layout, #pieces layout)


    in
     case (List.nth(piecesOf(layout), n))
      of VecTy(wid, w) => wid
    end

  (* return the n'th component of a vector layout as a vector type. *)
    fun nthVec (layout : vec_layout, n) = let
          fun get (_, _, []) = raise Fail(concat[
                  "nthVec(", VectorLayout.toString layout, ", ", Int.toString n, ")"
                ])
            | get (0, wid, w::_) = VecTy(wid, w)
            | get (n, wid, w::ws) = get (n-1, wid-w, ws)
          in
           (* get (n, #wid layout, #pieces layout) *)
	   List.nth(piecesOf(layout), n)
          end

  (* return the offset of the first component in the i'th piece of a vec_layout *)
    fun offsetOf (layout : vec_layout, i) = let
          fun add (0, _, off) = off
            | add (i, pw::pws, off) = add (i-1, pws, off+pw)
          in
            add (i, #pieces layout, 0)
          end

  (* is a vector padded? *)
    fun isPaddedVec (VecTy(wid, hwWid)) = (wid <> hwWid)
      | isPaddedVec _ = false

  (* are two type structurally equal? *)
    fun same (BoolTy, BoolTy) = true
      | same (StringTy, StringTy) = true
      | same (IntTy, IntTy) = true
      | same (VecTy(wid1, hwWid1), VecTy(wid2, hwWid2)) = (wid1 = wid2) andalso (hwWid1 = hwWid2)
      | same (TensorTy dd1, TensorTy dd2) = ListPair.allEq (op =) (dd1, dd2)
      | same (TensorRefTy dd1, TensorRefTy dd2) = ListPair.allEq (op =) (dd1, dd2)
      | same (TupleTy tys1, TupleTy tys2) = ListPair.allEq same (tys1, tys2)
      | same (SeqTy(ty1, NONE), SeqTy(ty2, NONE)) = same(ty1, ty2)
      | same (SeqTy(ty1, SOME n1), SeqTy(ty2, SOME n2)) = (n1 = n2) andalso same(ty1, ty2)
      | same (ImageTy info1, ImageTy info2) = ImageInfo.sameShape(info1, info2)
      | same (StrandIdTy n1, StrandIdTy n2) = Atom.same(n1, n2)
      | same (FemData(data1), FemData(data2)) = FemData.same(data1, data2)
      | same _ = false

  (* is a source type compatible with the destination type? *)
    fun match {src = TensorTy dd1, dst = TensorRefTy dd2} = ListPair.allEq (op =) (dd1, dd2)
      | match {src, dst} = same(src, dst)

    fun hash BoolTy = 0w1
      | hash StringTy = 0w2
      | hash IntTy = 0w3
      | hash (VecTy(wid, hwWid)) = 0w5 * Word.fromInt wid + 0w7 * Word.fromInt hwWid
      | hash (TensorTy dd) = List.foldl (fn (d, s) => 0w11 * Word.fromInt d + s) 0w13 dd
      | hash (TensorRefTy dd) = List.foldl (fn (d, s) => 0w13 * Word.fromInt d + s) 0w17 dd
      | hash (TupleTy tys) = List.foldl (fn (ty, s) => hash ty + s) 0w19 tys
      | hash (SeqTy(ty, NONE)) = hash ty + 0w23
      | hash (SeqTy(ty, SOME n)) = Word.fromInt n * hash ty + 0w29
      | hash (ImageTy info) = 0w37 * ImageInfo.hash info + 0w6
      | hash (StrandIdTy n) = 0w41 + Atom.hash n
      | hash (FemData(data)) = 0w43 + 0w47 * (FemData.hash data)

    fun toString BoolTy = "bool"
      | toString StringTy = "string"
      | toString IntTy = "int"
      | toString (VecTy(1,1)) = "real"
      | toString (VecTy(wid, hwWid)) = if (wid = hwWid)
          then "vec" ^ Int.toString wid
          else concat["vec", Int.toString wid, "{", Int.toString hwWid, "}"]
      | toString (TensorTy dd) = String.concat [
            "tensor[", String.concatWithMap "," Int.toString dd, "]"
          ]
      | toString (TensorRefTy dd) = String.concat [
            "&tensor[", String.concatWithMap "," Int.toString dd, "]"
          ]
      | toString (TupleTy tys) = String.concat [
            "(", String.concatWithMap " * " toString tys, ")"
          ]
      | toString (SeqTy(ty, NONE)) = toString ty ^ "[]"
      | toString (SeqTy(ty, SOME n)) = concat[toString ty, "[", Int.toString n, "]"]
      | toString (ImageTy info) = concat["image(", ImageInfo.toString info, ")"]
      | toString (StrandIdTy n) = concat["id(", Atom.toString n, ")"]
      | toString (FemData(data)) = concat["femData(", FemData.toString(data), ")"]

    structure Tbl = HashTableFn (
      struct
        type hash_key = t
        val hashVal = hash
        val sameKey = same
      end)

  end
