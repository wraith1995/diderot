(* mid-ir.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *
 * Middle-level version of the Diderot CFG IR
 *
 * Note: this file is generated from gen/ir/mid-ir.spec and gen/ir/mid-ir.in.
 *)

structure MidOps =
  struct

  (* required helper functions for types *)
    type ty = MidTypes.ty
    val samety = MidTypes.same
    val hashty = MidTypes.hash
    val tyToString = MidTypes.toString

  (* required helper functions for type lists *)
    type tys = ty list
    fun sametys (tys1, tys2) = ListPair.allEq samety (tys1, tys2)
    fun hashtys tys = List.foldl (fn (ty, s) => hashty ty + 0w3 * s) 0w0 tys
    fun tysToString tys = String.concat["[", String.concatWithMap "," tyToString tys, "]" ]

  (* required helper functions for the int type *)
    fun sameint (i1 : int, i2) = (i1 = i2)
    fun hashint i = Word.fromInt i
    fun intToString i = Int.toString i

  (* required helper functions for the string type *)
    fun samestring (s1 : string, s2) = (s1 = s2)
    val hashstring = HashString.hashString
    fun stringToString s = String.concat["\"", s, "\""]

  (* required helper functions for the shape type *)
    type shape = TensorShape.t
    val sameshape = TensorShape.same
    val hashshape = TensorShape.hash
    val shapeToString = TensorShape.toString

  (* required helper functions for the index control type type *)
    type idxctl = IndexCtl.t
    val sameidxctl = IndexCtl.same
    val hashidxctl = IndexCtl.hash
    val idxctlToString = IndexCtl.toString

    datatype rator
      = IAdd
      | ISub
      | IMul
      | IDiv
      | IMod
      | INeg
      | IAbs
      | LT of ty
      | LTE of ty
      | EQ of ty
      | NEQ of ty
      | GT of ty
      | GTE of ty
      | BAnd
      | BOr
      | BNot
      | Max of ty
      | Min of ty
      | TensorIndex of ty * shape
      | EigenVecs2x2
      | EigenVecs3x3
      | EigenVals2x2
      | EigenVals3x3
      | Zero of ty
      | Select of ty * int
      | Tuple of tys
      | Subscript of ty
      | MkDynamic of ty * int
      | Append of ty
      | Prepend of ty
      | Concat of ty
      | Range
      | Length of ty
      | SphereQuery of int * ty
      | Ceiling of int
      | Floor of int
      | Round of int
      | Trunc of int
      | IntToReal
      | RealToInt of int
      | NumStrands of StrandSets.t
      | Strands of ty * StrandSets.t
      | BuildPos of int
      | EvaluateBasis of BasisDataArray.t
      | Func3 of Stamp.t
      | EvalKernel of int * Kernel.t * int
      | Kernel of Kernel.t * int
      | Transform of ImageInfo.t
      | Translate of ImageInfo.t
      | LoadVoxels of ImageInfo.t * int
      | LoadVoxelsBase of ImageInfo.t * int
      | LoadVoxelsWithCtl of ImageInfo.t * int * idxctl
      | LoadVoxelsBaseWithCtl of ImageInfo.t * int * idxctl
      | Inside of ImageInfo.t * int
      | IndexInside of ImageInfo.t * int
      | ImageDim of ImageInfo.t * int
      | BorderCtlDefault of ImageInfo.t
      | BorderCtlClamp of ImageInfo.t
      | BorderCtlMirror of ImageInfo.t
      | BorderCtlWrap of ImageInfo.t
      | LoadSeq of ty * string
      | LoadImage of ty * string
      | LoadFem of ty
      | ExtractFemItem of ty * FemOpt.femOption
      | ExtractFemItem2 of ty * ty * FemOpt.femOption
      | ExtractFem of ty * ty
      | ExtractFemItemN of tys * ty * FemOpt.femOption * Stamp.t * string * tys * ty
      | KillAll
      | StabilizeAll
      | Print of tys
      | MathFn of MathFns.t
      | Check of int
      | Load of int * int * ty * int
      | Save of int * int * ty * int
      | intervalSimple
      | intervalMixed
      | intervalAffine
      | intervalToAffine
      | tensorToAffine
      | affineNative of ty * ty * ty
      | affineNative2 of ty
      | errors
      | lasterr
      | center
      | radius
      | minInterval
      | maxInterval
      | intersection
      | hull
      | extend
      | insideInterval
      | EvaluateBasisAff of BasisDataArray.t * int
      | twocomp

    fun resultArity IAdd = 1
      | resultArity ISub = 1
      | resultArity IMul = 1
      | resultArity IDiv = 1
      | resultArity IMod = 1
      | resultArity INeg = 1
      | resultArity IAbs = 1
      | resultArity (LT _) = 1
      | resultArity (LTE _) = 1
      | resultArity (EQ _) = 1
      | resultArity (NEQ _) = 1
      | resultArity (GT _) = 1
      | resultArity (GTE _) = 1
      | resultArity BAnd = 1
      | resultArity BOr = 1
      | resultArity BNot = 1
      | resultArity (Max _) = 1
      | resultArity (Min _) = 1
      | resultArity (TensorIndex _) = 1
      | resultArity EigenVecs2x2 = 2
      | resultArity EigenVecs3x3 = 2
      | resultArity EigenVals2x2 = 1
      | resultArity EigenVals3x3 = 1
      | resultArity (Zero _) = 1
      | resultArity (Select _) = 1
      | resultArity (Tuple _) = 1
      | resultArity (Subscript _) = 1
      | resultArity (MkDynamic _) = 1
      | resultArity (Append _) = 1
      | resultArity (Prepend _) = 1
      | resultArity (Concat _) = 1
      | resultArity Range = 1
      | resultArity (Length _) = 1
      | resultArity (SphereQuery _) = 1
      | resultArity (Ceiling _) = 1
      | resultArity (Floor _) = 1
      | resultArity (Round _) = 1
      | resultArity (Trunc _) = 1
      | resultArity IntToReal = 1
      | resultArity (RealToInt _) = 1
      | resultArity (NumStrands _) = 1
      | resultArity (Strands _) = 1
      | resultArity (BuildPos _) = 1
      | resultArity (EvaluateBasis _) = 1
      | resultArity (Func3 _) = 1
      | resultArity (EvalKernel _) = 1
      | resultArity (Kernel _) = 1
      | resultArity (Transform _) = 1
      | resultArity (Translate _) = 1
      | resultArity (LoadVoxels _) = 1
      | resultArity (LoadVoxelsBase _) = 1
      | resultArity (LoadVoxelsWithCtl _) = 1
      | resultArity (LoadVoxelsBaseWithCtl _) = 1
      | resultArity (Inside _) = 1
      | resultArity (IndexInside _) = 1
      | resultArity (ImageDim _) = 1
      | resultArity (BorderCtlDefault _) = 1
      | resultArity (BorderCtlClamp _) = 1
      | resultArity (BorderCtlMirror _) = 1
      | resultArity (BorderCtlWrap _) = 1
      | resultArity (LoadSeq _) = 1
      | resultArity (LoadImage _) = 1
      | resultArity (LoadFem _) = 1
      | resultArity (ExtractFemItem _) = 1
      | resultArity (ExtractFemItem2 _) = 1
      | resultArity (ExtractFem _) = 1
      | resultArity (ExtractFemItemN _) = 1
      | resultArity KillAll = 0
      | resultArity StabilizeAll = 0
      | resultArity (Print _) = 0
      | resultArity (MathFn _) = 1
      | resultArity (Check _) = 1
      | resultArity (Load _) = 1
      | resultArity (Save _) = 0
      | resultArity intervalSimple = 1
      | resultArity intervalMixed = 1
      | resultArity intervalAffine = 1
      | resultArity intervalToAffine = 1
      | resultArity tensorToAffine = 1
      | resultArity (affineNative _) = 1
      | resultArity (affineNative2 _) = 1
      | resultArity errors = 1
      | resultArity lasterr = 1
      | resultArity center = 1
      | resultArity radius = 1
      | resultArity minInterval = 1
      | resultArity maxInterval = 1
      | resultArity intersection = 1
      | resultArity hull = 1
      | resultArity extend = 1
      | resultArity insideInterval = 1
      | resultArity (EvaluateBasisAff _) = 1
      | resultArity twocomp = 3

    fun arity IAdd = 2
      | arity ISub = 2
      | arity IMul = 2
      | arity IDiv = 2
      | arity IMod = 2
      | arity INeg = 1
      | arity IAbs = 1
      | arity (LT _) = 2
      | arity (LTE _) = 2
      | arity (EQ _) = 2
      | arity (NEQ _) = 2
      | arity (GT _) = 2
      | arity (GTE _) = 2
      | arity BAnd = 2
      | arity BOr = 2
      | arity BNot = 1
      | arity (Max _) = 2
      | arity (Min _) = 2
      | arity (TensorIndex _) = 1
      | arity EigenVecs2x2 = 1
      | arity EigenVecs3x3 = 1
      | arity EigenVals2x2 = 1
      | arity EigenVals3x3 = 1
      | arity (Zero _) = 0
      | arity (Select _) = 1
      | arity (Tuple _) = ~1
      | arity (Subscript _) = 2
      | arity (MkDynamic _) = 1
      | arity (Append _) = 2
      | arity (Prepend _) = 2
      | arity (Concat _) = 2
      | arity Range = 2
      | arity (Length _) = 1
      | arity (SphereQuery _) = 2
      | arity (Ceiling _) = 1
      | arity (Floor _) = 1
      | arity (Round _) = 1
      | arity (Trunc _) = 1
      | arity IntToReal = 1
      | arity (RealToInt _) = 1
      | arity (NumStrands _) = 0
      | arity (Strands _) = 0
      | arity (BuildPos _) = 1
      | arity (EvaluateBasis _) = 1
      | arity (Func3 _) = 3
      | arity (EvalKernel _) = 1
      | arity (Kernel _) = 0
      | arity (Transform _) = 1
      | arity (Translate _) = 1
      | arity (LoadVoxels _) = 2
      | arity (LoadVoxelsBase _) = 2
      | arity (LoadVoxelsWithCtl _) = 2
      | arity (LoadVoxelsBaseWithCtl _) = 2
      | arity (Inside _) = 2
      | arity (IndexInside _) = 2
      | arity (ImageDim _) = 1
      | arity (BorderCtlDefault _) = 2
      | arity (BorderCtlClamp _) = 1
      | arity (BorderCtlMirror _) = 1
      | arity (BorderCtlWrap _) = 1
      | arity (LoadSeq _) = 0
      | arity (LoadImage _) = 0
      | arity (LoadFem _) = 2
      | arity (ExtractFemItem _) = 1
      | arity (ExtractFemItem2 _) = 2
      | arity (ExtractFem _) = 1
      | arity (ExtractFemItemN _) = ~1
      | arity KillAll = 0
      | arity StabilizeAll = 0
      | arity (Print _) = ~1
      | arity (MathFn _) = ~1
      | arity (Check _) = 1
      | arity (Load _) = 0
      | arity (Save _) = 2
      | arity intervalSimple = 1
      | arity intervalMixed = 2
      | arity intervalAffine = 1
      | arity intervalToAffine = 1
      | arity tensorToAffine = 1
      | arity (affineNative _) = 3
      | arity (affineNative2 _) = 2
      | arity errors = 1
      | arity lasterr = 1
      | arity center = 1
      | arity radius = 1
      | arity minInterval = 1
      | arity maxInterval = 1
      | arity intersection = 2
      | arity hull = 2
      | arity extend = 2
      | arity insideInterval = 2
      | arity (EvaluateBasisAff _) = 1
      | arity twocomp = 1

    fun isPure (MkDynamic _) = false
      | isPure (Append _) = false
      | isPure (Prepend _) = false
      | isPure (Concat _) = false
      | isPure KillAll = false
      | isPure StabilizeAll = false
      | isPure (Print _) = false
      | isPure (Save _) = false
      | isPure _ = true

    fun same (IAdd, IAdd) = true
      | same (ISub, ISub) = true
      | same (IMul, IMul) = true
      | same (IDiv, IDiv) = true
      | same (IMod, IMod) = true
      | same (INeg, INeg) = true
      | same (IAbs, IAbs) = true
      | same (LT(a0), LT(b0)) = samety(a0, b0)
      | same (LTE(a0), LTE(b0)) = samety(a0, b0)
      | same (EQ(a0), EQ(b0)) = samety(a0, b0)
      | same (NEQ(a0), NEQ(b0)) = samety(a0, b0)
      | same (GT(a0), GT(b0)) = samety(a0, b0)
      | same (GTE(a0), GTE(b0)) = samety(a0, b0)
      | same (BAnd, BAnd) = true
      | same (BOr, BOr) = true
      | same (BNot, BNot) = true
      | same (Max(a0), Max(b0)) = samety(a0, b0)
      | same (Min(a0), Min(b0)) = samety(a0, b0)
      | same (TensorIndex(a0,a1), TensorIndex(b0,b1)) = samety(a0, b0) andalso sameshape(a1, b1)
      | same (EigenVecs2x2, EigenVecs2x2) = true
      | same (EigenVecs3x3, EigenVecs3x3) = true
      | same (EigenVals2x2, EigenVals2x2) = true
      | same (EigenVals3x3, EigenVals3x3) = true
      | same (Zero(a0), Zero(b0)) = samety(a0, b0)
      | same (Select(a0,a1), Select(b0,b1)) = samety(a0, b0) andalso sameint(a1, b1)
      | same (Tuple(a0), Tuple(b0)) = sametys(a0, b0)
      | same (Subscript(a0), Subscript(b0)) = samety(a0, b0)
      | same (MkDynamic(a0,a1), MkDynamic(b0,b1)) = samety(a0, b0) andalso sameint(a1, b1)
      | same (Append(a0), Append(b0)) = samety(a0, b0)
      | same (Prepend(a0), Prepend(b0)) = samety(a0, b0)
      | same (Concat(a0), Concat(b0)) = samety(a0, b0)
      | same (Range, Range) = true
      | same (Length(a0), Length(b0)) = samety(a0, b0)
      | same (SphereQuery(a0,a1), SphereQuery(b0,b1)) = sameint(a0, b0) andalso samety(a1, b1)
      | same (Ceiling(a0), Ceiling(b0)) = sameint(a0, b0)
      | same (Floor(a0), Floor(b0)) = sameint(a0, b0)
      | same (Round(a0), Round(b0)) = sameint(a0, b0)
      | same (Trunc(a0), Trunc(b0)) = sameint(a0, b0)
      | same (IntToReal, IntToReal) = true
      | same (RealToInt(a0), RealToInt(b0)) = sameint(a0, b0)
      | same (NumStrands(a0), NumStrands(b0)) = StrandSets.same(a0, b0)
      | same (Strands(a0,a1), Strands(b0,b1)) = samety(a0, b0) andalso StrandSets.same(a1, b1)
      | same (BuildPos(a0), BuildPos(b0)) = sameint(a0, b0)
      | same (EvaluateBasis(a0), EvaluateBasis(b0)) = BasisDataArray.same(a0, b0)
      | same (Func3(a0), Func3(b0)) = Stamp.same(a0, b0)
      | same (EvalKernel(a0,a1,a2), EvalKernel(b0,b1,b2)) = sameint(a0, b0) andalso Kernel.same(a1, b1) andalso sameint(a2, b2)
      | same (Kernel(a0,a1), Kernel(b0,b1)) = Kernel.same(a0, b0) andalso sameint(a1, b1)
      | same (Transform(a0), Transform(b0)) = ImageInfo.same(a0, b0)
      | same (Translate(a0), Translate(b0)) = ImageInfo.same(a0, b0)
      | same (LoadVoxels(a0,a1), LoadVoxels(b0,b1)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1)
      | same (LoadVoxelsBase(a0,a1), LoadVoxelsBase(b0,b1)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1)
      | same (LoadVoxelsWithCtl(a0,a1,a2), LoadVoxelsWithCtl(b0,b1,b2)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1) andalso sameidxctl(a2, b2)
      | same (LoadVoxelsBaseWithCtl(a0,a1,a2), LoadVoxelsBaseWithCtl(b0,b1,b2)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1) andalso sameidxctl(a2, b2)
      | same (Inside(a0,a1), Inside(b0,b1)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1)
      | same (IndexInside(a0,a1), IndexInside(b0,b1)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1)
      | same (ImageDim(a0,a1), ImageDim(b0,b1)) = ImageInfo.same(a0, b0) andalso sameint(a1, b1)
      | same (BorderCtlDefault(a0), BorderCtlDefault(b0)) = ImageInfo.same(a0, b0)
      | same (BorderCtlClamp(a0), BorderCtlClamp(b0)) = ImageInfo.same(a0, b0)
      | same (BorderCtlMirror(a0), BorderCtlMirror(b0)) = ImageInfo.same(a0, b0)
      | same (BorderCtlWrap(a0), BorderCtlWrap(b0)) = ImageInfo.same(a0, b0)
      | same (LoadSeq(a0,a1), LoadSeq(b0,b1)) = samety(a0, b0) andalso samestring(a1, b1)
      | same (LoadImage(a0,a1), LoadImage(b0,b1)) = samety(a0, b0) andalso samestring(a1, b1)
      | same (LoadFem(a0), LoadFem(b0)) = samety(a0, b0)
      | same (ExtractFemItem(a0,a1), ExtractFemItem(b0,b1)) = samety(a0, b0) andalso FemOpt.same(a1, b1)
      | same (ExtractFemItem2(a0,a1,a2), ExtractFemItem2(b0,b1,b2)) = samety(a0, b0) andalso samety(a1, b1) andalso FemOpt.same(a2, b2)
      | same (ExtractFem(a0,a1), ExtractFem(b0,b1)) = samety(a0, b0) andalso samety(a1, b1)
      | same (ExtractFemItemN(a0,a1,a2,a3,a4,a5,a6), ExtractFemItemN(b0,b1,b2,b3,b4,b5,b6)) = sametys(a0, b0) andalso samety(a1, b1) andalso FemOpt.same(a2, b2) andalso Stamp.same(a3, b3) andalso samestring(a4, b4) andalso sametys(a5, b5) andalso samety(a6, b6)
      | same (KillAll, KillAll) = true
      | same (StabilizeAll, StabilizeAll) = true
      | same (Print(a0), Print(b0)) = sametys(a0, b0)
      | same (MathFn(a0), MathFn(b0)) = MathFns.same(a0, b0)
      | same (Check(a0), Check(b0)) = sameint(a0, b0)
      | same (Load(a0,a1,a2,a3), Load(b0,b1,b2,b3)) = sameint(a0, b0) andalso sameint(a1, b1) andalso samety(a2, b2) andalso sameint(a3, b3)
      | same (Save(a0,a1,a2,a3), Save(b0,b1,b2,b3)) = sameint(a0, b0) andalso sameint(a1, b1) andalso samety(a2, b2) andalso sameint(a3, b3)
      | same (intervalSimple, intervalSimple) = true
      | same (intervalMixed, intervalMixed) = true
      | same (intervalAffine, intervalAffine) = true
      | same (intervalToAffine, intervalToAffine) = true
      | same (tensorToAffine, tensorToAffine) = true
      | same (affineNative(a0,a1,a2), affineNative(b0,b1,b2)) = samety(a0, b0) andalso samety(a1, b1) andalso samety(a2, b2)
      | same (affineNative2(a0), affineNative2(b0)) = samety(a0, b0)
      | same (errors, errors) = true
      | same (lasterr, lasterr) = true
      | same (center, center) = true
      | same (radius, radius) = true
      | same (minInterval, minInterval) = true
      | same (maxInterval, maxInterval) = true
      | same (intersection, intersection) = true
      | same (hull, hull) = true
      | same (extend, extend) = true
      | same (insideInterval, insideInterval) = true
      | same (EvaluateBasisAff(a0,a1), EvaluateBasisAff(b0,b1)) = BasisDataArray.same(a0, b0) andalso sameint(a1, b1)
      | same (twocomp, twocomp) = true
      | same _ = false

    fun hash IAdd = 0w3
      | hash ISub = 0w5
      | hash IMul = 0w7
      | hash IDiv = 0w11
      | hash IMod = 0w13
      | hash INeg = 0w17
      | hash IAbs = 0w19
      | hash (LT(a0)) = 0w23 + hashty a0
      | hash (LTE(a0)) = 0w29 + hashty a0
      | hash (EQ(a0)) = 0w31 + hashty a0
      | hash (NEQ(a0)) = 0w37 + hashty a0
      | hash (GT(a0)) = 0w41 + hashty a0
      | hash (GTE(a0)) = 0w43 + hashty a0
      | hash BAnd = 0w47
      | hash BOr = 0w53
      | hash BNot = 0w59
      | hash (Max(a0)) = 0w61 + hashty a0
      | hash (Min(a0)) = 0w67 + hashty a0
      | hash (TensorIndex(a0,a1)) = 0w71 + hashty a0 + hashshape a1
      | hash EigenVecs2x2 = 0w73
      | hash EigenVecs3x3 = 0w79
      | hash EigenVals2x2 = 0w83
      | hash EigenVals3x3 = 0w89
      | hash (Zero(a0)) = 0w97 + hashty a0
      | hash (Select(a0,a1)) = 0w101 + hashty a0 + hashint a1
      | hash (Tuple(a0)) = 0w103 + hashtys a0
      | hash (Subscript(a0)) = 0w107 + hashty a0
      | hash (MkDynamic(a0,a1)) = 0w109 + hashty a0 + hashint a1
      | hash (Append(a0)) = 0w113 + hashty a0
      | hash (Prepend(a0)) = 0w127 + hashty a0
      | hash (Concat(a0)) = 0w131 + hashty a0
      | hash Range = 0w137
      | hash (Length(a0)) = 0w139 + hashty a0
      | hash (SphereQuery(a0,a1)) = 0w149 + hashint a0 + hashty a1
      | hash (Ceiling(a0)) = 0w151 + hashint a0
      | hash (Floor(a0)) = 0w157 + hashint a0
      | hash (Round(a0)) = 0w163 + hashint a0
      | hash (Trunc(a0)) = 0w167 + hashint a0
      | hash IntToReal = 0w173
      | hash (RealToInt(a0)) = 0w179 + hashint a0
      | hash (NumStrands(a0)) = 0w181 + StrandSets.hash a0
      | hash (Strands(a0,a1)) = 0w191 + hashty a0 + StrandSets.hash a1
      | hash (BuildPos(a0)) = 0w193 + hashint a0
      | hash (EvaluateBasis(a0)) = 0w197 + BasisDataArray.hash a0
      | hash (Func3(a0)) = 0w199 + Stamp.hash a0
      | hash (EvalKernel(a0,a1,a2)) = 0w211 + hashint a0 + Kernel.hash a1 + hashint a2
      | hash (Kernel(a0,a1)) = 0w223 + Kernel.hash a0 + hashint a1
      | hash (Transform(a0)) = 0w227 + ImageInfo.hash a0
      | hash (Translate(a0)) = 0w229 + ImageInfo.hash a0
      | hash (LoadVoxels(a0,a1)) = 0w233 + ImageInfo.hash a0 + hashint a1
      | hash (LoadVoxelsBase(a0,a1)) = 0w239 + ImageInfo.hash a0 + hashint a1
      | hash (LoadVoxelsWithCtl(a0,a1,a2)) = 0w241 + ImageInfo.hash a0 + hashint a1 + hashidxctl a2
      | hash (LoadVoxelsBaseWithCtl(a0,a1,a2)) = 0w251 + ImageInfo.hash a0 + hashint a1 + hashidxctl a2
      | hash (Inside(a0,a1)) = 0w257 + ImageInfo.hash a0 + hashint a1
      | hash (IndexInside(a0,a1)) = 0w263 + ImageInfo.hash a0 + hashint a1
      | hash (ImageDim(a0,a1)) = 0w269 + ImageInfo.hash a0 + hashint a1
      | hash (BorderCtlDefault(a0)) = 0w271 + ImageInfo.hash a0
      | hash (BorderCtlClamp(a0)) = 0w277 + ImageInfo.hash a0
      | hash (BorderCtlMirror(a0)) = 0w281 + ImageInfo.hash a0
      | hash (BorderCtlWrap(a0)) = 0w283 + ImageInfo.hash a0
      | hash (LoadSeq(a0,a1)) = 0w293 + hashty a0 + hashstring a1
      | hash (LoadImage(a0,a1)) = 0w307 + hashty a0 + hashstring a1
      | hash (LoadFem(a0)) = 0w311 + hashty a0
      | hash (ExtractFemItem(a0,a1)) = 0w313 + hashty a0 + FemOpt.hash a1
      | hash (ExtractFemItem2(a0,a1,a2)) = 0w317 + hashty a0 + hashty a1 + FemOpt.hash a2
      | hash (ExtractFem(a0,a1)) = 0w331 + hashty a0 + hashty a1
      | hash (ExtractFemItemN(a0,a1,a2,a3,a4,a5,a6)) = 0w337 + hashtys a0 + hashty a1 + FemOpt.hash a2 + Stamp.hash a3 + hashstring a4 + hashtys a5 + hashty a6
      | hash KillAll = 0w347
      | hash StabilizeAll = 0w349
      | hash (Print(a0)) = 0w353 + hashtys a0
      | hash (MathFn(a0)) = 0w359 + MathFns.hash a0
      | hash (Check(a0)) = 0w367 + hashint a0
      | hash (Load(a0,a1,a2,a3)) = 0w373 + hashint a0 + hashint a1 + hashty a2 + hashint a3
      | hash (Save(a0,a1,a2,a3)) = 0w379 + hashint a0 + hashint a1 + hashty a2 + hashint a3
      | hash intervalSimple = 0w383
      | hash intervalMixed = 0w389
      | hash intervalAffine = 0w397
      | hash intervalToAffine = 0w401
      | hash tensorToAffine = 0w409
      | hash (affineNative(a0,a1,a2)) = 0w419 + hashty a0 + hashty a1 + hashty a2
      | hash (affineNative2(a0)) = 0w421 + hashty a0
      | hash errors = 0w431
      | hash lasterr = 0w433
      | hash center = 0w439
      | hash radius = 0w443
      | hash minInterval = 0w449
      | hash maxInterval = 0w457
      | hash intersection = 0w461
      | hash hull = 0w463
      | hash extend = 0w467
      | hash insideInterval = 0w479
      | hash (EvaluateBasisAff(a0,a1)) = 0w487 + BasisDataArray.hash a0 + hashint a1
      | hash twocomp = 0w491

    fun toString IAdd = "IAdd"
      | toString ISub = "ISub"
      | toString IMul = "IMul"
      | toString IDiv = "IDiv"
      | toString IMod = "IMod"
      | toString INeg = "INeg"
      | toString IAbs = "IAbs"
      | toString (LT(a0)) = concat["LT<", tyToString a0, ">"]
      | toString (LTE(a0)) = concat["LTE<", tyToString a0, ">"]
      | toString (EQ(a0)) = concat["EQ<", tyToString a0, ">"]
      | toString (NEQ(a0)) = concat["NEQ<", tyToString a0, ">"]
      | toString (GT(a0)) = concat["GT<", tyToString a0, ">"]
      | toString (GTE(a0)) = concat["GTE<", tyToString a0, ">"]
      | toString BAnd = "BAnd"
      | toString BOr = "BOr"
      | toString BNot = "BNot"
      | toString (Max(a0)) = concat["Max<", tyToString a0, ">"]
      | toString (Min(a0)) = concat["Min<", tyToString a0, ">"]
      | toString (TensorIndex(a0,a1)) = concat["TensorIndex<", tyToString a0, ",", shapeToString a1, ">"]
      | toString EigenVecs2x2 = "EigenVecs2x2"
      | toString EigenVecs3x3 = "EigenVecs3x3"
      | toString EigenVals2x2 = "EigenVals2x2"
      | toString EigenVals3x3 = "EigenVals3x3"
      | toString (Zero(a0)) = concat["Zero<", tyToString a0, ">"]
      | toString (Select(a0,a1)) = concat["Select<", tyToString a0, ",", intToString a1, ">"]
      | toString (Tuple(a0)) = concat["Tuple<", tysToString a0, ">"]
      | toString (Subscript(a0)) = concat["Subscript<", tyToString a0, ">"]
      | toString (MkDynamic(a0,a1)) = concat["MkDynamic<", tyToString a0, ",", intToString a1, ">"]
      | toString (Append(a0)) = concat["Append<", tyToString a0, ">"]
      | toString (Prepend(a0)) = concat["Prepend<", tyToString a0, ">"]
      | toString (Concat(a0)) = concat["Concat<", tyToString a0, ">"]
      | toString Range = "Range"
      | toString (Length(a0)) = concat["Length<", tyToString a0, ">"]
      | toString (SphereQuery(a0,a1)) = concat["SphereQuery<", intToString a0, ",", tyToString a1, ">"]
      | toString (Ceiling(a0)) = concat["Ceiling<", intToString a0, ">"]
      | toString (Floor(a0)) = concat["Floor<", intToString a0, ">"]
      | toString (Round(a0)) = concat["Round<", intToString a0, ">"]
      | toString (Trunc(a0)) = concat["Trunc<", intToString a0, ">"]
      | toString IntToReal = "IntToReal"
      | toString (RealToInt(a0)) = concat["RealToInt<", intToString a0, ">"]
      | toString (NumStrands(a0)) = concat["NumStrands<", StrandSets.toString a0, ">"]
      | toString (Strands(a0,a1)) = concat["Strands<", tyToString a0, ",", StrandSets.toString a1, ">"]
      | toString (BuildPos(a0)) = concat["BuildPos<", intToString a0, ">"]
      | toString (EvaluateBasis(a0)) = concat["EvaluateBasis<", BasisDataArray.toString a0, ">"]
      | toString (Func3(a0)) = concat["Func3<", Stamp.toString a0, ">"]
      | toString (EvalKernel(a0,a1,a2)) = concat["EvalKernel<", intToString a0, ",", Kernel.toString a1, ",", intToString a2, ">"]
      | toString (Kernel(a0,a1)) = concat["Kernel<", Kernel.toString a0, ",", intToString a1, ">"]
      | toString (Transform(a0)) = concat["Transform<", ImageInfo.toString a0, ">"]
      | toString (Translate(a0)) = concat["Translate<", ImageInfo.toString a0, ">"]
      | toString (LoadVoxels(a0,a1)) = concat["LoadVoxels<", ImageInfo.toString a0, ",", intToString a1, ">"]
      | toString (LoadVoxelsBase(a0,a1)) = concat["LoadVoxelsBase<", ImageInfo.toString a0, ",", intToString a1, ">"]
      | toString (LoadVoxelsWithCtl(a0,a1,a2)) = concat["LoadVoxelsWithCtl<", ImageInfo.toString a0, ",", intToString a1, ",", idxctlToString a2, ">"]
      | toString (LoadVoxelsBaseWithCtl(a0,a1,a2)) = concat["LoadVoxelsBaseWithCtl<", ImageInfo.toString a0, ",", intToString a1, ",", idxctlToString a2, ">"]
      | toString (Inside(a0,a1)) = concat["Inside<", ImageInfo.toString a0, ",", intToString a1, ">"]
      | toString (IndexInside(a0,a1)) = concat["IndexInside<", ImageInfo.toString a0, ",", intToString a1, ">"]
      | toString (ImageDim(a0,a1)) = concat["ImageDim<", ImageInfo.toString a0, ",", intToString a1, ">"]
      | toString (BorderCtlDefault(a0)) = concat["BorderCtlDefault<", ImageInfo.toString a0, ">"]
      | toString (BorderCtlClamp(a0)) = concat["BorderCtlClamp<", ImageInfo.toString a0, ">"]
      | toString (BorderCtlMirror(a0)) = concat["BorderCtlMirror<", ImageInfo.toString a0, ">"]
      | toString (BorderCtlWrap(a0)) = concat["BorderCtlWrap<", ImageInfo.toString a0, ">"]
      | toString (LoadSeq(a0,a1)) = concat["LoadSeq<", tyToString a0, ",", stringToString a1, ">"]
      | toString (LoadImage(a0,a1)) = concat["LoadImage<", tyToString a0, ",", stringToString a1, ">"]
      | toString (LoadFem(a0)) = concat["LoadFem<", tyToString a0, ">"]
      | toString (ExtractFemItem(a0,a1)) = concat["ExtractFemItem<", tyToString a0, ",", FemOpt.toString a1, ">"]
      | toString (ExtractFemItem2(a0,a1,a2)) = concat["ExtractFemItem2<", tyToString a0, ",", tyToString a1, ",", FemOpt.toString a2, ">"]
      | toString (ExtractFem(a0,a1)) = concat["ExtractFem<", tyToString a0, ",", tyToString a1, ">"]
      | toString (ExtractFemItemN(a0,a1,a2,a3,a4,a5,a6)) = concat["ExtractFemItemN<", tysToString a0, ",", tyToString a1, ",", FemOpt.toString a2, ",", Stamp.toString a3, ",", stringToString a4, ",", tysToString a5, ",", tyToString a6, ">"]
      | toString KillAll = "KillAll"
      | toString StabilizeAll = "StabilizeAll"
      | toString (Print(a0)) = concat["Print<", tysToString a0, ">"]
      | toString (MathFn(a0)) = concat["MathFn<", MathFns.toString a0, ">"]
      | toString (Check(a0)) = concat["Check<", intToString a0, ">"]
      | toString (Load(a0,a1,a2,a3)) = concat["Load<", intToString a0, ",", intToString a1, ",", tyToString a2, ",", intToString a3, ">"]
      | toString (Save(a0,a1,a2,a3)) = concat["Save<", intToString a0, ",", intToString a1, ",", tyToString a2, ",", intToString a3, ">"]
      | toString intervalSimple = "intervalSimple"
      | toString intervalMixed = "intervalMixed"
      | toString intervalAffine = "intervalAffine"
      | toString intervalToAffine = "intervalToAffine"
      | toString tensorToAffine = "tensorToAffine"
      | toString (affineNative(a0,a1,a2)) = concat["affineNative<", tyToString a0, ",", tyToString a1, ",", tyToString a2, ">"]
      | toString (affineNative2(a0)) = concat["affineNative2<", tyToString a0, ">"]
      | toString errors = "errors"
      | toString lasterr = "lasterr"
      | toString center = "center"
      | toString radius = "radius"
      | toString minInterval = "minInterval"
      | toString maxInterval = "maxInterval"
      | toString intersection = "intersection"
      | toString hull = "hull"
      | toString extend = "extend"
      | toString insideInterval = "insideInterval"
      | toString (EvaluateBasisAff(a0,a1)) = concat["EvaluateBasisAff<", BasisDataArray.toString a0, ",", intToString a1, ">"]
      | toString twocomp = "twocomp"

  end

structure MidIR = SSAFn(
  val irName = "mid-ir"
  structure Ty = MidTypes
  structure Op = MidOps)

structure MidCensus = CensusFn(MidIR)
