(* mid-contract.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *
 * Expand intervals:
scalars = consts, scalar tensor
intervals = interval tensor
termains = intervals | scalars
ops = op1 terminas | op2 terminas x 2 | op3 terminas x 3 | opn terminas x end
ret = ops | E.Sum ops

ONE ISSUE: Make sure sizing and special ops
(case op1
      increasing, decreasing (neg, expr, sqrt, pow) -> obvious
      abs -> min(max(min, 0), max), max(abs(max), abs(min)) -- vectorized abs via or max(-xmin, max(xmax, -xmax)) 3 -- 4
      
      rest - two 2 vec rep, unwrap into interval functions, rewrap
      ignore: sgn

      op2:
      Sub: clear - divide into mins and maxes - make sure to do scalars
      div: make the division tensor, nan test (vectorized), and run multiplication - make sure to do scalars
      ignore: min, max

      op3: clamp ignored
      opn: ADD: clear, but scalars
      opn: PROD:
      --build the 2^n combinations for the n interval operations, add in the scalars, take mins and maxes as vectors.
      )

Sum: just substitute back in
SUM: any chances for vectorization we miss?

AA:
(case op1
      Just use interval versions! Avoid pow when you fix that and in probe

      op2:
      sub - clear, but again the filter
      div: convert to interval, use nan test, convert back to affine -- make div clearer (LAST)

      op3: ignore:

      opn:ADD: clear, but again the filter
      opn: PROD: for the operations, we lift and do all the non scalar ones. Then do the scalar ones. Then combine.
      
      )

SUM: you just subtitute back in and do it a few times. Easy and peasy.
SUM: any chances for vectorization we miss?

NOTE: what about scalar ops in the div, add, 
Plan: fix probe, power
Plan: ops to mmid
Plan: make util functions we will need and framework (we can use mk ops though)
Plan: Do interval, push the new functions down to low and vectorize and cxx
Plan: affine later.
Plan: Push and test
Plan:vectorization
PLan:fast-math
Plan:parser
.
 *)

structure MidIntervalExpand : sig

    val transform : MidIR.program -> MidIR.program

  end = struct

    structure IR = MidIR
    structure Op = MidOps
    structure Ty = MidTypes
    structure V = IR.Var
    structure ST = Stats
    structure E = Ein
    structure EU = EinUtil
    structure EP = EinPP
    (*given a tensor shape and an alpha, compute a new shape and alpha to represent the result after alpha is used in another op*)
    (*
\sum ijk...

*)
    fun alphaIndex(argshape, alpha) =
	let
	 val rank = List.length argshape
	 val newShape = List.mapPartial (fn E.V x => SOME(List.nth(argshape, x)) | E.C _ => NONE ) alpha
	 val newAlpha = List.filter (fn E.V _ => true | E.C _ => false) alpha
				     
	in
	 (newShape, newAlpha)
	end

    fun groupConversionS name avail (lst, s) = let
     val v::vst = lst
     val ty = IR.Var.ty v
     val Ty.TensorTy(shape', _) = ty
    in
     AvailRHS.assignCons(avail, name, lst, Ty.TensorTy(s::shape', NONE))
    end

		     
    fun cvtTy (Ty.TensorTy(ts, SOME j)) = if j = 1
					  then Ty.TensorTy(2 :: ts, NONE)
					  else if j > 1
					  then Ty.TensorTy(j :: ts, NONE)
					  else raise Fail "impossible Tensor Ty"
      | cvtTy (Ty.TensorTy(ts, NONE)) = Ty.TensorTy(ts, NONE)
      | cvtTy (Ty.TupleTy(ts)) = Ty.TupleTy(List.map cvtTy ts)
      | cvtTy (Ty.SeqTy(t, opt)) = Ty.SeqTy(cvtTy t, opt)
      | cvtTy (t) = t (*bool, string, int, img, strand, fem, kern*)

    structure Env = TranslateEnvFn (
     struct
     structure SrcIR = IR
     structure DstIR = IR
     val cvtTy = cvtTy
     end)		    
		     

  (********** Counters for statistics **********)
    val cntUnused               = ST.newCounter "mid-opt:unusedp"
    val firstCounter            = cntUnused
    val lastCounter             = cntUnused

    (*Plan: opts here, expands here, deal with new mid in mid to low, vectorization, c++, lifting, parsing, affine fix
     *)
    (*plan: probe and such - ops expand here -  further expansion (TOM)
Do assign, dops
     ops first because easy - tensor[j::sigma] -- low opt tensor index*)

    fun intervalCons(avail, a,b) =
	let
	 val tenTy = IR.Var.ty a
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val tenTy' = Ty.TensorTy(2::ts, NONE)
	in
	 AvailRHS.assignCons(avail, "intervalConst", [a,b], tenTy')
	end

    fun affine2Cons(avail, a,b) =
	let
	 val tenTy = IR.Var.ty a
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val tenTy' = Ty.TensorTy(2::ts, NONE)
	in
	 AvailRHS.assignCons(avail, "affine2Cons", [a,b], tenTy')
	end
	  

    fun affineCons(avail, v : IR.var, errs : IR.var list, errn : IR.var) =
	let
	 val ty = IR.Var.ty v
	 val tys = List.map (IR.Var.ty) (errn::errs)
	 val _ = if List.all (fn x => Ty.same(x, ty)) tys
		 then ()
		 else raise Fail "bad affine cons"
	 val Ty.TensorTy(shp, NONE) = ty
	 val total = 2 + List.length errs
	in
	 AvailRHS.assignCons(avail, "affineCons", v::(errs @ [errn]), Ty.TensorTy(total :: shp, NONE))
	end

    fun sliceOut(avail, outerIdx, src) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val outer::inner = ts
	 val _ = if 0 <= outerIdx andalso outerIdx <= outer
		 then ()
		 else raise Fail "bad interval expansion"
	 val _ = if List.length inner = 0
		 then raise Fail "bad interval expansion"
		 else ()

	 local
	  val null = ArrayNd.array'(0, inner)
	  val indexInfo = ArrayNd.fromArray'(ArrayNd.getInverseIndex(ArrayNd.computeInverseIndex(null)), inner)
	  fun tenIdxFun(shape) = AvailRHS.assignOp(avail, "intervalRebuildAcc", Ty.realTy, Op.TensorIndex(tenTy, outerIdx::shape), [src])
	  fun groupConversion(lst, s) = let
	   val v::vst = lst
	   val ty = IR.Var.ty v
	   val Ty.TensorTy(shape', _) = ty
	  in
	   AvailRHS.assignCons(avail, "intervalInterCons", lst, Ty.TensorTy(s::shape', NONE))
	  end
	 in
	 val sliced = ArrayNd.convertToTree(indexInfo, tenIdxFun, groupConversion)
	 end
	in
	 sliced
	end

    fun sliceOutEin(avail, outerIdx, src) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val outer::inner = ts
	 val _ = if 0 <= outerIdx andalso outerIdx <= outer
		 then ()
		 else raise Fail "bad interval expansion"
	 val _ = if List.length inner = 0
		 then raise Fail "bad interval expansion"
		 else ()

	 val params = [E.TEN(true, ts, NONE)]
	 val index = inner
	 val body = E.Tensor(0, E.C outer :: (List.tabulate(List.length inner, fn x => E.V x)))
	 val ein = E.EIN{params=params, index=index, body=body}

	in
	 AvailRHS.assignEin(avail, "intervalInterConsEin", Ty.TensorTy(inner, NONE), ein, [src])
	end


    fun makeAbsTen(avail, src) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val params = [E.TEN(true, ts, NONE)]
	 val index = ts
	 val body = E.Op1(E.Abs, E.Tensor(0, (List.tabulate(List.length ts, fn x => E.V x))))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "absTensor", tenTy, ein, [src])
	end
	  

    fun makeAbsTenAtIdx(avail, outerIdx, src) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val outer::inner = ts
	 val _ = if 0 <= outerIdx andalso outerIdx <= outer
		 then ()
		 else raise Fail "bad interval expansion"
	 val _ = if List.length inner = 0
		 then raise Fail "bad interval expansion"
		 else ()

	 val params = [E.TEN(true, ts, NONE)]
	 val index = inner
	 val body = E.Op1(E.Abs, E.Tensor(0, (E.C outerIdx) :: (List.tabulate(List.length inner, fn x => E.V x))))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "absTensorC", tenTy, ein, [src])
	end	  

    fun makeAddTen(avail, srcs) =
	let
	 val tys = List.map (IR.Var.ty) srcs
	 val Ty.TensorTy(ts, NONE) :: tys' = tys
	 val _ = if List.all (fn x => Ty.same(Ty.TensorTy(ts, NONE), x)) tys'
		 then ()
		 else raise Fail "bad tensor add"
	 val len = List.length ts
	 val tabed = List.tabulate(len, fn x => E.V x)
	 val params = List.tabulate(len, fn x => E.TEN(true, ts, NONE))
	 val tensors = List.tabulate(len, fn x => E.Tensor(x, tabed))
	 val index = ts
	 val body = E.Opn(E.Add, tensors)
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "addTensor", Ty.TensorTy(ts, NONE), ein, srcs)
	end

    fun makeSubTen(avail, src1, src2) =
	let
	 val srcs = [src1, src2]
	 val _ = if Ty.same(IR.Var.ty src1, IR.Var.ty src2)
		 then ()
		 else raise Fail "bad tensor sub"
	 val Ty.TensorTy(ts, NONE) = IR.Var.ty src1
	 val len = List.length ts
	 val tabed = List.tabulate(len, fn x => E.V x)
	 val params = [E.TEN(true, ts, NONE), E.TEN(true, ts, NONE)]
	 val index = ts
	 val body = E.Op2(E.Sub, E.Tensor(0, tabed), E.Tensor(1, tabed))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "subTensor", Ty.TensorTy(ts, NONE), ein, srcs)
	end
    fun makeIntervalWidth(avail, src) =
	let
	 val Ty.TensorTy(t::ts, NONE) = IR.Var.ty src
	 val _ = if t <> 2
		 then raise Fail "bad interval width"
		 else ()
	 val len = List.length ts
	 val tabed = List.tabulate(len, fn x => E.V x)
	 val params = [E.TEN(true, t::ts, NONE)]
	 val index = ts
	 val body = E.Op2(E.Sub, E.Tensor(0, E.C 1 :: tabed), E.Tensor(0, E.C 0 :: tabed))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "subTensor", Ty.TensorTy(ts, NONE), ein, [src])
	end

    fun makeIntervalMidPoint(avail, src) =
	let
	 val Ty.TensorTy(t::ts, NONE) = IR.Var.ty src
	 val _ = if t <> 2
		 then raise Fail "bad interval width"
		 else ()
	 val len = List.length ts
	 val tabed = List.tabulate(len, fn x => E.V x)
	 val params = [E.TEN(true, t::ts, NONE)]
	 val index = ts
	 val body = E.Opn(E.Add, [E.Tensor(0, E.C 0 :: tabed), E.Tensor(0, E.C 1 :: tabed)])
	 val ein = E.EIN{params=params, index=index, body=body}
	 val einVar = AvailRHS.assignEin(avail, "addIntervalTensor", Ty.TensorTy(ts, NONE), ein, [src])

	 val one = Rational.fromInt 1
	 val rational = Rational.fromInt 2
	 val dived =  (Rational.div(one, rational)) 
	 val params' = [E.TEN(true, ts, NONE)]	
	 val body' = E.Opn(E.Prod, [E.ConstR dived, E.Tensor(0, tabed)])
	 val ein' = E.EIN{params=params', index=index, body=body'}
	in
	 AvailRHS.assignEin(avail, "averageTensor", Ty.TensorTy(ts, NONE), ein', [einVar])
	end

    fun intervalZeros(avail, src) =
	let
	 val Ty.TensorTy(t::ts, NONE) = IR.Var.ty src
	 val _ = if t <> 2
		 then raise Fail "bad interval width"
		 else ()
	in
	 AvailRHS.makeRealLit(avail, ts, RealLit.zero true)
	end

    fun vmax2(b, avail, src1, src2) =
	let
	 val srcs = [src1, src2]
	 val _ = if Ty.same(IR.Var.ty src1, IR.Var.ty src2)
		 then ()
		 else raise Fail "bad tensor sub"
	 val Ty.TensorTy(ts, NONE) = IR.Var.ty src1
	 val len = List.length ts
	 val tabed = List.tabulate(len, fn x => E.V x)
	 val params = [E.TEN(true, ts, NONE), E.TEN(true, ts, NONE)]
	 val index = ts
	 val body = E.Op2(if b then E.Max else E.Min, E.Tensor(0, tabed), E.Tensor(1, tabed))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "maxTensor", Ty.TensorTy(ts, NONE), ein, srcs)
	end
    fun vmax2alpha(b, avail, src1, src2, index, alpha1, alpha2) =
	let
	 val srcs = [src1, src2]
	 val Ty.TensorTy(ts, NONE) = IR.Var.ty src1
	 val Ty.TensorTy(ts', NONE) = IR.Var.ty src2
	 val params = [E.TEN(true, ts, NONE), E.TEN(true, ts', NONE)]
	 val body = E.Op2(if b then E.Max else E.Min, E.Tensor(0, alpha1), E.Tensor(1, alpha2))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "maxTensor", Ty.TensorTy(ts, NONE), ein, srcs)
	end
    fun makeAbsTenAlpha(avail, src, alpha, index) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy
	 val params = [E.TEN(true, ts, NONE)]
	 val body = E.Op1(E.Abs, E.Tensor(0, alpha))
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "absTensor", tenTy, ein, [src])
	end
	  

    fun vmax(b, avail, srcs) =
	let
	 val tys = List.map (IR.Var.ty) srcs
	 val Ty.TensorTy(ts, NONE) :: tys' = tys

	 val _ = if List.all (fn x => Ty.same(Ty.TensorTy(ts, NONE), x)) tys'
		 then ()
		 else raise Fail "bad tensor max"
	 val _ = if List.length srcs < 1
		 then raise Fail "vmax needs at least 1 args"
		 else ()
	 val fs::fss = srcs
	in
	 List.foldr (fn (x, y) => vmax2(b, avail, x, y)) fs fss
	end

    fun minInterval(avail, arg) = sliceOutEin(avail, 0, arg)
    fun maxInterval(avail, arg) = sliceOutEin(avail, 1, arg)

    fun radius(avail, arg, k) =
	 if k = 1
	 then makeIntervalWidth(avail, arg)
	 else if k = 2
	 then sliceOutEin(avail, 1, arg)
	 else if k > 2
	 then let
	  val nume = k - 2
	  val errArgs = List.tabulate(nume, fn x => makeAbsTenAtIdx(avail, x + 1, arg))
	  val errornArg = sliceOutEin(avail, k - 1, arg)
	  val added = makeAddTen(avail, errArgs@[errornArg])
	 in
	  added
	 end
	 else raise Fail "impossible k"


    fun average(avail, args) =
	let
	 val added = makeAddTen(avail, args)
	 val num = List.length args
	 val one = Rational.fromInt 1
	 val rational = Rational.fromInt num
	 val dived =  (Rational.div(one, rational)) 

	 val Ty.TensorTy(shp, NONE) = IR.Var.ty added
	 val index = shp
	 val params = [E.TEN(true, shp, NONE)]
	 val body = E.Opn(E.Prod, [E.ConstR dived, E.Tensor(0, List.tabulate(List.length shp, fn x => E.V x))])
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "average", Ty.TensorTy(shp, NONE), ein, [added])

	end
    fun scatterPatternVecs(avail, arg) = (*run on ten or radius of interval*)
	let
	 val ty = IR.Var.ty arg
	 val Ty.TensorTy(shp, NONE) = ty
	 val size = (List.foldr (fn (x,y) => x * y) 1 shp)
	 val shapedArray = ArrayNd.getInverseIndex(ArrayNd.computeInverseIndex(ArrayNd.array'(0, shp)))
	 val knockOuts = List.tabulate(size, fn x =>
						List.tabulate(size,
							      fn y => if x = y
								      then
								       AvailRHS.assignOp(avail, "deltaii",
											 Ty.realTy,
											 Op.TensorIndex(ty, Array.sub(shapedArray, x)),
											 [arg])
								      else AvailRHS.makeRealLit(avail,
												[],
												RealLit.zero true)))
	 val knockOutArrays = List.map (fn x => ArrayNd.fromList'(x, shp)) knockOuts
	 fun groupConversion(lst, s) = let
	  val v::vst = lst
	  val ty = IR.Var.ty v
	  val Ty.TensorTy(shape', _) = ty
	 in
	  AvailRHS.assignCons(avail, "shapeSpliceCons", lst, Ty.TensorTy(s::shape', NONE))
	 end
	 val results = List.map (fn x => ArrayNd.convertToTree(x, fn y => y, groupConversion)) knockOutArrays
	in
	 results
	end
    (*intersection, hull; extend, inside - special*)
    (*rest: covered above + cons bs (or deltas maybe)*)
    (*Radius: get t_1, t_2,... -> add, add*)
    (*make vadd, make vsub
      make abs, make vmax
      intervalToAffine shape pattern
      tensorToAffine/co -> cons builder
    all destructors obvious
    Vcomp for inside, extend is cons, hull/intersection need vmax
    
    plan: this today + fix the eval basis thing or the cvtTy thing... - fix evalBasis here!
    plan: ein/basis/vectorized next 2 days
    plan: const lift monday and tuesday
    plan: parse on Wed (simple changes + the n-tree)
    plan: test next days
    
     *)
			

    fun handleOp(avail, y, oper, args, oldArgs, oldY) =
	let
	 fun finish(link) = AvailRHS.addAssignToList(avail, (y, IR.VAR link))
	 fun handleOp' (Op.intervalSimple, [a]) = intervalCons(avail, a, a)
	   | handleOp' (Op.intervalMixed, [a, b]) = intervalCons(avail, a, b)
	   | handleOp' (Op.intervalAffine, [a]) = (*builds interval via affine -> compute rad, center and then sub/add*)
	     let
	      val [a'] = oldArgs
	      val Ty.TensorTy(_, SOME k ) = IR.Var.ty a'
	      val rad = radius(avail, a, k)
	      val center = sliceOutEin(avail, 0, a)
	      val min = makeSubTen(avail, center, rad)
	      val max = makeAddTen(avail, [center, rad])
				  
	     in
	      intervalCons(avail, min, max)
	     end
	   | handleOp' (Op.intervalToAffine, [a]) =
	     let
	      val rad = makeIntervalWidth(avail, a)
	      val scattered = scatterPatternVecs(avail, rad)
	      val center = makeIntervalMidPoint(avail, a)
	      val zeros = intervalZeros(avail, a)
						
	     in
	      affineCons(avail, center, scattered, zeros)
	     end
	   | handleOp' (Op.affineNative(t1, t2, t3), [a, b, c]) =
	     let
	      val (innerTy, count) = (case t2
				       of Ty.SeqTy(t as Ty.TensorTy(r::rs, NONE), SOME d) => if (List.length rs) <> 0
											     then (t, d)
											     else raise Fail "bad error size"
					| Ty.TensorTy(r::ts, NONE) => if (List.length ts) <> 0
								      then (Ty.TensorTy(ts, NONE), r)
								      else raise Fail "bad error size"
					| _ => raise Fail "bad error types"
				     (* end case*))
	      val  _ = if Ty.same(t1, innerTy) andalso Ty.same(innerTy, t3)
		       then ()
		       else raise Fail "bad error types for affine cons"

				  
	      val tenargs = (case t2
			      of Ty.SeqTy(t as Ty.TensorTy(_, NONE), SOME d) =>
				 List.tabulate(count, fn x => AvailRHS.assignOp(avail, "affineSeqAcc", t, Op.Subscript(t2),
										[b, AvailRHS.makeIntLit(avail, x)]))
			       | Ty.TensorTy(t::ts, NONE) => List.tabulate(count, fn x => sliceOutEin(avail, x + 1, b))
			       | _ => raise Fail "impossible"
			    (* end case*))
	     in
	      affineCons(avail, a, tenargs, c)
	     end
	   | handleOp' (Op.affineNative2(t), [a, b]) = affine2Cons(avail, a, b)
	   | handleOp' (Op.errors, [a]) =
	     let
	      val [a'] = oldArgs
	      val Ty.TensorTy(shp, SOME k ) = IR.Var.ty a'
	     in
	      if k < 2
	      then raise Fail "bad error op: k <= 1"
	      else if k = 2 (*d=0*)
	      then AvailRHS.makeNan(avail, shp) (*raise Fail "bad error op: k = 2" (*QUESTION: should there be a shp zero here*)*)
	      else if k = 3 (*d = 1; just a scalar*)
	      then sliceOut(avail, 1, a) (*d >= 2*)
	      else AvailRHS.assignCons(avail, "errorCons", List.tabulate(k - 2,
									 fn x => sliceOut(avail, x + 1, a)),
				       Ty.TensorTy((k-2)::shp, NONE))
	     end
	   | handleOp' (Op.lasterr, [a]) =
	     let
	      val [a'] = oldArgs
	      val Ty.TensorTy(shp, SOME k ) = IR.Var.ty a'
	     in
	      if k < 2
	      then raise Fail "bad errorn op: k <= 1"
	      else sliceOut(avail, k - 1, a)
	     end
	   | handleOp' (Op.radius, [a]) =
	     let
	      val [a'] = oldArgs
	      val Ty.TensorTy(shp, SOME k ) = IR.Var.ty a'
	     in
	      radius(avail, a, k)
	     end
	   | handleOp' (Op.minInterval, [a]) = minInterval(avail, a)
	   | handleOp' (Op.maxInterval, [a]) = maxInterval(avail, a)
	   | handleOp' (Op.intersection, [a, b]) =
	     let
	   (*clampTTT, a < x < b*)
	   (*min_a min_b max_a*)
	   (*min_a max_b max_a*)
	   (*if B \subset_neq A then this will presreve b
	    if A \subset_neq B  then will give A
	    if the intersect A then B so min_b < min_a < max_b < max_a
	   or min_a < min_b < max_a < max_b then

           if they are seperate so either max_a < min_b < max_b or max_b < min_a

	   We now do clamp a_0... < b_i... < a_1... for i=0..1, ....
	    *)
	      val [t1, t2] = List.map (IR.Var.ty) [a, b]
	      val (Ty.TensorTy(t::ts, NONE), Ty.TensorTy(t'::ts', NONE)) = (t1, t2)
	      val _ = if Ty.same(t1, t2)
		      then ()
		      else raise Fail "bad intersection -not same shape"									     
	      val _ = if t=t' andalso t=2
		      then ()
		      else raise Fail "bad intersection-non-interval"

	      val tens = [E.TEN(true, t::ts, NONE), E.TEN(true, t::ts, NONE)]
	      val vs = (List.tabulate(List.length ts, fn x => E.V (x+1)))
	      val clamp = E.Op3(E.Clamp, E.Tensor(0, E.C 0 :: vs), E.Tensor(1, E.V 0 :: vs), E.Tensor(0, E.C 1 :: vs))
	      val ein1 = E.EIN{params=tens, index=t::ts, body=clamp}
	     in
	      AvailRHS.assignEin(avail, "intervalIntersection", Ty.TensorTy(t::ts, NONE), ein1, [a, b])
	     end
	   | handleOp' (Op.hull, [a, b]) =
	     let
	      val mins = [minInterval(avail, a), minInterval(avail, b)]
	      val maxes = [maxInterval(avail, a), maxInterval(avail, b)]
	      val min = vmax(false, avail, mins)
	      val max = vmax(true, avail, maxes)
	     in
	      intervalCons(avail, min, max)
	     end
	   | handleOp' (Op.extend, [a, new]) =
	     let
	      val ty = IR.Var.ty a
	      val Ty.TensorTy(k::shp, NONE) = ty
	      val start = sliceOut(avail, 0, a)
	      val endV = sliceOut(avail, k - 1, a)
              val rest = List.tabulate(k-2, fn x => sliceOut(avail, x + 1, a))
	     in
	      affineCons(avail, start, rest@[new], endV)
	     end
	   | handleOp' (Op.insideInterval, [ten, interval]) = AvailRHS.assignOp(avail, "insided", IR.Var.ty ten, Op.twocomp,
										[minInterval(avail, interval), ten,
										 maxInterval(avail, interval)])
	   | handleOp' (oper as Op.EvaluateBasis(bda), [arg]) =
	     (case IR.Var.ty oldY
	       of Ty.TensorTy(_, NONE) => AvailRHS.assignOp(avail, "plainbasis", IR.Var.ty y, oper, [arg])
		| Ty.TensorTy(_, SOME k) =>
		  let
		  in
		   AvailRHS.assignOp(avail, "intervalbasis", IR.Var.ty y, Op.EvaluateBasisAff(bda, k), [arg])
		  end
	     (* end case*))
	   | handleOp' (oper, vs) = AvailRHS.assignOp(avail, "temp", IR.Var.ty y, oper, vs)
	in
	 finish(handleOp'(oper, args))
	end


    fun verifyEinForm(body) =
	let
	 val (sx, body') = (case body
			     of E.Sum(sx, e) => (sx, e)
			      | _ => ([], body)
			   (* end case *))
	 fun highIr () = raise Fail "high ir term in mid ir"
	 fun verifyLeaf (E.Const _) = true
	   | verifyLeaf (E.ConstR _) = true
	   | verifyLeaf (E.Tensor _) = true
	   | verifyLeaf (E.Zero _) = true
	   | verifyLeaf (E.Delta _) = true
	   | verifyLeaf (E.Epsilon _) = true
	   | verifyLeaf (E.Eps2 _) = true
	   | verifyLeaf (E.Field _) = highIr ()
	   | verifyLeaf (E.Lift _) = highIr ()
	   | verifyLeaf (E.Identity _) = highIr ()
	   | verifyLeaf (E.Fem _) = highIr ()
	   | verifyLeaf (E.Conv _) = highIr ()
	   | verifyLeaf (E.Partial _) = highIr ()
	   | verifyLeaf (E.Apply _) = highIr ()
	   | verifyLeaf (E.Comp _) = highIr ()
	   | verifyLeaf (E.Probe _) = highIr ()
	   | verifyLeaf (E.OField _) = raise Fail "disallowed term"
	   | verifyLeaf (E.Value _) = true
	   | verifyLeaf (E.Img _) = true
	   | verifyLeaf (E.Krn _) = true
	   | verifyLeaf (E.Poly _) = raise Fail "disallowed term"
	   | verifyLeaf (E.If _) = raise Fail "disallowed term"
	   | verifyLeaf (E.Op1 _) = false
	   | verifyLeaf (E.Op2 _) = false
	   | verifyLeaf (E.Op3 _) = false
	   | verifyLeaf (E.Opn _) = false

	 fun leafs (E.Op1(_, b)) = [b]
	   | leafs (E.Op2(_, b1, b2)) = [b1, b2]
	   | leafs (E.Op3(_, b1, b2, b3)) = [b1, b2, b3]
	   | leafs (E.Opn(_, bs)) = bs
	   | leafs _ = raise Fail "impossible"
	in
	 if verifyLeaf body'
	 then (sx, true, body')
	 else if List.all verifyLeaf (leafs body')
	 then (sx, false, body')
	 else raise Fail "badly formed mid ein"
	end

    fun rankOfExp'(e, params) =
	let
	 fun maxi s = List.foldr (Int.max) 0 s
	 fun minAboveZero s = List.foldr (Int.min) 0 (List.filter (fn x => x > 0) s)
	 fun failOrMax(es, rs) =
	     if List.length es <> List.length rs
	     then raise Fail "internal failure in rank of exp"
	     else if (maxi rs) <> (minAboveZero rs)
	     then raise Fail ("internal rank error: typechecker or later created expressions of incompatible interval rank:\n "
			      ^ (String.concat(ListPair.map (fn (x,y) => "rank(" ^ (EinPP.expToString x) ^ ") = " ^ (Int.toString y) ^ "\n") (es, rs))))
	     else maxi rs

	 fun recur e' = rankOfExp'(e', params)
	 fun failOrMax'(es) = failOrMax(es, List.map recur es)

	 fun doTensor(idx) = if 0 > idx orelse idx >= List.length params
			     then raise Fail "bad tensor param"
			     else (case List.nth(params, idx)
				    of E.TEN(_, _, SOME j) => j
				     | E.TEN(_, _, NONE) => 0
				     | _ => raise Fail ("bad Tensor arg: " ^ (Int.toString idx))
				  (* end case *))


				    
	in
	 (case e
	   of E.Const _ => 0
	    | E.ConstR _ => 0
	    | E.Tensor(pid, _) => doTensor pid
	    | E.Zero _ => 0
	    | E.Delta _ => 0
	    | E.Epsilon _ => 0
	    | E.Eps2 _ => 0
	    | E.Field _ => raise Fail "Field in mid"
	    | E.Lift e => recur e
	    | E.Identity _ => raise Fail "id-field in mid"
	    | E.Conv _ => raise Fail "Convo in mid"
	    | E.Fem _ => raise Fail "Fem in mid"
	    | E.Partial _ => raise Fail "Partial in mid"
	    | E.Apply _ => raise Fail "Apply in mid"
	    | E.Comp _ => raise Fail "comp in mid"
	    | E.Probe(e, at) => recur at
	    | E.OField _ => raise Fail "ignore ofield"
	    | E.Value _ => 0
	    | E.Img _ => 0
	    | E.Krn _ => 0
	    | E.Poly _ => raise Fail "ignore Poly"
	    | E.If(_, a, b) => failOrMax'([a,b])
	    | E.Sum(sx, e') => recur e'
	    | E.Op1(u, e') => recur e'
	    | E.Op2(b, e1, e2) => failOrMax'([e1, e2])
	    | E.Op3(t, e1, e2, e3) => failOrMax'([e1, e2, e3])
	    | E.Opn(opn, es) => failOrMax' es
	 (*end case*))
	end

    fun leafPid (E.Tensor(pid, _)) = SOME(pid)
      | leafPid (E.Lift(e)) = leafPid e
      | leafPid (E.Const _) = NONE
      | leafPid (E.ConstR _) = NONE
      | leafPid (E.Zero _) = NONE
      | leafPid (E.Delta _) = NONE
      | leafPid (E.Epsilon _) = NONE
      | leafPid (E.Eps2 _) = NONE
      | leafPid _ = raise Fail "bad leaf"


    fun cvtParams (E.TEN(b, shp, NONE)) = E.TEN(b, shp, NONE)
      | cvtParams (E.TEN(b, shp, SOME k)) = if k = 1
					    then E.TEN(b, 2::shp, NONE)
					    else if k > 1
					    then E.TEN(b, k::shp, NONE)
					    else raise Fail "bad TEN"
      | cvtParams (E.FLD _) = raise Fail "high-ir param"
      | cvtParams (other) = other

    fun detectBadTen r (E.TEN(b, shp, SOME k)) = r = k
      | detectBadTen r _ = true


    fun toTwoVec(avail, arg) =
	let
	 val ty = IR.Var.ty arg
	 val Ty.TensorTy(s::innerShape, NONE) = ty
	 val _ = if s <> 2
		 then raise Fail "bad interval tensor"
		 else ()
	 val temp = ArrayNd.array'(0, innerShape)
	 val inverse = (ArrayNd.getInverseIndex o ArrayNd.computeInverseIndex)(temp)
	 val make2Var = List.map (fn x => let
				   val x' = x
				   val a = AvailRHS.assignOp(avail, "acc02", Ty.realTy, Op.TensorIndex(ty, 0::x'), [arg])
				   val b = AvailRHS.assignOp(avail, "acc12", Ty.realTy, Op.TensorIndex(ty, 1::x'), [arg])
				  in
				   AvailRHS.assignCons(avail, "acc2", [a,b], Ty.TensorTy([2],NONE))
				  end) (Array.toList(inverse))
									      
	in
	 make2Var
	end

    fun buildIntervalFromTwoVec(avail, innerShape, twoVecs) =
	let
	 val temp = ArrayNd.array'(0, innerShape)
	 val inverse = (ArrayNd.getInverseIndex o ArrayNd.computeInverseIndex)(temp)

	 fun build i = List.map (fn x => AvailRHS.assignOp(avail, "unacc0"^(Int.toString i), Ty.realTy, Op.TensorIndex(Ty.TensorTy([2], NONE), [i]), [x])) twoVecs
	 val twovec1 = build 0
	 val twovec2 = build 1
	 val indexInfo = ArrayNd.fromList'(List.@(twovec1, twovec2), 2::innerShape)
	 val ret = ArrayNd.convertToTree(indexInfo, (fn x => x), groupConversionS "unaccgroup" avail)
	in
	 ret
	end


    fun swapIntervalMinMax(avail, interval) =
	let
	 val ty as Ty.TensorTy(s::shp, NONE) = IR.Var.ty interval
	 val _ = if s <> 2
		 then raise Fail "bad interval to swapIntervalMinMax"
		 else ()
	 val temp = ArrayNd.array'(0, s::shp)
	 val inverse = (ArrayNd.getInverseIndex o ArrayNd.computeInverseIndex)(temp)
	 val inverseArray = ArrayNd.fromArray'(inverse, s::shp)
	 fun swapFn(i::idx) =
	     let
	      val i' = if i = 0
		       then 1
		       else if i = 1
		       then 0
		       else raise Fail "impossible"
				  
	     in
	      AvailRHS.assignOp(avail, "intervalSawpAcc", Ty.realTy, Op.TensorIndex(ty, i'::idx), [interval])
	     end
	in
	 ArrayNd.convertToTree(inverseArray, swapFn , groupConversionS "swapInterval" avail)
	end

    fun inc(E.V x) = E.V (x + 1)
      | inc (E.C j) = E.C j
    fun modAlpha alpha = E.V 0 :: (List.map inc alpha)
    fun modTensor(E.Tensor(pid, alpha)) = E.Tensor(0, E.V 0 :: (List.map inc alpha))
    fun modSumrange(xs) = List.map (fn (x,y,z) => (x + 1, y, z)) xs

    fun mkSum([], b) = b
      | mkSum (sx, b) = E.Sum(sx, b)
	  


    fun handleIntervalEin(y, oldY, params, index, sx, body, body', leafForm, args, oldArgs) =
	let
	 val avail = AvailRHS.new ()
	 fun rankOfExp(e) = let val i = rankOfExp'(e, params) in if i = 0
								 then NONE
								 else if i > 0
								 then SOME(i)
								 else raise Fail "bad interval rank" end

	 fun mkEin (index, param, b) = E.EIN{index=index, params=params, body=b}
	 val tests = List.all (detectBadTen 2) params
	 val _ = if tests
		 then ()
		 else raise Fail "bad params to interval ein"
	 val params' = List.map cvtParams params
	 val noInterval = List.all (detectBadTen (~1)) params
	 val overallRank = rankOfExp(body')


	 fun handleOp1 (E.PowInt(j), t as E.Tensor(pid, alpha)) =
	     if Int.mod(j, 2) = 1
	     then
	      AvailRHS.assignEin(avail, "tempIntervalPowOdd", IR.Var.ty y, mkEin(2::index, params',
									      mkSum(modSumrange(sx),
										    E.Op1(E.PowInt(j),
											  modTensor(t)))),
				 [List.nth(args, pid)])
	     else 
	     let
	      (*min(xmax, max(0, xmin))*)
	      (*max(abs(xmin), xmax)*)
	      val arg = List.nth(args, pid)
	      val E.TEN(_, argShape, _) = List.nth(params, pid)
	      val (argShape', alpha') = alphaIndex(argShape, alpha)

	      val xmin = minInterval(avail, arg)
	      val xmax = minInterval(avail, arg)
	      val minabs = makeAbsTenAlpha(avail, xmin, alpha, argShape')
	      val xmax' =  vmax2alpha(true, avail, minabs, xmax, argShape', alpha', alpha)

	      val zeros = AvailRHS.makeRealLit(avail, argShape', RealLit.zero true)
	      val maxzeromin = vmax2alpha(true, avail, zeros, xmin, argShape', alpha', alpha)
	      val xmin' = vmax2alpha(false, avail, xmax, maxzeromin, argShape', alpha', alpha')

	      val result = intervalCons(avail, xmin', xmax')				    
				       
	     in
	      AvailRHS.assignEin(avail, "tempIntervalPowEven", IR.Var.ty y, mkEin(2::index,
										  params',
										  mkSum(modSumrange(sx),
											(E.Tensor(0, modAlpha alpha')))), [result])
				
	     end
	   | handleOp1 (E.Abs, t as E.Tensor(pid, alpha)) =
	     let
	      
	     (*
max(abs(xmin), xmax)
abs(min(xmax, max(0, xmin))) 
	      *)

	      val arg = List.nth(args, pid)
	      val E.TEN(_, argShape, _) = List.nth(params, pid)
	      val (argShape', alpha') = alphaIndex(argShape, alpha)
	      val xmin = minInterval(avail, arg)
	      val xmax = minInterval(avail, arg)
	      val minabs = makeAbsTenAlpha(avail, xmin, alpha, argShape') 
	      val xmax' = vmax2alpha(true, avail, minabs, xmax, argShape', alpha', alpha) (*alpha' doesn't have any more constants*)

	      val zeros = AvailRHS.makeRealLit(avail, argShape', RealLit.zero true)
	      val zeromax= vmax2alpha(true, avail, zeros, xmin, argShape, alpha', alpha)
	      val minzeromax= vmax2alpha(false, avail, xmax, zeromax, argShape', alpha', alpha')
	      val xmin' = makeAbsTen(avail, minzeromax)
	      val result = intervalCons(avail, xmin', xmax')
	     in
	      AvailRHS.assignEin(avail, "tempIntervalAbs", IR.Var.ty y, mkEin(2::index,
									      params',
									      mkSum(modSumrange(sx),
										    (E.Tensor(0, modAlpha alpha')))), [result])
	     end
	   | handleOp1 (operator : E.unary, t as E.Tensor(pid, alpha)) =
	     let
	      val dec = (false, true)
	      val inc = (true, false)
	      val nei = (false, false)
	      val (increasing, decreasing) = (case operator
					       of E.Neg => dec
						| E.Exp => inc
						| E.Sqrt => inc
						| E.Cosine => nei
						| E.ArcCosine => nei
						| E.Sine => nei
						| E.ArcSine => nei
						| E.Tangent => nei
						| E.ArcTangent => nei
						| E.PowInt(j) => raise Fail "handle otherway"
						| E.Abs => raise Fail "handle otherway"
						| E.Sgn => raise Fail "no sgn of interval"
					     (* end case *))
	     in
	      if increasing
	      then AvailRHS.assignEin(avail, "tempIntervalInc", IR.Var.ty y, mkEin(2::index, params',
										   mkSum(modSumrange(sx),
											 E.Op1(operator, modTensor(t)))), [List.nth(args, pid)])
	      else if decreasing
	      then let
	       val arg = List.nth(args, pid)
	       val swaped = swapIntervalMinMax(avail, arg)
	      in
	       AvailRHS.assignEin(avail, "tempIntervalDec", IR.Var.ty y, mkEin(2::index,
									       params',
									       mkSum(modSumrange(sx),
										     E.Op1(operator, modTensor(t)))),
				  [swaped])
	      end
	      else
	       let
		val arg = List.nth(args, pid)
		val twoVecs = toTwoVec(avail, arg)
		val name = (case operator
			     of E.Cosine => "cosine"
			      | E.ArcCosine => "arccosine"
			      | E.Sine => "sine"
			      | E.ArcSine => "arcsine"
			      | E.Tangent => "tangent"
			      | E.ArcTangent => "arctangent"
			   (*end case*))
		val twoVecs' = List.map (fn x => AvailRHS.assignOp(avail, "twovecfunc", Ty.TensorTy([2], NONE),
								   Op.scalarIntervalFun(name), [x])) twoVecs
		val result = buildIntervalFromTwoVec(avail, index, twoVecs')
						    
	       in
		AvailRHS.assignEin(avail, "tempIntervalDec", IR.Var.ty y, mkEin(2::index,
										params',
										mkSum(modSumrange(sx),
										      modTensor(t))), [result])
	       end
	     end


	 (**)
	 (*doubleConst or doubleTensor or ...*)
	 (*(sx or []), duplicate leaf form -> \sum ..
	   op1 - increasing we do obvious - and insert in sum
           op1 - decreasing we seperate, swap, same
	   op1 - otherwise, we two vector -- how?
           ---split into two -  using inverse map -> (i, ijks)-> (a_0ijk, a_1ijk)
           ---apply to 2op to get a_i
           ---build back using tree func
           ---substitute in place.
	   ---DO this part in about 80 mins ish


	 
	   op2 - sub is a swap one and then do substitute into sub
	   op2 - div - nancheck and then simple product chase -- need to check for scalar stuff - swap stuff...??

	   opn - add is clear - split or no?
	   opn - prod - complciated - lots of avail

           JUST DO IT --affine is similar for one op probably and div... rest is pretty clean
	  *)
	in
	 (case overallRank
	   of NONE => (*pure scalar operation that we just need to duplicate the index of 2::rest and fix the indexes*)
	      let
	       val index' = 2::index
	       val body'' = EU.mapIndex(body, fn x => x + 1) 
	       val ein = E.EIN{params=params, index=index', body=body''}
	      in
	       if noInterval
	       then [IR.ASSGN(y, IR.EINAPP(ein, args))]
	       else raise Fail "impossible"
	      end
	    | SOME _ => (case body'
			  of E.Op1(unary, t) => (AvailRHS.addAssignToList(avail, (y, IR.VAR (handleOp1(unary, t))));
						 AvailRHS.getAssignments' avail)
			   | _ => raise Fail "NYI"
			(* end case *))
	 (* end case *))
	end

    fun handleEin(y, oldY, ein as E.EIN{params, index, body}, args, oldArgs) =
	let
	 val (sx, leafForm, body') = verifyEinForm(body)
	 val Ty.TensorTy(shp, iv) = IR.Var.ty oldY
	in
	 (case iv
	   of NONE => [IR.ASSGN(y, IR.EINAPP(ein, args))]
	    | SOME k => if k = 1
			then handleIntervalEin(y, oldY, params, index, sx, body, body', leafForm, args, oldArgs)
			else raise Fail "NYI"
	 (* end case*))
	end




    fun expand (env, (y, rhs)) = let
          fun assign rhs = [IR.ASSGN(Env.rename (env, y), rhs)]
          in
            case rhs
             of IR.GLOBAL x => assign (IR.GLOBAL(Env.renameGV(env, x)))
              | IR.STATE(NONE, fld) => assign (IR.STATE(NONE, Env.renameSV(env, fld)))
              | IR.STATE(SOME x, fld) =>
                  assign (IR.STATE(SOME(Env.rename(env, x)), Env.renameSV(env, fld)))
              | IR.VAR x => assign (IR.VAR(Env.rename(env, x)))
              | IR.LIT lit => assign (IR.LIT lit)
              | IR.OP(rator, args) =>
		let
		 val y' = Env.rename(env, y)
		 val args' = Env.renameList(env, args)
		 val avail = AvailRHS.new()
		 val _ = handleOp(avail, y', rator, args', args, y)
		in
		 (AvailRHS.getAssignments'(avail))
		end
	      (*Affine ops used -- all else pass on (BASIS or pass into low?)*)
                  (* List.map IR.ASSGN (expandOp (env, Env.rename (env, y), rator, args)) *)
              | IR.CONS(args, ty) => assign (IR.CONS(Env.renameList(env, args), cvtTy ty))
              | IR.SEQ(args, ty) => assign (IR.SEQ(Env.renameList(env, args), cvtTy ty))
              | IR.EINAPP(rator, args) => raise Fail "oops"
              | IR.APPLY(f, args) =>
                  assign (IR.APPLY(Env.renameFV(env, f), Env.renameList(env, args)))
              | _ => raise Fail("bogus rhs for ASSIGN: " ^ IR.RHS.toString rhs)
            (* end case *)
    end



  (* expand a IR multi-assignment to a IR CFG *)
    fun mexpand (env, (ys, rhs)) = let
          fun massign rhs = let
                val nd = IR.Node.mkMASSIGN(Env.renameList(env, ys), rhs)
                in
                  IR.CFG{entry=nd, exit=nd}
                end
          fun mkOP (rator, xs) = massign(IR.OP(rator, Env.renameList(env, xs)))
          in
            case rhs
             of IR.OP(Op.EigenVecs2x2, xs) => mkOP (Op.EigenVecs2x2, xs)
              | IR.OP(Op.EigenVecs3x3, xs) => mkOP (Op.EigenVecs3x3, xs)
              | IR.OP(Op.KillAll, []) => mkOP (Op.KillAll, [])
              | IR.OP(Op.StabilizeAll, []) => mkOP (Op.StabilizeAll, [])
              | IR.OP(Op.Print tys, xs) => mkOP (Op.Print(List.map cvtTy tys), xs)
	      | IR.OP(Op.Save args, xs) => mkOP(Op.Save args, xs)
              | IR.MAPREDUCE mrs => let
                  val mrs = List.map
                        (fn (r, f, xs) => (r, Env.renameFV(env, f), Env.renameList(env, xs)))
                          mrs
                  in
                    massign (IR.MAPREDUCE mrs)
                  end
              | _ => raise Fail("bogus rhs for MASSIGN: " ^ IR.RHS.toString rhs)
            (* end case *)
          end

    structure Trans =  TranslateFn (
      struct
        open Env
        val expand = IR.CFG.mkBlock o expand
        val mexpand = mexpand
      end)

    structure Promote = PromoteFn (IR)

    fun transform prog =
	let
	 val prog = Trans.translate prog
	in
	 MidCensus.init prog;
	 Promote.transform prog
	end

  end
