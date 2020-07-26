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
      abs -> min(max(min, 0), max), max(abs(max), abs(min)) -- vectorized
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
		 List.map (fn x => IR.ASSGN x) (AvailRHS.getAssignments(avail))
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
