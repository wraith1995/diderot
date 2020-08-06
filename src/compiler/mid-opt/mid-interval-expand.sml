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
    structure IMap = IntRedBlackMap
    (*given a tensor shape and an alpha, compute a new shape and alpha to represent the result after alpha is used in another op*)
    (*
\sum ijk...

    *)




    (**index, sx -> alpha_ijk..
       sum index, index index, constant index
       --constants are eliminated
       --sums are preserved 
       --index are permuted
\sum gamma E.Tensor(gamma..constants..alpha) 

--Remove constants
alpha, gamma -> E.Tensor(gamma..alpha)
\sum_krs -A_jkirs
B_ijkrs = -A_jkirs
\sum B_ijk
3 3 E.V 1 E.V 2 E.V 0
3 3 E.V 0 E.V 1 E.V 2
3 3 ... ... 

.*)


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

    fun sliceOutEinAlpha(avail, alpha, src) =
	let
	 val tenTy = IR.Var.ty src
	 val Ty.TensorTy(ts, NONE) = tenTy

	 val params = [E.TEN(true, ts, NONE)]
	 val index = ts
	 val body = E.Tensor(0, alpha)
	 val ein = E.EIN{params=params, index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "intervalInterConsEinAlpha", Ty.TensorTy(index, NONE), ein, [src])
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

    fun make2Prod(avail, arg1, arg2, index, alpha1, alpha2) =
	let
	 val Ty.TensorTy(ts1, NONE) = IR.Var.ty arg1
	 val Ty.TensorTy(ts2, NONE) = IR.Var.ty arg2
	 val body = E.Opn(E.Prod, [E.Tensor(0, alpha1), E.Tensor(1, alpha2)])
	 val ein = E.EIN{params=[E.TEN(true, ts1, NONE), E.TEN(true, ts2, NONE)], index=index, body=body}
	in
	 AvailRHS.assignEin(avail, "tensor2Mult", Ty.TensorTy(index, NONE), ein, [arg1, arg2])
	end
    fun makeVProd(avail, arg1, arg2)  =
	let
	 val Ty.TensorTy([dim], NONE) = IR.Var.ty arg1
	 val Ty.TensorTy([dim'], NONE) = IR.Var.ty arg1
	 val _ = if dim = dim'
		 then ()
		 else raise Fail "impossible: unequal dims in Vprod"
	 val alpha = [E.V 0]
	 val index = [dim]
	in
	 make2Prod(avail, arg1, arg2, index, alpha, alpha)
	end

    fun makeVScale(avail, scale, arg) =
	let
	 val _ = (case IR.Var.ty scale
		   of Ty.TensorTy([], NONE) => ()
		    | Ty.TensorTy([dim], NONE) =>  if dim = 1 orelse dim = 0
						   then ()
						   else raise Fail "can't scale by vector"
		    | _ => raise Fail "can't scale by vector"
		 (* end case *))
	 val Ty.TensorTy([dim], NONE) = IR.Var.ty arg
	 val params = [E.TEN(true, [], NONE), E.TEN(true, [dim], NONE)]
	 val args = [scale, arg]
	 val body = E.Opn(E.Prod, [E.Tensor(0, []), E.Tensor(0, [E.V 0])])
	 val ein = E.EIN{params=params, index=[dim], body=body}
	in
	 AvailRHS.assignEin(avail, "scale", Ty.TensorTy([dim], NONE), ein, args)
	end

    fun makeVIndex(avail, arg, idx) =
	let
	 val ty = IR.Var.ty arg
	 val Ty.TensorTy([dim], NONE) = ty
	 val _ = if 0 <= idx andalso idx < dim
		 then ()
		 else raise Fail ("bad idx in VIndex")
	 val opper = Op.TensorIndex(ty, [idx])
	in
	 AvailRHS.assignOp(avail, "acc", Ty.realTy, opper, [arg])
	end
    fun getVElements(avail, arg) =
	let
	 val ty = IR.Var.ty arg
	 val Ty.TensorTy([dim], NONE) = ty
	in
	 List.tabulate(dim, fn x => makeVIndex(avail, arg, x))
	end
		    
    fun makeCons(avail, args) =
	let
	 val len = List.length args
	 val _ = if len >= 2
		 then ()
		 else raise Fail "bad cons len"
	 val a::_ = args
	 val ty = IR.Var.ty a
	 val Ty.TensorTy(shp, NONE) = ty
	 val newShp = len::shp
	 val _ = if List.all (fn x => Ty.same(x, ty)) (List.map IR.Var.ty args)
		 then ()
		 else raise Fail "bad tys for cons"
	in
	 AvailRHS.assignCons(avail, "cons", args, Ty.TensorTy(newShp, NONE))
	end
    fun flatConsVectors(avail, args : IR.var list) =
	let
	 fun seperateThem(num, stuff, r::args : IR.var list) = let val vs = getVElements(avail, r)
						 in seperateThem(num+List.length vs, List.@(stuff, vs), args)
						 end
	   | seperateThem (num, stuff, []) = (num, stuff)

	 val (num, stuff) = seperateThem(0, [], args)
	in
	 AvailRHS.assignCons(avail, "consSplatter", stuff, Ty.TensorTy([num], NONE))
	end

    fun makeTensorProdFactorized(avail, args) =
	let
	 val tys = List.map (IR.Var.ty) args
	 val _ = if List.length args < 2
		 then raise Fail "bad tensor product"
		 else ()

	 val t::tys' = tys
	 val Ty.TensorTy([dim], NONE) = t
	 val _ = List.all (fn x => Ty.same(t, x)) tys'

	 fun recurse (old, []) = old
	   | recurse (old, r::rest) =
	     let
	      val seperated = List.tabulate(dim, fn x => makeVIndex(avail, r, x))
	      val scaled = List.map (fn x => makeVScale(avail, x, old)) seperated
	      val cons = flatConsVectors(avail, scaled)
	     in
	      recurse(cons, rest)
	     end
	 val a::args' = args

			  
	in
	 recurse(a, args)
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

    fun vminmax(b, avail, srcs, index, alpha) =
	let
	 fun foldf(x, y) = vmax2alpha(b, avail, x, y, index, alpha, alpha)
	 fun dumbFoldr f (x::y::xs) = dumbFoldr f (f(x,y)::xs)
	   | dumbFoldr f [x] = x
	   | dumbFoldr f [] = raise Fail "bad vminmax"
	in
	 dumbFoldr foldf srcs
	end

    fun vminmaxV(b, avail, srcs, index) =
	let
	 val alpha = List.tabulate(List.length(index), fn x => E.V x)
	 fun foldf(x, y) = vmax2alpha(b, avail, x, y, index, alpha, alpha)
	 fun dumbFoldr f (x::y::xs) = dumbFoldr f (f(x,y)::xs)
	   | dumbFoldr f [x] = x
	   | dumbFoldr f [] = raise Fail "bad vminmax"
	in
	 dumbFoldr foldf srcs
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

    fun normalizeEin(sx, index, body)=
	let
	 (*for all sx, we need to make them be index::len(sx)*)
	 val rank = List.length index
	 val sx' = List.map (fn (x,y,z) => x) sx
	 val sxRank = List.length sx
	 val rankPairs = (List.tabulate(sxRank, fn x => (List.nth(sx', x), x + rank)))
	 val rankMap = (List.foldr (fn ((x,y), a) => IMap.insert(a, x, y)) IMap.empty rankPairs)

	 val newSx = List.map (fn (x,y,z) => (Option.valOf(IMap.find(rankMap, x)), y, z)) sx
	 fun fixMu(E.V x) = if x < rank
			     then E.V x
			     else (case IMap.find (rankMap, x)
				    of SOME x' => E.V x'
				     | NONE => raise Fail "bad rank map"
				  (* end case *))
	   | fixMu (E.C x) = E.C x
	 fun fixAlpha(alpha) = List.map fixMu alpha
	 fun fixLeaf (E.Tensor(id, alpha)) = E.Tensor(id, fixAlpha alpha)
	   | fixLeaf (E.Zero(alpha)) = E.Zero(fixAlpha alpha)
	   | fixLeaf (E.Delta(mu1, mu2)) = E.Delta(fixMu mu1, fixMu mu2)
	   | fixLeaf (E.Epsilon(mu1, mu2, mu3)) = E.Epsilon(fixMu mu1, fixMu mu2, fixMu mu3)
	   | fixLeaf (E.Eps2(mu1, mu2)) = E.Eps2(fixMu mu1, fixMu mu2)
	   | fixLeaf (E.Img(pid, alpha, pos, i)) = E.Img(pid, fixAlpha alpha, pos, i)
	   | fixLeaf (E.Krn(pid, mus, j)) = E.Krn(pid, List.map (fn (x,y) => (fixMu x, fixMu y)) mus, j)
	   | fixLeaf (E.Poly _) = raise Fail "don't use Poly"

	 fun fixOp(E.Op1(u, e1)) = E.Op1(u, fixLeaf e1)
	   | fixOp (E.Op2(b, e1, e2)) = E.Op2(b, fixLeaf e1, fixLeaf e2)
	   | fixOp (E.Op3(tr, e1, e2, e3)) = E.Op3(tr, e1, e2, e3)
	   | fixOp (E.Opn(opn, es)) = E.Opn(opn, List.map fixLeaf es)
	in
	 if sxRank = 0
	 then (sx, body)
	 else (newSx, fixOp body)
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
      | leafPid (E.Img(pid, _, _, _)) = SOME(pid)
      | leafPid (E.Krn(pid, _, _)) = SOME(pid)
      | leafPid (E.ConstR _) = NONE
      | leafPid (E.Zero _) = NONE
      | leafPid (E.Delta _) = NONE
      | leafPid (E.Epsilon _) = NONE
      | leafPid (E.Eps2 _) = NONE
      | leafPid _ = raise Fail "bad leaf"

    fun mapPid f e =
	let
	 fun leafPid (E.Tensor(pid, a)) = E.Tensor(f pid, a)
	   | leafPid (E.Img(pid, a, b, c)) = E.Img(f pid, a, b, c)
	   | leafPid (E.Krn(pid, a, b)) = E.Krn(f pid, a, b)
	   | leafPid e = e
	in
	 leafPid e 
	end



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
    fun incAlpha alpha = List.map inc alpha
    fun modTensor(E.Tensor(pid, alpha)) = E.Tensor(0, E.V 0 :: (List.map inc alpha))
    fun modSumrange(xs) = List.map (fn (x,y,z) => (x + 1, y, z)) xs
    fun makeIvMask (E.TEN(_, _ ,SOME _)::ts) = true::makeIvMask(ts)
      | makeIvMask (_::ts) = false::makeIvMask(ts)
      | makeIvMask ([]) = []
    fun modLeaf ivmask (e as E.Const _) = e
      | modLeaf ivmask (e as E.Tensor(pid, alpha)) = E.Tensor(pid, if List.nth(ivmask, pid)
								 then modAlpha alpha
								 else incAlpha alpha)
      | modLeaf ivmask (E.Zero alpha) = E.Zero (incAlpha alpha)
      | modLeaf ivmask (E.Delta(mu1, mu2)) = E.Delta(inc mu1, inc mu2)
      | modLeaf ivmask (E.Epsilon(mu1, mu2, mu3)) = E.Epsilon(inc mu1, inc mu2, inc mu3)
      | modLeaf ivmask (E.Eps2(mu1, mu2)) = E.Eps2(inc mu1, inc mu2)
      | modLeaf ivmask _ = raise Fail "bad leaf"
    fun mkSum([], b) = b
      | mkSum (sx, b) = E.Sum(sx, b)

    fun argToTen arg = (case IR.Var.ty arg
			 of Ty.TensorTy(shp, j) => E.TEN(true, shp, j)
			  | _ => raise Fail "bad param")

    fun mkEin (index, params, b) = E.EIN{index=index, params=params, body=b}



    fun collectAlpha (E.Tensor(_, alpha)) = alpha
      | collectAlpha (E.Zero alpha) = alpha
      | collectAlpha (E.Delta(mu1, mu2)) = [mu1, mu2]
      | collectAlpha (E.Epsilon(mu1, mu2, mu3)) = [mu1, mu2, mu3]
      | collectAlpha (E.Eps2(mu1, mu2)) = [mu1, mu2]
      | collectAlpha (E.Img(_, alpha, _, _)) = alpha
      | collectAlpha _ = []

    fun collectAlphas ops = List.map collectAlpha ops

    fun alphaReplaceError(a, b) r = if List.length a = List.length b
				    then r
				    else raise Fail "bad alpha replace"
    fun uncollectAlpha (alpha, E.Tensor(pid, alpha')) =
	(alphaReplaceError(alpha, alpha') (E.Tensor(pid, alpha)))
      | uncollectAlpha (alpha, E.Zero(alpha')) = alphaReplaceError(alpha, alpha') (E.Zero(alpha))
      | uncollectAlpha ([mu1, mu2], E.Delta(_)) = E.Delta(mu1, mu2)
      | uncollectAlpha ([mu1, mu2, mu3], E.Epsilon _) = E.Epsilon(mu1, mu2, mu3)
      | uncollectAlpha ([mu1, mu2], E.Eps2 _) = E.Eps2(mu1, mu2)
      | uncollectAlpha (alpha, E.Img(a, alpha', b, c)) = (alphaReplaceError(alpha, alpha') (E.Img(a, alpha, b, c)))
      | uncollectAlpha ([], e) = e
      | uncollectAlpha _ = raise Fail "uncollect"
    fun uncollectAlphas (alphas, opss) = ListPair.map uncollectAlpha (alphas, opss)

    fun extractOpts(sx, index, opts) =
	let
	 (* build map from E.V x to size *)
	 val indexRank = List.length index
	 val sxRank = List.length sx
	 fun sizeFn x = if x < indexRank andalso x >= 0
			then List.nth(index, x)
			else if x < sxRank
			then let val (_, x, y) = List.nth(sx, x)
			     in 1 + (y - x) end (*size is inclusize...*)
			else raise Fail ("bad size index: " ^ (Int.toString x))

	 (* Find all relevant E.V xs used in our opts*)
	 val alphas = collectAlphas opts
	 val alphaIdx = List.mapPartial (fn E.V x => SOME(x) | E.C _ => NONE) (List.concat alphas)

	 (*Order E.V x in the order in which they are used in the original expression and remove duplicates*)
	 fun isolate [] = []
	   | isolate (x::xs) = x::isolate(List.filter (fn y => y <> x) xs)
	 fun sort data = ListMergeSort.sort (fn (x,y) => x > y) data
	 val alphaIdxSort = isolate (sort alphaIdx)
	 (* Calculate new index via the used E.V x*)
	 val index' = List.map sizeFn alphaIdxSort

	 (* build map from the old E.V x to the new E.V x (i.e with respect to index and index') and the inverse back*)
	 val alphaPairs = (List.tabulate(List.length alphaIdxSort, fn x => (List.nth(alphaIdxSort, x), x)))
	 val alphaIdxMap = (List.foldr (fn ((x,y), a) => IMap.insert(a, x, y)) IMap.empty alphaPairs) (*oldIdx -> idx of idx*)
	 val alphaIdxMapInv = (List.foldr (fn ((y,x), a) => IMap.insert(a, x, y)) IMap.empty alphaPairs) (* idx of idx -> oldIdx*)
	 fun mapReplace mapp idx = (case IMap.find(mapp, idx)
				     of NONE => raise Fail "bad map in split"
				      | SOME k => k
				   (* end case*))
	 fun mapReplace' mapp idx = (case idx
				      of E.V k => E.V(mapReplace mapp k)
				       | E.C k => E.C k
				    (* end case*))

	 (* Using the map from old E.V x to new E.V x, updated all used alphas to be used in the new *)
	 val alphas' = List.map (List.map (mapReplace' alphaIdxMap)) alphas
	 (* Update the opts accordingly *)
	 val opts' = uncollectAlphas(alphas', opts)
	 (*Using the inverse mod, build an alpha to access the resulting tensor in the old expression; no constants, obviously*)
	 val replaceAlpha = List.map (fn x => E.V x) (List.map (mapReplace alphaIdxMapInv) (List.tabulate(List.length index', fn x => x)))
	in
 	(*return opts with renamed alphas to operate within index' and 
	  then be indexed in the old expression via replaceAlpha*)
	 (opts, index', replaceAlpha)
	end

    fun splitOpt(sx, index, opts1, opts2) =
	(extractOpts(sx, index, opts1), extractOpts(sx, index, opts2), extractOpts(sx, index, opts1@opts2))

    fun partitionScalarInterval(es, params) = List.partition (fn x => rankOfExp'(x, params) = 0) es

    fun binarySequences(n) =
	let
	 fun bs(r) =
	     if r = 1
	     then [[0], [1]]
	     else let val many = bs(r - 1)
		  in List.concatMap (fn x => [0::x, 1::x]) many
		  end

	in
	 if n <=0
	 then raise Fail "bad bin seq"
	 else bs(n)
	end

    fun tensorsAndSequence(bnseq, tensors) =
	let
	 fun append(i::is, E.Tensor(pid, alpha)::ts) = E.Tensor(pid, E.C i :: alpha):: append(is, ts)
	in
	 append(bnseq, tensors)
	end
    fun tensorsAndSequences(bnseqs, tensors) = List.map (fn x => tensorsAndSequence(x, tensors)) bnseqs

    fun makeProducts(avail, args, sx, index, tensors) =
	let
	 fun run(E.Tensor(pid1, alpha1)::E.Tensor(pid2, alpha2)::rest, NONE) =
	     let
	      val arg1 = List.nth(args, pid1)
	      val arg2 = List.nth(args, pid2)
	      val ([E.Tensor(_, alpha1'), E.Tensor(_, alpha2')], index', rAlpha) =
		  extractOpts(sx, index, [E.Tensor(pid1, alpha1), E.Tensor(pid2, alpha2)])
	      
	      val ret = make2Prod(avail, arg1, arg2, index', alpha1', alpha2')
	     in
	      run(E.Tensor(pid2, rAlpha)::rest, SOME(ret))
	     end
	   | run (E.Tensor(pid1, alpha1)::E.Tensor(pid2, alpha2)::rest, SOME arg1) =
	     let
	      val arg2 = List.nth(args, pid2)
	      val (([E.Tensor(_, alpha1'), E.Tensor(_, alpha2')], index', rAlpha)) =
		  extractOpts(sx, index, [E.Tensor(pid1, alpha1), E.Tensor(pid2, alpha2)])
	      val ret = make2Prod(avail, arg1, arg2, index', alpha1', alpha2')
	     in
	      run(E.Tensor(pid2, rAlpha)::rest, SOME(ret))
	     end
	   | run (E.Tensor(pid, alpha1)::[], SOME(ret)) = ret
	   | run (E.Tensor(pid, alpha1)::[], NONE) = (*only possible if we start with one; basically just does selection*)
	     let
	      val arg = List.nth(args, pid)
	      val Ty.TensorTy(shp, NONE) = IR.Var.ty arg
	      val (([E.Tensor(_, alpha1')], index', rAlpha)) =
		  extractOpts(sx, index, [E.Tensor(pid, alpha1)])
	      val ein = E.EIN{params=[E.TEN(true, shp, NONE)], index=index', body=E.Tensor(0, alpha1')}
	     in
	      AvailRHS.assignEin(avail, "prodTwo1", Ty.TensorTy(index', NONE), ein, [arg])
	     end
	in
	 run(tensors, NONE)
	end
	  

    fun handleOpn (avail, y, sx, index, params, args) (opn, ts) =
	let
	 fun handleOpn' (E.Add, ts) =
	     let
	      val params' = List.map cvtParams params
	      val mask = makeIvMask params
	      val ts' = List.map (modLeaf mask) ts
	      val index' = 2::index
	      val opn = mkSum(modSumrange(sx), E.Opn(E.Add, ts'))
	      val ein = E.EIN{params=params', index=index', body=opn}
	     in
	      AvailRHS.assignEin(avail, "intervalAdd", Ty.TensorTy(index', NONE), ein, args)
	     end
	   | handleOpn' (E.Prod,ts) =
	     let
	      val (scalar, interval) = partitionScalarInterval(ts, params)
	      val (scalarT, intervalT, combinedT) = splitOpt(sx, index, scalar, interval)
	      val numInterval = List.length interval
	      val numScalar = List.length scalar
	      val binSeqs = binarySequences(numInterval)
	      val (intervalTOpts, intervalIndex, intervalAlpha) = intervalT
	      val intervalTOptsMin = tensorsAndSequences(binSeqs, intervalTOpts)
	      val (prodVars, prodAlpha, prodIndex) =
		   let
		    val prods = List.map (fn x => makeProducts(avail, args, sx, intervalIndex, x)) intervalTOptsMin
		   in
		    (prods, intervalAlpha, intervalIndex)
		   end

	      val (preMin, preMax) =
		  if numInterval = 1
		  then (List.nth(prodVars, 0), List.nth(prodVars, 1)) (*min and max already*)
		  else let
		   val min = vminmaxV(false, avail, prodVars, prodIndex)
		   val max = vminmaxV(true, avail, prodVars, prodIndex)
		  in
		   (min, max)
		  end
	      val (_, combinedIndex, combinedAlpha) = combinedT (* combined alpha is the alpha with respect to the original index/sum*)
	      val (sMin, sMax, sIndex, sAlpha) =
		  if numScalar = 0
		  then (preMin, preMax, combinedIndex, combinedAlpha)
		  else
		   let
		    (*gather and fix scalar ones*)
		    val (scalarTOpts, scalarIndex, scalarAlpha) = scalarT
		    fun rebuildParams(E.Tensor(pid, alpha)::rest, cnt, new, params, args) =
			let val param = List.nth(params, pid)
			    val arg = List.nth(args, pid)
			in
			 rebuildParams(rest, cnt + 1, E.Tensor(cnt, alpha)::new, param::params, arg::args)
			end
		      | rebuildParams (E.Img(pid, alpha, pos, j)::rest, cnt, new, params, args) =
			let val param = List.nth(params, pid)
			    val arg = List.nth(args, pid)
			in
			 rebuildParams(rest, cnt + 1, E.Img(cnt, alpha, pos, j)::new, param::params, arg::args)
			end
		      | rebuildParams (E.Krn(pid, mus, j)::rest, cnt, new, params, args) =
			let val param = List.nth(params, pid)
			    val arg = List.nth(args, pid)
			in
			 rebuildParams(rest, cnt + 1, E.Krn(cnt, mus, j)::new, param::params, arg::args)
			end
		      | rebuildParams ([], _, rest, params, args) =
			(List.rev rest, List.rev params, List.rev args)
		    val (newScalarTopts, newParams, newArgs) = rebuildParams(scalarTOpts, 1, [], [E.TEN(true, intervalIndex, NONE)], [])
		    val body = E.Opn(E.Prod, scalarTOpts@[E.Tensor(0, prodAlpha)])
		    val ein = E.EIN{
			 params = newParams,
			 index = combinedIndex,
			 body = body
			}
		    val min' = AvailRHS.assignEin(avail, "minWithScalar", Ty.TensorTy(combinedIndex, NONE), ein, preMin::newArgs)
		    val max' = AvailRHS.assignEin(avail, "maxWithScalar", Ty.TensorTy(combinedIndex, NONE), ein, preMax::newArgs)
		    (*FIXME: we could add optimizations here based on what the constants are...*)
		    val min'' = vmax2(false, avail, min', max')
		    val max'' = vmax2(true, avail, min', max')
		   in
		    (min'', max'', combinedIndex, combinedAlpha)
		   end
	      val einSum = mkSum(sx, E.Tensor(0, sAlpha))
	      val einParam = [E.TEN(true, combinedIndex, NONE)]
	      val einIndex = index (*combinedIndex will be part original index and part sx*)
	      val ein = E.EIN{params=einParam, index=einIndex, body=einSum}
	      val sum1 = AvailRHS.assignEin(avail, "sumMin", Ty.TensorTy(index, NONE), ein, [sMin])
	      val sum2 = AvailRHS.assignEin(avail, "sumMax", Ty.TensorTy(index, NONE), ein, [sMax])
	     in
	      intervalCons(avail, sum1, sum2)
	     end
	in
	 handleOpn' (opn, ts)
	end
    fun alphaAccessM(b, avail, arg, alpha, sx, index) =
	let
	 val Ty.TensorTy(argshp, NONE) = IR.Var.ty arg
	 val opts = [E.Tensor(0, E.C (if b then 1 else 0)::alpha)]
	 val ([E.Tensor(_, alpha')], index', replaceAlpha) = extractOpts(sx, index, opts)
	 val ein = Ein.EIN{params=[E.TEN(true, argshp, NONE)], index=index', body=E.Tensor(0, alpha')}
	in
	 (AvailRHS.assignEin(avail, "alphaAccM", Ty.TensorTy(index', NONE), ein, [arg]), index', replaceAlpha, sx)
	end

    fun makeExtractedDiv(avail, sx, index, opt, params, args) =
	let
	 val params' = List.map cvtParams params
	 val ([opt'], newIndex, replaceAlpha) = extractOpts(sx, index, [opt])
	 val oneLit = AvailRHS.makeRealLit(avail, newIndex, RealLit.one)
	 val default = List.tabulate(List.length newIndex, fn x => E.V x)
	 (**addAssign(t, "lit", IR.Ty.realTy, IR.LIT(Literal.Real (RealLit.one)))*)
	 val possiblePid = leafPid opt'
	 val (newArgs, opt'') = 
	     (case possiblePid
	       of NONE => ([oneLit], mapPid (fn x => 1) opt')
		| SOME pid => ([oneLit, List.nth(args, pid)], mapPid (fn x => 1) opt')
	     (* end case *))
	 val divEin = E.Op2(E.Div, E.Tensor(0, default), opt'')
	 val divParams = (case possiblePid
			   of NONE => [E.TEN(true, newIndex, NONE)]
			    | SOME pid => [E.TEN(true, newIndex, NONE), List.nth(params', pid)])
	 val ein = E.EIN{params=divParams, index=newIndex, body=divEin}
	 val ret = AvailRHS.assignEin(avail, "inverseEin", Ty.TensorTy(newIndex, NONE), ein, newArgs)
	in
	 (ret, newIndex, replaceAlpha)
	end
	  
    fun handleOp1s (avail, y, sx, index, params, args) (op1, t) =
	let
	 fun handleOp1 (E.PowInt(j), t as E.Tensor(pid, alpha)) =
	     if Int.mod(j, 2) = 1
	     then
	      AvailRHS.assignEin(avail, "tempIntervalPowOdd", IR.Var.ty y, mkEin(2::index, [cvtParams (List.nth(params, pid))],
										 mkSum(modSumrange(sx),
										       E.Op1(E.PowInt(j),
											     modTensor(t)))),
				 [List.nth(args, pid)])
	     else 
	      let
	       (*min(xmax, max(0, xmin))*)
	       (*max(abs(xmin), xmax)*)
	       val arg = List.nth(args, pid)
	       (* val xmin = minInterval(avail, arg) *)
	       (* val xmax = minInterval(avail, arg) *)
	       val (xminAcc, index', alpha', sx') = alphaAccessM(false, avail, arg, alpha, sx, index) (*index and alpha to subsittute in if you extractouside of the sum and opt*)
	       val newRank = (List.length index')
	       val simple = List.tabulate(newRank, fn x => E.V x) (*just do pairwise operations*)
	       val (xmaxAcc, _, _, _) = alphaAccessM(true, avail, arg, alpha, sx, index)
	       val minabs = makeAbsTenAlpha(avail, xminAcc, simple, index')
	       val xmax' =  vmax2alpha(true, avail, minabs, xmaxAcc, index', simple, simple)

	       val zeros = AvailRHS.makeRealLit(avail, index', RealLit.zero true)
	       val maxzeromin = vmax2alpha(true, avail, zeros, xminAcc, index', simple, simple)
	       val xmin' = vmax2alpha(false, avail, xmaxAcc, maxzeromin, index', simple, simple)
	       val result = intervalCons(avail, xmin', xmax')				    
	      in
	       AvailRHS.assignEin(avail, "tempIntervalPowEven", IR.Var.ty y, mkEin(2::index, [E.TEN(true, index', NONE)],
										   mkSum(modSumrange(sx'),
											 E.Op1(E.PowInt(j),
											       modTensor(E.Tensor(0, alpha'))))),
				  [result])				 
	      end
	   | handleOp1 (E.Abs, t as E.Tensor(pid, alpha)) =
	     let
	      
	      (*
max(abs(xmin), xmax)
abs(min(xmax, max(0, xmin))) 
	      *)

	      val arg = List.nth(args, pid)
	      val (xminAcc, index', alpha', sx') = alphaAccessM(false, avail, arg, alpha, sx, index)
	      val (xmaxAcc, _, _, _) = alphaAccessM(true, avail, arg, alpha, sx, index)

	      val newRank = (List.length index')
	      val simple = List.tabulate(newRank, fn x => E.V x) (*just do pairwise operations*)
	  
						  
	      val minabs = makeAbsTenAlpha(avail, xminAcc, simple, index') 
	      val xmax' = vmax2alpha(true, avail, minabs, xmaxAcc, index', simple, simple) 

	      val zeros = AvailRHS.makeRealLit(avail, index', RealLit.zero true)
	      val zeromax= vmax2alpha(true, avail, zeros, xminAcc, index', simple, simple)
	      val minzeromax= vmax2alpha(false, avail, xmaxAcc, zeromax, index', simple, simple)
	      val xmin' = makeAbsTen(avail, minzeromax)

	      val result = intervalCons(avail, xmin', xmax')
	     in
	      AvailRHS.assignEin(avail, "tempIntervalAbs", IR.Var.ty y, mkEin(2::index,
									      [argToTen result],
									      mkSum(sx',
										    modTensor(E.Tensor(0, alpha')))), [result])
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
	      then AvailRHS.assignEin(avail, "tempIntervalInc", IR.Var.ty y, mkEin(2::index, [argToTen (List.nth(args, pid))],
										   mkSum(modSumrange(sx),
											 E.Op1(operator, modTensor(t)))), [List.nth(args, pid)])
	      else if decreasing
	      then let
	       val arg = List.nth(args, pid)
	       val swaped = swapIntervalMinMax(avail, arg)
	      in
	       AvailRHS.assignEin(avail, "tempIntervalDec", IR.Var.ty y, mkEin(2::index,
									       [argToTen swaped],
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
										[argToTen result],
										mkSum(modSumrange(sx),
										      modTensor(t))), [result])
	       end
	     end

	in
	 handleOp1 (op1, t)
	end

    fun handleOp2s (avail, y, sx, index, params, args) (op1, t1, t2) =
	let
	 fun handleOp2 (E.Sub) = (*if second is an interval, swap it; then do as in add*)
	     let
	      val (t1p, pid1, t2id) = (case t1
						  of E.Tensor(pid, alpha) => (E.Tensor(0, alpha), SOME(pid), 1)
						   | E.Img(pid, alpha, pos, j) => (E.Img(0, alpha, pos, j), SOME(pid), 1)
						   | E.Krn(pid, mus, j) => (E.Krn(0, mus, j), SOME(pid), 1)
						   | _ => (t1, NONE, 0)
						(* end case*))
	      val (t2p, pid2, t2arg) =
		  (case t2
		    of E.Tensor(pid, alpha) =>
		       (case List.nth(params, pid)
			 of E.TEN(_, _, SOME 1) =>
			    let
			     val arg = List.nth(args, pid)
			     val swaped = swapIntervalMinMax(avail, arg)
			    in
			     (E.Tensor(t2id, alpha), SOME(pid), SOME(swaped))
			    end
			  | E.TEN(_, _, NONE) =>  (E.Tensor(t2id, alpha), SOME(pid),
						   SOME(List.nth(args, pid)))
		       (* end case*))
		     | E.Img(pid, alpha, pos, j) => (E.Img(t2id, alpha, pos, j), SOME(pid),
						     SOME(List.nth(args, pid)))
		     | E.Krn(pid, mus, j) => (E.Krn(t2id, mus, j), SOME(pid),
					      SOME(List.nth(args, pid)))
		     | _ => (t2, NONE, NONE)
		  (* end case *))

	      val oldPids = List.map (Option.valOf) (List.filter (Option.isSome) [pid1, pid2])
	      val params = List.map (fn x => List.nth(params, x)) oldPids
	      val mask = makeIvMask params
	      val newParams = List.map cvtParams params
	      val newLeafs = List.map (modLeaf mask) ([t1p, t2p])

	      val newArgs = (case (pid1, t2arg)
			      of (SOME j, SOME arg) => [List.nth(args, j), arg]
			       | (NONE, SOME arg) => [arg]
			       | (SOME j, NONE) => [List.nth(args, j)]
			       | (NONE, NONE) => raise Fail "no tensor argument for E.Sub interval is impossible!"
			    (*end case*))
	      val [e1, e2] = newLeafs

	      val newSub = E.Op2(E.Sub, e1, e2)
	      val body = mkSum(modSumrange(sx), newSub)
	      val ein = E.EIN{params = newParams, index=2::index, body = body}
	      val ret = AvailRHS.assignEin(avail, "intervalSub", Ty.TensorTy(2::index, NONE), ein, newArgs)
	     in
	      ret
	     end
	   | handleOp2 (E.Div) =
	     let
	      (*This code:
		  make a/b into a_ialpha/b_ibeta with i the interval var (depending on prescence of interval)
		  extract a, index, ialpha
		  extract 1/b, index, ibeta - this makes ibeta work
		  
		  for 1/b, if b is an interval, account for infentities as needed and min/max

		  access a so ialpha works

		  rebuild provide alpha as alpha and beta as beta, build corresponding tensor expression, build old style params.
		  This allows one to basically build the old ein with converted args ( the accesed 1/a and the 1/b accessed with some swaping)
		  This can be passed to the handleOpn function with the prod argument.
	       *)
	      val mask = makeIvMask params
	      val (t1', t2') = (modLeaf mask t1, modLeaf mask t2)
	      val sx' = modSumrange(sx)
	      val index' = 2 :: index
	      val ([t1''], numIndex, numAlpha) = extractOpts(sx', index', [t1'])
	      val (divVar, divIndex, divAlpha) = makeExtractedDiv(avail, sx', index', t2', params, args)
	      val bottomArg =
		  (case t2
		    of E.Tensor(pid, alpha) =>
		       (case List.nth(params, pid)
			 of E.TEN(_, shp, SOME 1) =>
			    let 
				val zero = AvailRHS.makeRealLit(avail, shp, RealLit.zero true)
				val posInf = AvailRHS.makeRealLit(avail, shp, RealLit.posInf)
				val negInf = AvailRHS.makeRealLit(avail, shp, RealLit.negInf)
				val min = sliceOutEin(avail, 0, divVar)
				val max = sliceOutEin(avail, 1, divVar)

				fun opj j = Op.zerotestselect j
				val op1 = AvailRHS.assignOp(avail, "minnan", Ty.TensorTy(divIndex, NONE), opj 0, [min, max, zero, negInf])
				val op2 = AvailRHS.assignOp(avail, "minnan", Ty.TensorTy(divIndex, NONE), opj 1, [min, max, zero, posInf])
				val min' = vmax2(false, avail, op1, op2)
				val max' = vmax2(true, avail, op1, op2)
				val consed = intervalCons(avail, min', max')
			    in
			     SOME(consed)
			    end
			  | E.TEN(_, _, NONE) => (NONE) (*just a scalar*)
			    
			  | E.TEN(_, _, _) => raise Fail "E.Div: / by affine with d > 1"
		       (* end case *))
		     | _ => NONE
		  (* end case *))

	      val divVar' = Option.getOpt(bottomArg, divVar)

	      val (numVar, topInterval) =
		  let
		   val pid = leafPid t1'
		   val t1''' = mapPid (fn x => 0) t1''
		   val (params, args, interval) = (case pid
						    of NONE => ([], [], false)
						     | SOME pid' =>
						       (case List.nth(params, pid')
							 of E.TEN(_, _, NONE) =>
							    ([cvtParams (List.nth(params, pid'))],
							     [List.nth(args, pid')], false)
							  | E.TEN(_, _, SOME 1) =>
							    ([cvtParams (List.nth(params, pid'))],
							     [List.nth(args, pid')], true)
							  | E.TEN(_, _, _) => raise Fail "affine div"
							  | _ => ([cvtParams (List.nth(params, pid'))],
								  [List.nth(args, pid')], false)
						       (* end case *))
					(* end case *))
		   val ein = E.EIN{params=params, index=numIndex, body=t1'''}
		  in
		   (AvailRHS.assignEin(avail, "numAcc", Ty.TensorTy(numIndex, NONE), ein, args), interval)
		  end

	      fun revertParamInterval(v) =
		  let
		   val Ty.TensorTy(s::shp, NONE) = IR.Var.ty v
		   val _ = if s = 2
			   then ()
			   else raise Fail "bad revertParamInterval"
		  in
		   E.TEN(true, shp, SOME 1)
		  end
	      fun grabScalarParam(v) =
		  let
		   val Ty.TensorTy(shp, NONE) = IR.Var.ty v
		  in
		   E.TEN(true, shp, NONE)
		  end
	      fun decAlpha alpha = List.map (fn E.V x => E.V (x - 1) | E.C j => E.C j) alpha
	      fun dropAlpha alpha = decAlpha (List.drop(alpha, 1))

	      val bottomInterval = Option.isSome bottomArg
	      val (args1, args2) =
		  (case (bottomInterval, topInterval)
		    of (false, false) => raise Fail "scalar div in interval"
		     | (true, true) => ((sx, index, [revertParamInterval numVar, revertParamInterval divVar'], [numVar, divVar']),
					(E.Prod, [E.Tensor(0, dropAlpha numAlpha), E.Tensor(1, dropAlpha divAlpha)]))
		     | (true, false) => ((sx, index, [revertParamInterval numVar, grabScalarParam divVar'], [numVar, divVar']),
					 (E.Prod, [E.Tensor(0, dropAlpha numAlpha), E.Tensor(1, decAlpha divAlpha)]))
		     | (false, true) => ((sx, index, [grabScalarParam numVar, revertParamInterval divVar'], [numVar, divVar']),
					 (E.Prod, [E.Tensor(0,  decAlpha numAlpha), E.Tensor(1, dropAlpha divAlpha)]))
		  (* end case *))
	      val args1' = (avail, y, #1 args1, #2 args1, #3 args1, #4 args1)
	     in
	      handleOpn args1' args2
	     end
				  
	   | handleOp2 (E.Max) = raise Fail "bad: no interval max allowed"
	   | handleOp2 (E.Min) = raise Fail "bad: no interval min allowed"
	in
	 handleOp2 op1
	end

    fun handleIntervalEin(y, oldY, params, index, sx, body, body', leafForm, args, oldArgs) =
	let
	 val avail = AvailRHS.new ()
	 fun rankOfExp(e) = let val i = rankOfExp'(e, params) in if i = 0
								 then NONE
								 else if i > 0
								 then SOME(i)
								 else raise Fail "bad interval rank" end


	 val tests = List.all (detectBadTen 2) params 
	 val _ = if tests
		 then ()
		 else raise Fail "bad params to interval ein"
	 val noInterval = List.all (detectBadTen (~1)) params
	 val overallRank = rankOfExp(body')




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
			  of E.Op1(unary, t) =>
			     let
			      val v = (handleOp1s (avail, y, sx, index, params, args) (unary, t))
			     in
			     (AvailRHS.addAssignToList(avail, (y, IR.VAR v));
			      AvailRHS.getAssignments' avail)
			     end
			   | E.Op2(binary, t1, t2) =>
			     let
			      val v = (handleOp2s (avail, y, sx, index, params, args) (binary, t1, t2))
			     in
			      (AvailRHS.addAssignToList(avail, (y, IR.VAR v));
			       AvailRHS.getAssignments' avail)
			     end
			   | E.Opn(opn, ts) =>
			     let
			      val v = handleOpn (avail, y, sx, index, params, args) (opn, ts)
			     in
			      (AvailRHS.addAssignToList(avail, (y, IR.VAR v));
			       AvailRHS.getAssignments' avail)
			     end
			   | _ => raise Fail ("invalid interval in:" ^ (EinPP.toString(E.EIN{params=params,index=index,body=body})))
			(* end case *))
	 (* end case *))
	end

    fun handleEin(y, oldY, ein as E.EIN{params, index, body}, args, oldArgs) =
	let
	 val (sx, leafForm, body') = verifyEinForm(body)
	 val (sx', body'') = normalizeEin(sx, index, body') (*ensure that summation indecies are from [len(index), len(index) + num)*)
	 val Ty.TensorTy(shp, iv) = IR.Var.ty oldY
	in
	 (case iv
	   of NONE => [IR.ASSGN(y, IR.EINAPP(ein, args))]
	    | SOME k => if k = 1
			then handleIntervalEin(y, oldY, params, index, sx', body, body'', leafForm, args, oldArgs)
			else raise Fail "NYI"
	 (* end case*))
	end
    structure BDA = BasisDataArray
    fun intervalBasisExpand(avail, shp, bda, arg, BDA.Unknown) =
	let
	 val (basisArrayData, meta) = BDA.explode bda
	 (*for the even, we need:*)
	 (*min(xmax, max(0, xmin))*)
	 (*max(abs(xmin), xmax)*)
	 val Ty.TensorTy([dim], NONE) = IR.Var.ty arg
	 val oddmin = sliceOutEin(avail, 0, arg)
	 val oddmax = sliceOutEin(avail, 1, arg)
	 val zeros = AvailRHS.makeRealLit(avail, [dim], RealLit.zero true)
	 val ones = AvailRHS.makeRealLit(avail, [dim], RealLit.one)

	 val zeromax = vmax2(true, avail, zeros, oddmin)
	 val evenmin = vmax2(false, avail, oddmax, zeromax)
			    
	 val absmin = makeAbsTen(avail, oddmin)
	 val evenmax = vmax2(true, avail, absmin, oddmax)

	 val dim' = BDA.domainDim bda
	 val minDegree = BDA.minDegree meta
	 val maxDegree = BDA.maxDegree meta
	 (*(0,1) -> (0,..., degree) -> (0,.., dim -1)*)
	 fun getPowers(degree, oddstart, evenstart) =
	     if degree = 0
	     then (ones, ones)::getPowers(degree + 1, oddstart, evenstart)
	     else if Int.mod(degree, 2) = 1 andalso degree <= maxDegree
	     then
	      let
	       val (omi, omx) = oddstart
	       val (newoddmin, newoddmax) = (makeVProd(avail, oddmin, makeVProd(avail, oddmin, omi)),
					     makeVProd(avail, oddmax, makeVProd(avail, oddmax, omx)))
	      in
	       (newoddmin, newoddmax) :: getPowers(degree + 1, (newoddmin, newoddmax), evenstart)
	      end
	     else if Int.mod(degree, 2) = 0 andalso degree <= maxDegree
	     then
	      let
	       val (emi, emx) = evenstart
	       val (newevenmin, newevenmax) = (makeVProd(avail, evenmin, makeVProd(avail, evenmin, emi)),
					       makeVProd(avail, evenmax, makeVProd(avail, evenmax, emx)))

	      in
	       (newevenmin, newevenmax) :: getPowers(degree + 1, oddstart, (newevenmin, newevenmax))
	      end
	     else []

	 val pows : (IR.var * IR.var) list = getPowers(0, (oddmin, oddmax), (evenmin, evenmax))
	 val (mins, maxes) = ListPair.unzip pows (*mins: degree 0 (x,y,z), degree 1 (x,y,z), ...*)
	 (*min: By var by degree list (x0 x1,...), (y0,...), ...*)
	 fun nth(v, x) = makeVIndex(avail, v, x)
	 val minVectors = List.foldr (fn (x, bydim) => List.tabulate(dim, fn z=> nth(x, z)::List.nth(bydim, z))) (List.tabulate(dim, fn _ => [])) mins
	 val maxVectors = List.foldr (fn (x, bydim) => List.tabulate(dim, fn z => nth(x, z)::List.nth(bydim, z))) (List.tabulate(dim, fn _ => [])) maxes (*size is dim; list of degree vars*)
	 val minVectorVars = List.map (fn x => AvailRHS.assignCons(avail, "pows", x, Ty.TensorTy([List.length x], NONE))) minVectors
	 val maxVectorVars = List.map (fn x => AvailRHS.assignCons(avail, "pows", x, Ty.TensorTy([List.length x], NONE))) maxVectors

				      (*for binary sequences, do tensor products:*)
	 val binseq = binarySequences(dim)
	 val acc = List.tabulate(dim, fn x => x)
	 fun pairVec x = if x = 0
			 then minVectorVars
			 else if x = 1 then maxVectorVars
			 else raise Fail "bad"
	 fun pairMap(x,y) = List.nth(pairVec x, y)
	 val tensorProductGroupings = List.map (fn x => ListPair.map pairMap (x, acc)) binseq

	 val tensorProducts = List.map (fn  x => makeTensorProdFactorized(avail, x)) tensorProductGroupings
	 val size = let
	  val (e::es) = tensorProducts
	  val Ty.TensorTy([dim], NONE) = IR.Var.ty e
	 in
	  dim
	 end

	 (*min/max ^*)
	 val tensorProductMin = vminmaxV(false, avail, tensorProducts, [size])
	 val tensorProductMax = vminmaxV(true, avail, tensorProducts, [size])

	 fun nthswap(b, idx) =
	     let
	      val (m, ma) = if b (*if swap*)
			    then (tensorProductMax, tensorProductMin)
			    else (tensorProductMin, tensorProductMax)
			       
	     in
	      (AvailRHS.assignOp(avail, "swapx", Ty.realTy, Op.TensorIndex(Ty.TensorTy([size], NONE), [idx]), [m]),
	       AvailRHS.assignOp(avail, "swapy", Ty.realTy, Op.TensorIndex(Ty.TensorTy([size], NONE), [idx]), [ma]))
	     end

	 fun swapperVar(reallist) =
	     let
	      val swappers = List.map (fn x => Real.<(x,0.0)) reallist (*swap if realit < 0*)
	      val skippers = List.map (fn x => Real.==(x,0.0)) reallist (*skip if realit = 0 *)
	      fun build((swap, skip)::rest, r::rest', c) =
		  if skip
		  then build(rest, rest', c + 1)
		  else (nthswap(swap, c), r) :: build(rest, rest', c + 1)

	      val correctOrder = build(ListPair.zip (swappers, skippers), reallist, 0)
	      val (minlist, maxlist, reallitlist) = List.foldr (fn (((a,b), c), (l1, l2, l3)) => (a::l1, b::l2, c::l3)) ([], [], []) correctOrder
	      val reallitvars = AvailRHS.makeRealLits(avail, reallitlist)
	     in
	      (flatConsVectors(avail, minlist), flatConsVectors(avail, maxlist), flatConsVectors(avail, reallitvars))
	     end
	 fun convertMin(basis) = (*FIXME NOT MIN/MAX*)
	     let
	      val monoArray = BasisData.explode basis
	      val monoList = ArrayNd.toList monoArray
	      val (minVec, maxVec, realVec) = swapperVar(monoList)
	      val vprodmin = makeVProd(avail, minVec, realVec)
	     in
	      vprodmin
	     end
	 fun convertMax(basis) = (*FIXME NOT MIN/MAX*)
	     let
	      val monoArray = BasisData.explode basis
	      val monoList = ArrayNd.toList monoArray
	      val (minVec, maxVec, realVec) = swapperVar(monoList)
	      val vprodmax = makeVProd(avail, maxVec, realVec)
	     in
	      vprodmax
	     end
	 val minCons = ArrayNd.convertToTree(basisArrayData, convertMin, (groupConversionS "consBasisMin" avail))
	 val maxCons = ArrayNd.convertToTree(basisArrayData, convertMax, (groupConversionS "consBasisMax" avail))
	in
	 intervalCons(avail, minCons, maxCons)
	end


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
		| Ty.TensorTy(shp, SOME 1) => intervalBasisExpand(avail, shp, bda, arg, BDA.analyzeBasis bda)
		| Ty.TensorTy(_, SOME k) =>
		  let
		   val () =()
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
		 (AvailRHS.getAssignments'(avail))
		end
	      (*Affine ops used -- all else pass on (BASIS or pass into low?)*)
                  (* List.map IR.ASSGN (expandOp (env, Env.rename (env, y), rator, args)) *)
              | IR.CONS(args, ty) => assign (IR.CONS(Env.renameList(env, args), cvtTy ty))
              | IR.SEQ(args, ty) => assign (IR.SEQ(Env.renameList(env, args), cvtTy ty))
              | IR.EINAPP(rator, args) =>
		let
		 val newY = Env.rename(env, y)
		 val newArgs = Env.renameList(env, args)
		in
		 handleEin(newY, y, rator, newArgs, args)
		end
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
