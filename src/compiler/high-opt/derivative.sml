(* derivative.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *
 *)

structure Derivative : sig

    val mkApply : Ein.ein_exp * Ein.ein_exp * Ein.param_kind list * Ein.index_id -> Ein.ein_exp option

  end  = struct

    structure E = Ein

    fun err str=raise Fail (String.concat["Ill-formed EIN Operator: ", str])

  (* Helper functions used to make EIN expressions *)
    fun mkAdd exps = E.Opn(E.Add, exps)
    fun mkSub (e1, e2) = E.Op2(E.Sub, e1, e2)
    fun mkProd exps = E.Opn(E.Prod, exps)
    fun mkDiv (e1, e2) = E.Op2(E.Div, e1, e2)
    fun mkNeg e = E.Op1(E.Neg, e)

    val zero = E.Const 0
    val one = E.Const 1

  (* Filters list of expressions in a product *)
    fun filterProd(args:Ein.ein_exp list) = (case EinFilter.mkProd args
           of SOME e => e
            | NONE => mkProd args
          (* end case *))

    fun rewriteProd [a] = a
      | rewriteProd exps = E.Opn(E.Prod, exps)

  (* Apply derivative to unary operatos *)
   fun applyOp1 (op1:Ein.unary, e1:Ein.ein_exp, dx:Ein.alpha)  = let
        val d0::dn = dx
        val px = E.Partial dx
        val inner = E.Apply(E.Partial[d0], e1)
        val square = mkProd [e1, e1]
        val e2 = mkDiv(one, E.Op1(E.Sqrt, mkSub(one, square)))
        fun iterDn e2 = if null dn then e2 else E.Apply(E.Partial dn, e2)
        in
          case op1
            (* distribute derivative *)
           of E.Neg => mkNeg(E.Apply(px, e1))
            (* derivative of exponent function *)
            | E.Exp => iterDn (mkProd [inner, E.Op1(E.Exp, e1)])
            (* chain rule *)
            | E.Sqrt => let
                val half = mkDiv (E.Const 1, E.Const 2)
                val e3 = mkDiv (inner, E.Op1(op1, e1))
                in
                  case dn
                   of [] => mkProd [half, e3]
                    | _ => mkProd [half, E.Apply(E.Partial dn, e3)]
                  (* end case *)
                end
            (* trigonometric identities *)
            | E.Cosine => iterDn (mkProd [mkNeg (E.Op1(E.Sine, e1)), inner])
            | E.ArcCosine => iterDn (mkProd [mkNeg e2, inner])
            | E.Sine => iterDn (mkProd [E.Op1(E.Cosine, e1), inner])
            | E.ArcSine => iterDn (mkProd [e2, inner])
            | E.Tangent =>
                iterDn (mkProd [mkDiv(one, mkProd[E.Op1(E.Cosine, e1), E.Op1(E.Cosine, e1)]), inner])
            | E.ArcTangent =>
                iterDn (mkProd [mkDiv(one, mkAdd[one, square]), inner])
            (* Power rule *)
            | E.PowInt n => iterDn (mkProd [E.Const n, E.Op1(E.PowInt(n-1), e1), inner])
            (* Absolute function *)
            | E.Abs => iterDn  (mkProd [inner, E.Op1(E.Sgn,e1)])
            (* Sign *)
            | E.Sgn => zero
          (* end case *)
        end

  (* Apply derivative to binary operators *)
    fun applyOp2 (op2:Ein.binary, e1:Ein.ein_exp, e2:Ein.ein_exp, dx:Ein.alpha) = let
          val (d0::dn) = dx
          val p0 = E.Partial [d0]
          (* Apply a single derivative to each term *)
          val inner1 = E.Apply(p0, e1)
          val inner2 = E.Apply(p0, e2)
          fun iterDn e2 = if null dn then e2 else E.Apply(E.Partial dn, e2)
          in
            case op2
             (* Distribute derivative over subtraction *)
             of E.Sub => mkSub (inner1, E.Apply(p0, e2))
              | E.Div => (case (e1, e2)
                   of (_, E.Const e2) => mkDiv (inner1, E.Const e2)
                    | _ => (case EinFilter.partitionField [e2]
                      of (_, []) => (case e1
                          of E.Const _ => zero        (* no fields *)
                           | _ => mkDiv (inner1, e2)  (* push derivative to numerator *)
                         (* end case *))
                       | (exp_T, exp_F) => let
                          val exp_F' = E.Apply(p0, rewriteProd exp_F)
                        (* get numerator*)
                          val num = (case e1
                                 of E.Const 1 => mkProd [E.Const ~1, exp_F']
                                  | E.Const c => mkNeg (mkProd [E.Const c, exp_F'])
                                  | _ => mkSub (mkProd (inner1 :: exp_F), mkProd[e1, exp_F'])
                                (* end case *))
                        (* denominator *)
                          val denom = mkProd (exp_T @ exp_F @ exp_F)
                        (* Quotient rule *)
                          in
                            iterDn (mkDiv (num, denom))
                          end
                     (* end case *))
                  (* end case *))
             | E.Max =>  let
                val comp = E.Compare(E.GT, e1, e2)
                val exp  = E.If(comp, inner1, inner2)
                in iterDn exp end
             | E.Min =>  let
                val comp = E.Compare(E.LT, e1, e2)
                val exp  = E.If(comp, inner1, inner2)
                in iterDn exp end
            (* end case *)
          end

  (* Apply derivative to addition and product *)
    fun applyOpn (E.Add:E.opn, es:Ein.ein_exp list, dx:Ein.alpha)  = mkAdd (List.map (fn a => E.Apply(E.Partial dx, a)) es)
      | applyOpn (E.Prod, es, d0::dn) = (
        (* filter out tensor and field terms *)
          case (EinFilter.partitionField es)
           of (_, []) => E.Const 0        (* no fields in expression *)
            | (exp_T, exp_F) => let
              (* Implements product rule for a single derivatie p0 and list of expressions *)
                fun prodAppPartial (p0,[e1]) = E.Apply(p0, e1)
                  | prodAppPartial (p0,e::es) = let
                     val l = prodAppPartial(p0,es)
                     val e2' = filterProd [e, l]
                     val e1' = filterProd (es @ [E.Apply(p0, e)])
                     in
                       mkAdd[e1', e2']
                     end
              (* Apply product rule to field terms then add tensor terms back in *)
                val expDeriv = filterProd (exp_T @ [prodAppPartial (E.Partial [d0], exp_F)])
              (* wrap remaining derivatives around the result *)
                val expFinal = if null dn then expDeriv else E.Apply(E.Partial dn, expDeriv)
                in
                  expFinal
                end
        (* end case*))

    (* rewrite free index (i) in  expression e to vk *)
    (* The idea of this function is that when I take the derivative of a composition e1 \circ e2, a new index, vk, is added     *)
    (* to the expresssion e1 as a result of taking a derivative of it. This index needs to be coordinated with outermost index of e2 for the sum to work *)
    (* Thus, this function finds the outermost index, E.V 0, and replaces all instances of it with vk *)
    (* We again note that implict in this working is that composition is restricted to e2 being a scalar or vector field*)

    fun rewriteIx(vk, e) = let
        (* the variable index inside a composition will start at 0*)
        fun change(E.V 0) = E.V vk
          | change ix     = ix
        in (case e
            of E.Const _ => e
            | E.Tensor(id2, alpha) => E.Tensor(id2, List.map (fn e => change(e)) alpha)
            | E.Zero (alpha) => E.Zero(List.map (fn e => change(e)) alpha)
            | E.Delta (i, j) => E.Delta(change i, change j)
            | E.Epsilon (i, j, k) => E.Epsilon(change i, change j, change k)
            | E.Eps2 (i, j) => E.Eps2(change i, change j)
            | E.Field(id2, alpha) => E.Field(id2, List.map (fn e => change(e)) alpha)
            | E.Lift e1  => E.Lift(rewriteIx(vk, e1))
	    (*I've corected this rule from previous code: *)
	    (* The previous version only tackled the alpha index *)
	    (* This version searches all indicies because it is possible this could be a gradient of a vector field*)
	    (* or a scalar field formed as say the trace of a jacobian: Sum_{i} F_{i}_dx_{i}, which would require us to rewrite both i-s! *)
	    (* The actual bug that made me see this was the case of a gradient of a scalar field: Dx_{i}, where the previous version ignored the i*)
	    | E.Conv(id2, alpha, h2, dxes) => E.Conv(id2, List.map (fn e => change(e)) alpha, h2, List.map (fn e => change(e)) dxes)
	    | E.Fem(fem, id1, id2, id3, alpha, dxes) => E.Fem(fem, id1, id2, id3, List.map (fn e => change(e)) alpha, List.map (fn e => change(e)) dxes)
            | E.Partial(alpha) => E.Partial (List.map (fn e => change(e)) alpha)
	    (* FIX ME:  I'm not sure about these 3: Comeback to check later.*)
	    (* I guess the main issue is: are the indicies of e1,e2 disjoint in these cases?*)
	    (* I think the answer is yes and so this is okay, but I'm not sure*)
	    (* Overlap in indicies here, as opposed to the convo case, would indicate an outter sum*)
	    (* meaning the index was not used in the shape and so not relevant to rewriteIx*)
            | E.Apply(e1, e2) => E.Apply(e1, rewriteIx(vk, e2)) (*implicit here is that e1 is not a E.Partial*)
            | E.Probe(e1, e2) => E.Probe(rewriteIx(vk, e1), e2)
            | E.Comp(e1, es) => E.Comp(rewriteIx(vk, e1), es)
	    (* These are fine:*)
            | E.Sum(sx, e1) => E.Sum(sx, rewriteIx(vk, e1))
            | E.Op1(op1, e1) => E.Op1(op1, rewriteIx(vk, e1))
            | E.Op2(op2, e1, e2) => E.Op2(op2, rewriteIx(vk, e1), rewriteIx(vk, e2))
            | E.Opn(opn, es) => E.Opn(opn, List.map (fn e => rewriteIx(vk, e)) es)
            |  _ => e
            (*end case*))
        end

    fun findDim(e, params) = (case e
            of E.Const _ => NONE
            | E.Tensor _ => NONE
            | E.Zero _ => NONE
            | E.Delta _ => NONE
            | E.Epsilon _ => NONE
            | E.Eps2 _ => NONE
            | E.Lift _ => NONE
            | E.Field(id, _) => let
             val E.IMG(d, _) = List.nth(params, id) 
            in SOME(d) end
				
            | E.Conv(id, _, _, _) => let
                val E.IMG(d, _) = List.nth(params, id)
            in SOME(d) end
	    | E.Fem(_,_, fem, _, _, _) =>
	      let
	       val E.FEM(data) = List.nth(params, fem)


	       val d = FemData.underlyingDim data
	      in
	       SOME(d)
	      end
            | E.Partial _ => NONE
            | E.Apply(e1, e2) => findDim(e2, params)
            | E.Probe(e1, e2) => findDim(e1, params)
            | E.Comp(e1, _) => findDim(e1, params)
            | E.Sum(_, e1) => findDim(e1, params)
            | E.Op1(_, e1) => findDim(e1, params)
            | E.Op2(_, e1, e2) => (case findDim(e1, params)
                of NONE => findDim(e2, params)
                |  e => e
              (*end case *))
            | E.Opn(_, es) => let
                fun iter([]) = NONE
                  | iter(e1::es) = (case findDim(e1, params)
                      of NONE => iter(es)
                      |  e => e
                    (*end case *))
                in iter(es) end
            |  _ => raise Fail(String.concat["find Dim does not handle: ", EinPP.expToString(e)])
		    
							 
        (* end case*))
        
  (* rewrite Apply nodes *)
    fun mkApply (E.Partial dx, e:Ein.ein_exp, params, sumX) = (case e
           of E.Const _ => SOME zero
            | E.ConstR _ => SOME zero
            | E.Tensor _ => err "Tensor without Lift"
            | E.Zero _ => err "Zero without Lift"
            | E.Delta _ => SOME zero
            | E.Epsilon _ => SOME zero
            | E.Eps2 _ => SOME zero
            | E.Field _ => NONE
            | E.Lift _ => SOME zero
            | E.Conv(v, alpha, h, d2) => SOME(E.Conv(v, alpha, h, d2@dx))
	    | E.Fem((E.Plain(bda, n, f)), id1, id2, id3, alpha, dxes)  => SOME(E.Fem(E.Plain( BasisDataArray.D(bda), n, f), id1, id2, id3, alpha, dxes@dx))
	    | E.Fem((E.Invert(bda, n, f)), id1, id2, id3, [mu1], [])  =>
	      let
	       (*QUESTION: What to do about 1D? THIS CODE fucks shit up*)

	       fun fact 0 = 1
		 | fact 1 = 1
		 | fact n = n * ( fact (n-1))
				  
	       val mu2 = List.hd dx
	       val dx' = List.tl dx
	       fun probeAtIndex(mu1, mu2) = E.Fem(E.Plain( BasisDataArray.D(bda), n, SOME(f)), id1, id2, id3, [mu1], [mu2])
	       val dim = BasisDataArray.domainDim bda
	       val _ = if dim <> 2 andalso dim <> 3
		       then raise Fail "One d and higher than 3d fields not yet supported"
		       else ()
	       val sumStart = sumX
	       val sumItter = List.tabulate(dim, fn x => x +  sumStart + 1)
	       val sumMus = List.map E.V sumItter
	       val sumRanges = List.map (fn x => (x, 0, dim - 1 )) sumItter (*QUESTION: dim - 1correct?*)
	       (*get mesh and get the dim*)
	       val accRanges = List.tabulate(dim, fn x => E.C x)

	       fun makeEps([a,b]) = E.Eps2(a,b)
		 | makeEps([a,b,c]) = E.Epsilon(a,b,c)
	       val eps = makeEps sumMus
	       val detAcces = List.map probeAtIndex (ListPair.zip(accRanges, sumMus))
	       val det = E.Op2(E.Div, E.Sum(sumRanges, E.Opn(E.Prod, eps::detAcces)), E.Const(fact dim))
	       (* fix sum indecies*)
	       val sumStartAdj = sumX + dim
	       val adjSumItter1 = List.tabulate(dim - 1, fn x => x + sumStartAdj + 1)
	       val adjSumItter2 = List.tabulate(dim - 1, fn x => x + sumStartAdj + dim + 1)
	       val adjSumMus1 = List.map E.V adjSumItter1
	       val adjSumMus2 = List.map E.V adjSumItter2
	       val adjSumRanges = List.map (fn x => (x, 0, dim - 1)) (List.@(adjSumItter1, adjSumItter2))
	       val finalSumX = sumStartAdj + 2 * (dim - 1) (*TODO: fix me*)
	       val eps1 = makeEps (List.@(adjSumMus1, [mu1]))
	       val eps2 = makeEps (List.@(adjSumMus2, [mu2]))
	       val adjAccess = List.map probeAtIndex (ListPair.zip(adjSumMus1, adjSumMus2))
	       val adj = E.Op2(E.Div, E.Sum(adjSumRanges, E.Opn(E.Prod, eps1::eps2::adjAccess)), E.Const(fact (dim - 1)))

	       val inv = E.Op2(E.Div, adj, det)
	      in
	       mkApply(E.Partial dx', inv, params, finalSumX)
	      end
            | E.Partial _ => err("Apply of Partial")
            | E.Apply(E.Partial d2, e2) => SOME(E.Apply(E.Partial(dx@d2), e2))
            | E.Apply _ => err "Apply of non-Partial expression"
            | E.Probe _ => err "Apply of Probe"
            | E.Value _ => err "Value used before expand"
            | E.Img _ => err "Probe used before expand"
            | E.Krn _ => err "Krn used before expand"
            | E.Comp(e1, [(e2, n)]) =>
	      let
	       (* I want to comment on why this works:*)
	       (* There are two key assumptions: *)
	       (* one: the type signature of compose means we only compose things with vector fields:*)
	       (* two: vector fields that are derivatives of other fields are treated as vector fields in eucldiean space and not on some manifold*)
	       (* meaning the transforms between types of indicies (or in our case world/index space) are absorbed into the expression for *)
	       (* a vector field by the act of pushing the derivatives to convo and taken into account properly when a prob of convo is expanded in high to mid.*)
	       (* With these assumptions, the chain rule reduces to the tensor-vector case every time.*)
	       val (d0::dn) = dx
               val vk = 100+sumX (* FIX ME: fresh index*)
               val e3 = E.Comp(E.Apply(E.Partial[E.V vk], e1), [(e2, n)])
               val e4 = E.Apply(E.Partial[d0], rewriteIx(vk, e2)) 
               val SOME(dim) = findDim(e1, params)
               val e5 = E.Sum([(vk, 0, dim-1)], E.Opn(E.Prod, [e3, e4]))
	       val e = if null dn then e5 else E.Apply(E.Partial dn, e5)
	      in SOME (e) end
            | E.Comp _ => err "unsupported differentiation of comp"
            | E.OField(ofld, e2, E.Partial alpha) => SOME(E.OField(ofld, e2, E.Partial(alpha @ dx)))
            | E.If(comp, e3, e4) => SOME(E.If(comp, E.Apply(E.Partial dx, e3), E.Apply(E.Partial dx, e4)))
            | E.Sum(sx, e1) => SOME(E.Sum(sx, E.Apply(E.Partial dx, e1)))
            | E.Op1(op1, e1) => SOME(applyOp1(op1, e1, dx))
            | E.Op2(op2, e1, e2) => SOME(applyOp2(op2, e1, e2, dx))
            | E.Op3(op3, e1, e2, e3) => SOME zero (*assume clamp is not lifted*)
            | E.Opn(opn, es) => SOME(applyOpn(opn, es, dx))
          (* end case *))

  end
