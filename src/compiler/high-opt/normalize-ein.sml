(* normalize-ein.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure NormalizeEin : sig

  (* normalize an Ein function; if there are no changes, then NONE is returned. *)
    val transform : HighIR.var list  -> Ein.ein -> Ein.ein option

  end = struct

    structure E = Ein
    structure ST = Stats
    structure IR = HighIR
    structure Ty = HighTypes

  (********** Counters for statistics **********)
    val cntNullSum              = ST.newCounter "high-opt:null-sum"
    val cntSumRewrite           = ST.newCounter "high-opt:sum-rewrite"
    val cntProbe                = ST.newCounter "high-opt:normalize-probe"
    val cntFilter               = ST.newCounter "high-opt:filter"
    val cntApplyPartial         = ST.newCounter "high-opt:apply-partial"
    val cntNegElim              = ST.newCounter "high-opt:neg-elim"
    val cntSubElim              = ST.newCounter "high-opt:sub-elim"
    val cntDivElim              = ST.newCounter "high-opt:div-elim"
    val cntDivDiv               = ST.newCounter "high-opt:div-div"
    val cntAddRewrite           = ST.newCounter "high-opt:add-rewrite"
    val cntSqrtElim             = ST.newCounter "high-opt:sqrt-elim"
    val cntEpsElim              = ST.newCounter "high-opt:eps-elim"
    val cntEpsToDeltas          = ST.newCounter "high-opt:eps-to-deltas"
    val cntNegDelta             = ST.newCounter "high-opt:neg-delta"
    val cntReduceDelta          = ST.newCounter "high-opt:reduce-delta"
    val cntIdentityProbe        = ST.newCounter "high-opt:identity-probe"
    val cntCompCancel           = ST.newCounter "high-opt:cntcompCancel"
    val cntInnerLoop            = ST.newCounter "high-opt:inner-loop"
    val firstCounter            = cntNullSum
    val lastCounter'            = cntCompCancel
    val lastCounter             = cntInnerLoop
    val cntRounds               = ST.newCounter "high-opt:normalize-round"



    fun err str = raise Fail(String.concat["Ill-formed EIN Operator",str])

    val zero = E.Const 0

    fun mkProd exps = E.Opn(E.Prod, exps)
    fun mkDiv (e1, e2) = E.Op2(E.Div, e1, e2)

  (* build a normalized summation *)
    fun mkSum ([], b) = (ST.tick cntNullSum; b)
      | mkSum (sx, b) = let
          fun return e = (ST.tick cntSumRewrite; e)
          in
            case b
             of E.Lift e => return (E.Lift(E.Sum(sx, e)))
              | E.Tensor(_, []) => return b
              | E.Zero [] => return b
              | E.Const _ => return b
              | E.ConstR _ => return b
              | E.Opn(E.Prod, es) => (case EinFilter.filterSca (sx, es)
                   of (true, e) => return e
                    | _ => E.Sum(sx, b)
                  (* end case *))
              | _ => E.Sum(sx, b)
            (* end case *)
          end

  (* build a normalized probe operation *)
    fun mkProbe (fld, x) = let
          fun return e = (ST.tick cntProbe; e)
          in
            case fld
             of E.Tensor _         => err "Tensor without Lift"
              | E.Lift e           => return e
              | E.Zero _           => return fld
              | E.Partial _        => err "Probe Partial"
              | E.Probe _          => err "Probe of a Probe"
              | E.Value _          => err "Value used before expand"
              | E.Img _            => err "Probe used before expand"
              | E.Krn _            => err "Krn used before expand"
              | E.Comp _           => E.Probe(fld, x) (* handled next stage*)
              | E.Epsilon _        => return fld
              | E.Eps2 _           => return fld
              | E.Const _          => return fld
              | E.Delta _          => return fld
              | E.If(E.Compare(op1, e1, e2), e3, e4) => let
                    val comp2 = E.Compare(op1, E.Probe(e1, x), E.Probe(e2, x))
                    in (ST.tick cntProbe; (E.If(comp2, E.Probe(e3, x), E.Probe(e4, x)))) end
              | E.If(E.Var id, e3, e4) 
                => (ST.tick cntProbe;  E.If(E.Var id, E.Probe(e3, x), E.Probe(e4, x)))
              | E.Sum(sx1, e)      => return (E.Sum(sx1, E.Probe(e, x)))
              | E.Op1(op1, e)      => return (E.Op1(op1, E.Probe(e, x)))
              | E.Op2(op2, e1, e2) => return (E.Op2(op2, E.Probe(e1, x), E.Probe(e2, x)))
              | E.Op3(op3, e1, e2, e3) =>
                  return (E.Op3(op3, E.Probe(e1, x), E.Probe(e2, x), E.Probe(e3, x)))
              | E.Opn(opn, [])     => err "Probe of empty operator"
              | E.Opn(opn, es)     => return (E.Opn(opn, List.map (fn e => E.Probe(e, x)) es))
              | _                  => E.Probe(fld, x)
            (* end case *)
          end

    (* rewrite expression with composition operation *)
    fun mkComp(F, es, x) = let
        fun return e = (ST.tick cntProbe; e)
        fun setInnerProbe e = E.Probe(E.Comp(e, es), x)
        val probe = setInnerProbe F
        in (case F
            of E.Tensor _        => err "Tensor without Lift"
            | E.Lift e           => return e
            | E.Zero _           => return F
            | E.Partial _        => err "Probe Partial"
            | E.Probe _          => err "Probe of a Probe"
            | E.Value _          => err "Value used before expand"
            | E.Img _            => err "Probe used before expand"
            | E.Krn _            => err "Krn used before expand"
            | E.Comp c           => probe (* handled next stage*)
            | E.Epsilon _        => return F
            | E.Eps2 _           => return F
            | E.Const _          => return F
            | E.Delta _          => return F
            | E.Sum(sx1, e)      => return (E.Sum(sx1, setInnerProbe e))
            | E.Op1(op1, e)      => return (E.Op1(op1, setInnerProbe e))
            | E.Op2(op2, e1, e2) => let
                val exp1 = setInnerProbe e1
                val exp2 = setInnerProbe e2
                val xexp = E.Op2(op2, exp1, exp2)
                in return xexp end
            | E.Opn(opn, [])     => err "Probe of empty operator"
            | E.Opn(opn, es1)     => let
                val exps =  List.map (fn e1 => E.Probe(E.Comp(e1, es), x)) es1
                val xexp = E.Opn(opn, exps)
                in return xexp end
            | _                  => probe
            (* end case *))
    end

    val basisTable : (BasisDataArray.t, BasisDataArray.t) HashTable.hash_table = HashTable.mkTable(BasisDataArray.hash, BasisDataArray.same) (256, Fail "basis D not found")
    fun getBasisDerivative(B : BasisDataArray.t) : BasisDataArray.t =
	(case HashTable.find basisTable B
	  of SOME(b) => b
	   | NONE =>
	     let val b = BasisDataArray.D(B) in (HashTable.insert basisTable (B,b); b) end
	(* end case*))



  (* rewrite body of EIN *)
    fun transform vars (ein as Ein.EIN{params, index, body}) = let
     (* val _ = print(String.concat["\ntransform: ", EinPP.expToString(body)]) *)
     val compDepth = ref 0
          fun filterProd args = (case EinFilter.mkProd args
                 of SOME e => (ST.tick cntFilter; e)
                  | NONE => mkProd args
                (* end case *))
          fun filterProdNoCnt args = (case EinFilter.mkProd args
                 of SOME e => e
                  | NONE => mkProd args
                (* end case *))
          val sumX = ref (length index)
	  val debug = true
          fun incSum() = (sumX:= (!sumX+2) handle ex => raise ex; if debug
								  then print("new sumX:"^(Int.toString (!sumX))^"\n")
								  else ())
         fun addSum((v, _, _)::sx) =
             ((sumX:= (!sumX)+v handle ex => (print("adding sum:"^(Int.toString v)^" to " ^ (Int.toString (!sumX))^"\n"); raise ex));
	      if debug
	      then print("raised sumX by "^ (Int.toString v)^" to "^(Int.toString (!sumX))^"\n")
	      else ()
	     )
          fun rewrite body = (case body
                 of E.Const _ => body
                  | E.ConstR _ => body
                  | E.Tensor _ => body
                  | E.Zero _ => body
                  | E.Delta(E.C(1), E.C(1)) => (ST.tick cntReduceDelta; E.Const 1)
                  | E.Delta(E.C(0), E.C(0)) => (ST.tick cntReduceDelta; E.Const 1)
                  | E.Delta(E.C(0), E.C(1)) => (ST.tick cntReduceDelta; E.Const 0)
                  | E.Delta(E.C(1), E.C(0)) => (ST.tick cntReduceDelta; E.Const 0)
                  | E.Delta _ => body
                  | E.Epsilon _ => body
                  | E.Eps2 _ => body
                (************** Field Terms **************)
                  | E.Field _ => body
                  | E.Lift e1 => E.Lift(rewrite e1)
                  | E.Conv _ => body
		  | E.Identity _ => body
		  | E.Fem _ => body
                  | E.Partial _ => body
                  | E.Apply(E.Partial [], e1) => e1
                  | E.Apply(E.Partial d1, e1) => let
                      val e1 = rewrite e1
                      in
                        (case (Derivative.mkApply(E.Partial d1, e1, params, sumX, getBasisDerivative) handle ex=>raise ex) 
                         of SOME e => (ST.tick cntApplyPartial; e)
                          | NONE => E.Apply(E.Partial d1, e1)
                        (* end case *)) handle ex => raise ex
                      end 
                  | E.Apply _ => err "Ill-formed Apply expression"
                (************** Field Terms **************)
                  | E.OField(ofld, e, alpha) => E.OField(ofld, rewrite e, alpha)
                  | E.Value _ => err "Value before Expand"
                  | E.Img _ => err "Img before Expand"
                  | E.Krn _ => err "Krn before Expand"
                  | E.Poly _ => err "Poly before Expand"
                  (************** Composition **************)
                  | E.Comp(E.Comp(a, es1), es2) => (ST.tick cntProbe; rewrite (E.Comp(a, es1@es2)))
                  | E.Comp(a, (E.Comp(b, es1), m)::es2) =>  (ST.tick cntProbe; rewrite (E.Comp(a, ((b, m)::es1)@es2)))
                  | E.Comp(e1, es)  =>
                    let

		     (* val _ = print(" Start comp " ^ (Int.toString(!compDepth)) ^ ":"^(EinPP.expToString(E.Comp(e1,es)))^"\n") *)
		     (* val _ = compDepth := (!compDepth) + 1 handle ex => raise ex *)
		    (* Rewrite e1 and es fully:*)
		     val needRewrite = ref false
		     fun getOpt(a,b) = (case a
					 of SOME(a') => (needRewrite := true; a')
					  | NONE => b)
                     val e1' = innerLoop(e1, ST.sum{from = firstCounter, to = lastCounter'}, false) handle ex => raise ex
		     val e1'' = getOpt(e1', e1)
                     val es' = List.map (fn (e2, n2)=> (innerLoop(e2, ST.sum{from = firstCounter, to = lastCounter'}, false), n2)) es handle ex => raise ex
		     val es'' = ListPair.map (fn ((e2, n2), (e2', n2')) => (getOpt(e2', e2), n2)) (es, es')

		     (* If any rewrites occur, we try to rewrite the whole thing further*)
		     val E.Comp(e1''', es''') = if !needRewrite
						then (case innerLoop(E.Comp(e1'', es''), ST.sum{from = firstCounter, to = lastCounter'}, false)
						       of SOME(a) => a
							| NONE => E.Comp(e1'', es'')
						     (* end case*))
						else E.Comp(e1'', es'')
		     val compRet = E.Comp(e1''', es''')

		     (*We now setup the cancellation: this function will, based on eIN, 
		       replace field ins in eOut with an identity or a comp!
		     This will return: 
		     -the resulting eOut
		     -a flag for if any comp replaces happened
		     -a flag for if any identity replaces happened
		     -a flag if eIn couldn't result in any cancells
		      *)
		     fun tryCancel(eOut : E.ein_exp, eIn : E.ein_exp, eInBind : E.index_bind list) =
		     	 let
		     	  val fail = ref false
		     	  val aRet = ref false
		     	  fun sucRet e = (aRet := true; e)
		     	  fun failRet e = (fail := true; E.Comp(e, [(eIn, eInBind)]))
		     	  fun vTy x = IR.Var.ty (List.nth(vars, x))
		     	  fun vS(x,y) = x=y orelse IR.Var.same(List.nth(vars, x), List.nth(vars, y))
		     	  fun findInvert(mesh, index, indexSource, dofSource) e = (*T^-1 \circ T_i*)
		     	      (case e
		     		of (E.Fem(E.Invert(_, _, NONE, w), index', indexSource', dofSource', [acc], [])) (*trf case - w is false*)
		     		   => if vS(index, index') andalso vS(indexSource', indexSource)
		     			 andalso vS(dofSource, dofSource') andalso (if w
										    then raise Fail "impossible"
										    else true)
		     		      then sucRet(E.Identity(FemData.meshDim mesh, acc, SOME(false, index)))
		     		      else failRet e
		     		 | (E.Fem(E.Invert(_, _, SOME _, w), index', indexSource', dofSource', [acc], [])) => (*Tinv or tF*)
		     		   if vS(indexSource', indexSource) andalso vS(dofSource, dofSource')
				      andalso (case (vTy index')
						of HighTypes.FemData(FemData.Mesh _) => true
											andalso if w
												then true
												else raise Fail "invalid global invert"
						 (*true except on boundary, which we ignore - w is true here.*)
						 | HighTypes.IntTy => vS(index', index)
								      andalso if w
									      then raise Fail "invalid local invert producing world space"
									      else true
						 (*T_j^-1 \circ T_i \neq id unless i=j - w is false here*)
						 | _ => raise Fail "impossible"
					      (*end case*))
		     		   then sucRet(E.Identity(FemData.meshDim mesh, acc, SOME(false, index)))
		     		   else failRet e
		     		 | _ => failRet e
		     	      (* end case*))
		     	  fun findPlain(mesh, index, indexSource, dofSource) e = (*T_i \circ T^-1*)
		     	      (case e
		     		of E.Fem(E.Plain(_, _, NONE), index', indexSource', dofSource', [acc], []) (*T_i - the only case*)
		     		   => if vS(index, index') andalso vS(indexSource', indexSource)
					 (*index test exludes Global inverse, which is impossible anyway as that only lives in F*)
		     			 andalso vS(dofSource, dofSource')
		     		      then sucRet(E.Identity(FemData.meshDim mesh, acc, SOME(true, index)))
		     		      else failRet e
		     		 | _ => failRet e
		     	      (* end case*))
		     	  val f : (Ein.ein_exp -> Ein.ein_exp) option =
		     	      (case eIn
		     		of E.Fem(femEin, index, indexSource, dofSource, [E.V _], []) =>
		     		   (case (femEin, vTy indexSource)
		     		     of (E.Plain(_, _, _), Ty.FemData(FemData.Mesh mesh))
		     			=> SOME(findInvert(mesh, index, indexSource, dofSource))
		     		      | (E.Invert(_, _, _, _), Ty.FemData(FemData.Mesh mesh)) =>
		     			SOME(findPlain(mesh, index, indexSource, dofSource))
		     		      | _ => NONE
		     		   (* end case*))
		     		 | _ => NONE
		     	      (* end ase*))
		     	  val eOut' = Option.map (fn f' => EinUtil.mapInNodes(eOut, f')) f
		     	  val possibleReplace = Bool.not (Option.isSome f)
		     	 in
		     	  (Option.getOpt(eOut', eOut), !fail, !aRet, possibleReplace)
		     	 end
		     (* (*We do the canellation directly with this function, marking it here:*) *)
		     (* val cancelRef = ref false *)
		     (* fun doCompLefts(e1::e2::es) = *)
		     (* 	 let *)
		     (* 	  val (e1e, e1b) = e1 *)
		     (* 	  val (e2e, e2b) = e2 *)
		     (* 	  val (e1e', anyFail, anyReplace, nothing) = tryCancel(e1e, e2e, e2b) *)
		     (* 	  val _ = if Bool.not anyFail andalso Bool.not anyReplace *)
		     (* 		  then raise Fail ("impossible:ill-formed comp: " ^ (EinPP.expToString(e1e)) ^ " and " ^ (EinPP.expToString(e2e)) ^"\n") *)
		     (* 		  else () *)
		     (* 	 in *)
		     (* 	  if nothing orelse anyFail (*Question: might be a good idea to allow any-fail*) *)
		     (* 	  then e1::doCompLefts(e2::es) *)
		     (* 	  else (cancelRef := true; ST.tick cntCompCancel; (e1e', e1b) :: doCompLefts(es)) *)
		     (* 	 end *)
		     (*   | doCompLefts ([e]) = [e] *)
		     (*   | doCompLefts ([]) = [] *)
		     (* val allEs = (e1''', []) :: es''' *)
		     (* val allEs' = doCompLefts(allEs) *)
			     
		     (* (*Now we prep to filter identities*) *)
		     (* val filterId = ref false *)
		     (* fun filterIdentities(alles) = *)
		     (* 	 let *)
		     (* 	  val n = List.length alles *)
		     (* 	  fun filterFn((E.Identity(_, E.V _, _), _)) = false *)
		     (* 	    | filterFn _ = true *)
		     (* 	  val alles' = List.filter filterFn alles *)
		     (* 	  val possibleId = List.find (Bool.not o filterFn) alles *)
		     (* 	  val _ = filterId := (n <> 1 andalso n <> (List.length alles')) *)
		     (* 	 in *)
		     (* 	  (case alles' *)
		     (* 	    of [] => (case possibleId *)
		     (* 		       of SOME(pid)=> [pid] *)
		     (* 			| _ => raise Fail "impossible:empy alles with no identity" *)
		     (* 		     (* end case *)) *)
		     (* 	     | es => es *)
		     (* 	  (* end case *)) *)
		     (* 	 end *)
		     (* val allEs'' = filterIdentities(allEs') *)
		     (* (*We produce the resulting ein here:*) *)
		     (* val compRet = (case allEs'' *)
		     (* 		     of [] => raise Fail "impossible: fail comp cancel somehow." *)
		     (* 		      | [(r,rbind)] => r *)
		     (* 		      | (r,rbind)::rs => E.Comp(r,rs) *)
		     (* 		   (*end case*)) *)


		     (* val compRet = if !cancelRef orelse !filterId *)
		     (* 		   then *)
		     (* 		    Option.getOpt(innerLoop(E.Comp(e1'', es''), ST.sum{from = firstCounter, to = lastCounter'}, false), *)
		     (* 				  compRet) *)
				     (* 		   else compRet *)
		     (* val _ = compDepth := (!compDepth) - 1 *)
		     (* val _ = print("End comp " ^ (Int.toString(!compDepth)) ^ ":" ^ (EinPP.expToString(compRet))^"\n") *)


                    in  compRet end
                  | E.Probe(E.Comp(e1, es), x)  =>
                    let
                    val e1' = rewrite e1
                    val es' = List.map (fn (e2, n2)=> (rewrite e2, n2)) es
                    in (case (rewrite(E.Comp(e1', es')))
                        of E.Comp(e1', es') => mkComp(e1', es', x)
                        | e => e
                        (*end case*))
                    end
		  | E.Probe(E.Identity(dim, mu1, _), e as E.Tensor(tid, [])) => (*Id(e2) -> Id * e2 \sum_i=(0,dim-1) delta_ij e_2_i...*)
		    let
		     (* val _ = print("inner:" ^ (EinPP.expToString e) ^ "\n") *)
		     val newSumRange = !sumX + 1 handle ex => raise ex
		     val _ = incSum()
		     val ret = mkSum([(newSumRange, 0, dim - 1)], E.Opn(E.Prod, [E.Delta(mu1, E.V newSumRange), E.Tensor(tid, [E.V newSumRange])]))
		     (* val new = print("inner':" ^ (EinPP.expToString ret)) *)
		    in
		     (ST.tick cntIdentityProbe; ret)
		    end
                  | E.Probe(e1, e2) => mkProbe(rewrite e1, rewrite e2)
                (************** Sum **************)
                  | E.If(E.Compare(op1, e1, e2), e3, e4) => E.If(E.Compare(op1, rewrite e1, rewrite e2), rewrite e3, rewrite e4)
                  | E.If(E.Var id, e3, e4) => E.If(E.Var id, rewrite e3, rewrite e4)
                  | E.Sum(sx, e) => (addSum(sx); mkSum (sx, rewrite e))
                (************* Algebraic Rewrites Op1 **************)
                  | E.Op1(E.Neg, E.Op1(E.Neg, e)) => (ST.tick cntNegElim; rewrite e)
                  | E.Op1(E.Neg, E.Const 0) => (ST.tick cntNegElim; zero)
                  | E.Op1(E.Neg, e1 as E.Zero _) => (ST.tick cntNegElim; e1)
                  | E.Op1(op1, e1) => E.Op1(op1, rewrite e1)
                (************* Algebraic Rewrites Op2 **************)
                  | E.Op2(E.Sub, E.Const 0, e2) => (ST.tick cntSubElim; E.Op1(E.Neg, rewrite e2))
                  | E.Op2(E.Sub, e1, E.Const 0) => (ST.tick cntSubElim; rewrite e1)
                  | E.Op2(E.Sub, e1 as E.Zero _, e2) => (ST.tick cntSubElim; e1)
                  | E.Op2(E.Sub, e1, e2 as E.Zero _) => (ST.tick cntSubElim; e2)
                  | E.Op2(E.Div, E.Const 0, e2) => (ST.tick cntDivElim; zero)
                  | E.Op2(E.Div, e1 as E.Zero _, e2) => (ST.tick cntDivElim; e1)
                  | E.Op2(E.Div, E.Op2(E.Div, a, b), E.Op2(E.Div, c, d)) => (
                      ST.tick cntDivDiv;
                      rewrite (mkDiv (mkProd[a, d], mkProd[b, c])))
                  | E.Op2(E.Div, E.Op2(E.Div, a, b), c) => (
                      ST.tick cntDivDiv;
                      rewrite (mkDiv (a, mkProd[b, c])))
                  | E.Op2(E.Div, a, E.Op2(E.Div, b, c)) => (
                      ST.tick cntDivDiv;
                      rewrite (mkDiv (mkProd[a, c], b)))
                 (************** min|max **************)
                  | E.Op2(E.Min, e1, e2)     =>
                    let
                        val comp = E.Compare(E.LT, e1, e2)
                        val exp  = E.If(comp, e1, e2)
                    in (ST.tick cntProbe; exp) end
                  | E.Op2(E.Max, e1, e2)     =>
                    let
                        val comp = E.Compare(E.GT, e1, e2)
                        val exp  = E.If(comp, e1, e2)
                    in (ST.tick cntProbe; exp) end
                  | E.Op2(op1, e1, e2) => E.Op2(op1, rewrite e1, rewrite e2)
                  | E.Op3(op3, e1, e2, e3) =>
                      E.Op3(op3, rewrite e1, rewrite e2, rewrite e3)
                (************* Algebraic Rewrites Opn **************)
                  | E.Opn(E.Add, es) => let
                      val es' = List.map rewrite es
                      in
                        case EinFilter.mkAdd es'
                         of SOME body' => (ST.tick cntAddRewrite; body')
                          | NONE => E.Opn(E.Add, es')
                      end
                (************* Product **************)
                  | E.Opn(E.Prod, []) => err "missing elements in product"
                  | E.Opn(E.Prod, [e1]) => rewrite e1
                  | E.Opn(E.Prod, [e1 as E.Op1(E.Sqrt, s1), e2 as E.Op1(E.Sqrt, s2)]) =>
                      if EinUtil.sameExp(s1, s2)
                        then (ST.tick cntSqrtElim; s1)
                        else filterProd [rewrite e1, rewrite e2]
                (************* Product EPS **************)
                  | E.Opn(E.Prod, (eps1 as E.Epsilon(i,j,k))::ps) => (case ps
                       of ((p1 as E.Apply(E.Partial (d as (_::_::_)), e)) :: es) => (
                            case (EpsUtil.matchEps ([i,j,k], d), es)
                             of (true, _) => (ST.tick cntEpsElim; zero)
                              | (_, []) => mkProd[eps1, rewrite p1]
                              | _ => filterProd [eps1, rewrite (mkProd (p1 :: es))]
                            (* end case *))
                        | ((p1 as E.Conv(_, _, _, (d as (_::_::_)))) :: es) => (
                            case (EpsUtil.matchEps ([i,j,k], d), es)
                             of (true, _) => (ST.tick cntEpsElim; E.Lift zero)
                              | (_, []) => mkProd[eps1, p1]
                              | _ => (case rewrite (mkProd(p1 :: es))
                                   of E.Opn(E.Prod, es') => filterProd (eps1 :: es')
                                    | e => filterProd [eps1, e]
                                  (* end case *))
                        (* end case *))
			(* Recreation of convo rule but for Fem's convo
			 TODO: shouldn't the idea of Convo as simply a wrapper for special computations be implemented? *)		
			| ((p1 as E.Fem(_, _, _,_,_, (d as (_::_::_)))) :: es) => (
                         case (EpsUtil.matchEps ([i,j,k], d), es)
                          of (true, _) => (ST.tick cntEpsElim; E.Lift zero)
                           | (_, []) => mkProd[eps1, p1]
                           | _ => (case rewrite (mkProd(p1 :: es))
                                    of E.Opn(E.Prod, es') => filterProd (eps1 :: es')
                                     | e => filterProd [eps1, e]
                                  (* end case *))
                            (* end case *))
                        | [E.Tensor(_, [i1, i2])] =>
                            if (j=i1 andalso k=i2)
                              then (ST.tick cntEpsElim; zero)
                              else body
                        | _  => (case EpsUtil.epsToDels (eps1::ps)
                             of (SOME(e, sx), NONE, rest) => (case sx
                                (* Changed to Deltas*)
                                   of [] => (
                                        ST.tick cntEpsToDeltas;
                                        E.Opn(E.Prod, e::rest))
                                    | _  => (
                                        ST.tick cntEpsToDeltas;
                                        E.Opn(E.Prod, E.Sum(sx, e)::rest))
                                  (* end case *))
                              | (SOME _ , _ , _) => raise Fail "not possible"
                              | (NONE, NONE, rest) => raise Fail "not possible"
                              | (NONE, SOME(epsAll), rest) => (case rest
                                   of [] => (body) (* empty eps-product and empty rest no change *)
                                    | [r1] => (E.Opn(E.Prod, epsAll@[rewrite r1]))
                                    | _ => (case rewrite(E.Opn(E.Prod, rest))
                                         of E.Opn(E.Prod, p) => (E.Opn(E.Prod, epsAll@p))
                                          | r2 => (E.Opn(E.Prod, epsAll@[r2]))
                                        (* end case *))
                                    (* end case *))
                            (* end case *))
                      (* end case *))
                  | E.Opn(E.Prod, E.Delta d::es) => (case es
                       of [E.Op1(E.Neg, e1)] => (
                            ST.tick cntNegDelta; E.Op1(E.Neg, mkProd[E.Delta d, e1]))
                        | _ => let
                            val (pre', eps, dels, post) = EinFilter.partitionGreek(E.Delta d::es)
                            in
                              case EpsUtil.reduceDelta(eps, dels, post)
                               of (false, _) => (case (rewrite(mkProd es))
                                     of E.Opn(E.Prod, p) => mkProd (E.Delta d::p)
                                      | e2 => mkProd [E.Delta d,  e2]
                                    (* end case*))
                              (*  | (_, E.Opn(E.Prod, p)) => (ST.tick cntReduceDelta; filterProd p)*)
                                | (_, a) => (ST.tick cntReduceDelta; a)
                              (* end case *)
                            end
                      (* end case *))
                (************* Product Generic **************)
                  | E.Opn(E.Prod, [e1 as E.Zero alpha, e2]) =>
                      if (EinFilter.isScalar e2)
                        then E.Zero alpha
                        else filterProd [rewrite e1, rewrite e2]
                  | E.Opn(E.Prod, [e1, e2 as E.Zero alpha]) =>
                      if (EinFilter.isScalar e1)
                        then E.Zero alpha
                        else filterProd [rewrite e1, rewrite e2]
                  | E.Opn(E.Prod, [e1, e2]) => filterProd [rewrite e1, rewrite e2]
                  | E.Opn(E.Prod, e1::es) => let
                      val e' = rewrite e1
                      val e2 = rewrite (mkProd es)
                      in
                        case e2
                         of E.Opn(Prod, p') => filterProd (e' :: p')
                          | _ => filterProd [e',e2]
                        (* end case *)
                      end
			     (* end case *))
	  and innerLoop(exp, total, changed) =
	      let
	       val exp' = rewrite exp handle ex => raise ex
	       val totalTicks = ST.sum{from = firstCounter, to = lastCounter'} handle ex => raise ex
	      in
	       if totalTicks > total
	       then innerLoop(exp', totalTicks, true) handle ex => raise ex
	       else if changed
	       then (ST.tick cntInnerLoop; SOME(exp'))
	       else NONE
	      end
(* (*DEBUG*)val start = ST.count cntRounds *)
          fun loop (body, total, changed) = let
	   val _ = ST.report()
                val body' = rewrite body
		val _ =print(String.concat["\n\n ==> X:", EinPP.expToString(body),"\n ==> Y:", EinPP.expToString(body'), "\n"])
                val totalTicks = ST.sum{from = firstCounter, to = lastCounter}
		val _ = print("Total Ticks:"^(Int.toString totalTicks) ^"\n")
                in
                  ST.tick cntRounds;
		  (* (*DEBUG*)if (ST.count cntRounds - start > 50) then raise Fail "too many steps" else (); *)
                  if (totalTicks > total) (* something changed *)
                    then loop(body', totalTicks, true) handle ex => raise ex (* keep going *)
                  else if changed (* nothing changed - if any changes happened at all, provide result*)
                    then SOME(Ein.EIN{params=params, index=index, body=body'})
                    else NONE
          end
	  val einRet = loop(body, ST.sum{from = firstCounter, to = lastCounter}, false) handle ex => raise ex
	  (* val _ = Option.app (fn x => print(String.concat(["ret:", EinPP.toString x, "\n"]))) einRet *)
          in
            einRet
          end

  end
