(* ein-util.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure EinUtil : sig

  (* Are Ein functions/expressions the same? *)
    val same : Ein.ein * Ein.ein -> bool
    val sameExp : Ein.ein_exp * Ein.ein_exp -> bool

  (* compute a hash for an Ein function *)
    val hash : Ein.ein -> Word.word

    val iterPP: Ein.ein_exp list -> Ein.ein_exp
    val iterAA: Ein.ein_exp list -> Ein.ein_exp

    (*Used to replace on the inputs to an ein_exp with another; used to do cancellations in Comp*)
    val mapInNodes : Ein.ein_exp * (Ein.ein_exp -> Ein.ein_exp) -> Ein.ein_exp
    val analyzeInNodes : Ein.ein_exp * (Ein.ein_exp -> 'a option) * ('a * 'a -> 'a) -> 'a option
    val detectInNodes : Ein.ein_exp * (Ein.ein_exp -> bool) -> bool
    val cleanIds : Ein.ein_exp -> Ein.ein_exp

    val collectIndicies : Ein.ein_exp -> int list

    (* val typeCheckTenShape : Ein.ein -> ((int list * int option) * bool)  (*shape, interval, matches internal shape*) *)
    (* val typeCheckTenShapeExp : Ein.ein_exp * (int -> ((int list * int option)) option) -> (int list * int option) *)
  end = struct

    structure E = Ein

    fun sameExp (e1, e2) = let
        fun sameIndex ([], []) = true
          | sameIndex ((E.V i)::ix, (E.V j)::jx) = (i = j) andalso sameIndex (ix, jx)
          | sameIndex ((E.C i)::ix, (E.C j)::jx) = (i = j) andalso sameIndex (ix, jx)
          | sameIndex _ = false
        fun sameKx ([], [])=true
          | sameKx ((E.V i,_)::ix, (E.V j,_)::jx) = (i = j) andalso sameKx(ix, jx)
          | sameKx ((E.C i,_)::ix, (E.C j,_)::jx) = (i = j) andalso sameKx(ix, jx)
          | sameKx _ = false
        fun sameSx ([], []) = true
          | sameSx ((i,_,_)::ix, (j,_,_)::jx) = (i = j) andalso sameSx(ix, jx)
          | sameSx _ = false
        fun sameOp1 (E.Neg, E.Neg) = true
          | sameOp1 (E.Exp, E.Exp) = true
          | sameOp1 (E.Sqrt, E.Sqrt) = true
          | sameOp1 (E.Cosine, E.Cosine) = true
          | sameOp1 (E.ArcCosine, E.ArcCosine) = true
          | sameOp1 (E.Sine, E.Sine) = true
          | sameOp1 (E.ArcSine, E.ArcSine) = true
          | sameOp1 (E.Tangent, E.Tangent) = true
          | sameOp1 (E.ArcTangent, E.ArcTangent) = true
          | sameOp1 (E.PowInt n1, E.PowInt n2) = (n1 = n2)
          | sameOp1 (E.Abs, E.Abs) = true
          | sameOp1 (E.Sgn, E.Sgn) = true
          | sameOp1 _ = false
        fun same (e1, e2) = (case (e1, e2)
               of (E.Const c1, E.Const c2) => (c1 = c2)
                | (E.ConstR r1, E.ConstR r2) => Rational.same(r1, r2)
                | (E.Tensor(id1, ix1), E.Tensor(id2, ix2)) =>
                    (id1 = id2) andalso sameIndex(ix1, ix2)
                | (E.Zero(ix1), E.Zero(ix2)) => sameIndex(ix1, ix2)
                | (E.Delta(ix1, jx1), E.Delta(ix2, jx2)) =>
                    (ix1 = ix2) andalso (jx1 = jx2)
                | (E.Epsilon(i1, j1, k1), E.Epsilon(i2, j2, k2)) =>
                    (i1 = i2) andalso (j1 = j2) andalso (k1 = k2)
                | (E.Eps2(i1, j1), E.Eps2(i2, j2)) => (i1 = i2) andalso (j1 = j2)
                | (E.Field(id1, ix1), E.Field(id2, ix2)) =>
                    (id1 = id2) andalso sameIndex(ix1, ix2)
                | (E.Lift e1, E.Lift e2) => same (e1, e2)
                | (E.Conv(fid1, alpha1, tid1, ix1), E.Conv(fid2, alpha2, tid2, ix2)) =>
                    (fid1 = fid2) andalso (tid1 = tid2) andalso
                    sameIndex (alpha1, alpha2) andalso sameIndex (ix1, ix2)
		| (E.Identity(d1, mu1, NONE), E.Identity(d2, mu2, NONE)) => d1=d2 andalso sameIndex([mu1], [mu2])
		| (E.Identity(d1, mu1, SOME(b,idx)), E.Identity(d2, mu2, SOME(b',idx'))) =>
		  d1=d2 andalso sameIndex([mu1], [mu2]) andalso b=b' andalso idx=idx'
		| (E.Fem(femEin1, cell1, index1, dof1, shape1, dxes1),E.Fem(femEin2, cell2, index2, dof2, shape2, dxes2)) =>
		  cell1=cell2 andalso index1=index2 andalso dof1=dof2
		  andalso sameIndex(shape1,shape2) andalso sameIndex(dxes1,dxes2)
		  andalso
		  (case (femEin1, femEin2)
		    of (E.Plain(bda1,n1, SOME(s1, s1')), E.Plain(bda2,n2, SOME(s2, s2'))) =>
		       (n1=n2)
		       andalso (BasisDataArray.same ( bda1, bda2))
		       andalso Stamp.same(s1,s2) andalso Stamp.same(s1', s2')
		    |  (E.Plain(bda1,n1, NONE), E.Plain(bda2,n2, NONE)) => 
		       (n1=n2)
		       andalso (BasisDataArray.same ( bda1, bda2))
		    | (E.Invert(bda1,n1, s1, w1), E.Invert(bda2,n2, s2, w2)) =>
		       (n1=n2)
		       andalso (BasisDataArray.same ( bda1, bda2))
		       andalso (w1=w2)
		       andalso (case (s1,s2)
				 of (SOME(s1',s1''),SOME(s2',s2'')) => Stamp.same(s1',s2') andalso Stamp.same(s1'', s2'')
				  | (NONE, NONE) => true
				  | _ => false (* end case*))
		    | _ => false
		  (*end case*))
		
		| (E.Partial ix, E.Partial jx) => sameIndex(ix, jx)
                | (E.Apply(e11, e12), E.Apply(e21, e22)) => same(e11, e21) andalso same(e12, e22)
                | (E.Comp(e11, es1), E.Comp(e21, es2)) =>
                    same(e11, e21) andalso sameSubEin(es1, es2)
                | (E.Probe(e11, e12), E.Probe(e21, e22)) => same(e11, e21) andalso same(e12, e22)
                | (E.Value i, E.Value j) => (i = j)
                | (E.Img(id1, ix1, pos1, s1), E.Img(id2, ix2, pos2, s2)) =>
                    (id1 = id2) andalso sameList(pos1, pos2) andalso sameIndex(ix1, ix2) andalso (s1 = s2)
                | (E.Krn(id1, ix1, dim1), E.Krn(id2, ix2, dim2)) =>
                    (id1 = id2) andalso sameKx(ix1, ix2) andalso (dim1 =  dim2)
                | (E.OField(E.CFExp (es1), e1, ix1), E.OField(E.CFExp (es2), e2, ix2)) =>
                    same(e1, e2) andalso same(ix1, ix2) andalso ListPair.allEq (op =) (es1, es2)
                | (E.Poly(e1, n1,alpha1), E.Poly(e2, n2, alpha2)) =>
                     same(e1, e2) andalso (n1 = n2) andalso sameIndex(alpha1, alpha2)
                | (E.Sum(c1, e1), E.Sum(c2, e2)) => sameSx(c1, c2) andalso same(e1, e2)
                | (E.Op1(op1, e1), E.Op1(op2, e2)) => sameOp1(op1, op2) andalso same(e1, e2)
                | (E.Op2(op1, e11, e12), E.Op2(op2, e21, e22)) =>
                    (op1 = op2) andalso same(e11, e21) andalso same(e12, e22)
                | (E.Op3(op1, e11, e12, e13), E.Op3(op2, e21, e22, e23)) =>
                    (op1 = op2) andalso same(e11, e21) andalso same(e12, e22) andalso same(e13, e23)
                | (E.Opn(op1, es1), E.Opn(op2, es2)) =>
                    (op1 = op2) andalso sameList(es1, es2)
                | _ => false
              (* end case *))
        and sameSubEin([], []) = true
            | sameSubEin ((e1 ,_)::es1, (e2, _)::es2) = same(e1, e2) andalso sameSubEin(es1, es2)
            | sameSubEin _ = false
        and sameList ([], []) = true
          | sameList (e1::es1, e2::es2) = same(e1, e2) andalso sameList(es1, es2)
          | sameList _ = false
        in
          same (e1, e2)
        end

    fun same (E.EIN{params=p1, index=ix1, body=e1}, E.EIN{params=p2, index=ix2, body=e2}) = let
          fun sameParam (E.TEN(i1, shp1, NONE), E.TEN(i2, shp2, NONE)) =
              (i1 = i2) andalso ListPair.allEq (op =) (shp1, shp2)
	    | sameParam (E.TEN(i1, shp1, SOME k1), E.TEN(i2, shp2, SOME k2)) =
	      (k1=k2) andalso i1=i2 andalso ListPair.allEq (op =) (shp1, shp2)
            | sameParam (E.FLD i1, E.FLD i2) = (i1 = i2)
            | sameParam (E.KRN, E.KRN) = true
            | sameParam (E.IMG(i1, shp1), E.IMG(i2, shp2)) =
              (i1 = i2) andalso ListPair.allEq (op =) (shp1, shp2)
	    | sameParam (E.INT, E.INT) = true
	    | sameParam (E.FEM(d1), E.FEM(d2)) = FemData.same(d1,d2)
            | sameParam _ = false
          in
            ListPair.allEq sameParam (p1, p2) andalso
              ListPair.allEq (op =) (ix1, ix2) andalso
                sameExp (e1, e2)
          end

    fun hash (Ein.EIN{body, ...}) = let
        fun hash' body = let
            fun hashInt i = Word.fromInt i
            fun iter [e] = hash' e
              | iter (e1::es) = hash' e1 + iter es
            fun iterS [(e,_)] = hash' e
              | iterS ((e1,_)::es) = hash' e1 + iterS es
            fun hashMu (E.C c) = hashInt c + 0w17
              | hashMu (E.V v) = hashInt v
            fun hashAlpha [] = 0w3
              | hashAlpha (e1::es) = hashMu e1 + hashAlpha es
            fun hashDels [] = 0w5
              | hashDels ((i, j)::es) = hashMu i + hashMu j + hashDels es
            in
              case body
               of E.Const i => hashInt i + 0w3
                | E.ConstR _ => 0w5
                | E.Tensor(_, alpha) => 0w7 + hashAlpha alpha
                | E.Zero(alpha) => 0w11 + hashAlpha alpha
                | E.Delta _ => 0w17
                | E.Epsilon _ => 0w19
                | E.Eps2 _ => 0w23
                | E.Field(_, alpha) => 0w29 + hashAlpha alpha
                | E.Lift e1 => 0w31 + hash' e1
                | E.Conv(_, alpha, _, dx) =>
                    0w37 + hashAlpha alpha + hashAlpha dx + hashInt(length dx)
                | E.Partial alpha => 0w41+hashAlpha alpha
                | E.Apply(e1, e2) => 0w43 + hash' e1 + hash' e2
                | E.Comp(e1, es) => 0w47 + hash' e1 + iterS es
                | E.Probe(e1, e2) => 0w53 + hash' e1 + hash' e2
                | E.Value _ => 0w59
                | E.Img (_, alpha, es, _) => 0w61 + hashAlpha alpha + iter es
                | E.Krn (_, dels, dim) => 0w67 + hashDels dels + hashInt dim
                | E.OField(ofld, e2, alpha) => 0w71 +hash' e2  + hash' alpha
                | E.Poly(e1, n1, alpha2) => 0w73 + hash' e1 + hashInt n1 + hashAlpha alpha2
                | E.If(comp, e3, e4) => 0w79+hash' e3 + hash' e4
                | E.Sum(c,e1) => 0w83 + hash' e1
                | E.Op1(e1,e2) => (case e1
                     of E.Cosine => 0w89 + hash' e2
                      | E.ArcCosine => 0w97 + hash' e2
                      | E.Sine => 0w101 + hash' e2
                      | E.ArcSine => 0w103 + hash' e2
                      | E.Tangent => 0w107 + hash' e2
                      | E.ArcTangent => 0w109 + hash' e2
                      | E.Neg => 0w113 + hash' e2
                      | E.Sqrt => 0w127 + hash' e2
                      | E.PowInt _ => 0w131 + hash' e2
                      | E.Exp => 0w137 + hash' e2
                      | E.Abs => 0w139 + hash' e2
                      | E.Sgn => 0w149 + hash' e2
                    (* end case *))
                | E.Op2(E.Sub, e1, e2) => 0w151 + hash' e1 + hash' e2
                | E.Op2(E.Div, e1, e2) => 0w157 + hash' e1 + hash' e2
                | E.Op2(E.Max, e1, e2) => 0w163 + hash' e1 + hash' e2
                | E.Op2(E.Min, e1, e2) => 0w167 + hash' e1 + hash' e2
                | E.Op3(E.Clamp, e1, e2, e3) => 0w173 + hash' e1 + hash' e2 + hash' e3
                | E.Opn(E.Add, es) => 0w179 + iter es
                | E.Opn(E.Prod, es) => 0w181 + iter es
		| E.Fem(femEin1,_,_,_,shape1,dxes1) =>
		  0w191 + hashAlpha shape1 + hashAlpha dxes1
		  + (case femEin1
		      of E.Plain(bda, n, SOME(s, s')) =>
			 0w193 + BasisDataArray.hash bda + Stamp.hash s + Stamp.hash s'
		       | E.Plain(bda, _, _) => 0w193 + BasisDataArray.hash bda
		       | E.Invert(bda,_,SOME(s,s'),w) => 0w197 + BasisDataArray.hash bda + Stamp.hash s + Stamp.hash s' + (if w
															   then 0w227
															   else 0w0)
		       | E.Invert(bda,_,NONE, w) => 0w197 + BasisDataArray.hash bda + (if w
										       then 0w227
										       else 0w0)
		    (*end case*))
		| E.Identity(d, mu, NONE) => 0w199 + hashAlpha [mu] + (Word.fromInt d)
		| E.Identity(d, mu, SOME(b, d')) => 0w199 + hashAlpha [mu] + (Word.fromInt d) + (Word.fromInt d') + (if b
														     then 0w211
														     else 0w223)
              (* end case *)
            end
        in
          hash' body
        end

    fun iterPP es = let
          fun iterP ([], [r]) = r
           | iterP ([], rest) = E.Opn(E.Prod, rest)
           | iterP (E.Const 0::es, rest) = E.Const(0)
           | iterP (E.Const 1::es, rest) = iterP(es, rest)
           | iterP (E.Delta(E.C c1, E.V v1)::E.Delta(E.C c2, E.V v2)::es, rest) =
             (* variable can't be 0 and 1 '*)
               if (c1 = c2)
                 then iterP (es, E.Delta(E.C c1, E.V v1)::E.Delta(E.C c2, E.V v2)::rest)
                 else E.Const(0)
           | iterP (E.Opn(E.Prod, ys)::es, rest) = iterP(ys@es, rest)
           | iterP (e1::es, rest)   = iterP(es, e1::rest)
         in
            iterP (es, [])
          end

    fun iterAA es = let
         fun iterA ([], []) = E.Const 0
           | iterA ([], [r]) = r
           | iterA ([], rest) = E.Opn(E.Add, rest)
           | iterA (E.Const 0::es, rest) = iterA(es, rest)
           | iterA (E.Opn(E.Add, ys)::es, rest) = iterA(ys@es, rest)
           | iterA (e1::es, rest) = iterA(es, e1::rest)
         in
            iterA (es, [])
          end



    fun analyzeInNodes(exp : E.ein_exp, f: E.ein_exp -> 'a option, disamg : 'a * 'a -> 'a) =
	let
	 fun combine(n,[]) = n
	   | combine (n, a::aa) = (case n
				    of NONE => combine(mapper a, aa)
				     | SOME(t) => (case mapper a
						    of NONE => combine(n, aa)
						     | SOME(t') => combine(SOME(disamg(t, t')), aa)))
	 and mapper (e as E.Const _) = NONE
	   | mapper (e as E.ConstR _ ) = NONE
	   | mapper (e as E.Tensor _) = NONE
	   | mapper (e as E.Zero _) = NONE
	   | mapper (e as E.Delta _) = NONE
	   | mapper (e as E.Epsilon _) = NONE
	   | mapper (e as E.Eps2 _) = NONE
	   | mapper (e as E.Field(i, alpha)) = raise Fail "field arguments impossible post arg substitution" 
	   | mapper (E.Lift(e)) = mapper e (*Question: can't only tensor exps be in here... probably safe*)
	   | mapper (e as E.Identity _) = f e
	   | mapper (e as E.Conv _) = f e
	   | mapper (e as E.Fem _ ) = f e
	   | mapper (e as E.Partial _) = raise Fail "can't map through Partials"
	   | mapper (e as E.Apply _) = raise Fail "can't mapp through apply partials"
	   | mapper (E.Comp(e1, e2s)) = let val (last, binds):: rest = List.rev e2s
					in mapper last end
	   | mapper (E.Probe(e1,e2)) = mapper e2
	   | mapper (E.OField _) = raise Fail "impossible: disallow OField"
	   | mapper (e as E.Value _) = NONE
	   | mapper (e as E.Img _) = NONE
	   | mapper (e as E.Krn _) = NONE
	   | mapper (e as E.Poly _) = raise Fail "impossible: disallow Poly"
	   | mapper (E.If(c, e1, e2)) = combine(NONE, [e1,e2])
	   | mapper (E.Sum(sx, e')) =  mapper e'
	   | mapper (E.Op1(u, e')) =  mapper e'
	   | mapper (E.Op2(u, e1, e2)) = combine(NONE, [e1,e2])
	   | mapper (E.Op3(u, e1, e2, e3)) = combine(NONE, [e1,e2, e3])
	   | mapper (E.Opn(u, es)) = combine(NONE, es)
	in
	 mapper exp
	end

    fun mapInNodes(exp : E.ein_exp, f : E.ein_exp -> E.ein_exp) =
	let
	 fun mapper (e as E.Const _) = e
	   | mapper (e as E.ConstR _ ) = e
	   | mapper (e as E.Tensor _) = e
	   | mapper (e as E.Zero _) = e
	   | mapper (e as E.Delta _) = e
	   | mapper (e as E.Epsilon _) = e
	   | mapper (e as E.Eps2 _) = e
	   | mapper (e as E.Field(i, alpha)) = raise Fail "field arguments impossible post arg substitution" 
	   | mapper (E.Lift(e)) = E.Lift(mapper e) (*Question: can't only tensor exps be in here... probably safe*)
	   | mapper (e as E.Identity _) = f e
	   | mapper (e as E.Conv _) = f e
	   | mapper (e as E.Fem _ ) = f e
	   | mapper (e as E.Partial _) = raise Fail "can't map through Partials"
	   | mapper (e as E.Apply _) = raise Fail "can't mapp through apply partials"
	   | mapper (E.Comp(e1, e2s)) = let val (last, binds):: rest = List.rev e2s
					  val e2s' = List.rev ((mapper last, binds) :: rest)
				      in E.Comp(e1, e2s') end
	   | mapper (E.Probe(e1,e2)) = E.Probe(e1, mapper e2)
	   | mapper (E.OField _) = raise Fail "impossible: disallow OField"
	   | mapper (e as E.Value _) = e
	   | mapper (e as E.Img _) = e
	   | mapper (e as E.Krn _) = e
	   | mapper (e as E.Poly _) = raise Fail "impossible: disallow Poly"
	   | mapper (E.If(c, e1, e2)) = E.If(c, mapper e1, mapper e2)
	   | mapper (E.Sum(sx, e')) = E.Sum(sx, mapper e')
	   | mapper (E.Op1(u, e')) = E.Op1(u, mapper e')
	   | mapper (E.Op2(u, e1, e2)) = E.Op2(u, mapper e1, mapper e2)
	   | mapper (E.Op3(u, e1, e2, e3)) = E.Op3(u, mapper e1, mapper e2, mapper e3)
	   | mapper (E.Opn(u, es)) = E.Opn(u, List.map mapper es)
	   (* | mapper _ = raise Fail "impossible" *)
				       
	in
	 mapper exp
	end

    fun detectInNodes(e, f) =
	let
	 fun mapper (e as E.Const _) = false
	   | mapper (e as E.ConstR _ ) = false
	   | mapper (e as E.Tensor _) = false
	   | mapper (e as E.Zero _) = false
	   | mapper (e as E.Delta _) = false
	   | mapper (e as E.Epsilon _) = false
	   | mapper (e as E.Eps2 _) = false
	   | mapper (e as E.Field(i, alpha)) = raise Fail "field arguments impossible post arg substitution" 
	   | mapper (E.Lift(e)) = mapper e (*Question: can't only tensor exps be in here... probably safe*)
	   | mapper (e as E.Identity _) = f e
	   | mapper (e as E.Conv _) = f e
	   | mapper (e as E.Fem _ ) = f e
	   | mapper (e as E.Partial _) = raise Fail "can't map through Partials"
	   | mapper (e as E.Apply _) = raise Fail "can't mapp through apply partials"
	   | mapper (E.Comp(e1, e2s)) = let val (last, binds):: rest = List.rev e2s
					in mapper last end
	   | mapper (E.Probe(e1,e2)) = mapper e2
	   | mapper (E.OField _) = raise Fail "impossible: disallow OField"
	   | mapper (e as E.Value _) = false
	   | mapper (e as E.Img _) = false
	   | mapper (e as E.Krn _) = false
	   | mapper (e as E.Poly _) = raise Fail "impossible: disallow Poly"
	   | mapper (E.If(c, e1, e2)) = mapper e1 orelse mapper e2
	   | mapper (E.Sum(sx, e')) =  mapper e'
	   | mapper (E.Op1(u, e')) =  mapper e'
	   | mapper (E.Op2(u, e1, e2)) = mapper e1 orelse mapper e2
	   | mapper (E.Op3(u, e1, e2, e3)) = mapper e1 orelse mapper e2 orelse mapper e3
	   | mapper (E.Opn(u, es)) = List.foldr (fn (x,y) => mapper x orelse y) false es
	in
	 mapper e
	end
    fun cleanIds(e) =
	let
	 fun mapper (e as E.Const _) = e
	   | mapper (e as E.ConstR _ ) = e
	   | mapper (e as E.Tensor _) = e
	   | mapper (e as E.Zero _) = e
	   | mapper (e as E.Delta _) = e
	   | mapper (e as E.Epsilon _) = e
	   | mapper (e as E.Eps2 _) = e
	   | mapper (e as E.Field(i, alpha)) = raise Fail "field arguments impossible post arg substitution" 
	   | mapper (E.Lift(e)) = E.Lift(mapper e) (*Question: can't only tensor exps be in here... probably safe*)
	   | mapper (e as E.Identity _) = e
	   | mapper (e as E.Conv _) = e
	   | mapper (e as E.Fem _ ) =  e
	   | mapper (e as E.Partial _) = raise Fail "can't map through Partials"
	   | mapper (e as E.Apply _) = raise Fail "can't mapp through apply partials"
	   | mapper (E.Comp(e1, e2s)) =
	     let
	      fun filterFn((E.Identity(_, E.V _, _), _)) = false
		| filterFn _ = true

	      val e2s' = List.filter filterFn (List.map (fn (x, y) => (mapper x, y)) e2s)

	     in
	      if (Bool.not o filterFn) (e1, [])
	      then (case e2s'
		     of [] => e1
		      | (e,_)::[] => e
		      | (e,_)::es => E.Comp(e, es)
		   (* end case*))
	      else (case e2s'
		     of [] => e1
		      | _ => E.Comp(e1, e2s')
		   (* end case*))
	     end
	   | mapper (E.Probe(e1,e2)) = E.Probe(mapper e1, mapper e2)
	   | mapper (E.OField _) = raise Fail "impossible: disallow OField"
	   | mapper (e as E.Value _) = e
	   | mapper (e as E.Img _) = e
	   | mapper (e as E.Krn _) = e
	   | mapper (e as E.Poly _) = raise Fail "impossible: disallow Poly"
	   | mapper (E.If(c, e1, e2)) = (E.If(c, mapper e1, mapper e2))
	   | mapper (E.Sum(sx, e')) =  E.Sum(sx, mapper e')
	   | mapper (E.Op1(u, e')) =  E.Op1(u, mapper e')
	   | mapper (E.Op2(u, e1, e2)) = (E.Op2(u, mapper e1, mapper e2))
	   | mapper (E.Op3(u, e1, e2, e3)) =(E.Op3(u, mapper e1, mapper e2, mapper e3))
	   | mapper (E.Opn(u, es)) = (E.Opn(u, List.map mapper es))
	in
	 mapper e
	end

    structure ISet = IntRedBlackSet
    fun collectIndicies(e) =
	let
	 val start = ISet.empty
	 fun addIdx (E.V e, set) = ISet.add(set, e)
	   | addIdx (E.C _, set) = set
	 fun addAlpha(set, alpha) = List.foldl addIdx set alpha
	 fun mapper (E.Const _, set) = set
	   | mapper (E.ConstR _, set) = set
	   | mapper (E.Tensor(_, alpha), set) = addAlpha(set, alpha)
	   | mapper (E.Zero(alpha), set) = addAlpha(set, alpha)
	   | mapper (E.Delta(a,b), set) = addAlpha(set, [a, b])
	   | mapper (E.Epsilon(a,b,c), set) = addAlpha(set, [a, b, c])
	   | mapper (E.Eps2(a,b), set) = addAlpha(set, [a, b])
	   | mapper (E.Field(_, alpha), set) = addAlpha(set, alpha)
	   | mapper (E.Lift e, set) = mapper(e, set)
	   | mapper (E.Identity(_, a, _), set) = addAlpha(set, [a])
	   | mapper (E.Conv(_, alpha, _, dx), set) = addAlpha(set, List.@(alpha, dx))
	   | mapper (E.Fem(_, _, _, _, alpha, dx), set) = addAlpha(set, List.@(alpha, dx))
	   | mapper (E.Partial(alpha), set) = addAlpha(set, alpha)
	   | mapper (E.Apply(e1, e2), set) = mapper(e2,  mapper(e1, set))
	   | mapper (E.Comp(e, es), set) =
	     let
	      val set1 = mapper(e, set)
	      val (es', bindsIGUess) = ListPair.unzip es (*CHECK ME: go over the use of this bind in compfloat*)
	      val set' = List.foldr (mapper) set1 es'
	     in
	      set'
	     end
	   | mapper (E.Probe(e1, e2), set) = mapper(e2, mapper(e1, set))
	   | mapper (E.OField _, set) = raise Fail "OField disallowed"
	   | mapper (E.Value _, set) = set
	   | mapper (E.Img(_, alpha, _, _), set) = raise Fail "unexpected Mid-IL term"
	   | mapper (E.Krn _, set) = raise Fail "unexpected Mid-IL term"
	   | mapper (E.Poly _, set) = raise Fail "Poly disallowed"
	   | mapper (E.If(_, e1, e2), set) = List.foldr mapper set [e1, e2]
	   | mapper (E.Sum(sx, e), set) =
	     let
	      val sx' = List.map (fn (x, _, _) => x ) sx
	      val set' = ISet.addList(set, sx')
	     in
	      mapper(e, set)
	     end
	   | mapper (E.Op1(_, e), set) = mapper(e, set)
	   | mapper (E.Op2(_, e1, e2), set) = List.foldr mapper set [e1, e2]
	   | mapper (E.Op3(_, e1, e2, e3), set) = List.foldr mapper set [e1, e2, e3]
	   | mapper (E.Opn(_, es), set) = List.foldl mapper set es
					  
	in
	 ISet.listItems(mapper(e, start))
	end

    (* fun typeCheckTenShapeExp (exp, f) = *)
    (* 	(case exp *)
    (* 	  of *)
    (* 	(* end case*)) *)
    (* fun typeCheckTenShape (E.EIN{params, index, body}) = *)
    (* 	let *)
    (* 	 val p = List.length params *)
    (* 	 fun f idx = if idx < p orelse p < 0 *)
    (* 		     then NONE *)
    (* 		     else (case List.nth(params, idx) *)
    (* 			    of E.TEN(_, shape, iv) => SOME((shape, iv)) *)
    (* 			     | E.FLD j => raise Fail "normalized ein has a FLD" *)
    (* 			     | E.KRN => raise Fail "E.Krn checked" *)
    (* 			     | E.IMG(j, shp) => raise Fail "E.IMG checked" *)
    (* 			     | E.FEM(d) => raise Fail "E.FEM checked" *)
    (* 			     | E.INT => raise Fail "Int checked" *)
    (* 			  (* end case *)) *)

    (* 	 val (shape, interval) = typeCheckTenShapeExp(body, f) *)
    (* 	 val sameShape = ListPair.all (fn (x,y) => x=y) (shape, index) *)
    (* 	in *)
    (* 	 ((shape, interval), sameShape) *)
    (* 	end *)
  end
