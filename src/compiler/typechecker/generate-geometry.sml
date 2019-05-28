structure PT = ParseTree
structure Ty = Types
structure TU = TypeUtil
structure E = Env
structure L = Literal
structure SS = Substring
structure CO = CheckOverload
structure FT = FemData
structure FO = FemOpt
structure FN = FemName
structure N = BasisNames
structure CF = CheckTypeFile
structure BV = BasisVars

structure FemGeometry : sig
	   val makeGeometryFuncs : Env.t * (Error.err_stream * Error.span) * Error.span * FemData.femType * 
				   CheckTypeFile.dimensionalObjects list * 
				   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) *
				   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) *
				   (AST.expr list -> AST.expr)
				   -> {cell:(Atom.atom * Var.t) list * AST.global_dcl list * 
					    (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list,
				       pos:(Atom.atom * Var.t) list * AST.global_dcl list * 
					   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list,
				       ref:(Atom.atom * Var.t) list * AST.global_dcl list * 
					   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list}
	  end = struct


val verbosity = false;

fun makePrinStatement(msg, vars, endMsg) = AST.S_Print((AST.E_Lit(Literal.String(msg)))::(vars@[AST.E_Lit(Literal.String(endMsg))]))
						      


fun makePrim'(var, args, argTys, resultTy) =
    let
     val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf var)
     val count1 = List.length domTy
     val count2 = List.length argTys
     val count3 = List.length args
     val _ = if count1 <> count2 orelse count1 <> count3
		orelse count2 <> count3
	     then raise Fail "makePrim!"
	     else ()
     val temp = Unify.matchArgs(domTy, args, argTys) handle exn => raise exn
     val _ = if Option.isSome(temp)
	     then ()
	     else print((String.concatWith "," (List.map TU.toString argTys))^"\n");
     val args' = Option.valOf(temp) 
     val _ = Unify.matchType(resultTy, rngTy)
    in
     AST.E_Prim(var, tyArgs, args', rngTy)
     handle exn => raise exn
    end
    handle exn => raise exn

fun makeAnd v1 v2 = makePrim'(BV.and_b, [v1,v2], [Ty.T_Bool, Ty.T_Bool], Ty.T_Bool)
fun makeAnds([v1,v2]) = makeAnd v1 v2
  | makeAnds (v::vs) = makeAnd v (makeAnds vs)


fun makeGeometryFuncs(env, cxt, span, meshData, geometry, inverse, forwardInfo, makeRefCellInsideFunc) = let
 (*Make types and associated data*)
 val FT.Mesh(mesh) = meshData
 val cellData = FT.MeshCell(mesh)
 val posData = FT.MeshPos(mesh)
 val refData = FT.RefCell(mesh)

 val meshName = FT.nameOf meshData

 val cellTy = Ty.T_Fem(cellData, SOME(meshName))
 val posTy = Ty.T_Fem(posData, SOME(meshName))
 val refTy = Ty.T_Fem(refData, SOME(meshName))
 val meshTy = Ty.T_Fem(meshData, NONE)   
 val dim = FT.meshDim mesh
 val (_, makeTransformFunc, _) = forwardInfo


 val invTransformFieldTy = Ty.T_Field{
      diff = Ty.DiffConst(NONE),
      dim = Ty.DimConst(dim),
      shape = Ty.Shape([Ty.DimConst(dim)])
     }
 val transformFieldTy = invTransformFieldTy

 val dTransformFieldTy = Ty.T_Field{
      diff = Ty.DiffConst(NONE),
      dim = Ty.DimConst(dim),
      shape = Ty.Shape([Ty.DimConst(dim),Ty.DimConst(dim)])
     }						 
 val vecTy = Ty.vecTy dim
 val vec2Ty = Ty.vecTy 2
 val matTy = Ty.matTy dim
 (*for the case where we need to make 4d vecs and mats to do transforms*)
 val vecTyExtra = Ty.vecTy (dim + 1) 
 val matTyExtra = Ty.matTy (dim + 1);

 val vec2SeqTy = Ty.T_Sequence(Ty.realTy, SOME(Ty.DimConst 2))

 val invertBasisVar = if dim = 2
		      then BV.fn_inv2_t
		      else if dim = 3
		      then BV.fn_inv3_t
		      else raise Fail "impossible"
				 
 fun makeRealExpr x = AST.E_Lit(Literal.Real(CF.realToRealLit x))
			       
 fun makeInvTransformFunc([cellExp], eval) =
     let
      val transform = makeTransformFunc([cellExp])
      val dTransform = makePrim'(BV.op_Dotimes, [transform], [transformFieldTy], dTransformFieldTy)
				
      val zero = AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(RealLit.zero true))), vecTy)

      val A = makePrim'(BV.op_probe, [dTransform, zero], [dTransformFieldTy, vecTy], matTy)
      val invA = makePrim'(invertBasisVar, [A], [matTy], matTy) (*FUCK*)
      val offset = makePrim'(BV.op_probe, [transform, zero], [transformFieldTy, vecTy], vecTy)
			    
      val offseted = makePrim'(BV.sub_tt, [eval, offset], [vecTy, vecTy], vecTy)
      val result = makePrim'(BV.op_inner_tt, [invA, offseted], [matTy, vecTy], vecTy)
     in
      (invA, result)
     end
   | makeInvTransformFunc (_) = raise Fail "impossible"			      
 (*if dim=2, line-line intersection; if dim=3 - line-plane intersection*)
 fun rayAtT(refPosExp, dPosExp, t) =
     let
      val scale = makePrim'(BV.mul_tr, [dPosExp, t], [vecTy, Ty.realTy], vecTy)
      val result = makePrim'(BV.add_tt, [refPosExp, scale], [vecTy, vecTy], vecTy)
     in
      result
     end

  fun buildInsideBool(refPosExp, dPosExp, test1) =
      let
       (*comptue new pos *)
       val adjustedTest = test1
       (* val scale = makePrim'(BV.mul_tr, [dPosExp, test1], [vecTy, Ty.realTy], vecTy) *)
       val newPos = rayAtT(refPosExp, dPosExp, adjustedTest)
       val insideTest = makeRefCellInsideFunc([newPos, newPos]);
       (*avoid the need to pass refPos because under our scheme we current forget it...*)
      in
       insideTest
      end
	
 fun kernAtT(normal, dScalar, refPosExp, dPosExp, t) =
     let
      val ray = rayAtT(refPosExp, dPosExp, t)
      val norDot = makePrim'(BV.op_inner_tt, [normal, ray], [vecTy, vecTy], Ty.realTy)
      val result = makePrim'(BV.sub_tt, [norDot, dScalar], [Ty.realTy, Ty.realTy], Ty.realTy)
     in
      result
     end

 fun gradKernAtT(normal, dScalar, refPosExp, dPosExp, t) =
     let
      val data = makePrim'(BV.op_inner_tt, [normal, dPosExp], [vecTy, vecTy], Ty.realTy)
     in
      data
     end

 fun newtonUpdate(normal, dScalar, refPosExp, dPosExp, t) =
     let
      val derv = gradKernAtT(normal, dScalar, refPosExp, dPosExp, t)
      val kern = kernAtT(normal, dScalar, refPosExp, dPosExp, t)
      val update = makePrim'(BV.div_rr, [kern, derv], [Ty.realTy, Ty.realTy], Ty.realTy)
      val updated = makePrim'(BV.sub_tt, [t, update], [Ty.realTy, Ty.realTy], Ty.realTy)
     in
      updated
     end
fun newtonUpdate'(normal, dScalar, refPosExp, dPosExp, t) =
     let
      val derv = gradKernAtT(normal, dScalar, refPosExp, dPosExp, t)
      val kern = kernAtT(normal, dScalar, refPosExp, dPosExp, t)
      val update = makePrim'(BV.div_rr, [kern, derv], [Ty.realTy, Ty.realTy], Ty.realTy)
      val updated = makePrim'(BV.sub_tt, [t, update], [Ty.realTy, Ty.realTy], Ty.realTy)
     in
      (updated, update)
     end

fun newtonLoopBlock(normal, dScalar, refPosExp, dPosExp, maxN, eps, t) =
    let
     val timeVar = Var.new(Atom.atom "t", span, Var.LocalVar, Ty.realTy)
     val timeExp = AST.E_Var(timeVar, span)
     val start = AST.S_Assign((timeVar, span), t)
			     
     val itterVar = Var.new(Atom.atom "i", span, Var.IterVar, Ty.T_Int)
     val itterExp = makePrim'(BV.range, [AST.E_Lit(Literal.intLit 0), maxN], [Ty.T_Int, Ty.T_Int],
			      Ty.T_Sequence(Ty.T_Int, NONE))
     val iter = (itterVar, itterExp)
     val (updated, update) = newtonUpdate'(normal, dScalar, refPosExp, dPosExp, timeExp)
     val absUpdate = makePrim'(BV.op_norm_t, [update], [Ty.realTy], Ty.realTy)
     val smallUpdate = makePrim'(BV.gte_rr, [eps, absUpdate], [Ty.realTy, Ty.realTy], Ty.T_Bool)
				
     val ifStatement = AST.S_IfThenElse(smallUpdate, AST.S_Return(timeExp), AST.S_Block([]))
     val forLoop = AST.S_Block([AST.S_Foreach(iter, AST.S_Block([
						   AST.S_Assign((timeVar,span), updated),
						   ifStatement])),
				AST.S_Return(AST.E_Lit(Literal.Real(RealLit.m_one)))])
    in
     forLoop
    end

 fun newtonUpdates(normal, dScalar, refPosExp, dPosExp, t, n) =
     let
      fun nUpdate(t0) = newtonUpdate(normal, dScalar, refPosExp, dPosExp, t0)
      fun nUp(t0, 0) = t0
	| nUp (t0, n0) = let val up = nUpdate(t0)
			 in nUp(up, n0 - 1)
			 end
     in
      nUp(t, n)
     end

 fun planePointDistSgn(normal, dScalar, refPosExp) =
     let
      val norDot = makePrim'(BV.op_inner_tt, [normal, refPosExp], [vecTy, vecTy], Ty.realTy)
      val result = makePrim'(BV.sub_tt, [dScalar, norDot], [Ty.realTy, Ty.realTy], Ty.realTy)
			    
     in
      result
     end
 fun planePointDist(normal, dScalar, refPosExp) =
     let
      val result = planePointDistSgn(normal, dScalar, refPosExp)
      val abs = makePrim'(BV.op_norm_t, [result], [Ty.realTy], Ty.realTy)
     in
      (abs, result)
     end				      
 fun lineLineIntersect(targetBasis, targetD, refBasis, refD, normal, d) =
     let
      val sub1 = makePrim'(BV.sub_tt, [refBasis, targetBasis], [vecTy, vecTy], vecTy)
      val cross1 = makePrim'(BV.op_cross2_tt, [targetD, refD], [vecTy, vecTy], Ty.realTy)
      val div1 = makePrim'(BV.div_tr, [refD, cross1], [vecTy, Ty.realTy], vecTy)
      val cross2 = makePrim'(BV.op_cross2_tt, [sub1, div1], [vecTy, vecTy], Ty.realTy)
      fun errorFunc t = planePointDistSgn(normal, d, rayAtT(refBasis, refD, t))
      fun refine(n,t) = newtonUpdates(normal, d, refBasis, refD, t, n)
			    
     in
      (cross2, cross1, buildInsideBool(targetBasis, targetD, cross2), refine, errorFunc) (*the first is the time and the second is in case we need to check for nan problems*)
     end

 fun linePlaneIntersect(refBasis, refD, normal, dScalar) =
     let
      val dot1 = makePrim'(BV.op_inner_tt, [normal, refBasis], [vecTy, vecTy], Ty.realTy)
      val num = makePrim'(BV.sub_tt, [dScalar, dot1], [Ty.realTy, Ty.realTy], Ty.realTy)
      val dot2 = makePrim'(BV.op_inner_tt, [normal, refD], [vecTy, vecTy], Ty.realTy)
      val result = makePrim'(BV.div_rr, [num, dot2], [Ty.realTy, Ty.realTy], Ty.realTy)
      fun errorFunc t = planePointDistSgn(normal, dScalar, rayAtT(refBasis, refD, t))
			    
      fun refine(n, t) = newtonUpdates(normal, dScalar, refBasis, refD, t, n)

     in
      (result, dot2, buildInsideBool(refBasis, refD, result), refine, errorFunc)
     end




 fun listPairCheck(v1, v2) = if List.length v1 = List.length v2
			     then ()
			     else raise Fail "pairCheck!"
       
 fun buildIntersectionTestInfo(2, geometry, refPosExp, dRefPosExp) =
     let
      fun lineIntersect(bma, a, n, d) = lineLineIntersect(refPosExp, dRefPosExp, a, bma, n, d)
      val lineParams = Option.valOf (List.find (fn CF.LineParam(xs) => true | _ => false) geometry)
		       handle exn => raise exn
      val CF.LineParam(lineData) = lineParams
      (*convert to realLits - thankfully only vectors*)

      val vecExprs = List.map (fn (x,y,(z1,z2)) => (AST.E_Tensor(List.map makeRealExpr x, vecTy),
						    AST.E_Tensor(List.map makeRealExpr y, vecTy))) lineData;
      val normalDVecs = List.map (fn (_, _, (z1, z2)) => (AST.E_Tensor(List.map makeRealExpr z1, vecTy), makeRealExpr z2)) lineData
      val _ = listPairCheck(vecExprs, normalDVecs)
      val combinedExprs = ListPair.map (fn ((x,y), (a,b)) => (x,y,a,b)) (vecExprs, normalDVecs) handle exn => raise exn
				 (*time, parallel param*)
      val intersectionExprs = List.map lineIntersect combinedExprs
      fun dist(n, d) = planePointDist(n, d, refPosExp)
				     
      (*dist from refPos to facet, distance signed*)
      val intersectionExprs' = List.map dist normalDVecs
      val _ = listPairCheck(intersectionExprs, intersectionExprs')
      val resultingTests = ListPair.map (fn ((x, y, j, r, f), (w, z)) => (x, y, j, r, f, w, z)) (intersectionExprs, intersectionExprs') handle exn => raise exn
					
				       
     in
      resultingTests handle exn => raise exn
     end
   | buildIntersectionTestInfo(3, geometry, refPosExp, dRefPosExp ) =
     let
      fun planeIntersect(d, normal) = linePlaneIntersect(refPosExp, dRefPosExp, normal, d)
      fun dist(d, n) = planePointDist(n, d, refPosExp)
      val planeParams = Option.valOf (List.find (fn CF.PlaneParam(_) => true | _ => false) geometry) handle exn => raise exn
      val CF.PlaneParam(xs, _) = planeParams
      val planeParamExprs = List.map (fn (x, ys) => (makeRealExpr x, AST.E_Tensor(List.map makeRealExpr ys, vecTy))) xs
      val intersectionExprs = List.map planeIntersect planeParamExprs
      val distTestExprs = List.map dist planeParamExprs
      val _ = listPairCheck(intersectionExprs, distTestExprs)
      val resultingTests = ListPair.map (fn ((x, y, j, r, f), (w, z)) => (x, y, j, r, f, w, z)) (intersectionExprs, distTestExprs)
			   handle exn => raise exn
				       
     in
      resultingTests handle exn => raise exn
     end
										      


 fun facetPrint(test1, test2, insideTestval, dist, signedDist, intExp, (r,d), error, off, dp, maybeInt, otherPossibleTime) =
     let
      val place = rayAtT(r, d, test1)
      fun mkStr x = AST.E_Lit(Literal.String(x))
      val extraPrints = (case maybeInt
			  of NONE => []
			   | SOME(fint) => [mkStr "\n original face: ", fint, mkStr " "])
     in
      AST.S_Print(List.@([mkStr "\ntest :", test1, mkStr " test2: ", test2, mkStr " insideTest: ", insideTestval,
			  mkStr " dist: ", dist, mkStr "signed dist: ", signedDist,
			  mkStr " face: ", intExp, mkStr " place: ", place,
		   mkStr " error: ", error,
		   mkStr " offset :", off,
		   mkStr " dp: ", dp,
		   mkStr " otherPossibleTime: ", otherPossibleTime,
		   mkStr "\n" ],extraPrints))
     end

 fun intersectionTesting(intersectionExprs, parallelTest, parallelTestNeps,
			 insideTest, preRefine, postRefine, printFlag, facetIntTest, max, (r, d)) =
     let
      fun coerce e = AST.E_Coerce({srcTy=Ty.T_Int, dstTy=Ty.realTy, e= e})
      val tempVar = Var.new (Atom.atom "time", span, AST.LocalVar, Ty.realTy)
      val faceReserveVar = Var.new (Atom.atom "faceTime", span, AST.LocalVar, Ty.realTy)
      val tempVar' = Var.new (Atom.atom "face", span, AST.LocalVar, Ty.T_Int)
      val tempExp = AST.E_Var(tempVar, span)
      val tempExp' = AST.E_Var(tempVar', span)
      val faceExp = AST.E_Var(faceReserveVar, span)
      val neg1 = AST.E_Lit(Literal.Int(IntLit.fromInt (~1)))
      val faceStart = AST.S_Decl(faceReserveVar, SOME(AST.E_Lit(Literal.Real(RealLit.negInf))))
      val tempStart = AST.S_Decl(tempVar, SOME(AST.E_Lit(Literal.Real(max))))

      val tempStart' = AST.S_Decl(tempVar', SOME(neg1))
				 
      val zero = AST.E_Lit(Literal.Real(RealLit.zero false))
      val eps = AST.E_Lit(Literal.Real(RealLit.fromDigits{isNeg = false, digits = [1], exp = IntInf.fromInt (~parallelTestNeps)}))
      val zeroEps = AST.E_Lit(Literal.Real(RealLit.fromDigits{isNeg = true, digits = [1], exp = IntInf.fromInt (~15)}))

      val timeEps = AST.E_Lit(Literal.Real(RealLit.fromDigits{isNeg = true, digits = [1], exp = IntInf.fromInt (~07)}))
      val normalEps = AST.E_Lit(Literal.Real(RealLit.fromDigits{isNeg = false, digits = [1], exp = IntInf.fromInt (~07)}))
      val distEps = AST.E_Lit(Literal.Real(RealLit.fromDigits{isNeg = false, digits = [1], exp = IntInf.fromInt (~07)}))

      fun build(isCloseToZero, eps) =
	  let
	   val absTest2 = makePrim'(BV.op_norm_t, [isCloseToZero], [Ty.realTy], Ty.realTy)
	   val antiNanInfTest = makePrim'(BV.gte_rr, [absTest2, eps], [Ty.realTy, Ty.realTy], Ty.T_Bool)
					 
	  in
	   antiNanInfTest
	  end
			 
      fun buildIf(test1, test2, insideTestVal, refine, errorAtT, dist, sgnDist, intExpr) =
	  let
	   (*compute special test, int*)
	   val testError = errorAtT test1
	   val absTestError = makePrim'(BV.op_norm_t, [testError], [Ty.realTy], Ty.realTy)
				       
	   val preRefinedTest1 = refine(preRefine, test1)

	   val posRefineTest2 = preRefinedTest1 (*refine(postRefine, preRefinedTest1)*)
	   val absTest2 = makePrim'(BV.op_norm_t, [test2], [Ty.realTy], Ty.realTy)
				   
	   val positiveTest = makePrim'(BV.gte_rr, [preRefinedTest1, timeEps], [Ty.realTy, Ty.realTy], Ty.T_Bool)
	   val newUpdateTest = makePrim'(BV.gt_rr, [tempExp, preRefinedTest1], [Ty.realTy, Ty.realTy], Ty.T_Bool)
	   val antiNanInfTest = makePrim'(BV.gte_rr, [absTest2, normalEps], [Ty.realTy, Ty.realTy], Ty.T_Bool)
	   val savePrint = if printFlag
			   then [makePrinStatement("Saving at this face!", [preRefinedTest1, newUpdateTest], "\n")]
			   else []
	   (*Checks if a stored entry facet is equal to this facet; prevents failure if the pos is slightly outside the ref*)
	   val optionalIntTest = (case facetIntTest
				   of SOME(oldFacetInt) => (fn b =>
							       let
								val test = makePrim'(BV.neq_ii, [oldFacetInt, intExpr], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
								val printBackup = if printFlag
									    then [makePrinStatement("Saving backup time: ", [faceExp], "\n")]
									    else []
							       in
								AST.S_IfThenElse(test, b, AST.S_Block([AST.S_Assign((faceReserveVar, span), test1)]@printBackup))
							       end)
				   |  NONE => (fn b => b))
	   val preTests = if parallelTest
			  then [positiveTest, newUpdateTest, antiNanInfTest]
			  else [positiveTest, newUpdateTest]
	   val test = makeAnds(preTests)
	   val proc = if insideTest
		      then (fn b => AST.S_IfThenElse(insideTestVal, b, AST.S_Block([])))
		      else (fn b => b)


	   (* makePrinStatement("Suc with:", [test1,AST.E_Lit(Literal.String(", ")), test2], "\n"), *)
	   val fin = AST.S_IfThenElse(test, (optionalIntTest o proc)
					      (AST.S_Block(savePrint@[
							   AST.S_Assign((tempVar, span), posRefineTest2),
							   AST.S_Assign((tempVar', span), intExpr)
					      ])),
				      AST.S_Block([]))
	   val print1 = facetPrint(preRefinedTest1, test2, insideTestVal, dist, sgnDist, intExpr, (r, d), absTestError, r, d, facetIntTest, faceExp)

	  in
	   if printFlag
	   then	AST.S_Block([print1, fin])
	   else fin
	  end
      val tests = List.length intersectionExprs
      val intLits = List.tabulate(tests, fn x => AST.E_Lit(Literal.intLit x))
      val _ = listPairCheck(intersectionExprs, intLits)
      val zip = ListPair.map (fn ((x,y,a,b,c,d,e), z) => (x,y,a,b,c,d,e, z)) (intersectionExprs, intLits)
		handle exn => raise exn
      val onPlaneVals = List.map (fn  (x,y,a,b,c,d,e, z) => y) zip
			     
      val ifs = List.map buildIf zip
      val nanTests = List.map (fn x => build(x, eps)) onPlaneVals
      val nanTestsWithFace = ListPair.zip (List.tabulate(List.length nanTests, fn x => x), nanTests)


      val timeReturn = makePrim'(BV.fn_max_r, [tempExp, zero], [Ty.realTy, Ty.realTy], Ty.realTy)
      val backupTimeReturn = makePrim'(BV.fn_max_r, [faceExp, zero], [Ty.realTy, Ty.realTy], Ty.realTy)
      (*build base + time*dpos       val end = rayAtT(r, d, timeReturn)
*)
      val workedTestFacet = makePrim'(BV.neq_ii, [tempExp', neg1], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
      val workedTestInside = buildInsideBool(r, d, timeReturn) (*use intersection time*)
      val workedTest = makeAnds([workedTestFacet, workedTestInside])

      (*make actual fail*)



      val failBlock = AST.S_Block([ (*TODO: should this be ~1? What is the exact standard for this function?*)
						  AST.S_Return(AST.E_Tensor([coerce neg1, coerce neg1], vec2Ty))
				 ])
      val failRet = AST.S_Return(AST.E_Tensor([coerce neg1, coerce neg1], vec2Ty))
      fun wrapReturn(res as AST.S_Return(e), msg) = if printFlag
						    then AST.S_Block([makePrinStatement("\nReturning via " ^ msg ^ " : ", [e], "\n"), res])
						    else res
      fun makeOnPlaneBackup([]) = failBlock
	| makeOnPlaneBackup ((i, nanTest)::rest) = failBlock
	  (* let *)
	  (*  val thisInt = AST.E_Lit(Literal.Int(IntLit.fromInt i)) *)
	  (*  val thisCellTest = (case facetIntTest *)
	  (* 			of SOME(umm) => makePrim'(BV.neq_ii, [umm, thisInt], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)) *)
	  (*  val test =  makeAnds([nanTest, thisCellTest]) *)
	  (*  val ret = AST.S_Return(AST.E_Tensor([zero, coerce thisInt], vec2Ty)) *)
	   
	  (* in *)
	  (*  AST.S_IfThenElse(test, ret, makeOnPlaneBackup rest) *)
      (* end *)


      val failBlock' = (case facetIntTest
			 of NONE => failBlock
			  | SOME(faceInt) =>
			    let
			     val antiNanFaceBackUp = makeOnPlaneBackup nanTestsWithFace
			     val facetTest = makePrim'(BV.equ_rr, [AST.E_Lit(Literal.Real(RealLit.negInf)), faceExp], [Ty.realTy, Ty.realTy], Ty.T_Bool)
			     val otherReturn = AST.S_Return(AST.E_Tensor([backupTimeReturn, coerce faceInt], vec2Ty))
			    in
			     AST.S_IfThenElse(facetTest, wrapReturn(failRet, "fail"), wrapReturn(otherReturn, "backup"))
			    end) (*TODO: This conditional is not needed... if facet is never tested I think...*)

      val ifReturn = AST.S_IfThenElse(workedTest,
				      AST.S_Block([
						  wrapReturn(AST.S_Return(AST.E_Tensor([timeReturn, coerce tempExp'], vec2Ty)), "standard")
						 ]),
				      failBlock'
				     )
      (*check if we are on a facet and die just in case*)

      val stms = [tempStart, tempStart',faceStart]@ifs@[ifReturn]

     in
      stms handle exn => raise exn
     end

			 
			 
 local
  (*_exit and exit; _enter and enter*)
  val refPosParam = Var.new(Atom.atom "refPos", span, AST.FunParam, vecTy)
  val dposParam = Var.new(Atom.atom "dPos", span, AST.FunParam, vecTy)
  val intPosParam = Var.new(Atom.atom "i", span, AST.FunParam, Ty.T_Int)
  val refPosExp = AST.E_Var(refPosParam, span)
  val dPosExp = AST.E_Var(dposParam, span)
  val intPosExp = AST.E_Var(intPosParam, span)
  val funType = Ty.T_Fun([vecTy, vecTy, Ty.T_Int], vec2Ty)
  val funAtom = Atom.atom "_exit"
  val funVar = Var.new (funAtom, span, Var.FunVar, funType)
		       
  val tests = buildIntersectionTestInfo(dim, geometry, refPosExp, dPosExp)


  val body = AST.S_Block(intersectionTesting(tests, true, 20, false, 0, 0, verbosity, SOME(intPosExp), RealLit.posInf, (refPosExp, dPosExp)))
  val result = ((funAtom, funVar), AST.D_Func(funVar, [refPosParam, dposParam, intPosParam], body))

		 
  val funAtom' = Atom.atom "_enter"
  val funType' = Ty.T_Fun([vecTy, vecTy], vec2Ty)
  val funVar' = Var.new (funAtom', span, Var.FunVar, funType')
	
  val body' = AST.S_Block(intersectionTesting(tests, true, 20, true, 0, 0, verbosity, NONE, RealLit.posInf, (refPosExp, dPosExp)))
  val result' = ((funAtom', funVar'), AST.D_Func(funVar', [refPosParam, dposParam], body'))


  val exitFuncTy = Ty.T_Fun([refTy, posTy, vecTy], Ty.realTy)
  val exitFuncName = Atom.atom (FemName.refExit)
  fun replaceExit ([re, pos, vec]) =
      let
       (*get refPos, call the function, etc*)
       val refPos = AST.E_ExtractFemItem(pos, vecTy, (FemOpt.RefPos, posData))
       val posEntryFacet = AST.E_ExtractFemItem(pos, Ty.T_Int, (FemOpt.PosEntryFacet, posData))
       val res = AST.E_Apply((funVar, span), [refPos, vec, posEntryFacet], vec2Ty)
       val ret = AST.E_Slice(res, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy )
      in
       ret

      end
  val funcResult = (exitFuncName, replaceExit, exitFuncTy)

  val enterFuncName = Atom.atom (FemName.refEnter)
  val enterFuncTy = Ty.T_Fun([refTy, vecTy, vecTy], Ty.realTy)
  fun replaceEnter ([re, vec1, vec2]) =
      let
       val res = AST.E_Apply((funVar', span), [vec1, vec2], vec2Ty)
       val ret = AST.E_Slice(res, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy )
      in
       ret
      end
    | replaceEnter (_) = raise Fail "typechecker error"

  val enterFuncResult = (enterFuncName, replaceEnter, enterFuncTy)

 in
 val hiddenExitFuncResult = result
 val hiddenEnterFuncResult = result'
 val enterFuncResult = enterFuncResult
 val exitFuncResult = funcResult
 end
 (*Several versions of cell connectivity:
		 1. Raw
		 2. Nearby cells only -> search through them
		 3. (cell, face) -> cell
		 4. (cell, face) -> (cell', face')
		For now, we have the 2nd and 4th one. 
		This is where you index into cells and the into face -> an int2 pair.
 
  *)
 fun faceConnectivityToCellFact(meshExp, cellExp, facetExp) =
     let
      val FT.Mesh(mesh) = meshData
      val FemData.RefCellData{numFaces, ...} = FemData.refCell mesh
      val cellOffset = AST.E_Lit(Literal.intLit (2 * numFaces))

      val cellIndexIn = makePrim'(BV.mul_ii, [cellExp, cellOffset], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
      val two = AST.E_Lit(Literal.intLit 2)
      val facetIndexIn = makePrim'(BV.mul_ii, [facetExp, two], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
      val offset = makePrim'(BV.add_ii, [cellIndexIn, facetIndexIn], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
			    
      val neg1 = AST.E_Lit(Literal.intLit (~1))
      val validCellTest = makePrim'(BV.equ_ii, [cellExp, neg1], [Ty.T_Int, Ty.T_Int],Ty.T_Bool)
      val validFacet = makePrim'(BV.equ_ii, [facetExp, neg1], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
      val combinedTest = makePrim'(BV.and_b, [validCellTest, validFacet], [Ty.T_Bool, Ty.T_Bool], Ty.T_Bool)
      val opt = (FemOpt.CellFaceCell, meshData)
      val retTy = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2))
      val femExp = AST.E_ExtractFemItemN([meshExp, offset], [meshTy, Ty.T_Int], retTy, opt, NONE)
      val printRet = makePrinStatement("temp:", [femExp, cellExp, facetExp], "\n");
      val ret = AST.S_Return(femExp)
      val badReturn = AST.S_Return(AST.E_Seq([neg1, neg1], retTy));

      val test = AST.S_IfThenElse(validFacet, badReturn, ret)
      (*Could insert test here*)
     in
      AST.S_Block(test::[])
     end

 local
  val facetIntParam = Var.new(Atom.atom "faceIdx", span, AST.FunParam, Ty.T_Int)
  val cellParam = Var.new (Atom.atom "cell", span, AST.FunParam, Ty.T_Int)
  val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)
  val params = [facetIntParam, cellParam, meshParam]
  val facetIntExp = AST.E_Var(facetIntParam, span)
  val cellExp = AST.E_Var(cellParam, span);
  val meshExp = AST.E_Var(meshParam, span);
  (* val buildNextCellFunction = buildIndexAnalysis(dim, geometry) *)
  (* val functionBody = buildNextCellFunction(meshExp, meshData, cellExp, facetIntExp) *)
  (* val nextCellAtom = Atom.atom "$nextCell2" *)
  (* val nextCellTy = Ty.T_Fun([Ty.T_Int, Ty.T_Int, meshTy], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2))) *)
  (* val nextCellFuncVar = Var.new (nextCellAtom, span, AST.FunVar, nextCellTy) *)
  (* val nextCellFunc = AST.D_Func(nextCellFuncVar, params, AST.S_Block(functionBody)) *)


  (* fun internalCellFuncReplace([seqExp, meshPosExp]) = *)
  (*     let *)
  (*      (*extract meshPos*) *)
  (*      val meshExp = AST.E_ExtractFem(meshPosExp, meshData) *)
  (*      val cellExpr = AST.E_ExtractFemItem(meshPosExp,Ty.T_Int, (FemOpt.CellIndex, meshData)) *)
  (*      val facetIdExpr = makePrim'(BV.floor, [AST.E_Slice(seqExp, [SOME(AST.E_Lit(Literal.intLit 1))], Ty.realTy)], [Ty.realTy], Ty.T_Int) *)
  (*     in *)
  (*      AST.E_Apply((nextCellFuncVar, span), [facetIdExpr, cellExpr, meshExp], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2))) *)
  (*     end *)

  val nextCellAtom4 = Atom.atom "nextCell4"
  val nextCellTy4 = Ty.T_Fun([Ty.T_Int, Ty.T_Int, meshTy], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))
  val nextCellFuncVar4 = Var.new (nextCellAtom4, span, AST.FunVar, nextCellTy4)
  val functionBody4 = faceConnectivityToCellFact(meshExp, cellExp, facetIntExp)
  val nextCellFunc4 = AST.D_Func(nextCellFuncVar4, params, AST.S_Block([functionBody4]))

 in
 (* val cellFunc = ((nextCellAtom,nextCellFuncVar), nextCellFunc) *)
 val cellFunc4 = ((nextCellAtom4, nextCellFuncVar4), nextCellFunc4)
 (* val internalCellFuncReplace = internalCellFuncReplace *)
 end

 (*let's think about validity here and what not 
		  - calling on an invalid cell means everything breaks; we do have an option of an invalid facet-then no need to do either of the above
-		  - should be caught here *)
 (*call first function - args go to this one - if facet valid, call next one - get new one, extract tensor, translate, return *)
 (*wait - if facet is invalid, it means soething went wrong in intersection.... - means in boundary*)
 (*on boundary and coordinates problem...*)
 (*correct - is-inside issue*)
 (*invalid facet is probably impossible because a line has to intersect one of them... assuming no nans or bs like that...put the if condition just in case and comment it in later.
		 also need to check if planes could be invalid - i.e if all three points are co-line - either sub-det is 0...*)
 (*do build new pos outside*)
 (*function: if invalidFacet 1 or invalidFacet 2... return invalid build...  *)
 (*inside valid -> dot or division -> to solve based on new pos -> *)

 (*build a function body; build replace that coordinates everything else.*)
 fun buildSolveBlock(dim, srcFacetExp, dstFacetExpr, meshExp, newCellExp, refPosExp, geometry) =
     let
      (*Check the src,dst,and cell for valid; if yes, we do this; otherwise, return an invalid one*)
      fun printSolveInfo(bigmatrix, face1, face2, selected) =
	  let
	   fun mkStr x = AST.E_Lit(Literal.String(x))
	  in
	   AST.S_Print([
		       mkStr "face ", face1, mkStr " to ", face2, mkStr "\n",
		       mkStr "got: ", selected, mkStr "\n"
		      ])
	  end
      
      fun buildSolveOperation(CF.LineParam(xs)) =
	  let
	   val correctIndexes = List.map (fn ([a,b],y,z) => if Real.<=(Real.abs(b),0.000001)
							    then AST.E_Lit(Literal.intLit 0)
							    else AST.E_Lit(Literal.intLit 1)) xs
	   val count = List.length correctIndexes
	   val selectTy = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst count))
	   val seqVecTy = Ty.T_Sequence(Ty.realTy, SOME(Ty.DimConst(2)));
	   val selectTy' = Ty.T_Sequence(seqVecTy, SOME(Ty.DimConst count))
	   val selectTy'' = Ty.T_Sequence(vecTy, SOME(Ty.DimConst count))

	   val (betaSeqs, betaTensor) = ListPair.unzip (List.map (fn (xs, _, _) => let val x = List.map makeRealExpr xs in 
										    (AST.E_Seq(x, seqVecTy), AST.E_Tensor(x, vecTy)) end ) xs)
	   val betas = AST.E_Seq(betaSeqs, selectTy')
	   val betasTensors = AST.E_Seq(betaTensor, selectTy'')
	   val (alphaSeqs, alphaTensor) = ListPair.unzip (List.map (fn (_, xs, _) => let val x = List.map makeRealExpr xs in 
										      (AST.E_Seq(x, seqVecTy), AST.E_Tensor(x, vecTy)) end ) xs)
	   val alphas = AST.E_Seq(alphaSeqs, selectTy')
	   val alphasTensor = AST.E_Seq(alphaTensor, selectTy'')
				       
	   val solveIndex = makePrim'(BV.subscript, [AST.E_Seq(correctIndexes, selectTy), srcFacetExp], [selectTy, Ty.T_Int], Ty.T_Int)
	   (*a_dst[i] + (b-a)_dst[i] t = refPosExp[i]*)
	   (*alpha + tbeta = gamma -> gamma-alpha/beta*)
	   val (alphaSeq, betaSeq, targetSeq) =
	       (
		 makePrim'(BV.subscript, [alphas, srcFacetExp], [selectTy', Ty.T_Int], seqVecTy),
		 makePrim'(BV.subscript, [betas, srcFacetExp], [selectTy', Ty.T_Int], seqVecTy),
		 AST.E_Seq([AST.E_Slice(refPosExp, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy),
			    AST.E_Slice(refPosExp, [SOME(AST.E_Lit(Literal.intLit 1))], Ty.realTy)], seqVecTy))
		 
	   val (alpha, beta, target) =
	       (
		 makePrim'(BV.subscript, [alphaSeq, solveIndex], [seqVecTy, Ty.T_Int], Ty.realTy),
		 makePrim'(BV.subscript, [betaSeq, solveIndex], [seqVecTy, Ty.T_Int], Ty.realTy),
		 makePrim'(BV.subscript, [targetSeq, solveIndex], [seqVecTy, Ty.T_Int], Ty.realTy)
	       )

	   val sub1 = makePrim'(BV.sub_tt, [target, alpha], [Ty.realTy, Ty.realTy], Ty.realTy)
	   val timeAlongLine = makePrim'(BV.div_rr, [sub1, beta], [Ty.realTy, Ty.realTy], Ty.realTy)

	   val (newAlpha, newBeta) = (makePrim'(BV.subscript, [alphasTensor, dstFacetExpr], [selectTy'', Ty.T_Int], vecTy),
				      makePrim'(BV.subscript, [betasTensors, dstFacetExpr], [selectTy'', Ty.T_Int], vecTy))
	   val newRefPos = makePrim'(BV.add_tt, [newAlpha, makePrim'(BV.mul_tr, [newBeta, timeAlongLine], [vecTy, Ty.realTy], vecTy)],
				     [vecTy, vecTy], vecTy)
	  in
	   (*meshExp,cellExpr, posExpr*)
	   (NONE, AST.E_ExtractFemItemN([meshExp, newCellExp, newRefPos, dstFacetExpr], [meshTy, Ty.T_Int, vecTy, Ty.T_Int], posTy, (FemOpt.RefBuild, meshData), NONE))
	  end
	| buildSolveOperation(CF.PlaneParam(_, xs)) =
	  let
	   val mat4 = Ty.matTy 4
	   val vec4 = Ty.vecTy 4

	   fun buildTensor(a) = AST.E_Tensor(List.map (fn x => AST.E_Tensor(List.map makeRealExpr x , vec4)) a, mat4) handle exn => raise exn
	   val count = List.length xs
	   val seqTy1 = Ty.T_Sequence(mat4, (SOME(Ty.DimConst count)))
	   val seqTy2 = Ty.T_Sequence(seqTy1, SOME(Ty.DimConst count))

	   fun buildSeq(xs') = AST.E_Seq(List.map (fn x => AST.E_Seq(List.map buildTensor x, seqTy1)) xs', seqTy2)	
	   val seqExp = buildSeq xs

	   val selectedSeq = makePrim'(BV.subscript, [seqExp, srcFacetExp], [seqTy2, Ty.T_Int], seqTy1) handle exn => raise exn

	   val selectedTensor = makePrim'(BV.subscript, [selectedSeq, dstFacetExpr], [seqTy1, Ty.T_Int], mat4) handle exn => raise exn

	   val vec4refPos = AST.E_Tensor(List.tabulate(3, fn x => AST.E_Slice(refPosExp, [SOME(AST.E_Lit(Literal.intLit x))], Ty.realTy))@[AST.E_Lit(Literal.Real(RealLit.one))],
					 vec4)
					
	   val resultVec4 = makePrim'(BV.op_inner_tt, [selectedTensor, vec4refPos], [mat4, vec4], vec4) handle exn => raise exn
	   val newRefPos = AST.E_Tensor(List.tabulate(3, fn x => AST.E_Slice(resultVec4, [SOME(AST.E_Lit(Literal.intLit x))], Ty.realTy)), vecTy)
	   val result =  AST.E_ExtractFemItemN([meshExp, newCellExp, newRefPos, dstFacetExpr], [meshTy, Ty.T_Int, vecTy, Ty.T_Int], posTy, (FemOpt.RefBuild, meshData), NONE)
	   val printStm = printSolveInfo(seqExp, srcFacetExp, dstFacetExpr, selectedTensor )

	  in
	   (SOME(AST.S_Block([])), result) handle exn => raise exn
	  end


      val geometryPass = if dim=2
			 then Option.valOf (List.find (fn CF.LineParam(_) => true | _ => false) geometry)
			 else if dim = 3
			 then Option.valOf (List.find (fn CF.PlaneParam(_) => true | _ => false) geometry) handle exn => raise exn
			 else raise Fail "ops"
      val solveCall = buildSolveOperation(geometryPass) handle exn => raise exn

     in
      solveCall
      handle exn => raise exn
     end
 (*note exposing this...for the world based ones we should be doing something else.... also, we could take the time argument!*)
 fun makeRefExitPosBody(meshExp, cellExp, refPosExp, dPosExp, timeAndFaceExp, dim, geometry, debug) =
     let
      (*first, build new expression*)
      val time = AST.E_Slice(timeAndFaceExp, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy)
      val srcFacet = makePrim'(BV.floor, [AST.E_Slice(timeAndFaceExp, [SOME(AST.E_Lit(Literal.intLit 1))], Ty.realTy)], [Ty.realTy], Ty.T_Int)
      val srcFacetPrint = makePrinStatement("SrcFacet: ", [srcFacet], "\n");
      val newVec = makePrim'(BV.add_tt, [makePrim'(BV.mul_rt, [time, dPosExp], [Ty.realTy, vecTy], vecTy), refPosExp], [vecTy, vecTy], vecTy)
      val ((_, nextCellFuncVar),_) = cellFunc4
      val nextCellAndFace = AST.E_Apply((nextCellFuncVar, span), [srcFacet, cellExp, meshExp], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))
      val newCell = makePrim'(BV.subscript, [nextCellAndFace, AST.E_Lit(Literal.intLit 0)], [Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)), Ty.T_Int], Ty.T_Int)
      val dstFace = makePrim'(BV.subscript, [nextCellAndFace, AST.E_Lit(Literal.intLit 1)], [Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)), Ty.T_Int], Ty.T_Int)
      val dstFacetPrint = makePrinStatement("dstFacet: ", [dstFace], "\n");
      val (printOption, solveBody) = buildSolveBlock(dim, srcFacet, dstFace, meshExp, newCell, newVec, geometry)
      (* srcFacetPrint,dstFacetPrint, *)
      val return = (case printOption
		     of SOME(pr) => AST.S_Block([pr, AST.S_Return(solveBody)])
		      | NONE => AST.S_Block([AST.S_Return(solveBody)])
		    (*end case*))
      val boundaryMeshPos = AST.E_ExtractFemItemN([meshExp, newVec], [meshTy, vecTy], posTy, (FemOpt.InvalidBuildBoundary, posData), NONE)
      val failReturnMesh = AST.S_Return(boundaryMeshPos)
      val condition = makePrim'(BV.neq_ii, [AST.E_Lit(Literal.intLit (~1)), newCell], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
      val stm = AST.S_IfThenElse(condition, return, failReturnMesh);
      (*TODO:This code does not check if the time is invalid; assumption is won't be called without time for now.*)
     in
      stm
     end
 local
  (*build complete function that depends on nothing:*)
  (**)
  val ((_, exitFuncVar), _) = hiddenExitFuncResult

  val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)
  val cellParam = Var.new (Atom.atom "cellId", span, AST.FunParam, Ty.T_Int)
  val refPosParam = Var.new (Atom.atom "refPos", span, AST.FunParam, vecTy)
  val dPosParam = Var.new (Atom.atom "dPos", span, AST.FunParam, vecTy)
  val timeAndFace = Var.new(Atom.atom "time", span, AST.FunParam, vec2SeqTy)

  fun mkVar v = AST.E_Var(v,span)
			 
  val meshExp = mkVar meshParam
  val cellExp = mkVar cellParam
  val refExp = mkVar refPosParam
  val dExp = mkVar dPosParam
  val tfExp = mkVar timeAndFace;
  
  val catchFacetIssue = true
  val hiddenExitPosAtom = Atom.atom "$exitPos"
  val hiddenExitPosTy = Ty.T_Fun([meshTy, Ty.T_Int, vecTy, vecTy, vec2SeqTy], posTy) (*call this in build world*)
  val hiddenExitPosVar = Var.new (hiddenExitPosAtom, span, Var.FunVar, hiddenExitPosTy)
  val bodyStm = makeRefExitPosBody(meshExp, cellExp, refExp, dExp, tfExp, dim, geometry, catchFacetIssue)
  val hiddeExitPosFunc = AST.D_Func(hiddenExitPosVar, [meshParam, cellParam, refPosParam, dPosParam, timeAndFace], bodyStm)
  val hiddenExitFuncPosResult = ((hiddenExitPosAtom, hiddenExitPosVar), hiddeExitPosFunc)

  val exitPos = Atom.atom (FemName.refExitPos)
  val exitPosTy = Ty.T_Fun([refTy, posTy, vecTy], posTy)
  fun exitPosReplace'([re, pos, vec]) =
      let
       val mesh = AST.E_ExtractFem(pos, meshData)
       val cell = AST.E_ExtractFemItem(pos, Ty.T_Int, (FemOpt.CellIndex, meshData))
       val refPos = AST.E_ExtractFemItem(pos, vecTy, (FemOpt.RefPos, meshData))
       val posEntryFacet = AST.E_ExtractFemItem(pos, Ty.T_Int, (FemOpt.PosEntryFacet, posData))
       val seq = AST.E_Apply((exitFuncVar, span), [refPos, vec, posEntryFacet], vec2SeqTy) (*ugg*)
       (*apply func*)
       val result = AST.E_Apply((hiddenExitPosVar,span), [mesh, cell, refPos, vec, seq], posTy)
      in
       (result, seq)
      end
    | exitPosReplace'(_) = raise Fail "impossible got past type checking."
  fun exitPosReplace(vs) = (fn (x,y) => x) (exitPosReplace'(vs))
  val refExitPosReplace = (exitPos, exitPosReplace, exitPosTy)

			    
			    
 in
 val seqRefPosExitFunc = exitPosReplace'
 val refPosExitFunc = hiddenExitFuncPosResult
 val refExitPosReplace = refExitPosReplace
 end

 (*build transforms between reference cells...*)
 (*modify others to produce correct division scheme.*)
 (*one function that selects transform and does it.*)
 (*emebds all results in the building of mesh pos - call first function, give results to second function, give results to third, which builds and gets inlined*)

 (*meshPos/meshCell based ones:*)
 (*Transform*)
 (* val transform = AST.E_FemField(meshExp, meshExp, SOME(cellIntExp), transformFieldTy, NONE) *)
 (* val dTransform = makePrim'(BV.op_Dotimes, [transform], [transformFieldTy], dTransformFieldTy) *)
 
 (*meshPos.exit/meshCell.exit*)
 (*issue: have world vs. have reference... make a toggle...*)
 fun buildWorldIntersectsExp (true, exitEnter, meshExp, cellExp, cellIntExp, refVec1Exp, vec1Exp, vec2Exp, facetInt) =
     let
      val transformField = makeTransformFunc([cellExp])
      val vec1Exp = if Option.isSome(refVec1Exp)
      		    then let val SOME(new') = refVec1Exp in new' end
      		    else vec1Exp
			   
      val (dInvAtZero, inverted) = makeInvTransformFunc([cellExp], vec1Exp)
      val trans1 =  if Option.isSome(refVec1Exp)
 		    then vec1Exp
 		    else inverted
			   
      val derv = dInvAtZero 
      val trans2 = makePrim'(BV.op_inner_tt, [derv, vec2Exp], [matTy, vecTy], vecTy) handle exn => raise exn
      val appFunc  = if exitEnter
		     then let val ((_, enterFuncVar), _) = hiddenEnterFuncResult
			  in enterFuncVar end
		     else let val ((_, exitFuncVar), _) = hiddenExitFuncResult
			  in exitFuncVar end
      val args = if exitEnter
		 then  [trans1, trans2]
		 else [trans1, trans2, facetInt]
			    
      val seqResult = AST.E_Apply((appFunc, span), args, vec2SeqTy)
      val ret = AST.E_Slice(seqResult, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy)
			   (*other exit call...*)
			   
     in
      (ret, (seqResult, trans1, trans2))
     end
   | buildWorldIntersectsExp (false, exitEnter, meshExp, cellExp, cellIntExp, refVec1Exp, vec1Exp, vec2Exp, facetInt) = raise Fail "later"



														    
 local
  (*build _exit, build exit func, enter func*)
  val FT.Mesh(mesh) = meshData
  val cellHiddenExit = Atom.atom "$cellExit" (*hidden*)
  val cellHiddenExitTy = Ty.T_Fun([meshTy, cellTy, Ty.T_Int, vecTy, vecTy, vecTy], vec2SeqTy)
  fun cellHiddenExitFun(b,v1, v2, v3, v4, v5, v6, v7) = buildWorldIntersectsExp(FemData.isAffine mesh,b, v1, v2, v3, v4, v5, v6, v7)


									   

  val posExit = Atom.atom (FemName.posExit)
  val posExitTy = Ty.T_Fun([posTy, vecTy], Ty.realTy)
  fun posExitFun([v1, v2]) =
      let
       val mesh = AST.E_ExtractFem(v1, meshData)
       val cellInt = AST.E_ExtractFemItem(v1, Ty.T_Int, (FO.CellIndex, posData))
       val cellExp =  AST.E_LoadFem(cellData, SOME(mesh), SOME(cellInt))
       val refPos = AST.E_ExtractFemItem(v1, vecTy, (FemOpt.RefPos, meshData))
       val transformField = makeTransformFunc([cellExp])
       val posEntryFacet = AST.E_ExtractFemItem(v1, Ty.T_Int, (FemOpt.PosEntryFacet, posData))
       val worldPos =  makePrim'(BV.op_probe, [transformField, refPos], [transformFieldTy, vecTy], vecTy) handle exn => raise exn (*TODO: fix worldPos*)
       val (result, _) = cellHiddenExitFun(false, mesh, cellExp, cellInt, SOME(refPos), worldPos, v2, posEntryFacet)

      in
       result
      end


  val cellEnter = Atom.atom (FemName.cellEnter)
  val cellEnterTy = Ty.T_Fun([cellTy, vecTy, vecTy], Ty.realTy)
  fun cellEnterFun([v1, v2, v3]) =
      let
       val mesh = AST.E_ExtractFem(v1, meshData)
       val cellExp = v1
       val cellInt = AST.E_ExtractFemItem(v1, Ty.T_Int, (FO.CellIndex, cellData))
       val refVecExp = NONE (*transform from v2*)
       val (result, _) = cellHiddenExitFun(true, mesh, cellExp, cellInt, NONE, v2, v3, AST.E_Lit(Literal.intLit (~1)))
      in
       result
      end



  val cellEnterPos = Atom.atom (FemName.cellEnterPos)
  val cellEnterPosTy = Ty.T_Fun([cellTy, vecTy, vecTy], posTy)
  fun cellEnterPosFun([v1, v2, v3]) =
      let
       val mesh = AST.E_ExtractFem(v1, meshData)
       val cellExp = v1
       val cellInt = AST.E_ExtractFemItem(v1, Ty.T_Int, (FO.CellIndex, cellData))
       val refVecExp = NONE (*transform from v2*)
       val (timeResult, (seq, refBPos, refDPos)) = cellHiddenExitFun(true, mesh, cellExp, cellInt, NONE, v2, v3, AST.E_Lit(Literal.intLit (~1)))
       val refPos = rayAtT(refBPos, refDPos, timeResult)
       val worldPos = rayAtT(v2, v3, timeResult)
       val facetIdExpr = makePrim'(BV.floor, [AST.E_Slice(seq, [SOME(AST.E_Lit(Literal.intLit 1))], Ty.realTy)], [Ty.realTy], Ty.T_Int)
       val posResult = AST.E_ExtractFemItemN([mesh, cellInt, refPos, worldPos, facetIdExpr], [meshTy, Ty.T_Int, vecTy, vecTy, Ty.T_Int], posTy, (FemOpt.AllBuild, posData), NONE)
      in
       posResult
      end
	

	
 in
 val cellHiddenResult = (cellHiddenExit,cellHiddenExitFun, cellHiddenExitTy)
 val posExitResult = (posExit, posExitFun, posExitTy)
 val cellEnterResult = (cellEnter, cellEnterFun, cellEnterTy)
 val cellEnterPosResult = (cellEnterPos, cellEnterPosFun, cellEnterPosTy)
 end

 fun buildWorld(true, posExp, vecExp) =
     let
      (*get the int, ref*)
      val refExitPosFunc = seqRefPosExitFunc
      val mesh = AST.E_ExtractFem(posExp, meshData)
      val cellInt = AST.E_ExtractFemItem(posExp, Ty.T_Int, (FO.CellIndex, posData))
      val cellExp =  AST.E_LoadFem(cellData, SOME(mesh), SOME(cellInt))
      val zeroExp = AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(RealLit.zero true))), vecTy)


      val transformField = makeTransformFunc([cellExp])
      val invertVar = if dim = 2
 		      then BV.fn_inv2_f
 		      else if dim = 3
 		      then BV.fn_inv3_f
 		      else raise Fail "dim <> 2 and dim <> 3"
      val dTransformVal = makePrim'(BV.op_Dotimes, [transformField], [invTransformFieldTy], dTransformFieldTy)
      val invertedField = makePrim'(invertVar, [dTransformVal], [dTransformFieldTy], dTransformFieldTy)
      (* val refPos = AST.E_ExtractFemItem(posExp, vecTy, (FemOpt.RefPos, meshData)) *)
      val probe = makePrim'(BV.op_probe, [invertedField, zeroExp], [dTransformFieldTy, vecTy], matTy)
      (* val (invA, _) = makeInvTransformFunc([cellExp], zeroExp) *)
      val newDpos = makePrim'(BV.op_inner_tt, [probe, vecExp], [matTy, vecTy], vecTy)

			     
      val mesh = AST.E_ExtractFem(cellExp, meshData)
      val (posResult, seqResult) = refExitPosFunc([mesh, posExp, newDpos]) (*TODO: Compute worldPos here too and percolate to simple.*)
     in
      posResult
     end
   | buildWorld (false, posExp, vecExp) =  raise Fail "not implemented"

 local
  val FT.Mesh(mesh) = meshData
  val posExitName = Atom.atom (FemName.posExitPos)
  val posExitTy = Ty.T_Fun([posTy, vecTy], posTy)
  fun posExitFun([v1,v2]) = buildWorld(FemData.isAffine mesh, v1,v2)
    | posExitFun (_) = raise Fail "typechecker error"
			     
 in
 val posExitPosResult = (posExitName, posExitFun, posExitTy)
 end
 (* call the above, send to per cell newton with a wider epsilon control? -- go back to meshCell and that refCell epsilon...*)
 (*other option is to: transform to reference small steps....*)
 (* 0.00001 -> get a new pos out via previous method...*)
 (*My schedule: deal with maps later if we need them - we probably don't for the moment.*)

 
 val refActualFuncs = [hiddenExitFuncResult, hiddenEnterFuncResult, cellFunc4, refPosExitFunc]
 val (refActualFuncsInfo, refActualFuncDcls) = ListPair.unzip refActualFuncs
 val refReplaceFuncs = [exitFuncResult, refExitPosReplace, enterFuncResult]

 val posReplaceFuncs = [posExitResult, posExitPosResult]
 val posActualFuncs = []
 val (posActualFuncInfo, posActualFuncDcls) = ListPair.unzip posActualFuncs
 val results = {ref = (refActualFuncsInfo, refActualFuncDcls, refReplaceFuncs),
		pos = (posActualFuncInfo,posActualFuncDcls,posReplaceFuncs),
		cell = ([],[],[cellEnterResult, cellEnterPosResult])}


in
 results
 handle exn => raise exn
end

end
