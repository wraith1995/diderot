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
fun makePrim(prim, argTys) =
    let
     val AST.E_Prim(var, metargs, exprs, result) = prim
     val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf var)
     val exprs' = Option.valOf(Unify.matchArgs(domTy, exprs, argTys))
     val _ = Unify.matchType(result, rngTy)

    in
     AST.E_Prim(var, tyArgs, exprs', result)
    end

fun makePrinStatement(msg, vars, endMsg) = AST.S_Print((AST.E_Lit(Literal.String(msg)))::(vars@[AST.E_Lit(Literal.String(endMsg))]))
						      


fun makePrim'(var, args, argTys, resultTy) =
    let
     val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf var)
     val temp = Unify.matchArgs(domTy, args, argTys)
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

structure FemGeometry : sig
	   val makeGeometryFuncs : Env.t * (Error.err_stream * Error.span) * Error.span * FemData.femType * 
				   CheckTypeFile.dimensionalObjects list * 
				   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) *
				   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty)  
				   -> {cell:(Atom.atom * Var.t) list * AST.global_dcl list * 
					    (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list,
				       pos:(Atom.atom * Var.t) list * AST.global_dcl list * 
					   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list,
				       ref:(Atom.atom * Var.t) list * AST.global_dcl list * 
					   (Atom.atom * (AST.expr list -> AST.expr) * Types.ty) list}
	  end = struct
fun makeGeometryFuncs(env, cxt, span, meshData, geometry, inverse, forwardInfo) = let
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
 val (_, makeInvTransformFunc, invTransformFieldTy) = inverse
 val (_, makeTransformFunc, transformFieldTy) = forwardInfo
						  

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

 (*if dim=2, line-line intersection; if dim=3 - line-plane intersection*)
 fun lineLineIntersect( targetBasis, targetD, refBasis, refD) =
     let
      val sub1 = makePrim'(BV.sub_tt, [refBasis, targetBasis], [vecTy, vecTy], vecTy)
      val cross1 = makePrim'(BV.op_cross2_tt, [targetD, refD], [vecTy, vecTy], Ty.realTy)
      val div1 = makePrim'(BV.div_tr, [refD, cross1], [vecTy, Ty.realTy], vecTy)
      val cross2 = makePrim'(BV.op_cross2_tt, [sub1, div1], [vecTy, vecTy], Ty.realTy)
			    
     in
      (cross2, cross1) (*the first is the time and the second is in case we need to check for nan problems*)
     end

 fun linePlaneIntersect(refBasis, refD, normal, dScalar) =
     let
      val dot1 = makePrim'(BV.op_inner_tt, [normal, refBasis], [vecTy, vecTy], Ty.realTy)
      val num = makePrim'(BV.sub_tt, [dScalar, dot1], [Ty.realTy, Ty.realTy], Ty.realTy)
      val dot2 = makePrim'(BV.op_inner_tt, [normal, refD], [vecTy, vecTy], Ty.realTy)
      val result = makePrim'(BV.div_rr, [num, dot2], [Ty.realTy, Ty.realTy], Ty.realTy)
			    
     in
      (result, dot2)
     end

 fun makeRealExpr x = AST.E_Lit(Literal.Real(CF.realToRealLit x))
 fun twoDimTests(refPosExp, dPosExp, geometry) =
     let
      fun lineIntersect(bma, a) = lineLineIntersect(refPosExp, dPosExp, a, bma)
      val lineParams = Option.valOf (List.find (fn CF.LineParam(xs) => true | _ => false) geometry) handle exn => raise exn
      val CF.LineParam(lineData) = lineParams
      (*convert to realLits - thankfully only vectors*)

      val vecExprs = List.map (fn (x,y,z) => (AST.E_Tensor(List.map makeRealExpr x, vecTy),
					      AST.E_Tensor(List.map makeRealExpr y, vecTy))) lineData
      val intersectionExprs = List.map lineIntersect vecExprs
     in
      intersectionExprs
     end

 fun threeDimTests(refPosExp, dPosExp, geometry) =
     let
      fun planeIntersect(d, normal) = linePlaneIntersect(refPosExp, dPosExp, normal, d)

      val planeParams = Option.valOf (List.find (fn CF.PlaneParam(_) => true | _ => false) geometry) handle exn => raise exn
      val CF.PlaneParam(xs, _) = planeParams
      val planeParamExprs = List.map (fn (x, ys) => (makeRealExpr x, AST.E_Tensor(List.map makeRealExpr ys, vecTy))) xs
      val intersectionExprs = List.map planeIntersect planeParamExprs
     in
      intersectionExprs handle exn => raise exn
     end
     handle exn => raise exn

 (*create local vars with +Inf, -1; Store compute in local var, if >=0 and <= current, update  *)
 fun intersectionTesting intersectionExprs =
     let
      val tempVar = Var.new (Atom.atom "time", span, AST.LocalVar, Ty.realTy)
      val tempVar' = Var.new (Atom.atom "face", span, AST.LocalVar, Ty.T_Int)
      val tempExp = AST.E_Var(tempVar, span)
      val tempExp' = AST.E_Var(tempVar', span)
      val tempStart = AST.S_Decl(tempVar, SOME(AST.E_Lit(Literal.Real(RealLit.posInf))))
      val neg1 = AST.E_Lit(Literal.Int(IntLit.fromInt (~1)))
      val tempStart' = AST.S_Decl(tempVar', SOME(neg1))
      val zero = AST.E_Lit(Literal.Real(RealLit.zero true))
      fun buildIf(test1, test2, intExpr) =
	  let
	   (*compute special test, int*)
	   val positiveTest = makePrim'(BV.gt_rr, [test1, zero], [Ty.realTy, Ty.realTy], Ty.T_Bool)
	   val newUpdateTest = makePrim'(BV.gt_rr, [tempExp, test1], [Ty.realTy, Ty.realTy], Ty.T_Bool)
	   val combined = makePrim'(BV.and_b, [positiveTest, newUpdateTest], [Ty.T_Bool, Ty.T_Bool], Ty.T_Bool)
	   val fin = AST.S_IfThenElse(combined, AST.S_Block([
							    AST.S_Assign((tempVar, span), test1),
							    AST.S_Assign((tempVar', span), intExpr)
							   ]),
				      AST.S_Block([]))
	  in
	   fin
	  end
	    
      val tests = List.length intersectionExprs
      val intLits = List.tabulate(tests, fn x => AST.E_Lit(Literal.intLit x))
      val zip = ListPair.map (fn ((x,y), z) => (x,y,z)) (intersectionExprs, intLits)
			     
      val ifs = List.map buildIf zip

      val workedTest = makePrim'(BV.neq_ii, [tempExp', neg1], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
      fun coerce e = AST.E_Coerce({srcTy=Ty.T_Int, dstTy=Ty.realTy, e= e})

      val ifReturn = AST.S_IfThenElse(workedTest,
				      AST.S_Block([
						  AST.S_Return(AST.E_Tensor([tempExp, coerce tempExp'], vec2Ty))
						 ]),
				      AST.S_Block([ (*TODO: should this be ~1? What is the exact standard for this function?*)
						  AST.S_Return(AST.E_Tensor([coerce neg1, coerce neg1], vec2Ty))
				     ]))

      val stms = [tempStart, tempStart']@ifs@[ifReturn]

     in
      stms  handle exn => raise exn
     end
     handle exn => raise exn
			 
			 
 local
  (*_exit and exit*)
  val refPosParam = Var.new(Atom.atom "refPos", span, AST.FunParam, vecTy)
  val dposParam = Var.new(Atom.atom "dPos", span, AST.FunParam, vecTy)
  val refPosExp = AST.E_Var(refPosParam, span)
  val dPosExp = AST.E_Var(dposParam, span)
  val funType = Ty.T_Fun([vecTy, vecTy], vec2Ty)
  val funAtom = Atom.atom "_exit"
  val funVar = Var.new (funAtom, span, Var.FunVar, funType)
  val tests = if dim = 2
	      then twoDimTests(refPosExp, dPosExp, geometry)
	      else if dim = 3
	      then threeDimTests(refPosExp, dPosExp, geometry) handle exn => raise exn
	      else
	       raise Fail "Dim ought to be 2 or 3 in check-global.sml; this should not have been called at this point."
  val body = AST.S_Block(intersectionTesting tests)
  val result = ((funAtom, funVar), AST.D_Func(funVar, [refPosParam, dposParam], body))
  (*make type, var*)
  (*build intersections based on dim*)
  (*build stms*)
  (*build function*)

  val exitFuncTy = Ty.T_Fun([refTy, posTy, vecTy], Ty.realTy)
  val exitFuncName = Atom.atom (FemName.refExit)
  fun replaceExit ([re, pos, vec]) =
      let
       (*get refPos, call the function, etc*)
       val refPos = AST.E_ExtractFemItem(pos, vecTy, (FemOpt.RefPos, meshData))
       val res = AST.E_Apply((funVar, span), [refPos, vec], vec2Ty)
       val ret = AST.E_Slice(res, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy )
      in
       ret

      end
  val funcResult = (exitFuncName, replaceExit, exitFuncTy)


  val enterFuncName = Atom.atom (FemName.refEnter)
  val enterFuncTy = Ty.T_Fun([refTy, vecTy, vecTy], Ty.realTy)
  fun replaceEnter ([re, vec1, vec2]) =
      let
       val res = AST.E_Apply((funVar, span), [vec1, vec2], vec2Ty)
       val ret = AST.E_Slice(res, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy )
      in
       ret
      end
    | replaceEnter (_) = raise Fail "typechecker error"

  val enterFuncResult = (enterFuncName, replaceEnter, enterFuncTy)
			  
 in
 val hiddenExitFuncResult = result
 val enterFuncResult = enterFuncResult
 val exitFuncResult = funcResult
 end

 fun buildIndexAnalysis(dim, geometry) =
     let
      val SOME(CF.Points(_,SOME(buildIndex))) = (List.find (fn CF.Points(_) => true | _ => false) geometry) handle exn => raise exn
      val SOME(CF.Higher(_, xs)) = (List.find (fn CF.Higher(dim', _) => dim'=dim -1 | _ => false) geometry) handle exn => raise exn
      fun buildNodeIndex(idx) = List.nth(buildIndex, idx)
      val planeIndecies : int list list = List.map (fn (x,y,z) => List.map (buildNodeIndex o IntInf.toInt) z) xs;


      (*first, build sizing *)
      (*second build giant array*)
      (*build expression index and Known index*)
      val sizes = List.map List.length planeIndecies
      val size = List.foldr (op+) 0 sizes
      val planeCount = List.length sizes
      val sizes' = List.map (fn x => List.take(sizes,x)) (List.tabulate(planeCount, fn x => x));
      val offsets = List.map (List.foldr Int.* 1) sizes'
      val offsetAndSize = ListPair.zip (offsets, sizes)
      fun buildIntLit x = AST.E_Lit(Literal.intLit x)
      val int2vecTy = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2))
      val offsetAndSizeExprTy = Ty.T_Sequence(int2vecTy, SOME(Ty.DimConst planeCount))
      val offsetAndSizeExpr = AST.E_Seq(List.map (fn (x,y) => AST.E_Seq([buildIntLit x, buildIntLit y], int2vecTy)) offsetAndSize, offsetAndSizeExprTy)

      fun flat xs = List.foldr op@ [] xs
      val planeIndeciesTy =  Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst size))
      val planeIndeciesExpr = AST.E_Seq(List.map buildIntLit (flat planeIndecies), planeIndeciesTy)
      (*We need to partion planeSizes to things of the same size*)
      val planeIndeciesAndIndex : (int list * int) list = ListPair.zip (planeIndecies, List.tabulate(planeCount, fn x => x));
      val maxSize = List.foldr (Int.max) 0 sizes
      val planeIndeciesBySize = List.filter (fn (y,x) => not (List.null x)) (List.tabulate(maxSize+1, fn i => (i, List.filter (fn (x,y) => List.length x = i) planeIndeciesAndIndex)) : (int * (int list * int) list) list)
      (*for each size, we have a lits of facets with their indecies into an array of nodes and there integer id*)
      (*{nodes[allIndecies[offset[face][0]+i]] : i in (0, offset[face][1]}*)
					    
      fun getUnknownNodeAccesses faceExpr (nodesExpr, nodesTy) =
	  let
	   val offsetInto = makePrim'(BV.subscript, [offsetAndSizeExpr, faceExpr], [offsetAndSizeExprTy, Ty.T_Int], int2vecTy)
	   val offset = makePrim'(BV.subscript, [offsetInto, AST.E_Lit(Literal.intLit 0)], [int2vecTy, Ty.T_Int], Ty.T_Int)
	   val size = makePrim'(BV.subscript, [offsetInto, AST.E_Lit(Literal.intLit 1)], [int2vecTy, Ty.T_Int], Ty.T_Int)
			       
	   val itterVar = Var.new (Atom.atom "i", span, Var.IterVar, Ty.T_Int)
	   val itterVarRange = makePrim'(BV.range, [AST.E_Lit(Literal.intLit 0), size], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
	   val access1 = makePrim'(BV.add_ii, [offset, AST.E_Var(itterVar, span)],
				   [Ty.T_Int, Ty.T_Int], Ty.T_Int)
	   val access2 = makePrim'(BV.subscript, [planeIndeciesExpr, access1], [planeIndeciesTy, Ty.T_Int], Ty.T_Int)
	   (*todo: can this be sub?*)
	   val fin = makePrim'(BV.dynSubscript, [nodesExpr, access2], [nodesTy, Ty.T_Int], Ty.T_Int)
	   val comprehension = AST.E_Comprehension(access2, (itterVar, itterVarRange), Ty.T_Sequence(Ty.T_Int, NONE))
	  in
	   comprehension
	  end


      fun getKnownNodeAccess faceInt (nodesExpr, nodesTy) =
	  let
	   (*since we know the face in advance, we just grab the nodes from known indecies*)
	   val faceIdxes = List.nth(planeIndecies, faceInt)
	   val faceIdxesCount = List.nth(sizes, faceInt)
	   val faceIntLits = List.map buildIntLit faceIdxes
	   val nodeAcceses = List.map (fn x => makePrim'(BV.subscript, [nodesExpr, x], [nodesTy, Ty.T_Int], Ty.T_Int)) faceIntLits
	   val ty = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst faceIdxesCount))
	   val result = AST.E_Seq(nodeAcceses, ty)
				 
	  in
	   (result, ty)
	  end

      fun comparisionStm count cellIntExp faceInt unknownFaceNodesExp (knownFaceNodesExp, ty) =
	  let

	   (*compare nodes from the first cell to the current cell; we've already loaded both at this point *)
	   
	   fun makeIntEquality(i) = makePrim'(BV.equ_ii,
					      [
						makePrim'(BV.subscript, [knownFaceNodesExp, AST.E_Lit(Literal.intLit i)],
							  [ty, Ty.T_Int], Ty.T_Int),
						makePrim'(BV.dynSubscript, [unknownFaceNodesExp, AST.E_Lit(Literal.intLit i)],
							  [Ty.T_Sequence(Ty.T_Int, NONE), Ty.T_Int], Ty.T_Int)
							 
					      ], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
	   val tests = List.map makeIntEquality (List.tabulate(count, fn x => x))
	   val test = makeAnds(tests);
	   val ifStm = AST.S_IfThenElse(test,
					AST.S_Return(AST.E_Seq([cellIntExp, AST.E_Lit(Literal.intLit faceInt)], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))),
					AST.S_Block([]))
				       (*return new cell and new facet*)
				       

	  in
	   ifStm
	  end

      fun perSizeComparisionLoop(meshExp, meshData, faceSize, cellSeq, faceInts : int list, unknownFaceNodesExp) =
	  (*for each cell in c: nodes <- loadNodes(c) ifstm1, ifstm2,*)
	  let
	   val FT.Mesh(mesh) = meshData
	   val spaceDim = FT.meshMapDim mesh
	   val itterVar = Var.new (Atom.atom "cellInt", span, Var.IterVar, Ty.T_Int)
	   val itterVarExp = AST.E_Var(itterVar, span)
	   val loadMeshIndecies = (AST.E_ExtractFemItem2(meshExp, itterVarExp, Ty.T_Int, Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst spaceDim)), (FemOpt.ExtractIndices, meshData)), Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst spaceDim)))
	   fun buildFaceNodes(i : int)= (getKnownNodeAccess i loadMeshIndecies, i)
	   val knownNodeAcc = List.map buildFaceNodes faceInts
	   fun buildComps (exp, i) = comparisionStm faceSize itterVarExp i unknownFaceNodesExp exp
	   val comps = AST.S_Block(List.map buildComps knownNodeAcc)
	   val loop = AST.S_Foreach((itterVar, cellSeq), comps)
				   
	  in
	   loop
	  end



      fun buildFunction(meshExp, meshData, cellIntExp, faceIntExp) =
	  let
	   (*todo: lift some defs here up?*)
	   (*todo: check if cell and facet are valid.*)
	   val FT.Mesh(mesh) = meshData
	   val spaceDim = FT.meshMapDim mesh
	   (*get nearby cells*)
	   val nearbyCellVar = Var.new (Atom.atom "nearbyCells", span, Var.LocalVar, Ty.T_Sequence(Ty.T_Int, NONE))
	   val nearbyCellVarExp = AST.E_Var(nearbyCellVar, span)
	   val nearbyCellVarAssign = AST.S_Assign((nearbyCellVar,span), AST.E_ExtractFemItem2(meshExp, cellIntExp, Ty.T_Int, Ty.T_Sequence(Ty.T_Int, NONE), (FemOpt.CellConnectivity, meshData)));
	   
	   (*get the relevant nodes*)
	   val nodesTy = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst spaceDim))
	   val generalNodesVar = Var.new (Atom.atom "cellNodes", span, Var.LocalVar, nodesTy)
	   val generalNodesExp = AST.E_Var(generalNodesVar,span)
	   val generalNodesAssign = AST.S_Assign((generalNodesVar,span), AST.E_ExtractFemItem2(meshExp, cellIntExp, Ty.T_Int, nodesTy, (FemOpt.ExtractIndices, meshData)));
	   (*get the face-node sequence*)
	   val faceNodes = getUnknownNodeAccesses faceIntExp (generalNodesExp, nodesTy);
	   val numUnknownKnownsVar = Var.new (Atom.atom "numFacetNodes", span, Var.LocalVar, Ty.T_Int)
	   val numUknownKnownsAssign = AST.S_Assign((numUnknownKnownsVar, span), makePrim'(BV.fn_length, [faceNodes], [Ty.T_Sequence(Ty.T_Int, NONE)], Ty.T_Int))
	   fun makeSizeTest i = makePrim'(BV.equ_ii, [AST.E_Lit(Literal.intLit i), AST.E_Var(numUnknownKnownsVar, span)], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)

	   (*for each size, partition and build comparision loop*)

	   fun buildLoopsByFaceSize(size, faceIntsIntPair) =
	       let
		val ret = perSizeComparisionLoop(meshExp, meshData, size, nearbyCellVarExp, (List.map (fn (x,y) => y) faceIntsIntPair) , faceNodes)
	       in
		AST.S_IfThenElse(makeSizeTest size,
				 ret,
				 AST.S_Block([]))
	       end
	   val loops = List.map buildLoopsByFaceSize planeIndeciesBySize
	   val failRet = AST.S_Return(AST.E_Seq([AST.E_Lit(Literal.intLit (~1)), AST.E_Lit(Literal.intLit (~1))], int2vecTy))				 
				     (*build return at the end.*)
	  in
	   nearbyCellVarAssign::numUknownKnownsAssign::loops@[failRet]
	  end
     in
      buildFunction
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
      (*Could insert test here*)
     in
      AST.S_Block(printRet::ret::[])
     end

 local
  val facetIntParam = Var.new(Atom.atom "faceIdx", span, AST.FunParam, Ty.T_Int)
  val cellParam = Var.new (Atom.atom "cell", span, AST.FunParam, Ty.T_Int)
  val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)
  val params = [facetIntParam, cellParam, meshParam]
  val facetIntExp = AST.E_Var(facetIntParam, span)
  val cellExp = AST.E_Var(cellParam, span);
  val meshExp = AST.E_Var(meshParam, span);
  val buildNextCellFunction = buildIndexAnalysis(dim, geometry)
  val functionBody = buildNextCellFunction(meshExp, meshData, cellExp, facetIntExp)
  val nextCellAtom = Atom.atom "$nextCell2"
  val nextCellTy = Ty.T_Fun([Ty.T_Int, Ty.T_Int, meshTy], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))
  val nextCellFuncVar = Var.new (nextCellAtom, span, AST.FunVar, nextCellTy)
  val nextCellFunc = AST.D_Func(nextCellFuncVar, params, AST.S_Block(functionBody))


  fun internalCellFuncReplace([seqExp, meshPosExp]) =
      let
       (*extract meshPos*)
       val meshExp = AST.E_ExtractFem(meshPosExp, meshData)
       val cellExpr = AST.E_ExtractFemItem(meshPosExp,Ty.T_Int, (FemOpt.CellIndex, meshData))
       val facetIdExpr = makePrim'(BV.floor, [AST.E_Slice(seqExp, [SOME(AST.E_Lit(Literal.intLit 1))], Ty.realTy)], [Ty.realTy], Ty.T_Int)
      in
       AST.E_Apply((nextCellFuncVar, span), [facetIdExpr, cellExpr, meshExp], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))
      end

  val nextCellAtom4 = Atom.atom "nextCell4"
  val nextCellTy4 = Ty.T_Fun([Ty.T_Int, Ty.T_Int, meshTy], Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst 2)))
  val nextCellFuncVar4 = Var.new (nextCellAtom4, span, AST.FunVar, nextCellTy4)
  val functionBody4 =faceConnectivityToCellFact(meshExp, cellExp, facetIntExp)
  val nextCellFunc4 = AST.D_Func(nextCellFuncVar4, params, AST.S_Block([functionBody4]))


				



 in
 val cellFunc = ((nextCellAtom,nextCellFuncVar), nextCellFunc)
 val cellFunc4 = ((nextCellAtom4, nextCellFuncVar4), nextCellFunc4)
 val internalCellFuncReplace = internalCellFuncReplace
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
	   AST.E_ExtractFemItemN([meshExp, newCellExp, newRefPos], [meshTy, Ty.T_Int, vecTy], posTy, (FemOpt.RefBuild, meshData), NONE)
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
	   val result =  AST.E_ExtractFemItemN([meshExp, newCellExp, newRefPos], [meshTy, Ty.T_Int, vecTy], posTy, (FemOpt.RefBuild, meshData), NONE)

	  in
	   result handle exn => raise exn
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
      val solveBody = buildSolveBlock(dim, srcFacet, dstFace, meshExp, newCell, newVec, geometry)
      val return = AST.S_Block([srcFacetPrint,dstFacetPrint, AST.S_Return(solveBody)])
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

  val exitPos = Atom.atom "exitPos"
  val exitPosTy = Ty.T_Fun([refTy, posTy, vecTy], posTy)
  fun exitPosReplace([re, pos, vec]) =
      let
       val mesh = AST.E_ExtractFem(pos, meshData)
       val cell = AST.E_ExtractFemItem(pos, Ty.T_Int, (FemOpt.CellIndex, meshData))
       val refPos = AST.E_ExtractFemItem(pos, vecTy, (FemOpt.RefPos, meshData))
       val seq = AST.E_Apply((exitFuncVar, span), [refPos, vec], vec2SeqTy) (*ugg*)
       (*apply func*)
       val result = AST.E_Apply((hiddenExitPosVar,span), [mesh, cell, refPos, vec, seq], posTy)
      in
       result
      end
    | exitPosReplace(_) = raise Fail "impossible got past type checking."

  val refExitPosReplace = (exitPos, exitPosReplace, exitPosTy)

			    
			    
 in
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
 fun buildWorldIntersectsExp (true, meshExp, cellExp, cellIntExp, refVec1Exp, vec1Exp, vec2Exp) =
     let
      val invTransformField = makeInvTransformFunc([cellExp])
      val transformField = makeTransformFunc([cellExp])
      val vec1Exp = if Option.isSome(refVec1Exp)
      		    then let val SOME(new') = refVec1Exp in new' end
      		    else vec1Exp
			   
      val dInvTransformField = if Option.isSome(refVec1Exp)
 			       then
 				let
 				 val temp = makePrim'(BV.op_Dotimes, [transformField], [invTransformFieldTy], dTransformFieldTy)
 				 val inverse = if dim = 2
 					       then BV.fn_inv2_f
 					       else if dim = 3
 					       then BV.fn_inv3_f
 					       else raise Fail "dim <> 2 and dim <> 3"
							  
 				 val invTemp = makePrim'(inverse, [temp], [dTransformFieldTy], dTransformFieldTy)
							
 				in
 				 invTemp
 				end
 			       else
 				makePrim'(BV.op_Dotimes, [invTransformField], [invTransformFieldTy], dTransformFieldTy)
					 
					 
      val trans1 =  if Option.isSome(refVec1Exp)
 		    then vec1Exp
 		    else makePrim'(BV.op_probe, [invTransformField, vec1Exp], [invTransformFieldTy, vecTy], vecTy)
				  
      val derv = makePrim'(BV.op_probe, [dInvTransformField, vec1Exp], [dTransformFieldTy, vecTy], matTy)
      val trans2 = makePrim'(BV.op_inner_tt, [derv, vec2Exp], [matTy, vecTy], vecTy)

      val ((_, exitFuncVar), _) = hiddenExitFuncResult
      val seqResult = AST.E_Apply((exitFuncVar, span), [trans1, trans2], vec2SeqTy)
				 (* val ret = AST.E_Slice(seqResult, [SOME(AST.E_Lit(Literal.intLit 0))], Ty.realTy) *)
				 (*other exit call...*)
     in
      seqResult
     end
   | buildWorldIntersectsExp (false, meshExp, cellExp, cellIntExp, refVec1Exp, vec1Exp, vec2Exp) = raise Fail "later"



													 
 local
  (*build _exit, build exit func, enter func*)
  val FT.Mesh(mesh) = meshData
  val cellHiddenExit = Atom.atom "$cellExit" (*hidden*)
  val cellHiddenExitTy = Ty.T_Fun([meshTy, cellTy, Ty.T_Int, vecTy, vecTy, vecTy], vec2SeqTy)
  fun cellHiddenExitFun(v1, v2, v3, v4, v5, v6) = buildWorldIntersectsExp(FemData.isAffine mesh, v1, v2, v3, v4, v5, v6)


									 

  val posExit = Atom.atom (FemName.posExit)
  val posExitTy = Ty.T_Fun([posTy, vecTy], Ty.realTy)
  fun posExitFun([v1, v2]) =
      let
       val mesh = AST.E_ExtractFem(v1, meshData)
       val cellInt = AST.E_ExtractFemItem(v1, Ty.T_Int, (FO.CellIndex, posData))
       val cellExp =  AST.E_LoadFem(cellData, SOME(mesh), SOME(cellInt))
       val refPos = AST.E_ExtractFemItem(v1, vecTy, (FemOpt.RefPos, meshData))
       val transformField = makeTransformFunc([cellExp])
       val worldPos =  makePrim'(BV.op_probe, [transformField, refPos], [transformFieldTy, vecTy], vecTy) (*TODO: fix worldPos*)
       val result = cellHiddenExitFun(mesh, cellExp, cellInt, SOME(refPos), worldPos, v2)

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
       val result = cellHiddenExitFun(mesh, cellExp, cellInt, NONE, v2, v3)
      in
       result
      end
	

	
 in
 val cellHiddenResult = (cellHiddenExit,cellHiddenExitFun, cellHiddenExitTy)
 val posExitResult = (posExit, posExitFun, posExitTy)
 val cellEnterResult = (cellEnter, cellEnterFun, cellEnterTy)
 end

 fun buildWorld(true, posExp, vecExp) =
     let
      (*get the int, ref*)
      val (_, refExitPosFunc, _) = refExitPosReplace
      val mesh = AST.E_ExtractFem(posExp, meshData)
      val cellInt = AST.E_ExtractFemItem(posExp, Ty.T_Int, (FO.CellIndex, posData))
      val cellExp =  AST.E_LoadFem(cellData, SOME(mesh), SOME(cellInt))
				  
      val invTransformField = makeInvTransformFunc([cellExp])
      val transformField = makeTransformFunc([cellExp])
      val invertVar = if dim = 2
 		      then BV.fn_inv2_f
 		      else if dim = 3
 		      then BV.fn_inv3_f
 		      else raise Fail "dim <> 2 and dim <> 3"
      val dTransformVal = makePrim'(BV.op_Dotimes, [transformField], [invTransformFieldTy], dTransformFieldTy)
      val invertedField = makePrim'(invertVar, [dTransformVal], [dTransformFieldTy], dTransformFieldTy)
      val refPos = AST.E_ExtractFemItem(posExp, vecTy, (FemOpt.RefPos, meshData))
      val probe = makePrim'(BV.op_probe, [invertedField, refPos], [dTransformFieldTy, vecTy], matTy)
      val newDpos = makePrim'(BV.op_inner_tt, [probe, vecExp], [matTy, vecTy], vecTy)

				   
      val mesh = AST.E_ExtractFem(cellExp, meshData)
      val result = refExitPosFunc([mesh, posExp, newDpos])
     in
      result
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

							  
 val refActualFuncs = [hiddenExitFuncResult, cellFunc,cellFunc4, refPosExitFunc]
 val (refActualFuncsInfo, refActualFuncDcls) = ListPair.unzip refActualFuncs
 val refReplaceFuncs = [exitFuncResult, refExitPosReplace, enterFuncResult]

 val posReplaceFuncs = [posExitResult, posExitPosResult]
 val posActualFuncs = []
 val (posActualFuncInfo, posActualFuncDcls) = ListPair.unzip posActualFuncs
 val results = {ref = (refActualFuncsInfo, refActualFuncDcls, refReplaceFuncs),
		pos = (posActualFuncInfo,posActualFuncDcls,posReplaceFuncs),
		cell = ([],[],[cellEnterResult])}


in
 results
end

end
