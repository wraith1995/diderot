(* check-globals.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)

structure GenFemAst  =
struct
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

val err = TypeError.error

datatype token = datatype TypeError.token

fun validateFemType (cxt, env, tyName) femTyDef fileinfo =
    let
     (*TODO: parse the file and make these correct...*)
     (*TODO : check against the file in the following.*)
     val parsedJson = CF.loadJson(fileinfo, cxt)

	

    in
     (case (femTyDef, parsedJson)
       of (PT.T_Mesh, SOME(parsed)) =>
	  let
	   val (tempMesh, consts, geometry, meshOpts) = CF.parseMesh(env, cxt, tyName, parsed)
	  in
	   (Option.mapPartial (fn x => SOME(x, NONE)) tempMesh, consts, geometry, meshOpts)
	  end
	| (PT.T_Space(mesh, shape), SOME(parsed)) =>
	  let
	   (*validate shape*)
	   val shapeTy = CheckExpr.checkShape(env, cxt, shape)
	   val shapeOpt = (case shapeTy
			    of Ty.Shape(dims) =>
			       (List.map (fn x =>
					     (case x
					       of Ty.DimConst(d) => SOME(d)
						| _ => NONE)) dims)
				 
			     | _ => [NONE]
			  (*end case*))
	   val shape' = List.foldr (fn (x, y) => (case (x, y)
						   of (SOME(d), SOME(y')) => SOME(d :: y')
						    | _ => NONE)) (SOME([])) shapeOpt
 	   (* validate shape against file? *)
	   val meshType = Option.mapPartial (FT.extractMesh) (Option.map (fn (x,y) => x) (Option.mapPartial TU.extractFemType
													    (Option.mapPartial
													       (fn x => SOME(TypeEnv.findDef x))
													       (E.findTypeEnv(env, mesh)))))


	   val spaceType = Option.map (fn (shape'', mesh) => CF.parseSpace(env, cxt, tyName, [], shape'', mesh, parsed)) (case (shape',meshType)
															   of (SOME(a),SOME(b)) => SOME(a,b)
															    | _ => NONE (* end case*))
				      
	  in
	   (case (spaceType, meshType)
	     of (SOME(SOME(space''), consts), SOME(meshType')) => (SOME(space'', SOME(mesh)), consts, [], [])
	      (*NOTE: these errors could be improved.*)
	      | (_, SOME(_)) => ((err (cxt, [
				       S "Declared a function space type ",
				       A(tyName), S" with an invalid space "])); (NONE, [],[], []))
	      | (_, NONE) => ((err (cxt, [
				    S "Declared a function space type ",
				    A(tyName), S" with an underlying type that is not a mesh or not defined", A(mesh)])); (NONE, [], [], []))
	   (*end case*))
	  end

	| (PT.T_Func(space), SOME(parsed)) =>
	  let
	   val spaceType = Option.mapPartial (FT.extractSpace) (Option.map (fn (x,y) => x) (Option.mapPartial TU.extractFemType (Option.mapPartial
																   (fn x => SOME(TypeEnv.findDef x))
																   (E.findTypeEnv(env, space)))))

	   val funcType = Option.map (fn s => CF.parseFunc(env, cxt, tyName, NONE, s, parsed)) spaceType
	  in
	   (case (funcType,spaceType)
	     of (SOME(SOME(func), consts),SOME(space')) => (SOME(func, SOME(space)), consts, [], [])
	      | (SOME(NONE,_), SOME(_)) => ((err (cxt, [
						  S "Declared a femFunction type ",
						  A(tyName), S" with an invalid definition in the file or incompatability with the space.", A(space)])); (NONE, [],[], []))
	      |  (_, NONE) =>  ((err (cxt, [
				      S "Declared a femFunction type ",
				      A(tyName), S" with an underlying type that is not a space or not defined", A(space)])); (NONE, [],[], []))
	      (* end case *)
	      | _ => raise Fail "cast not analyzed")
	  end
	| (_, NONE) => (err (cxt, [S "Unable to parse file for declared fem type:", A(tyName) ]); (NONE, [],[], []))
	| _ => raise Fail ("Non fem type passed to validateFemType: " ^ Atom.toString(tyName))
     (* end case *))
    end

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
fun makeRefCellInside(env, cxt, span, dim, refCellInfo, posVar, meshVar, insert, meshData, meshTy) =
    let
     val insideVec = Ty.vecTy dim
     val FemData.RefCellData({ty=refCellClass, eps,...}) = refCellInfo
     val epsExpr = AST.E_Lit(Literal.Real(eps))
     val one = AST.E_Lit(Literal.Real(RealLit.one))
     val onePlusEps = makePrim'(BV.add_tt, [one, epsExpr], [Ty.realTy, Ty.realTy], Ty.realTy)
     val posVar' = posVar
     val meshVar' = meshVar
     fun greaterThanTest top bot =
	 let
	  val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.gt_rr)
	 in
	  makePrim(AST.E_Prim(BV.gt_rr, tyArgs, [top, bot], Ty.T_Bool), [Ty.realTy, Ty.realTy])
	 end
     fun subs vecExpr i = AST.E_Slice(vecExpr, [SOME(AST.E_Lit(Literal.intLit i))], Ty.realTy)


				     
				     
    in
     (case refCellClass
       of FemData.KSimplex(dim) =>
	  let
	   val epsilonVec = AST.E_Tensor(List.tabulate(dim, fn x => epsExpr),
					 insideVec)
	   val oneVec = AST.E_Tensor(List.tabulate(dim, fn x => one), insideVec)
	   val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.add_tt)
	   val adjustedPos = makePrim(AST.E_Prim(BV.add_tt, tyArgs, [epsilonVec, posVar'], insideVec), [insideVec, insideVec])
	   val adjustedPosAcc = List.tabulate(dim,
					      fn x => subs adjustedPos x)
	   val posTests = List.map (fn x => greaterThanTest x (AST.E_Lit(Literal.Real(RealLit.zero true)))) adjustedPosAcc
	   val (tyArgs', Ty.T_Fun(domTy', rngTy')) = TU.instantiate(Var.typeOf BV.op_inner_tt)
	   val sumVar = makePrim(AST.E_Prim(BV.op_inner_tt, tyArgs', [oneVec, posVar'], Ty.realTy), [insideVec, insideVec])
	   val sumTest = greaterThanTest onePlusEps sumVar
	   val endResult = makeAnds (sumTest::posTests)
				    
	  in
	   endResult
	  end
	| FemData.KCube(_) =>
	  let
	   val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_norm_t)
	   val norm = makePrim(AST.E_Prim(BV.op_norm_t, tyArgs, [posVar'], Ty.realTy), [insideVec])
	   val endTest = greaterThanTest onePlusEps norm
	  in
	   endTest
	  end
	| FemData.Other(_, _) =>
	  (case insert
	    of SOME(file) =>
	       let
		val _ = ()
	       in
		AST.E_ExtractFemItemN([meshVar, posVar, epsExpr],
				      [meshTy, insideVec, Ty.realTy],
				      Ty.T_Bool,
				      (FemOpt.InsideInsert(Atom.atom file), meshData), NONE)
	       end
	     | NONE => raise Fail "insert error"))
    end
(*TODO: make parameterized newton inverse*)
(*TODO: make meshPosSearch*)
fun itterStartPos(cxt, span, refCellClass, vecTy, dim, start) =
    (case refCellClass
      of FemData.KCube(dim) => AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(RealLit.zero false))), vecTy)
       | FemData.KSimplex(dim) =>
	 let
	  val digits = List.tabulate(100, fn x => 3)
	  val exp = IntInf.fromInt (0)
	  val third = RealLit.fromDigits({isNeg = false, digits = digits, exp = exp})
	 in
	  AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(third))), vecTy)
	 end
       | FemData.Other(dim, _) => (case start
				       (*todo: check length*)
				 of SOME(s) => AST.E_Tensor((List.map (fn x => AST.E_Lit(Literal.Real(x))) s),
							    Ty.T_Tensor(Ty.Shape([Ty.DimConst dim])))
				  | NONE => raise Fail "insert an error here")
    (* end case*))

(*build n derivative functions of the transform*)
(*also refcell insideCheck...*)
(*with these functions, we can do more interesting things...*)
      
fun dnT(env, cxt, span, mesh, m) =
    let
     val FT.Mesh(meshData) = mesh
     val dim = FT.meshDim meshData
     val range = dim::(List.tabulate(m, fn x => dim))
     val dimC = Ty.DimConst(dim)
     val dimT = Ty.T_Tensor(Ty.Shape([dimC]))
     val meshTy = Ty.T_Fem(mesh, NONE)
     val meshPosData = FT.MeshPos(meshData)
     val meshPosTy = Ty.T_Fem(meshPosData, SOME(FT.nameOf mesh))
     fun rangeC n = Ty.Shape(List.take(List.map Ty.DimConst range, n + 1))
     fun resultFieldTy n = Ty.T_Field{diff=Ty.DiffConst(NONE),
				      dim = dimC,
				      shape = rangeC n}
     fun resultTensorTy n = Ty.T_Tensor(rangeC n)
     fun functionTy n = Ty.T_Fun([dimT,
				  Ty.T_Int,
				  meshTy
				 ], resultTensorTy n)
     val params = [Var.new (Atom.atom "pos", span, AST.FunParam, dimT),
		   Var.new (Atom.atom "cell", span, AST.FunParam, Ty.T_Int),
		   Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)]
     val [vecIn, cellIn, meshIn] = params
				     
     fun functionName n = Atom.atom ((Atom.toString (FT.nameOf mesh)) ^ "_transform_" ^ (Int.toString n))
     fun functionVar n = Var.new(functionName n , span, AST.FunVar, functionTy n)

     val baseTransform = AST.E_FemField(AST.E_Var(meshIn, span), AST.E_Var(meshIn, span), SOME(AST.E_Var(cellIn, span)), resultFieldTy 0, FemOpt.Transform, NONE)

				       
				       
				       
     fun makeTransform 0 = baseTransform
       | makeTransform n = makePrim'(BV.op_Dotimes,
				     [ makeTransform (n-1)],
				     [resultFieldTy (n-1)],
				     (resultFieldTy (n)))

     fun probeTransform n = makePrim'(BV.op_probe, [makeTransform n, AST.E_Var(vecIn, span)],
				      [resultFieldTy n, dimT],
				      (resultFieldTy n))


     fun transformFunction n =
	 let
	  val at = functionName n
	  val var = functionVar n
	 in
	  (at,
	   var,
	   AST.D_Func(var, params, AST.S_Block([AST.S_Return((probeTransform n))])))
	 end
	   

	   
     val results = List.tabulate(m+1, transformFunction)
     val (_, firstFun, _) = List.nth(results, 0)
     val buildWorldPosFunTy = Ty.T_Fun([meshPosTy],meshPosTy)
     val meshPosParam = Var.new(Atom.atom"pos", span, AST.FunParam, meshPosTy)



     val dumbWordPos = Ty.T_Fun([meshPosTy], dimT)
     val dumbWorldPosName = FemData.functionNameMake mesh (FemName.worldPos)
     val dumWorldPosOutSideName = Atom.atom (FemName.worldPos)
     val dumbWorldPosFunVar = Var.new(dumWorldPosOutSideName, span, AST.FunVar, dumbWordPos)
     val dumbBody = AST.S_Return(AST.E_Apply(
				   (firstFun, span),
				   [AST.E_ExtractFemItem(AST.E_Var(meshPosParam, span),
							 dimT,
							 (FemOpt.RefPos, meshPosData)),
				    AST.E_ExtractFemItem(AST.E_Var(meshPosParam, span),
							 Ty.T_Int,
							 (FemOpt.CellIndex, meshPosData)),
				    AST.E_ExtractFem(AST.E_Var(meshPosParam, span),mesh)], dimT))
     val dumbFun = AST.D_Func(dumbWorldPosFunVar,[meshPosParam],dumbBody)
     val dumbSpecialResult = (dumWorldPosOutSideName, dumbWorldPosFunVar, dumbFun)
			       
    in
     (dumbSpecialResult, results)
    end
fun makeInvVar(dim) =
    if dim = 2
    then BV.fn_inv2_f
    else if dim = 3
    then BV.fn_inv3_f
    else raise Fail "invalid dim;"
(*TODO: this resgtriction should be left as I should be able to do arbitrary inverses of a matrix*)
(*TODO: a \ operator would be nice.*)


fun makeMeshPosSearch(env, cxt, span, refCellClass, meshData, newtonTol, newtonAttempts, contraction, killAfterTol, posExpr, meshExpr, insideFunction, start) =
    let

     (*set up types*)
     val (FT.Mesh(mesh)) = meshData
     val meshTy = Ty.T_Fem(meshData, NONE)
     val dim = FemData.meshDim mesh
     val insideVec = Ty.vecTy dim
     val insideMat = Ty.matTy dim
     val dimConst = Ty.DimConst(dim)
     val inf = Ty.DiffConst(NONE)
     val transformFieldTy = Ty.T_Field({diff=inf, dim = dimConst, shape=Ty.Shape([dimConst])})
     val dTransformFieldTy = Ty.T_Field({diff=inf, dim = dimConst, shape=Ty.Shape([dimConst, dimConst])})
     val meshPosData = FT.MeshPos(mesh)
     val meshPosType = Ty.T_Fem(meshPosData, SOME(FT.nameOf meshData))
     val zero = AST.E_Lit(Literal.intLit 0)
     val one = AST.E_Lit(Literal.intLit 1)
			
     val initStms = []

     (*setup vars, constant expressions, inserts*)
     val startPosition = itterStartPos(cxt, span, refCellClass, insideVec, dim, start)
     val invVar = makeInvVar(dim)
     val tolExpr = newtonTol
     val maxNewtonExpr = newtonAttempts

     val numCellExpr' = AST.E_ExtractFemItem(meshExpr, Ty.T_Int, (FemOpt.NumCell, meshData))
     val numCellExpr = makePrim'(BV.sub_ii, [numCellExpr', one], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
     val dummExp = AST.E_ExtractFemItem(meshExpr, Ty.T_Int, (FemOpt.StartCell, meshData))
     val numCellVar = Var.new(Atom.atom "numCell", span, Var.LocalVar, Ty.T_Int)
     val numCellAssign = AST.S_Decl((numCellVar, SOME(numCellExpr)))
     val numCellVarExpr = AST.E_Var((numCellVar, span))
     val initStms = numCellAssign::initStms
				     
     val maxNewtonExpr' = makePrim'(BV.add_ii, [maxNewtonExpr, one], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
     val maxTotal = makePrim'(BV.mul_ii, [maxNewtonExpr', numCellVarExpr], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
     val itters = makePrim'(BV.range, [zero, maxTotal], [Ty.T_Int, Ty.T_Int], Ty.T_Sequence(Ty.T_Int, NONE))
     val itterVar = Var.new(Atom.atom "itter", span, AST.IterVar, Ty.T_Int)
     val itter = (itterVar, itters)
		   

     (*setup the vars that we use in the itteration.*)
     val cellItterVar = Var.new(Atom.atom "cellInt", span, Var.LocalVar, Ty.T_Int)
     val newtonItterVar = Var.new(Atom.atom "newtonInt", span, Var.LocalVar, Ty.T_Int)
     val newtonReset =  AST.S_Decl((newtonItterVar, SOME(zero)))
     val initStms = newtonReset::AST.S_Decl((cellItterVar,  SOME(dummExp)))::initStms;
     val cellItterVarExp = AST.E_Var((cellItterVar, span))
     val newtonItterVarExp = AST.E_Var((newtonItterVar, span))

     val newPosVar = Var.new(Atom.atom "xn", span, Var.LocalVar, insideVec)
     val newPosInit = AST.S_Decl((newPosVar, SOME(startPosition)))
     val newPosReset = AST.S_Assign((newPosVar, span), startPosition)
     val newPosVarExpr = AST.E_Var((newPosVar, span))
     val initStms = newPosInit::initStms
     val updateVar = Var.new(Atom.atom "delta", span, Var.LocalVar, insideVec)


     (*setup the fields*)	
     val transformField = AST.E_FemField(meshExpr, meshExpr, SOME(cellItterVarExp), transformFieldTy, FemOpt.Transform, NONE)
     val transformFieldModPos = makePrim'(BV.sub_ft, [transformField, posExpr], [transformFieldTy, insideVec], transformFieldTy)
     val dTransformField = makePrim'(BV.op_Dotimes,[transformField], [transformFieldTy], dTransformFieldTy)
     val invDTransformField = makePrim'(invVar, [dTransformField], [dTransformFieldTy], dTransformFieldTy)	    
     (*setup the key tests:*)
     val tolExpr' = makePrim'(BV.mul_rr, [tolExpr, tolExpr], [Ty.realTy, Ty.realTy], Ty.realTy)
     val deltaNorm = makePrim'(BV.op_inner_tt, [AST.E_Var((updateVar, span)), AST.E_Var((updateVar, span))], [insideVec, insideVec], Ty.realTy)
     val deltaNormTest = makePrim'(BV.gte_rr, [tolExpr', deltaNorm], [Ty.realTy, Ty.realTy], Ty.T_Bool)
				  
     val insideTest = insideFunction([meshExpr, newPosVarExpr])
     val makeMeshPos = AST.E_ExtractFemItemN([meshExpr, cellItterVarExp,newPosVarExpr, posExpr], [meshTy, Ty.T_Int, insideVec, insideVec], meshPosType, (FemOpt.AllBuild, meshPosData), NONE)
     val invalidMeshPos = AST.E_ExtractFemItemN([meshExpr], [ meshTy], meshPosType, (FemOpt.InvalidBuild, meshPosData), NONE)

     val succesReturn = AST.S_Return(makeMeshPos)
     val failReturn = AST.S_Return(invalidMeshPos)

     val incC = AST.S_Assign((cellItterVar, span),
			     makePrim'(BV.add_ii,
				       [cellItterVarExp, one],
				       [Ty.T_Int, Ty.T_Int], Ty.T_Int))
     val incN = AST.S_Assign((newtonItterVar, span),
			     makePrim'(BV.add_ii,
				       [newtonItterVarExp, one],
				       [Ty.T_Int, Ty.T_Int], Ty.T_Int))

     val nTest = makePrim'(BV.gte_ii,
			   [newtonItterVarExp, maxNewtonExpr],
			   [Ty.T_Int, Ty.T_Int],
			   Ty.T_Bool)
     val cTest = makePrim'(BV.gte_ii,
			   [cellItterVarExp, numCellVarExpr],
			   [Ty.T_Int, Ty.T_Int],
			   Ty.T_Bool)
     val emptyBlock = AST.S_Block([])


     val ifStm = AST.S_IfThenElse(
	  deltaNormTest,
	  AST.S_Block([
		      AST.S_IfThenElse(insideTest,
				       AST.S_Block([ succesReturn]),
				       emptyBlock)
		     ]),
	  emptyBlock)
     val incrementNewton = incN

     val ifStm' = AST.S_IfThenElse(
	  deltaNormTest,
	  AST.S_Block([
		      AST.S_IfThenElse(insideTest,
				       AST.S_Block([succesReturn]),
				       emptyBlock)
		     ]), emptyBlock)
     val cellNTest = AST.S_IfThenElse(
	  nTest,
	  AST.S_Block([newtonReset,
		       AST.S_IfThenElse(cTest,
					AST.S_Block([failReturn]),
					AST.S_Block([incC]))]), emptyBlock)

     (*incrementN, icrementC, *)
     (*reoriegnet: extraSetup, updateDetlaStm, updateCurrentPosStm*)
     (*if data, buildFirstForloop: needs mesh and original pos*)
     (*if conservative, build second capable*)

     val (setupVars, updateDeltaStm, updateCurrentPosStm) =
	 if contraction
	 then
	  let
	   val probeD = makePrim'(BV.op_probe, [invDTransformField, startPosition], [dTransformFieldTy, insideVec], insideMat)
	   val matVar = Var.new(Atom.atom "A", span, AST.LocalVar, insideMat)
	   val probeDAssign = AST.S_Decl((matVar, SOME(probeD)))
	   val dotFields = makePrim'(BV.op_inner_tf,  [AST.E_Var(matVar, span), transformFieldModPos], [insideMat, transformFieldTy], transformFieldTy)
	   val probeUpdate = makePrim'(BV.op_probe, [dotFields, newPosVarExpr], [transformFieldTy, insideVec], insideVec)
	   val updateDeltaStm = AST.S_Decl((updateVar, SOME(probeUpdate)))
	   val updateCurrentPosExpr = makePrim'(BV.sub_tt, [AST.E_Var(newPosVar,span), AST.E_Var(updateVar, span)], [insideVec, insideVec], insideVec)
	   val updateCurrentPosStm = AST.S_Assign((newPosVar,span), updateCurrentPosExpr)
	  in
	   ([probeDAssign],updateDeltaStm,updateCurrentPosStm)
	  end
	 else
	  let
	   val dotFields' = makePrim'(BV.op_inner_ff,  [invDTransformField, transformFieldModPos], [dTransformFieldTy, transformFieldTy], transformFieldTy) 
	   val probeUpdate = makePrim'(BV.op_probe, [dotFields', newPosVarExpr], [transformFieldTy, insideVec], insideVec) 
	   val updateDeltaStm = AST.S_Assign((updateVar, span), probeUpdate)
	   val updateCurrentPosExpr = makePrim'(BV.sub_tt, [AST.E_Var(newPosVar,span), AST.E_Var(updateVar, span)], [insideVec, insideVec], insideVec)
	   val updateCurrentPosStm = AST.S_Assign((newPosVar,span), updateCurrentPosExpr)

	  in
	   ([], updateDeltaStm, updateCurrentPosStm)
	  end
     val currentlyAt = makePrinStatement("New RefPos Estimate is: ", [newPosVarExpr], "\n")
     val currentDelta = makePrinStatement("New Delta Estimate is: ", [AST.E_Var((updateVar, span)), deltaNorm], "\n")

     val secondForLoop = AST.S_Foreach(itter, AST.S_Block(setupVars@[
							  updateDeltaStm,
							  updateCurrentPosStm,
							  ifStm,
							  incrementNewton,
							  cellNTest]))
     val loopStms =
	 (case FemData.meshAccInsert mesh
	   of NONE => [secondForLoop]
	    | SOME(insert, cons) =>
	      let
	       val opt = (FemOpt.NearbyCellQuery(insert), FemData.Mesh(mesh))
	       val intSeq = Ty.T_Sequence(Ty.T_Int, NONE)
	       val cellExpr = AST.E_ExtractFemItem2(meshExpr,posExpr,insideVec,intSeq, opt);
	       val printSeq = makePrinStatement("Cells:", [cellExpr],"\n")
	       val currentlyAt = makePrinStatement("RefPos Estimate is: ", [newPosVarExpr], "\n")
	       val newCellsVar = Var.new(Atom.atom "yayCells", span, AST.LocalVar, Ty.T_Sequence(Ty.T_Int, NONE))
	       val newAssign = AST.S_Decl((newCellsVar, SOME(cellExpr)))

	       val itterVar = Var.new(Atom.atom "cellItter", span, AST.IterVar, Ty.T_Int)
	       val newtonItterVar = Var.new(Atom.atom "newtonItter", span, AST.IterVar, Ty.T_Int)
	       val itter1 = (itterVar, AST.E_Var(newCellsVar, span))
	       val itter2 = (newtonItterVar, makePrim'(BV.range, [zero, maxNewtonExpr], [Ty.T_Int, Ty.T_Int], Ty.T_Sequence(Ty.T_Int, NONE)))

	       val tempAssignment = AST.S_Assign((cellItterVar, span), AST.E_Var(itterVar,span))
	       val tempAssignment' = AST.S_Assign((cellItterVar, span), dummExp)
						 
	       (*TODO: rewrite as above maybe - interesting to compare these two*)
	       val loop = AST.S_Foreach(
		    itter1,
		    AST.S_Block((newPosReset::tempAssignment::setupVars)@[
				AST.S_Foreach(itter2,
					      AST.S_Block([
							  updateDeltaStm,
							  updateCurrentPosStm,
							  ifStm'
			       ]))]
		   ))
	      in
	       if cons
	       then [newAssign, loop]
	       else [newAssign, loop, tempAssignment', secondForLoop]
	      end
	 (*end case*))

     val bodyStms = loopStms@[failReturn]
     val bodyStms' = initStms@bodyStms
     val body = AST.S_Block(bodyStms')

			   
			   
    in
     body
    end

fun makeNewtonInversesBody(env, cxt, span, refCellClass, meshData, newtonTol, newtonAttempts, contraction, killAfterTol, posExpr, cellIntExpr, meshExpr, insideFunction, start, meshTy) =
    let
     (*Setup useful types*)
     val dim = FemData.meshDim meshData
     val insideVec = Ty.vecTy dim
     val insideMat = Ty.matTy dim
     val dimConst = Ty.DimConst(dim)
     val inf = Ty.DiffConst(NONE)
     val transformFieldTy = Ty.T_Field({diff=inf, dim = dimConst, shape=Ty.Shape([dimConst])})
     val dTransformFieldTy = Ty.T_Field({diff=inf, dim = dimConst, shape=Ty.Shape([dimConst, dimConst])}) 

     (*setup useful vars and expressions*)
     val tolExpr = newtonTol
     val startPosition = itterStartPos(cxt, span, refCellClass, insideVec, dim, start)
     val invVar = makeInvVar(dim)
     val newPosVar = Var.new(Atom.atom "xn", span, AST.LocalVar, insideVec)
     val newPosVarExpr = AST.E_Var((newPosVar, span))
     val initPos = AST.S_Decl((newPosVar), SOME(startPosition))
     val updateVar = Var.new(Atom.atom "delta", span, AST.LocalVar, insideVec)
     val itterVar = Var.new(Atom.atom "itter", span, AST.IterVar, Ty.T_Int)
     val itterRange = makePrim'(BV.range,
				[AST.E_Lit(Literal.intLit 0), newtonAttempts],
				[Ty.T_Int, Ty.T_Int], Ty.T_Sequence(Ty.T_Int, NONE)) (*TODO: shouldn't we put a SOME here.*)
			       
     val itter = (itterVar, itterRange)
		   

     (*TODO: should we use newton or contraction - if contraction, lift derivative out...*)
     val transformField = AST.E_FemField(meshExpr, meshExpr, SOME(cellIntExpr), transformFieldTy, FemOpt.Transform, NONE)
     val transformFieldModPos = makePrim'(BV.sub_ft, [transformField, posExpr], [transformFieldTy, insideVec], transformFieldTy)
     val dTransformField = makePrim'(BV.op_Dotimes,[transformField], [transformFieldTy], dTransformFieldTy)
     val invDTransformField = makePrim'(invVar, [dTransformField], [dTransformFieldTy], dTransformFieldTy)
				       
     val deltaNorm = makePrim'(BV.op_inner_tt, [AST.E_Var((updateVar, span)),AST.E_Var((updateVar, span))], [insideVec,insideVec], Ty.realTy)
     val tolExpr' = makePrim'(BV.mul_rr, [tolExpr, tolExpr], [Ty.realTy, Ty.realTy], Ty.realTy)
     val deltaNormTest = makePrim'(BV.gte_rr, [tolExpr', deltaNorm], [Ty.realTy, Ty.realTy], Ty.T_Bool)
     val insideTest = insideFunction([meshExpr, newPosVarExpr])
     (* val combinedTest = AST.Andalso(deltaNormTest, insideTest) (*note: we might want to combined them...... it might be better to do delta test then inside test????*) *)
     val meshPosData = FT.MeshPos(meshData)
     val meshPosTy = Ty.T_Fem(meshPosData, SOME(FemData.nameOf (FT.Mesh(meshData))))
     val makeMeshPos = AST.E_ExtractFemItemN([meshExpr, cellIntExpr, newPosVarExpr, posExpr], [meshTy, Ty.T_Int, insideVec, insideVec], meshPosTy, (FemOpt.AllBuild, meshPosData), NONE)
     val invalidMeshPos = AST.E_ExtractFemItemN([meshExpr], [ meshTy], meshPosTy, (FemOpt.InvalidBuild, meshPosData), NONE)
     val succesIntermediate = Var.new (Atom.atom "dump", span, Var.LocalVar, meshPosTy)
     val sia = AST.S_Assign((succesIntermediate,span), (makeMeshPos))
     val failReturn = AST.S_Return(invalidMeshPos)
     val succesReturn = AST.S_Return(makeMeshPos)
     val rightBefore = makePrinStatement("Hello\n",[newPosVarExpr],"\n");
     val ifStm = if true
		 then AST.S_IfThenElse(deltaNormTest, (*we might want to change this test to be this and that, rather than a nested if.*)
				       AST.S_Block([AST.S_IfThenElse(insideTest,
								     AST.S_Block([succesReturn]), (*converged and arrived*)
								     AST.S_Block([failReturn]))]), (*converged and failed; somehow*)
				       AST.S_Block([]))
		 else AST.S_IfThenElse(makeAnd deltaNormTest insideTest, AST.S_Block([succesReturn]), AST.S_Block([failReturn]))

				      
     val completeBody =
	 if contraction
	 then
	  let

	   val probeD = makePrim'(BV.op_probe, [invDTransformField, startPosition], [dTransformFieldTy, insideVec], insideMat)
	   val matVar = Var.new(Atom.atom "A", span, AST.LocalVar, insideMat)
	   val probeDAssign = AST.S_Decl((matVar, SOME(probeD)))
	   val dotFields = makePrim'(BV.op_inner_tf,  [AST.E_Var(matVar, span), transformFieldModPos], [insideMat, transformFieldTy], transformFieldTy)
	   val probeUpdate = makePrim'(BV.op_probe, [dotFields, newPosVarExpr], [transformFieldTy, insideVec], insideVec)
	   val updateDeltaStm = AST.S_Decl(updateVar, SOME(probeUpdate))
	   val updateCurrentPosExpr = makePrim'(BV.sub_tt, [AST.E_Var(newPosVar,span), AST.E_Var(updateVar, span)], [insideVec, insideVec], insideVec)
	   val updateCurrentPosStm = AST.S_Assign((newPosVar,span), updateCurrentPosExpr)

	   val forLoop = AST.S_Foreach(itter, AST.S_Block([
							  updateDeltaStm,
							  updateCurrentPosStm,
							  ifStm]))
	  in
	   AST.S_Block( [
		       initPos,
		       probeDAssign,
		       forLoop,
		       failReturn])
	  end
	 else
	  let

	   val dotFields' = makePrim'(BV.op_inner_ff,  [invDTransformField, transformFieldModPos], [dTransformFieldTy, transformFieldTy], transformFieldTy) 
				     
	   val probeUpdate = makePrim'(BV.op_probe, [dotFields', newPosVarExpr], [transformFieldTy, insideVec], insideVec) (**)
	   val updateDeltaStm = AST.S_Assign((updateVar, span), probeUpdate)
	   val updateCurrentPosExpr = makePrim'(BV.sub_tt, [AST.E_Var(newPosVar,span), AST.E_Var(updateVar, span)], [insideVec, insideVec], insideVec)
	   val updateCurrentPosStm = AST.S_Assign((newPosVar,span), updateCurrentPosExpr)

	   val forLoop = AST.S_Foreach(itter, AST.S_Block( [
							  updateDeltaStm,
							  updateCurrentPosStm,
							  ifStm]))
	  in
	   AST.S_Block([
		       initPos,
		       forLoop,
		       failReturn])
	  end
	    
    in
     completeBody
    end
fun makeFunctionRegistration(name, argTys, resultTy, f) = let
 val atomName = Atom.atom name
 val argCount = List.length argTys
 val functionTy = Ty.T_Fun(argTys, resultTy)
 val funcCall = (fn x => if List.length x = argCount
			 then f x
			 else raise Fail ("impossible number of args to "^ name ^"; this should not have gotten pass the type checker.")
		)
in
 (atomName, funcCall, functionTy)
end

fun makeFemMethods (cxt, env, tyName, span) femTyDef file cellAccData femInfo =
    (case (femInfo, femTyDef)
      of ((func as FT.Func(f), n), PT.T_Func(space)) =>
	 let
	  val funcName = FT.nameOf func

	  val spaceFunc = Atom.atom (FemName.space)
	  val space' = FT.spaceOf(FT.Func(f))
	  val spaceName = FT.nameOf space'
	  val funcTy' = Ty.T_Fem(func, SOME(spaceName))
	  val funcTy = Ty.T_Named(funcName, funcTy')
	  val mesh = FT.meshOf func
	  val FT.Mesh(meshData) = mesh
	  val meshName = FT.nameOf mesh
				   
	  val spaceType = let
	   val m = FT.meshOf(space')
	   val meshName = FT.nameOf(m)
	  in Ty.T_Fem(space', SOME(meshName)) end
	  val spaceFuncTy = Ty.T_Fun([Ty.T_Named(tyName,Ty.T_Fem(femInfo))],Ty.T_Named(space, spaceType))
	  val spaceFuncVar = Var.new (spaceFunc, span, AST.FunVar, spaceFuncTy)
	  val param = Var.new (Atom.atom "arg0", span, AST.FunParam, Ty.T_Named(tyName,Ty.T_Fem(femInfo)))
			      
	  val spaceFuncBody = AST.S_Block([AST.S_Return(AST.E_ExtractFem(AST.E_Var(param,span),space'))])
	  val spaceFun = AST.D_Func(spaceFuncVar, [param], spaceFuncBody)

	  val funcArgTy = Ty.T_Named(tyName, Ty.T_Fem(femInfo))
	  val funcCellTy = Ty.T_Fem(func, SOME(funcName))
	  val FT.Mesh(meshData) = FT.meshOf func
	  val meshCell = FT.MeshCell(meshData)
	  val funcCellVal = FT.FuncCell(f)
	  val meshCellTy = Ty.T_Fem(meshCell, SOME(meshName))
	  val funcCellTy = Ty.T_Fem(funcCellVal, SOME(funcName))
				   
	  val funcCell = Atom.atom (FemName.funcCell)
	  val funcCellType = Ty.T_Fun([funcArgTy, meshCellTy], funcCellTy)
	  val funcCellFuncVar = Var.new (funcCell, span, AST.FunVar, funcCellType)

	  val param = Var.new (Atom.atom "arg0", span, AST.FunParam, funcArgTy)
	  val param' = Var.new (Atom.atom "arg1", span, AST.FunParam, meshCellTy)
	  (*TODO: optional debug information could go here.*)
	  val funcCellBody = AST.S_Return(AST.E_LoadFem(funcCellVal,
							SOME(AST.E_Var(param, span)),
							SOME(AST.E_ExtractFemItem(
							       AST.E_Var(param', span),
							       Ty.T_Int,
							       (FO.CellIndex, meshCell)))))
	  val funcCellFunc = AST.D_Func(funcCellFuncVar,
					[param, param'],
					funcCellBody);
	  local
	   
	   val SOME(meshTyEnv) = Env.findTypeEnv(env, meshName)
	   val SOME(findPosVar) = TypeEnv.findMethod(meshTyEnv, Atom.atom (FemName.meshPos))

	   val dim = FT.meshDim meshData
	   val rangeShape = FT.rangeShape f
	   val diff = Ty.DiffConst(NONE)
	   val dimConst = Ty.DimConst(dim)
	   val rangeShape = Ty.Shape(List.map Ty.DimConst rangeShape)
	   val fieldTy2 = Ty.T_Field({diff = diff,
				      dim = dimConst,
				      shape = rangeShape})
	   val fieldTy1 = Ty.T_Field({diff = diff,
				      dim = dimConst,
				      shape = Ty.Shape([dimConst])})

				    
	   val fieldFunTy = Ty.T_Fun([funcTy], fieldTy2)
	   val fieldAtom = Atom.atom (FemName.field)

	   val badCell = AST.E_Lit(Literal.intLit (~1))
	   fun fieldMakerFun([funExpr]) =
	       let
		val spaceExpr = AST.E_ExtractFem(funExpr, space')
		val meshExpr = AST.E_ExtractFem(spaceExpr, mesh)
		val field2 = AST.E_FemField(funExpr, spaceExpr, SOME(meshExpr),
					    fieldTy2, FemOpt.RefField,
					    SOME(findPosVar, span))
					   
		val field1 = AST.E_FemField(meshExpr, meshExpr, SOME(meshExpr),
					    fieldTy1, FemOpt.InvTransform,
					    SOME(findPosVar, span))
		val field = makePrim'(BV.comp, [field2, field1],
				      [fieldTy2, fieldTy1], fieldTy2)
	       in
		field
	       end
	     | fieldMakerFun _ = raise Fail "impossible: invalid number of args go through typechecker"


	  in
	  val result = (fieldAtom, fieldMakerFun, fieldFunTy)
	  end

	  local
	   val invalidCellAtom = Atom.atom (FemName.invalidCell)
	   val invalidCellTy = Ty.T_Fun([funcTy], funcCellTy)
	   fun invalidCellFun([v1]) = AST.E_LoadFem(FT.FuncCell(f), SOME(v1), SOME(AST.E_Lit(Literal.intLit (~1))))
	     | invalidCellFun (_) = raise Fail "typechecker error."

	  in
	  val invalidCellReplace = (invalidCellAtom, invalidCellFun, invalidCellTy)
	  end

	 in
	  ([(funcCell, funcCellFuncVar), (spaceFunc, spaceFuncVar)],
	   [funcCellFunc, spaceFun],
	   [result, invalidCellReplace],
	   [])
	 end
       | ((FT.Space(s), n), PT.T_Space(meshName)) =>
	 let
	  val meshFunc = Atom.atom "domain"
	  val mesh = FT.meshOf(FT.Space(s))
	  val meshName = FT.nameOf(mesh)
	  val meshType = Ty.T_Fem(mesh, SOME(meshName))
	  val meshFuncTy = Ty.T_Fun([Ty.T_Named(tyName, Ty.T_Fem(femInfo))], Ty.T_Named(meshName, meshType))
	  val meshFuncVar = Var.new(meshFunc, span, AST.FunVar, meshFuncTy)
	  val param = Var.new (Atom.atom "arg0", span, AST.FunParam, Ty.T_Named(tyName,Ty.T_Fem(femInfo)))
	  val meshFuncBody =  AST.S_Block([AST.S_Return(AST.E_ExtractFem(AST.E_Var(param,span),mesh))])
	  val meshFun = AST.D_Func(meshFuncVar, [param], meshFuncBody)
	 in
	  ([(meshFunc, meshFuncVar)], [meshFun],[],[])
	 end
       | ((meshData as FT.Mesh(m), n), PT.T_Mesh) =>
	 let
	  val meshName = FT.nameOf meshData
	  val meshArgTy = Ty.T_Named(tyName, Ty.T_Fem(femInfo))
	  val meshCellTy = Ty.T_Fem(FT.MeshCell(m), SOME(tyName))
	  val meshPosTy = Ty.T_Fem(FT.MeshPos(m), SOME(tyName))
	  val numCell = Atom.atom "numCell"
	  val meshType = Ty.T_Fem(FT.Mesh(m), SOME(tyName))
	  val numCellFuncTy = Ty.T_Fun([meshArgTy],Ty.T_Int)
	  val numCellFuncVar = Var.new (numCell, span, AST.FunVar, numCellFuncTy)
	  val param = Var.new (Atom.atom "arg0", span, AST.FunParam, meshArgTy)
	  val numCellBody = AST.S_Block([AST.S_Return(AST.E_ExtractFemItem(AST.E_Var(param, span), Ty.T_Int, (FemOpt.NumCell, FT.Mesh(m))))])
	  val numCellFun = AST.D_Func(numCellFuncVar, [param], numCellBody)


	  (*NOTE: get list of cells -> we are going to use a extract*)
	  val cells' = Atom.atom "0cells"
	  val cells = Atom.atom "cells"
	  val cellType = Ty.T_Fem(FT.MeshCell(m), SOME(meshName))
	  val cellSeqType = Ty.T_Sequence(cellType, NONE)
	  val cellTypeFuncTy = Ty.T_Fun([meshArgTy], cellSeqType)
	  val cellFuncVar' = Var.new (cells, span, AST.FunVar, cellTypeFuncTy)
	  val param' = Var.new (Atom.atom "arg0", span, AST.FunParam, meshArgTy)

	  val itterVar = Var.new (Atom.atom "i", span, AST.IterVar, Ty.T_Int)
	  val iter = (itterVar, AST.E_Prim(BasisVars.range, [],
			 		   [AST.E_Lit(Literal.Int(IntLit.fromInt 0)),
			 		    AST.E_ExtractFemItem(AST.E_Var(param', span), Ty.T_Int, (FemOpt.NumCell, FT.Mesh(m)))],
			 		   Ty.T_Sequence(Ty.T_Int, NONE)))
	  val comp  = AST.E_Comprehension(AST.E_LoadFem(FT.MeshCell(m), SOME(AST.E_Var(param', span)), SOME(AST.E_Var(itterVar, span))), iter, cellSeqType)
	  val cellsBody = AST.S_Block([AST.S_Return(comp)])
	  val cellsFun' =  AST.D_Func(cellFuncVar', [param'], cellsBody)

	  val cellFuncVar = Var.new (cells, span, AST.FunVar, cellTypeFuncTy)
	  val cellsFun = AST.D_Func(cellFuncVar, [param'],
				    AST.S_Block([AST.S_Return(AST.E_ExtractFemItem(AST.E_Var(param',span), cellSeqType, ( FemOpt.Cells, FT.Mesh(m))))]))


	  fun makeCells vars = (case vars
				 of [v] => AST.E_ExtractFemItem(v, cellSeqType, (FemOpt.Cells, FT.Mesh(m)))
				  | _ => raise Fail "impossible argument to method got pass type checking")

	  val refCellName = Atom.atom "refcell"
	  val refCell = FT.RefCell(m)
	  val refCellTy = Ty.T_Fem(refCell, SOME(meshName))
	  val refCellEnvName = FT.envNameOf refCell
	  val refCellFuncTy = Ty.T_Fun([meshArgTy], refCellTy)
	  val refCellFuncBody = AST.S_Block([AST.S_Return(AST.E_ExtractFemItem(AST.E_Var(param, span), refCellTy, (FemOpt.RefCell, FT.Mesh(m))))])
	  val refCellFuncVar = Var.new(refCellName, span, AST.FunVar, refCellFuncTy)
	  val refCellFunc = AST.D_Func(refCellFuncVar, [param], refCellFuncBody)



	  val meshPosFuncAtom = Atom.atom FemName.meshPos
	  val meshPos = FT.MeshPos(m)
	  val meshPosTy = Ty.T_Fem(meshPos, SOME(meshName))
	  val refCellData = FemData.refCell m
	  val FemData.RefCellData({ty=refCellClass, eps, newtonControl,insideInsert,...}) = refCellData
	  val {contraction, itters, newtonTol, killAfterTol, start} = newtonControl
	  val newtonTolExp = AST.E_Lit(Literal.Real(newtonTol))
	  val newtonAttempts = AST.E_Lit(Literal.intLit itters)

	  val mapDim = FemData.meshDim m
	  val vecTy = Ty.vecTy mapDim

	  val posParam = Var.new (Atom.atom "pos", span, AST.FunParam, vecTy)
	  val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshType)
	  val posExpr = AST.E_Var(posParam,span)
	  val meshExpr = AST.E_Var(meshParam, span)
	  val meshPosfuncTy = Ty.T_Fun([meshArgTy, vecTy], meshPosTy)
	  val meshPosFuncVar = Var.new (meshPosFuncAtom, span, AST.FunVar, meshPosfuncTy)

	  val meshInsideAtom = Atom.atom FemName.isInsideMesh
	  val meshInsideType = Ty.T_Fun([cellType], Ty.T_Bool)
	  fun meshInsideFunc([v1, v2]) =
	      let
	       val temp1 = AST.E_Apply((meshPosFuncVar, span), [v1, v2], meshPosTy)
	       val ret = makePrim'(BV.neq_ii, [AST.E_ExtractFemItem(temp1, Ty.T_Int,  (FO.CellIndex, meshPos)), AST.E_Lit(Literal.intLit (~1))], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
	      in
	       ret
	      end
	    | meshInsideFunc (_) = raise Fail "impossible type checker error"


	  fun makeRefCellInsideFunc vars = (case vars
					     of [v1,v2] => makeRefCellInside(env, cxt, span, mapDim, refCellData, v2, v1, insideInsert, meshData, Ty.T_Fem(meshData, NONE))
					      | _ => raise Fail "impossble got through type checker"
					   (* end case*))
	  val meshPosBody = makeMeshPosSearch(env, cxt, span,
					      refCellClass,
					      FT.Mesh(m),
					      newtonTolExp,newtonAttempts,
					      contraction, killAfterTol,
					      posExpr, meshExpr,
					      makeRefCellInsideFunc,
					      start)

	  val meshPosFunc = AST.D_Func(meshPosFuncVar, [meshParam, posParam], meshPosBody)


	  (*Equality Function:N.op_equ, N.op_neq*)
	  local
	   val cell1 = Var.new(Atom.atom "cell1", span, Var.FunParam, meshCellTy)
	   val cell2 = Var.new(Atom.atom "cell1", span, Var.FunParam, meshCellTy)
	   fun cellInt x = AST.E_ExtractFemItem(AST.E_Var(x,span), Ty.T_Int, (FemOpt.CellIndex, FT.Mesh(m)))
	   val result1 = makePrim'(BV.equ_ii, [cellInt cell1, cellInt cell2], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
	   val result2 = makePrim'(BV.neq_ii, [cellInt cell1, cellInt cell2], [Ty.T_Int, Ty.T_Int], Ty.T_Int)
	   val ret1 = AST.S_Return(result1)
	   val ret2 = AST.S_Return(result2)

	   val atom1 = BasisNames.op_equ
	   val atom2 = BasisNames.op_neq

	   val funTy = Ty.T_Fun([meshCellTy, meshCellTy], Ty.T_Bool)
	   val funVar1 = Var.new(atom1, span, Var.FunVar, funTy)
	   val funVar2 = Var.new(atom2, span, Var.FunVar, funTy)

	   val func1 = AST.D_Func(funVar1, [cell1, cell2], ret1)
	   val func2 = AST.D_Func(funVar2, [cell1, cell2], ret2)

				 

				 
	  in
	  val (eqFunc, nEqFunc) = (func1, func2)
	  val eqFuncPair = (atom1, funVar1)
	  val nEqFuncPair = (atom2, funVar2)
	  end


	  local
	   (*TODO: Invalid stuff*)
	   val invalidPos = Atom.atom (FemName.invalidPos)
	   val invalidPosTy = Ty.T_Fun([meshArgTy], meshPosTy)
	   fun invalidPosFun([v1]) = AST.E_ExtractFemItemN([v1], [meshArgTy], meshPosTy, (FemOpt.InvalidBuild, FT.MeshPos(m)), NONE)
	     | invalidPosFun(_)  = raise Fail "typechecker error."

	   val invalidCell = Atom.atom (FemName.invalidCell)
	   val invalidCellTy = Ty.T_Fun([meshArgTy], meshCellTy)
	   fun invalidCellFun([v1]) = AST.E_LoadFem(FT.MeshCell(m), SOME(v1), SOME(AST.E_Lit(Literal.intLit (~1))))
	     | invalidCellFun (_) = raise Fail "typechecker error."
				       
	  in
	  val invalidReplace = (invalidPos, invalidPosFun, invalidPosTy)
	  val invalidCellReplace = (invalidCell, invalidCellFun, invalidCellTy)
	  end

	  


	 in
	  ([(numCell, numCellFuncVar), (cells', cellFuncVar'), (refCellName, refCellFuncVar), (meshPosFuncAtom, meshPosFuncVar), eqFuncPair, nEqFuncPair],
	   [numCellFun, cellsFun', cellsFun, refCellFunc, meshPosFunc, eqFunc, nEqFunc],
	   [(cells, makeCells, cellTypeFuncTy), (meshInsideAtom, meshInsideFunc, meshInsideType), invalidReplace, invalidCellReplace],
	   [eqFuncPair, nEqFuncPair])
	 end

    (*end case*))



fun makeGeometryFuncsWrap(env, cxt, span, meshData, geometry, inverseInfo, forwardInfo, makeRefCellInsideFunc) =
    let
     val FT.Mesh(mesh) = meshData
     val dim = FT.meshDim mesh
     (*check points exist*)
     (*check correct dim wrappers exist*)
     val test1 = List.exists (fn CF.Points(_) => true | _ => false ) geometry
     val test2 = List.exists (fn CF.Higher(dim', _) => dim'  = dim - 1 | _ => false) geometry
     val test3 = if dim = 2
		 then List.exists (fn CF.LineParam(_) => true | _ => false) geometry
		 else
		  if dim = 3
		  then List.exists (fn CF.PlaneParam(_) => true | _ => false) geometry
		  else (false)
     val test = test3 andalso test2 andalso test1
    in
     if test
     then FemGeometry.makeGeometryFuncs(env, cxt, span, meshData, geometry, inverseInfo, forwardInfo, makeRefCellInsideFunc)
	  handle exn => raise exn
     else {ref = ([],[],[]), pos = ([],[],[]), cell = ([],[],[])}
    end

fun makeMeshPos(env, cxt, span, meshData, transformVars, dumb, extraFuns, extraReplace) = let
 val results = []
 val meshName = FT.nameOf meshData
 val FT.Mesh(mesh) = meshData
 val dim = FT.meshDim mesh
 val vecTy = Ty.vecTy dim
 val cellData = FT.cellOf meshData
 val cellName = FT.envNameOf cellData
 val cellTy = Ty.T_Fem(cellData, SOME(meshName))
 val meshPosData = FT.MeshPos(mesh)
 val meshPosTy = Ty.T_Fem(meshPosData, SOME(meshName))
 val meshPosName = FT.envNameOf meshPosData
 val firstFunc::_ = transformVars

 fun femOpt opt = (opt, meshPosData)

 fun valid([v]) = makePrim'(BV.neq_ii, [AST.E_ExtractFemItem(v, Ty.T_Int, femOpt FO.CellIndex), AST.E_Lit(Literal.intLit (~1))], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
 val results = makeFunctionRegistration(FemName.valid, [meshPosTy], Ty.T_Bool, valid)::results

 fun cell([v]) =
     let
      val cellInt = AST.E_ExtractFemItem(v, Ty.T_Int, femOpt FO.CellIndex)
      val mesh = AST.E_ExtractFem(v, meshData)
     in
      AST.E_LoadFem(cellData, SOME(mesh), SOME(cellInt)) 
     end 
 val results = makeFunctionRegistration(FemName.cell, [meshPosTy], cellTy, cell)::results

 fun refPos([v]) = AST.E_ExtractFemItem(v, vecTy, femOpt FO.RefPos) (*go for it*)
 val results = makeFunctionRegistration(FemName.refPos, [meshPosTy], vecTy, refPos)::results

 val hiddenFuns = results
 val methods = [(fn (x,y,z) => (x,y)) dumb]
 val constants = []
(*TODO: FIX: THE NAMING OF FUNCTIONS*)
(*TODO: BUILD THE TRANSFORM FUNCTION AND DERIVATIVES IN NEWTON STYLE*)

		   
		   
in
 Env.insertNamedType(env, cxt, meshPosName, meshPosTy, constants, extraFuns@methods, extraReplace@hiddenFuns)
end

fun makeDescendentFemTypes (cxt, tyName, span) geometry cellAccData (env, femType) =
    let
     val constants = []
     val methods = []
     val methodRefs = []
     val name = FT.nameOf femType


    in
     (case femType
       of FT.Mesh(m) =>
	  let
	   val meshData = femType
	   val meshName = FT.nameOf femType
	   val meshTy = Ty.T_Fem(femType, NONE)
	   val mapDim = FT.meshDim m
	   val dimSizeVec = Ty.vecTy mapDim
	   val meshMapDim = FT.meshMapDim m
	   val shape = FT.dataShapeOf femType

	   val cellData = FT.cellOf femType			  
	   val cellName = FT.envNameOf (cellData)
	   val cellTy = Ty.T_Fem(cellData, SOME(meshName))
				
	   val refCellData = FT.RefCell(m)
	   val refCellName = FT.envNameOf refCellData
	   val refCellTy = Ty.T_Fem(refCellData, SOME(meshName))
	   val refCellInfo = FT.refCell m

	   val posData = FT.MeshPos(m)
	   val posTy = Ty.T_Fem(posData, SOME(meshName))

	   val transformDofs = Atom.atom "transformDofs"
	   val transformDofsTy = Ty.T_Tensor(Ty.Shape((List.map Ty.DimConst shape))) (*careful*)
	   val transformsDofsFunTy = Ty.T_Fun([cellTy], transformDofsTy)
	   fun makeTransformDofs vars = (case vars
					  of [v] =>
					     let
					      val getMesh =  AST.E_ExtractFem(v, femType)
					      val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, FT.MeshCell(m)))
					      val getTensor = AST.E_ExtractFemItemN([getMesh, getMesh ,getCellInt],
										    [Ty.T_Fem(femType, NONE), Ty.T_Fem(femType, NONE), Ty.T_Int],
										    transformDofsTy, (FemOpt.ExtractDofs, FT.Mesh(m)), NONE)
					     in
					      getTensor
					     end
					   | _ => raise Fail "impossible argument to method got passed type checking"
					(* end case*))

	   val transform = Atom.atom FemName.transform
	   val transform' = Atom.atom FemName.transform'
	   val dimConst = Ty.DimConst(mapDim)
	   val inf = Ty.DiffConst(NONE)
	   val transformFieldTy = Ty.T_Field({diff=inf, dim = dimConst, shape=Ty.Shape([dimConst])})
	   val transformFieldFunc = Ty.T_Fun([cellTy],transformFieldTy)
	   fun makeTransformFieldFunc vars = (case vars
					       of [v] =>
						  let
						   val getMesh = AST.E_ExtractFem(v, femType)
						   val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, FT.MeshCell(m)))
						  in
						   AST.E_FemField(getMesh, getMesh, SOME(getCellInt), transformFieldTy, FemOpt.Transform, NONE)
						  end
						| _ => raise Fail "impossible argument to method got passed type checking"
					     (* end case*))




	   val refCellInside = Atom.atom (FemName.refCellIsInside)
	   val refCellInsideTy = Ty.T_Fun([refCellTy, dimSizeVec], Ty.T_Bool)
	   val FemData.RefCellData({ty=refCellClass, eps, newtonControl,insideInsert,...}) = refCellInfo (*move this outside? Eh.*)
	   fun makeRefCellInsideFunc vars = (case vars
					      of [v1,v2] => makeRefCellInside(env, cxt, span, mapDim, refCellInfo, v2, v1, insideInsert, meshData, meshTy)
					       | _ => raise Fail "impossble got through type checker"
					    (* end case*))

	   fun mkVar (var) = AST.E_Var(var,span)
	   val (newtonResult) = let 
	    val hiddenFuncAtom = (FT.functionNameMake femType (FemName.hiddenNewtonInverse)) (*TODO: inconsistent scheme for hiding this...*)



	    val {contraction, itters, newtonTol, killAfterTol, start} = newtonControl

	    val meshData = m
	    val posTy = Ty.vecTy mapDim
	    val newtonTol = AST.E_Lit(Literal.Real(newtonTol)) (*TODO fix me with ref cell info*)
	    val newtonAttempts = AST.E_Lit(Literal.intLit itters) (*TODO fix me with ref cell info*)


	    val posParam = Var.new (Atom.atom "pos", span, AST.FunParam, posTy)
	    val cellIntParam = Var.new (Atom.atom "cellInt", span, AST.FunParam, Ty.T_Int)
	    val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)
	    val meshPosTy = Ty.T_Fem(FT.MeshPos(m), SOME(FT.nameOf (FT.Mesh(m))))
	    val funTy = Ty.T_Fun([Ty.vecTy mapDim, Ty.T_Int, meshTy], meshPosTy)
	    val hiddenFuncVar = Var.new (hiddenFuncAtom, span, AST.FunVar, funTy)
					
	    val posExpr = mkVar posParam
	    val cellIntExpr = mkVar cellIntParam
	    val meshExpr = mkVar meshParam
				 
	    val insideFunc = makeRefCellInsideFunc

	    val body = makeNewtonInversesBody(env, cxt, span, refCellClass, meshData, newtonTol, newtonAttempts, contraction, killAfterTol, posExpr, cellIntExpr, meshExpr, insideFunc, start, meshTy)
	    val newtonFunc = AST.D_Func(hiddenFuncVar, [posParam, cellIntParam, meshParam], body)


				       
	   in
	    ((hiddenFuncAtom, hiddenFuncVar, newtonFunc))
	   end

	   val (hiddenFuncAtom, hiddenFuncVar, newtonFunc) = newtonResult

	   val (invTransformSpec, invTransformSpec') =
	       let
		val invTransformField = Atom.atom FemName.invTransform
		val invTransformField' = Atom.atom FemName.invTransform'
		val invTransformFieldFuncTy = Ty.T_Fun([cellTy],transformFieldTy)
		fun makeInvTransformFieldFunc vars = (case vars
						       of [v] =>
							  let
							   val getMesh = AST.E_ExtractFem(v, femType)
							   val getCellInt = AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, FT.MeshCell(m)))
							  in
							   AST.E_FemField(getMesh, getMesh, SOME(getCellInt), transformFieldTy, FemOpt.InvTransform, SOME(hiddenFuncVar, span))
							  end
							| _ => raise Fail "impossible for this to be get past typechecker"
						     )
	       in
		((invTransformField, makeInvTransformFieldFunc, invTransformFieldFuncTy),
		 (invTransformField', makeInvTransformFieldFunc, invTransformFieldFuncTy))
	       end

	   local
	    val insideCellAtom = Atom.atom (FemName.isInsideMeshCell)
	    val insideCellFunTy = Ty.T_Fun([cellTy, dimSizeVec], Ty.T_Bool)
	    fun makeInsideCell([v1,v2]) =
		let
		 val getPos = AST.E_Apply((hiddenFuncVar,span), [v2, AST.E_ExtractFemItem(v1, Ty.T_Int, (FemOpt.CellIndex, cellData)), AST.E_ExtractFem(v1, meshData)], posTy)
		 val getValid = makePrim'(BV.neq_ii, [AST.E_ExtractFemItem(getPos, Ty.T_Int, (FemOpt.CellIndex, posData)), AST.E_Lit(Literal.Real(RealLit.m_one))], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
		in
		 getValid
		end
	      | makeInsideCell _ = raise Fail "typechecker error"
	   in
	   val insideCellTriple = (insideCellAtom, makeInsideCell, insideCellFunTy)
	   end

	   local
	    val isValidAtom = Atom.atom (FemName.isValidCell)
	    val isValidFunTy = Ty.T_Fun([cellTy], Ty.T_Bool)
	    fun makeValidCheck([v1]) =
		let
		 val int =  AST.E_ExtractFemItem(v1, Ty.T_Int, (FemOpt.CellIndex, cellData))
		 val test = makePrim'(BV.gte_ii, [int, AST.E_Lit(Literal.intLit 0)], [Ty.T_Int, Ty.T_Int], Ty.T_Bool)
		in
		 test
		end
	      | makeValidCheck _ = raise Fail "typechecker error"
	   in
	   val validCellTriple = (isValidAtom, makeValidCheck, isValidFunTy)
	   end

	   local
	    val refMeshPosAtom = Atom.atom (FemName.refMeshPos)
	    val refMeshFunTy = Ty.T_Fun([cellTy, dimSizeVec], posTy)
	    fun refMeshPosBuild([v1, v2]) =
		let
		 val mesh = AST.E_ExtractFem(v1, meshData)
		 val int = AST.E_ExtractFemItem(v1, Ty.T_Int, (FemOpt.CellIndex, cellData))
		 val intNeg = AST.E_Lit(Literal.intLit (~1))
		 val result = AST.E_ExtractFemItemN([mesh, int, v2, intNeg], [meshTy, Ty.T_Int, dimSizeVec, Ty.T_Int], posTy, (FemOpt.RefBuild, meshData), NONE)
		in
		 result
		end
	      | refMeshPosBuild _ = raise Fail "typechecker error"
	   in
	   val refMeshPosTriple = (refMeshPosAtom, refMeshPosBuild, refMeshFunTy)
	   end
	   val transformSpec = (transform, makeTransformFieldFunc, transformFieldFunc)
	   val transformSpec' = (transform', makeTransformFieldFunc, transformFieldFunc)
	   val transformFuncs = [(transformDofs, makeTransformDofs, transformsDofsFunTy), transformSpec, transformSpec',
				 invTransformSpec, invTransformSpec', insideCellTriple,validCellTriple, refMeshPosTriple]
	   val cellMethods = [(hiddenFuncAtom, hiddenFuncVar)]
	   (* TODO: The naming in this area of the compiler is really inconsistent*)

	   val (dumb, dnTFuncs) = dnT(env, cxt, span, FT.Mesh(m),  3) (*maybe more?*)
	   val vars = List.map (fn (x,y,z) => y) dnTFuncs
	   val dTFuncDcls = List.map (fn (x,y,z) => z) dnTFuncs
	   val dumbDef = (fn (x,y,z) => z) dumb
	   val refCellTransformFuncs = [(refCellInside, makeRefCellInsideFunc, refCellInsideTy)] (*pass to geometry.....yay!*)
	   val geometryInfo = makeGeometryFuncsWrap(env, cxt, span, meshData, geometry, invTransformSpec, transformSpec, makeRefCellInsideFunc)

	   val (extraCellFuncInfo, extraCellFuncDcls, extraCellReplaceFuncs) = #cell geometryInfo

	   local
	    fun cellInt x = AST.E_ExtractFemItem(x, Ty.T_Int, (FemOpt.CellIndex, FT.Mesh(m)))
	    fun meshVal x=  AST.E_ExtractFem(x, FT.Mesh(m))
	    fun makeOne(name, ty, discard) =
		let
		 val funTy = Ty.T_Fun([cellTy], ty)
		 fun doit([v1]) = AST.E_ExtractFemItem2(meshVal v1, cellInt v1, Ty.T_Int, ty, (FemOpt.CellData name, FT.Mesh(m)))
		   | doit _ = raise Fail "typechecker error"


		in
		 (name, doit, funTy)
		end
	    val funs = List.map makeOne cellAccData
	   in
	   val dataAccReplaceFuns = funs
	   end



	   val env' = Env.insertNamedType(env, cxt, cellName, cellTy, constants, extraCellFuncInfo@cellMethods, extraCellReplaceFuncs@transformFuncs@dataAccReplaceFuns)


	   val (refActualFuncI, refActualFuncDcls, refReplaceFunc) = #ref geometryInfo
	   val refCellTransformFuncs' = refReplaceFunc@refCellTransformFuncs


	   val env'' = Env.insertNamedType(env', cxt, refCellName, refCellTy, constants, refActualFuncI, refCellTransformFuncs')

	   (*get the appropraite transform vars*)
	   val (posActualFuncI , posActualFuncDcls, posReplaceFunc) = #pos geometryInfo

	   val env''' = makeMeshPos(env, cxt, span, FT.Mesh(m), vars, dumb, posActualFuncI, posReplaceFunc)
				   
	  in
	   (env''', extraCellFuncDcls@posActualFuncDcls@refActualFuncDcls@(dumbDef::newtonFunc::dTFuncDcls))
	  end
	| FT.Func(f) =>
	  let
	   val funcName = FT.nameOf femType
	   val cellName = FT.envNameOf (FT.FuncCell(f))
	   val SOME(space) = FT.dependencyOf femType
	   val mesh = FT.meshOf femType

	   val FT.Space(spaceData) = space
	   val FT.Mesh(meshData) = mesh
	   val meshCellEnvName = FT.envNameOf (FT.MeshCell(meshData))
	   val mapDim = FT.meshDim meshData

	   val funcCellData = FT.cellOf femType
	   val cellTy = Ty.T_Fem(funcCellData, SOME(funcName))
	   val funcDataShape = FT.dataShapeOf femType
	   val funcRangeShape = FT.dataRangeShapeOf femType
	   val dofCount = FT.spaceDim spaceData
	   val dofIndexSeq = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst(dofCount)))

	   val functionDofs = Atom.atom FN.funcDofs
	   val funcDofsTy = Ty.T_Tensor(Ty.Shape (List.map Ty.DimConst funcDataShape))
	   val functionDofsFunTy = Ty.T_Fun([cellTy], funcDofsTy)
	   fun makeFunctionDofs vars = (case vars
					 of [v] =>
					    let
					     val getFunc = AST.E_ExtractFem(v, femType)
					     val depFemType = Option.valOf (FemData.dependencyOf femType)
					     val getSpace = AST.E_ExtractFem(getFunc, depFemType)
					     (* val getSpace = AST.E_ExtractFem(getFunc,space) *)
					     (*might want to extract it now.*)
					     val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
					     (* val getIndexi = AST.E_ExtractFemItem2(getSpace, (getCellInt), Ty.T_Int, dofIndexSeq, (FemOpt.ExtractIndices, space)) *)
										 (* val getTensor = AST.E_ExtractFemItem2(getFunc, (getIndexi), dofIndexSeq, funcDofsTy, (FemOpt.ExtractDofsSeq, femType)) *)
					     val getTensor = AST.E_ExtractFemItemN([getFunc, getSpace, getCellInt],
										   [Ty.T_Fem(femType, NONE), Ty.T_Fem(depFemType, NONE), Ty.T_Int],
										   funcDofsTy, (FemOpt.ExtractDofs, femType), NONE)
										  
					    in
					     getTensor
					    end
					  | _ => raise Fail "impossible argument to method got passed type checking")


	   val refField = Atom.atom FemName.refField
	   val refField' = Atom.atom FemName.refField'
	   val dimConst = Ty.DimConst(mapDim)
	   val rangeShape = Ty.Shape((List.map Ty.DimConst funcRangeShape))
	   val inf = Ty.DiffConst(NONE)
	   val refFieldTy = Ty.T_Field({diff = inf, dim = dimConst, shape=rangeShape})
	   val refFieldFunc =  Ty.T_Fun([cellTy], refFieldTy)
	   fun makeRefFieldFunc vars = (case vars
					 of [v] =>
					    let
					     val getFunc = AST.E_ExtractFem(v, femType)
					     val getSpace = AST.E_ExtractFem(getFunc, space)
					     val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
					    in
					     AST.E_FemField(getFunc,getSpace, SOME(getCellInt),refFieldTy,FO.RefField, NONE )
					    end
					  | _ => raise Fail "impossible argument to method got passed type checking")

	   val transformRefFieldFunc = Atom.atom (FemName.trf)
	   val transformRefFieldFunc' = Atom.atom (FemName.trf')
	   val transformFieldTy = Ty.T_Field({diff = inf, dim = dimConst, shape = Ty.Shape([dimConst])})
	   fun makeTransformRefFieldFunc vars = (case vars
						  of [v] =>
						     let
						      val getFunc = AST.E_ExtractFem(v, femType)
						      val getSpace = AST.E_ExtractFem(getFunc, space)
						      val getMesh = AST.E_ExtractFem(getSpace, mesh)
						      val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
						      (*get mesh cell env, get the function, make the fields*)
						      (* val SOME(meshCellEnv) = Env.findTypeEnv(env, meshCellEnvName) *)
						      (* val SOME(newtonMethod) = TypeEnv.findMethod(meshCellEnv, FT.functionNameMake mesh (FemName.hiddenNewtonInverse )) *)
						      (*get the basis*)
						      val field1 = AST.E_FemField(getFunc,getSpace, SOME(getCellInt), refFieldTy, FO.RefField, NONE )
						      val field2 = AST.E_FemField(getMesh, getMesh, SOME(getCellInt), transformFieldTy, FemOpt.InvTransform, NONE)
						      val field = makePrim'(BV.comp, [field1, field2], [refFieldTy, transformFieldTy], refFieldTy)
						     in
						      field
						     end
						   | _ => raise Fail "impossible that this got past the typechecker"
						)

	   val finvt = Atom.atom (FemName.finvt)
	   fun makeFinvtFunc vars = (case vars
				      of [v] =>
					 let
					  val getFunc = AST.E_ExtractFem(v, femType)
					  val getSpace = AST.E_ExtractFem(getFunc, space)
					  val getMesh = AST.E_ExtractFem(getSpace, mesh)
					  val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
					  (*get mesh cell env, get the function, make the fields*)
					  val SOME(meshCellEnv) = Env.findTypeEnv(env, meshCellEnvName)
					  val SOME(newtonMethod) = TypeEnv.findMethod(meshCellEnv, FT.functionNameMake mesh (FemName.hiddenNewtonInverse))
					  (*get the basis*)
					  val field1 = AST.E_FemField(getFunc,getSpace, SOME(getCellInt), refFieldTy, FO.RefField, NONE )
					  val field2 = AST.E_FemField(getMesh, getMesh, SOME(getCellInt), transformFieldTy, FemOpt.InvTransform, SOME(newtonMethod, span))
					  val field = makePrim'(BV.comp, [field1, field2], [refFieldTy, transformFieldTy], refFieldTy)
					 in
					  field
					 end
				       | _ => raise Fail "impossible that this got past the typechecker"
				    )


	   val transformFuncs = [(functionDofs, makeFunctionDofs, functionDofsFunTy),
				 (refField, makeRefFieldFunc, refFieldFunc), (refField', makeRefFieldFunc, refFieldFunc),
				 (transformRefFieldFunc, makeTransformRefFieldFunc, refFieldFunc),
				 (transformRefFieldFunc', makeTransformRefFieldFunc, refFieldFunc),
				 (finvt, makeFinvtFunc, refFieldFunc)
				]
	   val env' = Env.insertNamedType(env, cxt, cellName, cellTy, constants, methods, transformFuncs)
	  in
	   (env', [])
	  end

	| _ => (env, [])
     (* end case *))
    end
      

end
  
