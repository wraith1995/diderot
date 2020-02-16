(* gen-tys-and-ops.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *)

structure GenTysAndOps : sig

    val gen : CodeGenEnv.t * CollectInfo.t -> {
            preWorld : CLang.decl list,         (* decls to place before world outside of namespace *)
            postWorld : CLang.decl list         (* decls to place after world inside of namespace *)
          }

  end = struct

    structure IR = TreeIR
    structure Ty = TreeTypes
    structure CL = CLang
    structure RN = CxxNames
    structure Env = CodeGenEnv
    structure FT = FemData

    val zero = RealLit.zero false

    fun mkReturn exp = CL.mkReturn(SOME exp)
    fun mkInt i = CL.mkInt(IntInf.fromInt i)
    fun mkFunc (ty, name, params, body) = CL.D_Func(["inline"], ty, [], name, params, body)
   (* make a constructor function's prototype and out-of-line definition *)
    fun mkConstr (cls, params, inits) = (
            CL.D_Constr([], [], cls, params, NONE),
            CL.D_Constr(["inline"], [CL.SC_Type(CL.T_Named cls)], cls, params, SOME(inits, CL.mkBlock[]))
          )
   (* make a member function prototype and out-of-line definition *)
    fun mkMemberFn (cls, ty, f, params, body) = (
            CL.D_Proto([], ty, f, params),
            CL.D_Func(["inline"], ty, [CL.SC_Type(CL.T_Named cls)], f, params, body)
          )

    fun trType namespace env ty = TypeToCxx.trQType(env, namespace, ty)

    fun sequenceTy (elemTy, sz) =
          CL.T_Template("diderot::array", [elemTy, CL.T_Named(Int.toString sz)])

    fun genTyDecl env = let
          val (realTy, realTyName, realTySz) = if #double(Env.target env)
					       then (CL.double, "double", 8)
					       else (CL.float, "float", 4)
	  val (intTy, intTyName, intTySize) = if #longint(Env.target env)
					      then (CL.int64, "int64_t", 8)
					      else (CL.int32, "int32_t", 4)
	  (*generate parse function for tuple type*)
	  fun generateTupleTyParse(tys, name) =
	      let
	       val tplname = "tpl_" ^ name
	       val parsename = tplname ^ "_parse"
	       val destoryname = tplname ^ "_destory"
	       val paramsParse = [
		CL.PARAM([], CL.voidPtr, "ptr"),
		CL.PARAM([], CL.charPtr, "str"),
		CL.PARAM([], CL.charPtr, "err")
	       ]
	       val paramsDestory = [CL.PARAM([], CL.voidPtr, "ptr")]
				     (*need to make sure input types properly converted and generated!*)
				     (*assume [] are only on teem types*)
				     (*init struct*)
				     (*create prase stuff*)
				     (*try to parse each one and check*)
				     (*if fail -> parse and use error*)
				     (*if you get a nrrd, we are pasing strings*)
				     (*build struct*)

				     
	      in
	      end
						      
          (* generate the type and member declarations for a recorded type *)

	  fun generateFemType(ty) =
	      let
	       val name = Atom.toString (FT.nameOf ty)
	       val underDim = FT.underlyingDim ty
	      in
	      (case ty
		of FT.Mesh(m) => (
		 let
		  val FemData.RefCellData({insideInsert=insideTest, ...}) =  (FemData.refCell m)
		  
		  val indexMap = CL.D_Var([],CL.T_Ptr(intTy),[], "indexMap", NONE)
		  val coordinatesMap = CL.D_Var([],CL.T_Ptr(realTy),[], "coordMap", NONE)
		  val dimConst = CL.D_Var([], intTy,[], "dim", NONE)
		  val mapDimConst = CL.D_Var([], intTy,[], "mapDim", NONE)
		  val numCells = CL.D_Var([], intTy,[], "numCells", NONE)
		  val indexPtr = CL.D_Var([], CL.voidPtr,[], "index", NONE)
		  val conTy = CL.T_Ptr(intTy)
		  val con = CL.D_Var([], conTy, [], "con", NONE)

		  val extra =
		      let
		       val results = [con]
		       val results = if (Option.isSome (FemData.meshAccInsert(m)))
				     then indexPtr::results
				     else indexPtr::results
		      in
		       results
		      end

		  local
		   val extraData = FT.getCellDataTypes m
		   fun sortData ss = ListMergeSort.sort (fn ((a1, s1), (a2, s2)) => (Atom.toString a1) > (Atom.toString a2)) ss;
		   val extraData' = sortData extraData
		   fun convertTy(FT.String) = CL.charPtr
		     | convertTy(FT.Bool) = CL.boolTy
		     | convertTy (FT.Int) = intTy
		     | convertTy (FT.Real) = realTy
		     | convertTy (FT.Tensor(_)) = realTy
		     | convertTy (FT.Array(ty, _)) = convertTy ty
		   fun makeParam (a, ty) = CL.D_Var([], CL.T_Ptr(convertTy ty), [], Atom.toString a, NONE)
		  in
		  val extrasParams = List.map makeParam extraData'
		  end


		  (* We need to add a load function via string*)
		  (*load the file*)
		  (*ugg*)
		  (*val load = CL.D_Func([],) *)
		  val block = CL.S_Block([CL.S_Comment(["No something with the nrrd"])])
		  val load = CL.mkFuncDcl(CL.T_Named(name),"operator=", [CL.PARAM([],CL.T_Named("std::string"), "file")], block)

		  val insideTests = (case insideTest
				      of SOME(file) =>
					 let
					  val ins = TextIO.openIn file
					  (*this should be centrealized somewhere:*)
					  fun loop ins = 
					      case TextIO.inputLine ins of 
						  SOME line => line :: loop ins 
						| NONE      => []
					  val lines = loop ins before TextIO.closeIn ins
					  val insert = String.concatWith "\n" lines
					 in
					  [CL.D_Verbatim(["bool ", name, "_inside(", realTyName, " * x, ", realTyName, " eps", "){\n", insert, "\n}"])]
					 end
				       | NONE => [])


		  val fields = [indexMap, coordinatesMap, dimConst, mapDimConst, numCells]@extra@extrasParams@[load]


														(*extrac*)

		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []} (*add public declerations*)
		 in
		  (class, insideTests)
		 end)
		| FT.Space(s) => (
		 let
		  val spaceTy = CL.T_Named(name)
		  val indexMap = CL.D_Var([],CL.T_Ptr(intTy),[], "indexMap", NONE)
		  val meshName = Atom.toString (FT.nameOf (FT.meshOf ty))
		  val mesh = CL.D_Var([],CL.T_Named(meshName),[], "mesh", NONE)
		  val block = CL.S_Block([CL.S_Comment(["No something with the nrrd"])])
		  val load = CL.mkFuncDcl(CL.T_Named(name),"operator=", [CL.PARAM([], CL.T_Named("std::string"), "file")], block)

		  val loadFemBlock = CL.S_Block([CL.S_Decl([],spaceTy, "space", SOME(CL.I_Exp (CL.E_UnOp(CL.%*, CL.E_Var("this"))))),
						 CL.S_Exp(CL.mkAssignOp(CL.E_Select(CL.mkVar "space", "mesh"), CL.$=, CL.mkVar "mesh")),
						 CL.mkReturn (SOME(CL.mkVar "space"))])
		  val loadFem = CL.mkFuncDcl(CL.T_Named(name),"loadFem", [CL.PARAM([],CL.T_Named(meshName), "mesh")], loadFemBlock)
					 
		  val fields = [indexMap, mesh, load, loadFem]
		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []}
		 in
		  
		  (class, [])
		 end
		)
		| FT.Func(f) => (
		 let
		  val funcTy = CL.T_Named(name)
		  val coords = CL.D_Var([],CL.T_Ptr(realTy),[], "coordMap", NONE)
		  val spaceName = Atom.toString (FT.nameOf (FT.spaceOf ty))
		  val space =  CL.D_Var([],CL.T_Named(spaceName),[], "space", NONE)
		  val block = CL.S_Block([CL.S_Comment(["No something with the nrrd"])])
		  val load = CL.mkFuncDcl(CL.T_Named(name),"operator=", [CL.PARAM([],CL.T_Named("std::string"), "file")], block)
		  val loadFemBlock = CL.S_Block([CL.S_Decl([],funcTy, "func", SOME(CL.I_Exp (CL.E_UnOp(CL.%*, CL.E_Var("this"))))),
						 CL.S_Exp(CL.mkAssignOp(CL.E_Select(CL.mkVar "func", "space"), CL.$=, CL.mkVar "space")),
						 CL.mkReturn (SOME(CL.mkVar "func"))])
		  val loadFem = CL.mkFuncDcl(CL.T_Named(name),"loadFem", [CL.PARAM([],CL.T_Named(spaceName), "space")], loadFemBlock)
					 
		  val fields = [coords, space, load, loadFem]

		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []}
		 in
		  (class, [])
		 end
		)
		| FT.MeshCell(mesh) => (
		 let
		  val meshCellTy = CL.T_Named(name)
		  val meshCellRefTy = CL.T_Named(name ^ " & ")
		  val int = CL.D_Var([], intTy ,[], "cell", NONE)
		  val meshName = Atom.toString (FT.nameOf (FT.meshOf ty))
		  val meshTy = CL.T_Named(meshName)
		  val meshPtrTy = CL.T_Ptr(meshTy)
		  val meshPtr =  CL.D_Var([],meshTy,[], "mesh", NONE)
		  val osty = CL.T_Named("std::ostream &") (*attr doesn't support ty & *)
		  val printer = CL.mkFuncDcl(osty, "operator<<", [CL.PARAM([], osty, "os"),
								  CL.PARAM(["const"], meshCellRefTy, "cell")], CL.S_Block([CL.S_Return(SOME(CL.E_BinOp(CL.E_Var("os"),
																		CL.#<<,
																		   CL.E_Select(CL.E_Var("cell"),"cell"))))]))
		  val printerPrime = CL.mkFuncDcl(osty, "operator<<", [CL.PARAM([], osty, "os")] ,  CL.S_Block([CL.S_Return(SOME(CL.E_BinOp(CL.E_Var("os"),
																		CL.#<<,
																		   CL.mkIndirect(CL.E_Var("this"),"cell"))))]))
		  val builder = CL.mkFuncDcl(meshCellTy, "makeFem",
					     [CL.PARAM([], meshTy, "mesh"),
					      CL.PARAM([], intTy, "cellInt")],
					     CL.S_Block([
							CL.S_Decl([], meshCellTy, "cell", NONE),
							CL.S_Exp(
							 CL.E_AssignOp(CL.E_Select(CL.E_Var("cell"), "cell"), CL.$=, CL.E_Var("cellInt"))),
							CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var("cell"), "mesh"), CL.$=,
															 (CL.E_Var("mesh")))),
							CL.S_Return(SOME(CL.E_Var("cell")))
					    ]))

		  val copy_to = CL.mkFuncDcl(CL.charPtr, "copy_to", [CL.PARAM(["const"],meshCellTy, "dumb"), CL.PARAM([], CL.charPtr, "cp")],
					     CL.S_Block([
							CL.S_Decl([], CL.size_t,
								  "nbytes",
								  SOME(CL.I_Exp(CL.E_Sizeof(intTy)))),
							CL.mkCallExp(CL.E_Var("std::memcpy"),
								     [
								       CL.E_Var("cp"),
								       CL.mkAddrOf (CL.E_Grp((CL.E_Select(CL.E_Var("dumb"), "cell")))),
								       CL.E_Var("nbytes")
								    ]),
							CL.S_Return(SOME(CL.mkBinOp(CL.E_Var("cp"), CL.#+, CL.E_Var("nbytes"))))
					    ]))

		  val fields = [int, meshPtr, printerPrime]
		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []}
		 in
		  (class, [printer, builder, copy_to])
		 end
		)
		| FT.FuncCell(func) => (
		 let
		  val funcCellTy = CL.T_Named(name)
		  val funcCellRefTy = CL.T_Named(name ^ " & ")
		  val int = CL.D_Var([], intTy ,[], "cell", NONE)
		  val funcName = Atom.toString ((FT.nameOf o Option.valOf o FT.dependencyOf) ty)
		  val funcTy = CL.T_Named(funcName)
		  val funcPtrTy = CL.T_Ptr(funcTy)
		  val funcPtr = CL.D_Var([],funcTy,[], "func", NONE)
					
		  val osty = CL.T_Named("std::ostream &") (*attr doesn't support ty & *)
		  val printer = CL.mkFuncDcl(osty, "operator<<", [CL.PARAM([], osty, "os"),
								  CL.PARAM(["const"], funcCellRefTy, "cell")],
					     CL.S_Block(
					      [CL.S_Return(SOME(CL.E_BinOp(CL.E_Var("os"),CL.#<<, CL.E_Select(CL.E_Var("cell"),"cell"))))]))

		  val builder = CL.mkFuncDcl(funcCellTy, "makeFem",
					     [CL.PARAM([], funcTy, "func"),
					      CL.PARAM([], intTy, "cellInt")],
					     CL.S_Block([
							CL.S_Decl([], funcCellTy, "cell", NONE),
							CL.S_Exp(
							 CL.E_AssignOp(CL.E_Select(CL.E_Var("cell"), "cell"), CL.$=, CL.E_Var("cellInt"))),
							CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var("cell"), "func"), CL.$=,(CL.E_Var("func")))),
							CL.S_Return(SOME(CL.E_Var("cell")))
					    ]))
		  val copy_to = CL.mkFuncDcl(CL.charPtr, "copy_to", [CL.PARAM(["const"], funcCellTy, "dumb"), CL.PARAM([], CL.charPtr, "cp")],
					     CL.S_Block([
							CL.S_Decl([], CL.size_t,
								  "nbytes",
								  SOME(CL.I_Exp(CL.E_Sizeof(intTy)))),
							CL.mkCallExp(CL.E_Var("std::memcpy"),
								     [
								       CL.E_Var("cp"),
								       CL.mkAddrOf (CL.E_Grp((CL.E_Select(CL.E_Var("dumb"), "cell")))),
								       CL.E_Var("nbytes")
								    ]),
							CL.S_Return(SOME(CL.mkBinOp(CL.E_Var("cp"), CL.#+, CL.E_Var("nbytes"))))
					    ]))
		  val fields = [int, funcPtr]
		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []}
		 in
		  (class, [printer, builder, copy_to])
		 end)
		| FT.MeshPos(mesh) => (
		 let
		  val meshPosTy = CL.T_Named(name)
		  val meshPosRefTy = CL.T_Named(name ^ " & ")
		  val meshPosPointerTy = CL.T_Ptr meshPosTy
		  val meshName = Atom.toString (FT.nameOf (FT.meshOf ty))
		  val meshTy = CL.T_Named(meshName);

		  val dim = FT.meshDim mesh
		  val tensorTy = RN.tensorRefTy [dim] (*TODO: maybe tensor ref or vec???? - figure out later.*)
		  val tensorDataTy = RN.tensorTy [dim] (*to figure out sizing for copies*)

		  val blockedTensorDataTy = RN.tensorRefTy [2*dim + 1]
		  (*params:*)
		  val meshParam = CL.D_Var([], meshTy, [], "mesh", NONE)
		  val intParam = CL.D_Var([], intTy, [], "cell", NONE)
		  val refPosParam = CL.D_Var([], tensorDataTy, [], "refPos", NONE)
		  val worldPosParam = CL.D_Var([], tensorDataTy, [], "worldPos", NONE)
		  val worldPosCompute = CL.D_Var([], CL.boolTy, [], "wpc", NONE)
		  val valid = CL.D_Var([], CL.boolTy, [], "valid", NONE)
		  val faceParam = CL.D_Var([], intTy, [], "face", NONE)

		  val osty = CL.T_Named("std::ostream &")



		  (*TODO: fix meshPos printing*)
		  val printer = CL.mkFuncDcl(osty, "operator<<", [CL.PARAM([], osty, "os"), CL.PARAM(["const"], meshPosRefTy, "pos")], 
					     CL.S_Block([
							CL.S_Return(SOME(CL.E_BinOp(CL.E_Var("os"),
										    CL.#<<, CL.E_Select(CL.E_Var("pos"),"cell")
										       )))]))

		  fun assgn src name = CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var(src),name), CL.$=, CL.E_Var(name)))
		  fun assgnTensor src name = CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var(src),name), CL.$=, CL.E_Select(CL.E_Var(name), "_data")))
					   
		  fun init name = (name, CL.I_Exp(CL.E_Var(name)))
		  fun init' name exp = (name, CL.I_Exp(exp))

		  fun noAssgn name = CL.S_Exp(CL.E_AssignOp(CL.E_Var(name), CL.$=, CL.E_Var(name)))

		  fun assgn' name expr = CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var("pos"), name), CL.$=, expr))

		  fun assgnTrue name = CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var("pos"), name), CL.$=, CL.E_Bool(true)))
		  fun assgnFalse name = CL.S_Exp(CL.E_AssignOp(CL.E_Select(CL.E_Var("pos"), name), CL.$=, CL.E_Bool(false)))
					       

		  fun nullPtrStart name = CL.E_Apply(CL.E_Var(name),[CL.E_Var("nullptr")])
		  fun varStart name = CL.E_Apply(CL.E_Var(name), [CL.E_Var(name)])
		  fun expStart name exp = CL.E_Apply(CL.E_Var(name), [exp])

		  fun simpleInit name = (name, CL.I_Exp(CL.E_Var(name)))
		  (*should we check for > 0 or < 0 and length >*)
		  val makeFem = CL.mkFuncDcl(meshPosTy, "makeFem", [CL.PARAM([], meshTy, "mesh"), CL.PARAM([], blockedTensorDataTy, "data")],
					     CL.S_Block([
							CL.S_Decl([], realTy, "cell", SOME(CL.I_Exp(CL.E_Subscript(CL.mkVar "data", CL.E_Int(IntLit.fromInt 0, intTy))))),
							CL.S_Decl([], tensorTy, "ref", SOME(CL.I_Exp(CL.mkAddrOf  (CL.E_Subscript(CL.mkVar "data", CL.E_Int(IntLit.fromInt 1, intTy)))))),
							CL.S_Decl([], tensorTy, "world", SOME(CL.I_Exp(CL.mkAddrOf (CL.E_Subscript(CL.mkVar "data", CL.E_Int(IntLit.fromInt 4, intTy)))))),
							CL.S_Decl([], intTy, "cellInt", SOME(CL.I_Exp(CL.E_XCast("static_cast", intTy, CL.mkVar "cell")))),
							CL.S_Decl([], CL.boolTy, "test", SOME(CL.I_Exp(CL.E_BinOp(CL.E_Int(IntLit.fromInt 0, intTy),
															    CL.#<,
															       CL.mkVar "cellInt")))),
							CL.S_Return(SOME(CL.E_Apply(CL.mkVar "allBuild", [CL.mkVar "mesh", CL.mkVar "cellInt", CL.mkVar "ref", CL.mkVar "world", CL.mkVar "test", CL.mkVar "test" ])))
							

							
							
					    ]))


		  val allBuild = CL.mkFuncDcl(meshPosTy, "allBuild", [CL.PARAM([], meshTy, "mesh"),
					       CL.PARAM([], intTy, "cell"),
					       CL.PARAM([], tensorTy, "refPos"),
					       CL.PARAM([], tensorTy, "worldPos"),
					       CL.PARAM([], CL.boolTy, "wpc"),
					       CL.PARAM([], CL.boolTy, "valid")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 assgn "pos" "cell",
							 
							 assgnTensor "pos" "refPos",
							 assgnTensor "pos" "worldPos",
							 assgn "pos" "wpc",
							 assgn "pos" "valid",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							])
					     )
		  val allBuild' = CL.mkFuncDcl(meshPosTy, "allBuild", [CL.PARAM([], meshTy, "mesh"),
								      CL.PARAM([], intTy, "cell"),
								      CL.PARAM([], tensorTy, "refPos"),
								      CL.PARAM([], tensorTy, "worldPos"),
								      CL.PARAM([], CL.boolTy, "wpc"),
								      CL.PARAM([], CL.boolTy, "valid"),
								      CL.PARAM([], intTy, "face")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 assgn "pos" "cell",
							 assgn "pos" "face",
							 assgnTensor "pos" "refPos",
							 assgnTensor "pos" "worldPos",
							 assgn "pos" "wpc",
							 assgn "pos" "valid",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							])
					     )
		  val refBuild = CL.mkFuncDcl(meshPosTy, "refBuild", [CL.PARAM([], meshTy, "mesh"),
								      CL.PARAM([], intTy, "cell"),
								      CL.PARAM([], tensorTy, "refPos")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 assgn "pos" "cell",
							 assgn "pos" "refPos",
							 assgnFalse "wpc",
							 assgnTrue "valid",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							]))
		  val refBuild' = CL.mkFuncDcl(meshPosTy, "refBuild", [CL.PARAM([], meshTy, "mesh"),
								      CL.PARAM([], intTy, "cell"),
								      CL.PARAM([], tensorTy, "refPos"),
								      CL.PARAM([], intTy, "face")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 assgn "pos" "cell",
							 assgn "pos" "refPos",
							 assgnFalse "wpc",
							 assgnTrue "valid",
							 assgn "pos" "face",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							]))
		  val invalidBuild = CL.mkFuncDcl(meshPosTy, "invalidBuild", [CL.PARAM([], meshTy, "mesh")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							])
					     )

		  val invalidBuild' = CL.mkFuncDcl(meshPosTy, "invalidBuild", [CL.PARAM([], meshTy, "mesh"),
									       CL.PARAM([], tensorTy, "refPos")],
					      CL.S_Block([
							 CL.S_Decl([], meshPosTy, "pos",NONE),
							 assgn "pos" "mesh",
							 assgn "pos" "refPos",
							 CL.S_Return(SOME(CL.E_Var("pos")))
							])
					     )
		      

					    
		  val genBuild = CL.D_Constr([],[], name,
					      [CL.PARAM([], meshTy, "mesh"),
					       CL.PARAM([], intTy, "cell"),
					       CL.PARAM([], tensorTy, "refPos"),
					       CL.PARAM([], tensorTy, "worldPos"),
					       CL.PARAM([], CL.boolTy, "wpc"),
					       CL.PARAM([], CL.boolTy, "valid"),
					       CL.PARAM([], intTy, "face")],
					      SOME([varStart "mesh",
						    varStart "cell",						    
						    varStart "wpc",
						    varStart "valid",
						    varStart "face"],
						   (CL.S_Block([
							       noAssgn "mesh",
							       noAssgn "cell",
							       noAssgn "refPos",
							       noAssgn "worldPos",
							       noAssgn "wpc",
							       noAssgn "valid",
							       noAssgn "face"])))
					      
					    )

		  val defaultBuild =  CL.D_Constr (
                       [], [], name, [], SOME([
					      expStart "cell" (CL.E_Int(IntLit.fromInt (~1), intTy)),
					      expStart "face" (CL.E_Int(IntLit.fromInt (~1), intTy)),
					      expStart "wpc" (CL.E_Bool(false)),
					      expStart "valid" (CL.E_Bool(false))], CL.S_Block([])))

		  val defaultBuild' =  CL.D_Constr (
                       [], [], name,  [CL.PARAM([], meshTy, "mesh")],
		       SOME([
			    expStart "cell" (CL.E_Int(IntLit.fromInt (~1), intTy)),
			    expStart "face" (CL.E_Int(IntLit.fromInt (~1), intTy)),
			    expStart "wpc" (CL.E_Bool(false)),
			    expStart "valid" (CL.E_Bool(false)),
			    varStart "mesh"], CL.S_Block([])))
						   
		  (* val invalidBuild = CL.D_Constr([],[], name, *)
		  (* 				  [CL.PARAM([], meshTy, "mesh")], *)
		  (* 				  [simpleInit "mesh", *)
		  (* 				   crapInit "cell" (CL.E_Int(IntLit.fromInt (~1), intTy)), *)
		  (* 				   crapInit "valid" (CL.E_Bool(false)), *)
		  (* 				   crapInit "wpc" (CL.E_Bool(false)), *)
		  (* 				   crapInit "refPos" (CL.E_Var("nullPtr")), *)
		  (* 				   crapInit "worldPos" (CL.E_Var("nullPtr")) *)
		  (* 				  ] *)
		  (* 				 ) *)

		  val arrayTy = CL.T_Array(realTy, SOME(dim*2))
		  val thisName = "dumb"
		  val copy_to = CL.D_Func([], CL.charPtr,[], "copy_to ", [CL.PARAM(["const"], meshPosTy, thisName), CL.PARAM([], CL.charPtr, "cp")],
					     CL.S_Block([
							CL.S_Decl([], CL.size_t,
								  "nbytes1",
								  SOME(CL.I_Exp(CL.E_Sizeof(tensorDataTy)))),
							CL.S_Decl([], CL.size_t,
								  "nbytes2",
								  SOME(CL.I_Exp(CL.E_Sizeof(intTy)))),
							
							CL.mkCallExp(CL.E_Var("std::memcpy"),
								     [
								       CL.E_Var("cp"),
								       CL.mkAddrOf (CL.E_Grp((CL.E_Select(CL.E_Var(thisName), "cell")))),
								       CL.E_Var("nbytes2")
								    ]),
							CL.S_Exp(CL.E_AssignOp(CL.E_Var("cp"), CL.+=, CL.E_Var("nbytes2"))),
							CL.mkIfThenElse (
							 CL.E_Select(CL.E_Var(thisName), "valid"),
							 CL.S_Block([
								    CL.mkCallExp(CL.E_Var("std::memcpy"),
										 [
										   CL.E_Var("cp"),
										   CL.E_Select((CL.E_Select(CL.E_Var(thisName), "refPos")), "_data"),
										   CL.E_Var("nbytes1")
										]),
								    CL.S_Exp(CL.E_AssignOp(CL.E_Var("cp"), CL.+=, CL.E_Var("nbytes1"))),
								    CL.mkCallExp(CL.E_Var("std::memcpy"),
										 [
										   CL.E_Var("cp"),
										   CL.E_Select((CL.E_Select(CL.E_Var(thisName), "worldPos")), "_data"),
										   CL.E_Var("nbytes1")
										]),
								    CL.S_Exp(CL.E_AssignOp(CL.E_Var("cp"), CL.+=, CL.E_Var("nbytes1")))]),
							 CL.S_Block([
								    
								    CL.S_Decl(["const"], arrayTy, "dumb", SOME(CL.I_Array([]))),
								    CL.S_Decl([], CL.size_t, "nbytes3", SOME(CL.I_Exp(CL.E_Sizeof(arrayTy)))),
								    CL.mkCallExp(CL.E_Var("std::memcpy"),
										 [
										   CL.E_Var("cp"),
										   CL.mkAddrOf(CL.E_Var("dumb")),
										   CL.E_Var("nbytes3")
										]),
								    CL.S_Exp(CL.E_AssignOp(CL.E_Var("cp"), CL.+=, CL.E_Var("nbytes3")))
								   ])

							),
						
							CL.S_Return(SOME(CL.E_Var("cp")))
					    ]))
		  val fields = [meshParam, intParam, refPosParam, worldPosParam, worldPosCompute, valid, faceParam, genBuild, defaultBuild,defaultBuild']
		  val class = CL.D_ClassDef{name=name, args=NONE, from=NONE,public = fields , protected = [], private = []}
		  
		 in
		  (class, [printer, allBuild, allBuild', invalidBuild, invalidBuild', refBuild, refBuild', copy_to, makeFem])
		 end
		)
		| _ => raise Fail "Illegal femtype"
		
		
		  
		
	      (* end case*))
	       end
          fun genDecl (ty, (tyDcls, fnDefs)) = (case ty
                 of Ty.VecTy(w, pw) => let
                      val cTyName = RN.vecTyName w
                      val cTy = CL.T_Named cTyName
                      val typedefDcl = CL.D_Verbatim[concat[
                              "typedef ", realTyName, " ", cTyName,
                              " __attribute__ ((vector_size (",
                              Int.toString(realTySz * pw), ")));"
                            ]]
                      in
                        (typedefDcl :: tyDcls, fnDefs)
                      end
                  | Ty.TensorRefTy shape => let
                      val name = RN.tensorRefStruct shape
                      val baseCls = concat[
                              "diderot::tensor_ref<", realTyName, ",",
                              Int.toString(List.foldl Int.* 1 shape), ">"
                            ]
                      fun mkConstr' (paramTy, paramId, arg) = mkConstr (
                            name,
                            [CL.PARAM([], paramTy, paramId)],
                            [CL.mkApply(baseCls, [arg])])
                    (* constructor from float/double pointer *)
                      val (constrProto1, constrDef1) = mkConstr' (
                            CL.constPtrTy realTy, "src", CL.mkVar "src")
                    (* constructor from tensor struct *)
                      val (constrProto2, constrDef2) = mkConstr' (
                            CL.T_Named(concat["struct ", RN.tensorStruct shape, " const &"]),
                            "ten", CL.mkSelect(CL.mkVar "ten", "_data"))
                    (* copy constructor *)
                      val (constrProto3, constrDef3) = mkConstr' (
                            CL.T_Named(name ^ " const &"), "ten",
                            CL.mkSelect(CL.mkVar "ten", "_data"))
                      val thisData = CL.mkIndirect(CL.mkVar "this", "_data")
                    (* last vector as tensor_ref *)
                      val lastDcl = (case shape
                            of [] => raise Fail "unexpected TensorRef[]"
                             | [_] => []
                             | _::dd => let
                                val d = List.last dd
                                in [
                                  CL.D_Func([], RN.tensorRefTy[d], [], "last",
                                    [CL.PARAM([], CL.uint32, "i")],
                                    CL.mkReturn(
                                      SOME(CL.mkAddrOf(CL.mkSubscript(thisData, CL.mkVar "i")))))
                                ] end
                            (* end case *))
                      val members = constrProto1 :: constrProto2 :: constrProto3 :: lastDcl
                      val structDcl = CL.D_ClassDef{
                              name = name,
                              args = NONE,
                              from = SOME("public " ^ baseCls),
                              public = members,
                              protected = [],
                              private = []
                            }
                      in
                        (structDcl :: tyDcls, constrDef1 :: constrDef2 :: constrDef3 :: fnDefs)
                      end
                  | Ty.TensorTy shape => let
                      val len = List.foldl Int.* 1 shape
                      val name = RN.tensorStruct shape
                      val baseCls = concat[
                              "diderot::tensor<", realTyName, ",",
                              Int.toString(List.foldl Int.* 1 shape), ">"
                            ]
                      fun mkConstr (paramTy, paramId, arg) = CL.D_Constr (
                            [], [], name,
                            [CL.PARAM([], paramTy, paramId)],
                            SOME([CL.mkApply(baseCls, [arg])], CL.mkBlock[]))
                    (* default constructor *)
                      val constrDcl1 = CL.D_Constr (
                            [], [], name, [], SOME([CL.mkApply(baseCls, [])], CL.mkBlock[]))
                    (* constructor from initializer list *)
                      val constrDcl2 = mkConstr (
                            CL.T_Template("std::initializer_list", [realTy]), "const & il",
                            CL.mkVar "il")
                    (* constructor from float/double pointer *)
                      val constrDcl3 = mkConstr (
                            CL.constPtrTy realTy, "src", CL.mkVar "src")
                    (* copy constructor *)
                      val constrDcl4 = mkConstr (
                            CL.T_Named(name ^ " const &"), "ten",
                            CL.mkSelect(CL.mkVar "ten", "_data"))
                    (* destructor *)
                      val destrDcl = CL.D_Destr([], [], name, SOME(CL.mkBlock[]))
                      val thisData = CL.mkIndirect(CL.mkVar "this", "_data")
                      val returnThis = CL.mkReturn(SOME(CL.mkUnOp(CL.%*, CL.mkVar "this")))
                    (* assignment from Tensor *)
                      val (assignProto1, assignDef1) = mkMemberFn(name,
                              CL.T_Named(name ^ " &"), "operator=",
                              [CL.PARAM([], CL.T_Named name, "const & src")],
                              CL.mkBlock[
                                  CL.mkCall("this->copy", [CL.mkSelect(CL.mkVar "src", "_data")]),
                                  returnThis
                                ])
                    (* assignment from TensorRef *)
                      val (assignProto2, assignDef2) = mkMemberFn(name,
                              CL.T_Named(name ^ " &"), "operator=",
                              [CL.PARAM([], CL.T_Named(RN.tensorRefStruct shape), "const & src")],
                              CL.mkBlock[
                                  CL.mkCall("this->copy", [CL.mkSelect(CL.mkVar "src", "_data")]),
                                  returnThis
                                ])
                    (* assignment from initializer list *)
                      val (assignProto3, assignDef3) = mkMemberFn(name,
                              CL.T_Named(name ^ " &"), "operator=",
                              [CL.PARAM([], CL.T_Template("std::initializer_list", [realTy]), "const & il")],
                              CL.mkBlock[
                                  CL.mkCall("this->copy", [CL.mkVar "il"]),
                                  returnThis
                                ])
                    (* assignment from array *)
                      val (assignProto4, assignDef4) = mkMemberFn(name,
                              CL.T_Named(name ^ " &"), "operator=",
                              [CL.PARAM([], CL.constPtrTy realTy, "src")],
                              CL.mkBlock[
                                  CL.mkCall("this->copy", [CL.mkVar "src"]),
                                  returnThis
                                ])
                    (* last vector as tensor_ref *)
                      val lastDcl = (case shape
                            of [] => raise Fail "unexpected TensorTy[]"
                             | [_] => []
                             | _::dd => let
                                val d = List.last dd
                                in [
                                  CL.D_Func([], RN.tensorRefTy[d], [], "last",
                                      [CL.PARAM([], CL.uint32, "i")],
                                      CL.mkReturn(
                                        SOME(CL.mkAddrOf(CL.mkSubscript(thisData, CL.mkVar "i")))))
                                ] end
                            (* end case *))
                      val structDcl = CL.D_ClassDef{
                              name = name,
                              args = NONE,
                              from = SOME("public " ^ baseCls),
                              public =
                                  constrDcl1 :: constrDcl2 :: constrDcl3 :: constrDcl4 ::
                                  destrDcl ::
                                  assignProto1 :: assignProto2 :: assignProto3 :: assignProto4 ::
                                  lastDcl,
                              protected = [],
                              private = []
                            }
                      val fnDefs = assignDef1 :: assignDef2 :: assignDef3 :: assignDef4 :: fnDefs
                      in
                        (structDcl :: tyDcls, fnDefs)
                      end
                  | Ty.TupleTy tys =>
		    let
		     val name = CodeGenUtil.tupleName(TreeTypes.TupleTy(tys))
		     val members = List.tabulate(List.length tys, fn x =>
								     let
								      val nty = trType TypeToCxx.NSDiderot env (List.nth(tys, x))
								      val nname = "t_"^(Int.toString x)
								     in
								      CL.D_Var([], nty, [], nname, NONE)
								     end)
		     val osty = CL.T_Named("std::ostream &")
		     (* Printer *)
		     val acceses = [CL.E_Str(")")]@(List.concat(List.rev(List.drop(List.rev(List.tabulate (List.length(tys) + 1, fn x => [ CL.E_Str(","), CL.E_Select(CL.E_Var("obj"), "t_"^(Int.toString (List.length(tys) - x - 1)))])), 1))))@[CL.E_Str("(")]

		     val print = List.foldr (fn (x,y) => CL.E_BinOp(y, CL.#<<, x)) (CL.E_Var("os")) acceses
		     val printer = CL.mkFuncDcl(osty, "operator<<",
						[CL.PARAM([], osty, "os"),
						 CL.PARAM(["const"], CL.T_Named(name), "obj")], 
						CL.S_Block([
							   CL.S_Return(
							    SOME(print))]))

		     (*dumb constructor*)
		     val params = List.tabulate (List.length tys,
						 fn idx =>
						    CL.PARAM([],
							     trType TypeToCxx.NSDiderot env (List.nth(tys, idx)), "p_" ^ (Int.toString idx)
						))
		     val inits = List.tabulate (List.length tys, fn idx =>
								    let
								     val v = CL.E_Var("t_" ^ (Int.toString idx))
								     val vp = CL.E_Var("p_" ^ (Int.toString idx))
								    in CL.E_Apply(v,[vp]) end)
						
		     val constr =  CL.D_Constr(["inline"], [CL.SC_Type(CL.T_Named name)], name, params, SOME(inits, CL.mkBlock([])))

					      (*tuple parse*)
					      
					      
						
		    in
		     (*(CL.D_StructDef(SOME(name), members, NONE) :: tyDcls,  constr :: printer :: fnDefs)*)
		     ( CL.D_ClassDef{name = name, args=NONE, from=NONE,
				   public = members@[constr],
				   protected = [],
				   private = []} :: tyDcls , printer :: fnDefs)
		    end
(* TODO: QUESTION: IS this really needed? I think this is handled elsewhere
                  | Ty.SeqTy(ty, NONE) =>
                  | Ty.SeqTy(ty, SOME n) =>
*)
                  
		  | Ty.FemData(data) =>
		    let
		     val (newTy, newFnDefs) = generateFemType data
		    in
		     (newTy :: tyDcls, List.@ (newFnDefs, fnDefs))
		    end
		  | ty => (tyDcls, fnDefs)
                (* end case *))
          in
            genDecl
          end

  (* generate an instance of the dynamic-sequence trait.  This struct has the following
   * components (see src/lib/include/diderot/dynseq.hxx)
   *
   *    template <>
   *    struct dynseq_traits<T> {
   *        using value_type = T;
   *        using base_type = BASE_TY;
   *        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
   *        static const uint32_t values_per_elem = N;
   *    };
   *    const __details::load_fn_ptr< dynseq_traits< T >::base_type > *dynseq_traits< T >::load_fn_tbl
   *        = LOAD_FN;
   *
   * where
   *    T is the C++ element type of the dynamic sequence (e.g., "tensor_2" for "vec2[]")
   *    BASE_TY is the type of values that comprise an element (e.g., "float" for "vec2[]")
   *    N is the number of values per sequence element (e.g., 1 for "int[]" and 9 for "tensor[3,3][]")
   *    LOAD_FN is the nrrd load function for the base type (e.g., nrrdDLoad if BASE_TY is double)
   *)
    fun genSeqTrait env = let
          val ns = #namespace(Env.target env)
          val realTy = Env.realTy env
          val trType = trType TypeToCxx.NSDiderot env
          fun trait ({argTy, baseTy, elemTy, nValsPerElem}, dcls) = let
              (* the name of the teem function table for the given base type *)
                val loadTbl = (case baseTy
                       of Ty.BoolTy => "nrrdILoad"
                        | Ty.IntTy => "nrrdILoad"
                        | Ty.VecTy(1, 1) => if #double(Env.target env)
					    then "nrrdDLoad"
					    else "nrrdFLoad"
                        | ty => raise Fail("genSeqTrait.loadFn: unexpected type " ^ Ty.toString ty)
                      (* end case *))
                val loadTblTy = CL.constPtrTy(CL.T_Named "__details::load_fn_ptr<base_type>")
                val seqTy = CL.T_Template("dynseq_traits", [argTy])
                val scope = CL.SC_Type seqTy
                in
                  CL.D_Template([], CL.D_ClassDef{
                      name = "dynseq_traits",
                      args = SOME[argTy],
                      from = NONE,
                      public = [
                          CL.D_Typedef("value_type", elemTy),
                          CL.D_Typedef("base_type", trType baseTy),
                          CL.D_Var(
                            ["static"],
                            CL.constPtrTy(CL.T_Named "__details::load_fn_ptr<base_type>"),
                            [], "load_fn_tbl", NONE),
                          CL.D_Var(
                            ["static", "const"], CL.uint32, [], "values_per_elem",
                            SOME(CL.I_Exp(mkInt nValsPerElem)))
                        ],
                      protected = [],
                      private = []
                    }) ::
                  CL.D_Var(
                    ["const"],
                    CL.T_Ptr(
                      CL.T_Template("__details::load_fn_ptr", [CL.T_Member(seqTy, "base_type")])),
                    [scope], "load_fn_tbl",
                    SOME(CL.I_Exp(CL.mkVar loadTbl))) ::
                  dcls
                end
          fun elemTy (Ty.SeqTy(ty, _)) = elemTy ty
            | elemTy (Ty.TensorTy[]) = Ty.realTy
            | elemTy ty = ty
          fun genTrait (ty, dcls) = (case ty
                 of Ty.SeqTy(argTy, NONE) => let
                      val argTy = trType argTy
                    (* for sequences of scalar values, we set nDims to 0 so that it matches the
                     * format of a nrrd, where the dimension is not represented.
                     *)
                      fun scalarSeqTrait ty = trait ({
                                argTy = argTy, baseTy = ty, elemTy = argTy, nValsPerElem = 1
                              },
						     dcls)

		      fun femSeqTrait ty = trait ({argTy = argTy,
						baseTy = ty,
						elemTy = argTy,
						nValsPerElem = 1
					       }, dcls)
					      
                      in
                        case elemTy ty
                         of ty as Ty.TensorTy(shp as _::_) => trait ({
                                  argTy = argTy, baseTy = Ty.realTy,
                                  elemTy = argTy,
                                  nValsPerElem = List.foldl Int.* 1 shp
                                },
                              dcls)
                          | ty as Ty.BoolTy => scalarSeqTrait ty
                          | ty as Ty.IntTy => scalarSeqTrait ty
                          | ty as Ty.VecTy(1, 1) => scalarSeqTrait ty
(* QUESTION: strands map to uint32_t and do not support loading; do we need a trait? *)
                          | ty as Ty.StrandIdTy _ => dcls
			  | ty as Ty.FemData(FemData.MeshCell(_)) => femSeqTrait Ty.IntTy
			  | ty as Ty.FemData(FemData.FuncCell(_)) => femSeqTrait Ty.IntTy
			  | ty as Ty.FemData(FemData.MeshPos(_)) => femSeqTrait Ty.IntTy
                          | ty => raise Fail("unexpected dynamic sequence of " ^ Ty.toString ty)
                        (* end case *)
                      end
                  | _ => dcls
                (* end case *))
          in
            genTrait
          end

    datatype operation = datatype CollectInfo.operation

    val ostreamRef = CL.T_Named "std::ostream&"

    fun output (e, e') = CL.mkBinOp(e, CL.#<<, e')

  (* generate code for the expression "e << s", where "s" is string literal *)
    fun outString (CL.E_BinOp(e, CL.#<<, CL.E_Str s1), s2) =
          output (e, CL.mkStr(s1 ^ String.toCString s2))
      | outString (e, s) = output (e, CL.mkStr(String.toCString s))

  (* generate a printing function for tensors with the given shape *)
    fun genTensorPrinter shape = let
          fun ten i = CL.mkSubscript(CL.mkSelect(CL.mkVar "ten", "_data"), mkInt i)
          fun prefix (true, lhs) = lhs
            | prefix (false, lhs) = outString(lhs, ",")
          fun lp (isFirst, lhs, i, [d]) = let
                fun lp' (_, lhs, i, 0) = (i, outString(lhs, "]"))
                  | lp' (isFirst, lhs, i, n) =
                      lp' (false, output (prefix (isFirst, lhs), ten i), i+1, n-1)
                in
                  lp' (true, outString(lhs, "["), i, d)
                end
            | lp (isFirst, lhs, i, d::dd) = let
                fun lp' (_, lhs, i, 0) = (i, outString(lhs, "]"))
                  | lp' (isFirst, lhs, i, n) = let
                      val (i, lhs) = lp (true, prefix (isFirst, lhs), i, dd)
                      in
                        lp' (false, lhs, i, n-1)
                      end
                in
                  lp' (true, outString(lhs, "["), i, d)
                end
          val params = [
                  CL.PARAM([], ostreamRef, "outs"),
                  CL.PARAM([], RN.tensorRefTy shape, "const & ten")
                ]
          val (_, exp) = lp (true, CL.mkVar "outs", 0, shape)
          in
            CL.D_Func(["static"], ostreamRef, [], "operator<<", params, mkReturn exp)
          end

  (* generate a printing function for fixed-size sequences *)
    fun genSeqPrinter (env, elemTy, size) = let
          val elemTy' = trType TypeToCxx.NSDiderot env elemTy
          val seqTy = sequenceTy (elemTy', size)
          val params = [
                  CL.PARAM([], ostreamRef, "outs"),
                  CL.PARAM([], seqTy, "const & seq")
                ]
          val outsV = CL.mkVar "outs"
          val seqV = CL.mkVar "seq"
          val body = [
                  CL.mkExpStm(outString(outsV, "}")),
                  mkReturn outsV
                ]
        (* get a sequence element *)
          fun getElem ix = CL.mkSubscript(seqV, ix)
          val body = if (size > 0)
                then CL.mkExpStm(output(outsV, getElem(CL.mkInt 0))) ::
                  CL.mkFor(
                    CL.int32, [("i", CL.mkInt 1)],
                    CL.mkBinOp(CL.mkVar "i", CL.#<, CL.mkDispatch(seqV, "size", [])),
                    [CL.mkUnOp(CL.%++, CL.mkVar "i")],
                    CL.mkExpStm(output(outString(outsV, ","), getElem(CL.mkVar "i"))))
                  :: body
                else body
          val body = CL.mkExpStm(outString(outsV, "{")) :: body
          in
            CL.D_Func(["static"], ostreamRef, [], "operator<<", params, CL.mkBlock body)
          end

  (* builds AST for the expression "(x <= lo) ? lo : (hi <= x) ? hi : x;" *)
    fun mkClampExp (lo, hi, x) =
          CL.mkCond(CL.mkBinOp(x, CL.#<=, lo), lo,
            CL.mkCond(CL.mkBinOp(hi, CL.#<=, x), hi,
              x))

    fun expandFrag (env, frag) =
          CL.verbatimDcl [frag] [("REALTY", if #double(Env.target env) then "double" else "float")]

    fun imageTy (realTy, d) =
          CL.T_Template(RN.qImageTyName d, [realTy, CL.T_Named "TY", CL.T_Named "VOXSZ"])

    fun doOp env (rator, dcls) = let
     val realTy = Env.realTy env
     val intTy = Env.intTy env
          fun mkVec (w, pw, f) = CL.mkVec(
                RN.vecTy w,
                List.tabulate(pw, fn i => if i < w then f i else CL.mkFlt(zero, realTy)))
          fun mkVMap (ty, name, f, w, pw) = let
                fun f' i = CL.mkApply(f, [CL.mkSubscript(CL.mkVar "v", mkInt i)])
                in
                  mkFunc(ty, name, [CL.PARAM([], ty, "v")], mkReturn (mkVec (w, pw, f')))
                end
          val dcl = (case rator
                 of Print(Ty.TensorRefTy shape) => genTensorPrinter shape
                  | Print(Ty.TupleTy tys) => CL.D_Verbatim[] (* no printer needed*)
                  | Print(Ty.SeqTy(ty, NONE)) => CL.D_Verbatim[] (* no printer needed *)
                  | Print(Ty.SeqTy(ty, SOME n)) => genSeqPrinter (env, ty, n)
                  | Print ty => CL.D_Verbatim[] (* no printer needed *)
                  | RClamp => let
                      val params = [
                              CL.PARAM([], realTy, "lo"),
                              CL.PARAM([], realTy, "hi"),
                              CL.PARAM([], realTy, "x")
                            ]
                      in
                        mkFunc(realTy, "clamp", params,
                          mkReturn(mkClampExp (CL.mkVar "lo", CL.mkVar "hi", CL.mkVar "x")))
                      end
                  | RLerp => let
                      val params = [
                              CL.PARAM([], realTy, "a"),
                              CL.PARAM([], realTy, "b"),
                              CL.PARAM([], realTy, "t")
                            ]
                      in
                        mkFunc(realTy, "lerp", params,
                          mkReturn (
                            CL.mkBinOp(
                              CL.mkVar "a",
                              CL.#+,
                              CL.mkBinOp(
                                CL.mkVar "t",
                                CL.#*,
                              CL.mkBinOp(CL.mkVar "b", CL.#-, CL.mkVar "a")))))
                      end
                  | VScale(w, pw) => let
                      val cTy = RN.vecTy w
                      in
                        mkFunc(cTy, RN.vscale w,
                          [CL.PARAM([], realTy, "s"), CL.PARAM([], cTy, "v")],
                          mkReturn(
                            CL.mkBinOp(mkVec(w, pw, fn _ => CL.mkVar "s"), CL.#*, CL.mkVar "v")))
                      end
                  | VSum(w, pw) => let
                      val name = RN.vsum w
                      val params = [CL.PARAM([], RN.vecTy w, "v")]
                      fun mkSum 0 = CL.mkSubscript(CL.mkVar "v", mkInt 0)
                        | mkSum i = CL.mkBinOp(mkSum(i-1), CL.#+, CL.mkSubscript(CL.mkVar "v", mkInt i))
                      in
                        mkFunc(realTy, name, params, mkReturn(mkSum(w-1)))
                      end
                  | VDot(w, pw) => let
                      val name = RN.vdot w
                      val vTy = RN.vecTy w
                      val params = [CL.PARAM([], vTy, "u"), CL.PARAM([], vTy, "v")]
                      fun mkSum 0 = CL.mkSubscript(CL.mkVar "w", mkInt 0)
                        | mkSum i = CL.mkBinOp(mkSum(i-1), CL.#+, CL.mkSubscript(CL.mkVar "w", mkInt i))
                      in
                        mkFunc(realTy, name, params,
                          CL.mkBlock[
                              CL.mkDeclInit(vTy, "w", CL.mkBinOp(CL.mkVar "u", CL.#*, CL.mkVar "v")),
                              mkReturn(mkSum(w-1))
                            ])
                      end
                  | VCeiling(w, pw) => mkVMap (RN.vecTy w, RN.vceiling w, "diderot::ceiling", w, pw)
                  | VFloor(w, pw) => mkVMap (RN.vecTy w, RN.vfloor w, "diderot::floor", w, pw)
                  | VRound(w, pw) => mkVMap (RN.vecTy w, RN.vround w, "diderot::round", w, pw)
                  | VTrunc(w, pw) => mkVMap (RN.vecTy w, RN.vtrunc w, "diderot::trunc", w, pw)
                  | VToInt(layout as {wid, pieces, ...}) => let
                      val intTy = Env.intTy env
                      val seqTy = sequenceTy (CL.intTy, wid)
                      fun mkItem i =
                            CL.I_Exp(CL.mkCons(intTy, [CL.mkSubscript(CL.mkVar "src", mkInt i)]))
                      val vParamTys = Ty.piecesOf layout
                      val vParams = List.mapi
                            (fn (i, Ty.VecTy(w, _)) => CL.PARAM([], RN.vecTy w, "v"^Int.toString i))
                              vParamTys
                      val initItems = let
                            fun doPiece (i, Ty.VecTy(w, _)) = let
                                  val src = CL.mkVar("v"^Int.toString i)
                                  fun mkItem j =
                                        CL.I_Exp(CL.mkCons(intTy, [CL.mkSubscript(src, mkInt j)]))
                                  in
                                    List.tabulate (w, mkItem)
                                  end
                            in
                              List.concat(List.mapi doPiece vParamTys)
                            end
                      in
                        mkFunc(seqTy, RN.vtoi wid,
                          vParams,
                          CL.mkBlock[
                              CL.mkDecl(seqTy, "res", SOME(CL.I_Exps initItems)),
                              mkReturn(CL.mkVar "res")
                            ])
                      end
                  | VLoad(w, pw) => let
                      val name = RN.vload w
                      val cTy = RN.vecTy w
                      fun arg i = CL.mkSubscript(CL.mkVar "vp", mkInt i)
                      in
                        mkFunc(cTy, name,
                          [CL.PARAM(["const"], CL.T_Ptr realTy, "vp")],
                          mkReturn(mkVec (w, pw, arg)))
                      end
                  | VCons(w, pw) => let
                      val name = RN.vcons w
                      val cTy = RN.vecTy w
                      val params = List.tabulate(w, fn i => CL.PARAM([], realTy, "r"^Int.toString i))
                      fun arg i = CL.mkVar("r"^Int.toString i)
                      in
                        mkFunc(cTy, name, params, mkReturn(mkVec (w, pw, arg)))
                      end
                  | VPack layout => let
                      val name = RN.vpack (#wid layout)
                      val vParamTys = Ty.piecesOf layout
                      val vParams = List.mapi
                            (fn (i, Ty.VecTy(w, _)) => CL.PARAM([], RN.vecTy w, "v"^Int.toString i))
                              vParamTys
                      val dstTy = RN.tensorTy[#wid layout]
                      fun mkAssign (i, v, j) =
                            CL.mkAssign(
                              CL.mkSubscript(CL.mkSelect(CL.mkVar "dst", "_data"), mkInt i),
                              CL.mkSubscript(v, mkInt j))
                      fun mkAssignsForPiece (dstStart, pieceIdx, wid, stms) = let
                            val piece = CL.mkVar("v"^Int.toString pieceIdx)
                            fun mk (j, stms) = if (j < wid)
                                  then mk (j+1, mkAssign (dstStart+j, piece, j) :: stms)
                                  else stms
                            in
                              mk (0, stms)
                            end
                      fun mkAssigns (_, [], _, stms) = CL.mkBlock(List.rev stms)
                        | mkAssigns (i, Ty.VecTy(w, _)::tys, offset, stms) =
                            mkAssigns (i+1, tys, offset+w, mkAssignsForPiece(offset, i, w, stms))
                      in
                        mkFunc(CL.voidTy, name,
                          CL.PARAM([], dstTy, "&dst") :: vParams,
                          mkAssigns (0, vParamTys, 0, []))
                      end
                  | TensorCopy shp => CL.D_Verbatim[]
(*
                  | TensorCopy shp => let
                      val name = RN.tensorCopy shp
                      val dim = List.foldl Int.* 1 shp
                      val dstTy = CL.T_Array(realTy, SOME dim)
                      in
                        mkFunc(CL.voidTy, name,
                          [CL.PARAM([], dstTy, "dst"), CL.PARAM([], CL.constPtrTy realTy, "src")],
                          CL.mkCall("std::memcpy", [
                              CL.mkVar "dst", CL.mkVar "src", CL.mkSizeof dstTy
                            ]))
                      end
*)
                  | Transform d => let
                      val e = CL.mkDispatch(CL.mkVar "img", "world2image", [])
                      val (resTy, e) = if (d = 1)
                            then (realTy, e)
                            else let val ty = RN.tensorRefTy[d, d]
                              in (ty, CL.mkCons(ty, [e])) end
                      in
                        CL.D_Template([CL.TypeParam "TY", CL.ConstParam(CL.intTy, "VOXSZ")],
                          mkFunc(resTy, "world2image",
                            [CL.PARAM([], imageTy (realTy, d), "const & img")],
                            CL.mkReturn(SOME e)))
                      end
                  | Translate d => let
                      val e = CL.mkDispatch(CL.mkVar "img", "translate", [])
                      val (resTy, e) = if (d = 1)
                            then (realTy, e)
                            else let val ty = RN.tensorRefTy[d]
                              in (ty, CL.mkCons(ty, [e])) end
                      in
                        CL.D_Template([CL.TypeParam "TY", CL.ConstParam(CL.intTy, "VOXSZ")],
                          mkFunc(resTy, "translate",
                            [CL.PARAM([], imageTy (realTy, d), "const & img")],
                            CL.mkReturn(SOME e)))
                      end
                  | Inside(layout, s) => let
                      val dim = #wid layout
                      val vTys = List.map
                            (fn ty => TypeToCxx.trType (env, ty))
                              (TreeTypes.piecesOf layout)
                      val xs = List.mapi (fn (i, ty) => "x"^Int.toString i) vTys
                      val vParams =
                            ListPair.map (fn (ty, x) => CL.PARAM([], ty, x)) (vTys, xs)
                    (* make the tests `(x < img.size(i)-s)` and `((s-1) < x)` *)
                      fun mkTests (x, i) = [
                              CL.mkBinOp(x, CL.#<,
                                CL.mkBinOp(
                                  CL.mkDispatch(CL.mkVar "img", "size", [mkInt i]),
                                  CL.#-, mkInt s)),
                              CL.mkBinOp(mkInt(s-1), CL.#<, x)
                            ]
                    (* build the test expression from the pieces *)
                      fun mkExps (i, w, v::vr, pw::pwr, tests) =
                            if (i < dim)
                              then if (w < pw)
                                then let
                                  val x = if (pw = 1)
                                        then CL.mkVar v
                                        else CL.mkSubscript(CL.mkVar v, mkInt w)
                                  in
                                    mkExps (i+1, w+1, v::vr, pw::pwr,
                                      mkTests(x, w) @ tests)
                                  end
                                else mkExps (i, pw, vr, pwr, tests)
                              else List.rev tests
                        | mkExps _ = raise Fail "inconsistent"
                      val (t1::tr) = mkExps (0, 0, xs, #pieces layout, [])
                      val exp = List.foldr
                            (fn (e1, e2) => CL.mkBinOp(e2, CL.#&&, e1))
                              t1 tr
                      in
                        CL.D_Template([CL.TypeParam "TY", CL.ConstParam(CL.intTy, "VOXSZ")],
                          mkFunc(CL.boolTy, RN.inside(dim, s),
                            vParams @ [CL.PARAM([], imageTy (realTy, dim), "img")],
                            mkReturn exp))
                      end
                  | EigenVals2x2 => expandFrag (env, CxxFragments.eigenvals2x2)
                  | EigenVals3x3 => expandFrag (env, CxxFragments.eigenvals3x3)
                  | EigenVecs2x2 => expandFrag (env, CxxFragments.eigenvecs2x2)
                  | EigenVecs3x3 => expandFrag (env, CxxFragments.eigenvecs3x3)
                  | SphereQuery(d, s) => let
                      val seqTy = CL.T_Template("diderot::dynseq", [CL.T_Named "strand_array::sid_t"])
                      val (posParam, posExp) = if d > 1
                            then (
                                CL.PARAM([], TypeToCxx.trType(env, Ty.TensorRefTy[d]), "_pos"),
                                CL.mkDispatch(CL.mkVar "_pos", "base", [])
                              )
                            else (CL.PARAM([], realTy, "_pos"), CL.mkAddrOf(CL.mkVar "_pos"))
                      val wrldV = CL.mkVar RN.worldVar
                      val selfParam = CL.PARAM([], CL.constPtrTy(RN.strandTy s), RN.selfVar)
                      val radiusParam = CL.PARAM([], realTy, "radius")
                      in
                        mkFunc(
                          seqTy, "sphere_query",
                          [RN.worldParam, selfParam, posParam, radiusParam],
                          CL.mkReturn(SOME(
                            CL.mkIndirectDispatch(
                              CL.mkIndirect(wrldV, "_tree"),
                              "sphere_query",
                              [CL.mkVar RN.selfVar, posExp, CL.mkVar "radius"]))))
                      end
                  | RIfWrap => let
                      val t1 = "float IfWrap (bool  e1, float e3, float e4)"
                      val t2 = "\n\t{"
                      val t3 = " \n\t\t if(e1){return e3;}"
                      val t4 = "\n\t\t    else{return e4;}"
                      val t5 = "\n\t}"
                      val es = [t1,t2,t3,t4,t5]
                  in CL.D_Verbatim es end
		  | MeshGeometryQueryInsert(file, ty) =>
		    let
		     val file = TextIO.openIn file
		     fun loop ins = 
			 (case TextIO.inputLine ins of 
			     SOME line => line :: loop ins 

			   | NONE      => []
			 (*end case*))
		     val lines = loop file
		     val (rawReal, rawInt) = (Env.rawRealTy env, Env.rawIntTy env)
		     val (rawRealStr, rawIntTy) = (PrintAsCxx.rawTyName rawReal,
						   PrintAsCxx.rawTyName rawInt
						  )
		     val subs = [("meshTy", ty), ("realTy", rawRealStr), ("intTy", rawIntTy)]
		     val stm = CL.S_Verbatim(List.map (StringSubst.expand subs) lines)
		     val seqTy = CL.T_Template("diderot::dynseq", [intTy])
		     val fin = CL.D_Func([], seqTy, [], "mesh_geom_"^ty,
					 [CL.PARAM([], CL.voidPtr, "index"),
					  CL.PARAM([], CL.T_Ptr(CL.T_Named(ty)), "mesh"),
					  CL.PARAM(["const"], CL.T_Ptr(realTy), "data")],
					 stm)
		    in
		     fin
		    end
                (* end case *))
          in
            dcl :: dcls
          end

    val firstTy = CL.D_Comment["***** Begin synthesized types *****"]
    val lastTy = CL.D_Comment["***** End synthesized types *****"]
    val noDclsTy = CL.D_Comment["***** No synthesized types *****"]

    val firstOp = CL.D_Comment["***** Begin synthesized operations *****"]
    val lastOp = CL.D_Comment["***** End synthesized operations *****"]
    val noDclsOp = CL.D_Comment["***** No synthesized operations *****"]

    fun gen (env, info) = let
          val spec = Env.target env
          val genTrait = genSeqTrait env
          val genTyDecl = genTyDecl env
          val opDcls = List.foldl (doOp env) [] (CollectInfo.listOps info)
          val tys = CollectInfo.listTypes info
          val (tyDcls, fnDefs) = List.foldr genTyDecl ([], []) tys
          val dcls = tyDcls @ fnDefs
          val traitDcls = List.foldl genTrait [] tys
          val preDcls = if List.null dcls andalso List.null traitDcls
                then [noDclsTy]
                else let
                  val res = [lastTy]
                  val res = if List.null traitDcls
                        then res
                        else CL.D_Namespace("diderot", traitDcls) :: res
                  val res = if List.null dcls
                        then res
                        else CL.D_Namespace(#namespace(Env.target env), dcls) :: res
                  in
                    firstTy :: res
                  end
          val postDcls = if List.null opDcls
                then [noDclsOp]
                else firstOp :: opDcls @ [lastOp]
          in
            {preWorld = preDcls, postWorld = postDcls}
          end

  end
