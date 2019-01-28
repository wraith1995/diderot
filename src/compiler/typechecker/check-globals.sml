(* check-globals.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2017 The University of Chicago
 * All rights reserved.
 *)

structure CheckGlobals : sig

  (* type check the global declarations of a program.  We partition the result
   * into constants (which can only depend on other constants), inputs (which
   * can depend on contants), and the other globals (functions and variables).
   *)
    val check : CmdLineConstants.t * Env.t * Env.context * ParseTree.global_dcl list -> {
            const_dcls : AST.var_dcl list,
            input_dcls : (AST.var_dcl * string option) list,
            other_dcls : AST.global_dcl list,
            env : Env.t
          }

  end = struct

    structure PT = ParseTree
    structure Ty = Types
    structure TU = TypeUtil
    structure E = Env
    structure L = Literal
    structure SS = Substring
    structure CO = CheckOverload
    structure FT = FemData
    structure N = BasisNames
    structure CF = CheckTypeFile

		     
    val err = TypeError.error

    datatype token = datatype TypeError.token

  (* tagged union of the different kinds of global declarations *)
    datatype dcl_kind
      = CONST of AST.var_dcl
      | INPUT of (AST.var_dcl * string option)
      | OTHER of AST.global_dcl
      | OTHERS of AST.global_dcl list
      | ERROR
      | IGNORE

  (* check the rhs initialization of a 'const' or 'input' declaration.  Return
   * (v, e), where v is the constant value and e is the AST version of the rhs.
   *)
    fun chkRHS (env, cxt, isInput, x', e) = let
          val (e', ty') = CheckExpr.check (env, cxt, e)
          val v = (case CheckConst.eval (cxt, isInput, e')
                 of SOME v => v
                  | NONE => ConstExpr.Expr e' (* error has already been reported *)
                (* end case *))
          val e' = ConstExpr.valueToExpr v
          val lhsTy = Var.monoTypeOf x'
          in
            case Util.coerceType (lhsTy, (e', ty'))
             of SOME e' => (v, e')
              | NONE => (
                  err (cxt, [
                      S "definition of ", V x', S " has wrong type\n",
                      S "  expected: ", TY lhsTy, S "\n",
                      S "  found:    ", TY ty'
                    ]);
                  (v, e'))
            (* end case *)
          end

    val bogusCExp = ConstExpr.Int 0
    val bogusExp = AST.E_Lit(L.Int 0)

    fun chkDcl (defs, env, cxt, dcl) = (case dcl
           of PT.GD_Mark{span, tree} => chkDcl (defs, env, (#1 cxt, span), tree)
            | PT.GD_Const(ty, {span, tree=x}, optDefn) => let
                val ty = CheckType.check (env, cxt, ty)
                val x' = Var.new(x, span, Var.ConstVar, ty)
                val override = (case CmdLineConstants.getConst defs x
                       of SOME def => (case ConstExpr.fromString (ty, def)
                             of NONE => (
                                  Error.error (#1 cxt, [
                                      "invalid command-line definition of constant '",
                                      Atom.toString x, "'"
                                    ]);
                                  NONE)
                              | someV => someV
                            (* end case *))
                        | NONE => NONE)
                val (v, e') = (case (optDefn, override)
                       of (SOME e, NONE) => chkRHS (env, cxt, false, x', e)
                        | (SOME e, SOME def) => (
                            ignore (chkRHS (env, cxt, false, x', e));
                            (def, ConstExpr.valueToExpr def))
                        | (NONE, SOME def) => (def, ConstExpr.valueToExpr def)
                        | (NONE, NONE) => (
                            err (cxt, [
                                S "must supply r.h.s. or command-line definition for constant '",
                                V x', S "'"
                              ]);
                            (bogusCExp, bogusExp))
                      (* end case *))
                in
                (* check that const variables have valid types *)
                  if not(TU.isValueType ty)
                    then err (cxt, [S "const variable ", V x', S " has invalid type ", TY ty])
                    else ();
                  E.checkForRedef (env, cxt, x);
                  E.recordProp (env, Properties.HasConsts);
                  ConstExpr.define(x', v);
                  (CONST(x', SOME e'), E.insertGlobal(env, cxt, x, x'))
                end
            | PT.GD_Input(ty, {span, tree=x}, optDesc, optDefn) => let
                val ty = CheckType.check (env, cxt, ty)
                val x' = Var.new(x, span, Var.InputVar, ty)
		val _ = print(Var.nameOf x')
                val rhs = (case optDefn
			    of NONE => (case ty
					 of Ty.T_Named(_,Ty.T_Fem(data, _)) =>
					    (case data
					      of FemData.Mesh(_) => (print("agg!\n");SOME(AST.E_LoadFem(data, NONE, NONE)))
					       | _ => (print("agg\n");NONE))
					  | _ => (print("smesh\n");NONE))
                             | SOME e => (print("aha!"); SOME(#2 (chkRHS (env, cxt, true, x', e))))
                      (* end case *))
                in
                (* check that input variables have valid types *)
                  if not(TU.isValueType ty orelse TU.isImageType ty orelse TU.isInputFemType ty)
                    then err (cxt, [S "input variable ", V x', S " has invalid type ", TY ty])
                    else ();
                  E.checkForRedef (env, cxt, x);
                  E.recordProp (env, Properties.HasGlobals);
                  E.recordProp (env, Properties.HasInputs);
                  (INPUT((x', rhs), optDesc), E.insertGlobal(env, cxt, x, x'))
                end
            | PT.GD_Var varDcl => let
                val (x, x', optDefn) = CheckStmt.checkVarDecl (env, cxt, Var.GlobalVar, varDcl)
                in
                  E.checkForRedef (env, cxt, x);
                  E.recordProp (env, Properties.HasGlobals);
                  (OTHER(AST.D_Var(x', optDefn)), E.insertGlobal(env, cxt, x, x'))
                end
            | PT.GD_Func(ty, {span, tree=f}, params, body) => let
                val ty' = CheckType.check(env, cxt, ty)
                val env' = E.functionScope (env, ty', f)
                val (params', env') = CheckParams.check (env', cxt, Var.FunParam, params)
                val body' = (case body
                       of PT.FB_Expr e => let
                            val eTy = CheckExpr.check (env', cxt, e)
                            in
                              case Util.coerceType(ty', eTy)
                               of SOME e' => AST.S_Return e'
                                | NONE => (
                                    err (cxt, [
                                        S "type of function body does not match return type\n",
                                        S "  expected: ", TY ty', S "\n",
                                        S "  found:    ", TY(#2 eTy)
                                      ]);
                                    AST.S_Block[])
                              (* end case *)
                            end
                        | PT.FB_Stmt s => CheckStmt.check(env', cxt, s)
                      (* end case *))
                val fnTy = Ty.T_Fun(List.map Var.monoTypeOf params', ty')
                val f' = Var.new (f, span, AST.FunVar, fnTy)
                in
(* QUESTION: should we also check for redefinition of basis functions? *)
                  E.checkForRedef (env, cxt, f);
                  (OTHER(AST.D_Func(f', params', body')), E.insertFunc(env, cxt, f, f'))
            end
	    | PT.GD_Overloading({span, tree=f}) =>
	      (case E.findFunc(env, f) 
		of E.PrimFun ([]) =>
		   (E.insertOverload(env, cxt, f); (IGNORE, env)) (*Added overload*)
		|  E.UserFun(g) => 
		   (*For now, I'm mandating overload declerations come before any defintions; This could be changed easily.*)
		   (err (cxt, [
			 S "Declared a function ", A(f), S" to be overloaded after the function was defined on ", LN (Error.location(#1 cxt ,Var.locationOf g))]); 
		    (ERROR, env))
		|  E.OverloadUserFun(g::gs)=>
		   (err (cxt, [
			 S "Declared a function ", A(f), S" to be overloaded multiple times, including here and on ", LN (Error.location(#1 cxt ,Var.locationOf g))]);
		    (ERROR, env))
		| _ =>
		  (err (cxt, [
			S "The symbol", A(f), S" is already an overloaded symbol in the basis yet is declared to be overloaded!"]);
		   (ERROR, env))
	      (* end case *) )
	    | PT.GD_Overload(ty, ({span, tree=ff}, isOp), params, body) =>
	      let
	       (* Setup this as if we have a function at first. *)
	       val ty' = CheckType.check(env, cxt, ty)
               val env' = E.functionScope (env, ty', ff)
               val (params', env') = CheckParams.check (env', cxt, Var.FunParam, params)
               val body' = (case body
			     of PT.FB_Expr e => let
                              val eTy = CheckExpr.check (env', cxt, e)
                             in
                              case Util.coerceType(ty', eTy)
                               of SOME e' => AST.S_Return e'
                                | NONE => (
                                 err (cxt, [
                                      S "type of overload function body does not match return type\n",
                                      S "  expected: ", TY ty', S "\n",
                                      S "  found:    ", TY(#2 eTy)
                                     ]);
                                 AST.S_Block[])
					    (* end case *)
                             end
                              | PT.FB_Stmt s => CheckStmt.check(env', cxt, s)
			   (* end case *))
			     
	       val paramTys = List.map Var.monoTypeOf params'
	       (*We interrupt this to catch the unary operators*)
	       val isBinOp = List.length paramTys = 2
	       val isUnaryOp = List.length paramTys = 1
	       (*QUESTION: any other unary operators to watch out for?*)
	       val f = if not isOp then ff
			   else if isUnaryOp then if Atom.same(ff, N.op_sub)
						  then N.op_neg
						  else ff
			   else ff

	       val fnTy = Ty.T_Fun(paramTys, ty')
               val f' = Var.new (ff, span, AST.FunVar, fnTy)
	       (* Now, we setup utilities to handle checks specific to an overload *)
	       (* If we can't find an overload, we will call this error.*)
	       fun noOverLoadErr r = (err (cxt, [
					 S "Declared an overloaded function ", A(r),
					 S" that has not beeen declared to be overloadedable"]);
				      (ERROR, env))
	       (*QUESTION: should use the overload mechanic to find all duplicates?*)
	       fun findConflictingOverload funcs =
		   let
		    val possibleOverload = CO.chkOverload(cxt, f, paramTys, [], funcs)
							 (*Note: ctx, f', and [] don't matter here unless BV.pow is involved, which is impossible right now*)
		    
		   in
		    case possibleOverload
		     of CO.Overload(AST.E_Prim(overloadName, _, _, _), ty) => [overloadName]
		      | CO.Overload(AST.E_Apply((overloadName, _), _ , _), ty) => [overloadName]
		      | CO.Overload(_) => raise Fail "Impossible: overload table holds non-prim and non-apply"
		      | _ => []
		   end
	       (* This produces a value that let's use identify when we need to worry about the special cases of :,\cdot,\otimes*)
	       val specialProductCheck = (case CO.checkSpecialProduct(paramTys, f)
					   of SOME(NONE, name) => true (* an operator that can be overloaded but is not in the basis  *)
					   |  SOME(s, name) => (err (cxt, [
							       S "Declared an overload for a builtin operator", A(f),
							       S "that is defined via the type ", TY(fnTy),
							       S "that is compatible with definition of ", A(name), S "."
							      ]); true) 
					    | NONE => false (* if it can be overloaded, it is in the basis *)
					 (*end case*))
	       fun overloadExtraToErrors (x::xs) = List.@ ([S "The symbol ", V(x), S "that has type", TY (Var.monoTypeOf x), S "and is defined on", LN (Error.location(#1 cxt ,Var.locationOf x)), S"\n"], (overloadExtraToErrors xs))
		 | overloadExtraToErrors [] = []
	       fun overloadPrimExtraToErrors (x::xs) = List.@ ([S "The symbol ", V(x), S "that has a type defined in the documentation", S "and might be defined on", LN (Error.location(#1 cxt ,Var.locationOf x)), S"\n"], (overloadExtraToErrors xs))
		 | overloadPrimExtraToErrors [] = []

						    
	      in
	       (case E.findFunc(env, f) 
		 of  E.UserFun(g) =>  noOverLoadErr f (* no overload declared*)
		   | E.OverloadUserFun(funcs) =>
		     if isOp
		     then (err (cxt, [
				S "Declared an operator overload ",
				A(f), S" for a user function"]);
			   (ERROR, env)) 
		     else
		      (case findConflictingOverload funcs 
			of [] => (OTHER(AST.D_Func(f', params', body')), E.insertOverloadFunc(env, cxt, f, f' :: funcs))
			 | possible =>
			   let
			    val overloadMsgs = overloadExtraToErrors possible
			   in
			   (err (cxt, List.@ ([
				 S "Declared a function overload ", A(f),
				 S " has a type", TY(fnTy) , S " that conflicts with previous definitions of this function: "
				], overloadMsgs)); 
			    (ERROR, env))
			   end
		      (* end case *))
		   | E.PrimFun (funcs) =>
		     if not isOp
		     then (err (cxt, [
				S "Declared a function overload '",
				A(f), S"' for an operator"]);
			   (ERROR, env))
		     else
		      (case (funcs, specialProductCheck)
			of ([], false) => noOverLoadErr f
			 | (_, _) => (
			  (case findConflictingOverload funcs 
			    of [] =>
				(OTHER(AST.D_Func(f', params', body')), E.insertOverloadPrim(env, cxt, f, f' :: funcs))
			     | others =>
			       let
				val overloadErros = overloadPrimExtraToErrors others
			       in
				
			       (err (cxt, List.@ ([
				     S "Declared an operator overload '", A(f),
				     S " has a type", TY(fnTy) , S " that conflicts with a previous definitions: " 
				    ], overloadErros));
				     (ERROR, env))
				end
			 (* end case *)))
		      (*end case*))
	       (* end case *))
	      end
	    | PT.GD_Type(tyDefBind, {tree=tyName, span=span}, file) =>
	      let
	       (*Note: I'm going to reuse the span here a lot in places where it hopefully shouldn't matter.*)
	       (*The zeroeth step is to determine if we are building an fem or non-fem based type.*)
	       (*FIX ME: reorganize this code*)
	       fun reformatConstants cs =  List.map (fn (x,y,z) => (Atom.atom z,y)) cs
	       val tyDef = (case tyDefBind
			     of PT.T_Mark({span, tree}) => tree
			      | _ => tyDefBind
			   (* end case*))
				
	       val isFemType = (case tyDef
				 of PT.T_Mesh => true
				  | PT.T_Space(_, _) => true
				  | PT.T_Func(_) =>  true 
				  (*TODO: Add Sequences if you want???*)
				  | _ =>   false)
	       (* TODO: Add file parsing; for now, we do dummy info.*)
	       (* This needs to be moved somewhere else... I think we need it in the type utilities*)
	       (* At the least, we need to isolate the key functionality *)
	       fun validateFemType femTyDef fileinfo =
		   let
		    (*TODO: parse the file and make these correct...*)
		    val parsedJson = CF.loadJson(fileinfo, cxt)

		    (*TODO : check against the file in the following.*)

		   in
		    (case (femTyDef, parsedJson)
		      of (PT.T_Mesh, SOME(parsed)) =>
			 let
			  val (tempMesh, consts) = CheckTypeFile.parseMesh(env, cxt, tyName, parsed)
			 in
			   (Option.mapPartial (fn x => SOME(x, NONE)) tempMesh, consts)
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
			    of (SOME(SOME(space''), consts), SOME(meshType')) => (SOME(space'', SOME(mesh)), consts)
			     (*NOTE: these errors could be improved.*)
			     | (NONE, SOME(_)) => ((err (cxt, [
				S "Declared a function space type ",
				A(tyName), S" with an invalid space "])); (NONE, []))
			    | (_, NONE) => ((err (cxt, [
				S "Declared a function space type ",
				A(tyName), S" with an underlying type that is not a mesh or not defined", A(mesh)])); (NONE, []))
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
			    of (SOME(SOME(func), consts),SOME(space')) => (SOME(func, SOME(space)), consts)
			     | (SOME(NONE,_), SOME(_)) => ((err (cxt, [
				S "Declared a femFunction type ",
				A(tyName), S" with an invalid definition in the file or incompatability with the space.", A(space)])); (NONE, []))
			    |  (_, NONE) =>  ((err (cxt, [
				S "Declared a femFunction type ",
				A(tyName), S" with an underlying type that is not a space or not defined", A(space)])); (NONE, []))
			  (* end case *))
			 end
		       | (_, NONE) => (err (cxt, [S "Unable to parse file for declared fem type:", A(tyName) ]); (NONE, []))
		       | _ => raise Fail ("Non fem type passed to validateFemType: " ^ Atom.toString(tyName))
		    (* end case *))
		   end

	       fun makeFemMethods femTyDef file femInfo =
		   (case (femInfo, femTyDef)
		     of ((FT.Func(f), n), PT.T_Func(space)) =>
			let

			 val spaceFunc = Atom.atom "space"
			 val space' = FT.spaceOf(FT.Func(f))
			 val spaceType = let
			  val m = FT.meshOf(space')
			  val meshName = FT.nameOf(m)
			 in Ty.T_Fem(space', SOME(meshName)) end
			 val spaceFuncTy = Ty.T_Fun([Ty.T_Named(tyName,Ty.T_Fem(femInfo))],Ty.T_Named(space, spaceType))
			 val spaceFuncVar = Var.new (spaceFunc, span, AST.FunVar, spaceFuncTy)
			 val param = Var.new (Atom.atom "arg0", span, AST.FunParam, Ty.T_Named(tyName,Ty.T_Fem(femInfo)))
					     
			 val spaceFuncBody = AST.S_Block([AST.S_Return(AST.E_ExtractFem(AST.E_Var(param,span),space'))])
			 val spaceFun = AST.D_Func(spaceFuncVar, [param], spaceFuncBody)
			 
			in
			 ([(spaceFunc, spaceFuncVar)],[spaceFun],[])
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
			 ([(meshFunc, meshFuncVar)], [meshFun],[])
			end
		      | ((FT.Mesh(m), n), PT.T_Mesh) =>
			let
			 val meshArgTy = Ty.T_Named(tyName, Ty.T_Fem(femInfo))
			 val numCell = Atom.atom "numCell"
			 val meshType = Ty.T_Fem(FT.Mesh(m), SOME(tyName))
			 val numCellFuncTy = Ty.T_Fun([meshArgTy],Ty.T_Int)
			 val numCellFuncVar = Var.new (numCell, span, AST.FunVar, numCellFuncTy)
			 val param = Var.new (Atom.atom "arg0", span, AST.FunParam, meshArgTy)
			 val numCellBody = AST.S_Block([AST.S_Return(AST.E_ExtractFemItem(AST.E_Var(param, span), Ty.T_Int, (FemOpt.NumCell, FT.Mesh(m))))])
			 val numCellFun = AST.D_Func(numCellFuncVar, [param], numCellBody)


			 (*get list of cells -> we are going to use a extract*)
			 val cells' = Atom.atom "0cells"
			 val cells = Atom.atom "cells"
			 val cellType = Ty.T_Fem(FT.MeshCell(m), NONE)
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
			 				(*AST.E_ExtractFemItem(AST.E_Var(param', span)*)
			 val cellsBody = AST.S_Block([AST.S_Return(comp)])
			 val cellsFun' =  AST.D_Func(cellFuncVar', [param'], cellsBody)

			 val cellFuncVar = Var.new (cells, span, AST.FunVar, cellTypeFuncTy)
			 val cellsFun = AST.D_Func(cellFuncVar, [param'],
						   AST.S_Block([AST.S_Return(AST.E_ExtractFemItem(AST.E_Var(param',span), cellSeqType, ( FemOpt.Cells, FT.Mesh(m))))]))


			 fun makeCells vars = (case vars
						of [v] => AST.E_ExtractFemItem(v, cellSeqType, (FemOpt.Cells, FT.Mesh(m)))
						 | _ => raise Fail "impossible argument to method got pass type checking")



			in
			 ([(numCell, numCellFuncVar), (cells', cellFuncVar')], [numCellFun, cellsFun', cellsFun], [(cells, makeCells, cellTypeFuncTy)])
			end
		   (*end case*))

	       fun makeDescendentFemTypes(env, femType) =
		   let
		    val constants = []
		    val methods = []
		    val methodRefs = []
		    val name = FT.nameOf femType


		   in
		    (case femType
		      of FT.Mesh(m) =>
			 let
			  val mapDim = FT.meshDim m
			  val meshMapDim = FT.meshMapDim m
			  val shape = FT.dataShapeOf femType
			  
			  val cellName = Atom.atom ("cell(" ^ (Atom.toString name) ^ ")" )
						   (*pos name*)
			  val cellData = FT.cellOf femType
			  val cellTy = Ty.T_Fem(cellData, NONE) (* no name*)

			  val transformDofs = Atom.atom "transformDofs"
			  val transformDofsTy = Ty.T_Tensor(Ty.Shape((List.map Ty.DimConst shape))) (*careful*)
			  val transformsDofsFunTy = Ty.T_Fun([cellTy], transformDofsTy)
			  fun makeTransformDofs vars = (case vars
							 of [v] =>
							    let
							     val getMesh =  AST.E_ExtractFem(v, femType)
							     val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, FT.MeshCell(m)))
							     val getTensor = AST.E_ExtractFemItem2(getMesh, getCellInt, Ty.T_Int, transformDofsTy, (FemOpt.ExtractDofs, FT.Mesh(m)))
							    in
							     getTensor
							    end
							  | _ => raise Fail "impossible argument to method got passed type checking"
						       (* end case*))

			  val transform = Atom.atom FemName.transform
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
								  AST.E_FemField(getMesh, SOME(getCellInt), transformFieldTy, FemOpt.Transform, NONE)
								 end
							       | _ => raise Fail "impossible argument to method got passed type checking"
							    (* end case*))


			  val transformFuncs = [(transformDofs, makeTransformDofs, transformsDofsFunTy), (transform, makeTransformFieldFunc, transformFieldFunc)]
			  (* TODO: The naming in this area of the compiler is really inconsistent*)

			  val env' = Env.insertNamedType(env, cxt, cellName, cellTy, constants, methods, transformFuncs)
							 
			 in
			  env'
			 end
		       | _ => env
		    (* end case *))
		   end


	       fun makeFemType femTyDef file =
		   let
		    val (femType, constants) : (FT.femType * Atom.atom option) option *  (Types.ty * ConstExpr.t * string) list = validateFemType femTyDef file
		    val (methodRefs, methods, newVars) = Option.getOpt (Option.map (makeFemMethods femTyDef file) femType, ([],[],[]))
		    val constants' = reformatConstants constants
		    val envI =  Option.map (fn x => Env.insertNamedType(env, cxt, tyName, Ty.T_Fem(x), constants', methodRefs, newVars)) femType

		    val argsOption = (case (envI, femType)
				       of (SOME(envI'), SOME((femTy', _))) => SOME((envI', femTy'))
					| _ => NONE)
		    val env' = Option.map makeDescendentFemTypes argsOption
		   in
		    (case env'
		      of SOME(env'') => (OTHERS(methods), env'')
		       | NONE => (ERROR, env)
		    (* end case *))
		    
		   end

	       fun makeBasicNamedType astTyDef file = let
		(*NOTE: error handling here is borked*)
		(*Step 1: Check the type definition*)
		val _ = print("ummm 1\n");
		val astTy = CheckType.check(env, cxt, astTyDef)
		val _ = print("umm 2\n");
		(* We enforce that this must be value type for now:*)
		(*bad way to do this:*)
		val isValType  = if TU.isValueType astTy
				 then true
				 else ((err (cxt, [
					     S "Declared a type", A(tyName), S " with a definition ", TY(astTy), S " that is not a value type or a fem type."
				       ])); false)
		(*Step 2: Build the constructor and add to the env*)
		val thisTy = Types.T_Named(tyName, astTy)
		val conFnTy = Ty.T_Fun([astTy], thisTy)
		val funVar = Var.new (tyName, span, AST.FunVar, conFnTy)
		val param = Var.new (Atom.atom "arg0", span, AST.FunParam, astTy)
		val body = AST.S_Block([AST.S_Return(AST.E_Coerce{srcTy=astTy, dstTy=thisTy, e=AST.E_Var(param,span)})])
		val astConFun = AST.D_Func(funVar, [param] , body)
		val _ = (E.checkForRedef (env, cxt, tyName); E.insertFunc(env, cxt, tyName, funVar)) (*finish creating constructor and adding to env*)

			  
		(*Step 3: Add methods: For now, we just add .unpack*)
		val unpack = Atom.atom "unpack"
		val deConFnTy = Ty.T_Fun([thisTy], astTy)
		val unpackVar =  Var.new (unpack, span, AST.FunVar, deConFnTy)
		val param' = Var.new (Atom.atom "arg0", span, AST.FunParam, thisTy)
		val body' = AST.S_Block([AST.S_Return(AST.E_Coerce{srcTy = thisTy, dstTy=astTy, e=AST.E_Var(param',span)})])
		val astDeConFun = AST.D_Func(unpackVar, [param'], body')
		(* Do Some fem or type specific stuff ... *)
		val methods = [(unpack, unpackVar)]
		val funDcls = [astConFun, astDeConFun]

		(*Step 4: Add constants*)
		val _ = print("Hummm!\n")
		val constants = reformatConstants (Option.valOf ((Option.map (fn x => CF.parseConstants(env, cxt, tyName, x))) (Option.mapPartial (fn x => CF.loadJson(x, cxt)) file)) handle exn => [])
		val _ = print("Hummm1!\n")
		val env' = Env.insertNamedType(env, cxt, tyName, thisTy, constants, methods, [])
	       in
		if isValType
		then (OTHERS(funDcls), env')
		else (ERROR, env)
	       end

	      in
	       if isFemType
	       then (case file
		      of SOME(file') => makeFemType tyDef file'
		       | _ => (err (cxt, [S "FemType declared without file"]); (ERROR, env)))
	       else makeBasicNamedType tyDef file

					  
	      end
	      
	
	
            | PT.GD_FieldFunc(ty, bindF, bindX, body) => (case CheckType.check(env, cxt, ty)
                 of ty' as Ty.T_Field{diff, dim, shape} => let
                      val f = Var.new(#tree bindF, #span bindF, Var.GlobalVar, ty')
                      val xTy = Ty.T_Tensor(Ty.Shape[dim])
                      val resTy = Ty.T_Tensor shape
(* QUESTION: should we check the arity of dim? *)
                      val x = Var.new(#tree bindX, #span bindX, Var.FunParam, xTy)
                      val env' = Env.insertLocal(
(* QUESTION: should we have a new kind of scope for these field functions? *)
                            Env.functionScope(env, resTy, #tree bindF),
                            cxt, #tree bindX, x)
                      val bodyTy = CheckExpr.check (env', cxt, body)
                      val _ = CheckFieldFunc.checkPT (env', cxt, body)
                      in
                        case Util.coerceType (resTy, bodyTy)
                         of SOME e => 
                            ( CheckFieldFunc.checkAST (env, cxt, e);                          
                             (
                                OTHER(AST.D_DiffFunc(f, [x], e)),
                                Env.insertGlobal(env, cxt, #tree bindF, f)
                              )
                            )
                          | NONE => (
                              err (cxt, [
                                  S "type of r.h.s. definition for '", A(#tree bindF),
                                  S "' does not match declaration",
                                  S "  expected: ", TY resTy, S "\n",
                                  S "  found:    ", TY(#2 bodyTy)
                                ]);
                              (ERROR, env))
                        (* end case *)
                      end
                  | ty' => (
                      err (cxt, [S "expected field type for '", A(#tree bindF), S "'"]);
                      (ERROR, env))
                (* end case *))
          (* end case *))

    fun check (defs, env, cxt, globs) = let
          fun chk (env, [], cdecls, idecls, gdecls) = {
                  const_dcls = List.rev cdecls,
                  input_dcls = List.rev idecls,
                  other_dcls = List.rev gdecls,
                  env = env
                }
            | chk (env, dcl::dcls, cdecls, idecls, gdecls) = (
                case chkDcl (defs, env, cxt, dcl)
                 of (CONST dcl, env) => chk (env, dcls, dcl::cdecls, idecls, gdecls)
                  | (INPUT dcl, env) => chk (env, dcls, cdecls, dcl::idecls, gdecls)
                  | (OTHER dcl, env) => chk (env, dcls, cdecls, idecls, dcl::gdecls)
		  | (OTHERS dcls', env) => chk (env, dcls, cdecls, idecls, dcls'@gdecls)
                  | (ERROR, env) => chk (env, dcls, cdecls, idecls, gdecls)
		  | (IGNORE, env) => chk (env, dcls, cdecls, idecls, gdecls)
                (* end case *))
          in
            chk (env, globs, [], [], [])
          end

  end
