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
    structure FO = FemOpt
    structure FN = FemName
    structure N = BasisNames
    structure CF = CheckTypeFile
    structure BV = BasisVars


		     
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
			     | (_, SOME(_)) => ((err (cxt, [
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
		     of ((func as FT.Func(f), n), PT.T_Func(space)) =>
			let
			 val funcName = FT.nameOf func
			 val spaceFunc = Atom.atom "space"
			 val space' = FT.spaceOf(FT.Func(f))
			 val spaceName = FT.nameOf space'
			 val mesh = FT.meshOf func
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

			 (* OH GOD THIS NEEDS TO BE CLEANED UP*)
			 val funcArgTy = Ty.T_Named(tyName, Ty.T_Fem(femInfo))
			 val funcCellTy = Ty.T_Fem(func, SOME(funcName))
			 val FT.Mesh(meshData) = FT.meshOf func
			 val meshCell = FT.MeshCell(meshData)
			 val funcCellVal = FT.FuncCell(f)
			 val meshCellTy = Ty.T_Fem(meshCell, SOME(meshName))
			 val funcCellTy = Ty.T_Fem(funcCellVal, SOME(funcName))
						  
			 val funcCell = Atom.atom "funcCell"
			 val funcCellType = Ty.T_Fun([funcArgTy, meshCellTy], funcCellTy)
			 val funcCellFuncVar = Var.new (funcCell, span, AST.FunVar, funcCellType)

			 val param = Var.new (Atom.atom "arg0", span, AST.FunParam, funcArgTy)
			 val param' = Var.new (Atom.atom "arg1", span, AST.FunParam, meshCellTy)
			 (* optional debug information could go here.*)
			 val funcCellBody = AST.S_Return(AST.E_LoadFem(funcCellVal,
								       SOME(AST.E_Var(param, span)),
								       SOME(AST.E_ExtractFemItem(
									      AST.E_Var(param', span),
									      Ty.T_Int,
									      (FO.CellIndex, meshCell)))))
			 val funcCellFunc = AST.D_Func(funcCellFuncVar,
						       [param, param'],
						       funcCellBody)


			in
			 ([(funcCell, funcCellFuncVar), (spaceFunc, spaceFuncVar)],
			  [funcCellFunc, spaceFun],
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
			 ([(meshFunc, meshFuncVar)], [meshFun],[])
			end
		      | ((meshData as FT.Mesh(m), n), PT.T_Mesh) =>
			let
			 val meshName = FT.nameOf meshData
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
			 				(*AST.E_ExtractFemItem(AST.E_Var(param', span)*)
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


			in
			 ([(numCell, numCellFuncVar), (cells', cellFuncVar'), (refCellName, refCellFuncVar)], [numCellFun, cellsFun', cellsFun,refCellFunc], [(cells, makeCells, cellTypeFuncTy)])
			end

		   (*end case*))

	       fun makePrim(prim, argTys) =
		   let
		    val AST.E_Prim(var, metargs, exprs, result) = prim
		    val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf var)
		    val exprs' = Option.valOf(Unify.matchArgs(domTy, exprs, argTys))
		    val _ = Unify.matchType(result, rngTy)

		   in
		    AST.E_Prim(var, tyArgs, exprs', result)
		   end

	       fun makePrim'(var, args, argTys, resultTy) =
		   let
		    val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf var)
		    val args' = Option.valOf(Unify.matchArgs(domTy, args, argTys))
		    val _ = Unify.matchType(resultTy, rngTy)
		   in
		    AST.E_Prim(var, tyArgs, args', rngTy)
		   end
	       fun makeAnd v1 v2 = makePrim'(BV.and_b, [v1,v2], [Ty.T_Bool, Ty.T_Bool], Ty.T_Bool)
	       fun makeAnds([v1,v2]) = makeAnd v1 v2
		 | makeAnds (v::vs) = makeAnd v (makeAnds vs)
	       fun makeRefCellInside(env, cxt, span, dim, refCellInfo, posVar, meshVar) =
		   let
		    val insideVec = Ty.vecTy dim
		    val FemData.RefCellData({ty=refCellClass, eps,...}) = refCellInfo
		    val epsExpr = AST.E_Lit(Literal.Real(eps))
		    val one = AST.E_Lit(Literal.Real(RealLit.one))
		    val onePlusEps = one (*TODO: FIX EPSILON and add me*)
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
		      of FemData.KSimplex =>
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
		       | FemData.KCube =>
			 let
			  val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_norm_t)
			  val norm = makePrim(AST.E_Prim(BV.op_norm_t, tyArgs, [posVar'], Ty.realTy), [insideVec])
			  val endTest = greaterThanTest onePlusEps norm
			 in
			  endTest
			 end)
		   end

	       fun makeNewtonInversesBody(env, cxt, span, refCellClass, meshData, newtonTol, newtonAttempts, contraction, killAfterTol, posExpr, cellIntExpr, meshExpr, insideFunction) =
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
		    val startPosition = (case refCellClass
					  of FemData.KCube => AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(RealLit.zero false))), insideVec)
					   | FemData.KSimplex =>
					     let
					      val digits = List.tabulate(100, fn x => 3)
					      val exp = IntInf.fromInt (1)
					      val third = RealLit.fromDigits({isNeg = false, digits = digits, exp = exp})
					     in
					      AST.E_Tensor(List.tabulate(dim, fn x => AST.E_Lit(Literal.Real(third))), insideVec)
					     end
					(* end case*))
		    val invVar = if dim = 2
				 then BV.fn_inv2_f
				 else if dim = 3
				 then BV.fn_inv3_f
				 else raise Fail "invalid dim;" (*TODO: this resgtriction should be left as I should be able to do arbitrary inverses of a matrix*)

		    val newPosVar = Var.new(Atom.atom "xn", span, AST.LocalVar, insideVec)
		    val newPosVarExpr = AST.E_Var((newPosVar, span))
		    val initPos = AST.S_Assign((newPosVar, span), startPosition)
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
						      
		    val deltaNorm = makePrim'(BV.op_norm_t, [AST.E_Var((updateVar, span))], [insideVec], Ty.realTy)
		    val deltaNormTest = makePrim'(BV.gte_rr, [tolExpr, deltaNorm], [Ty.realTy, Ty.realTy], Ty.T_Bool)
		    val insideTest = insideFunction([meshExpr, newPosVarExpr])
		    (* val combinedTest = AST.Andalso(deltaNormTest, insideTest) (*note: we might want to combined them...... it might be better to do delta test then inside test????*) *)
		    val succesReturn = AST.S_Return(newPosVarExpr) (*TODO: return the itterated thing*)
		    val nans = makePrim'(BV.nan, [], [], insideVec) (*return nans to signal failure*)
		    val failReturn = AST.S_Return(nans)
		    val ifStm = if true
				then AST.S_IfThenElse(deltaNormTest, (*we might want to change this test to be this and that, rather than a nested if.*)
						 AST.S_Block([AST.S_IfThenElse(insideTest,
									       AST.S_Block([succesReturn]), (*converged and arrived*)
									       AST.S_Block([failReturn]))]), (*converged and failed; somehow*)
						 AST.S_Block([]))
				else AST.S_IfThenElse(makeAnd deltaNormTest insideTest, AST.S_Block([succesReturn]), AST.S_Block([failReturn]))
		    val _ = print("The contraction is:"^(Bool.toString contraction) ^ "\n")
						
		    val completeBody =
			if contraction
			then
			 let
			  val _ = print("Contraction yes\n")
			  val probeD = makePrim'(BV.op_probe, [invDTransformField, startPosition], [dTransformFieldTy, insideVec], insideMat)
			  val matVar = Var.new(Atom.atom "A", span, AST.LocalVar, insideMat)
			  val probeDAssign = AST.S_Assign((matVar,span), probeD)
			  val dotFields = makePrim'(BV.op_inner_tf,  [AST.E_Var(matVar, span), transformFieldModPos], [insideMat, transformFieldTy], transformFieldTy)
			  val probeUpdate = makePrim'(BV.op_probe, [dotFields, newPosVarExpr], [transformFieldTy, insideVec], insideVec)
			  val updateDeltaStm = AST.S_Assign((updateVar, span), probeUpdate)
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
			  val _ = print("Contraction not")
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
								  AST.E_FemField(getMesh,getMesh, SOME(getCellInt), transformFieldTy, FemOpt.Transform, NONE)
								 end
							       | _ => raise Fail "impossible argument to method got passed type checking"
							    (* end case*))




			  val refCellInside = Atom.atom (FemName.refCellIsInside)
			  val refCellInsideTy = Ty.T_Fun([refCellTy, dimSizeVec], Ty.T_Bool)
			  fun makeRefCellInsideFunc vars = (case vars
							     of [v1,v2] => makeRefCellInside(env, cxt, span, mapDim, refCellInfo, v2, v1)
							      | _ => raise Fail "impossble got through type checker"
							   (* end case*))

			  fun mkVar (var) = AST.E_Var(var,span)
			  val (hiddenFuncAtom, hiddenFuncVar, newtonFunc) = let (*TODO: USE let a lot here too refactor this code... god this code is ugly*)
			   val hiddenFuncAtom = Atom.atom  (FemName.hiddenNewtonInverse)


			   val FemData.RefCellData({ty=refCellClass, eps, newtonControl}) = refCellInfo (*move this outside? Eh.*)
			   val {contraction, itters, newtonTol, killAfterTol} = newtonControl
			   val _ = print("The contraction is:"^(Bool.toString contraction) ^ "\n")
			   val meshData = m
			   val posTy = Ty.vecTy mapDim
			   val newtonTol = AST.E_Lit(Literal.Real(newtonTol)) (*TODO fix me with ref cell info*)
			   val newtonAttempts = AST.E_Lit(Literal.intLit itters) (*TODO fix me with ref cell info*)


			   val posParam = Var.new (Atom.atom "pos", span, AST.FunParam, posTy)
			   val cellIntParam = Var.new (Atom.atom "cellInt", span, AST.FunParam, Ty.T_Int)
			   val meshParam = Var.new (Atom.atom "mesh", span, AST.FunParam, meshTy)

			   val funTy = Ty.T_Fun([Ty.vecTy mapDim, Ty.T_Int, meshTy], posTy)
			   val hiddenFuncVar = Var.new (hiddenFuncAtom, span, AST.FunVar, funTy)
						       
			   val posExpr = mkVar posParam
			   val cellIntExpr = mkVar cellIntParam
			   val meshExpr = mkVar meshParam
						
			   val insideFunc = makeRefCellInsideFunc

			   val body = makeNewtonInversesBody(env, cxt, span, refCellClass, meshData, newtonTol, newtonAttempts, contraction, killAfterTol, posExpr, cellIntExpr, meshExpr, insideFunc)
			   val newtonFunc = AST.D_Func(hiddenFuncVar, [posParam, cellIntParam, meshParam], body)
			  in
			   (hiddenFuncAtom, hiddenFuncVar, newtonFunc)
			  end

			  val invTransformSpec =
			      let
			       val invTransformField = Atom.atom FemName.invTransform
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
			       (invTransformField, makeInvTransformFieldFunc, invTransformFieldFuncTy)
			      end

			  val transformFuncs = [(transformDofs, makeTransformDofs, transformsDofsFunTy), (transform, makeTransformFieldFunc, transformFieldFunc), invTransformSpec]
			  val cellMethods = [(hiddenFuncAtom, hiddenFuncVar)]
			  (* TODO: The naming in this area of the compiler is really inconsistent*)

			  val env' = Env.insertNamedType(env, cxt, cellName, cellTy, constants, cellMethods, transformFuncs)

			  
			 

			  val refCellTransformFuncs = [(refCellInside, makeRefCellInsideFunc, refCellInsideTy)]
			  val env'' = Env.insertNamedType(env', cxt, refCellName, refCellTy, constants, methods, refCellTransformFuncs)			
							 
			 in
			  (env'', [newtonFunc])
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
							    val getSpace = AST.E_ExtractFem(getFunc,space)
							    (*might want to extract it now.*)
							    val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
							    val getIndexi = AST.E_ExtractFemItem2(getSpace, (getCellInt), Ty.T_Int, dofIndexSeq, (FemOpt.ExtractIndices, space))
							    val getTensor = AST.E_ExtractFemItem2(getFunc, (getIndexi), dofIndexSeq, funcDofsTy, (FemOpt.ExtractDofsSeq, femType))
							   in
							    getTensor
							   end
							 | _ => raise Fail "impossible argument to method got passed type checking")


			  val refField = Atom.atom FemName.refField
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
						
							    val _ = print("Called!")
							   in
							    AST.E_FemField(getFunc,getSpace, SOME(getCellInt),refFieldTy,FO.RefField, NONE )
							   end
							 | _ => raise Fail "impossible argument to method got passed type checking")

			  val transformRefFieldFunc = Atom.atom (FemName.trf)
			  val transformFieldTy = Ty.T_Field({diff = inf, dim = dimConst, shape = Ty.Shape([dimConst])})
			  fun makeTransformRefFieldFunc vars = (case vars
								 of [v] =>
								    let
								     val getFunc = AST.E_ExtractFem(v, femType)
								     val getSpace = AST.E_ExtractFem(getFunc, space)
								     val getMesh = AST.E_ExtractFem(getSpace, mesh)
								     val getCellInt= AST.E_ExtractFemItem(v, Ty.T_Int, (FemOpt.CellIndex, funcCellData))
								     (*get mesh cell env, get the function, make the fields*)
								     val SOME(meshCellEnv) = Env.findTypeEnv(env, meshCellEnvName)
								     val SOME(newtonMethod) = TypeEnv.findMethod(meshCellEnv, Atom.atom FemName.hiddenNewtonInverse)
								     (*get the basis*)
								     val field1 = AST.E_FemField(getFunc,getSpace, SOME(getCellInt), refFieldTy, FO.RefField, NONE )
								     val field2 = AST.E_FemField(getMesh, getMesh, SOME(getCellInt), transformFieldTy, FemOpt.InvTransform, NONE)
								     val field = makePrim'(BV.comp, [field1, field2], [refFieldTy, transformFieldTy], refFieldTy)
								    in
								     field
								    end
								 | _ => raise Fail "impossible that this got past the typechecker"
							     )


			  val transformFuncs = [(functionDofs, makeFunctionDofs, functionDofsFunTy),
						(refField, makeRefFieldFunc, refFieldFunc),
						(transformRefFieldFunc, makeTransformRefFieldFunc, refFieldFunc)]
			  val env' = Env.insertNamedType(env, cxt, cellName, cellTy, constants, methods, transformFuncs)
			 in
			  (env', [])
			 end

		       | _ => (env, [])
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
		      of SOME(env'', methods') => (OTHERS(List.@(methods', methods)), env'')
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
