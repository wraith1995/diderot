(* check-expr.sml
 *
 * The typechecker for expressions.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure CheckExpr : sig

  (* type check an expression *)
    val check : Env.t * Env.context * ParseTree.expr -> (AST.expr * Types.ty)


  (* type check a list of expressions *)
    val checkList : Env.t * Env.context * ParseTree.expr list -> (AST.expr list * Types.ty list)

  (* type check an iteration expression (i.e., "x 'in' expr"), returning the iterator
   * and the environment extended with a binding for x.
   *)
    val checkIter : Env.t * Env.context * ParseTree.iterator -> ((AST.var * AST.expr) * Env.t)
								  
  (* type check a dimension that is given by a constant expression *)
    val checkDim : Env.t * Env.context * ParseTree.expr -> IntLit.t option

  (* type check a tensor shape, where the dimensions are given by constant expressions *)
    val checkShape : Env.t * Env.context * ParseTree.expr list -> Types.shape

  (* `resolveOverload (cxt, rator, tys, args, candidates)` resolves the application of
   * the overloaded operator `rator` to `args`, where `tys` are the types of the arguments
   * and `candidates` is the list of candidate definitions.
   *)
    val resolveOverload : Env.context * Atom.atom * Types.ty list * AST.expr list * Var.t list
          -> (AST.expr * Types.ty)

  end = struct

    structure PT = ParseTree
    structure L = Literal
    structure E = Env
    structure Ty = Types
    structure BV = BasisVars
    structure TU = TypeUtil
    structure CO = CheckOverload
    structure TE = TypeEnv
    structure FT = FemData
  (* an expression to return when there is a type error *)
    val bogusExp = AST.E_Lit(L.Int 0)
    val bogusExpTy = (bogusExp, Ty.T_Error)

    fun err arg = (TypeError.error arg; bogusExpTy)
    val warn = TypeError.warning

    datatype token = datatype TypeError.token

  (* mark a variable use with its location *)
    fun useVar (cxt : Env.context, x) = (x, #2 cxt)

  (* strip any marks that enclose an expression and return the span and the expression *)
    fun stripMark (_, PT.E_Mark{span, tree}) = stripMark(span, tree)
      | stripMark (span, e) = (span, e)

  (* resolve overloading: we use a simple scheme that selects the first operator in the
   * list that matches the argument types.
   *)
    fun resolveOverload (cxt, rator, argTys, args, candidates) =
	(case CO.chkOverload(cxt, rator, argTys, args, candidates)
	  of
	     CO.Overload(res) => res
	   | CO.Error(res) => err res
	   | CO.NoOverload => err(cxt, [
                                      S "Given an operator or function", A(rator), S" with argument types", TYS(argTys),
				      S " was unable to locate an overloaded operator."
                                    ])

	(* end case*))
  (* check the type of a literal *)
    fun checkLit lit = (AST.E_Lit lit, TypeOf.literal lit)

  (* check the well-formedness of a spatial query `e`, which has already been typechecked *)
    fun checkSpatialQuery (env, cxt, e, tyArgs, rngTy) = (case Env.strandTy env
           of SOME(strandTy, sEnv) => (case StrandEnv.findPosVar sEnv
					of SOME p => let
		      val (posTyVar, Ty.TYPE tv) = (case tyArgs
						    of [a] => (NONE, a)
						     | [a,b] => (SOME(a), b)
						     | _ => raise Fail "impossible")
						    
		      (*FIXME: add proper type errors here.*)
                      fun result dim = (
                            StrandEnv.recordSpaceDim (sEnv, dim);
                            (e, rngTy))
		      val posType = TU.prune (Var.monoTypeOf p)
                      in
                      (* instantiate the query's type to the strand type *)
                       ignore (Unify.matchType (Ty.T_Var tv, strandTy));

		       (case posTyVar
			 of NONE => ()
			  | SOME(Ty.TYPE tv') => ignore (Unify.matchType (Ty.T_Var tv', posType)));
                      (* check that the strand supports spatial queries *)
                        case StrandEnv.getSpaceDim sEnv
                         of SOME _ => (e, rngTy)  (* we have already processed a spatial query *)
                          | NONE => (
                              Env.recordProp (env, Properties.StrandCommunication);
                            (* check the type of the position; should be 1D, 2D, or 3D *)
                              case posType
                               of Ty.T_Tensor(Ty.Shape[]) => result 1
                                | Ty.T_Tensor(Ty.Shape[Ty.DimConst 2]) => result 2
                                | Ty.T_Tensor(Ty.Shape[Ty.DimConst 3]) => result 3
				| Ty.T_Fem(FT.MeshPos(meshData),_) => result (FemData.meshDim meshData)
                                | ty => err(cxt, [
                                      S "'expected one of real, vec2, vec3, or a meshPosition on a mesh of dim 2 or 3 for type of 'pos',\n",
                                      S "  but found: ", TY ty
                                    ])
                              (* end case *))
                        (* end case *)
                      end
                  | NONE => err(cxt, [
                        S "spatial queries require defining a 'pos' variable of suitable type"
                      ])
                (* end case *))
            | NONE => err(cxt, [
                  S "spatial queries are only allowed inside strands"
                ])
          (* end case *))

    (* check the type of an expression *)
    fun check (env, cxt, e) =
	(case e
           of PT.E_Mark m => check (E.withEnvAndContext (env, cxt, m))
            | PT.E_Cond(e1, cond, e2) => let
                val eTy1 = check (env, cxt, e1)
                val eTy2 = check (env, cxt, e2)
                in
                  case checkAndPrune(env, cxt, cond)
                   of (cond', Ty.T_Bool) => (case Util.coerceType2(eTy1, eTy2)
                         of SOME(e1', e2', ty) =>
                              if TU.isValueType ty
                                then (AST.E_Cond(cond', e1', e2', ty), ty)
                                else err (cxt, [
                                    S "result of conditional expression must be value type,\n",
                                    S "  but found ", TY ty
                                  ])
                          | NONE => err (cxt, [
                              S "types do not match in conditional expression\n",
                              S "  true branch:  ", TY(#2 eTy1), S "\n",
                              S "  false branch: ", TY(#2 eTy2)
                            ])
                        (* end case *))
                    | (_, Ty.T_Error) => bogusExpTy
                    | (_, ty') => err (cxt, [S "expected bool type, but found ", TY ty'])
                  (* end case *)
                end
            | PT.E_Range(e1, e2) => (case (check (env, cxt, e1), check (env, cxt, e2))
                 of ((e1', Ty.T_Int), (e2', Ty.T_Int)) => let
                      val resTy = Ty.T_Sequence(Ty.T_Int, NONE)
                      in
                        (AST.E_Prim(BV.range, [], [e1', e2'], resTy), resTy)
                      end
                  | ((_, Ty.T_Int), (_, ty2)) =>
                      err (cxt, [S "expected type 'int' on rhs of '..', but found ", TY ty2])
                  | ((_, ty1), (_, Ty.T_Int)) =>
                      err (cxt, [S "expected type 'int' on lhs of '..', but found ", TY ty1])
                  | ((_, ty1), (_, ty2)) => err (cxt, [
                        S "arguments of '..' must have type 'int', found ",
                        TY ty1, S " and ", TY ty2
                      ])
                (* end case *))
            | PT.E_OrElse(e1, e2) => checkCondOp (env, cxt, e1, "||", e2, AST.E_Orelse)
            | PT.E_AndAlso(e1, e2) => checkCondOp (env, cxt, e1, "&&", e2, AST.E_Andalso)
            | PT.E_BinOp(e1, rator, e2) => let
                val (e1', ty1) = check (env, cxt, e1)
                val (e2', ty2) = check (env, cxt, e2)
		fun resolveBinOp env rator =
		    (case Env.findFunc (env, rator)
                      of Env.PrimFun[rator] => let
                       val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf rator)
                      in
                       (case Unify.matchArgs(domTy, [e1', e2'], [ty1, ty2])
                        of SOME args => CO.chkPrim (AST.E_Prim(rator, tyArgs, args, rngTy), rngTy)
                         | NONE => err (cxt, [
					S "type error for binary operator ", V rator, S "\n",
					S "  expected: ", TYS domTy, S "\n",
					S "  found:    ", TYS[ty1, ty2]
                                       ])
		       (* end case *))
                      end
                       | Env.PrimFun ovldList =>
                         resolveOverload (cxt, rator, [ty1, ty2], [e1', e2'], ovldList)
                       | _ => raise Fail "impossible" (*Note: Not allowing things like a `lerp` b via parsing restrictions *)
                    (* end case *))
		fun checkSpecial checker =
		    (case checker (cxt, e1', ty1, e2', ty2)
		       of CO.Overload(result) => result
			| CO.Error((result)) => err result
			| CO.NoOverload =>  resolveBinOp env rator 
		    (*end case*))
                in
                  if Atom.same(rator, BasisNames.op_dot)
                  then checkSpecial CO.chkInnerProduct 
                  else if Atom.same(rator, BasisNames.op_outer)
                    then checkSpecial CO.chkOuterProduct
                  else if Atom.same(rator, BasisNames.op_colon)
                    then checkSpecial CO.chkColonProduct 
                  else resolveBinOp env rator
                end
            | PT.E_UnaryOp(rator, e) => let
                val eTy = check(env, cxt, e)
                in
                  case Env.findFunc (env, rator)  
                   of Env.PrimFun[rator] => let 
                        val (tyArgs, Ty.T_Fun([domTy], rngTy)) = TU.instantiate(Var.typeOf rator)
                        in
                          case Util.coerceType (domTy, eTy)
                           of SOME e' => CO.chkPrim  (AST.E_Prim(rator, tyArgs, [e'], rngTy), rngTy)
                            | NONE => err (cxt, [
                                  S "type error for unary operator ", V rator, S "\n",
                                  S "  expected: ", TY domTy, S "\n",
                                  S "  found:    ", TY (#2 eTy)
                                ])
                          (* end case *)
                        end
                    | Env.PrimFun ovldList => resolveOverload (cxt, rator, [#2 eTy], [#1 eTy], ovldList)
                    | _ => raise Fail "impossible"
                  (* end case *)
                end
            | PT.E_Apply(e, args) => let
                val (args, tys) = checkList (env, cxt, args)
                fun appTyError (f, paramTys, argTys) = err(cxt, [
                        S "type error in application of ", V f, S "\n",
                        S "  expected: ", TYS paramTys, S "\n",
                        S "  found:    ", TYS argTys
                      ])
                fun checkPrimApp f = if Var.isPrim f
                      then (case TU.instantiate(Var.typeOf f)
                         of (tyArgs, Ty.T_Fun(domTy, rngTy)) => (
                              case Unify.matchArgs (domTy, args, tys)
                               of SOME args => CO.chkPrim (AST.E_Prim(f, tyArgs, args, rngTy), rngTy)
                                | NONE => appTyError (f, domTy, tys)
                              (* end case *))
                          | _ => err(cxt, [S "application of non-function/field ", V f])
                        (* end case *))
                      else raise Fail "unexpected user function"
              (* check the application of a user-defined function *)
                fun checkFunApp (cxt, f) args' tys' = if Var.isPrim f
                      then raise Fail "unexpected primitive function"
                      else (case Var.monoTypeOf f
                         of Ty.T_Fun(domTy, rngTy) => (
                              case Unify.matchArgs (domTy, args', tys')
                               of SOME args => (AST.E_Apply(useVar(cxt, f), args, rngTy), rngTy)
                                | NONE => appTyError (f, domTy, tys')
                              (* end case *))
                          | _ => err(cxt, [S "application of non-function/field ", V f])
                           (* end case *))

                fun checkMethodApp (cxt, f) args' arg tys' ty = if Var.isPrim f
                      then raise Fail "unexpected primitive function"
                      else (case Var.monoTypeOf f
                         of Ty.T_Fun(domTy, rngTy) => (
                              case Unify.matchArgs (domTy, arg::args', ty::tys')
                               of SOME args => (AST.E_Apply(useVar(cxt, f), args, rngTy), rngTy)
                                | NONE => err(cxt, [
					      S "type error in application of method,", V f, S ", defined in type ", TY(ty),
					      S "  expected: ", TYS domTy, S "\n",
					      S "  found:    ", TYS (ty::tys')
                      ])
                              (* end case *))
                          | _ => err(cxt, [S "application of non-function/field ", V f])
                           (* end case *))
		fun checkFunTyMatchSingle cxt field funTy args argTys funcApp =
		    let
		     val Ty.T_Fun(domTys, rngTy) = funTy
		     (* check if this is a field*)
		    in
		     case Unify.matchArgs (domTys, args, argTys)
		      of SOME(args) => (funcApp args, rngTy)
		       | NONE => err(cxt, [
				     S "type error in application of method ", A field,
				     S "  expected: ", TYS domTys, S "\n",
				     S "  found:    ", TYS argTys
				    ])
		    end
		fun checkFemConstruction cxt loadFem rngTy domTy args argTys=
		    (case Unify.matchArgs (domTy, args, argTys)
		      of SOME(args) => (case (Env.currentScope env, FemData.baseFem (#1 loadFem))
					 of (Env.InputScope, true) => (AST.E_LoadFem(loadFem), rngTy)
					  | (_, false) => (AST.E_LoadFem(loadFem), rngTy)
					  | (_, true) => err(cxt, [S "cannot construct base FEM data outside of input init!"])
				       (*end case*))
		       | NONE => err(cxt, [
				     S "type error in application of fem creator ",
				     S "  expected: ", TYS domTy, S "\n",
				     S "  found:    ", TYS argTys
				    ])
		    (* end case*))
                fun checkFieldApp (e1', ty1) = (case (args, tys)
                       of ([e2'], [ty2]) => let
                            val (tyArgs, Ty.T_Fun([fldTy, domTy], rngTy)) =
                                  TU.instantiate(Var.typeOf BV.op_probe)
                            fun tyError () = err (cxt, [
                                    S "type error for field application\n",
                                    S "  expected: ", TYS[fldTy, domTy], S "\n",
                                    S "  found:    ", TYS[ty1, ty2]
						 ])
			    fun innerTyError (exp, fon) = err (cxt, [S "type error for field application of probe to pos\n",
								     S " expected: ", TYS exp, S "\n",
								     S " found: ", TYS fon]
							      )


			    val preCompose = (case TU.prune ty2
					       of Ty.T_Tensor _ => NONE
						| Ty.T_Fem(data as FemData.MeshPos(m), a) =>
						  let
						   val mesh = FemData.Mesh(m)
						   val meshTy = Ty.T_Fem(mesh, a)
						   val meshData = AST.E_ExtractFem(e2', mesh)
						   val cellInt = AST.E_ExtractFemItem(e2', Ty.T_Int, (FemOpt.CellIndex, data))
						   val dim = FemData.underlyingDim data
						   val transformFieldTy = Ty.vecField(NONE, dim)
						   val vecTy = Ty.vecTy dim
						   val femT = AST.E_FemField(meshData, meshData, SOME(cellInt), transformFieldTy,
									     FemOpt.Transform, NONE)
						   val (tyArgs', Ty.T_Fun([fldTy1, fldTy2], fldTy3)) =
						       TU.instantiate(Var.typeOf BV.comp)
						   fun tye () = innerTyError([fldTy1, ty2], [ty1, ty2])
								     
						   val (origField, _) = if Unify.equalType(ty1, fldTy1)
									then (case Util.coerceType(fldTy1, (e1', ty1))
									       of SOME(e1') => (e1', fldTy1)
										| NONE => tye ()
									     (* end case *))
									else tye ()
						   val (otherField, _) = if Unify.equalType(transformFieldTy, fldTy2)
									 then (case Util.coerceType(fldTy2, (femT, transformFieldTy))
										of SOME(ofield) => (ofield, fldTy2)
										 | NONE => bogusExpTy (*error above*)
									      (* end case*))
									 else bogusExpTy (*error above*)
						   val compose = AST.E_Prim(BV.comp, tyArgs', [origField, otherField], fldTy3)

						   val getRefPos = AST.E_ExtractFemItem(e2', vecTy, (FemOpt.RefPos, data))

						   val (doProbe, retTy) =
						       if Unify.equalType(fldTy, fldTy3)
						       then (case Util.coerceType(domTy, (getRefPos, vecTy))
							      of SOME(getRefPos') => (AST.E_Prim(BV.op_probe, tyArgs,
												[compose, getRefPos'],
												rngTy), rngTy)
							       | NONE => bogusExpTy (*error above*)
									   
							    (* end case *))
						       else bogusExpTy (*error above*)
						  in
						   SOME(doProbe, retTy)
						  end
						| _ => NONE
					     (* end case *))
                       in
			if Option.isSome preCompose
			then Option.valOf preCompose
			else
                         if Unify.equalType(fldTy, ty1)
                         then (case Util.coerceType(domTy, (e2', ty2))
                                of SOME e2' => (AST.E_Prim(BV.op_probe, tyArgs, [e1', e2'], rngTy), rngTy)
                                 | NONE => tyError()
                              (* end case *))
                         else tyError()
                       end
					      
                        | _ => err(cxt, [S "badly formed field application"])
					       (* end case *))
                in
                 case stripMark(#2 cxt, e)
		  of (span, PT.E_Select(exp, field)) => (*QUESTION: Is it possible the select is hidden behind something?*)
		     (
		       let
			val (astExp, astTy) = check (env, cxt, exp)
			val tyName = (case astTy
				       of Ty.T_Named(tyName, tyDef) => SOME(tyName)
					| Ty.T_Fem(data, _) => SOME( FT.envNameOf data)
					| _ => (err (cxt, [TY(astTy), S " is not a named type and so has no members."]); NONE)
				     (*end case*))

			val tyEnv = (Option.mapPartial (fn name => case Env.findTypeEnv(env, name)
								    of SOME(s) => SOME(s)
								     | NONE =>  (err (cxt, [S "There is no type named", A(name)]); NONE) )
						       tyName)
			val method : ((AST.var option) * ( (((AST.expr list) -> (AST.expr) ) * Ty.ty) option)) option
			    = Option.mapPartial (fn env => case TypeEnv.findMethod(env, field)
							    of SOME(s) => SOME ((SOME(s), NONE))
							     | NONE =>
							       (case TypeEnv.findHiddenVar(env, field)
								 of SOME(s) => SOME(NONE, SOME(s))
								  | NONE => ( (err (cxt, [S"The type named", A (TypeEnv.findName env),
										     S " does not have a member named",
										     A(field)])); NONE))) tyEnv
									
		       in
			case method 
			 of SOME(SOME(s), _) => checkMethodApp((#1 cxt, span), s) (args) astExp ( tys) astTy (*Note: This could result in a cryptic error message*)
			  | SOME(NONE, SOME(f, Ty.T_Fun([v], Ty.T_Field(fieldParams)))) => (*This is the fields case*)
			    let
			     (*Then we check if this is a valid application of the field.*)
			     val fieldExp = f [astExp]
			    in
			     (case tyName
			       of NONE => bogusExpTy (* THIS IS REDUNDANT*)
				| SOME(_) => checkFieldApp(fieldExp,Ty.T_Field(fieldParams)))
			    end
			  | SOME(NONE, SOME(f, funTy)) => checkFunTyMatchSingle cxt field funTy (astExp :: args) (astTy :: tys) f
			 | NONE => bogusExpTy
		       end
		     )
                    | (span, PT.E_Var f) => (case Env.findVar (env, f)
                         of SOME f' => checkFieldApp (
                              AST.E_Var(useVar((#1 cxt, span), f')),
                              Var.monoTypeOf f')
                          | NONE => (case Env.findFunc (env, f)
				      of Env.PrimFun[] => (
				       (case Env.findTypeEnv(env, f)
					 of NONE => err(cxt, [S "unknown function ", A f])
					  | SOME(tenv) =>
					    let (*args stored in args; the name f is the name of the type that we are clearing.*)
					     val tyDf = TE.findDef tenv
								   (*is FEM -> determine the underlying type if any that I need to match?*)
								   (*What I need: constructor type name : f, constructor tyargs : d, acceptable constructor ty args*)
					     (* I want to match the exact named type or the definition via the files???*)
					     val (data, tyArgs,name) = (case tyDf
								 of Ty.T_Fem(data,NONE) => (data, [],[])
								  | Ty.T_Fem(data, SOME(name)) =>
								    let
								     val tenv' = Option.valOf (Env.findTypeEnv(env, name))
								     val lowerTyDef = TE.findDef tenv'
								    in
								     (data, [Ty.T_Named(name, lowerTyDef)], [name])
								    end
								  | _ => raise Fail "non-FEM constructor")
					     val namedArg1 = (case args
							      of [] => NONE
							       | a::_ => SOME(a))
					     val namedArg2 = (case args
							      of [] => NONE
							       | [a,b] => SOME(b)
							       | _ => NONE)

					     val rngTy = Ty.T_Named(f, tyDf)
								 
					    in
					     checkFemConstruction cxt (data, namedArg1, namedArg2) rngTy tyArgs args tys
					    end
				       (* end case *))
			       )
                                | Env.PrimFun[f'] => checkPrimApp f' (*NOTE: This needs to change if any more operators like tensor, dot, colon product are added.*)
                                | Env.PrimFun ovldList => (
                                    case resolveOverload ((#1 cxt, span), f, tys, args, ovldList)
                                     of (e' as AST.E_Prim(f', tyArgs, _, _), rngTy) =>
(* NOTE: if/when we switch to matching type patterns (instead of unification),
 * we can use a "Self" type pattern to handle spatial queries.
 *)
                                          if Basis.isSpatialQueryOp f'
                                            then checkSpatialQuery (env, cxt, e', tyArgs, rngTy)
                                          else (e', rngTy)
				      | (e' as AST.E_Apply(_,_,_), rngTy) => (e', rngTy)
                                      | badResult => if Atom.same(Atom.atom "sphere", f)
						     then (print("SPHERE!\n");badResult)
						     else badResult
                                (* end case *))
                                | Env.UserFun f' => checkFunApp((#1 cxt, span), f') args tys
				| Env.OverloadUserFun(ovldList) => resolveOverload ((#1 cxt, span), f, tys, args, ovldList)
				    (* end case *))
					    (* end case*))
                    | _ => checkFieldApp (check (env, cxt, e))
                  (* end case *)
                end
            | PT.E_Subscript(e, indices) => let
                fun expectedTensor ty = err(cxt, [
                        S "expected tensor type for slicing, but found ", TY ty
                      ])
                fun chkIndex e = let
                      val eTy as (_, ty) = check(env, cxt, e)
                      in
                        if Unify.equalType(ty, Ty.T_Int)
                          then eTy
                          else err (cxt, [
                              S "expected type 'int' for index, but found ", TY ty
				   ])
                end
                (* for tensor/field/tuple slicing/indexing, the indices must be constant expressions *)
                fun chkConstIndex NONE = NONE
                  | chkConstIndex (SOME e) = (case chkIndex e
					       of (_, Ty.T_Error) => SOME bogusExp
						| (e', _) => (case CheckConst.eval (cxt, false, e')
										   (* FIXME: should check that index is in range for type! *)
							       of SOME cexp => SOME(ConstExpr.valueToExpr cexp)
								| NONE => SOME e' (* use e' to preserve variable uses *)
							     (* end case *))
					     (* end case *))				   
                val (e', ty) = check(env, cxt, e)
                in
                  case (TU.pruneHead ty, indices)
                   of (Ty.T_Error, _) => (
                        List.app (ignore o Option.map chkIndex) indices;
                        bogusExpTy)
                    | (ty1 as Ty.T_Sequence(elemTy, optDim), [SOME e2]) => let
                        val (e2', ty2) = chkIndex e2
                        val rator = if isSome optDim
                              then BV.subscript
                              else BV.dynSubscript
                        val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf rator)
                        in
                          if Unify.equalTypes(domTy, [ty1, ty2])
                            then let
                              val exp = AST.E_Prim(rator, tyArgs, [e', e2'], rngTy)
                              in
                                (exp, rngTy)
                              end
                            else raise Fail "unexpected unification failure"
                    end
		    | (ty1 as Ty.T_Tuple(tys), e2) => let
		     val size = List.length tys
		     val indicies = List.map chkConstIndex e2
		     val indexLength = List.length e2
		    in
		     if indexLength <> 1
		     then err (cxt, [S "Type error in slice operation\n",
				     S "Tuple types expect exactly one index."])
		     else
		      (case indicies
			of [SOME(AST.E_Lit(L.Int i))] => (let val idx = IntInf.toInt i
							      val _ = print("poop:" ^ (Int.toString idx))
							  in if 0 <= idx andalso idx < size
							     then (AST.E_Slice(e', indicies, List.nth(tys, idx)), List.nth(tys, idx))
							     else err (cxt, [S "Type error in slice operation\n",
									     S "Tuple type slices expect exactly one index in range.\n"])
							  end)
			 | _ => err (cxt, [S "Type error in slice operation\n",
					   S "Tuple type slices expect integers."])
		      (*end case*))
		    end
		      
                    | (ty as Ty.T_Sequence _, [NONE]) => expectedTensor ty
                    | (ty as Ty.T_Sequence _, _) => expectedTensor ty
                    | (ty, _) => let

                        val indices' = List.map chkConstIndex indices
                        val order = List.length indices'
                        (* val expectedTy = TU.mkTensorTy order*)
(* QUESTION: perhaps we should lift this case up above (i.e., one case for tensor and on for fields *)
                        val expectedTy = (case ty
                               of Ty.T_Field{diff, dim, shape=s as Ty.Shape(d2::dd2)} =>
                                    Ty.T_Field{diff=diff, dim=dim, shape=s}
                                | Ty.T_Tensor shape => TU.mkTensorTy order
                                | Ty.T_Field _ => raise Fail "unknown field type"
                                | ty => raise Fail("unexpected type for subscript: " ^ TU.toString ty)
                              (* end case *))
                        val resultTy = TU.slice(expectedTy, List.map Option.isSome indices')
                        in
                          if Unify.equalType(ty, expectedTy)
                            then (AST.E_Slice(e', indices', resultTy), resultTy)
                            else err (cxt, [
                                S "type error in slice operation\n",
                                S "  expected: ", S(Int.toString order), S "-order tensor\n",
                                S "  found:    ", TY ty
                              ])
                        end
                  (* end case *)
                end
            | PT.E_Select(e, field) => (case stripMark(#2 cxt, e)
                 of (_, PT.E_Var x) => (case E.findStrand (env, x)
                       of SOME _ => if E.inGlobalBlock env
                            then (case E.findSetFn (env, field)
                               of SOME setFn => let
                                    val (mvs, ty) = TU.instantiate (Var.typeOf setFn)
                                    val resTy = Ty.T_Sequence(Ty.T_Strand x, NONE)
                                    in
(* QUESTION: does it make sense to allow strand sets outside of reductions? *)
                                      E.recordProp (env, Properties.StrandSets);
                                      if Unify.equalType(ty, Ty.T_Fun([], resTy))
                                        then (AST.E_Prim(setFn, mvs, [], resTy), resTy)
                                        else raise Fail "impossible"
                                    end
                                | _ => err (cxt, [
                                      S "unknown strand-set specifier ", A field
                                    ])
                              (* end case *))
                            else err (cxt, [
                                S "illegal strand set specification in ",
                                S(E.scopeToString(E.currentScope env))
                              ])
                        | _ => (case Env.findTypeEnv(env, x)
				 of SOME(tyEnv) =>
				    let
				     val const = TypeEnv.findConstant(tyEnv, field)
				     (*Add all the other possible options here*)
				    in
				     (case const
				       of SOME(const') => (ConstExpr.valueToExpr const', ConstExpr.typeOfConst const')
					| NONE => (err (cxt, [
							S "Type Enviroment, ", A x, S ", does not have a type const", A field, S "."
						       ]);bogusExpTy)
				     (* end case *))
				    end
				   | _ => checkSelect (env, cxt, e, field)
			       (* end case*))
                      (* end case *))
                  | _ => checkSelect (env, cxt, e, field)
                (* end case *))
            | PT.E_Real e => (case checkAndPrune (env, cxt, e)
                 of (e', Ty.T_Int) =>
                      (AST.E_Prim(BV.i2r, [], [e'], Ty.realTy), Ty.realTy)
                  | (e', Ty.T_Error) => bogusExpTy
                  | (_, ty) => err(cxt, [
                        S "argument of 'real' must have type 'int', but found ",
                        TY ty
                      ])
                (* end case *))
            | PT.E_LoadSeq nrrd => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) = TU.instantiate(Var.typeOf(BV.fn_load_sequence))
                in
                  case chkStringConstExpr (env, cxt, nrrd)
                   of SOME nrrd => (AST.E_LoadNrrd(tyArgs, nrrd, rngTy), rngTy)
                    | NONE => (bogusExp, rngTy)
                  (* end case *)
                end
            | PT.E_LoadImage nrrd => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) = TU.instantiate(Var.typeOf(BV.fn_load_image))
                in
                  case chkStringConstExpr (env, cxt, nrrd)
                   of SOME nrrd => (AST.E_LoadNrrd(tyArgs, nrrd, rngTy), rngTy)
                    | NONE => (bogusExp, rngTy)
                  (* end case *)
                end
            | PT.E_Var x => (case E.findVar (env, x)
                 of SOME x' => (AST.E_Var(useVar(cxt, x')), Var.monoTypeOf x')
                  | NONE => (case E.findKernel (env, x)
                       of SOME h => (AST.E_Kernel h, TypeOf.kernel h)
                        | NONE => err(cxt, [S "undeclared variable ", A x])
                      (* end case *))
                (* end case *))
            | PT.E_Kernel(kern, dim) => let
                val k' = Int.fromLarge dim handle Overflow => 1073741823
                fun mkExp (e, k, ty) = if (k = k')
                      then (e, ty)
                      else let
                        val ty' = Ty.T_Kernel(Ty.DiffConst(SOME k'))
                        in
                          (AST.E_Coerce{srcTy = ty, dstTy = ty', e = e}, ty')
                        end
                in
                  case E.findVar (env, kern)
                   of SOME h => (case Var.monoTypeOf h
                         of ty as Ty.T_Kernel(Ty.DiffConst(SOME k)) =>
                              mkExp (AST.E_Var(useVar(cxt, h)), k, ty)
                          | _ => err(cxt, [
                                S "expected kernel, but found ", S(Var.kindToString h)
                              ])
                        (* end case *))
                    | NONE => (case E.findKernel (env, kern)
                         of SOME h => let
                              val k = Kernel.continuity h
                              in
                                mkExp (AST.E_Kernel h, k, TypeOf.kernel h)
                              end
                          | NONE => err(cxt, [S "unknown kernel ", A kern])
                        (* end case *))
                  (* end case *)
                end
            | PT.E_Lit lit => checkLit lit
            | PT.E_Id d => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) =
                      TU.instantiate(Var.typeOf(BV.identity))
                in
                  if Unify.equalType(Ty.T_Tensor(checkShape(env, cxt, [d, d])), rngTy)
                    then (AST.E_Prim(BV.identity, tyArgs, [], rngTy), rngTy)
                    else raise Fail "impossible"
            end
	    (*identityFld*)
	    | PT.E_IdFld d => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) =
                    TU.instantiate(Var.typeOf(BV.identityFld))
				  (*checkShape(env, cxt, [d, d])*)
		val dimCheck = (case checkDim(env, cxt, d)
				 of SOME(i) => Ty.DimConst (IntInf.toInt i)
				  | NONE => Ty.DimConst (~1) (*error later*)
			       (*end case*))
		val shapeCheck = checkShape(env, cxt, [d])
                in
                 if Unify.equalType(Ty.T_Field{diff=Ty.DiffConst NONE, dim=dimCheck, shape=shapeCheck}, rngTy)
                 then (AST.E_Prim(BV.identityFld, tyArgs, [], rngTy), rngTy)
                 else raise Fail "impossible"
            end
            | PT.E_Zero dd => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) =
                      TU.instantiate(Var.typeOf(BV.zero))
                in
                  if Unify.equalType(Ty.T_Tensor(checkShape(env, cxt, dd)), rngTy)
                    then (AST.E_Prim(BV.zero, tyArgs, [], rngTy), rngTy)
                    else raise Fail "impossible"
                end
            | PT.E_NaN dd => let
                val (tyArgs, Ty.T_Fun(_, rngTy)) =
                      TU.instantiate(Var.typeOf(BV.nan))
                in
                  if Unify.equalType(Ty.T_Tensor(checkShape(env, cxt, dd)), rngTy)
                    then (AST.E_Prim(BV.nan, tyArgs, [], rngTy), rngTy)
                    else raise Fail "impossible"
                end
            | PT.E_Sequence exps => (case checkList (env, cxt, exps)
                 of ([], _) => let
(* FIXME: metavar should have kind for concrete types here! *)
                      val ty = Ty.T_Sequence(Ty.T_Var(MetaVar.newTyVar()), SOME(Ty.DimConst 0))
                      in
                        (AST.E_Seq([], ty), ty)
                      end
                  | (args, tys) => (case Util.coerceTypes(List.map TU.pruneHead tys)
                       of SOME ty => if TU.isValueType ty
                            then let
                              fun doExp eTy = valOf(Util.coerceType (ty, eTy))
                              val resTy = Ty.T_Sequence(ty, SOME(Ty.DimConst(List.length args)))
                              val args = ListPair.map doExp (args, tys) handle exn => raise exn
                              in
                                (AST.E_Seq(args, resTy), resTy) handle exn => raise exn
				     end
                            else err(cxt, [S "sequence expression of non-value argument type"])
                        | NONE => err(cxt, [S "arguments of sequence expression must have same type"])
                      (* end case *))
                (* end case *))
            | PT.E_SeqComp comp => chkComprehension (env, cxt, comp)
	    | PT.E_Tuple args => let
	     val (args, tys) = checkList (env, cxt, args) (*Length args >= 2 b/c parser*)
	     val badType = List.find (Bool.not o TU.isValueOrStrandType) tys
	    in
	     (case badType
	       of NONE => (AST.E_Tuple(args, tys), Ty.T_Tuple(tys))
		| SOME(ty) => (print(TU.toString ty); err(cxt, [S "invalid element type for tuple: ", TY ty]))
	     (*end case*))
	    end
            | PT.E_Cons args => let
              (* Note that we are guaranteed that args is non-empty *)
                val (args, tys) = checkList (env, cxt, args)
              (* extract the first non-error type in tys *)
                val ty = (case List.find (fn Ty.T_Error => false | _ => true) tys
                       of NONE => Ty.T_Error
                        | SOME ty => ty
                      (* end case *))
              (* process the arguments checking that they all have the expected type *)
                fun chkArgs (ty, shape) = let
                      val Ty.Shape dd = TU.pruneShape shape (* NOTE: this may fail if we allow user polymorphism *)
                      val resTy = Ty.T_Tensor(Ty.Shape(Ty.DimConst(List.length args) :: dd))
                      fun chkArgs (arg::args, argTy::tys, args') = (
                            case Util.coerceType(ty, (arg, argTy))
                             of SOME arg' => chkArgs (args, tys, arg'::args')
                              | NONE => (
                                  TypeError.error(cxt, [
                                      S "arguments of tensor construction must have same type"
                                      (* FIXME: add types to error message *)
                                    ]);
                                  chkArgs (args, tys, bogusExp::args'))
                            (* end case *))
                        | chkArgs (_, _, args') = (AST.E_Tensor(List.rev args', resTy), resTy)
                      in
                        chkArgs (args, tys, [])
                      end
                fun chkArgsF (ty, diff, dim, shape) = let
                      val Ty.Shape dd = TU.pruneShape shape
                      val resTy = Ty.T_Field{diff=diff ,dim=dim, shape=Ty.Shape(Ty.DimConst(List.length args) :: dd)}
                      fun chkArgsF (arg::args, argTy::tys, args') = (case Util.coerceType(ty, (arg, argTy))
                             of SOME arg' => chkArgsF (args, tys, arg'::args')
                              | NONE => (
                                  TypeError.error(cxt, [
                                      S "arguments of tensor construction must have same type"
                                      (* FIXME: add types to error message *)
                                    ]);
                                  chkArgsF (args, tys, bogusExp::args'))
                            (* end case *))
                        | chkArgsF (_, _, args') = (AST.E_Field(List.rev args', resTy), resTy)
                      in
                        chkArgsF (args, tys, [])
                      end
                in
                  case TU.pruneHead ty
                   of Ty.T_Int => chkArgs(Ty.realTy, Ty.Shape[]) (* coerce integers to reals *)
                    | ty as Ty.T_Tensor shape => chkArgs(ty, shape)
                    | ty as Ty.T_Field{diff, dim, shape} => chkArgsF(ty, diff, dim, shape)
                    | _ => err(cxt, [S "Invalid argument type for tensor construction"])
                  (* end case *)
                end
        (* end case *))

  (* typecheck and the prune the result *)
    and checkAndPrune (env, cxt, e) = let
          val (e, ty) = check (env, cxt, e)
          in
            (e, TU.prune ty)
          end

  (* check a conditional operator (e.g., || or &&) *)
    and checkCondOp (env, cxt, e1, rator, e2, mk) = (
          case (check(env, cxt, e1), check(env, cxt, e2))
           of ((e1', Ty.T_Bool), (e2', Ty.T_Bool)) => (mk(e1', e2'), Ty.T_Bool)
            | ((_, Ty.T_Bool), (_, ty2)) =>
                err (cxt, [S "expected type 'bool' on rhs of '", S rator, S "', but found ", TY ty2])
            | ((_, ty1), (_, Ty.T_Bool)) =>
                err (cxt, [S "expected type 'bool' on lhs of '", S rator, S "', but found ", TY ty1])
            | ((_, ty1), (_, ty2)) => err (cxt, [
                  S "arguments of '", S rator, S "' must have type 'bool', but found ",
                  TY ty1, S " and ", TY ty2
                ])
          (* end case *))

  (* check a field select that is _not_ a strand-set *)
    and checkSelect (env, cxt, e, field) = (case checkAndPrune (env, cxt, e)
           of (e', Ty.T_Strand strand) => (case Env.findStrand(env, strand)
                 of SOME sEnv => (case StrandEnv.findStateVar(sEnv, field)
                       of SOME x' => let
                            val ty = Var.monoTypeOf x'
                            in
                              (AST.E_Select(e', useVar(cxt, x')), ty)
                            end
                        | NONE => err(cxt, [
                              S "strand ", A strand,
                              S " does not have state variable ", A field
                            ])
                      (* end case *))
                  | NONE => err(cxt, [S "unknown strand ", A strand])
					  (* end case *))
	    | (e' , Ty.T_Fem(data, _)) =>
	      let
	       val tyName = SOME(FT.envNameOf data ) handle exn => NONE
	       val tyEnv = (Option.mapPartial (fn name => case Env.findTypeEnv(env, name)
								    of SOME(s) => SOME(s)
								     | NONE =>  (err (cxt, [S "There is no type named", A(name)]); NONE) ) tyName)
	       val method = Option.mapPartial (fn env => (case TypeEnv.findHiddenVar(env, field)
								 of SOME(s) => SOME(s)
								  | NONE => ( (err (cxt, [S"The type named", A (TypeEnv.findName env),
										     S " does not have a member named",
										     A(field)])); NONE))) tyEnv
			    
	      in
	       (case method (*I think that any aditional typechecking here is redundant*)
		 of SOME(f, Ty.T_Fun([v], a)) =>
		    let
		     val fieldExpr = f [e']
		     val _ = print("found\n")		       
		    in
		     (fieldExpr, a)
		    end
		    
		  | NONE => bogusExpTy
	       (* end case *))
	      end
	    | (e', namedTy as Ty.T_Named(name, ty')) =>
	      let
	       (* There in one case here: 1 argument functions s*)
	       (* NOTE TO SELF: THIS checks a named type value*)
	       (* I'm only implementing constants for now.*)
	       val tyEnv = Env.findTypeEnv(env, name)
	       val const = Option.mapPartial (fn x => TypeEnv.findConstant(x, field)) tyEnv
	       val method =  Option.mapPartial (fn x => TypeEnv.findMethod(x, field)) tyEnv
	       val func = Option.mapPartial (fn x => TypeEnv.findHiddenVar(x, field)) tyEnv
	       (*Var.monoTypeOf f*)
	       fun checkFunApp(f, Ty.T_Fun([dom], rng)) = (case Unify.matchArgs([dom], [e'], [namedTy])
							    of SOME(args) => (AST.E_Apply(useVar(cxt, f), [e'], rng), rng)
							     | NONE => raise Fail "IS THIS ACTUALlY POSSIBLE???"
							  (*end case*))
		 | checkFunApp _ = ((err (cxt, [S"Application of function ",
						A (field), S ", a member of the type ",
						A(name), S "requires arguments"]));bogusExpTy)

	       fun checkFunApp'(f, Ty.T_Fun([dom], rng)) =
		   (case Unify.matchArgs([dom], [e'], [namedTy])
		     of SOME(args) => (f [e'], rng)
		      | NONE => ((err (cxt, [S"Application of function ", A (field), S ", a member of the type ", A(name),
					     S "requires argument of type ", TY(dom), S" but got an argument of type",
					     TY(namedTy), S"!"]));bogusExpTy)
		   (*end case*))
		 | checkFunApp' _ = ((err (cxt, [S"Application of function ", A (field),
						 S ", a member of the type ", A(name),
						 S "requires arguments"]));bogusExpTy)
				      
					    
	      in
	       (case (tyEnv, const, method, func)
	    	 of (_, SOME(const'), _, _) =>
		    (ConstExpr.valueToExpr const', ConstExpr.typeOfConst const')
		  | (_, NONE, SOME(method), _) => checkFunApp(method, Var.monoTypeOf method)
		  | (_, NONE, NONE, SOME(func', fTy)) => checkFunApp'(func', fTy)
	    	  | (SOME(_), NONE, NONE, NONE) => (err(cxt, [S "In selection of ", A field, S " from  the type env, ", A(name), S ", the field could not be located!"]); bogusExpTy)
		  | _ => (err (cxt, [S "Trying to select a field, ", A field, S ", from type, ", A(name), S ", but the type does not exist!"]);bogusExpTy)
	       (* end case*))
	      end
            | (_, Ty.T_Error) => bogusExpTy
            | (_, ty) => err (cxt, [
                  S "expected strand type or named type, but found ", TY ty,
                  S " in selection of ", A field
                ])
          (* end case *))

    and chkComprehension (env, cxt, PT.COMP_Mark m) =
          chkComprehension(E.withEnvAndContext(env, cxt, m))
      | chkComprehension (env, cxt, PT.COMP_Comprehension(e, [iter])) = let
          val (iter', env') = checkIter (E.blockScope env, cxt, iter)
          val (e', ty) = check (env', cxt, e)
          val resTy = Ty.T_Sequence(ty, NONE)
          in
            case iter'
             of (x, AST.E_Prim(f, _, [], _)) => if Basis.isStrandSet f
                  then (
                    Env.recordProp (env, Properties.GlobalReduce);
                    if not(Env.inGlobalBlock env)
                      then err (cxt, [
                          S "strand comprehension outside of global initialization or update"
                        ])
                    else if Env.inLoop env
                      then err (cxt, [
                          S "strand comprehension inside loop"
                        ])
                      else (AST.E_ParallelMap(e', x, f, resTy), resTy))
                  else (AST.E_Comprehension(e', iter', resTy), resTy)
              | _ => (AST.E_Comprehension(e', iter', resTy), resTy)
            (* end case *)
          end
      | chkComprehension _ = raise Fail "impossible"

    and checkIter (env, cxt, PT.I_Mark m) = checkIter (E.withEnvAndContext (env, cxt, m))
      | checkIter (env, cxt, PT.I_Iterator({span, tree=x}, e)) = (
          case checkAndPrune (env, cxt, e)
           of (e', ty as Ty.T_Sequence(elemTy, _)) => let
                val x' = Var.new(x, span, Var.LocalVar, elemTy)
                in
                  ((x', e'), E.insertLocal(env, cxt, x, x'))
                end
            | (e', ty) => let
                val x' = Var.new(x, span, Var.IterVar, Ty.T_Error)
                in
                  if TU.isErrorType ty
                    then ()
                    else TypeError.error (cxt, [
                        S "expected sequence type in iteration, but found '", TY ty, S "'"
                      ]);
                  ((x', bogusExp), E.insertLocal(env, cxt, x, x'))
                end
          (* end case *))

  (* typecheck a list of expressions returning a list of AST expressions and a list
   * of the types of the expressions.
   *)
    and checkList (env, cxt, exprs) = let
          fun chk (e, (es, tys)) = let
                val (e, ty) = checkAndPrune (env, cxt, e)
                in
                  (e::es, ty::tys)
                end
          in
            List.foldr chk ([], []) exprs
          end

  (* check a string that is specified as a constant expression *)
    and chkStringConstExpr (env, cxt, PT.E_Mark m) =
          chkStringConstExpr (E.withEnvAndContext (env, cxt, m))
      | chkStringConstExpr (env, cxt, e) = (case checkAndPrune (env, cxt, e)
           of (e', Ty.T_String) => (case CheckConst.eval (cxt, false, e')
                 of SOME(ConstExpr.String s) => SOME s
                  | SOME(ConstExpr.Expr e) => raise Fail "FIXME"
                  | NONE => NONE
                  | _ => raise Fail "impossible: wrong type for constant expr"
                (* end case *))
            | (_, Ty.T_Error) => NONE
            | (_, ty) => (
                TypeError.error (cxt, [
                    S "expected constant expression of type 'string', but found '",
                    TY ty, S "'"
                  ]);
                NONE)
          (* end case *))

  (* check a dimension that is given by a constant expression *)
    and checkDim (env, cxt, dim) = (case checkAndPrune (env, cxt, dim)
           of (e', Ty.T_Int) => (case CheckConst.eval (cxt, false, e')
                 of SOME(ConstExpr.Int d) => SOME d
                  | SOME(ConstExpr.Expr e) => (
                      TypeError.error (cxt, [S "unable to evaluate constant dimension expression"]);
                      NONE)
                  | NONE => NONE
                  | _ => raise Fail "impossible: wrong type for constant expr"
                (* end case *))
            | (_, Ty.T_Error) => NONE
            | (_, ty) => (
                TypeError.error (cxt, [
                    S "expected constant expression of type 'int', but found '",
                    TY ty, S "'"
                  ]);
                NONE)
          (* end case *))

  (* check a tensor shape, where the dimensions are given by constant expressions *)
    and checkShape (env, cxt, shape) = let
          fun checkDim' e = (case checkDim (env, cxt, e)
                 of SOME d => (
                      if (d <= 1)
                        then TypeError.error (cxt, [
                            S "invalid tensor-shape dimension; must be > 1, but found ",
                            S (IntLit.toString d)
                          ])
                        else ();
                      Ty.DimConst(IntInf.toInt d))
                  | NONE => Ty.DimConst ~1 
			    (* end case *))
	  fun normalizeShape(x::xs) =
	      (case checkAndPrune (env, cxt, x)
		of (e', Ty.T_Int) => x :: normalizeShape(xs)
		 | (e', Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst(n)))) =>
		   let
		    val newDims = List.tabulate(n, fn y=> PT.E_Subscript(x, ([SOME(PT.E_Lit(Literal.intLit y))])))
		   in
		    List.@(newDims,normalizeShape(xs))
		   end
		 | _ => x::normalizeShape(xs)
	      (* end case*))
	    | normalizeShape [] = []

          in
            Ty.Shape(List.map checkDim' (normalizeShape(shape)))
          end

  end
