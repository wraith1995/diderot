(* check-overload.sml
 *
 * A part of the type checker for figuring out if overloads can be resolved, conflict, or if there is some internal error in the code that prevents checking.
 * Includes various utilities that can be used to deal with special cases around products.

 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure CheckOverload : sig

	   datatype overload
	     = Overload of (AST.expr * Types.ty)
	     | NoOverload
	     | Error of Env.context * TypeError.token list

	   (* Finds the first possible overload *)
	   val chkOverload : Env.context * Atom.atom * Types.ty list * AST.expr list * Var.t list -> overload

	   (* Determines if a binary operator needs to be checked for overloads outside of the basis *)
	   val checkSpecialProduct : Types.ty list * Atom.atom -> Types.shape option option

	   (* Determines if a call matches a specific binary operator*)
	   val chkInnerProduct : Env.context * AST.expr * Types.ty * AST.expr * Types.ty -> overload
	   val chkOuterProduct : Env.context * AST.expr * Types.ty * AST.expr * Types.ty -> overload
	   val chkColonProduct : Env.context * AST.expr * Types.ty * AST.expr * Types.ty -> overload

	   (* Converts an application of a primitive to that of a function if needed.*)
	   val chkPrim : AST.expr * Types.ty ->  AST.expr * Types.ty
	  end = struct 

    structure PT = ParseTree
    structure L = Literal
    structure E = Env
    structure Ty = Types
    structure BV = BasisVars
    structure TU = TypeUtil

    datatype overload
      = Overload of (AST.expr * Types.ty)
      | NoOverload
      | Error of (Error.err_stream * Error.span) * TypeError.token list


  (* an expression to return when there is a type error *)
    val bogusExp = AST.E_Lit(L.Int 0)
    val bogusExpTy = (bogusExp, Ty.T_Error)
    val (S, A, TYS) = (TypeError.S, TypeError.A, TypeError.TYS)

    fun err arg = (TypeError.error arg; bogusExpTy)

    fun chkPrim e =
	(case e
	  of (AST.E_Prim(f', _, expList, fnTy), rngTy) =>
	     if  Var.kindOf f' = Var.BasisVar
	     then e
	     else if Var.kindOf f' = Var.FunVar
	     then let val span = Var.locationOf f' in (AST.E_Apply((f', span), expList, fnTy), rngTy) end
	     else raise Fail "impossible"
	   | _ => raise Fail "impossible"
	(* end case *))
		    
		    
    fun overloadResult a = Overload (chkPrim a)
    fun internalError a = Error a

(* resolve overloading: we use a simple scheme that selects the first operator in the
   * list that matches the argument types.
   *)
    fun chkOverload (_, rator, _, _, []) = NoOverload
      | chkOverload (cxt, rator, argTys, args, candidates) = let
(* FIXME: we could be more efficient by just checking for coercion matchs the first pass
 * and remembering those that are not pure EQ matches.
 *)
        (* build the result *)
          fun done (rator, tyArgs, args, rngTy) = if Var.same(rator, BV.pow_si)
                then let (* check that the second argument is a constant expression *)
                  val [e1, e2] = args
                  in
                    case CheckConst.eval (cxt, false, e2)
                     of SOME e2' =>
                          overloadResult (AST.E_Prim(rator, tyArgs, [e1, ConstExpr.valueToExpr e2'], rngTy), rngTy)
                      | NONE => internalError (cxt, [
                            S "constant-integer exponent is required when lhs of '^' is a field"
                          ])
                  end
                else overloadResult (AST.E_Prim(rator, tyArgs, args, rngTy), rngTy)
        (* try to match candidates while allowing type coercions *)
          fun tryMatchCandidates [] = internalError (cxt, [
                  S "unable to resolve overloaded operator ", A rator, S "\n",
                  S "  argument type is: ", TYS argTys, S "\n"
                ])
            | tryMatchCandidates (x::xs) = let
                val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf x)
                in
                  case Unify.tryMatchArgs (domTy, args, argTys)
                   of SOME args' => done(x, tyArgs, args', rngTy)
                    | NONE => tryMatchCandidates xs
                  (* end case *)
                end
        (* try to match candidates without type coercions *)
          fun tryCandidates [] = tryMatchCandidates candidates
            | tryCandidates (x::xs) = let
                val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf x)
                in
                  if Unify.tryEqualTypes(domTy, argTys)
                    then done(x, tyArgs, args, rngTy)
                    else tryCandidates xs
                end
          in
           tryCandidates candidates
      end

							       

    fun mkChkProductParams params shapeChk =
	if List.length params = 2
	then (case (List.map TU.prune params)
	       of [Ty.T_Tensor s1, Ty.T_Tensor s2] => shapeChk(s1, s2)
		| _ => NONE)
	else NONE							       


    fun chkInnerProductShape (Ty.Shape(dd1 as _::_), Ty.Shape(d2::dd2)) = let
     val (dd1, d1) = let
      fun splitLast (prefix, [d]) = (List.rev prefix, d)
        | splitLast (prefix, d::dd) = splitLast (d::prefix, dd)
        | splitLast (_, []) = raise Fail "impossible"
     in
      splitLast ([], dd1)
     end
    in
     if Unify.equalDim(d1, d2)
     then SOME(Ty.Shape(dd1@dd2))
     else NONE
    end
      | chkInnerProductShape _ = NONE

    fun chkOuterProductShape (Ty.Shape dd1, Ty.Shape dd2) = SOME(Ty.Shape(dd1@dd2))
      | chkOuterProductShape _ = NONE
				   
    fun chkColonProductShape (Ty.Shape(dd1 as _::_::_), Ty.Shape(d21::d22::dd2)) = let
     val (dd1, d11, d12) = let
      fun splitLast2 (prefix, [d1, d2]) = (List.rev prefix, d1, d2)
        | splitLast2 (prefix, d::dd) = splitLast2 (d::prefix, dd)
        | splitLast2 (_, []) = raise Fail "impossible"
     in
      splitLast2 ([], dd1)
     end
    in
     if Unify.equalDim(d11, d21) andalso Unify.equalDim(d12, d22)
     then SOME(Ty.Shape(dd1@dd2))
     else NONE
    end
      | chkColonProductShape _ = NONE

    fun chkInnerProductParams params = mkChkProductParams params chkInnerProductShape
    fun chkOuterProductParams params = mkChkProductParams params chkOuterProductShape
    fun chkColonProductParams params = mkChkProductParams params chkColonProductShape

    fun checkSpecialProduct(params, special) =
	if Atom.same(special, BasisNames.op_dot)
	then SOME(chkInnerProductParams params) 
	else if Atom.same(special, BasisNames.op_outer)
	then SOME(chkOuterProductParams params)
	else if Atom.same(special, BasisNames.op_colon)
	then SOME(chkColonProductParams params)
	else NONE


  (* type check a dot product, which has the constraint:
   *     ALL[sigma1, d1, sigma2] . tensor[sigma1, d1] * tensor[d1, sigma2] -> tensor[sigma1, sigma2]
   * and similarly for fields.
   *)
    fun chkInnerProduct (cxt, e1, ty1, e2, ty2) = let
     (* check the shape of the two arguments to verify that the inner constraint matches *)
     fun error () = internalError (cxt, [
			 S "type error for arguments of binary operator '•'\n",
			 S "  found: ", TYS[ty1, ty2], S "\n"
			])
    in
     case (TU.prune ty1, TU.prune ty2)
            (* tensor * tensor inner product *)
      of (Ty.T_Tensor s1, Ty.T_Tensor s2) => (case chkInnerProductShape(s1, s2)
					       of SOME shp => let
						val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_inner_tt)
						val resTy = Ty.T_Tensor shp
					       in
						if Unify.equalTypes(domTy, [ty1, ty2]) andalso Unify.equalType(rngTy, resTy)
						then overloadResult (AST.E_Prim(BV.op_inner_tt, tyArgs, [e1, e2], rngTy), rngTy)
						else error()
					       end
						| NONE => NoOverload
					     (* end case *))
       (* tensor * field inner product *)
       | (Ty.T_Tensor s1, Ty.T_Field{diff, dim, shape=s2}) => (case chkInnerProductShape(s1, s2)
								of SOME shp => let
								 val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_inner_tf)
								 val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
								in
								 if Unify.equalTypes(domTy, [ty1, ty2])
								    andalso Unify.equalType(rngTy, resTy)
								 then overloadResult (AST.E_Prim(BV.op_inner_tf, tyArgs, [e1, e2], rngTy), rngTy)
								 else error()
								end
								 | NONE => error()
							      (* end case *))
       (* field * tensor inner product *)
       | (Ty.T_Field{diff, dim, shape=s1}, Ty.T_Tensor s2) => (case chkInnerProductShape(s1, s2)
								of SOME shp => let
								 val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_inner_ft)
								 val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
								in
								 if Unify.equalTypes(domTy, [ty1, ty2])
								    andalso Unify.equalType(rngTy, resTy)
								 then overloadResult (AST.E_Prim(BV.op_inner_ft, tyArgs, [e1, e2], rngTy), rngTy)
								 else error()
								end
								 | NONE => error()
							      (* end case *))
       (* field * field inner product *)
       | (Ty.T_Field{diff=k1, dim=dim1, shape=s1}, Ty.T_Field{diff=k2, dim=dim2, shape=s2}) => (
        case chkInnerProductShape(s1, s2)
         of SOME shp => let
          val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_inner_ff)
          val resTy = Ty.T_Field{diff=k1, dim=dim1, shape=shp}
         in
	  (* FIXME: the resulting differentiation should be the minimum of k1 and k2 *)
          if Unify.equalDim(dim1, dim2)
             andalso Unify.equalTypes(domTy, [ty1, ty2])
             andalso Unify.equalType(rngTy, resTy)
          then overloadResult (AST.E_Prim(BV.op_inner_ff, tyArgs, [e1, e2], rngTy), rngTy)
          else error()
         end
          | NONE => error()
       (* end case *))
       | (ty1, ty2) => NoOverload
			    (* end case *)
    end

  (* type check an outer product, which has the constraint:
   *     ALL[sigma1, sigma2] . tensor[sigma1] * tensor[sigma2] -> tensor[sigma1, sigma2]
   * and similarly for fields.
   *)						    
    fun chkOuterProduct (cxt, e1, ty1, e2, ty2) = let
     (*QUESTION: CAN THIS ACTUALLY BE CALLED?*)
     fun shapeError () = internalError (cxt, [
			      S "unable to determine result shape of outer product\n",
			      S "  found: ", TYS[ty1, ty2], S "\n"
			     ])
     fun error () = internalError (cxt, [
			 S "type error for arguments of binary operator \"⊗\"\n",
			 S "  found: ", TYS[ty1, ty2], S "\n"
			])
    in
     case (TU.prune ty1, TU.prune ty2)
            (* tensor * tensor outer product *)
      of (Ty.T_Tensor s1, Ty.T_Tensor s2) => (case chkOuterProductShape(s1, s2)
					       of SOME shp => let
						val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_outer_tt)
						val resTy = Ty.T_Tensor shp
					       in
						if Unify.equalTypes(domTy, [ty1, ty2])
						   andalso Unify.equalType(rngTy, resTy)
						then overloadResult (AST.E_Prim(BV.op_outer_tt, tyArgs, [e1, e2], rngTy), rngTy)
						else error()
					       end
						| NONE => shapeError()
					     (* end case *))
       (* field * tensor outer product *)
       | (Ty.T_Field{diff, dim, shape=s1}, Ty.T_Tensor s2) => (case chkOuterProductShape(s1, s2)
								of SOME shp => let
								 val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_outer_ft)
								 val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
								in
								 if Unify.equalTypes(domTy, [ty1, ty2]) andalso Unify.equalType(rngTy, resTy)
								 then overloadResult (AST.E_Prim(BV.op_outer_ft, tyArgs, [e1, e2], rngTy), rngTy)
								 else error()
								end
								 | NONE => shapeError()
							      (* end case *))
       (* tensor * field outer product *)
       | (Ty.T_Tensor s1, Ty.T_Field{diff=diff, dim=dim, shape=s2}) => (case chkOuterProductShape(s1, s2)
									 of SOME shp => let
									  val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_outer_tf)
									  val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
									 in
									  if Unify.equalTypes(domTy, [ty1, ty2]) andalso Unify.equalType(rngTy, resTy)
									  then overloadResult (AST.E_Prim(BV.op_outer_tf, tyArgs, [e1, e2], rngTy), rngTy)
									  else error()
									 end
									  | NONE => shapeError()
								       (* end case *))
       (* field * field outer product *)
       | (Ty.T_Field{diff=k1, dim=dim1, shape=s1}, Ty.T_Field{diff=k2, dim=dim2, shape=s2}) => (
        case chkOuterProductShape(s1, s2)
         of SOME shp => let
          val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_outer_ff)
          val resTy = Ty.T_Field{diff=k1, dim=dim1, shape=shp}
         in
	  (* FIXME: the resulting differentiation should be the minimum of k1 and k2 *)
          if Unify.equalDim(dim1, dim2)
             andalso Unify.equalTypes(domTy, [ty1, ty2])
             andalso Unify.equalType(rngTy, resTy)
          then overloadResult (AST.E_Prim(BV.op_outer_ff, tyArgs, [e1, e2], rngTy), rngTy)
          else error()
         end
          | NONE => shapeError()
       (* end case *))
       | _ => NoOverload
		   (* end case *)
    end

  (* type check a colon product, which has the constraint:
   *     ALL[sigma1, d1, d2, sigma2] . tensor[sigma1, d1, d2] * tensor[d2, d1, sigma2] -> tensor[sigma1, sigma2]
   * and similarly for fields.
   *)
    fun chkColonProduct (cxt, e1, ty1, e2, ty2) = let
          fun error () = internalError (cxt, [
                  S "type error for arguments of binary operator \":\"\n",
                  S "  found: ", TYS[ty1, ty2], S "\n"
                ])
          in
            case (TU.prune ty1, TU.prune ty2)
            (* tensor * tensor colon product *)
             of (Ty.T_Tensor s1, Ty.T_Tensor s2) => (case chkColonProductShape(s1, s2)
                   of SOME shp => let
                        val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_colon_tt)
                        val resTy = Ty.T_Tensor shp
                        in
                          if Unify.equalTypes(domTy, [ty1, ty2])
                          andalso Unify.equalType(rngTy, resTy)
                            then overloadResult (AST.E_Prim(BV.op_colon_tt, tyArgs, [e1, e2], rngTy), rngTy)
                            else error()
                        end
                    | NONE => NoOverload
                  (* end case *))
            (* field * tensor colon product *)
              | (Ty.T_Field{diff, dim, shape=s1}, Ty.T_Tensor s2) => (case chkColonProductShape(s1, s2)
                   of SOME shp => let
                        val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_colon_ft)
                        val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
                        in
                          if Unify.equalTypes(domTy, [ty1, ty2]) andalso Unify.equalType(rngTy, resTy)
                            then overloadResult (AST.E_Prim(BV.op_colon_ft, tyArgs, [e1, e2], rngTy), rngTy)
                            else error()
                        end
                    | NONE => error()
                  (* end case *))
            (* tensor * field colon product *)
              | (Ty.T_Tensor s1, Ty.T_Field{diff=diff, dim=dim, shape=s2}) => (case chkColonProductShape(s1, s2)
                   of SOME shp => let
                        val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_colon_tf)
                        val resTy = Ty.T_Field{diff=diff, dim=dim, shape=shp}
                        in
                          if Unify.equalTypes(domTy, [ty1, ty2]) andalso Unify.equalType(rngTy, resTy)
                            then overloadResult (AST.E_Prim(BV.op_colon_tf, tyArgs, [e1, e2], rngTy), rngTy)
                            else error()
                        end
                    | NONE => error()
                  (* end case *))
            (* field * field colon product *)
              | (Ty.T_Field{diff=k1, dim=dim1, shape=s1}, Ty.T_Field{diff=k2, dim=dim2, shape=s2}) => (
                  case chkColonProductShape(s1, s2)
                   of SOME shp => let
                        val (tyArgs, Ty.T_Fun(domTy, rngTy)) = TU.instantiate(Var.typeOf BV.op_colon_ff)
                        val resTy = Ty.T_Field{diff=k1, dim=dim1, shape=shp}
                        in
(* FIXME: the resulting differentiation should be the minimum of k1 and k2 *)
                          if Unify.equalDim(dim1, dim2)
                          andalso Unify.equalTypes(domTy, [ty1, ty2])
                          andalso Unify.equalType(rngTy, resTy)
                            then overloadResult (AST.E_Prim(BV.op_colon_ff, tyArgs, [e1, e2], rngTy), rngTy)
                            else error()
                        end
                    | NONE => error()
                  (* end case *))
              | (ty1, ty2) => NoOverload
            (* end case *)
          end

						    
						    
						    
end
