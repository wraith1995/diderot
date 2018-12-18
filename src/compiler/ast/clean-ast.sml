(* clean-ast.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *
 * Removes all named types from the parse tree and global env. Replaces them with their definitions. This is valid by the construction of the new type feature. 
 *)

structure CleanAst : sig
	   (*Cleans the ast *)
	   val cleanAST : AST.program -> AST.program
	   (*Cleans a type.*)
	   val cleanTy : Types.ty -> Types.ty
	   (*Cleans an AST.var*)
	   val varCleanType : Var.t -> Var.t
end = struct

structure Ty = Types



fun cleanTy ty =
    (case ty
      of Ty.T_Named(name, def) => def
       | Ty.T_Fun(tys, ty) => Ty.T_Fun(List.map cleanTy tys, cleanTy ty)
       | Ty.T_Sequence(ty, dim) => Ty.T_Sequence(cleanTy ty, dim)
       | Ty.T_Var(meta) => Ty.T_Var(cleanMetaTy meta)
       | _ => ty)
and cleanMetaTy (res as Ty.TV {id,bind}) =
    (case !bind
      of NONE => res
       | SOME(ty) => let val _ = bind := SOME (cleanTy ty) in  res end)

val varCleanType = Var.cleanType cleanTy

fun cleanAST (AST.Program({props, const_dcls,
			  input_dcls, globals,
			  globInit, strand,
			  create, start, update})) =
	 AST.Program({props=props, const_dcls = List.map cleanVarDcl const_dcls,
		      input_dcls= List.map (fn (x,y) => (cleanVarDcl x,y)) input_dcls,
		      globals=List.map cleanGlobalDcl globals,
		      globInit= Option.map cleanStm globInit,
		      strand=strand,
		      create=cleanCreate create,
		      start= Option.map cleanStm start,
		      update=Option.map cleanStm update})

	  
    and cleanStrand Strand ({name, params, spatialDim, state,
    			     stateInit, startM, updateM, stabilizeM}) =
    	{
    	  name = name, params = List.map varCleanType params,
	  spaitalDim = spatialDim,
	  state = List.map cleanVarDcl state,
	  stateInit = Option.map cleanStm stateInit,
	  startM =  Option.map cleanStm startM,
	  updateM =  cleanStm updateM,
	  stabilizeM =  Option.map cleanStm stabilizeM
    	}

    and cleanCreate (Create.Create({dim, code})) = Create.Create{dim=dim, code = cleanStm code}
    and cleanVarDcl (v,e) = (varCleanType v, Option.map cleanExp e)
    and cleanGlobalDcl gdcl =
	(case gdcl
	  of AST.D_Var v => AST.D_Var (cleanVarDcl v)
	   | AST.D_Func( v, vs,s) => AST.D_Func(varCleanType v, List.map varCleanType vs, cleanStm s)
	   | AST.D_DiffFunc(v, vs, s) => AST.D_DiffFunc(varCleanType v, List.map varCleanType vs, cleanExp s)
	(* end case*))
    and cleanStm stm =
	(case stm
	  of AST.S_Block(stms) => AST.S_Block (List.map cleanStm stms)
	   | AST.S_Decl(( v,e)) => AST.S_Decl((varCleanType v,Option.map cleanExp e))
	   | AST.S_IfThenElse(e1,s1,s2) => AST.S_IfThenElse(cleanExp e1,cleanStm s1, cleanStm s2)
	   | AST.S_Foreach(( v,e),s) => AST.S_Foreach((varCleanType v,cleanExp e), cleanStm s)
	   | AST.S_Assign((v,e), e1) => AST.S_Assign((varCleanType v, e), cleanExp e1)
	   | AST.S_New(a, exps) => AST.S_New(a, List.map cleanExp exps)
	   | _ => stm
	(*end case*))
    and cleanExp exp =
	(case exp
	  of AST.E_Prim( v, tymeta, exps, ty) => AST.E_Prim(varCleanType v, tymeta, List.map cleanExp exps, cleanTy ty) (*NOTE:  List.map cleanMetaTy*)
	   | AST.E_Select(e, (v,s)) => AST.E_Select(cleanExp e, (varCleanType v,s))
	   | AST.E_Apply( (v,s), exps, ty) => AST.E_Apply((varCleanType v, s), List.map cleanExp exps, cleanTy ty)
	   | AST.E_Comprehension(e, ( v, eit), ty) => AST.E_Comprehension(cleanExp e, (varCleanType v, cleanExp eit), cleanTy ty)
	   | AST.E_ParallelMap(e, v1, v2,ty) => AST.E_ParallelMap(cleanExp e, varCleanType v1, varCleanType v2, cleanTy ty)
	   | AST.E_Tensor(exps, ty) => AST.E_Tensor(List.map cleanExp exps, cleanTy ty)
	   | AST.E_Field(exps, ty) => AST.E_Field(List.map cleanExp exps, cleanTy ty)
	   | AST.E_Seq(exps, ty) => AST.E_Seq(List.map cleanExp exps, cleanTy ty)
	   | AST.E_Slice(exp, exps, ty) => AST.E_Slice(cleanExp exp, List.map (Option.map cleanExp) exps, cleanTy ty)
	   | AST.E_Cond(e1,e2,e3, ty) => AST.E_Cond(cleanExp e1,cleanExp e2,cleanExp e3, cleanTy ty)
	   | AST.E_Orelse(e1,e2) => AST.E_Orelse(cleanExp e1,cleanExp e2)
	   | AST.E_Andalso(e1,e2) => AST.E_Andalso(cleanExp e1,cleanExp e2)
	   | AST.E_LoadNrrd(s) => AST.E_LoadNrrd(s) (*Question: Should we even need to worry about this?*)
	   | AST.E_Coerce{srcTy,dstTy, e} =>
	     let
	      val srcTy' = cleanTy srcTy
	      val dstTy' = cleanTy dstTy
	      val e' = cleanExp e
	      val res =AST.E_Coerce {srcTy=srcTy', dstTy=dstTy', e=e'}

	     in
	      (case (srcTy, dstTy)
		of (Ty.T_Named(_), _) => e'
		 | (_, Ty.T_Named(_)) => e'
		 | _ => res
	      )

	     end
	   | AST.E_LoadFem(data,someExp, someExp') => AST.E_LoadFem(data, Option.map cleanExp someExp, Option.map cleanExp someExp')
	   | AST.E_ExtractFem(e,data) => AST.E_ExtractFem(cleanExp e, data)
	   | AST.E_ExtractFemItem(e,ty,opt) => AST.E_ExtractFemItem(cleanExp e, ty, opt)
	   | _ => exp
	(* end case*))


	 
end
