(* fem-to-low.sml
 *
 * This module translates fem-opts that load sequences/arrays of DOFs to operations that directly build these objects; analogous to load-voxels. 

 * TODO: In some cases, we have an option of just directly copying a sequence and shouldn't do the translation; we should figure out when? 
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure FemToLow : sig
	   val expandOp : LowIR.var * LowIR.var * LowIR.var * LowTypes.ty * LowTypes.ty * FemOpt.femOption -> LowIR.assign list
	  end = struct


    structure IR = LowIR
    structure V = IR.Var
    structure Ty = LowTypes
    structure Op = LowOps
    structure MTy = MidTypes
    structure FT = FemData
    structure FO = FemOpt

    (* I'm going to shamelessly borrow code from load-voxels.sml; we can factor out commmon stuff later.*)

    fun assignOp (avail, pre, ty, opss, args) =
        AvailRHS.addAssign(avail, pre, ty, IR.OP(opss, args))
    fun assignCons (_, _, []) = raise Fail "empty cons"
      | assignCons (avail, pre, args as x::_) = let
       val ty = Ty.TensorTy(List.length args :: Ty.tensorShape(V.ty x))
      in
       AvailRHS.addAssign(avail, pre, ty, IR.CONS(args, ty))
      end
    fun imul' (avail, pre, a, b) = assignOp (avail, pre, Ty.IntTy, Op.IMul, [a, b])
    fun iadd' (avail, pre, a, b) = assignOp (avail, pre, Ty.IntTy, Op.IAdd, [a, b])
    fun ilit' (avail, pre, n) =
          AvailRHS.addAssign (avail, pre, Ty.IntTy, IR.LIT(Literal.Int(IntLit.fromInt n)))
    fun imul (avail, a, b) = imul' (avail, "mulRes", a, b)
    fun iadd (avail, a, b) = iadd' (avail, "addRes", a, b)
    fun ilit (avail, n) = ilit' (avail, "ilit", n)

    fun consFunc(vars, ty) = (case ty
			       of Ty.TensorTy[] => (case vars
						     of [v] => IR.VAR(v)
						      | _ => raise Fail "impossible")
				| Ty.TensorTy[_] => IR.CONS(vars, ty)
				| _ => raise Fail "impossible")
    (*Note: it might be better to not lower the index transfer and make it a mem copy:
     1. Avail rewrites the indexing to occur when we access tensors, which is bad for locality
     2. It might sometimes be more efficient to memcpy*)			  
    fun expandExtractIndex(lhs, seqTy, seqOpt, indexSource, indexStart) =
	let
	 val avail = AvailRHS.new ()
	 val Ty.SeqTy(Ty.IntTy, SOME n ) = seqTy
	 val newOpt = Op.ExtractFemItem2(Ty.IntTy, Ty.IntTy, seqOpt)
	 val indexStart' = imul(avail, indexStart, ilit(avail,n))
	 fun loadIndex idx = assignOp(avail, "idx", Ty.IntTy, newOpt, [indexSource, iadd(avail, indexStart', ilit(avail, idx))])
	 val idxes = List.tabulate (n, loadIndex)
	 val result = AvailRHS.addAssign (avail, "seq", seqTy, IR.SEQ(idxes, seqTy))
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result));
	 List.rev (AvailRHS.getAssignments avail)
	end

    fun expandDofsChoice(lhs, seqTy, loadOpt, loadSource, indexSeq, loadDims ) =
	let
	 val avail = AvailRHS.new ()
	 val perDofShape = FO.findTargetShape loadOpt
	 val dofSize = List.foldr (fn (x,y) => x*y) 1 perDofShape
				  
	 val _ = print ("[" ^ (String.concatWith ", " (List.map Int.toString perDofShape)) ^ "]\n")
	 val perDofShapeRev = List.rev perDofShape
	 val dim = List.length perDofShapeRev
	 val _ = if dim >= loadDims
		 then ()
		 else raise Fail "trying to load portions of a dof that are large than the dof"
	 val loadingShape = List.rev(List.take(perDofShapeRev, loadDims))
	 val loadTensorTy = Ty.TensorTy( loadingShape)
	 val accessShape = List.rev(List.drop(perDofShapeRev, loadDims))
	 val accessTensorTy = Ty.TensorTy( accessShape)
	 val Ty.SeqTy(Ty.IntTy, SOME numDofs ) = seqTy
	 val shift = List.foldr (fn (x,y) => x*y) 1 loadingShape (*size of the tensors that we are loading*)
	 val count = List.foldr (fn (x,y) => x*y) 1  accessShape (*size of the tensor that we load into*)

	 val newOpt = Op.ExtractFemItem2(Ty.IntTy, loadTensorTy, loadOpt)
	 fun loadTensor index = assignOp(avail, "dof_load", loadTensorTy, newOpt, [loadSource, index])
	 (*build an array of ops, get cons out, argument via index*)
	 fun buildDof idx =
	     let
	      fun getIndex index =  iadd(avail, ilit(avail, shift * index), imul(avail, ilit(avail,  dofSize), assignOp (avail, "idx", Ty.IntTy, Op.Subscript Ty.IntTy, [indexSeq, ilit (avail, idx)])))

	      val indicies = List.tabulate(count, getIndex)
	      val accessArray = ArrayNd.fromList'(indicies, accessShape)
	      val convert = loadTensor 
	      fun group(vars, n) =
		  let
		   val v::vs = vars
		   val ty = IR.Var.ty(v)
		   val Ty.TensorTy(shape') = ty

		  in
		   assignCons(avail, "interCons", vars)
		  end

	     in
	      ArrayNd.convertToTree(accessArray, convert, group) handle exn => raise exn
	     end

	 val result = assignCons(avail, "dofs", List.tabulate(numDofs, buildDof))
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result));
	 List.rev (AvailRHS.getAssignments avail)
	end


    fun expandDofsMemCpy(lhs, seqTy, loadOpt, loadSource, indexSeq ) =
	let
	 val avail = AvailRHS.new ()
	 val perDofShape = FO.findTargetShape loadOpt
	 val Ty.SeqTy(Ty.IntTy, SOME n ) = seqTy
	 val shift = List.foldr (fn (x,y) => x*y) 1 perDofShape
	 val newTensorTy = Ty.TensorTy(perDofShape)
	 val newOpt = Op.ExtractFemItem2(Ty.IntTy, newTensorTy, loadOpt)
	 fun getIndex idx =  imul(avail, ilit(avail, shift), assignOp (avail, "idx", Ty.IntTy, Op.Subscript Ty.IntTy, [indexSeq, ilit (avail, idx)]))
	 val indicies = List.tabulate(n, getIndex)
	 fun loadTensor idx = assignOp(avail, "dof", newTensorTy, newOpt, [loadSource, idx])
	 val result = assignCons(avail, "dofs", List.map  loadTensor indicies )
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result));
	 List.rev (AvailRHS.getAssignments avail)
	end

    fun expandDofs(arg as (lhs, seqTy, loadOpt, loadSource, indexSeq)) =
	(case IR.Var.ty lhs
	  of Ty.TensorTy([d]) => expandDofsMemCpy arg
	  |  Ty.TensorTy(_)=> expandDofsChoice(lhs, seqTy, loadOpt, loadSource, indexSeq, 0)
	 (*loads scalars one by one*)
	 (*TODO: is there a better way to make this choice?*)
	 (*end case*))
	


    fun expandOp(lhs, arg1, arg2, inputTy, outputTy, (opt, d)) =
	(case opt
	  of FemOpt.ExtractIndex => expandExtractIndex(lhs, outputTy, (opt, d), arg1, arg2)
	   | FemOpt.ExtractDof => expandDofs(lhs, inputTy, (opt, d), arg1, arg2)
	   | _ => raise Fail "impossible")
end
