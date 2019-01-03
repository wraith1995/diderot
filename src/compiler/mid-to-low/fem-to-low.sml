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


						  
    fun expandExtractIndex(lhs, seqTy, seqOpt, indexSource, indexStart) =
	let
	 val avail = AvailRHS.new ()
	 val Ty.SeqTy(Ty.IntTy, SOME n ) = seqTy
	 val newOpt = Op.ExtractFemItem2(Ty.IntTy, Ty.IntTy, seqOpt)
	 fun loadIndex idx = assignOp(avail, "idx", Ty.IntTy, newOpt, [indexSource, iadd(avail, indexStart, ilit(avail, idx))])
	 val idxes = List.tabulate (n, loadIndex)
	 val result = AvailRHS.addAssign (avail, "seq", seqTy, IR.SEQ(idxes, seqTy))
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result));
	 List.rev (AvailRHS.getAssignments avail)
	end

    fun expandDofs(lhs, seqTy, loadOpt, loadSource, indexSeq ) =
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

    fun expandOp(lhs, arg1, arg2, inputTy, outputTy, (opt, d)) =
	(case opt
	  of FemOpt.ExtractIndex => expandExtractIndex(lhs, outputTy, (opt, d), arg1, arg2)
	   | FemOpt.ExtractDof => expandDofs(lhs, inputTy, (opt, d), arg1, arg2)
	   | _ => raise Fail "impossible")
end
