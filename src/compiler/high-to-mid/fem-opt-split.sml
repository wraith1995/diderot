(* fem-to-low.sml
 *
 * This module handles the translation of fem ops to other fem ops.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure FemOptSplit : sig
	   val indexedDataLowering : MidIR.var * FemOpt.femOption * MidIR.var * MidIR.var -> MidIR.assign list
	   val indexedDataLoweringAvail : AvailRHS.t * MidIR.var * FemOpt.femOption * MidIR.var * MidIR.var * MidIR.var -> unit
	  end = struct


    structure IR = MidIR
    structure V = IR.Var
    structure Ty = MidTypes
    structure Op = MidOps
    structure FT = FemData
    structure FO = FemOpt

   fun assignOp (avail, pre, ty, opss, args) =
          AvailRHS.addAssign(avail, pre, ty, IR.OP(opss, args))


    fun imul' (avail, pre, a, b) = assignOp (avail, pre, Ty.IntTy, Op.IMul, [a, b])
    fun iadd' (avail, pre, a, b) = assignOp (avail, pre, Ty.IntTy, Op.IAdd, [a, b])
    fun ilit' (avail, pre, n) =
          AvailRHS.addAssign (avail, pre, Ty.IntTy, IR.LIT(Literal.Int(IntLit.fromInt n)))
    fun imul (avail, a, b) = imul' (avail, "mulRes", a, b)
    fun iadd (avail, a, b) = iadd' (avail, "addRes", a, b)
    fun ilit (avail, n) = ilit' (avail, "ilit", n)


    fun indexArrayLoad	(avail, len, opt, indexSource, indexVar) = 
	let
	 (* extract fem data*)
	 val targetType = Ty.SeqTy(Ty.IntTy, SOME(len))
	 val newOpt = Op.ExtractFemItem2(Ty.IntTy, targetType, opt) (*given a int, extracts that index from the data*)
	 val result = assignOp(avail, "femIndicies", targetType, newOpt, [indexSource, indexVar])
	in
	 (targetType, result)
	end

    fun findSources(avail, source, data) =
	(case FT.dependencyOf(data)
	  of NONE => (source, source)
	   | SOME(a) =>
	     let
	      val (indexSourceTy, currentTy) = (Ty.FemData(a), Ty.FemData(data))
	      val indexSource = assignOp(avail, "indexSource", indexSourceTy, Op.ExtractFem(indexSourceTy, currentTy), [source])
	     in
	      (source, indexSource)
	     end)

    fun dofLoad (avail, shape, opt, dofSource, seqVar, seqTy) =
	let
	 val targetType = Ty.TensorTy(shape)
	 val newOpt = Op.ExtractFemItem2(seqTy, targetType, opt)
	 val result = assignOp(avail, "femDofs", targetType, newOpt, [dofSource, seqVar])
	in
	 result
	end
    fun indexedDataLoweringAvail(avail, lhs, opt, dataArg, indexSrcVar, indexVar) =
	let
	 val SOME(indexingOpt, dataGetOpt) = FO.splitDataOpt opt
	 val (_, data) = opt
	 val (dataSource, indexSource) = (dataArg, indexSrcVar) 
	 val len = FO.findIndexLength indexingOpt
	 val resultShape = FO.findTargetShape dataGetOpt

	 val (indexTy, indexLoad) = indexArrayLoad(avail, len, indexingOpt, indexSource, indexVar)
	 val result = dofLoad(avail, resultShape, dataGetOpt, dataSource, indexLoad, indexTy)
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result))
	end

    fun indexedDataLowering (lhs, opt, dataArg, indexVar) =
	let
	 val avail = AvailRHS.new ()
	 val SOME(indexingOpt, dataGetOpt) = FO.splitDataOpt opt
	 val (_, data) = opt
	 val (dataSource, indexSource) = findSources(avail, dataArg, data)
	 val len = FO.findIndexLength indexingOpt
	 val resultShape = FO.findTargetShape dataGetOpt
					      
	 val (indexTy, indexLoad) = indexArrayLoad(avail, len, indexingOpt, indexSource, indexVar)
	 val result = dofLoad(avail, resultShape, dataGetOpt, dataSource, indexLoad, indexTy)
	in
	 AvailRHS.addAssignToList (avail, (lhs, IR.VAR result));
	 List.rev (AvailRHS.getAssignments avail)
	end


		     
end
