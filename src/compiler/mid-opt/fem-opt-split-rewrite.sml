(* fem-to-low.sml
 *
 * This module handles the translation of fem ops to other fem ops.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)


structure FemOptSplitRewrite : sig
	   val transform : MidIR.program -> MidIR.program
	  end = struct
structure IR = MidIR
structure Op = MidOps
structure Ty = MidTypes


fun doRewrite(lhs, IR.OP(Op.ExtractFemItemN(tys, outTy, opt as (FemOpt.ExtractDofs, data), _, _, _, _), [srcVar1, srcVar2, intVar])) =
    let
     val avail = AvailRHS.new ()
     val _ = FemOptSplit.indexedDataLoweringAvail(avail, lhs, opt, srcVar1, srcVar2, intVar)
     val stmts = List.rev (AvailRHS.getAssignments avail)
    in
     SOME(stmts)
    end
  | doRewrite _ = NONE
structure ST = Stats
val cntUnused               = ST.newCounter "mid-opt:unused'"
structure UnusedElim = UnusedElimFn (
 structure IR = IR
 val cntUnused = cntUnused)

		    
structure Rewrite = RewriteFn (
 struct
 structure IR = IR
 val doAssign = doRewrite

 fun doMAssign _ = NONE
 fun elimUnusedVars _ = false
 end)
val transform = Rewrite.transform		    

end
