(* mid-contract.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *
 * Expand intervals:
scalars = consts, scalar tensor
intervals = interval tensor
termains = intervals | scalars
ops = op1 terminas | op2 terminas x 2 | op3 terminas x 3 | opn terminas x end
ret = ops | E.Sum ops

ONE ISSUE: Make sure sizing and special ops
(case op1
      increasing, decreasing (neg, expr, sqrt, pow) -> obvious
      abs -> min(max(min, 0), max), max(abs(max), abs(min)) -- vectorized
      rest - two 2 vec rep, unwrap into interval functions, rewrap
      ignore: sgn

      op2:
      Sub: clear - divide into mins and maxes - make sure to do scalars
      div: make the division tensor, nan test (vectorized), and run multiplication - make sure to do scalars
      ignore: min, max

      op3: clamp ignored
      opn: ADD: clear, but scalars
      opn: PROD:
      --build the 2^n combinations for the n interval operations, add in the scalars, take mins and maxes as vectors.
      )

Sum: just substitute back in
SUM: any chances for vectorization we miss?

AA:
(case op1
      Just use interval versions! Avoid pow when you fix that and in probe

      op2:
      sub - clear, but again the filter
      div: convert to interval, use nan test, convert back to affine -- make div clearer (LAST)

      op3: ignore:

      opn:ADD: clear, but again the filter
      opn: PROD: for the operations, we lift and do all the non scalar ones. Then do the scalar ones. Then combine.
      
      )

SUM: you just subtitute back in and do it a few times. Easy and peasy.
SUM: any chances for vectorization we miss?

NOTE: what about scalar ops in the div, add, 
Plan: fix probe, power
Plan: ops to mmid
Plan: make util functions we will need and framework (we can use mk ops though)
Plan: Do interval, push the new functions down to low and vectorize and cxx
Plan: affine later.
Plan: Push and test
PLan:fast-math
Plan:parser
.
 *)

structure MidIntervalExpand : sig

    val transform : MidIR.program -> MidIR.program

  end = struct

    structure IR = MidIR
    structure Op = MidOps
    structure Ty = MidTypes
    structure V = IR.Var
    structure ST = Stats

  (********** Counters for statistics **********)
    val cntUnused               = ST.newCounter "mid-opt:unused"
    val firstCounter            = cntUnused
    val lastCounter             = cntUnused

    structure UnusedElim = UnusedElimFn (
        structure IR = IR
        val cntUnused = cntUnused)

    fun useCount (IR.V{useCnt, ...}) = !useCnt

  (* adjust a variable's use count *)
    fun incUse (IR.V{useCnt, ...}) = (useCnt := !useCnt + 1)
    fun decUse (IR.V{useCnt, ...}) = (useCnt := !useCnt - 1)
    fun use x = (incUse x; x)

    fun getRHSOpt x = (case V.getDef x
           of IR.OP arg => SOME arg
            | _ => NONE
          (* end case *))

    fun doAssign (lhs, IR.OP rhs) = raise Fail "umm"
      | doAssign _ = raise Fail "umm"




    structure Rewrite = RewriteFn (
      struct
        structure IR = IR
        val doAssign = doAssign
        fun doMAssign _ = NONE
        val elimUnusedVars = UnusedElim.reduce
      end)

    val transform = Rewrite.transform

  end
