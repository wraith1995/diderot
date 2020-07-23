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
Plan:vectorization
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
    structure E = Ein
    structure EU = EinUtil
    structure EP = EinPP

    fun cvtTy (Ty.TensorTy(ts, SOME j)) = if j = 1
					  then Ty.TensorTy(2 :: ts, NONE)
					  else if j > 1
					  then Ty.TensorTy(j :: ts, NONE)
					  else raise Fail "impossible Tensor Ty"
      | cvtTy (Ty.TensorTy(ts, NONE)) = Ty.TensorTy(ts, NONE)
      | cvtTy (Ty.TupleTy(ts)) = Ty.TupleTy(List.map cvtTy ts)
      | cvtTy (Ty.SeqTy(t, opt)) = Ty.SeqTy(cvtTy t, opt)
      | cvtTy (t) = t (*bool, string, int, img, strand, fem, kern*)

    structure Env = TranslateEnvFn (
     struct
     structure SrcIR = IR
     structure DstIR = IR
     val cvtTy = cvtTy
     end)		    
		     

  (********** Counters for statistics **********)
    val cntUnused               = ST.newCounter "mid-opt:unusedp"
    val firstCounter            = cntUnused
    val lastCounter             = cntUnused

    (*Plan: opts here, expands here, deal with new mid in mid to low, vectorization, c++, lifting, parsing, affine fix
    State vars? Globals? Might need a new rewrite - rewrite and change. Translate is a good idea.
     *)
			
    fun doAssign (lhs, IR.OP rhs) = raise Fail "umm"
      | doAssign _ = raise Fail "umm"




    fun expand (env, (y, rhs)) = let
          fun assign rhs = [IR.ASSGN(Env.rename (env, y), rhs)]
          in
            case rhs
             of IR.GLOBAL x => assign (IR.GLOBAL(Env.renameGV(env, x)))
              | IR.STATE(NONE, fld) => assign (IR.STATE(NONE, Env.renameSV(env, fld)))
              | IR.STATE(SOME x, fld) =>
                  assign (IR.STATE(SOME(Env.rename(env, x)), Env.renameSV(env, fld)))
              | IR.VAR x => assign (IR.VAR(Env.rename(env, x)))
              | IR.LIT lit => assign (IR.LIT lit)
              | IR.OP(rator, args) => raise Fail "oops"
                  (* List.map IR.ASSGN (expandOp (env, Env.rename (env, y), rator, args)) *)
              | IR.CONS(args, ty) => assign (IR.CONS(Env.renameList(env, args), cvtTy ty))
              | IR.SEQ(args, ty) => assign (IR.SEQ(Env.renameList(env, args), cvtTy ty))
              | IR.EINAPP(rator, args) => raise Fail "oops"
              | IR.APPLY(f, args) =>
                  assign (IR.APPLY(Env.renameFV(env, f), Env.renameList(env, args)))
              | _ => raise Fail("bogus rhs for ASSIGN: " ^ IR.RHS.toString rhs)
            (* end case *)
    end



  (* expand a IR multi-assignment to a IR CFG *)
    fun mexpand (env, (ys, rhs)) = let
          fun massign rhs = let
                val nd = IR.Node.mkMASSIGN(Env.renameList(env, ys), rhs)
                in
                  IR.CFG{entry=nd, exit=nd}
                end
          fun mkOP (rator, xs) = massign(IR.OP(rator, Env.renameList(env, xs)))
          in
            case rhs
             of IR.OP(Op.EigenVecs2x2, xs) => mkOP (Op.EigenVecs2x2, xs)
              | IR.OP(Op.EigenVecs3x3, xs) => mkOP (Op.EigenVecs3x3, xs)
              | IR.OP(Op.KillAll, []) => mkOP (Op.KillAll, [])
              | IR.OP(Op.StabilizeAll, []) => mkOP (Op.StabilizeAll, [])
              | IR.OP(Op.Print tys, xs) => mkOP (Op.Print(List.map cvtTy tys), xs)
	      | IR.OP(Op.Save args, xs) => mkOP(Op.Save args, xs)
              | IR.MAPREDUCE mrs => let
                  val mrs = List.map
                        (fn (r, f, xs) => (r, Env.renameFV(env, f), Env.renameList(env, xs)))
                          mrs
                  in
                    massign (IR.MAPREDUCE mrs)
                  end
              | _ => raise Fail("bogus rhs for MASSIGN: " ^ IR.RHS.toString rhs)
            (* end case *)
          end

    structure Trans =  TranslateFn (
      struct
        open Env
        val expand = IR.CFG.mkBlock o expand
        val mexpand = mexpand
      end)

structure Promote = PromoteFn (IR)

    fun transform prog =
	let
	 val prog = Trans.translate prog
	in
	 MidCensus.init prog;
	 Promote.transform prog
	end

  end
