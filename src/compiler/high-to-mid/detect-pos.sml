(* detect-pos.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)

structure DetectPos : sig
	   val removePos : Env.t * MidIR.assign -> MidIR.assign

	  end = struct
structure IR = MidIR
structure HIR = HighIR
structure Op = MidOps
structure HOp = HighOps
structure V = IR.Var
structure HV = HIR.Var
structure ST = Stats
structure E = Ein
structure Ty = MidTypes
structure HTy = HighTypes
(********** Counters for statistics **********)
val cntUnused = ST.newCounter "high-to-mid:pos-replace"


			      
fun detectTransform(e : Ein.ein, args : IR.var list) : (IR.var * IR.var * IR.var) option =
    let
     val E.EIN{params, index, body} = e
     fun vTy x = IR.Var.ty (List.nth(args, x))
     fun vS(x,y) = x=y orelse IR.Var.same(List.nth(args, x), List.nth(args, y))
    in
     (case body
       of E.Probe(E.Fem(E.Plain(_, _, NONE), index, idxSrc, dofSrc, [E.V _], []), E.Tensor(pid, _))
	  => (case (vTy index, vTy idxSrc, vTy dofSrc)
	       of (Ty.IntTy, Ty.FemData(FemData.Mesh _), Ty.FemData(FemData.Mesh _))
		  => if vS(idxSrc, dofSrc)
		     then
		      SOME(List.nth(args, index), List.nth(args, idxSrc), List.nth(args, pid))
		     else NONE
		| _ => NONE
	     (* end case *))
	| _ => NONE
     (* end case *))
    end

(*Fix me: check indSrc for general case with 4 tuple instead of 3-tuple*)
fun findTupleBinding env (indM : IR.var, indSrcM : IR.var, posVarM : IR.var) : HIR.var option =
    let
     val fstOpt = Option.map (fn (x,y) => x)
     val ind = fstOpt (Env.invertVar(env, indM))
     val indSrc = fstOpt(Env.invertVar(env, indSrcM))
     val posVar = fstOpt(Env.invertVar(env, posVarM))
     val indRhs = Option.map HV.getDef ind
     val posRhs = Option.map HV.getDef posVar
     fun validateTupleArgs j ((ty, i)) = j=i andalso (case ty
						       of HTy.TupleTy([HTy.TensorTy _, HTy.IntTy, HTy.IntTy]) => true
							| _ => false)
				   
    in
     (case (indRhs, posRhs)
       of (SOME(HIR.OP(HOp.Select arg1, [v1])), SOME(HIR.OP(HOp.Select arg2, [v2])))
	  => if HV.same(v1, v2) andalso (validateTupleArgs 1 arg1) andalso (validateTupleArgs 0 arg2)
	     then SOME(v1)
	     else NONE
	| _ => NONE
     (* end case*))
    end
      
fun findStateVarBinding (env : Env.t) (posTuple : HIR.var) =
    (case HV.getDef posTuple
      of HIR.STATE(strand, stateVar) =>
	 (case Env.getPosVars(env)
	   of SOME((srcPos, dstPos), (srcPos', dstPos')) =>
	      if IR.StateVar.same(Env.renameSV(env, stateVar), dstPos)
	      then SOME(Option.map (fn x => Env.rename(env, x)) strand, dstPos')
	      else NONE
	 )
       | _ => NONE
    (* end case *))

fun removePos(env, assign as (v, IR.EINAPP(ein, args))) =
    let
     val findTupleBinding = findTupleBinding env
     val findStateVarBinding = findStateVarBinding env
     val detected = detectTransform(ein, args)
     val tplBind = Option.mapPartial findTupleBinding detected
     val stateVar = Option.mapPartial findStateVarBinding tplBind
    in
     (case stateVar
       of NONE => assign
	| SOME(strand, sv) => (v, IR.STATE(strand, sv))
     (* end case *))
    end
  | removePos(env, e) = e
end
