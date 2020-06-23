(* simplify-fem-exp.sml
 *
 *
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)


structure SimplifyFemExp : sig
	   val transform : Simple.program -> Simple.program
	  end = struct

  structure S = Simple
  structure Ty = SimpleTypes
  structure V = SimpleVar
  structure F = SimpleFunc


fun cvtBlockWorld(S.Block{code, ...}) =
    let
     fun doMapReduce(S.MapReduce{result, reduction, mapf, args, source, domain}) =
	 let
	  val S.Func{f, params,body} = mapf
	  val mapf' = S.Func{f=f, params=params, body = cvtBlockWorld body}
	 in
	  S.MapReduce{result=result,
		      reduction=reduction,
		      mapf = mapf',
		      args = args,
		      source = source,
		      domain = domain}
	 end
     fun mapStms (S.S_Var(v, SOME(S.E_ExtractFemItem(pos, t, (FemOpt.WorldPos, m))))) = Util.cvtWorldPos(pos, v, false)
       | mapStms (S.S_Assign(v, S.E_ExtractFemItem(pos, t, (FemOpt.WorldPos, m)))) = Util.cvtWorldPos(pos, v, true)
       | mapStms (S.S_MapReduce(mrs)) = [S.S_MapReduce(List.map doMapReduce mrs)]
       | mapStms (S.S_IfThenElse(v,b1,b2)) = [S.S_IfThenElse(v, cvtBlockWorld b1, cvtBlockWorld b2)]
       | mapStms (S.S_Foreach(v1, v2, b)) = [S.S_Foreach(v1, v2, cvtBlockWorld b)]
       | mapStms (s) = [s]

    in
     S.Block{code = List.concatMap mapStms code, props=PropList.newHolder ()}
    end

fun transformWorld prog =
    let
     val S.Program{props, consts, inputs, constInit, globals,
		   globInit, funcs, strand, create, start, update} = prog

     val globInit' = cvtBlockWorld globInit
     val start' = Option.map cvtBlockWorld start
     val update' = Option.map cvtBlockWorld update
     val create' = Create.map cvtBlockWorld create
     val funcs' = List.map (fn S.Func{f, params, body} => S.Func{f=f, params=params, body=cvtBlockWorld body}) funcs

     val S.Strand{name, params, spatialDim, state, stateInit, startM, updateM, stabilizeM} = strand
     val strand' = S.Strand{name=name, params=params, spatialDim=spatialDim, state=state,
			    stateInit=cvtBlockWorld stateInit, startM = Option.map cvtBlockWorld startM,
			    updateM = cvtBlockWorld updateM, stabilizeM = Option.map cvtBlockWorld stabilizeM}
								       
    in
     S.Program{props=props, consts=consts, inputs=inputs, globals=globals, constInit=constInit, (*no _pos things*)
	       globInit=globInit', start=start', update=update', create=create', funcs = funcs', strand=strand'
	      }
    end
(*compose to allow other optimizations of this nature:*)
val transform = transformWorld
end
