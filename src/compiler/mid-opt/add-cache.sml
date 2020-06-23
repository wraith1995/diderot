(* add-cache.sml
 *
 * Add cache for dofs
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)

structure AddCache : sig
	   val translate : MidIR.program -> MidIR.program
	  end = struct

structure SrcIR = MidIR
structure SrcTy = MidTypes
structure DstIR = MidIR
structure DstTy = MidTypes

structure Env = TranslateEnvFn (
 struct
 structure SrcIR = SrcIR
 structure DstIR = DstIR
 val cvtTy = fn x => x
 end)		    


structure Trans = TranslateScopedFn(
 struct
 open Env
 type scope = bool
 val expand = fn (x,y,z) => raise Fail "umm"
 val mexpand = fn (x,y,z) => raise Fail "umm"

 val constInitScope = false
 val funcScope = false
 val globInitScope = false
 val globStart = false
 val globUpdate = false
 val createScope = false


 val strandInit = false
 val strandStart = false
 val strandUpScope = true
 val strandStabScope = false
 end
)


fun translate prog =
    let
     val prog = Trans.translate prog
    in
     (MidCensus.init prog;
     prog)
    end
end
