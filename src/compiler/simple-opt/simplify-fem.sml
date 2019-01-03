(* simplify-fem.sml
 *
 *
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)


structure SimplifyFem : sig

    val transform : Simple.program -> Simple.program

  end = struct

structure S = Simple
structure TY = SimpleTypes
structure V = SimpleVar
structure VSet = SimpleVar.Set
structure VMap = SimpleVar.Map
structure FT = FemData
structure FO = FemOpt



    fun transform prog = prog
end
