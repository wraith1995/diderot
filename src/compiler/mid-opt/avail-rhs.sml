(* avail-rhs.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)

(* tracking available LowIR rhs expressions *)
structure AvailRHS = AvailRHSFn (
    val phase = "mid-opt"
    structure IR = MidIR)
