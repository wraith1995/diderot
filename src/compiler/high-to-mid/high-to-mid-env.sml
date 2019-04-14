(* high-to-mid-env.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)
structure SrcIR = HighIR
structure SrcTy = HighTypes
structure DstIR = MidIR
structure DstTy = MidTypes

structure HME = struct

fun cvtTy SrcTy.BoolTy = DstTy.BoolTy
      | cvtTy SrcTy.StringTy = DstTy.StringTy
      | cvtTy SrcTy.IntTy = DstTy.intTy
      | cvtTy (SrcTy.TensorTy dd) = DstTy.tensorTy dd
      | cvtTy (SrcTy.TupleTy tys) = DstTy.TupleTy(List.map cvtTy tys)
      | cvtTy (SrcTy.SeqTy(ty, n)) = DstTy.SeqTy(cvtTy ty, n)
      | cvtTy (SrcTy.ImageTy info) = DstTy.ImageTy info
      | cvtTy (SrcTy.StrandTy n) = DstTy.StrandTy n
      | cvtTy SrcTy.KernelTy = DstTy.KernelTy
    (* we replace Field operations by 0, so the types are mapped to int *)
      | cvtTy SrcTy.FieldTy = DstTy.intTy
      | cvtTy (SrcTy.FemData(data)) = DstTy.FemData(data)

  (* instantiate the translation environment *)
    
end

structure Env = TranslateEnvFn (
 struct
 structure SrcIR = SrcIR
 structure DstIR = DstIR
 val cvtTy = HME.cvtTy
 end)
