(* fem-opt.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * Information about meshes, function spaces, and fuctions in function spaces.
 *)


structure FemOpt : sig

	   datatype femOpts = NumCell | ExtractDofs | ExtractDofsSeq | ExtractIndices |  ExtractIndex  | ExtractDof | Cells | CellIndex  | PromoteCell

	   type femOption = femOpts * FemData.femType
					    
	   val same : femOption * femOption -> bool
	   val hash : femOption -> word
	   val toString : femOption -> string
	   (*Finds the number of dofs an operator might be extracting*)
	   val findIndexLength : femOption -> int
	   (*Takes an operator and splits it into two (high-to-mid)*)
	   val splitDataOpt : femOption -> (femOption * femOption) option
	   (*Takes an operator and lowers it down from the mid-ir to the low-ir*)
	   val lowerDataOpt : femOption -> femOption option
	   (*Find the shape of a tensor that is the result a fem opt*)
	   val findTargetShape : femOption -> int list

	   datatype femField =  Transform | RefField | InvTransform | TRefField | Field						  
	   val fieldString : femField -> string
									     
	  end = struct
datatype femOpts = Cells
		 | ExtractDofs
		 | ExtractIndices |  ExtractDofsSeq
		 | NumCell  | ExtractIndex | ExtractDof (* primitize all*)
		 | CellIndex | PromoteCell (* primitize, cells*)

datatype femField =  Transform | RefField | InvTransform | TRefField | Field


fun fieldString (Transform) = "Transform"
  | fieldString (RefField) = "refField"
  | fieldString (InvTransform) = "InvTransform"
  | fieldString (TRefField) = "TRefField"
  | fieldString (Field) = "Field"

type femOption = femOpts * FemData.femType
structure FT = FemData
(* add airity*)			 
fun toStringOpt v =
    (case v
      of NumCell => "NumCell"
       | Cells => "Cells"  
       | CellIndex => "CellIndex" 
       | ExtractDofs => "ExtractDofs"
       | ExtractIndex => "ExtractIndex"
       | ExtractDof => "ExtractDof"
       | ExtractIndices => "ExtractIndices"
       | ExtractDofsSeq => "ExtractDofsSeq"
       | PromoteCell => "PromoteCell"
    (* end case*))

fun toString (v, d) = toStringOpt(v) ^ "(" ^ (FT.toString d) ^ ")"
         

fun arity (NumCell) = 1
  | arity (Cells) = 1
  | arity (CellIndex) = 1
  | arity (ExtractDofs) = 2
  | arity (ExtractIndices) = 2
  | arity (ExtractDofsSeq) = 2
  | arity (ExtractIndex) = 2
  | arity (ExtractDof) = 2
  | arity (PromoteCell) = 2



fun hash (NumCell, d) = 0w1 + FT.hash d
  | hash (Cells, d) = 0w2 + FT.hash d
  | hash (CellIndex, d) = 0w3 + FT.hash d
  | hash (ExtractIndex, d) = 0w5 + FT.hash d
  | hash (ExtractDof, d) = 0w7 + FT.hash d
  | hash (ExtractDofs, d) = 0w11 + FT.hash d
  | hash (ExtractIndices, d) = 0w13 + FT.hash d
  | hash (ExtractDofsSeq, d) = 0w17 + FT.hash d
  | hash (PromoteCell, d) = 0w19 + FT.hash d
		     
fun same ((v1, d1),(v2, d2)) = FT.same(d1,d2) andalso
    (case (v1,v2)
      of (NumCell, NumCell) => true
       | (Cells, Cells) => true
       | (CellIndex, CellIndex) => true
       | (ExtractDofs, ExtractDofs) => true
       | (ExtractIndex, ExtractIndex) => true
       | (ExtractIndices, ExtractIndices) => true
       | (ExtractDof, ExtractDof) => true
       | (ExtractDofsSeq,ExtractDofsSeq) => true
       | (PromoteCell, PromoteCell) => true
       | _ => false)


fun findIndexLength (opt, data) =
    (case (opt, data)
      of (ExtractDofs, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDof, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDofsSeq, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractIndices, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractIndex, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDofs, FT.Space(s)) => FT.spaceDim s
       | (ExtractIndex, FT.Space(s)) => FT.spaceDim s
       | _ => raise Fail "impossible optiion and data combination"
    (* end case *) )


fun findTargetShape (opt, data) =
     (case (opt, data)
       of (ExtractDofs, FT.Mesh(m)) => FT.dataShapeOf data
	| (ExtractDofsSeq, FT.Mesh(m)) => FT.dataShapeOf data
	| (ExtractDofs, FT.Space(s)) => FT.dataShapeOf data
	| (ExtractDofsSeq, FT.Space(s)) => FT.dataShapeOf data
	| (ExtractDof, FT.Mesh(m)) => [FT.meshDim m]
	| (ExtractDof, FT.Space(s)) => FT.dataRangeShapeOf data
	| _ => raise Fail "impossible optiion and data combination"
    (* end case *) )


fun splitDataOpt (opt, data) =
    (case (opt, data)
      of (ExtractDofs, FT.Mesh(m)) => SOME((ExtractIndices, data), (ExtractDofsSeq, data))
       | _ => raise Fail "impossible optiion and data combination"
    (* end case *))

      
fun lowerDataOpt (v,d) = (case v
			   of ExtractIndices => SOME(ExtractIndex, d)
			    | ExtractDofsSeq => SOME(ExtractDof, d)
			    | _ => NONE
			 (* end case *))


      
end
