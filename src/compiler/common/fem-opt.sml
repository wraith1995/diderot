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
datatype femOpts = Cells | RefCell
		 | ExtractDofs
		 | ExtractIndices |  ExtractDofsSeq
		 | NumCell  | ExtractIndex | ExtractDof (* primitize all*)
		 | CellIndex | PromoteCell (* primitize, cells*)
		 (*mesh pos operations:*)
		 | Valid | RefPos | WorldPos of Atom.atom option * Stamp.t option | UWorldPos
		 | RefBuild | WorldBuild of Atom.atom option * Stamp.t option | AllBuild
		 | InvalidBuild | WorldTest

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

	   datatype femField =  Transform | RefField | InvTransform | IdTransform | Field						  
	   val fieldString : femField -> string
									     
	  end = struct
datatype femOpts = Cells | RefCell
		 | ExtractDofs
		 | ExtractIndices |  ExtractDofsSeq
		 | NumCell  | ExtractIndex | ExtractDof (* primitize all*)
		 | CellIndex | PromoteCell (* primitize, cells*)
		 (*mesh pos operations:*)
		 | Valid | RefPos | WorldPos of Atom.atom option * Stamp.t option | UWorldPos
		 | RefBuild | WorldBuild of Atom.atom option * Stamp.t option | AllBuild
		 | InvalidBuild | WorldTest (*internal use only*)
				 

datatype femField =  Transform | RefField | InvTransform | IdTransform | Field


fun fieldString (Transform) = "Transform"
  | fieldString (RefField) = "refField"
  | fieldString (InvTransform) = "InvTransform"
  | fieldString (IdTransform) = "IdTransform"
  | fieldString (Field) = "Field"

type femOption = femOpts * FemData.femType
structure FT = FemData

fun mapAS f g e (a,b) = (case (a,b)
			  of (SOME(a'), NONE) => f a'
			  |  (NONE, SOME(b')) => g b'
			  | _ => e
			(*end case*))
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
       | RefCell => "RefCell"
       | Valid => "Valid"
       | RefPos => "RefPos"
       | WorldPos r => "WorldPos(" ^ (mapAS Atom.toString Stamp.toString "NONE" r) ^ ")"
       | RefBuild => "RefBuild"
       | WorldBuild r => "WorldBuild(" ^ ((mapAS Atom.toString Stamp.toString "Error") r) ^ ")"
       | InvalidBuild => "InvalidBuild"
       | AllBuild => "AllBuild"
       | WorldTest => "WorldTest"
       | UWorldPos => "UWorldPos"
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
  | arity (RefCell) = 1
  | arity (Valid) = 1
  | arity (RefPos) = 1
  | arity (WorldPos _) = 1
  | arity (RefBuild) = 2
  | arity (WorldBuild _) = 2
  | arity (InvalidBuild) = 1
  | arity (WorldTest) = 1
  | arity (UWorldPos) = 1
  | arity (AllBuild) = 4


fun hash (NumCell, d) = 0w1 + FT.hash d
  | hash (Cells, d) = 0w2 + FT.hash d
  | hash (CellIndex, d) = 0w3 + FT.hash d
  | hash (ExtractIndex, d) = 0w5 + FT.hash d
  | hash (ExtractDof, d) = 0w7 + FT.hash d
  | hash (ExtractDofs, d) = 0w11 + FT.hash d
  | hash (ExtractIndices, d) = 0w13 + FT.hash d
  | hash (ExtractDofsSeq, d) = 0w17 + FT.hash d
  | hash (PromoteCell, d) = 0w19 + FT.hash d
  | hash (RefCell, d) = 0w23 + FT.hash d
  | hash (Valid, d) =  0w29 + FT.hash d
  | hash (RefPos, d) = 0w31 + FT.hash d
  | hash (RefBuild, d) = 0w37 + FT.hash d
  | hash (InvalidBuild, d) = 0w53 + FT.hash d
  | hash (WorldBuild r, d) = 0w41 + FT.hash d + 0w43 * (mapAS Atom.hash Stamp.hash 0w47 r)
  | hash (WorldPos r, d) = 0w53 + FT.hash d + 0w59 * (mapAS Atom.hash Stamp.hash 0w61 r)
  | hash (WorldTest, d) = 0w67 + FT.hash d
  | hash (UWorldPos, d) = 0w71 + FT.hash d
  | hash (AllBuild, d) = 0w73 + FT.hash d

fun sameR ((a1,s1), (a2,s2)) = (case (a1, a2)
				 of (SOME(a1'), SOME(a2')) => Atom.same(a1', a2')
				  | (SOME(_), NONE) => false
				  | (NONE, SOME(_)) => false
				  | _ =>  (case (s1, s2)
					    of (SOME(s1'), SOME(s2')) => Stamp.same(s1', s2')
					     | (NONE, NONE) => true
					     | _ => false
					  (*end case*) )
			       (*end case*))
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
       | (RefCell,RefCell) => true
       | (Valid, Valid) => true
       | (RefPos, RefPos) => true
       | (RefBuild, RefBuild) => true
       | (InvalidBuild, InvalidBuild) => true
       | (WorldTest, WorldTest) => true
       | (WorldBuild r1, WorldBuild r2) => sameR(r1, r2)
       | (WorldPos r1, WorldPos r2) => sameR(r1, r2)
       | (UWorldPos, UWorldPos) => true
       | (AllBuild, AllBuild) => true
       | _ => false
    (*end case*))


fun findIndexLength (opt, data) =
    let
     val spaceDimGet = FT.spaceDim o Option.valOf o FT.extractSpace o Option.valOf o FT.dependencyOf
    in
    (case (opt, data)
      of (ExtractDofs, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDof, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDofsSeq, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractIndices, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractIndices, FT.Space(s)) => FT.spaceDim s
       | (ExtractIndex, FT.Mesh(m)) => FT.meshMapDim m
       | (ExtractDofs, FT.Space(s)) => FT.spaceDim s
       | (ExtractIndex, FT.Space(s)) => FT.spaceDim s
       | (ExtractDofsSeq, FT.Func(f)) => spaceDimGet data
       | (ExtractDofs, FT.Func(f)) => spaceDimGet data
       | (ExtractDof, FT.Func(f)) =>  spaceDimGet data
       | _ => raise Fail ("impossible optiion and data combination: " ^ (toString(opt, data)))
    (* end case *) )
    end

fun findTargetShape (opt, data) =
    
     (case (opt, data)
       of (ExtractDofs, FT.Mesh(m)) => FT.dataShapeOf data
	| (ExtractDofsSeq, FT.Mesh(m)) => FT.dataShapeOf data
	| (ExtractDofs, FT.Func(s)) => FT.dataShapeOf data
	| (ExtractDof, FT.Func(s)) => FT.dataRangeShapeOf data
	| (ExtractDofsSeq, FT.Space(s)) => FT.dataShapeOf data
	| (ExtractDof, FT.Mesh(m)) => [FT.meshDim m]
	| (ExtractDof, FT.Space(s)) => FT.dataRangeShapeOf data
	| (ExtractDofsSeq, FT.Func(f)) => FT.dataShapeOf data
	| _ => raise Fail ("impossible option and data combination: "  ^ (toString(opt, data)))
    (* end case *) )


fun splitDataOpt (opt, data) =
    (case (opt, data)
      of (ExtractDofs, FT.Mesh(m)) =>
	 SOME( (ExtractIndices, data), (ExtractDofsSeq, data))
       | (ExtractDofs, FT.Func(f)) =>
	 (SOME( (ExtractIndices, (Option.valOf o FT.dependencyOf) data), (ExtractDofsSeq,data)))
       | _ => NONE
    (* end case *))

      
fun lowerDataOpt (v,d) = (case v
			   of ExtractIndices => SOME(ExtractIndex, d)
			    | ExtractDofsSeq => SOME(ExtractDof, d)
			    | _ => NONE
			 (* end case *))


      
end
