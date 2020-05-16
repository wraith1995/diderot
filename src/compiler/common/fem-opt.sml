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
			      | ExtractDofs (*PST->High*)
			      | ExtractIndices |  ExtractDofsSeq (*Mid*)
			      | NumCell  | ExtractIndex | ExtractDof (* LOW->*)
			      | StartCell
			      | CellIndex (* primitize, cells*)
			      (*mesh pos operations:*)
			      | RefPos 
			      | RefBuild | WorldBuild of Atom.atom option * Stamp.t option | AllBuild
			      | NewWorld
			      | InvalidBuild | NearbyCellQuery of Atom.atom
			      | InvalidBuildBoundary | CellConnectivity | CellFaceCell
			      | InsideInsert of Atom.atom | PosEntryFacet
			      | CellData of Atom.atom

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

	   datatype femField =  Transform | RefField | InvTransform
			
	   val fieldString : femField -> string

	   datatype stages = PST | AST | SIMPLE | HIGH | MID | LOW | TREE | CXX
	   (*Range of allowed stages*)
	   val stageRange : femOption -> stages * stages

					   
									     
	  end = struct
datatype femOpts = Cells | RefCell
		 | ExtractDofs
		 | ExtractIndices |  ExtractDofsSeq
		 | NumCell  | StartCell | ExtractIndex | ExtractDof (* primitize all*)
		 | CellIndex (* primitize, cells*)
		 (*mesh pos operations:*)
		 | RefPos
		 | RefBuild | WorldBuild of Atom.atom option * Stamp.t option | AllBuild | NewWorld (*function to do inverse, but not needed*)
		 | InvalidBuild | InvalidBuildBoundary  (*internal use only*)
		 | NearbyCellQuery of Atom.atom
		 | CellConnectivity | CellFaceCell
		 | InsideInsert of Atom.atom | PosEntryFacet
		 | CellData of Atom.atom

		     
datatype femField =  Transform | RefField | InvTransform 


fun fieldString (Transform) = "Transform"
  | fieldString (RefField) = "refField"
  | fieldString (InvTransform) = "InvTransform"

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
       | RefCell => "RefCell"
       | RefPos => "RefPos"
       | RefBuild => "RefBuild"
       | WorldBuild r => "WorldBuild(" ^ ((mapAS Atom.toString Stamp.toString "Error") r) ^ ")"
       | InvalidBuild => "InvalidBuild"
       | InvalidBuildBoundary => "InvalidBuildBoundary"
       | AllBuild => "AllBuild"
       | NewWorld => "NewWorld"
       | NearbyCellQuery(a) => "NearbyCell(File="^(Atom.toString a)^")"
       | CellConnectivity => "CellConnectivity"
       | CellFaceCell => "CellFaceCell"
       | InsideInsert(a) => "InsideInsert(File="^(Atom.toString a) ^")"
       | StartCell => "StartCell"
       | PosEntryFacet => "PosEntryFacet"
       | (CellData(a)) => "CellData("^(Atom.toString a)^")"
    (* end case*))

fun toString (v, d) = toStringOpt(v) ^ "(" ^ (FT.toString d) ^ ")"
         

fun arity (NumCell) = 1
  | arity (StartCell) = 1
  | arity (Cells) = 1
  | arity (CellIndex) = 1
  | arity (ExtractDofs) = 2
  | arity (ExtractIndices) = 2
  | arity (ExtractDofsSeq) = 2
  | arity (ExtractIndex) = 2
  | arity (ExtractDof) = 2
  | arity (RefCell) = 1
  | arity (RefPos) = 1
  | arity (RefBuild) = 3
  | arity (WorldBuild _) = 2
  | arity (InvalidBuild) = 1
  | arity (InvalidBuildBoundary) = 2
  | arity (AllBuild) = 4
  | arity (NewWorld) = 2
  | arity (NearbyCellQuery(_)) = 2
  | arity (CellConnectivity) = 2
  | arity (CellFaceCell) = 3
  | arity (InsideInsert(_)) = 2
  | arity (PosEntryFacet) = 1
  | arity (CellData(_)) = 2


fun hash (NumCell, d) = 0w1 + FT.hash d
  | hash (Cells, d) = 0w2 + FT.hash d
  | hash (CellIndex, d) = 0w3 + FT.hash d
  | hash (ExtractIndex, d) = 0w5 + FT.hash d
  | hash (ExtractDof, d) = 0w7 + FT.hash d
  | hash (ExtractDofs, d) = 0w11 + FT.hash d
  | hash (ExtractIndices, d) = 0w13 + FT.hash d
  | hash (ExtractDofsSeq, d) = 0w17 + FT.hash d
  | hash (RefCell, d) = 0w23 + FT.hash d
  | hash (RefPos, d) = 0w31 + FT.hash d
  | hash (RefBuild, d) = 0w37 + FT.hash d
  | hash (InvalidBuild, d) = 0w53 + FT.hash d
  | hash (WorldBuild r, d) = 0w41 + FT.hash d + 0w43 * (mapAS Atom.hash Stamp.hash 0w47 r)
  | hash (AllBuild, d) = 0w73 + FT.hash d
  | hash (NewWorld, d) = 0w79 + FT.hash d
  | hash (NearbyCellQuery(a), d) = 0w83 + (Atom.hash a) * 0w89 + FT.hash d
  | hash (CellConnectivity, d) = 0w89 + FT.hash d
  | hash (InvalidBuildBoundary, d) = 0w97 + FT.hash d
  | hash (CellFaceCell, d) = 0w101 + FT.hash d
  | hash (InsideInsert(a),d) = 0w103 + (Atom.hash a) * 0w107 +  FT.hash d
  | hash (StartCell, d) = 0w107 + FT.hash d
  | hash (PosEntryFacet, d) = 0w109 + FT.hash d
  | hash (CellData(a), d) = 0w113 + 0w127 * (Atom.hash a) + FT.hash d

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
       | (RefCell,RefCell) => true
       | (RefPos, RefPos) => true
       | (RefBuild, RefBuild) => true
       | (InvalidBuild, InvalidBuild) => true
       | (WorldBuild r1, WorldBuild r2) => sameR(r1, r2)
       | (AllBuild, AllBuild) => true
       | (NewWorld, NewWorld) => true
       | (NearbyCellQuery(a1), NearbyCellQuery(a2)) => Atom.same(a1, a2)
       | (CellConnectivity, CellConnectivity) => true
       | (InvalidBuildBoundary,InvalidBuildBoundary) => true
       | (CellFaceCell, CellFaceCell) => true
       | (InsideInsert(a1), InsideInsert(a2)) => Atom.same(a1,a2)
       | (StartCell, StartCell) => true
       | (PosEntryFacet, PosEntryFacet) => true
       | (CellData(a1), CellData(a2)) => Atom.same(a1, a2)
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


datatype stages = PST | AST | SIMPLE | HIGH | MID | LOW | TREE | CXX
							     
fun stageRange (Cells, _) = (PST, SIMPLE)
  | stageRange (RefCell, _) = (PST, AST)
  | stageRange (NumCell, _) = (PST, CXX)
  | stageRange (StartCell, _) = (AST, CXX)
  | stageRange (CellIndex, _) = (AST, CXX)
  | stageRange (CellData _, _) = (AST, CXX)
  | stageRange (ExtractDofs, _) = (PST, HIGH)
  (*That these two exist at MID is a bit confusing for cache optimization*)
  | stageRange (ExtractIndices, _) = (MID, MID)
  | stageRange (ExtractDofsSeq, _) = (MID, MID)
  | stageRange (ExtractIndex, _) = (LOW, CXX)
  | stageRange (ExtractDof, _) = (LOW, CXX)
  | stageRange (RefPos, _) = (PST, CXX)
  | stageRange _ = raise Fail "NYI"
      
end
