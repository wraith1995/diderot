(* fem-data.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * Information about meshes, function spaces, and fuctions in function spaces.
 *)

structure FemData : sig
	   datatype refCellType = KSimplex of int | KCube of int | Other of int
	   val fromStr : string * int -> refCellType option
	   datatype refCell = RefCellData of {ty : refCellType,
					      eps : RealLit.t,
					      insideInsert : string option,
					      newtonControl :
					      {contraction : bool, itters: int, newtonTol : RealLit.t, killAfterTol : bool, start : RealLit.t list option}, numFaces : int}
					     
	   type mesh
	   val meshDim : mesh -> int
	   val meshBasis : mesh -> BasisDataArray.t
	   val refCell : mesh -> refCell
	   val meshAccInsert : mesh -> (Atom.atom * bool) option
		  
	   type space
	   val spaceShape : space -> int list
	   val spaceDim : space -> int
	   val spaceBasis : space -> BasisDataArray.t
	   type func
	   val rangeShape : func -> int list

	   datatype femType = Mesh of mesh
			    | Space of space
			    | Func of func
			    | RefCell of mesh
			    | MeshCell of mesh
			    | FuncCell of func
			    | MeshPos of mesh

	   (* Extract the name of a type of fem data as provided by the surface *)
	   val nameOf : femType -> Atom.atom
	   val envNameOf : femType -> Atom.atom
	   val functionNameMake : femType  -> string -> Atom.atom
	   val functionNameMake' : femType  -> string -> Stamp.t -> Atom.atom

	   (* Extract femData that helps to define other fem data *)
	   val meshOf : femType -> femType
	   val spaceOf : femType -> femType
	   val extractMesh : femType -> mesh option
	   val extractSpace : femType -> space option
	   val dependencyOf : femType -> femType option
	   val cellOf : femType -> femType

	   (* Extract dimensional constants*)
				      
	   val underlyingDim : femType -> int
	   val dataShapeOf : femType -> int list
	   val dataRangeShapeOf : femType -> int list
	   val meshMapDim : mesh -> int



	   (* various utilities for type checking, printing, value numbering. *)
	   val same : femType * femType -> bool
	   val toString : femType -> string
	   val hash : femType -> word

	   (* various utilities for figuring out when we can use some sort of fem data*)
	   val validInput : femType -> bool
	   val isValue : femType -> bool
	   val baseFem : femType -> bool


	   (* placeholders for creating fem data *)

	   val mkMesh : int * int * int  * BasisDataArray.t * Atom.atom * refCell * (Atom.atom * bool) option -> femType
	   val mkSpace : int * int list * mesh * BasisDataArray.t * Atom.atom -> femType
	   val mkFunc : space * int list * Atom.atom -> femType

	  end = struct

structure BD = BasisData
datatype refCellType = KSimplex of int | KCube of int | Other of int
fun fromStr(str,dim) = if str = "simplex"
		       then SOME(KSimplex(dim))
		       else if str = "cube"
		       then SOME(KCube(dim))
		       else if str = "other"
		       then SOME(Other(dim))
		       else NONE
			      
datatype refCell = RefCellData of {ty : refCellType,
				   eps : RealLit.t,
				   insideInsert : string option,
				   newtonControl :
				   {contraction : bool, itters: int, newtonTol : RealLit.t, killAfterTol : bool, start : RealLit.t list option}, numFaces : int}

fun getCellType(RefCellData{ty,...}) = ty
fun getCellEps(RefCellData{eps,...}) = eps
				     
fun sameRefCell(a,b) = RealLit.same(getCellEps(a), getCellEps(b)) andalso
		       (case (getCellType a, getCellType b)
			 of (KSimplex(d1), KSimplex(d2)) => d1=d2
			  | (KCube(d1), KCube(d2)) => d1=d2
			  | _ => false
		       (*end case*))
fun hashRefCell(a) = 0w3 * RealLit.hash(getCellEps(a)) + 
		     0w5 * (case (getCellType a)
			     of KSimplex(d1) => 0w3*(Word.fromInt d1)
			      | KCube(d2) => 0w5*(Word.fromInt d2)
			      | Other(d2) => 0w7*(Word.fromInt d2)
			   (*end case*))
(*refcells are `simplex, gen triangle, quad, *)

datatype mesh = Mesh' of {
	  dim : int,
	  mappingDim : int,
	  degree : int, (* default degree; I don't think we use this; should be removed*)
	  basis : BasisDataArray.t,
	  refCell : refCell,
	  name : Atom.atom,
	  insert : (Atom.atom * bool) option
	 }

fun meshDim (Mesh'({dim, mappingDim, basis, name, ...})) = dim
fun meshMapDim (Mesh'{dim, mappingDim, basis, name,...}) = mappingDim
fun meshBasis (Mesh'{dim, mappingDim, basis, name, ...}) = basis
fun meshAccInsert(Mesh'({insert,...})) = insert
fun refCell(Mesh'{refCell,...}) = refCell
				    
datatype space = Space' of {
	  dim : int,
	  shape : int list,
	  basis : BasisDataArray.t,
	  mesh : mesh,
	  name : Atom.atom
	 }

fun spaceDim (Space'{dim, shape, basis, mesh, name, ...}) = dim
fun spaceShape (Space'{dim, shape, basis, mesh, name, ...}) = shape
fun spaceBasis (Space'{dim, shape, basis, mesh, name, ...}) = basis
fun spaceMesh (Space'{dim, shape, basis, mesh, name, ...}) = mesh


datatype func = Func' of {
	  space : space,
	  shape : int list,
	  name : Atom.atom
	 }
fun funcSpace (Func'{space, name, ...}) = space
fun funcShape (Func'{shape, name, ...}) = shape
val rangeShape = funcShape


datatype femType = Mesh of mesh
		 | Space of space
		 | Func of func
		 | RefCell of mesh
		 | MeshCell of mesh
		 | FuncCell of func
		 | MeshPos of mesh


fun meshOf f =
    (case f
      of Mesh(m) => f
       | Space(Space'{mesh,...}) => Mesh(mesh)
       | Func(Func'{space, ...}) => meshOf (Space(space))
       | RefCell(m) => Mesh(m)
       | MeshCell(m) => Mesh(m)
       | FuncCell(s) => meshOf (Func(s))
       | MeshPos(m) => Mesh(m)
    (* end case *))

fun spaceOf f =
    (case f
      of Space(_) => f
       | Func(Func'{space, ...}) => Space(space)
       | FuncCell(Func'{space, ...}) => Space(space)
       | _ => raise Fail "No space."
    (* end case*))

fun nameOf data =
    (case data
      of Mesh(Mesh'{name,...}) => name
       | Space(Space'{name,...}) => name
       | Func(Func'{name,...}) => name
       | MeshCell(Mesh'{name, ...}) => Atom.atom ("mesh_cell_"^(Atom.toString name))
       | FuncCell(Func'{name, ...}) => Atom.atom ("func_cell_"^(Atom.toString name))
       | MeshPos(Mesh'{name, ...}) => Atom.atom ("mesh_pos_" ^(Atom.toString name))
       | RefCell(Mesh'{name,...}) =>  Atom.atom ("ref_cell_" ^(Atom.toString name))
    (* end case*)
    )

fun envNameOf data = 
 (case data
   of MeshCell(Mesh'{name, ...}) => Atom.atom ("call("^(Atom.toString name)^")")
    | FuncCell(Func'{name, ...}) => Atom.atom ("cell("^(Atom.toString name)^")")
    | MeshPos(Mesh'{name, ...}) => Atom.atom ("pos(" ^(Atom.toString name)^")")
    | RefCell(Mesh'{name,...}) =>  Atom.atom ("refCell(" ^(Atom.toString name)^")")
    | _ => raise Fail "Asked for env name on base fem type" (* end case*))
					     
fun functionNameMake data name = Atom.atom ("_" ^ (Atom.toString (nameOf data)) ^ "_" ^ name)
fun functionNameMake' data name stamp = Atom.atom ("_" ^ (Stamp.toString(stamp)) ^ "_" ^ (Atom.toString (nameOf data)) ^ "_" ^ name)
fun cellOf data =
    (case data
      of Mesh(m) => MeshCell(m)
       | Func(f) => FuncCell(f)
       | _ => raise Fail "impossible to turn this femdata into a cell type")
(* fun trueNameOf data = *)
(*     (case data *)
(*       of Mesh(m) => Atom.atom (String.concat(["mesh_", Int.toString (meshDim(m)), "_", Int.toString (meshMapDim(m)),]))) *)

fun underlyingDim data =
    (case data
      of Mesh(m) => meshDim m
       | Space(Space'{mesh,...}) => meshDim mesh
       | Func(Func'{space,...}) => let val Space'{dim,shape,basis,mesh,name} = space in meshDim mesh end
       | RefCell(m) => meshDim m
       | MeshCell(m) => meshDim m
       | FuncCell(f) => underlyingDim (Func(f))
       | MeshPos(m) => meshDim m
    (* end case*)
    )
fun dataShapeOf data =
    (case data
      of Mesh(m) => [meshMapDim m, meshDim m]
       | MeshCell(m) => [meshMapDim m, meshDim m] 
       | MeshPos(m) => [meshMapDim m, meshDim m]
       | Space(Space'{shape,dim,...}) => dim::shape
       | Func(Func'{space,...}) => dataShapeOf (Space(space)) (*TODO: This an error*)
       | FuncCell(f) => dataShapeOf (Func(f))

    (*end case*))

fun dataRangeShapeOf data =
    (case data
      of Mesh(m) => [meshDim m]
       | Func(Func'{shape,...}) => shape
       | _ => raise Fail "impossible")


fun sameMesh(m1, m2)
    = ((meshDim m1) = meshDim(m2))
      andalso (meshMapDim(m1) = meshMapDim(m2))
      andalso (BasisDataArray.same ((meshBasis m1), (meshBasis m2)))
      andalso (sameRefCell(refCell(m1), refCell(m2)))
      andalso (case (meshAccInsert(m1), meshAccInsert(m2))
		of (SOME(a1,b1), SOME(a2,b2)) => Atom.same(a1,a2) andalso b1=b2
		|  (NONE, NONE) => true
		| _ => false
	      (*end case*))

fun sameSpace(s1, s2)
    = (spaceDim(s1) = spaceDim(s2)) andalso (BasisDataArray.same (spaceBasis s1, spaceBasis s2)) andalso sameMesh((spaceMesh s1), spaceMesh(s2))
      andalso ((List.length (spaceShape s1)) = (List.length (spaceShape s2))) andalso (ListPair.all (fn (x,y) => x=y) ((spaceShape s1), (spaceShape s2)))
fun sameFunc(f1, f2) = sameSpace(funcSpace f1, funcSpace f2) andalso funcShape(f1)=funcShape(f2)
				
				

fun same(t1, t2) =
    (case (t1,t2)
      of (Mesh(m1), Mesh(m2)) => sameMesh(m1,m2)
       | (Space(s1), Space(s2)) => sameSpace(s1, s2)
       | (Func(f1), Func(f2)) => sameFunc(f1, f2)
       | (RefCell(m1), RefCell(m2)) => sameMesh(m1,m2)
       | (MeshCell(m1), MeshCell(m2)) => sameMesh(m1,m2)
       | (FuncCell(f1), FuncCell(f2)) => sameFunc(f1, f2)
       | (MeshPos(m1), MeshPos(m2)) => sameMesh(m1,m2)
       | _ => false
    )

fun hash ty =
    (case ty
      of Mesh((Mesh'{dim, degree, mappingDim, basis, name,refCell, insert}))
	 => 0w1 + 0w3 * (Word.fromInt dim) + 0w5 * (Word.fromInt mappingDim) + BasisDataArray.hash basis + 0w7 * hashRefCell(refCell) +
	    (case insert
	      of SOME(a1,b1) => 0w47 * Atom.hash a1 * (if b1
						       then 0w2
						       else 0w3)
	       | _ => 0w0
	    (*end case*)
	    )
       | Space(Space'({dim, shape, basis, mesh, name}))
       	 => 0w13 + 0w17 * (Word.fromInt dim) + BasisDataArray.hash basis + (hash (Mesh(mesh)))
       | Func(Func'{space, name, shape}) => 0w29 + hash (Space(space)) + 0w31 * (List.foldr (fn (x,y) => y  + 0w41 * (Word.fromInt x)) 0w37 (shape)) (*foldr consistnecy of hashing*)
       | RefCell(mesh) => 0w31 + (hash (Mesh(mesh)))
       | FuncCell(func)  => 0w37 + (hash (Func(func)))
       | MeshPos(mesh) => 0w41 + (hash (Mesh(mesh)))
       | MeshCell(mesh) => 0w43 + (hash (Mesh(mesh)))
    (* end case*))
      
      
fun validInput ty =
    (case ty
      of RefCell(_) => false
       | MeshPos(_) => false
       | _ => true)

fun isValue ty =
    (case ty
      of MeshCell(_) => true
       | FuncCell(_) => true
       | MeshPos(_) => true
       | _ => false)

fun toString ty =
    (case ty
      of Mesh(_) => "Mesh"
       | Space(_) => "Space"
       | Func(_) => "FemFunc"
       | RefCell(_) => "RefCell"
       | MeshCell(_) => "MeshCell"
       | FuncCell(_) => "FuncCell"
       | MeshPos(_) => "MeshPos"
    (* end case*))
fun baseFem ty =
    let
     val _ = ()
    in
    (case ty
      of Mesh(_) => true
       | Space(_) => true
       | Func(_) => true
       | _ => false)
    end      

fun extractMesh ty =
    (case ty
      of Mesh(m) => SOME(m)
       | _ => NONE
    (*end case*))
fun extractSpace ty =
    (case ty
      of Space(s) => SOME(s)
       | _ => NONE)

fun dependencyOf ty =
    (case ty
      of Mesh(_) => NONE
       | Space(Space'{mesh, ...}) => SOME(Mesh(mesh))
       | Func(Func'{space,...}) => SOME(Space(space))
       | RefCell(m) => SOME(Mesh(m))
       | MeshCell(m) => SOME(Mesh(m))
       | FuncCell(f) => SOME(Func(f))
       | MeshPos(m) => SOME(Mesh(m))
    (* end case *))

(* fem types that are valid inputs.*)
fun isInputFemType ty =
    (case ty
      of RefCell(_) => false
       | MeshPos(_) => false
       | _ => true)
      
fun mkMesh(dim, mDim, maxDegree, basis, name, refCell, insert) = Mesh (Mesh' {dim = dim, mappingDim = mDim, degree = maxDegree, basis = basis, name = name, refCell=refCell, insert=insert})

fun mkSpace(dim, shape, mesh, basis, name) = Space(Space' {dim = dim, shape = shape, basis= basis, mesh = mesh, name = name})

fun mkFunc (space, shape, name) = Func(Func' {space = space, name = name, shape=shape})

end
		  
