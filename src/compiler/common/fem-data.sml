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

	   type mesh
	   val meshDim : mesh -> int
		  
	   type space
	   type func

	   datatype femType = Mesh of mesh
			    | Space of space
			    | Func of func
			    | RefCell of mesh
			    | MeshCell of mesh
			    | FuncCell of func
			    | MeshPos of mesh

	   (* Extract the name of a type of fem data as provided by the surface *)
	   val nameOf : femType -> Atom.atom

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
	   val spaceDim : space -> int


	   (* various utilities for type checking, printing, value numbering. *)
	   val same : femType * femType -> bool
	   val toString : femType -> string
	   val hash : femType -> word

	   (* various utilities for figuring out when we can use some sort of fem data*)
	   val validInput : femType -> bool
	   val isValue : femType -> bool
	   val baseFem : femType -> bool


	   (* placeholders for creating fem data *)
	   val defaultSpace : mesh * Atom.atom -> femType   
	   val mkMesh : int * int * int  * BasisData.t list * Atom.atom -> femType
	   val mkSpace : int * int list * mesh * Atom.atom -> femType
	   val mkFunc : space * Atom.atom -> femType

	  end = struct

structure BD = BasisData

datatype mesh = Mesh' of {
	  dim : int,
	  mappingDim : int,
	  degree : int,
	  basis : BasisData.t list,
	  name : Atom.atom
	 }

fun meshDim (Mesh'({dim, mappingDim, basis, name, ...})) = dim
fun meshMapDim (Mesh'{dim, mappingDim, basis, name,...}) = mappingDim
fun meshBasis (Mesh'{dim, mappingDim, basis, name, ...}) = basis
datatype space = Space' of {
	  dim : int,
	  shape : int list,
	  basis : BasisData.t list,
	  mesh : mesh,
	  name : Atom.atom
	 }

fun spaceDim (Space'{dim, shape, basis, mesh, name, ...}) = dim
fun spaceShape (Space'{dim, shape, basis, mesh, name, ...}) = shape
fun spaceBasis (Space'{dim, shape, basis, mesh, name, ...}) = basis
fun spaceMesh (Space'{dim, shape, basis, mesh, name, ...}) = mesh


datatype func = Func' of {
	  space : space,
	  name : Atom.atom
	 }
fun funcSpace (Func'{space, name, ...}) = space



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
       | _ => raise Fail "Asked for the name of ref cell fem type" (*doesn't make sense*)
    (* end case*)
    )

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

    (*end case*))

fun dataRangeShapeOf data =
    (case data
      of Mesh(m) => [meshDim m]
       | Func(Func'{space,...}) => spaceShape space
       | _ => raise Fail "impossible")


fun sameMesh(m1, m2)
    = ((meshDim m1) = meshDim(m2)) andalso (meshMapDim(m1) = meshMapDim(m2)) andalso (ListPair.all BD.same ((meshBasis m1), (meshBasis m2)))
fun sameSpace(s1, s2)
    = (spaceDim(s1) = spaceDim(s2)) andalso (ListPair.all BD.same (spaceBasis s1, spaceBasis s2)) andalso sameMesh((spaceMesh s1), spaceMesh(s2))
      andalso ((List.length (spaceShape s1)) = (List.length (spaceShape s2))) andalso (ListPair.all (fn (x,y) => x=y) ((spaceShape s1), (spaceShape s2)))
fun sameFunc(f1, f2) = sameSpace(funcSpace f1, funcSpace f2)
				
				

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
      of Mesh((Mesh'{dim, degree, mappingDim, basis, name}))
	 => 0w1 + 0w3 * (Word.fromInt dim) + 0w5 * (Word.fromInt mappingDim) + (List.foldl (fn (d, s) => 0w11 * BD.hash d + s) 0w7 basis)
       | Space(Space'({dim, shape, basis, mesh, name}))
       	 => 0w13 + 0w17 * (Word.fromInt dim) + (List.foldl (fn (d, s) => 0w23 * BD.hash d + s) 0w19 basis) + (hash (Mesh(mesh)))
       | Func(Func'{space, name}) => 0w29 + hash (Space(space))
       | RefCell(mesh) => 0w31 + (hash (Mesh(mesh)))
       | FuncCell(func)  => 0w37 + (hash (Func(func)))
       | MeshPos(mesh) => 0w41 + (hash (Mesh(mesh)))
       | MeshCell(mesh) => 0w43 + (hash (Mesh(mesh)))
    (* end case*))
      
      
fun validInput ty =
    (case ty
      of Mesh(_) => true
       | Space(_) => true
       | Func(_) => true
       | _ => false)
fun isValue ty =
    (case ty
      of MeshCell(_) => true
       | FuncCell(_) => true
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
     val _ = print((toString ty) ^ "\n");
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
      
fun mkMesh(dim, mDim, maxDegree, basis, name) = Mesh (Mesh' {dim = dim, mappingDim = mDim, degree = maxDegree, basis = basis, name = name})


fun defaultSpace (m as Mesh'{dim,...}, name) = Space (Space' {dim = dim, shape = [], basis = [], mesh = m, name = name})

fun mkSpace(dim, shape, mesh, name) = Space(Space' {dim = dim, shape = shape, basis=[], mesh = mesh, name = name})

fun mkFunc (space, name) = Func(Func' {space = space, name = name})

end
		  
