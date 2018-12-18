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
	   type space
	   type func

	   datatype femType = Mesh of mesh
			    | Space of space
			    | Func of func
			    | RefCell of mesh
			    | MeshCell of mesh
			    | FuncCell of func
			    | MeshPos of mesh

(*valid input apply - valid output apply -> use these for all the damn utils*)
	   val nameOf : femType -> Atom.atom
	   val meshOf : femType -> femType
	   val spaceOf : femType -> femType
	   val underlyingDim : femType -> int
	   val meshMapDim : mesh -> int
	   val same : femType * femType -> bool
	   val validInput : femType -> bool
	   val extractMesh : femType -> mesh option
	   val extractSpace : femType -> space option
	   val femPP : femType -> string
	   val hash : femType -> word
	   val defaultSpace : mesh * Atom.atom -> femType   
	   val mkMesh : int * int * Atom.atom -> femType
	   val mkSpace : int * int list * mesh * Atom.atom -> femType
	   val mkFunc : space * Atom.atom -> femType

	  end = struct

structure BD = BasisData

datatype mesh = Mesh' of {
	  dim : int,
	  mappingDim : int,
	  basis : BasisData.t list,
	  name : Atom.atom
	 }

fun meshDim (Mesh'({dim, mappingDim, basis, name})) = dim
fun meshMapDim (Mesh'{dim, mappingDim, basis, name}) = mappingDim
fun meshBasis (Mesh'{dim, mappingDim, basis, name}) = basis
datatype space = Space' of {
	  dim : int,
	  shape : int list,
	  basis : BasisData.t list,
	  mesh : mesh,
	  name : Atom.atom
	 }

fun spaceDim (Space'{dim, shape, basis, mesh, name}) = dim
fun spaceShape (Space'{dim, shape, basis, mesh, name}) = shape
fun spaceBasis (Space'{dim, shape, basis, mesh, name}) = basis
fun spaceMesh (Space'{dim, shape, basis, mesh, name}) = mesh


datatype func = Func' of {
	  space : space,
	  name : Atom.atom
	 }
fun funcSpace (Func'{space, name}) = space


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
       | _ => raise Fail "Asked for the name of a non data fem type" (*doesn't make sense*)
    (* end case*)
    )

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
      of Mesh((Mesh'{dim, mappingDim, basis, name}))
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
fun femPP ty =
    (case ty
      of Mesh(_) => "Mesh"
       | Space(_) => "Space"
       | Func(_) => "FemFunc"
       | RefCell(_) => "RefCell"
       | MeshCell(_) => "MeshCell"
       | FuncCell(_) => "FuncCell"
       | MeshPos(_) => "MeshPos"
    (* end case*))

fun extractMesh ty =
    (case ty
      of Mesh(m) => SOME(m)
       | _ => NONE
    (*end case*))
fun extractSpace ty =
    (case ty
      of Space(s) => SOME(s)
       | _ => NONE)

(* fem types that are valid inputs.*)
fun isInputFemType ty =
    (case ty
      of RefCell(_) => false
       | MeshPos(_) => false
       | _ => true)
			  
fun mkMesh(dim, mDim, name) = Mesh (Mesh' {dim = dim, mappingDim = mDim, basis = [BD.empty], name = name})


fun defaultSpace (m as Mesh'{dim,...}, name) = Space (Space' {dim = dim, shape = [], basis = [], mesh = m, name = name})

fun mkSpace(dim, shape, mesh, name) = Space(Space' {dim = dim, shape = shape, basis=[BD.empty], mesh = mesh, name = name})

fun mkFunc (space, name) = Func(Func' {space = space, name = name})

end
		  
