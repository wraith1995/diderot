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
	   val same : femType * femType -> bool
	   val validInput : femType -> bool
	   val extractMesh : femType -> mesh option
	   val extractSpace : femType -> space option
	   val femPP : femType -> string
	   val hash : femType -> word
	   val defaultSpace : mesh -> femType   
	   val mkMesh : int * int -> femType
	   val mkSpace : int * int list * mesh -> femType
	   val mkFunc : space -> femType

	  end = struct

structure BD = BasisData

datatype mesh = Mesh' of {
	  dim : int,
	  mappingDim : int,
	  basis : BasisData.t list
	 }

fun meshDim (Mesh'({dim, mappingDim, basis})) = dim
fun meshMapDim (Mesh'{dim, mappingDim, basis}) = mappingDim
fun meshBasis (Mesh'{dim, mappingDim, basis}) = basis
datatype space = Space' of {
	  dim : int,
	  shape : int list,
	  basis : BasisData.t list,
	  mesh : mesh
	 }

fun spaceDim (Space'{dim, shape, basis, mesh}) = dim
fun spaceShape (Space'{dim, shape, basis, mesh}) = shape
fun spaceBasis (Space'{dim, shape, basis, mesh}) = basis
fun spaceMesh (Space'{dim, shape, basis, mesh}) = mesh


datatype func = Func' of {
	  space : space
	 }
fun funcSpace (Func'{space}) = space


datatype femType = Mesh of mesh
		 | Space of space
		 | Func of func
		 | RefCell of mesh
		 | MeshCell of mesh
		 | FuncCell of func
		 | MeshPos of mesh


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
      of Mesh((Mesh'{dim, mappingDim, basis}))
	 => 0w1 + 0w3 * (Word.fromInt dim) + 0w5 * (Word.fromInt mappingDim) + (List.foldl (fn (d, s) => 0w11 * BD.hash d + s) 0w7 basis)
       | Space(Space'({dim, shape, basis, mesh}))
       	 => 0w13 + 0w17 * (Word.fromInt dim) + (List.foldl (fn (d, s) => 0w23 * BD.hash d + s) 0w19 basis) + (hash (Mesh(mesh)))
       | Func(Func'{space})
       	 => 0w29 + hash (Space(space))
       | RefCell(mesh)
       	 => 0w31 + (hash (Mesh(mesh)))
       | FuncCell(func)
       	 => 0w37 + (hash (Func(func)))
       | MeshPos(mesh)
       	 => 0w41 + (hash (Mesh(mesh)))
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
			  
fun mkMesh(dim, mDim) = Mesh (Mesh' {dim = dim, mappingDim = mDim, basis = [BD.empty]})


fun defaultSpace (m as Mesh'{dim,...}) = Space (Space' {dim = dim, shape = [], basis = [], mesh = m})

fun mkSpace(dim, shape, mesh) = Space(Space' {dim = dim, shape = shape, basis=[BD.empty], mesh = mesh})

fun mkFunc (space) = Func(Func' {space = space})

end
		  
