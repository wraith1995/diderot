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


	   val extractMesh : femType -> mesh option
	   val extractSpace : femType -> space option
					       
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

datatype space = Space' of {
	  dim : int,
	  shape : int list,
	  basis : BasisData.t list,
	  mesh : mesh
	 }

datatype func = Func' of {
	  space : space
	 }


datatype femType = Mesh of mesh
		 | Space of space
		 | Func of func
		 | RefCell of mesh
		 | MeshCell of mesh
		 | FuncCell of func
		 | MeshPos of mesh

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

fun mkSpace(dim, shape, mesh) = Space(Space' {dim = dim, shape = shape, basis=[BD.empty], mesh = mesh})

fun mkFunc (space) = Func(Func' {space = space})

end
		  
