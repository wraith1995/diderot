(* types.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *
 * Internal representation of Diderot types.  These are the types produced
 * by the type checker.
 *)

structure Types =
  struct

  (* union type for meta variables and meta-variable bindings *)
    datatype ('ty, 'diff, 'shape, 'dim, 'interval) kind
      = TYPE of 'ty             (* ranges over types *)
      | DIFF of 'diff           (* ranges of differentiation levels (0, 1, ...) *)
      | SHAPE of 'shape         (* ranges over tensor shapes *)
      | DIM of 'dim             (* ranges over dimensions (1, 2, ...) *)
      | INTERVAL of 'interval

    datatype ty
      = T_Var of ty_var         (* meta variable of kind type *)
      | T_Bool
      | T_Int
      | T_String
      | T_Sequence of ty * dim option
      | T_Tuple of ty list
                                (* sequence of values; the dimension is none for *)
                                (* dynamic-length sequences *)
      | T_Strand of Atom.atom   (* named type (i.e., a strand) *)
      | T_Named of Atom.atom * ty    (* other named type *)
    (* convolution kernel; argument is number of levels of differentiation *)
      | T_Kernel of diff
    (* scalars, vectors, matrices, etc.; argument is tensor shape *)
      | T_Tensor of shape * interval
    (* data sets from NRRD *)
      | T_Image of {
          dim : dim,            (* 2D or 3D data set *)
          shape : shape         (* tensor shape; order is length of list *)
        }
    (* continuous field reconstructed from a data set *)
      | T_Field of {
          diff : diff,          (* number of levels of differentiation supported *)
          dim : dim,            (* dimension of domain (2D or 3D field) *)
          shape : shape         (* shape of tensors in range; order is length of list *)
        }
      | T_Fun of ty list * ty
      | T_Error                 (* for when there is a type error *)
      | T_Fem of FemData.femType * Atom.atom option (*fem type and a type name that it uses *)


  (* meta variable of kind type *)
    and ty_var = TV of {
          id : Stamp.t,
          bind : ty option ref
        }

  (* levels of differentiation *)
    and diff
      = DiffConst of int option         (* i (i >= 0); NONE means infinity differentiation *)
      | DiffVar of (diff_var * int)     (* d + i *)

  (* differentiation meta variable *)
    and diff_var = DfV of {
          id : Stamp.t,
          bound : int ref,              (* lower bound of differentiation *)
          bind : diff option ref        (* unification binding *)
        }

  (* tensor shapes *)
    and shape
	= Shape of dim list
      | ShapeVar of shape_var
      | ShapeExt of shape * dim         (* extension of shape (i.e., for D operator) *)

  (* shape meta variable *)
    and shape_var = SV of {
          id : Stamp.t,
          bind : shape option ref       (* unification binding *)
        }

  (* dimensions *)
    and dim
      = DimConst of int                 (* i *)
      | DimVar of dim_var

  (* dimension meta variable *)
    and dim_var = DV of {
          id : Stamp.t,
          bind : dim option ref         (* unification binding *)
        }
    and interval (* interval marker for a tensor: 0 for standard, 1 for interval, d = mid + d -2 errors + error*)
	= IC of int
        | MaxVar of interval_var * interval_var (* for bin ops where same types are added or one type + scalar *)
        | AddVar of interval_var * int (* D+j for various operations *)
    and interval_var = IV of {
	 id: Stamp.t,
	 bind : interval option ref
	}
    type meta_var = (ty_var, diff_var, shape_var, dim_var, interval_var) kind
    type var_bind = (ty, diff, shape, dim, interval) kind

    type scheme = meta_var list * ty

  (* useful types *)
    val realTy' = T_Tensor(Shape[], IC(0))
    val realTy = fn x => T_Tensor(Shape[], x)
    fun vecTy' n = T_Tensor(Shape[DimConst n], IC(0))
    fun matTy' n = T_Tensor(Shape[DimConst n, DimConst n], IC(0))
    fun vecTy n =  fn x => T_Tensor(Shape[DimConst n], x)
    fun matTy n = fn x => T_Tensor(Shape[DimConst n, DimConst n], x)

    val vec2Ty' = vecTy' 2
    val vec3Ty' = vecTy' 3
    val vec4Ty' = vecTy' 4
    val mat2Ty' = matTy' 2
    val mat3Ty' = matTy' 3
    val mat4Ty' = matTy' 4

    fun vecField(diffInt, dimInt) = T_Field {diff=DiffConst(diffInt),
					     dim=DimConst(dimInt),
					     shape=Shape([DimConst(dimInt)])}

  (* smart constructors for building normalized forms *)
    fun shapeExt (Shape dd, d) = Shape(dd @ [d])
      | shapeExt (shape, d) = ShapeExt(shape, d)

  end
  
  
