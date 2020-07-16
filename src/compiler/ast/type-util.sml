(* type-util.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure TypeUtil : sig

  (* constructor for building a tensor type of known order, but unknown
   * dimensions.
   *)
    val mkTensorTy : int -> Types.ty
    val mkTensorTy' : int * int -> Types.ty

  (* constructor for building a fixed-size sequence type of unknown size *)
    val mkSequenceTy : Types.ty -> Types.ty

  (* function to compute the slice of a tensor type based on a boolean
   * mask.  The value true in the mask means that the corresponding
   * dimension is being indexed, while false means that it is being
   * copied.
   *)
    val slice : Types.ty * bool list -> Types.ty

(* FIXME: make terminology be consistent between documentation and implementation w.r.t.
 * the kinds of types (e.g., concrete type, value type, ...)
 *)
  (* returns true if the type is a value type; i.e., a basic value type (bool, int,
   * string, or tensor), or a sequence of values.
   *)
    val isValueType : Types.ty -> bool

  (* returns true if the type is a value type, or a strand type, or a sequence of such types *)
    val isValueOrStrandType : Types.ty -> bool

  (* return true if the type is an image type *)
    val isImageType : Types.ty -> bool

  (* return true if the type is T_Error *)
    val isErrorType : Types.ty -> bool

    (* returns true if this is a fem type that we can input*)
    val isInputFemType : Types.ty -> bool
    val extractFemType : Types.ty -> (FemData.femType * Atom.atom option) option
    val femDatas : Types.ty -> (FemData.femType * Atom.atom option) list
    val hasFem : Types.ty -> bool
    val normalilzeFemInputTy : Types.ty -> Types.ty * Types.ty

  (* return the range (return type) of a function type *)
    val rngOf : Types.ty -> Types.ty

  (* prune out instantiated meta variables from a type.  We also normalize
   * tensor shapes (i.e., remove 1s).
   *)
    val prune : Types.ty -> Types.ty
    val pruneDiff : Types.diff -> Types.diff
    val pruneShape : Types.shape -> Types.shape
    val pruneInterval : Types.interval -> Types.interval
    val pruneInterval' : (Types.shape * Types.shape -> bool) -> Types.interval -> Types.interval
    val pruneDim : Types.dim -> Types.dim

  (* prune the head of a type *)
    val pruneHead : Types.ty -> Types.ty

  (* resolve meta variables to their instantiations (or else variable) *)
    val resolve : Types.ty_var -> Types.ty
    val resolveDiff : Types.diff_var -> Types.diff
    val resolveShape : Types.shape_var -> Types.shape
    val resolveDim : Types.dim_var -> Types.dim
    val resolveInterval : Types.interval_var -> Types.interval
    val resolveVar : Types.meta_var -> Types.var_bind
					 

  (* equality testing *)
    (* val sameDim : Types.dim * Types.dim -> bool *)

  (* string representations of types, etc *)
    val toString : Types.ty -> string
    val diffToString : Types.diff -> string
    val shapeToString : Types.shape -> string
    val dimToString : Types.dim -> string
    val intervalToString : Types.interval -> string

  (* convert to fully resolved monomorphic forms *)
    val monoDim : Types.dim -> int
    val monoShape : Types.shape -> int list
    val monoDiff : Types.diff -> int option
    val monoInterval : Types.interval -> int

  (* instantiate a type scheme, returning the argument meta variables and the resulting type.
   * Note that we assume that the scheme is closed.
   *)
    val instantiate : Types.scheme -> (Types.meta_var list * Types.ty)

  end = struct

    structure Ty = Types
    structure MV = MetaVar

  (* constructor for building a tensor type of known order, but unknown
   * dimensions.
   *)
    fun mkTensorTy order =
          Ty.T_Tensor(
           Ty.Shape(List.tabulate(order, fn _ => Ty.DimVar(MetaVar.newDimVar()))),
	   Ty.AddVar(MetaVar.newIntervalVar (), 0))
    fun mkTensorTy' (order, interval) = if interval < 0 then raise Fail "impossible"
				    else
          Ty.T_Tensor(
           Ty.Shape(List.tabulate(order, fn _ => Ty.DimVar(MetaVar.newDimVar()))),
	   Ty.AddVar(MetaVar.newIntervalVar (), 0))


    fun mkSequenceTy ty = Ty.T_Sequence(ty, SOME(Ty.DimVar(MetaVar.newDimVar())))

  (* prune out instantiated meta variables from a type.  We also normalize
   * tensor dimensions (i.e., remove 1s).
   *)
    fun prune ty = (case ty
           of (ty as Ty.T_Var(Ty.TV{bind, ...})) => (case !bind
                 of NONE => ty
                  | SOME ty => prune ty
						    (* end case *))
	    | Ty.T_Named(name, ty') =>  Ty.T_Named(name, prune ty')
            | Ty.T_Sequence(ty, NONE) => Ty.T_Sequence(prune ty, NONE)
            | Ty.T_Sequence(ty, SOME dim) => Ty.T_Sequence(prune ty, SOME(pruneDim dim))
	    | Ty.T_Tuple(tys) => Ty.T_Tuple(List.map prune tys)
            | (Ty.T_Kernel diff) => Ty.T_Kernel(pruneDiff diff)
            | (Ty.T_Tensor (shape, iv)) => Ty.T_Tensor(pruneShape shape, pruneInterval iv)
            | (Ty.T_Image{dim, shape}) => Ty.T_Image{
                  dim = pruneDim dim,
                  shape = pruneShape shape
                }
            | (Ty.T_Field{diff, dim, shape}) => Ty.T_Field{
                  diff = pruneDiff diff,
                  dim = pruneDim dim,
                  shape = pruneShape shape
                }
            | (Ty.T_Fun(tys1, ty2)) => Ty.T_Fun(List.map prune tys1, prune ty2)
            | ty => ty
          (* end case *))

    and pruneDiff (Ty.DiffVar(Ty.DfV{bind=ref(SOME diff), ...}, i)) = (
          case pruneDiff diff
           of Ty.DiffVar(dv, i') => Ty.DiffVar(dv, i+i')
            | Ty.DiffConst NONE => Ty.DiffConst NONE
            | Ty.DiffConst(SOME i') => Ty.DiffConst(SOME(i+i'))
          (* end case *))
      | pruneDiff diff = diff

    and pruneDim dim =
	(case dim
	  of Ty.DimVar(Ty.DV{bind=ref(SOME dim), ...}) => pruneDim dim
	   | Ty.IntervalDim(iv) =>
	     let
	      val Ty.IV{bind, id} = iv
	     in
	      (case !bind
		of NONE => dim
		 | SOME(iv') =>
		   (case pruneInterval iv'
		     of Ty.IC j => Ty.DimConst j
		      | iv'' => (bind := SOME(iv''); dim)))
	     end
           | dim => dim
        (* end case *))

    and filterDim dim = (case pruneDim dim
           of Ty.DimConst 1 => NONE
            | dim => SOME dim
          (* end case *))

    and pruneShape shape = (case shape
           of Ty.Shape dd => Ty.Shape(List.mapPartial filterDim dd)
            | Ty.ShapeVar(Ty.SV{bind=ref(SOME shape), ...}) => pruneShape shape
            | Ty.ShapeExt(shape, dim) =>
	      (case filterDim dim
                of SOME dim => Ty.shapeExt(pruneShape shape, dim)
                 | NONE => pruneShape shape
	      (* end case *))
            | Ty.ShapeExt'(dim, shape) =>
	      (case filterDim dim
                of SOME dim => Ty.ShapeExt'(dim, pruneShape shape)
                 | NONE => pruneShape shape
              (* end case *))
            | _ => shape
			   (* end case *))
			     
    and pruneInterval interval =
	(case interval
	  of Ty.IC j => Ty.IC j
	   | Ty.MaxVar(iv1, iv2) =>
	     (case (pruneInterval iv1, pruneInterval iv2)
	       of (a as Ty.IC(iv1'), b as Ty.IC(iv2')) => 
		  if iv1' = 0 andalso iv2' = 0
		  then Ty.IC 0
		  else if (iv1' = 0) andalso (iv2' > 0)
		  then Ty.IC iv2'
		  else if (iv2' = 0) andalso (iv1' > 0)
		  then Ty.IC iv1'
		  else if (iv1' > 0) andalso (iv1' = iv2')
		  then Ty.IC iv1'
		  else Ty.MaxVar(a, b)
		| (Ty.IC 0, b) => b
		| (a, Ty.IC 0) => a
		| (a as Ty.AddVar(iv1,r), b as Ty.AddVar(iv2, r')) =>
		  if MV.sameIntervalVar(iv1, iv2) andalso r=r'
		  then Ty.AddVar(iv1,r)
		  else Ty.MaxVar (a,b)
		| (iv1', iv2') => Ty.MaxVar(iv1', iv2')
	     (* end case *))
	   | Ty.AddVar(Ty.IV{id=id, bind=ref(SOME(iv))}, j) =>
	     (case pruneInterval iv
	       of Ty.IC(iv') => Ty.IC(iv' + j)
		| Ty.AddVar(iv', k) => Ty.AddVar(iv', j+k)
		| iv' => interval
	     (* end case *))
	   | Ty.ShapeSize(shape) =>
	     (case pruneShape shape
	       of Ty.Shape ds =>
		  let
		   val size = List.foldr (fn (Ty.DimConst j, k) => j *k | _ => 0) 1 ds
					 
		  in
		   if size = 0
		   then Ty.ShapeSize (Ty.Shape ds)
		   else Ty.IC (size + 2)
		  end
		| shape' => Ty.ShapeSize(shape'))
	   | _ => interval
	(* end case *))
    and pruneInterval' equalShape interval  =
	(case pruneInterval interval
	  of Ty.MaxVar(Ty.ShapeSize(shape1), Ty.ShapeSize(shape2)) =>
	     if equalShape(shape1, shape2)
	     then Ty.ShapeSize(shape1)
	     else interval
	   | _ => interval
	(* end case *))

  (* resolve meta variables to their instantiations (or else variable) *)
    fun resolve (tv as Ty.TV{bind, ...}) = (case !bind
           of NONE => Ty.T_Var tv
            | SOME ty => prune ty
          (* end case *))

    fun resolveDiff (dv as Ty.DfV{bind, ...}) = (case !bind
           of NONE => Ty.DiffVar(dv, 0)
            | SOME diff => pruneDiff diff
          (* end case *))

    fun resolveInterval (iv as Ty.IV{bind, ...}) =
	(case !bind
	  of NONE => Ty.AddVar(iv, 0)
	   | SOME iv => pruneInterval iv
	(* end case*))

    fun resolveShape (sv as Ty.SV{bind, ...}) = (case !bind
           of NONE => Ty.ShapeVar sv
            | SOME shape => pruneShape shape
          (* end case *))

    fun resolveDim (dv as Ty.DV{bind, ...}) = (case !bind
           of NONE => Ty.DimVar dv
            | SOME dim => pruneDim dim
          (* end case *))

    fun resolveVar (Ty.TYPE tv) = Ty.TYPE(resolve tv)
      | resolveVar (Ty.DIFF dv) = Ty.DIFF(resolveDiff dv)
      | resolveVar (Ty.SHAPE sv) = Ty.SHAPE(resolveShape sv)
      | resolveVar (Ty.DIM d) = Ty.DIM(resolveDim d)
      | resolveVar (Ty.INTERVAL iv) = Ty.INTERVAL(resolveInterval iv)

  (* prune the head of a type *)
    fun pruneHead ty = let
          fun prune' (ty as Ty.T_Var(Ty.TV{bind, ...})) = (case !bind
                 of NONE => ty
                  | SOME ty => prune' ty
                (* end case *))
            | prune' (Ty.T_Sequence(ty, NONE)) = Ty.T_Sequence(ty, NONE)
            | prune' (Ty.T_Sequence(ty, SOME dim)) = Ty.T_Sequence(ty, SOME(pruneDim dim))
	    | prune' (Ty.T_Tuple(tys)) = Ty.T_Tuple(List.map prune' tys)
            | prune' (Ty.T_Kernel diff) = Ty.T_Kernel(pruneDiff diff)
            | prune' (Ty.T_Tensor (shape, interval)) = Ty.T_Tensor(pruneShape shape, pruneInterval interval)
            | prune' (Ty.T_Image{dim, shape}) = Ty.T_Image{
                  dim = pruneDim dim,
                  shape = pruneShape shape
                }
            | prune' (Ty.T_Field{diff, dim, shape}) = Ty.T_Field{
                  diff = pruneDiff diff,
                  dim = pruneDim dim,
                  shape = pruneShape shape
              }
	    | prune' (Ty.T_Named(name, def)) = Ty.T_Named(name,  prune' def)
            | prune' ty = ty
          in
            prune' ty
          end

  (* helper function for isValueType and isValueOrStrandType; checks for fixed-size types
   * inside a dynamic sequence (i.e., no nested dynamic sequences).
   *)
    fun isFixedSize (allowStrand, ty) = (case ty
           of Ty.T_Bool => true
            | Ty.T_Int => true
            | Ty.T_String => true
            | Ty.T_Sequence(ty, SOME _) => isFixedSize (allowStrand, ty)
            | Ty.T_Strand _ => allowStrand
            | Ty.T_Tensor _ => true
	    | Ty.T_Tuple(tys) => List.all (fn x => isFixedSize(allowStrand, x)) tys
            | Ty.T_Error => true
	    | Ty.T_Named(_, def) => isFixedSize (allowStrand, def)
	    | Ty.T_Fem(_,_) => true
            | _ => false
          (* end case *))

  (* returns true if the type is a value type; i.e., a basic value type (bool, int,
   * string, or tensor), or a sequence of values.
   *)
    fun isValueType ty = (case prune ty
           of Ty.T_Bool => true
            | Ty.T_Int => true
            | Ty.T_String => true
            | Ty.T_Sequence(ty, SOME _) => isValueType ty
            | Ty.T_Sequence(ty, NONE) => isFixedSize (false, ty)
            | Ty.T_Tensor _ => true
	    | Ty.T_Tuple(tys) => (List.all isValueType tys)
			 (*andalso (List.all (fn x => isFixedSize (false, x)) tys)*)
            | Ty.T_Error => true
	    | Ty.T_Named(_, ty') => isValueType ty'
	    | Ty.T_Fem(data,_) => FemData.isValue data
            | _ => false
          (* end case *))

  (* returns true if the type is a value type, or a strand type, or a sequence of such types *)
    fun isValueOrStrandType ty = (case prune ty
           of Ty.T_Bool => true
            | Ty.T_Int => true
            | Ty.T_String => true
            | Ty.T_Sequence(ty, SOME _) => isValueOrStrandType ty
            | Ty.T_Sequence(ty, NONE) => isFixedSize (true, ty)
	    | Ty.T_Tuple(tys) => List.all isValueOrStrandType tys
            | Ty.T_Strand _ => true
	    | Ty.T_Fem(data, _) => FemData.isValue data
	    | Ty.T_Named(_, ty') => isValueOrStrandType ty'
            | Ty.T_Tensor _ => true
            | Ty.T_Error => true
            | _ => false
          (* end case *))

  (* returns true if the type is an ImageTy *)
    fun isImageType ty = (case prune ty
           of Ty.T_Image _ => true
            | Ty.T_Error => true
            | _ => false
          (* end case *))

    fun isErrorType ty = (case prune ty
			   of Ty.T_Error => true
			    | _ => false
			 (* end case *))

    fun extractFemType ty =
	(case ty
	  of Ty.T_Fem(v) => SOME(v)
	   | Ty.T_Named(_, ty') => extractFemType ty'
	   | _ => NONE
	(* end case *))

    fun femDatas ty =
	(case prune ty
	  of Ty.T_Sequence(ty', _) => femDatas ty'
	   | Ty.T_Tuple(tys) => List.concatMap femDatas tys
	   | Ty.T_Named(_, ty') => femDatas ty'
	   | Ty.T_Fem(data) => [data]
	   | _ => []
	(*end case*))
    fun hasFem ty =
	(case prune ty
	  of Ty.T_Sequence(ty', _) => hasFem ty'
	   | Ty.T_Tuple tys => List.exists hasFem tys
	   | Ty.T_Named(_, ty') => hasFem ty'
	   | Ty.T_Fem(data) => true
	   | _ => false)
    fun p(a) = (a,a)
    fun normalilzeFemInputTy (Ty.T_Var(_)) = raise Fail "failed type clean"
      | normalilzeFemInputTy (t as Ty.T_Bool) = p(t)
      | normalilzeFemInputTy (t as Ty.T_Int) = p(t)
      | normalilzeFemInputTy (t as Ty.T_String) = p(t)
      | normalilzeFemInputTy (t as Ty.T_Sequence(ty, opt)) = let val (l,r) = normalilzeFemInputTy(ty)
							     in (Ty.T_Sequence(l, opt), Ty.T_Sequence(r, opt))
							     end
      | normalilzeFemInputTy (t as Ty.T_Tuple(tys)) = let val lrs = (List.map normalilzeFemInputTy tys)
							  val (lefts, rights) = ListPair.unzip lrs
						      in
						       (Ty.T_Tuple(lefts), Ty.T_Tuple(rights))
						      end
      | normalilzeFemInputTy (Ty.T_Named(_, ty)) = normalilzeFemInputTy(ty)
      | normalilzeFemInputTy (Ty.T_Kernel _) = raise Fail "invalid io type"
      | normalilzeFemInputTy (t as Ty.T_Tensor _) = p(t)
      | normalilzeFemInputTy (Ty.T_Image _) = raise Fail "invalid io type"
      | normalilzeFemInputTy (Ty.T_Field _) = raise Fail "invalid io type"
      | normalilzeFemInputTy (Ty.T_Fun _) = raise Fail "invalid io type"
      | normalilzeFemInputTy (t as Ty.T_Error) = p(t)
      | normalilzeFemInputTy (t as Ty.T_Fem(data, _)) =
	(case data
	  of FemData.Mesh _ => raise Fail "covered seperately"
	   | FemData.Space _ => raise Fail "covered seperately"
	   | FemData.Func _ => raise Fail "covered seperately"
	   | FemData.RefCell _ => raise Fail "invalid io type"
	   | FemData.MeshCell _ => (Ty.T_Int, t)
	   | FemData.FuncCell _ => (Ty.T_Int, t)
	   | FemData.MeshPos m =>
	     let val ten = Ty.T_Tensor(Ty.Shape([Ty.DimConst(FemData.underlyingDim data)]), Ty.IC(0))
	     in (Ty.T_Tuple([ten, Ty.T_Int, Ty.T_Int]), t)
	     end (*TODO: Order correct with data.*)
	(*end case*))
						  

			   
    (*mesh, func, space must be alone; all other fem types can appear wherever*)
    fun isInputFemType ty =
	(case prune ty
	  of Ty.T_Fem(data, _) => FemData.validInput data
	   | Ty.T_Named(_, ty') => isInputFemType ty'
	   | _ => false
	(*end case*))
    (*Question: Do we need this function?*)
  (* equality testing *)
    (* fun sameDim (Ty.DimConst d1, Ty.DimConst d2) = (d1 = d2) *)
    (*   | sameDim (Ty.DimVar v1, Ty.DimVar v2) = MetaVar.sameDimVar(v1, v2) *)
    (*   | sameDim (Ty.IntervalDim iv1, Ty.IntervalDim iv2) = MetaVar.sameIntervalVar(iv1, iv2) *)
    (*   | sameDim _ = false *)

    fun listToString fmt sep items = String.concatWith sep (List.map fmt items)

    fun diffToString diff = (case pruneDiff diff
           of Ty.DiffConst NONE => ""  (* QUESTION: should we do something else here? *)
            | Ty.DiffConst(SOME n) => Int.toString n
            | Ty.DiffVar(dv, 0) => MV.diffVarToString dv
            | Ty.DiffVar(dv, i) => if i < 0
                then String.concat["(", MV.diffVarToString dv, "-", Int.toString(~i), ")"]
                else String.concat["(", MV.diffVarToString dv, "+", Int.toString i, ")"]
			    (* end case *))


    fun shapeToString shape = (case pruneShape shape
           of Ty.Shape shape => concat["[", listToString dimToString "," shape, "]"]
            | Ty.ShapeVar sv => concat["[", MV.shapeVarToString sv, "]"]
            | Ty.ShapeExt(shape, d) => let
                fun toS (Ty.Shape shape) = (listToString dimToString "," shape) ^ ","
                  | toS (Ty.ShapeVar sv) = MV.shapeVarToString sv ^ ";"
                  | toS (Ty.ShapeExt(shape, d)) = concat[toS shape, dimToString d, ","]
		  | toS (Ty.ShapeExt'(d,shape)) = concat[dimToString d, toS shape, ","]
                in
                  concat["[", toS shape, dimToString d, "]"]
            end
	    | Ty.ShapeExt'(d, shape) => let
             fun toS (Ty.Shape shape) = (listToString dimToString "," shape) ^ ","
               | toS (Ty.ShapeVar sv) = MV.shapeVarToString sv ^ ";"
               | toS (Ty.ShapeExt(shape, d)) = concat[toS shape, dimToString d, ","]
	       | toS (Ty.ShapeExt'(d,shape)) = concat[dimToString d, toS shape, ","]
            in
             concat["[", dimToString d, ",", toS shape, "]"]
            end
          (* end case *))

    and dimToString dim = (case pruneDim dim
           of Ty.DimConst n => Int.toString n
            | Ty.DimVar v => MV.dimVarToString v
	    | Ty.IntervalDim iv => MV.intervalVarToString iv
			  (* end case *))
    and intervalToString iv =
	(case pruneInterval iv
	  of Ty.IC j => if j=0
			then ""
			else if j = 1
			then "interval "
			else if j > 1
			then "affine["^(Int.toString j)^"] " 
			else "invalid["^(Int.toString j)^"] " 
	   | Ty.MaxVar(a,b) => "affineMax(" ^ (intervalToString a) ^ ", " ^ (intervalToString b) ^") "
	   | Ty.AddVar(iv, j) => "intervalAdd(" ^ (MV.intervalVarToString iv) ^ ", " ^ (Int.toString j) ^ ") "
	   | Ty.ShapeSize(shape) => "affine[|" ^ (shapeToString shape) ^"| + 1]"
	(* end case *))

    fun toString ty = (case pruneHead ty
           of Ty.T_Var(Ty.TV{bind=ref(SOME ty), ...}) => toString ty
            | Ty.T_Var tv => MV.tyVarToString tv
            | Ty.T_Bool => "bool"
            | Ty.T_Int => "int"
            | Ty.T_String => "string"
            | Ty.T_Sequence(ty, NONE) => concat[toString ty, "[]"]
            | Ty.T_Sequence(ty, SOME dim) => concat[toString ty, "[", dimToString dim, "]"]
	    | Ty.T_Tuple(tys) => concat ["(", concat(List.map toString tys), ")"]
            | Ty.T_Strand id => Atom.toString id
	    | Ty.T_Named (id, ty') => (Atom.toString id ) ^"( using " ^ toString ty' ^ ")"
            | Ty.T_Kernel(Ty.DiffConst NONE) => raise Fail "unexpected infinite kernel"
            | Ty.T_Kernel diff => "kernel#" ^ diffToString diff
            | Ty.T_Tensor(Ty.Shape[], iv) => intervalToString(iv) ^ "real"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 2], iv) => intervalToString(iv) ^ "vec2"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 3], iv) => intervalToString(iv) ^ "vec3"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 4], iv) => intervalToString(iv) ^ "vec4"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 2, Ty.DimConst 2], iv) => intervalToString(iv) ^ "mat2"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 3, Ty.DimConst 3], iv) => intervalToString(iv) ^ "mat3"
            | Ty.T_Tensor(Ty.Shape[Ty.DimConst 4, Ty.DimConst 4], iv) => intervalToString(iv) ^ "mat4"
            | Ty.T_Tensor (shape, iv) => intervalToString(iv) ^ "tensor" ^ shapeToString shape
            | Ty.T_Image{dim, shape} => concat[
                  "image(", dimToString dim, ")", shapeToString shape
                ]
            | Ty.T_Field{diff=(Ty.DiffConst NONE), dim, shape} => concat[
                  "field", "(", dimToString dim, ")", shapeToString shape
                ]
            | Ty.T_Field{diff, dim, shape} => concat[
                  "field#", diffToString diff, "(", dimToString dim,
                  ")", shapeToString shape
                ]
            | Ty.T_Fun(tys1, ty2) => let
                fun tysToString [] = "()"
                  | tysToString [ty] = toString ty
                  | tysToString tys = String.concat[
                        "(", listToString toString " * " tys, ")"
                      ]
                in
                  String.concat[tysToString tys1, " -> ", toString ty2]
                end
            | Ty.T_Error => "<error-type>"
	    | Ty.T_Fem(data, a) => let
	     val tyDep = (case a
			   of SOME(atom) => Atom.toString atom
			    | NONE => "NONE")
	    in
	     "FemType: "^ (FemData.toString data) ^ " with type var:" ^ tyDep
	    end
          (* end case *))

  (* return the range (return type) of a function type *)
    fun rngOf (Ty.T_Fun(_, ty)) = ty
      | rngOf ty = raise Fail(concat["TypeUtil.rngOf(", toString ty, ")"])

    fun slice (Ty.T_Tensor(Ty.Shape l, iv), mask) = let
          fun f (d, true, dd) = dd
            | f (d, false, dd) = d::dd
          in
            Ty.T_Tensor(Ty.Shape(ListPair.foldr f [] (l, mask)), iv)
          end
      | slice (Ty.T_Field{shape as Ty.Shape l,diff,dim}, mask) = let
          fun f (d, true, dd) = dd
            | f (d, false, dd) = d::dd
          in
            Ty.T_Field{diff=diff, dim=dim, shape= Ty.Shape (ListPair.foldr f [] (l, mask))}
      end
      | slice (ty, _) = raise Fail(concat["slice(", toString ty, ", _)"])

  (* convert to fully resolved monomorphic forms *)
    fun monoDim dim = (case pruneDim dim
           of Ty.DimConst d => d
            | dim => raise Fail(concat["dim ", dimToString dim, " is not constant"])
          (* end case *))

    fun monoShape shp = (case pruneShape shp
           of Ty.Shape shp => List.map monoDim shp
            | shp => raise Fail(concat["shape ", shapeToString shp, " is not constant"])
          (* end case *))

    fun monoDiff diff = (case pruneDiff diff
           of Ty.DiffConst k => k
            | diff => raise Fail(concat["diff ", diffToString diff, " is not constant"])
			(* end case *))

    fun monoInterval iv =
	(case pruneInterval iv
	  of Ty.IC j => j
	   | _ => raise Fail(concat["interval ", intervalToString iv, " is not constant"])
	(* end case*))

  (* instantiate a type scheme, returning the argument meta variables and the resulting type.
   * Note that we assume that the scheme is closed.
   *)
    fun instantiate ([], ty) = ([], ty)
      | instantiate (mvs, ty) = let
          fun instantiateVar (mv, (mvs, env)) = let
                val mv' = MV.copy mv
                in
                  (mv'::mvs, MV.Map.insert(env, mv, mv'))
                end
          val (mvs, env) = List.foldr instantiateVar ([], MV.Map.empty) mvs
          fun iDiff (Ty.DiffVar(k, i)) = (case MV.Map.find(env, Ty.DIFF k)
                 of SOME(Ty.DIFF k) => Ty.DiffVar(k, i)
                  | _ => raise Fail "impossible"
                (* end case *))
            | iDiff diff = diff
          fun iDim (Ty.DimVar dv) = (case MV.Map.find(env, Ty.DIM dv)
				      of SOME(Ty.DIM dv) => Ty.DimVar dv
				       | _ => raise Fail "impossible"
				    (* end case *))
	    | iDim (Ty.IntervalDim iv) = (case MV.Map.find(env, Ty.INTERVAL iv)
					   of SOME(Ty.INTERVAL iv') => Ty.IntervalDim iv'
					    | _ => raise Fail "impossible"
					 (* end case *))
            | iDim dim = dim
          fun iShape (Ty.ShapeVar sv) = (case MV.Map.find(env, Ty.SHAPE sv)
                 of SOME(Ty.SHAPE sv) => Ty.ShapeVar sv
                  | _ => raise Fail "impossible"
                (* end case *))
            | iShape (Ty.ShapeExt(shape, dim)) = Ty.ShapeExt(iShape shape, iDim dim)
	    | iShape (Ty.ShapeExt'(dim, shape)) = Ty.ShapeExt'(iDim dim, iShape shape)
            | iShape (Ty.Shape dims) = Ty.Shape(List.map iDim dims)

	  fun iInterval (Ty.AddVar(iv, j)) =
	      (case MV.Map.find(env, Ty.INTERVAL iv)
		of SOME(Ty.INTERVAL iv') => Ty.AddVar(iv', j)
		 | _ => raise Fail "impossible"
	      (* end case *))
	    | iInterval (Ty.IC(j)) = Ty.IC(j)
	    | iInterval (Ty.MaxVar(iv1, iv2)) = Ty.MaxVar(iInterval iv1, iInterval iv2)
	    | iInterval (Ty.ShapeSize(shape)) = Ty.ShapeSize(iShape shape)
							 
          fun ity (Ty.T_Var tv) = (case MV.Map.find(env, Ty.TYPE tv)
                 of SOME(Ty.TYPE tv) => Ty.T_Var tv
                  | _ => raise Fail "impossible"
                (* end case *))
            | ity Ty.T_Bool = Ty.T_Bool
            | ity Ty.T_Int = Ty.T_Int
            | ity Ty.T_String = Ty.T_String
            | ity (Ty.T_Sequence(ty, NONE)) = Ty.T_Sequence(ity ty, NONE)
            | ity (Ty.T_Sequence(ty, SOME d)) = Ty.T_Sequence(ity ty, SOME(iDim d))
	    | ity (Ty.T_Tuple(tys)) = Ty.T_Tuple(List.map ity tys)
            | ity (ty as Ty.T_Strand _) = ty
            | ity (Ty.T_Kernel k) = Ty.T_Kernel(iDiff k)
            | ity (Ty.T_Tensor (shape, iv)) = Ty.T_Tensor(iShape shape, iInterval iv)
            | ity (Ty.T_Image{dim, shape}) = Ty.T_Image{dim=iDim dim, shape=iShape shape}
            | ity (Ty.T_Field{diff, dim, shape}) =
                Ty.T_Field{diff=iDiff diff, dim=iDim dim, shape=iShape shape}
            | ity (Ty.T_Fun(dom, rng)) = Ty.T_Fun(List.map ity dom, ity rng)
	    | ity (ty as Ty.T_Named(_)) = ty
	    | ity (ty as Ty.T_Fem(_)) = ty
            | ity Ty.T_Error = Ty.T_Error
          in
            (mvs, ity ty)
          end

  end
