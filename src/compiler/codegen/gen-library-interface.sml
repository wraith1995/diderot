(* gen-library-interface.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *
 * Generate the C header file for the library produced by the compiler.
 *
 * The format of header file is:
 *
 *      HEAD
 *      decls for input variables
 *      BODY
 *      decls for output variables
 *      FOOT
 *)

structure GenLibraryInterface : sig

    val gen : {
            env : CodeGenEnv.t,                 (* target information *)
            subs : (string * string) list,      (* substitutions for fragment expansion *)
            rt : string option,                 (* fragment with extra runtime system hooks *)
            inputs : TreeGlobalVar.t Inputs.input list,
            outputs : OutputUtil.output_info list
          } -> unit

  (* translate an APIType to the C type used to represent it in the library API *)
    val toCType : CodeGenEnv.t * APITypes.t -> CLang.ty

  (* standard names and types used in the library interface *)
    val inputDesc      : TargetSpec.t * string -> string
    val inputGet       : TargetSpec.t * string -> string
    val inputGetSize   : TargetSpec.t * string -> string
    val inputSet       : TargetSpec.t * string -> string
    val inputSetByName : TargetSpec.t * string -> string
    val snapshotGet    : TargetSpec.t * string -> string
    val outputGet      : TargetSpec.t * string -> string
    val worldTy        : TargetSpec.t -> CLang.ty

  end = struct

    structure Ty = APITypes
    structure CL = CLang
    structure Spec = TargetSpec
    structure Env = CodeGenEnv
    structure Out = CodeOutput

    val openOut = Out.openOut {ext = "h", ppDecl = PrintAsC.output}

    val nrrdPtrTy = CL.T_Ptr(CL.T_Named "Nrrd")
    val stringTy = CL.constPtrTy CL.charTy

  (* translate an APIType to the C type used to represent it in the external API *)
    fun toCType (env, ty) = (case ty
           of Ty.IntTy => Env.intTy env
            | Ty.BoolTy => Env.boolTy env
            | Ty.TensorTy[] => Env.realTy env
            | Ty.TensorTy dd => CL.T_Array(Env.realTy env, SOME(List.foldl Int.* 1 dd))
            | Ty.StringTy => CL.constPtrTy CL.charTy
            | Ty.ImageTy(dim, shp) => CL.constPtrTy CL.charTy
	    | Ty.FemData(data) => (case data
				    of FemData.Mesh(_) => CL.voidPtr
				     | FemData.Space(_) => CL.voidPtr
				     | FemData.Func(_) => CL.voidPtr
				     | FemData.RefCell(_) => raise Fail "RefCell should have disappeared"
				     | FemData.MeshCell(_) => raise Fail "undecided"
				     | FemData.FuncCell(_) => raise Fail "undecided"
				     | FemData.MeshPos(_) => raise Fail "undecided")
            | Ty.SeqTy(ty, NONE) => raise Fail "unexpected dynamic SeqTy"
            | Ty.SeqTy(ty, SOME n) => CL.T_Array(toCType(env, ty), SOME n)
			    (* end case *))
    (*copied from gen-inputs.sml -> TODO: fix code duplication between these*)
    fun trTypeAlt (env, ty) = (case ty
			     of Ty.IntTy => Env.intTy env
			      | Ty.BoolTy => Env.boolTy env
			      | Ty.TensorTy[] => Env.realTy env
			      | Ty.TensorTy dd => CL.T_Array(Env.realTy env, SOME(List.foldl Int.* 1 dd))
			      | Ty.StringTy => CL.constPtrTy CL.charTy
			      | Ty.ImageTy(dim, _) => CL.T_Ptr(CL.T_Named "Nrrd")
			      | Ty.SeqTy(ty, NONE) => CL.T_Ptr(CL.T_Named "Nrrd")
			      | Ty.SeqTy(ty, SOME n) => CL.T_Array(trTypeAlt(env, ty), SOME n)
			   (* end case *))

    fun mkSymbol base = let
          fun tr c = if Char.isAlpha c then Char.toUpper c
                else if Char.isDigit c then c
                else #"_"
          in
            String.concat["_", CharVector.map tr base, "_H_"]
          end

  (* standard names and types used in the library interface *)
    val inputDesc      = TargetSpec.qualifyCId' "input_desc"
    val inputGet       = TargetSpec.qualifyCId' "input_get"
    val inputGetSize   = TargetSpec.qualifyCId' "input_get_size"
    val inputSet       = TargetSpec.qualifyCId' "input_set"
    val inputSetByName = TargetSpec.qualifyCId' "input_set_by_name"
    val snapshotGet    = TargetSpec.qualifyCId' "snapshot"
    val outputGet      = TargetSpec.qualifyCId' "output_get"

    fun worldTy spec = CLang.T_Ptr(CLang.T_Named(TargetSpec.qualifyCId "world_t" spec))

    fun gen {env, subs, rt, inputs, outputs} = let
          val spec = Env.target env
          val baseName = OS.Path.joinDirFile{
                  dir = #outDir spec,
                  file = #outBase spec
                }
        (* the world pointer type *)
          val worldPtrTy = worldTy spec

				   
        (* create decls for an input variable *)
          fun mkInputDecls (Inputs.INP{var, name, ty,  desc, init}) = let
                val wrldParam = CL.PARAM([], worldPtrTy, "wrld")
              (* create a description declaration for the input variable *)
                val descDcl = (case desc
                       of NONE => []
                        | SOME desc => [
                              CL.D_Var(["extern"], CL.T_Ptr(CL.T_Named "const char"),
                                [], inputDesc(spec, name), NONE)
                            ]
                      (* end case *))
                val getDcl =  if not(Ty.hasTupleEsq ty)
                      then let
                        val name = inputGet(spec, name)
                      (* convert the input type to a by-reference C type *)
                        val outTy = (case ty
                               of Ty.BoolTy => CL.T_Ptr(toCType(env, ty))
                                | Ty.StringTy => CL.T_Ptr CL.charPtr
                                | Ty.IntTy => CL.T_Ptr(toCType(env, ty))
                                | Ty.TensorTy[] => CL.T_Ptr(toCType(env, ty))
                                | Ty.TensorTy _=> toCType(env, ty)
                                | Ty.SeqTy(_, SOME _) => toCType(env, ty)
(* QUESTION: for images and sequences, it is not clear what the get function should return.
* For now, we just return 0
*)
                                | Ty.SeqTy(_, NONE) => CL.T_Ptr CL.voidPtr
                                | Ty.ImageTy _ => CL.T_Ptr CL.voidPtr
				| Ty.FemData(data) => CL.T_Ptr (toCType(env, ty)) (*QUESTION: Is this even possible? Only for cells.*)
                              (* end case *))
                        val dcls = if Ty.hasDynamicSize ty
                              then [
                                  CL.D_Proto(
                                    [], CL.T_Named "size_t", inputGetSize(spec, name),
                                    [wrldParam])
                                ]
                              else []
                        val dcls =
                              CL.D_Proto([], CL.voidTy, name, [wrldParam, CL.PARAM([], outTy, "v")])
                                :: dcls
                        in
                          dcls
                        end
                      else []
                val setDcl =
		    let
		     val useOld = not(Ty.hasTupleEsq ty) orelse (Ty.hasFem ty)
		     fun loadPrototypes (params) = [
                      CL.D_Proto(
                       [], Env.boolTy env, inputSetByName(spec, name),
                       List.@([wrldParam, CL.PARAM([], stringTy, "s")], params)),
                      CL.D_Proto(
                       [], Env.boolTy env, inputSet(spec, name),
                       List.@([wrldParam, CL.PARAM([], nrrdPtrTy, "nin")], params))
                     ]
		     (* prototypes for setting an image or dynamic sequence from a nrrd *)

		    in
		     if useOld
		     then (case ty
                            of Ty.ImageTy _ => loadPrototypes([])
                             | Ty.SeqTy(elemTy, NONE) => if APITypes.hasFem elemTy
							 then loadPrototypes([CL.PARAM([], CL.voidPtr, "data1")])
							 else loadPrototypes([])
									    
                             | _ => [
                              CL.D_Proto(
                               [], Env.boolTy env, inputSet(spec, name),
                               [wrldParam, CL.PARAM([], toCType(env, ty), "v")])
                             ]
                          (* end case *))
		     else let
			  val outTys = Ty.toOutputAbleTypes(ty)
			  val paramCTypes = List.map (fn x => trTypeAlt(env, x)) outTys
			  val tabs = List.tabulate(List.length outTys, fn x => x)
			  val params = ListPair.map (fn (idx, v) => CL.PARAM([], v, "v_" ^ (Int.toString idx))) (tabs, paramCTypes)
			 in
			  [CL.D_Proto(
			   [], Env.boolTy env, inputSet(spec, name),
			   wrldParam::params
			  )]
			 end
                    end
                in
                  descDcl @ getDcl @ setDcl
                end
        (* create a decl for an output variable *)
          (* fun mkGetDecl snapshot {ty = Ty.SeqTy(_, NONE), name, kind} = [ *)
          (*         CL.D_Proto( *)
          (*           [], Env.boolTy env, (if snapshot then snapshotGet else outputGet)(spec, name), *)
          (*           [CL.PARAM([], worldPtrTy, "wrld"), CL.PARAM([], nrrdPtrTy, "lengths"), CL.PARAM([], nrrdPtrTy, "data")]) *)
          (*       ] *)
	  fun mkGetDecl snapshot {name,ty, kind} =let
	   val outputTys = Ty.toOutputAbleTypes ty
	   fun outputParamsFold ((ty, n), params) =
	       (case ty
		 of Ty.SeqTy(_, NONE) => params@[CL.PARAM([], nrrdPtrTy, "lengths_" ^ (Int.toString n)),
						 CL.PARAM([], nrrdPtrTy, "data_" ^ (Int.toString n))]
		  | _ => params@[CL.PARAM([], nrrdPtrTy, "data_" ^ (Int.toString n))]
	       (*end case*))
			    
	   val params = List.foldr outputParamsFold [] (List.tabulate(List.length outputTys, fn x=> (List.nth(outputTys,x), x)))
	  in
	   [
             CL.D_Proto(
              [], Env.boolTy env, (if snapshot then snapshotGet else outputGet)(spec, name),
              [CL.PARAM([], worldPtrTy, "wrld")]@params)
           ]
	  end
          val outS = openOut baseName
          val placeholders =
                ("H_FILE", Out.filename outS) ::
                ("H_DEFINE", mkSymbol(#outBase spec)) ::
                subs
          val prDecl = Out.decl outS
          val prFrag = Out.fragment placeholders outS
          val prDecls = List.app prDecl
    in
     prFrag Fragments.libHHead;
     prDecl (CL.D_Verbatim ["\n/**** Functions etc. for input variables ****/\n"]);
     List.app (fn input => prDecls (mkInputDecls input)) inputs;
     case rt of SOME rt => prFrag rt | _ => ();
     prFrag Fragments.libHBody;
     if (Spec.isParallel spec)
     then prFrag Fragments.libHParExtras
     else ();
     prDecl (CL.D_Verbatim ["\n/**** Getters for output values ****/\n"]);
     if (#snapshot spec)
     then List.app (fn output => prDecls (mkGetDecl true output)) outputs
     else ();
     List.app (fn output => prDecls (mkGetDecl false output)) outputs;
     prFrag Fragments.libHFoot;
     Out.closeOut outS
    end

  end
