(* global-env.sml
 *
 * The global environment is an imperative environment that tracks global variables,
 * program properties, and the name of the strand.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)

structure GlobalEnv : sig

    type t

  (* create a new (empty) global environment *)
    val new : unit -> t

  (* functions are either user-defined or pre-defined *)
    datatype fun_def
      = PrimFun of AST.var list                 (* possibly overloaded builtin function *)
      | UserFun of AST.var                      (* user-defined function *)
      | OverloadUserFun of AST.var list         (* overloaded user-defined function*)
  (* lookup a strand *)
    val findStrand : t * Atom.atom -> StrandEnv.t option

  (* lookup a function definition *)
    val findFunc : t * Atom.atom -> fun_def

  (* lookup a global variable *)
    val findVar : t * Atom.atom -> AST.var option

  (* lookup a kernel *)
    val findKernel : t * Atom.atom -> Kernel.t option

    val findTypeEnv : t * Atom.atom -> TypeEnv.t option

  (* add the strand to the environment; this function returns the strand
   * environment that we use to track strand state variables.
   *)
    val insertStrand : t *  StrandEnv.t -> unit

  (* insert a function binding *)
    val insertFunc : t * Atom.atom * fun_def -> unit

  (* insert a global variable binding *)
    val insertVar : t * Atom.atom * AST.var -> unit

  (* insert a kernel binding*)
    val insertKernel : t * Atom.atom * Kernel.t -> unit

    val insertNamedType : t * Atom.atom * Types.ty * (Atom.atom * ConstExpr.t ) list  * (Atom.atom * AST.var ) list -> unit

  (* tracking program features *)
    val recordProp : t * Properties.t -> unit
    val hasProp : t * Properties.t -> bool
    val properties : t -> Properties.t list

    val cleanEnv : t -> t

  end = struct

    structure ATbl = AtomTable
    structure AMap = AtomMap
    structure TE = TypeEnv
    structure TY = Types

  (* functions are either user-defined or pre-defined *)
    datatype fun_def
      = PrimFun of AST.var list                 (* possibly overloaded builtin function *)
      | UserFun of AST.var                      (* user-defined function *)
      | OverloadUserFun of AST.var list         (* overloaded user-defined function*)


  (* global environment holds global variables and strands *)
    datatype t = GE of {
        strand : StrandEnv.t option ref,        (* the strand *)
        fEnv : fun_def ATbl.hash_table,         (* functions, which may be overloaded *)
        vEnv : AST.var ATbl.hash_table,         (* global variable bindings *)
        kEnv : Kernel.t ATbl.hash_table,        (* kernel bindings *)
	tEnv :  TE.t ATbl.hash_table,            (* type environment *)
        props : Properties.t list ref           (* record program properties *)
      }

    fun new () = GE{
            strand = ref NONE,
            fEnv = ATbl.mkTable(128, Fail "global function env"),
            vEnv = ATbl.mkTable(128, Fail "global variable env"),
            kEnv = ATbl.mkTable(16, Fail "kernel env"),
	    tEnv = ATbl.mkTable(16, Fail "type env"),
            props = ref[]
          }

  (* add the strand to the environment; this function returns the strand
   * environment that we use to track strand state variables.
   *)
    fun insertStrand (GE{strand, ...}, sEnv) = (case !strand
           of NONE => strand := SOME sEnv
            | _ => raise Fail "impossible: redefinition of strand"
          (* end case *))

    fun insertFunc (GE{fEnv, ...}, f, f') = ATbl.insert fEnv (f, f')

    fun insertVar (GE{vEnv, ...}, x, x') = ATbl.insert vEnv (x, x')

    fun insertKernel (GE{kEnv, ...}, k, k') = ATbl.insert kEnv (k, k')

    fun insertNamedType (GE{tEnv, ...}, name, def, constants, methods) =
	let
	 val newType = TE.new(name, def)
	 val _ = List.app (fn (n,c) => TE.insertConstant(newType,n,c)) constants
	 val _ = List.app (fn (n,c) => TE.insertMethod(newType,n,c)) methods
				
	in
	 (ATbl.insert tEnv (name, newType))
	end

  (* lookup a strand *)
    fun findStrand (GE{strand, ...}, name) = (case !strand
           of NONE => NONE
            | someEnv as SOME sEnv => if Atom.same(name, StrandEnv.strandName sEnv)
                then someEnv
                else NONE
          (* end case *))

    fun findFunc (GE{fEnv, ...}, x) = (case ATbl.find fEnv x
           of NONE => PrimFun[]
            | SOME fdef => fdef
          (* end case *))

    fun findVar (GE{vEnv, ...}, x) = ATbl.find vEnv x

    fun findKernel (GE{kEnv, ...}, x) = ATbl.find kEnv x

    fun findTypeEnv (GE{tEnv, ...}, x) = ATbl.find tEnv x

  (* tracking program features *)
    fun recordProp (GE{props, ...}, p) = if (Properties.hasProp p (!props))
          then ()
          else props := p :: !props

    fun hasProp (GE{props, ...}, p) = Properties.hasProp p (!props)

    fun properties (GE{props, ...}) = !props

    
    fun cleanAstVarTable table = ATbl.map CleanAst.varCleanType table
    fun cleanFunTable table = ATbl.map (fn x => case x
						 of PrimFun(vars) => PrimFun(List.map CleanAst.varCleanType vars)
						  | UserFun(var) => UserFun(CleanAst.varCleanType var)
						  | OverloadUserFun(vars) => OverloadUserFun(List.map CleanAst.varCleanType vars)
				       (*end case*)
				       ) table					
					

    fun cleanEnv (GE{strand, fEnv, vEnv, kEnv, tEnv, props}) =
	GE{
	 strand = ref (Option.map StrandEnv.cleanStrand (!strand)),
	 fEnv = cleanFunTable fEnv,
	 vEnv = cleanAstVarTable vEnv,
	 kEnv = kEnv,
	 tEnv = tEnv,
	 props = props
	}

  end
