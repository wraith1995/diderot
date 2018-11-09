(* type-env.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)


structure TypeEnv : sig
	   type t
	   (* Create a new enviroment for a type*)
	   val new : Atom.atom * Types.ty -> t
	   val insertMethod : t * Atom.atom * AST.var -> unit
	   val insertConstant : t * Atom.atom * ConstExpr.t -> unit
	   val findName : t -> Atom.atom
	   val findDef : t -> Types.ty

	   val findMethod : t * Atom.atom -> AST.var option
	   val findConstant : t * Atom.atom -> ConstExpr.t option
		  
	  end =  struct


    structure ATbl = AtomTable
    structure AMap = AtomMap		   
    structure Ty = Types
		     
    datatype t = TE of {
	      name : Atom.atom,
	      def : Ty.ty,
	      methods: AST.var ATbl.hash_table,
	      constants: ConstExpr.t ATbl.hash_table
	     }
    fun new (name, def) = TE {
	 name = name,
	 def = def,
	 methods = ATbl.mkTable(16, Fail ("Method env for " ^ Atom.toString name)),
	 constants = ATbl.mkTable(16, Fail ("Constant env for " ^ Atom.toString name))
	}
			     
    fun insertMethod (TE{methods, ...}, name, var) = ATbl.insert methods (name, var)
    fun insertConstant (TE{constants, ...}, name, var) = ATbl.insert constants (name, var)

    fun findName (TE{name, ...}) = name
				   
    fun findDef (TE{def, ...}) = def
				 
    fun findMethod (TE{methods, ...}, name) = ATbl.find methods name
    fun findConstant (TE{constants, ...}, name) = ATbl.find constants name
					
			 
			
end
