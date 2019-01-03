(* type-env.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)


structure TypeEnv : sig
	   type t
	   (* Create a new enviroment for a named type*)
	   val new : Atom.atom * Types.ty -> t
					       
	   (* Inserts a new method into a type env*)
	   val insertMethod : t * Atom.atom * AST.var -> unit

	   (* Inserts a new constant into a type env*)							   
	   val insertConstant : t * Atom.atom * ConstExpr.t -> unit

	   val insertHiddenVar : t * Atom.atom * ((AST.expr list) -> (AST.expr )) * Types.ty -> unit
							   
	   val findHiddenVar : t * Atom.atom -> (((AST.expr list) -> AST.expr ) * Types.ty) option

	   (* Finds the name of a type associated to an env*)
	   val findName : t -> Atom.atom

	   (* Finds the def of a type associated to an env*)
	   val findDef : t -> Types.ty
				
	   (* Finds a given method in a type env*)
	   val findMethod : t * Atom.atom -> AST.var option

	   (* Finds a give constant in a type env*)
	   val findConstant : t * Atom.atom -> ConstExpr.t option
end =  struct


    structure ATbl = AtomTable
    structure AMap = AtomMap		   
    structure Ty = Types
		     
    datatype t = TE of {
	      name : Atom.atom,
	      def : Ty.ty,
	      methods: AST.var ATbl.hash_table,
	      vars : (((AST.expr list)-> AST.expr ) * Types.ty) ATbl.hash_table,
	      constants: ConstExpr.t ATbl.hash_table
	     }
    fun new (name, def) = TE {
	 name = name,
	 def = def,
	 methods = ATbl.mkTable(16, Fail ("Method env for " ^ Atom.toString name)),
	 vars = ATbl.mkTable(16, Fail ("Var env for " ^ Atom.toString name)),
	 constants = ATbl.mkTable(16, Fail ("Constant env for " ^ Atom.toString name))
	}
			     
    fun insertMethod (TE{methods, ...}, name, var) = ATbl.insert methods (name, var)
    fun insertConstant (TE{constants, ...}, name, var) = ATbl.insert constants (name, var)
    fun insertHiddenVar (TE{vars, ...}, name, var, ty) = ATbl.insert vars (name, (var, ty))

    fun findName (TE{name, ...}) = name
				   
    fun findDef (TE{def, ...}) = def
				 
    fun findMethod (TE{methods, ...}, name) = ATbl.find methods name
    fun findConstant (TE{constants, ...}, name) = ATbl.find constants name
    fun findHiddenVar (TE{vars, ...}, name) = ATbl.find vars name
					
			 
			
end
