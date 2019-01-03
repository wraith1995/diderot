(* inputs.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure Inputs : sig

  (* input initialization *)
    datatype 'var input_init
      = NoDefault                       (* uninitialized input *)
      | ConstExpr                       (* initialized to constant.  The actual initialization
                                         * is handled in code.
                                         *)
      | LoadSeq of string               (* load a sequence from a file *)
      | Proxy of string * ImageInfo.t   (* input image specified by proxy *)
      | Image of ImageInfo.t            (* input image w/o proxy *)
      | LoadFem of FemData.femType * ('var option)

    datatype 'var input = INP of {
        var : 'var,                     (* the global input variable *)
        name : string,                  (* the source-language input name *)
        ty : APITypes.t,                (* source-language type *)
        desc : string option,           (* the optional description *)
        init : 'var input_init               (* the default initialization (or NoDefault) *)
      }

    val varOf : 'a input -> 'a

    val same : 'a input * 'a input -> bool

    val hash : 'a input -> word

    val initToString : 'var input_init -> string

    val toString : 'a input -> string

    val imageInfo : 'a input -> ImageInfo.t option

  (* return true if the initialization specifies that there is a default value *)
    val isDefault : 'var input_init -> bool

  (* return info about a proxy image (or NONE) *)
    val proxy : 'var input_init -> ImageInfo.t option

  (* type conversion *)
    val map : ('a -> 'b) -> 'a input -> 'b input

  end = struct

    datatype 'var input_init
      = NoDefault                       (* uninitialized input *)
      | ConstExpr                       (* initialized to constant.  The actual initialization
                                         * is handled in the constInit block.
                                         *)
      | LoadSeq of string               (* load a sequence from a file *)
      | Proxy of string * ImageInfo.t   (* input image specified by proxy *)
      | Image of ImageInfo.t            (* input image w/o proxy *)
      | LoadFem of FemData.femType * ('var option)

    datatype 'var input = INP of {
        var : 'var,                     (* the global input variable *)
        name : string,                  (* the source-language input name *)
        ty : APITypes.t,                (* source-language type *)
        desc : string option,           (* the optional description *)
        init : 'var input_init               (* the default initialization (or NoDefault) *)
      }

    fun varOf (INP{var, ...}) = var

    fun same (INP{name=n1, ...}, INP{name=n2, ...}) = (n1 = n2)

    fun hash (INP{name, ...}) = HashString.hashString name

    fun initToString NoDefault = "<no-default>"
      | initToString ConstExpr = "<default>"
      | initToString (LoadSeq name) = String.concat["load(\"", name, "\")"]
      | initToString (Proxy(name, _)) = String.concat["image(\"", name, "\")"]
      | initToString (Image info) = ImageInfo.toString info
      | initToString (LoadFem(data, var)) =
	let
	 val varStr = (case var
			of SOME(v) => "SOME(var)"
			 | NONE => "NONE"
		      (* end case*))
	in
	 String.concat["loadFem(", (FemData.toString data), ", ",varStr, ")" ]
	end

    fun toString (INP{name, desc=NONE, init=NoDefault, ...}) = name
      | toString (INP{name, desc=SOME desc, init=NoDefault, ...}) =
          String.concat[name, "(\"", String.toString desc, "\")"]
      | toString (INP{name, desc=NONE, init, ...}) =
          String.concat[name, " = ", initToString init]
      | toString (INP{name, desc=SOME desc, init, ...}) =
          String.concat[name, "(\"", String.toString desc, "\") = ", initToString init]

    fun imageInfo (INP{init=Proxy(_, info), ...}) = SOME info
      | imageInfo (INP{init=Image info, ...}) = SOME info
      | imageInfo _ = NONE

    fun proxy (Proxy(_, info)) = SOME info
      | proxy _ = NONE

    fun initConver( init : 'a input_init , f : 'a -> 'b ) : 'b input_init =
	(case init
	  of LoadFem(data, SOME(v)) => LoadFem(data, SOME( f v))
	   | LoadFem(data, NONE) => LoadFem(data, NONE)
	   | NoDefault => NoDefault
	   | Proxy(a,b) => Proxy(a,b)
	   | Image(a) => Image(a)
	   | LoadSeq(s) => LoadSeq(s)
	   | ConstExpr => ConstExpr
	(* end case*))

  (* type conversion *)
    fun map f (INP{var, name, ty, desc, init}) =
	let
	 val init' = initConver(init, f)
	in
         INP{var = f var, name = name, ty = ty, desc = desc, init = init'}
	end

    fun isDefault NoDefault = false
      | isDefault (Image _) = false
      | isDefault _ = true

  end
