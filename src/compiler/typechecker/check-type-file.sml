(* check-tye-file.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)
structure CheckTypeFile  : sig
	   type constant = (Types.ty * ConstExpr.t * string)
	   val loadJson : string *  Env.context -> JSON.value option
	   val matchType : string -> (Types.ty * int list * int list) option
	   val parseConstant : Env.t * Env.context * Atom.atom * Types.ty option *  JSON.value  -> constant
	   val parseConstants : Env.t * Env.context * Atom.atom * JSON.value -> constant list
												      
	   val parseMesh : Env.t * Env.context * Atom.atom * JSON.value -> (FemData.femType option * constant list)
												
           (* val parseSpace : Env.t * Env.context * Atom.atom * FemData.space * JSON.value -> FemData.space * constant list *)
	   (* val parseFunc : Env.t * Env.context * Atom.atom * FemData.space * JSON.value -> FemData.func   * constant list *)
											    



	  end = struct
type constant = (Types.ty * ConstExpr.t * string)


structure PT = ParseTree
structure FT = FemData
structure FN = FemName
structure BD = BasisData
structure Ty = Types
structure L = Literal
structure CE = ConstExpr
structure J = JSON
structure JP = JSONParser
structure JPP = JSONPrinter
structure JU = JSONUtil
structure RE = RegExpFn(structure P=AwkSyntax
			structure E=BackTrackEngine)
structure MT = MatchTree
structure V = Vector




		
val tensorTy = "tensor((\\[([0-9]+,)*[0-9]*\\]))"
val otherbase = "string|int|bool"
val seqDims = "(\\[[0-9]+\\])+"

val tensor = RE.compileString tensorTy
val otherBase = RE.compileString otherbase
val sequenceDim = RE.compileString seqDims


fun tryTensor string =
    let
     val matchTensor = StringCvt.scanString (RE.find tensor) string
    in
     (case matchTensor
       of SOME(tree) =>
	  let
	   val {len, pos}: {len:int, pos:StringCvt.cs} = MatchTree.nth(tree,2) (* 2 must exist if tensor is matched*)
	   val next = String.extract(string, pos+len, NONE)
	   val sub = String.substring(string, pos+1, len - 1)
	   val size = String.size sub
	   val tokens = String.tokens (fn x => x = #"," orelse x = #"]" orelse x = #"[") sub
	   val tokens' = List.filter (fn x => x <> "") tokens
	   val ints = List.map (Option.valOf o Int.fromString) tokens'
	  in
	   if size = 0
	   then SOME([],next)
	   else SOME(ints, next)
	  end
	 
	| NONE => NONE
		  handle exn => NONE
     (* end case *))
    end
fun trySeqDims string =
    let
     val matchTensor = StringCvt.scanString (RE.find sequenceDim) string
    in
     (case matchTensor
       of SOME(tree) =>
	  let
	   val {len, pos}: {len:int, pos:StringCvt.cs} = MatchTree.nth(tree,0) (* 2 must exist if tensor is matched*)
	   val next = String.extract(string, pos+len, NONE)
	   val sub = String.substring(string, pos, len)
	   val tokens = String.tokens (fn x => x = #"," orelse x = #"]" orelse x = #"[") sub
	   val ints = List.map (Option.valOf o Int.fromString) tokens      
	  in
	   SOME(ints, next)
	  end
	| NONE => NONE
		  handle exn => NONE
     (* end case*))
    end

fun tryOtherStart string =
    let
     val matchTensor = StringCvt.scanString (RE.find otherBase) string
    in
     (case matchTensor
       of SOME(tree) =>
	  let
	   val {len, pos}: {len:int, pos:StringCvt.cs} = MatchTree.nth(tree,0) (* 2 must exist if tensor is matched*)
	   val next = String.extract(string, pos+len, NONE)
	   val sub = String.substring(string, pos, len)
	  in
	   SOME(sub, next)
	  end
	| NONE => NONE
		  handle exn => NONE
     (* end case*))
    end

fun matchType string =
    let
     val baseMatch = tryOtherStart string
     val tensorMatch = tryTensor string
     val baseTyAndRest = (case (baseMatch, tensorMatch)
			   of (SOME(sub, next), _) => SOME(next,
							   if sub = "bool"
							   then Ty.T_Bool
							   else if sub = "string"
							   then Ty.T_String
							   else  Ty.T_Int, [])

			    | (NONE, SOME(ints, next)) => SOME(next, Ty.T_Tensor(Ty.Shape(List.map Ty.DimConst ints)), ints)
			    | (NONE, NONE) => NONE
			 (* end case*))
     fun mkSeq(ty, []) = ty
       | mkSeq (ty, x::xs) = mkSeq(Ty.T_Sequence(ty, SOME(Ty.DimConst(x))), xs)
    in
     (case baseTyAndRest
       of NONE => NONE (* TODO: we are chosing not to validate the rest of the string! We should check there is nothing left!*)
	| SOME(next, ty, tensorInts) => (case trySeqDims(next)
					  of SOME(ints, rest) => SOME(mkSeq(ty, ints), ints, tensorInts) (* [1][2][3] -> [1,2,3]*)
					   | NONE => SOME(ty, [], tensorInts)
					(* end case *))
     (* end case *))
    end

(* we define a constant as an object with three fields: name, type, value*)
(* to hande arrays of arbitrary nested depth, we need to parse directly to AST*)
(* Why isn't this part of the Standard Basis?! Source: https://stackoverflow.com/questions/30201666/convert-array-to-list-in-sml *)
fun arrayToList  arr = V.foldr (op ::) [] arr
fun handleArray [] (x::xs) valueWrap json =
    let
     val vec = JU.asArray json
     val values = V.map (handleArray [] xs valueWrap) vec
    in
     PT.E_Cons(arrayToList values)
    end
  | handleArray (y::ys) zs valueWrap json =
    let
     val vec = JU.asArray json
     val values = V.map (handleArray ys zs valueWrap) vec
    in
     PT.E_Sequence(arrayToList values)
    end
  | handleArray [] [] valueWrap json = valueWrap json

fun oE f x = SOME(f x) handle exn => NONE
fun realToRealLit x =
    let
     fun mkReal ss = let
      val (isNeg, rest) = (case Substring.getc ss
			    of SOME(#"-", r) => (true, r)
			     | SOME(#"+", r) => (false, r)
			     | _ => (false, ss)
			  (* end case *))
      val (whole, rest) = Substring.splitl Char.isDigit rest
      val rest = Substring.triml 1 rest (* remove "." *)
      val (frac, rest) = Substring.splitl Char.isDigit rest
      val exp = if Substring.isEmpty rest
		then 0
		else let
                 val rest = Substring.triml 1 rest (* remove "e" or "E" *)
		in
                 #1(valOf(IntInf.scan StringCvt.DEC Substring.getc rest))
		end
     in
      (RealLit.real{
         isNeg = isNeg,
         whole = Substring.string whole,
         frac = Substring.string frac,
         exp = exp
      })
     end
    in
     mkReal (Substring.extract (Real.toString x, 0, NONE))
    end
fun makeString x = oE (PT.E_Lit o Literal.String o JU.asString ) x
fun makeBool x = oE (PT.E_Lit o Literal.Bool  o JU.asBool) x
fun makeInt x = oE (PT.E_Lit o Literal.Int o JU.asIntInf ) x
fun makeReal x = oE (PT.E_Lit o Literal.Real o realToRealLit o JU.asNumber) x
fun tryAll x = (Option.valOf o Option.join) (List.find Option.isSome (List.map (fn g => g x) [makeString, makeBool, makeInt, makeReal]))
					     
val bogusExp = AST.E_Lit(L.Int 0)
val bogusExpTy = (Ty.T_Error, bogusExp)

fun mkConstExpr cxt expr = Option.valOf (CheckConst.eval(cxt, false, expr))
fun extractIntConst cxt (ty, cexpr, name) =
    (case cexpr
      of ConstExpr.Int(i) => SOME(IntInf.toInt i)
       | _ => NONE)
		   
fun bogusExpTy' cxt = (Ty.T_Error, mkConstExpr cxt bogusExp, "")
		   
fun err arg = (TypeError.error arg; bogusExpTy)
		     
val err = TypeError.error

datatype token = datatype TypeError.token		

fun loadJson(fileName, cxt) = SOME(JP.parseFile fileName)
			      handle exn =>
				     (TypeError.error (cxt,
						       [S "When parsing file, ", S fileName, S ", an exception was raised: ", S(exnName exn),
							S" with a message: ", S(exnMessage exn)]); NONE)

fun typeEquality(ty, ty') =
    (case Unify.matchType(ty, ty')
      of Unify.EQ => true
       | _ => false
    (* end case *))




		
(* todo: this function needs to be refactored to provide better errors *)
fun parseConstant(env, cxt, tyName, optionalSpecTy, json) =
    let
     val findField = JU.findField json
     val name = Option.valOf (Option.map JU.asString (findField "name"))
     val tyString = Option.map JU.asString (findField "type")
     val value = Option.valOf (findField "value")
     val typeVal = Option.map matchType tyString
			      
    in
     (case Option.join typeVal
       of SOME(ty, seqInts, tensorInts) =>
	  let
	   val resultExpr = handleArray seqInts tensorInts tryAll value
	   val (astExpr, ty') = CheckExpr.check(env, cxt, resultExpr)
	   val result = (case optionalSpecTy
			  of NONE => (ty, mkConstExpr cxt astExpr, name)
			   | SOME(ty'') => (case Unify.matchType(ty, ty'')
					     of Unify.EQ => (ty, mkConstExpr cxt astExpr, name)
					      | _ => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
								 A tyName, S ", has type of ", TY(ty), S "but it is contrained to have type ",
								 TY ty'', S"."]) ;bogusExpTy' cxt)
					   (*end case*))
			(*end case *))
	  in
	   (case Unify.matchType(ty, ty')
	     of Unify.EQ => result
	      | _ => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
				 A tyName, S ", has type of ", TY(ty), S "but the actual type is ",
				 TY ty', S"."]) ;bogusExpTy' cxt)
	   (* end case*))
	  end 
	| NONE => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
			      A tyName, S " could not be parsed!"
		       ]) ;bogusExpTy' cxt)
		  handle exn => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
					    A tyName, S " could not be parsed as an exception occured!"
				     ]) ;bogusExpTy' cxt)
     (*end case*))
    end

fun parseConstants(env, cxt, tyName, json) =
    let
     val findField = JU.findField json
     val constants = findField "constants"
     
    in
     (case constants
       of SOME(a) =>  JU.arrayMap (fn x => parseConstant(env, cxt, tyName, NONE, x) ) a
	| NONE => (err (cxt, [S ", In the definition of", A tyName, S "there is no constants object"
				     ]) ;[])
				      (* end case *))

    end
fun findField x y = (JU.findField y x) handle exn => NONE
fun optionList list = SOME(List.map Option.valOf list) handle exn => NONE
fun optionTuple (x,y) = SOME((Option.valOf x, Option.valOf y)) handle exn => NONE
							   
fun constantExists(name, ty) =
    let
     fun finder (ty', expr, name') = name = name' andalso
				     (case ty
				       of  NONE => true
					| SOME(ty'') => typeEquality(ty'', ty'))
	 
    in
     List.find finder
    end

fun parseScalarBasis(env, cxt, tyName, json, dim, degree, spaceDim) =
    let
     fun maker(x : real Array.array ) : BD.t =  Option.valOf (BD.makeBasisFunc(x, dim, degree))
     fun realArray(x : JSON.value) : real Array.array = Array.fromList (JU.arrayMap (JU.asNumber) x)
     val f = maker o realArray
    in
     JU.arrayMap f json
     handle exn => (err (cxt, [S"Unable to parse polynomial basis for type ", A tyName, S".",
			S"An exception occured in the parsing:",S (exnMessage(exn)), S"."]);
     		   [BD.empty(dim, degree)])
    end
(* I don't like this.... maybe turn it into specs and process this.*)
fun parseMesh(env, cxt, tyName, json) =
    let
     val mesh = findField "mesh" json

     val constantsField = Option.mapPartial (findField "constants") mesh

     val constant = Option.map (fn x => parseConstants(env, cxt, tyName, x)) (mesh)

     val dimConstCheck = Option.valOf (Option.map (constantExists(FN.dim, SOME(Ty.T_Int))) constant)
     val meshMapConstCheck = Option.valOf (Option.map (constantExists(FN.meshMapDim, SOME(Ty.T_Int))) constant)
     val degreeConstCheck = Option.valOf (Option.map (constantExists(FN.maxDegree, SOME(Ty.T_Int))) constant)
     val dim : int option = Option.mapPartial (extractIntConst cxt) dimConstCheck
     val spaceDim : int option = Option.mapPartial (extractIntConst cxt) meshMapConstCheck
     val degree : int option = Option.mapPartial (extractIntConst cxt) degreeConstCheck
				   
     val combined = optionList [dimConstCheck, meshMapConstCheck]
     val transformShapeConst = Option.map (fn ([(ty, expr, name), (ty', expr', name')])
					      =>
						let
						 val ty'' = (Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst(2))))
						in
						 [(ty'', CE.Seq([expr', expr], ty''), FN.tds)]
						end
					  ) combined
     val newConstants = Option.map (List.@) (optionTuple(constant, transformShapeConst))


     val scalarBasis = Option.mapPartial (findField "basis") mesh 
     val parsedBasis = Option.map
			 (fn (json', dim', degree', spaceDim') =>
			     ArrayNd.fromList(parseScalarBasis(env, cxt, tyName, json', dim', degree', spaceDim')))
			 (case (scalarBasis, dim, degree, spaceDim)
			   of (SOME(a), SOME(b), SOME(c), SOME(d)) => SOME((a,b,c,d))
			    | _ => NONE (*end case *))

     val meshVal = Option.map (fn (x,y,z,basis) => FT.mkMesh(x,y,z,
							     BasisDataArray.makeUniform(basis,x),
							     tyName))
			      (case (dim, spaceDim, degree, parsedBasis)
				of (SOME(a), SOME(b), SOME(c), SOME(d)) => SOME((a,b,c,d))
				 | _ => NONE (* end case*))
     val newConstants' = Option.getOpt(newConstants, [])
    in
     (meshVal, newConstants')
    end
end
