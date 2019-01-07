(* check-tye-file.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)
structure CheckTypeFile  : sig


	   val matchType : string -> (Types.ty * int list) option
	   val parseConstant : Env.t * Env.context * Atom.atom * JSON.value * Types.ty option -> Types.ty * ConstExpr.t * string
	   val parseConstants : Env.t * Env.context * Atom.atom * JSON.value -> (Types.ty * ConstExpr.t * string) list
												      
	   val parseMesh : Env.t * Env.context * Atom.atom * JSON.value -> (FemData.mesh option * (Types.ty * ConstExpr.t * string) list)
															 (* grab specific constants; parse basis can be general*)
	   (* val parseSpace : Env.t * Env.context * Atom.atom * FemData.space * JSON.value -> FemData.space *)
	   (* val parseFunc : Env.t * Env.context * Atom.atom * FemData.space * JSON.value -> FemData.func *)
											    

									(* step 2: make fem spec, and parse to femData *)
									(* step 3: drop idea of result, Other, etc -> just make parseX where X=mesh,...,..,..*)
									(* step 4: reviese type checker to handle constants, checking of fem data against file info*)
									(* step 5: infustructure for transform sampling*)

	  end = struct

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

(* fun extractMatchFromTree(string, idx, tree, upOffset : int, downOffset : int) = *)
(*     (case MatchTree.nth(tree, idx) *)
(*      of NONE => NONE *)
(*       | SOME{pos, len} => SOME(String.substring(string, pos + upOffset, len - downOffset), *)
(* 			     pos + upOffset + len - downOffset) *)
(*     (* end case*)) *)

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
       | mkSeq (ty, x::xs) = mkSeq(Ty.T_Sequence(ty, SOME(Ty.DimConst(x)), xs))
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
     val values = V.map (handleArray xs valueWrap) vec
    in
     PT.E_Tensor(arrayToList values)
    end
  | handleArray (y::ys) zs valueWrap json =
    let
     val vec = JU.asArray json
     val values = V.map (handeArray ys zs valueWrap) vec
    in
     PT.E_Sequence(arrayToList values)
    end
  | handleArray [] [] valueWrap json = valueWrap json

fun oE f x = SOME(f x) handle exn => NONE
fun makeString x = PT.E_Lit (Literal.String (oE JU.asString x))
fun makeBool x = PT.E_Lit (Literal.Bool (oE JU.asBool x))
fun makeInt x = PT.E_Lit (Literal.Int (oE JU.asIntInf x))
fun makeReal x = PT.E_Lit (Literal.Real (oE JU.asNumber x))
fun tryAll x = (Option.valof o Option.valof) (List.find Option.isSome (List.map (fn g => g x) [makeString, makeBool, makeInt, makeReal]))
					     
val bogusExp = AST.E_Lit(L.Int 0)
val bogusExpTy = (Ty.T_Error, bogusExp)
val bogusExpTy' = (Ty.T_Error, bogusExp, "")
		   
fun err arg = (TypeError.error arg; bogusExpTy)

fun typeEquality ty ty' =
    (case Unify.matchType(ty, ty')
      of Unify.EQ => true
       | _ => false
    (* end case *))

fun mkConstExpr cxt expr = CheckConst.CheckConst(cxt, false, expr)
fun extractIntConst cxt expr =
    (case mkConstExpr cxt expr
      of SOME(ConstExpr.Int(i)) => SOME(IntInf.toInt i)
       | NONE => NONE)



		
(* todo: this function needs to be refactored to provide better errors *)
fun parseConstant(env, cxt, tyName, optionalSpecTy, json) =
    let
     val findField = JU.findField json
     val name = JU.asString (findField "name")
     val tyString = JU.asString (findField "type")
     val value = JU.asString (findField "value")
     val typeVal = matchType tyString
    in
     (case typeVal
       of SOME(ty, seqInts, tensorInts) =>
	  let
	   val resultExpr = handleArray seqInts tensorInts tryAll json
	   val (astExpr, ty') = CheckExpr.check(env, cxt, resultExpr)
	   val result = (case optionalSpecTy
			  of NONE => (mkConstExpr astExpr, ty)
			   | SOME(ty'') => (case Unify.matchType(ty, ty'')
					     of Unify.EQ => (ty, astExpr, name)
					      | _ => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
									A tyName, S ", has type of ", TY(ty), S "but it is contrained to have type ",
									TY ty'', S"."]) ;bogusExpTy')
					   (*end case*))
			(*end case *))
	  in
	   (case Unify.matchType(ty, ty')
	     of Unify.EQ => result
	      | _ => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
				 A tyName, S ", has type of ", TY(ty), S "but the actual type is ",
				 TY ty', S"."]) ;bogusExpTy')
	   (* end case*))
	  end 
	| NONE => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
			      A tyName, S " could not be parsed!"
		       ]) ;bogusExpTy')
		  handle exn => (err (cxt, [S "declared type of constant", S name, S ", in definition of",
					    A tyName, S " could not be parsed!"
				     ]) ;bogusExpTy')
     (*end case*))
    end

fun parseConstants(env, cxt, tyName, json) =
    let
     val findField = oE JU.findField json
     val constants = findField "constants"
     
    in
     (case constants
       of SOME(a) =>  JU.arrayMap (fn x => parseConstant(env, cxt, tyName, NONE, x) ) a
	| NONE => (err (cxt, [S ", In  the definition of", A tyName, S "there is no constants object"
				     ]) ;bogusExpTy)
				      (* end case *))

    end
fun findField x y = SOME(JU.findField y x) handle exn => NONE
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
     fun maker x =  BD.makeBasis(x, dim, degree)
     fun realArray x = JU.arrayMap (JU.asNumber) x
     fun f = maker o realArray
    in
     JU.arrayMap f json
     handle exn => (err (cxt, [S"Unable to parse polynomial basis for type ", A tyName, S"."]);
		   [BD.empty(dim, degree)])
    end
fun parseMesh(env, cxt, tyName, json) =
    let
     val mesh = findField "mesh" json
     val constantsField = Option.mapPartial (findField "constants") mesh
     val constant = Option.map (fn x => parseConstants(env, cxt, tyName, x)) constantsField
			       
     val dimConstCheck = (Option.join o Option.map) (constantExists(FN.dim, Ty.T_Int)) constant
     val meshMapConstCheck = (Option.join o Option.map) (constantExists(FN.meshMapDim, Ty.T_Int)) constant
     val degreeConstCheck = (Option.join o Option.map) (constantExists(FN.maxDegree, Ty.T_Int)) constant
					  
     val dim = extractIntConst cxt dimConstCheck
     val spaceDim = extractIntConst cxt meshMapConstCheck
     val degree = degreeConstCheck cxt degreeConstCheck
				   
     val combined = optionList [dimConstCheck, meshMapConstCheck]
     val transformShapeConst = Option.map (fn ([(ty, expr, name), (ty', expr', name')])
					      =>
						let
						 val ty'' = (Ty.T_Sequence(Ty.T_Int, 2))
						in
						 [(ty'', CE.Seq([expr', expr], ty''), FN.tds)]
						end
					  ) combined
     val newConstants = Option.map (List.@) (optionTuple(constants, transformShapeConst))

     (* now we get the polynomial are we get ints...*)
     val scalarBasis = findField "poly" json
     val parsedBasis = Option.map
			 (fn [json, dim, degree, spaceDim] =>
			     parseScalarBasis(env, cxt, tyName, json, dim, degree, spaceDim))
			 (optionList ([scalarBasis, dim, degree, spaceDim]))

     val meshVal = Option.map (fn [x,y,basis] => FD.mkMesh(x,y,basis)) (optionList ([dim, spaceDim,parsedBasis]))
     val newConstants' = Option.getOpt(newConstants, [])
			      


    in
     (meshVal, newConstants')
    end
    
      
    


end
