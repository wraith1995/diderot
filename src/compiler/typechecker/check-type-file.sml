(* check-tye-file.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)
structure CheckTypeFile  : sig
	   type constants = (string * Types.ty * AST.expr ) list
	   datatype other = M of FT.mesh | S of FT.space | F of FT.func
	   datatype result = C of constants | Fem of constants * other

	   val matchType : string -> (Types.ty * int list) option
	   val parseConstant : JSON.value -> Types.ty * ParseTree.expr
	   val parseConstants : JSON.value -> Types.ty * ParseTree.expr list
									(* step 1: send env, error context, and remake parse constant*)
									(* step 2: make fem spec, and parse to femData *)
									(* step 3: drop idea of result, Other, etc -> just make parseX where X=mesh,...,..,..*)
									(* step 4: reviese type checker to handle constants, checking of fem data against file info*)
									(* step 5: infustructure for transform sampling*)

	  end = struct

structure PT = ParseTree
structure FT = FemData
structure Ty = Types
structure L = Literal
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
fun makeString x = SOME(PT.E_Lit (Literal.String (JU.asString x))) handle exn => NONE
fun makeBool x = SOME(PT.E_Lit (Literal.Bool (JU.asBool x))) handle exn => NONE
fun makeInt x = SOME(PT.E_Lit (Literal.Int (JU.asIntInf x))) handle exn => NONE
fun makeReal x = SOME(PT.E_Lit (Literal.Real (JU.asNumber x))) handle exn => NONE
fun tryAll x = (Option.valof o Option.valof) (List.find Option.isSome (List.map (fn g => g x) [makeString, makeBool, makeInt, makeReal]))
			    
(*real, int, string, bool to ast*)
fun parseConstant json =
    let
     val findField = JU.findField json
     val name = findField "name"
     val tyString = findField "type"
     val value = findField "value"

     val typeVal = matchType tyString
    in
     (case typeVal
       of SOME(ty, seqInts, tensorInts) => SOME(ty, handleArray seqInts tensorInts tryAll json) 
	| NONE => NONE
		  handle exn => NONE
     (* end case *) )
     
    end
    
      
    

		 (*figure out depth and build parser*)
		 (*check-expr*)
		 (*do this in a list*)
		 (* check for mesh key*)
type constants = (string * Ty.ty * AST.expr ) list
datatype other = M of FT.mesh | S of FT.space | F of FT.func
datatype result = C of constants | Fem of constants * other
end
