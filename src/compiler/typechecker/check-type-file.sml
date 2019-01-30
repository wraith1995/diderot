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
												
           val parseSpace : Env.t * Env.context * Atom.atom * int list * int list * FemData.mesh * JSON.value -> FemData.femType option * constant list
																	  (* basis function shape, range shape*)
	   val parseFunc : Env.t * Env.context * Atom.atom * int list option * FemData.space * JSON.value -> FemData.femType option * constant list
															       (*function space shape*)
															       
											    



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
val otherbase = "string|int|bool|real"
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
     fun empty x = (x,[])
     val baseTyAndRest = (case (baseMatch, tensorMatch)
			   of (SOME(sub, next), _) => SOME(next,
							   if sub = "bool"
							   then Ty.T_Bool
							   else if sub = "string"
							   then Ty.T_String
							   else if sub = "int"
							   then Ty.T_Int
							   else Ty.realTy, [])

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
fun optionInt x = oE (IntInf.toInt o JU.asIntInf) x
fun optionListInt x = oE (JU.arrayMap (IntInf.toInt o JU.asIntInf)) x
fun optionRealLit x = oE ( realToRealLit o JU.asNumber) x
					     
val bogusExp = AST.E_Lit(L.Int 0)
val bogusExpTy = (Ty.T_Error, bogusExp)

fun mkConstExpr cxt expr = Option.valOf (CheckConst.eval(cxt, false, expr))
fun extractIntConst cxt (ty, cexpr, name) =
    (case cexpr
      of ConstExpr.Int(i) => SOME(IntInf.toInt i)
       | _ => NONE)

fun extractShapeConst cxt (ty, cexpr, name) =
    (case cexpr
      of ConstExpr.Seq(seq, Ty.T_Sequence(Ty.T_Int,_)) =>
	 SOME(List.map (Option.valOf o (fn x => extractIntConst cxt (Ty.T_Int, x, ""))) seq)
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
(*an object with coeffs and optional var list  or just an array
if no object, do array
if object:
-> if no var, replace it with dim
-> for each one, look for a coeff list and parse according to global
-> if there is one, override it.
*)
fun parseScalarBasis(env, cxt, tyName, json, dim, degree, spaceDim) =
    let
     fun maker(x : real Array.array ) : BD.t =  Option.valOf (BD.makeBasisFunc(x, dim, degree))
     fun realArray(x : JSON.value) : real Array.array = Array.fromList (JU.arrayMap (JU.asNumber) x)
     val f = maker o realArray
     fun arrayVersion x = let
      val array = ArrayNd.fromList(JU.arrayMap f x)
     in SOME(BasisDataArray.makeUniform(array, dim))  end
			  handle exn => NONE
			    

     fun parseAsObject(x : JSON.value, dim, degree) =
	 let
	  val vars = Option.mapPartial optionListInt (findField "vars" x)
	  val vars' = Option.getOpt(vars, List.tabulate(dim, fn x => x))
	  val dim' = List.length vars'
	  val polys = Option.map (JU.arrayMap (fn x => x )) (findField "polys" x)
	  val deg = Option.mapPartial optionInt ((findField "degree") x)
	  val degree' = Option.getOpt(deg, degree )
				     (* TODO: insert compatability checks*)
				     
	 in
	  (case polys
	    of SOME(ps) => SOME(ps, vars', dim', degree')
	     | _ => NONE)
	 
	 end
	 handle exn =>  NONE
     fun parseObjectPoly(x : JSON.value, dim, degree) =
	 let
	  val coeffs = (findField "coeffs") x
	  val parsedCoeffs = (Option.map realArray coeffs) handle exn => NONE
	  val deg = Option.mapPartial optionInt ((findField "degree") x)
	  val degree' = Option.getOpt(deg, degree)
	 in
	  (case parsedCoeffs
	    of SOME(array) => BD.makeBasisFunc(array, dim, degree')
	     | _ => NONE)
	 end
	 handle exn => NONE
			 
     fun parsePolys(xs, dim, degree) =
	 let
	  fun tryCoeff x =  let val _ = () in SOME(f x) end
			    handle exn => NONE
	  fun tryObj x = parseObjectPoly(x, dim, degree)
			 handle exn => NONE
	  fun try x = let

	  in
	   (case tryCoeff x
	     of SOME(b) => b
	      | NONE => Option.valOf ( tryObj x) ) end
	  val polys = ArrayNd.fromList (List.map try xs )
	 in
	  SOME(polys) handle exn => NONE
	 end
	 handle exn => NONE
	   

     (*parse array first*)
     (*parse object*)
     (*for each poly, parse array and then object*)
 
     fun parseBasis(x, dim, degree) = 
	 (case arrayVersion x
	   of SOME(basis) => SOME(basis) (* yuck*)
	    | NONE =>
	      ((case parseAsObject(x, dim, degree)
		of NONE => (print("umm!!!!!!!!!!!!!!!!!\n"); NONE)
		 | SOME(polys, vars, dim', degree') =>
		   (case parsePolys(polys, dim', degree')
		     of SOME(basis) => SOME(BasisDataArray.makeVarSubset(basis,vars))
		      | NONE => (print("umm!!!!!!!!!!!!!!!!!\n"); NONE)))))
	      

	   
     val leftOver = ArrayNd.fromList (List.tabulate(1, (fn _ => BD.empty(dim, degree))))
     val error = BasisDataArray.makeUniform(leftOver, 1)



    in
     Option.valOf (parseBasis(json, dim, degree))
     handle exn => (err (cxt, [S"Unable to parse polynomial basis for type ", A tyName, S".",
			S"An exception occured in the parsing:",S (exnMessage(exn)), S"."]);
     		   error)
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
			     (parseScalarBasis(env, cxt, tyName, json', dim', degree', spaceDim')))
			 (case (scalarBasis, dim, degree, spaceDim)
			   of (SOME(a), SOME(b), SOME(c), SOME(d)) => SOME((a,b,c,d))
			    | _ => NONE (*end case *))

     val meshVal = Option.map (fn (x,y,z,basis) => FT.mkMesh(x,y,z,
							     basis,
							     tyName))
			      (case (dim, spaceDim, degree, parsedBasis)
				of (SOME(a), SOME(b), SOME(c), SOME(d)) => SOME((a,b,c,d))
				 | _ => NONE (* end case*))
     val newConstants' = Option.getOpt(newConstants, [])
    in
     (meshVal, newConstants')
    end


fun parseSpace(env, cxt, tyName, basisFunctionShape : int list, rangeShape : int list, meshData, json) =
    let

     (*we will note use basisFunctionShape, but I'm keeping it here for now.*)
     val _ = if List.length basisFunctionShape = 0
	     then ()
	     else raise Fail "illegal operation"
			
     val space = findField "space" json

     val constantsField = Option.mapPartial (findField "constants") space
     val constant = Option.map (fn x => parseConstants(env, cxt, tyName, x)) (space) 
     val rangeShapeLength = List.length rangeShape
     val rangeTy = Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst rangeShapeLength))

     val dimConstCheck = Option.valOf (Option.map (constantExists(FN.dim, SOME(Ty.T_Int))) constant) handle exn => raise exn
     val spaceMapDim = Option.valOf (Option.map (constantExists(FN.spaceMapDim, SOME(Ty.T_Int))) constant) handle exn => raise exn
     val degreeConstCheck = Option.valOf (Option.map (constantExists(FN.maxDegree, SOME(Ty.T_Int))) constant) handle exn => raise exn
     val rangeConstCheck =  Option.join (Option.map (constantExists(FN.rangeShape, SOME(rangeTy))) constant) handle exn => raise exn
     
					 

     val dim : int option = Option.mapPartial (extractIntConst cxt) dimConstCheck
     val spaceDim : int option = Option.mapPartial (extractIntConst cxt) spaceMapDim
     val degree : int option = Option.mapPartial (extractIntConst cxt) degreeConstCheck
     val rangeShape : int list option = Option.mapPartial (extractShapeConst cxt) rangeConstCheck
     val rangeShape = (case rangeShape
			of SOME([1]) => SOME([])
			 | _ => rangeShape)
     val combined = optionList [spaceMapDim, rangeConstCheck]
     val spaceDofShapeConst = Option.map (fn ([(ty, expr, name), (ty', expr', name')])
					     =>
					       let
						val ty'' = (Ty.T_Sequence(Ty.T_Int, SOME(Ty.DimConst(1+rangeShapeLength))))
						val CE.Seq(lst, _) = expr'
					       in
						[(ty'', CE.Seq(expr::lst, ty''), FN.sds)]
					       end
					 ) combined
     val newConstants = Option.map (List.@) (optionTuple(constant, spaceDofShapeConst))

     val scalarBasis = Option.mapPartial (findField "basis") space
     val parsedBasis = Option.map
			 (fn (json', dim', degree', spaceDim') =>
			     (parseScalarBasis(env, cxt, tyName, json', dim', degree', spaceDim')))
			 (case (scalarBasis, dim, degree, spaceDim)
			   of (SOME(a), SOME(b), SOME(c), SOME(d)) => SOME((a,b,c,d))
			    | _ => NONE (*end case *))

     val spaceVal = Option.map (fn (x,y,basis) => FT.mkSpace(x,y, meshData,
							     basis,
							     tyName))
			      (case (spaceDim, rangeShape, parsedBasis)
				of (SOME(a), SOME(b), SOME(c)) => SOME((a,b,c))
				 | _ => NONE (* end case*))
			      
     val newConstants' : constant list = Option.getOpt(newConstants, [])




    in
     (spaceVal, newConstants') handle exn => raise exn
    end
fun parseFunc(env, cxt, tyName, expectedRangeShape, space, json) =
    let
     val _ = (case expectedRangeShape
	       of NONE => ()
		| _ => raise Fail "illegal operation")
     val func = findField "function" json

     val constantsField = Option.mapPartial (findField "constants") func
     val constant = Option.map (fn x => parseConstants(env, cxt, tyName, x)) func
     val _ = Option.app (fn x => print(Int.toString(List.length(x)))) constant


     (*these could actually be parsed and used but for now we block them to avoid confusion*)
     val optionalRange = (Option.mapPartial (constantExists(FN.rangeShape, NONE)) constant)
     val _ = (case optionalRange
	       of (NONE) => ()
		| SOME(_) => raise Fail "illegal constant assignment")

     val spaceShape = FT.spaceShape space
     val spaceDim = FT.spaceDim space
     val spaceDofShape = spaceDim :: spaceShape

     fun makeShapeConst(seq) =
	 let
	  val ty = Ty.T_Sequence(Ty.T_Int, SOME((Ty.DimConst o List.length) seq))
	 in
	  (ConstExpr.Seq(List.map (ConstExpr.Int o IntLit.fromInt) seq, ty), ty)
	 end
     val (spaceShapeTy, spaceShapeConst) = makeShapeConst spaceShape
     val (dofShapeTy, dofShapeConst) = makeShapeConst spaceDofShape
				       

     val spaceRangeConst = (spaceShapeConst, spaceShapeTy, FN.rangeShape)
     val spaceDofConst = (dofShapeConst, dofShapeTy, FN.fds)
     val constants' = spaceDofConst :: spaceRangeConst :: (Option.getOpt(constant, []))
     val funcType = SOME(FT.mkFunc(space, spaceShape, tyName))
				       

				  

	       
			   
			   
    in
     (funcType, constants')
    end
end
