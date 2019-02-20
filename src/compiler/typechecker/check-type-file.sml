(* check-tye-file.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 *)
structure CheckTypeFile  : sig
	   type constant = (Types.ty * ConstExpr.t * string)
	   datatype dimensionalObjects = Points of (int * real list) list * int list option
				       | Higher of int * ((int * IntInf.int list * IntInf.int list) list)
				       | Mapping of int * string
				       | LineParam of (real list * real list * real) list
				       (*b-a, a, |b-a|*)
				       | PlaneParam of (real * real list) list * real list list list list (*d, normal, matrix of matrices*)
	   
	   val loadJson : string *  Env.context -> JSON.value option
	   val matchType : string -> (Types.ty * int list * int list) option
	   val parseConstant : Env.t * Env.context * Atom.atom * Types.ty option *  JSON.value  -> constant
	   val parseConstants : Env.t * Env.context * Atom.atom * JSON.value -> constant list
												      
	   val parseMesh : Env.t * Env.context * Atom.atom * JSON.value -> (FemData.femType option * constant list * dimensionalObjects list)
												
           val parseSpace : Env.t * Env.context * Atom.atom * int list * int list * FemData.mesh * JSON.value -> FemData.femType option * constant list
																	  (* basis function shape, range shape*)
	   val parseFunc : Env.t * Env.context * Atom.atom * int list option * FemData.space * JSON.value -> FemData.femType option * constant list
															       (*function space shape*)
															       
	   val realToRealLit : real -> RealLit.t								    


	  end = struct
type constant = (Types.ty * ConstExpr.t * string)
datatype dimensionalObjects = Points of (int * real list) list * int list option
			    | Higher of int * ((int * IntInf.int list * IntInf.int list) list)
			    | Mapping of int * string
			    | LineParam of (real list * real list * real) list
			    (*b-a, a*)
			    | PlaneParam of (real * real list) list * real list list list list (*d, normal, matrix of matrices*)

(*index of a point and its coords 
| dim of objects and index to higher dimensional objects with addtional objects with index to points for figuring out the plane!= *)

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
     val preProc  = String.implode o (List.map (fn #"~" => #"-" | a => a)) o String.explode
     fun mkReal ss = let
      val (isNeg, rest) = (case Substring.getc ss
			    of SOME(#"-", r) => (true, r)
			     | SOME(#"+", r) => (false, r)
			     | _ => (false, ss)
			  (* end case *))
      val (whole, rest) = Substring.splitl Char.isDigit rest
      val rest = (case Substring.getc rest
		   of SOME(#".", _) => Substring.triml 1 rest (* remove "." if it exists*)
		    | _ => rest
		  (* end case*))
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
     mkReal (Substring.extract (preProc(Real.toString x), 0, NONE))
    end
fun makeString x = oE (PT.E_Lit o Literal.String o JU.asString ) x
fun makeBool x = oE (PT.E_Lit o Literal.Bool  o JU.asBool) x
fun makeInt x = oE (PT.E_Lit o Literal.Int o JU.asIntInf ) x
fun makeReal x = oE (PT.E_Lit o Literal.Real o realToRealLit o JU.asNumber) x
fun tryAll x = (Option.valOf o Option.join) (List.find Option.isSome (List.map (fn g => g x) [makeString, makeBool, makeInt, makeReal]))
fun optionInt x = oE (IntInf.toInt o JU.asIntInf) x
fun optionListInt x = oE (JU.arrayMap (IntInf.toInt o JU.asIntInf)) x
fun optionRealLit x = oE ( realToRealLit o JU.asNumber) x
fun optionBool x = oE (JU.asBool) x
					     
val bogusExp = AST.E_Lit(L.Int 0)
val bogusExpTy = (Ty.T_Error, bogusExp)
fun mkConstExpr cxt expr = Option.valOf (CheckConst.eval(cxt, false, expr))
fun bogusExpTy' cxt = (Ty.T_Error, mkConstExpr cxt bogusExp, "")
			
fun err arg = (TypeError.error arg; bogusExpTy)
		
val err = TypeError.error

datatype token = datatype TypeError.token		


fun makeParseError(err, cxt, field, msg, default) = (err (cxt, [S "Unable to parse field, ", S field, S" with message:", S msg]); default)


fun extractIntConst cxt (ty, cexpr, name) =
    (case cexpr
      of ConstExpr.Int(i) => SOME(IntInf.toInt i)
       | _ => NONE)

fun extractShapeConst cxt (ty, cexpr, name) =
    (case cexpr
      of ConstExpr.Seq(seq, Ty.T_Sequence(Ty.T_Int,_)) =>
	 SOME(List.map (Option.valOf o (fn x => extractIntConst cxt (Ty.T_Int, x, ""))) seq)
       | _ => NONE)
		   

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

fun warnIfNo(env, cxt, field1, field2, tyName) json =
    (case json
      of NONE => NONE
       | SOME(json') => (case findField field1 json'
			  of SOME(a) => SOME(a)
			   | NONE => (TypeError.warning(cxt,
							[S "Expected field, ", S(field1),
							 S", in parsing of another field, ", S(field2),
							 S", in the pasing of type", A(tyName), S"."]);
				      NONE)
			(*end case*))
    (*end case*))

fun errIfNo(env, cxt, field1, field2, tyName) json =
    (case json
      of NONE => NONE
       | SOME(json') => (case findField field1 json'
			  of SOME(a) => SOME(a)
			   | NONE => (TypeError.error(cxt,
							[S "Expected field, ", S(field1),
							 S", in parsing of another field, ", S(field2),
							 S", in the pasing of type", A(tyName), S"."]);
				      NONE)
			(*end case*))
     (*end case*))       


(*NOTE: this is only needed because there is no RealLit.t or Rational.t arithemtic...*)
fun ugg(f, [], _)= []
  | ugg(f : 'x * 'y -> 'z, (x::xs) : 'x list, (y::ys) : 'y list) =
    let
     val temp1 = (List.map (fn y => f(x,y)) ys)
     val temp2 = ugg(f,xs, (y::ys))
    in
     temp1::temp2
    end
fun subtractSeg([x,y]) = (ListPair.map (Real.-) (y,x), x)
  | subtractSeg _ = raise Fail "invalid usbtract"
fun analyzeGeometry1(cxt, points, higher1) =
    let
     val Points(xs, _) = points
     val vecs = List.map (fn (x,y) => y) xs
     fun getVector i = List.nth(vecs, i)
     val Higher(1, xs) = higher1
     val segs = List.map (fn (x,y,z) => z) xs
     val segs' = List.map (List.map IntInf.toInt) segs
     val pairSegs' = List.map (subtractSeg o (List.map getVector)) segs'
			      (*To do *)
     fun norm xs = Math.sqrt(List.foldr (fn (x,y) => x*x + y) 0.0 xs)
     val pairSegs'' = List.map (fn (x,y)=> (x, y, norm x)) pairSegs'
     val testNorms = (List.exists (fn (x,y,z) => Real.<=(z, 0.00001)) pairSegs'')
     val _ = if testNorms
	     then TypeError.warning(cxt, [S"A line defined in a json file is probably degenerate!"])
	     else ()
    in
     LineParam(pairSegs'')
    end
fun subtract(a,b) = (ListPair.map (Real.-)(a,b))
fun crossProduct ([a1,a2,a3] : real list, [b1, b2, b3]) =
    [a2*b3 - a3*b2, a1*b3 - a3*b1, a1*b2 - a2*b1]
fun planeNormal([p1,p2,p3]) =
    let
     val n =  crossProduct(subtract(p2,p1), subtract(p3,p1))
     val norm =  Math.sqrt(List.foldr (op+) 0.0 (List.map (fn x => x*x) n))
     val n' = List.map (fn x => x/norm) n
     val dot = List.foldr (op+) 0.0 (ListPair.map Real.* (n', p1))
     val p4 = ListPair.map Real.+ (n', p1)
    in
     (n',dot, p4) (*plane of the form x\in P iff n*x = dot*)
    end

fun buildTansform(pl1 : real list list, pl2 : real list list) =
    let
     val [p1,p2,p3,p4] = pl1 (*from i.e inverted*)
     val [q1,q2,q3,q4] = pl2 (*to i.e not inverted*)
     val [a11,a21,a31] = p1
     val [a12,a22,a32] = p2
     val [a13,a23,a33] = p3
     val [a14,a24,a34] = p4
     val [b11,b21,b31] = q1
     val [b12,b22,b32] = q2
     val [b13,b23,b33] = q3
     val [b14,b24,b34] = q4

     val result = [[((~(a23*a32) + a22*a33)*b11)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a23*a31 - a21*a33)*b12)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a22*a31) + a21*a32)*b13)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a13*a32 - a12*a33)*b11)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a13*a31) + a11*a33)*b12)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a12*a31 - a11*a32)*b13)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((~(a13*a22) + a12*a23)*b11)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a13*a21 - a11*a23)*b12)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a12*a21) + a11*a22)*b13)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34)*b11)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34)*b12)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34)*b13)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + b14],
   [((~(a23*a32) + a22*a33)*b21)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a23*a31 - a21*a33)*b22)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a22*a31) + a21*a32)*b23)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a13*a32 - a12*a33)*b21)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a13*a31) + a11*a33)*b22)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a12*a31 - a11*a32)*b23)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((~(a13*a22) + a12*a23)*b21)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a13*a21 - a11*a23)*b22)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a12*a21) + a11*a22)*b23)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34)*b21)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34)*b22)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34)*b23)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + b24],
   [((~(a23*a32) + a22*a33)*b31)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a23*a31 - a21*a33)*b32)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a22*a31) + a21*a32)*b33)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a13*a32 - a12*a33)*b31)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a13*a31) + a11*a33)*b32)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a12*a31 - a11*a32)*b33)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((~(a13*a22) + a12*a23)*b31)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a13*a21 - a11*a23)*b32)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a12*a21) + a11*a22)*b33)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33),
    ((a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34)*b31)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((~(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34)*b32)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + 
     ((a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34)*b33)/(~(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33) + b34],[0.0,0.0,0.0,1.0]]
    in
     result
    end
					  
fun analyzeGeometry2(points, higher2) =
    let
     val Points(xs, _) = points
     val vecs = List.map (fn (x,y) => y) xs
     fun getVector i = List.nth(vecs, i)

     val Higher(2, xs) = higher2
     val planes = List.map (fn (x,y,z) => (List.map (getVector o IntInf.toInt)) z) xs
     val numPlanes = List.length planes
     val computedInfo = List.map planeNormal planes
     val extraPoints = List.map (fn (x,y,z) => z) computedInfo
     val dNormPairs = List.map (fn (x,y,z) => (y, x)) computedInfo
     val affineDefs = ListPair.map (fn ([x,y,z], a) => [x,y,z,a]) (planes, extraPoints)
     (*for x in affineDef, for y in affineDef: *)

     val transforms : real list list list list = ugg(buildTansform, affineDefs, affineDefs) (*[x][y] maps from facet x to facet y*)
						    (*todo: build *)
    in
     PlaneParam(dNormPairs, transforms)
    end


      
fun parseGeometryOfRefCell(env, cxt, dim, json, meshName) =
    let
     val geometry = errIfNo(env, cxt, "geometry", "refCell", meshName) json
     val verticies = errIfNo(env, cxt, "vertices", "geometry", meshName) geometry
     val vertToNode = errIfNo(env, cxt, "toNode", "geometry", meshName) geometry


     fun parseToNode(SOME(json)) =
	 let
	  val array = (JU.arrayMap (IntInf.toInt o JU.asIntInf)) json
	 in
	  SOME(array)
	  handle exn => NONE
	 end
       | parseToNode (NONE) = NONE
     val vertToNode' = parseToNode(vertToNode)

     val zero = List.tabulate(dim, fn x => 0.0)
     fun parsePoint (idx, jsonPoint) =
	 let
	  val array = JU.arrayMap (JU.asNumber) jsonPoint
	  val _ = if List.length array = dim
		  then ()
		  else raise (Fail "...")
			     
	 in
	  (idx,array)
	  handle exn => ((TypeError.error(cxt,[S "Unable to process point ",
					       S(Int.toString idx),
					       S" for mesh ", A(meshName),
					       S". Chechk that all verticies have the correct dim and are arrays of reals."]); (idx, zero)))
	 end

     fun parseVerticies'(SOME(verts)) = Vector.foldr (op ::) [] (JU.asArray(verts))
       | parseVerticies'(NONE) = []
				 handle exn => ((TypeError.error(cxt,[S "Unable to process verticies for mesh named", A(meshName), S"."]); []))
											

     fun parseVerticies(verts) = let val len = List.length verts
				 in
				  Points(ListPair.map parsePoint (List.tabulate(len, fn x =>x ), verts), vertToNode')
				 end
     val verts = parseVerticies(parseVerticies'(verticies))

     val higherObjects = List.tabulate(dim - 1 , fn x => (x+1, errIfNo(env, cxt, "object"^(Int.toString (x+1)), "geometry", meshName) geometry))
     fun parserHigherEntry(subDim, json, index) =
	 let
	  val name = "object"^(Int.toString subDim)
	  val lowerObjects = errIfNo(env, cxt, "entity", name, meshName) (SOME(json))
	  val planePoints = errIfNo(env, cxt, "plane", name, meshName) (SOME(json))

	  fun arrayOfInts(json, ty) = JU.arrayMap (JU.asIntInf) json
				      handle exn => (TypeError.error(cxt, [S "Unable to process ", S ty, S" of ", S(name),
									   S(" entry "^(Int.toString index)),
									   S" in ", A(meshName), S"."]); [])

	  val entities = Option.getOpt(Option.map (fn x => arrayOfInts(x, "entity")) lowerObjects, [])
	  val planeDef =  Option.getOpt(Option.map (fn x => arrayOfInts(x, "plane")) planePoints, List.tabulate(subDim + 1, fn x => IntInf.fromInt 0))
	  (*points to define a subdimensional object is that dim*)
	  val _ = if List.length planeDef = subDim + 1 (*TODO: get indexs hight.*)
		  then ()
		  else TypeError.error(cxt, [S"Unable to process ", S(name), S (" entry" ^(Int.toString index)),
					     S" because plane has wrong number of points."])


	 in
	  (index, entities, planeDef)
	 end
	 handle exn => raise exn

     fun parseHigher(subDim, SOME(json)) =
	 let
	  val _ = print("Ummm:"^(Int.toString subDim))
	  val listOfJson =  Vector.foldr (op ::) [] (JU.asArray(json))
	  val count = List.length listOfJson

	  val results = ListPair.map (fn (x,y) => parserHigherEntry(subDim, x, y)) (listOfJson, List.tabulate(count, fn x => x))
	 in
	  Higher(subDim, results)
	 end
       | parseHigher (subDim, NONE) = Higher(subDim, [])

     val higher = List.map parseHigher higherObjects

     fun mappingInserts(json) =
	 let
	  val possibleFields = List.tabulate(dim - 1, fn x => ("map"^(Int.toString (x+1)), x))
	  val grabField = fn x => (Option.map JU.asString) (JU.findField json x)
	  val possible = (List.filter (fn (x,y) =>Option.isSome y)) (List.map (fn (f,i) => (i, grabField f)) possibleFields)
	  val possible' = List.map (fn (x,SOME(y)) => Mapping(x,y)) possible
	 in
	  possible'
	  handle exn => []
	 end
     val maps = Option.getOpt((Option.map mappingInserts json),[])
	   (*do basis version*)
			     (*tabulate over it...*)

     val parameterized = if dim = 2
			 then
			  (case List.find(fn Higher(1, _) => true | _ => false) higher
			    of SOME(h) => analyzeGeometry1(cxt, verts, h)
			     | _ => raise Fail "case not covered but insert C")
			 else if dim = 3
			 then (case List.find (fn Higher(2, _) => true | _ => false) higher
				of SOME(h) => analyzeGeometry2(verts, h)
				 | _ => raise Fail "case not covered but insert C")
			 else raise Fail "invalid dim"
     val result = parameterized::verts::(List.@(higher, maps))
    in
     result
    end

      
      
fun paresRefCell(env, cxt, refCellJson, dim, machinePres, meshName) =
    let
     val tyOption = (findField "type") refCellJson
     val epsilonOption = (findField "epsilon") refCellJson
     val epsilon = Option.mapPartial (optionRealLit) epsilonOption
     val eps = (case epsilon
		 of SOME(e) => e
		  | _ => makeParseError(err,cxt, "epsilon", "field doesn't exist or isn't a real number.",  realToRealLit machinePres))
     val cellClass = Option.mapPartial (oE JU.asString) tyOption
     val cell = Option.mapPartial (fn  x => FemData.fromStr(x, dim)) cellClass

     val newtonParamsJson = (findField "newtonParams") refCellJson
     val newtonParamsDict = (case newtonParamsJson
			      of NONE => {contraction = true, itters = 16, killAfterTol = true, newtonTol = realToRealLit 0.001}
			       | SOME(json) =>
				 let
				  val contraction = Option.getOpt(Option.mapPartial (optionBool) (findField "contraction" json), true)
				  val _ = print("Contraction test:" ^ (Bool.toString contraction) ^ "\n")
				  val itters = Option.getOpt(Option.mapPartial (optionInt) (findField "itters" json), 16)
				  val killAfterTol = Option.getOpt(Option.mapPartial (optionBool) (findField "killAfterTol" json), true)
				  val newtonTol = Option.getOpt(Option.mapPartial (optionRealLit) (findField "newtonTol" json), realToRealLit 0.001)

				 in
				  {contraction=contraction, itters=itters, killAfterTol=killAfterTol, newtonTol=newtonTol}
				 end
			    (*end case*))
     val _ = print("epsilon: " ^ (RealLit.toString eps) ^ "\n");
     val geometry = parseGeometryOfRefCell(env, cxt, dim, SOME(refCellJson), meshName)
    in
     (case cell
       of SOME(cellVal) => SOME(FemData.RefCellData({ty=cellVal, eps = eps, newtonControl = newtonParamsDict}), geometry)
	| NONE => makeParseError(err, cxt, "cell:type", "field doesn't exsist, isn't a string, or string is stupid", NONE)
				
     (*end case*))
    end


    
      
fun parseAccelerate(env, cxt, json, meshName) =
    let
     val acc = findField "accelerate" json

     val insert = warnIfNo(env, cxt, "insert", "accelerate", meshName) acc
     val includes = warnIfNo(env, cxt, "includes", "accelerate", meshName) acc
     val linkDirs =  warnIfNo(env, cxt, "linkDirs", "accelerate", meshName) acc
     val includeDirs =  warnIfNo(env, cxt, "includeDirs", "accelerate", meshName) acc
     val libs =  warnIfNo(env, cxt, "libs", "accelerate", meshName) acc;
     val conservative =  warnIfNo(env, cxt, "conservative", "accelerate", meshName) acc

     fun getStrings(json) = (case json
			      of SOME(json') => JU.arrayMap (JU.asString) json'			
			       | NONE => []
			    )
			    handle exn => []
     val includes' = getStrings(includes)
     val linkDirs' = getStrings(linkDirs)
     val includeDirs' = getStrings(includeDirs)
     val libs' = getStrings(libs)
     (*store properties*)
     val newProp = Properties.NeedsExtraLibs(includes', includeDirs', libs', linkDirs')
     val _ = Env.recordProp(env, newProp);
     (*check if the env exists*)
     val parseInsert = Option.map JU.asString insert handle exn => NONE
     val conservative' = (JU.asBool) (Option.valOf conservative) handle exn => false
    in
     (case parseInsert
       of SOME(file) => if OS.FileSys.access(file, [])
			then SOME(Atom.atom file, conservative')
			else
			 (TypeError.warning(cxt, [S"File ", S(file), S" does not exist."]);
			  NONE))
    end
      
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

     val refCellField = Option.mapPartial (findField "refCell") mesh
     val refCellData = Option.mapPartial (fn (x,y) => paresRefCell(env, cxt,x,y, 0.0000001, tyName))
					 (optionTuple(refCellField, dim))

     val dimMethods = Option.getOpt(Option.map (fn (x,y) => y) refCellData, [])

     val optionalAcceleration = Option.mapPartial (fn  x => parseAccelerate(env, cxt, x, tyName)) mesh
     val meshVal = Option.map (fn (x,y,z,basis, cell) => FT.mkMesh(x,y,z,
							     basis,
							     tyName, cell, optionalAcceleration))
			      (case (dim, spaceDim, degree, parsedBasis, refCellData)
				of (SOME(a), SOME(b), SOME(c), SOME(d), SOME(e, _)) => SOME((a,b,c,d,e))
				 | _ => NONE (* end case*))
     val newConstants' = Option.getOpt(newConstants, [])
    in
     (meshVal, newConstants', dimMethods)
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
