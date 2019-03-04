(* eval-basis.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *)

structure EvalBasis : sig

	   val expand : LowIR.var * BasisDataArray.t * LowIR.var -> (LowIR.var * LowIR.rhs) list
	  end  = struct
val dumb = true
structure IR = LowIR
structure V = IR.Var
structure BD = BasisData
structure Op = LowOps
structure AN = ArrayNd
structure Ty = LowTypes
structure BDA = BasisDataArray
datatype EvalMethod = Mono

fun getLiteral x = (case V.getDef x
		    of IR.LIT arg => SOME arg
		     | _ => NONE
		  (* end case *))			

fun isOne x = (case getLiteral x
		of NONE => false
		 | SOME((Literal.Real(r))) => RealLit.same(r, RealLit.one)
		 | SOME((Literal.Int(i))) => IntLit.same(i, IntLit.fromInt 1)
		 | _ => false
	      (*end case*))
fun isNegOne x = (case getLiteral x
		   of NONE => false
		    | SOME((Literal.Real(r))) => RealLit.same(r, RealLit.m_one)
		    | SOME((Literal.Int(i))) => IntLit.same(i, IntLit.fromInt (~1))
		    | _ => false
		 (*end case*))

fun isZero x = (case getLiteral x
		 of SOME((Literal.Real(r))) => RealLit.same(r, RealLit.zero true)
						     orelse
						     RealLit.same(r, RealLit.zero false)
		  | SOME((Literal.Int(i))) => IntLit.same(i, IntLit.fromInt 0)
		  | _ => false
	       (*end case*))

fun consFunc(vars, ty) = (case ty
			   of Ty.TensorTy[] => (case vars
					      of [v] => IR.VAR(v)
					       | _ => raise Fail "impossible")
			    | Ty.TensorTy[1] => raise Fail "impossible"
			    | Ty.TensorTy(_) => IR.CONS(vars, ty)
			    | _ => raise Fail "impossible")
fun multFunc(v1, v2) = (case (IR.Var.ty(v1), IR.Var.ty(v1))
			 of (Ty.TensorTy[], Ty.TensorTy[]) =>  IR.OP(Op.RMul, [v1,v2])
			  | (Ty.TensorTy[d1], Ty.TensorTy[d2]) =>
			    if d1=d2 andalso (d1<>1) andalso (d2 <> 1)
			    then IR.OP(Op.VDot(d1), [v1,v2])
			    else raise Fail "impossible"

			  | _ => raise Fail "impossible"
		       )
(*USE THIS UNTIL THE VECTORIZATION IS FIXED:*)
fun sumVars(avail, [v]) = AvailRHS.addAssign(avail, "sum", Ty.realTy, IR.VAR(v))
  | sumVars(avail, [v1, v2]) = AvailRHS.addAssign(avail, "sum", Ty.realTy, IR.OP(Op.RAdd, [v1,v2]))
  | sumVars(avail, v::vs) =
    let
     val sum = sumVars(avail, vs)
    in
     AvailRHS.addAssign(avail, "sum", Ty.realTy, IR.OP(Op.RAdd, [v, sum]))
    end


fun power(avail, var, 0) =
    let
     val one = IR.LIT(Literal.Real(RealLit.one))
    in
     AvailRHS.addAssign(avail, "one", Ty.realTy, one)
    end
  | power (avail, var, 1) = AvailRHS.addAssign(avail, "iden", Ty.realTy, IR.VAR(var))
  | power(avail, var, n) =
    if isZero var
    then var
    else if isOne var
    then var
    else if isNegOne var
    then (if Int.mod(n, 2) = 1 then var (*neg as 2n+1*)
	  else power(avail, var, 0) (*pos as 2n*))
    else
    let
     val prev = power(avail, var, n - 1)
     val prod = IR.OP(Op.RMul, [prev, var])
    in
     AvailRHS.addAssign(avail, "prod"^(Int.toString n), Ty.realTy, prod)
    end
(*filter ones*)
(*filter *)
fun multList(avail, xs) =
    let
     fun multList'(avail, true, [x]) = x
       | multList' (avail, false, [x]) =
	 let
	  val neg = IR.OP(Op.RNeg, [x])
	  val negVar = AvailRHS.addAssign(avail, "neg", Ty.realTy, neg)
	 in
	  negVar
	 end
       | multList'(avail, pos, x::xs) =
	 if isZero x
	 then x
	 else if isOne x
	 then multList'(avail, pos, xs)
	 else if isNegOne x
	 then multList'(avail, Bool.not(pos), xs)
	 else
	  let
	   val y = multList'(avail, pos, xs)
	   val prod = IR.OP(Op.RMul, [x,y])
	   val prodVar = AvailRHS.addAssign(avail, "prod", Ty.realTy, prod)
	  in
	   prodVar
	  end
    in
     multList'(avail, true, xs)
    end
fun extractVars(avail, var, idxes) =
    let
     val d = List.length idxes
     val accessOp = List.map ( fn x => Op.VIndex(d,x)) idxes
    in
     (*need to use TensorIndex or something else like that*)
     List.map
       (fn x => AvailRHS.addAssign(avail, "varAcc", Ty.realTy, IR.OP(x, [var])))
       accessOp
    end

      
fun powerIndex(avail, vars, idx) =
    let
     val powers = List.map (fn (x,y) => power(avail, x,y)) (ListPair.zip(vars, idx))
    in
     multList(avail, powers)
    end

(*We need this function because lack of choice about RealLit and real at various points*)
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

(*See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.330.7430&rep=rep1&type=pdf algorithim 1*)
fun countOccurance(x : int, nonZeroTerms : (real * int list) list) : int=
    List.foldr (op+) 0 (List.map (fn (_, idxes) => List.nth(idxes,x)) nonZeroTerms)

fun partitionById(id : int, nonZeroTerms) =
    let
     val has = List.filter (fn (z, x) => List.nth(x, id) > 0) nonZeroTerms
     val notHas = List.filter (fn (z, x)=> List.nth(x, id) = 0) nonZeroTerms
    in
     (has, notHas)
    end

fun decId(l, ls) =
    let
     fun decId'(y, c, x::xs) = if y = c
			       then (x-1)::xs
			       else x::decId'(y, c + 1, xs)
       | decId'(y, c, []) = raise Fail "imossible"
					 
    in
     decId'(l, 0, ls)
    end

fun constants(ls) = List.filter (fn (_,x) => List.all (fn z => z=0) x) ls
fun notConstants(ls) = List.filter (fn (_,x) => List.exists (fn z => z>0) x) ls

fun evalFunctionLazyHorner'(avail, [(r,c)], vars) =
    let
     val pows = powerIndex(avail, vars, c)
     val coeff =  AvailRHS.addAssign(avail, "extra", Ty.realTy, IR.LIT(Literal.Real(realToRealLit(r))))
     val result = multList(avail, [pows, coeff])
    in
     result
    end
  | evalFunctionLazyHorner'(avail, [], vars) = AvailRHS.addAssign(avail, "extra", Ty.realTy, IR.LIT(Literal.Real(RealLit.zero true)))
  | evalFunctionLazyHorner'(avail, nonZeroTerms , vars)  = 
    let
     (*check if count is 0 or 1*)
     val varCount = List.length vars
     val _ = print("vars:"^(Int.toString varCount) ^ "\n")
     val idxTerms = List.map (fn (x,y) => y) nonZeroTerms
     val _ = print("Idxes len:"^(String.concatWith "," (List.map (Int.toString o List.length) idxTerms))^"\n")
     val varCounts = List.tabulate(varCount, fn x => (x, countOccurance(x, nonZeroTerms))) handle exn => raise exn
     val maxPair = (List.foldr (fn ((newId,newCount), (oldId, oldCount)) =>
				  if newCount > oldCount
				  then (newId, newCount)
				  else (oldId, oldCount))) (~1,~1) varCounts
	
     val (maxId, _) = maxPair
     val _ = if maxId = ~1
	     then raise Fail "ops"
	     else ()
     val maxIdVar = List.nth(vars, maxId)

     val (hasId, notHasId) = partitionById(maxId, nonZeroTerms)
					  
     val hasIdDec = List.map (fn (r,x) => (r,decId(maxId, x))) hasId
			     
     val constants = constants(notHasId)
     val nonConstantNotId = notConstants(notHasId)
     (*I'm not sure it is possible or constants to have len > 1.*)
     val _ = if List.length constants > 1
     	     then print("odd constants somehow\n")
	     else ()
     val constants' = List.map (fn (r, _) => r) constants
     val extra = List.foldr (Real.+) 0.0 constants' 
     val extraVar = AvailRHS.addAssign(avail, "extra", Ty.realTy, IR.LIT(Literal.Real(realToRealLit(extra))))
     (*check the sizes here before doing this:*)
     val a0 = evalFunctionLazyHorner'(avail, nonConstantNotId, vars)
     val a1 = evalFunctionLazyHorner'(avail, hasIdDec, vars)
     val xa1 = AvailRHS.addAssign(avail, "xa1", Ty.realTy, IR.OP(Op.RMul, [maxIdVar, a1]))
     val sum = if List.length nonConstantNotId > 0
	       then AvailRHS.addAssign(avail, "a0pxa1", Ty.realTy, IR.OP(Op.RAdd, [a0, xa1]))
	       else xa1
     val sum' = if List.length constants > 0
		then AvailRHS.addAssign(avail, "sa0pxa1", Ty.realTy, IR.OP(Op.RAdd, [extraVar, sum]))
		else sum
    in
     sum'
    end



fun evalFunctionDumb(avail, basisFunc, vars) =
    let
     val nonZeroTerms = BD.nonZeroTerms basisFunc
     val count = List.length nonZeroTerms
     fun foldrFunc((x,y), (lst1, lst2)) =
	 let
	  val realLit = IR.LIT(Literal.Real(realToRealLit(x)))
	  val assign = AvailRHS.addAssign(avail, "coeff", Ty.realTy, realLit)
	  val pow = powerIndex(avail, vars, y)				 
	 in
	  (assign ::lst1, pow :: lst2)
	 end
     val result =
	 if count = 0
	 then AvailRHS.addAssign(avail, "basisEval", Ty.realTy, IR.LIT(Literal.Real(RealLit.zero(false))))
	 else if count = 1
	 then
	  let
	   val ([c], [p]) = List.foldr foldrFunc ([],[]) nonZeroTerms
	   val resultOP = multFunc(c, p)
	   val next = AvailRHS.addAssign(avail, "basisEval", Ty.realTy, resultOP)
	  in
	   next
	  end
	 else
	  if dumb
	  then
	   let
	    val (coeffs, pows) = List.foldr foldrFunc ([],[]) nonZeroTerms
	    val combined = ListPair.zip (coeffs, pows)
	    val ops = List.map multFunc combined
	    val mults = List.map (fn x => AvailRHS.addAssign(avail, "mult", Ty.realTy, x)) ops
	    val sum = sumVars(avail, mults)
	    val next = IR.Var.new("intermediate", Ty.realTy)
	   in
	    (AvailRHS.addAssignToList(avail, (next, IR.VAR(sum))); next)
	   end
	  else
	  let
	   val (coeffs, pows) = List.foldr foldrFunc ([],[]) nonZeroTerms
					   
	   val coeffVar = AvailRHS.addAssign(avail, "coeffs", Ty.vecTy(count), consFunc(coeffs, Ty.vecTy(count)))
	   val powVar = AvailRHS.addAssign(avail, "powers", Ty.vecTy(count), consFunc(pows, Ty.vecTy(count)))
	   val resultOP = multFunc(coeffVar, powVar) (*vectorization*)
	   val next = AvailRHS.addAssign(avail, "basisEval", Ty.realTy, resultOP)
	   val next' = IR.Var.new("intermediate", Ty.realTy)
	  in
	   (AvailRHS.addAssignToList(avail, (next', IR.VAR(next))); next')
	  end
	    
    in
     result
    end
fun evalFunctionLazyHorner(avail, basisFunc, vars) =
    let
     val nonZeroTerms = BD.nonZeroTerms basisFunc
    in
      evalFunctionLazyHorner'(avail, nonZeroTerms, vars)
    end

fun evalBasisDumb(avail, basisArray, pos) =
    let
     val affineTest = BDA.isAffine basisArray
     val (basisArrayData, meta) = BDA.explode basisArray
     val varAcc = BDA.vars(meta)
     val vars = extractVars(avail, pos, varAcc)
     fun convert(func) = if true (*affineTest*)
			 then evalFunctionDumb(avail, func, vars)
			 else evalFunctionLazyHorner(avail, func, vars)
     fun group(vars, n) =
	 let
	  val v::vs = vars
	  val ty = IR.Var.ty(v)
	  val Ty.TensorTy(shape) = ty
	  val cons = consFunc(vars, Ty.TensorTy(n::shape))
	  
			    
	 in
	  AvailRHS.addAssign(avail, "intermediateCons", Ty.TensorTy(n::shape), cons)
	 end

    in
     AN.convertToTree(basisArrayData, convert, group)
    end


fun expand(lhs, data, pos) =
    let
     val avail = AvailRHS.new()
     val int = (case BDA.analyzeBasis data
		 of Unknown => evalBasisDumb(avail, data, pos))
     val newVar = IR.VAR(int)
     val _ = AvailRHS.addAssignToList(avail, (lhs, newVar))
			
    in
     List.rev (AvailRHS.getAssignments(avail))
    end

end
