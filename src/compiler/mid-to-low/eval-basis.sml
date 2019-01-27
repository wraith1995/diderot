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

structure IR = LowIR
structure BD = BasisData
structure Op = LowOps
structure AN = ArrayNd
structure Ty = LowTypes
structure BDA = BasisDataArray

datatype EvalMethod = Mono

fun consFunc(vars, ty) = (case ty
			   of Ty.TensorTy[] => (case vars
					      of [v] => IR.VAR(v)
					       | _ => raise Fail "impossible")
			    | Ty.TensorTy(_) => IR.CONS(vars, ty)
			    | _ => raise Fail "impossible")
fun multFunc(v1, v2) = (case (IR.Var.ty(v1), IR.Var.ty(v1))
			 of (Ty.TensorTy[], Ty.TensorTy[]) =>  IR.OP(Op.RMul, [v1,v2])
			  | (Ty.TensorTy[d1], Ty.TensorTy[d2]) => if d1=d2
								  then IR.OP(Op.VDot(d1), [v1,v2])
								  else raise Fail "impossible"
			  | _ => raise Fail "impossible"
		       )


fun power(avail, var, 0) =
    let
     val one = IR.LIT(Literal.Real(RealLit.one))
    in
     AvailRHS.addAssign(avail, "one", Ty.realTy, one)
    end
  | power (avail, var, 1) = AvailRHS.addAssign(avail, "one", Ty.realTy, IR.VAR(var))
  | power(avail, var, n) =
    let
     val prev = power(avail, var, n - 1)
     val prod = IR.OP(Op.RMul, [prev, var])
    in
     AvailRHS.addAssign(avail, "prod"^(Int.toString n), Ty.realTy, prod)
    end

fun multList(avail, [x]) = x
  | multList(avail, x::xs) =
    let
     val y = multList(avail, xs)
     val prod = IR.OP(Op.RMul, [x,y])
     val prodVar = AvailRHS.addAssign(avail, "prod", Ty.realTy, prod)
    in
     prodVar
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

fun evalBasisDumb(avail, basisArray, pos) =
    let
     val (basisArrayData, meta) = BDA.explode basisArray
     val varAcc = BDA.vars(meta)
     val vars = extractVars(avail, pos, varAcc)
     fun convert(func) = evalFunctionDumb(avail, func, vars)
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
