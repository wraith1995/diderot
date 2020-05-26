(* basis-data.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * Information about a basis for a function space.
*)


structure BasisData : sig
	   type t

	   val empty : int * int -> t
	   val same : t * t -> bool
	   val hash : t -> word
	   val toString : t -> string

	   val isMono : t -> bool
	   val isAffine : t -> bool
	   val extractMonoTerm : t * int list -> real
	   val nonZeroTerms : t -> (real * int list) list
	   val dim : t -> int
	   val degree : t -> int


	   val makeBasisFunc : real Array.array * int * int -> t option
	   val zero : int -> t

	   val dx : t * int -> t
	   (*not implemented:*)
	   val scale : t * real -> t
	   val sum : t * t -> t
	   val cast : t * int -> t
	   val product : t* t -> t

			     
			 
	  end = struct

structure R = Real
type r = Real.real
structure A = ArrayNd

datatype Rep = Array of r A.ArrayNd (*to add sparse representation later*)

(* why I can't find this in the interger spec is beyond me *)
fun power (x, 0) = 1  
  | power (x, n) = x * power(x, n-1)		

fun arrayToList arr = Array.foldr (op ::) [] arr
fun dimDegreeShape(dim, degree) = List.tabulate(dim, fn _ => degree + 1)
					       
datatype t = BasisFunc of {
	  dim : int, (* number of vars*)
	  degree : int, (*height sum of powers*)
	  mapDim : int, (* number of coefficiens used*)
	  strForm : Atom.atom, (* string representation used*)
	  coeffs : Rep (*storage of coefficients*)
	 }
fun dimDomain(BasisFunc{dim,...}) = dim
fun empty(dim, degree) = BasisFunc {dim = dim,
				    degree = degree,
				    mapDim = 0,
				    coeffs = Array (A.fromArray'(Array.tabulate(power(degree, dim), fn x => 0.0), dimDegreeShape(dim,degree))),
				    strForm = Atom.atom "0.0"}
fun toString(BasisFunc{strForm, ...}) = Atom.toString strForm
fun getString(BasisFunc{strForm, ...}) = strForm
fun same(b1,b2) = Atom.same(getString(b1), getString(b2))
fun hash(BasisFunc{strForm,  ...}) = Atom.hash(strForm)


fun numNonZero(Array(coeffs')) = let
 fun nonZero x = if Real.==(x, 0.0)
		 then 0
		 else 1 
			
in
 (ArrayNd.foldr (fn (x,y) => y + (nonZero x)) 0 coeffs')
end
					      
fun isMono(BasisFunc{coeffs=Array(coeffs'),...}) = let
 fun nonZero x = if Real.==(x, 0.0)
		 then 0
		 else 1 
			
in
 (ArrayNd.foldr (fn (x,y) => y + (nonZero x)) 0 coeffs') <= 1
end

fun isAffine(BasisFunc{coeffs=Array(coeffs'), ...}) =
    let
     fun zeroIfNonAffine (idx, x, y) =
	 if y
	 then if Real.==(x, 0.0)
	      then true
	      else List.all (fn x => x <= 1) idx
	 else false
						  
				      
    in
     ArrayNd.foldri' zeroIfNonAffine true coeffs'
    end


fun extractMonoTerm(BasisFunc{coeffs=Array(coeffs'),...}, idx) = if A.indexInside(coeffs', idx)
								 then A.sub'(coeffs', idx)
								 else 0.0
fun nonZeroTerms(BasisFunc{coeffs=Array(coeffs'),...}) =
    let
     fun foldrFunc(inlist, a, lst) =
	 if Real.==(a,0.0)
	 then lst
	 else (a,inlist) :: lst
    in
     ArrayNd.foldri' foldrFunc [] coeffs'
    end
fun dim(BasisFunc{dim,...}) = dim
fun degree(BasisFunc{degree,...}) = degree

(* The following bits of code are all dedicated to our storage and manipulation of polynomials.
   The idea is that we will store them in the simplest possible way, easily extensible to all dimensions and degrees and so on.
   If you, have three paramters: 
   dim -> the number of variables
   degree -> the maximal power that any variable might be taken too
   mapDim -> the exact number of coefficients needed to specify the polynomial.

  With this, we can store the whole thing in an array of size (degree +1 )^dim; we treat this an array of size [degree][degree]...
  Thus, for a polynomial, coeff * x^a * y^b * z^c -> array[a][b][c] = coeff.
  We observe that this will allow us to do things like use a basis of only quadratic terms or other weird shit like that.
  There is some waste, but we are going to keep things clean.
  To see how C-style N-d arrays work, refer to: https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/

  Based on this sytem, we will take derivatives, tensor products, sums, direct sums, and produce string representations. This is all that we need to do. 

*)


					      
(*Given an index [a,b,c,d,..] of length dim, we translate to a single integer index so that array[idx] = array[a][b]...*)
(* fun printIdxes([]) = "" *)
(*   | printIdxes (x::xs) = "["^ (String.concatWith "," (List.map (Int.toString) x)) ^"] " ^ printIdxes(xs) *)
      
(* Given a single integer index, idx, to an N-d array of size (degree+1)^dim, we wish to compute the index [a,b,c,d,...] so indexScheme([a,b,c,d,...], dim, degree)=idx *)
(* This function builds a list so that list[idx] = [a,b,c,d,...]*)
(* We do this in the following manner: 
   We wish to find [a,b,c,d...,last] so that a*(degree+1)^(dim - 1) + b (degree+1)^(dim - 2) + c (degree+1)^(dim - 3) + ... last = idx
   If I take the mod and div of idx by (degree+1)^(dim - 1), we get a*(degree+1)^(dim - 1) + rem = idx.
   If we apply the same to rem, we can get b. So we just continue on like this.
 *)

(*In order to hash these, we need a consistent string rep of basis functions. 
  This function creates the monomial form of an arbitrary basis function.
  It uses the above machinery. And obeys the following rules:
  Vars are of the form x0, x1, ... x(dim-1)
  Terms are printed coeff * (x0)^a * (x1)^b * ... 
  The poly is printed term + term + term.
  The exception is the zero poly, which is just "0.0".
 *)

fun makeMonoTerms(reals, dim, degree) =
    let
     val vars = List.tabulate(dim, (fn x => "(x"^Int.toString(x)^")"))
     val idxes = arrayToList(A.getInverseIndex(reals))

     fun makePowers idx = 
	 let
	  val powers = List.map (fn (x,y) => x ^ "^" ^ (Int.toString(y))) (ListPair.zip (vars, idx))
	  val product = String.concatWith " * " powers
	 in
	  product
	 end
     fun getIndex idx =
	 let
	  val test = A.sub'(reals, idx)
		     handle exn => raise exn
	 in
	  if Real.==(test, 0.0)
	  then NONE
	  else SOME(Real.toString(test) ^ " * " ^ makePowers(idx))
	 end
     val monoTerms = List.map (Option.valOf) (List.filter Option.isSome (List.map getIndex idxes))
     val string = String.concatWith " + " monoTerms
    in
     string
    end

fun actualDegree(array) =
    ArrayNd.foldri' (fn (idx,x,y) => if Real.==(x,0.0)
				    then y
				    else Int.max(List.foldr (op+) 0 idx, y)) 0 array
      
fun makeBasisFunc(reals, dim, degree) =
    let
     val shape = dimDegreeShape(dim, degree)
     val reals' = A.fromArray'(reals, shape)
     val monoRep = Atom.atom (makeMonoTerms(reals', dim, degree))
     val degree' = actualDegree(reals')
    in
     SOME(BasisFunc({dim = dim, degree = degree', mapDim = numNonZero(Array(reals')), coeffs = Array(reals'), strForm = monoRep}))
     handle exn => (print(exnMessage(exn)); NONE)
    end

fun zero dim = BasisFunc({dim = dim, degree = 0,mapDim = 1, coeffs = Array(A.fromList [0.0]), strForm = Atom.atom "0.0"})
      
fun dx(t as BasisFunc({dim, degree, coeffs=Array(coeffs'), strForm, mapDim,...}), varIndex) =
    let
     val degree' = Int.max(degree - 1, 0)
     val newShape = dimDegreeShape(dim, degree')
     val start = List.tabulate(dim, fn _ => 0)
     val newArray = ArrayNd.array'(0.0, newShape)(* ArrayNd.subregion(coeffs', start, newShape) *)
     fun printList(a) = "[" ^ (String.concatWith "," (List.map Int.toString a)) ^ "]"
     fun allSameBut(idx1, idx2, i) = List.foldri (fn (j,(x,y), n) => n andalso (i=j orelse x=y))
						 true (ListPair.zip(idx1,idx2))
     fun modNew(idx, a) =
	 let
	  
	  val power = List.nth(idx, varIndex) + 1
	  val realPower = Real.fromInt power
	  fun foldrFunc(index, a', b') = if List.nth(index, varIndex) = power
					    andalso allSameBut(index,idx, varIndex)
				       then a'+b'
				       else b'
	  (* grab everyone with idx_{varIndex} = power*)
	 in
	   realPower * (ArrayNd.foldri' foldrFunc 0.0 coeffs')
	 end
     val _ = if  Atom.same(strForm, Atom.atom "0.0")
	     then  ()
	     else ArrayNd.modifyi' modNew newArray
     val strForm' = Atom.atom (makeMonoTerms(newArray, dim, degree'))
     val mapDim' = numNonZero(Array(newArray))
     (* val _ = print("We started with poly:"^(Atom.toString strForm)^"\n") *)
     (* val _ = print("We took the derivative:"^(Int.toString varIndex)^"\n") *)
     (* val _ = print("We got the poly:"^(Atom.toString strForm')^"\n") *)
     (* val _ = print("It has degree:"^(Int.toString degree')^"\n") *)

    in
     if degree = 0  (*fix bad checking of this.*)
     then zero dim
     else
      let
       val _ = ()
      in
       BasisFunc{dim = dim, degree = degree', coeffs = Array(newArray), mapDim = mapDim', strForm = strForm'}
      end
    end



fun scale(t,s) = t
fun sum(t1,t2) = t1
fun cast(t, i) = t
fun product(t1,t2) = t1


    

end
