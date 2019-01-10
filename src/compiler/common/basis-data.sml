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

	   val makeBasis : real Array.array * int * int * int -> t option

	   val dx : t * int -> t
			     
			 
	  end = struct

structure R = Real
type r = Real.real
structure A = Array

(* why I can't find this in the interger spec is beyond me *)
fun power (x, 0) = 1  
  | power (x, n) = x * power(x, n-1)		

datatype t = BasisFunc of {
	  dim : int,
	  degree : int,
	  mapDim : int,
	  strForm : Atom.atom,
	  coeffs : r Array.array 
	 }
fun empty(dim, degree) = BasisFunc {dim = dim,
				    degree = degree,
				    mapDim = 0,
				    coeffs = A.tabulate(power(degree, dim), fn x => 0.0),
				    strForm = Atom.atom "0.0"}
fun toString(BasisFunc{strForm, ...}) = Atom.toString strForm
fun getString(BasisFunc{strForm, ...}) = strForm
fun same(b1,b2) = Atom.same(getString(b1), getString(b2))
fun hash(BasisFunc{strForm,  ...}) = Atom.hash(strForm)

(* The following bits of code are all dedicated to our storage and manipulation of polynomials.
   The idea is that we will store them in the simplest possible way, easily extensible to all dimensions and degrees and so on.
   If you, have three paramters: 
   dim -> the number of variables
   degree -> the maximal power that any variable might be taken too
   mapDim -> the exact number of coefficients needed to specify a polynomial.

  With this, we can store the whole thing in an array of size (degree +1 )^dim; we treat this an array of size [degree][degree]...
  Thus, for a polynomial, coeff * x^a * y^b * z^c -> array[a][b][c] = coeff.
  We observe that this will allow us to do things like use a basis of only quadratic terms or other weird shit like that.
  There is some waste, but we are going to keep things clean.
  To see how C-style N-d arrays work, refer to: https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/

  Based on this sytem, we will take derivatives, tensor products, sums, direct sums, and produce string representations. This is all that we need to do. 

*)


					      
(*Given an index [a,b,c,d,..] of length dim, we translate to a single integer index so that array[idx] = array[a][b]...*)
fun indexScheme(index, dim, degree) =
    let

     val degreePowers = List.tabulate(dim, (fn x => power(degree + 1, x))) (* compute (degree+1)^0,..., (degree + 1)^(dim - 1)*)
     val intemediate = ListPair.zip (index, degreePowers) (*pair these with a member of the index*)
     val index' = List.foldr (fn ((a,b),y) => a*b+y) 0 intemediate (* compute the desired sum ... *)
    in
     index'
    end
(* fun printIdxes([]) = "" *)
(*   | printIdxes (x::xs) = "["^ (String.concatWith "," (List.map (Int.toString) x)) ^"] " ^ printIdxes(xs) *)
      
(* Given a single integer index, idx, to an N-d array of size (degree+1)^dim, we wish to compute the index [a,b,c,d,...] so indexScheme([a,b,c,d,...], dim, degree)=idx *)
(* This function builds a list so that list[idx] = [a,b,c,d,...]*)
(* We do this in the following manner: 
   We wish to find [a,b,c,d...,last] so that a*(degree+1)^(dim - 1) + b (degree+1)^(dim - 2) + c (degree+1)^(dim - 3) + ... last = idx
   If I take the mod and div of idx by (degree+1)^(dim - 1), we get a*(degree+1)^(dim - 1) + rem = idx.
   If we apply the same to rem, we can get b. So we just continue on like this.
 *)
fun makeIndexes(dim, degree) =
    let
     val modIndexes = List.rev (List.tabulate(dim, fn x => power(degree + 1, x)))

     fun modDiv(num, divider) = (Int.mod(num,divider), Int.div(num, divider))
     fun invert'(num, x::xs, ys) = if x = 1
				   then num::ys
				   else let val (mod',div') = modDiv(num, x)
					in invert'(mod', xs, div'::ys) end

     fun idxes x = invert'(x, modIndexes, [])
     val tab = List.tabulate(power(degree + 1, dim), fn x => idxes x)
    in
     tab
    end

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
     val idxes = makeIndexes(dim, degree)

     fun makePowers idx = 
	 let
	  val powers = List.map (fn (x,y) => x ^ "^" ^ (Int.toString(y))) (ListPair.zip (vars, idx))
	  val product = String.concatWith " * " powers
	 in
	  product
	 end
     fun getIndex idx =
	 let
	  val test = A.sub(reals, indexScheme(idx, dim, degree))
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
      
fun makeBasis(reals, dim, degree, mapDim) =
    let
     val monoRep = Atom.atom (makeMonoTerms(reals, dim, degree))
     val _ = print("Imported basis func with mono: " ^ (Atom.toString monoRep) ^ "\n")
    in
     SOME(BasisFunc({dim = dim, degree = degree, mapDim = mapDim, coeffs = reals, strForm = monoRep}))
     handle exn => (print(exnMessage(exn)); NONE)
    end

fun zero dim mapDim = BasisFunc({dim = dim, degree = 0,mapDim = mapDim, coeffs = A.fromList [0.0], strForm = Atom.atom "0.0"})
      
fun dx(BasisFunc({dim, degree, coeffs, strForm, mapDim,...}), varIndex) =
    let
     val len = A.length coeffs
     val len' = power(Int.min(degree - 1, 0), dim)
     (* val primeArray = A.array(, real) *)
     (* val newArray = A.copy({src = coeffs, dst = primeArray, di = 0 }) *)
     val modIndexes = List.rev (List.tabulate(degree, fn x => power(degree + 1, x)))
     fun idxes x = Array.fromList (List.map (fn d => Int.mod(x, d))  modIndexes)
     fun getVarIdx(idx : int A.array) : int = A.sub( idx, varIndex)
     fun setIdx idx =
	 let
	  val index = idxes idx (* translate array index to triplet index*)
	  val power : int = getVarIdx index
				      (* get power required *)
				      (*create the lower index*)
				      (* set *)
	 in

	  ()
	 end
    in
     if dim - varIndex < 0
     then raise Fail "impossible derivative"
     else if degree = 1
     then zero dim 1
     else 
      BasisFunc({dim = dim, degree = degree - 1,mapDim = mapDim, coeffs = coeffs, strForm = strForm })
    end


    

end
