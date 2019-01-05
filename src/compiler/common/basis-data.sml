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

	   val empty : t
	   val same : t * t -> bool
	   val hash : t -> word

	   val makeBasis : real Array.array * int * int -> t

	   val dx : t * int -> t
			     
			 
	  end = struct
structure R = Real
type r = Real.real
structure A = Array


datatype t = BasisFunc of {
	  dim : int,
	  degree : int,
	  strForm : Atom.atom,
	  coeffs : r Array.array 
	 }
val empty = BasisFunc {dim = 0, degree = 0, coeffs = A.fromList [], strForm = Atom.atom "0.0"}

fun toString(BasisFunc{strForm, ...}) = strForm

fun same(b1,b2) = Atom.same(toString(b1), toString(b2))
fun hash(BasisFunc{strForm,  ...}) = Atom.hash(strForm)


					      

(* we use an index scheme for our polynomials*)
(* We store them in the monomial basis where our array is degree^dim and indexed as such*)


(* translates an index x^a y^b z^c -> index: recall how multi-dimensional indeicies work: https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/*)
(* why I can't find this in the interger spec is beyond me *)
fun power (x, 0) = 1  
  | power (x, n) = x * power(x, n-1)
fun indexScheme(index, dim, degree) =
    let

     val degreePowers = List.tabulate(dim, (fn x => power(degree, x)))
     val intemediate = ListPair.zip (index, degreePowers)
     val index = let val ps = List.map (fn (x,y) => x*y) intemediate
			   in List.foldr (fn (x,y) => x+y) 0 ps end
    in
     index
    end


      
fun makeIndexes(dim, degree) =
    let
     val modIndexes = List.rev (List.tabulate(degree, fn x => power(degree, x)))
     fun idxes x = List.map (fn d => Int.mod(x, d)) modIndexes
     val tab = List.tabulate(power(degree, dim), fn x => x)
    in
     List.map idxes tab
    end

fun makeMonoTerms(reals, dim, degree) =
    let
     val vars = List.tabulate(dim, (fn x => "x"^Int.toString(x)))
     val idxes = makeIndexes(dim, degree)
     fun makePowers idx =
	 let
	  val powers = List.map (fn (x,y) => x ^ (Int.toString(y))) (ListPair.zip (vars, idx))
	  val product = String.concatWith "*" powers
	 in
	  product
	 end
     fun getIndex idx =
	 let
	  val test = A.sub(reals, indexScheme(idx, dim, degree))
			  
	 in
	  if Real.==(test, 0.0)
	  then ""
	  else Real.toString(test) ^ "*" ^ makePowers(idx)
	 end
     val monoTerms = List.filter (fn x => x <> "") (List.map getIndex idxes)
     val string = String.concatWith "+" monoTerms
     
    in
     string
    end
      (*1d, 2d, and 3d versions*)
      (*1d, 2d, 3d index*)
      (**)

fun makeBasis(reals, dim, degree) =
    let
     val monoRep = Atom.atom (makeMonoTerms(reals, dim, degree))
    in
     BasisFunc({dim = dim, degree = degree, coeffs = reals, strForm = monoRep})
    end

fun zero dim = BasisFunc({dim = dim, degree = 0, coeffs = A.fromList [0.0], strForm = Atom.atom "0.0"})
      
fun dx(BasisFunc({dim, degree, coeffs, strForm}), varIndex) =
    let
     val len = A.length coeffs
     val len' = power(Int.min(degree - 1, 0), dim)
     (* val primeArray = A.array(, real) *)
     (* val newArray = A.copy({src = coeffs, dst = primeArray, di = 0 }) *)
     val modIndexes = List.rev (List.tabulate(degree, fn x => power(degree, x)))
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
     then zero dim
     else 
      BasisFunc({dim = dim, degree = degree - 1, coeffs = coeffs, strForm = strForm })
    end


    

end
