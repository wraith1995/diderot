(* basis-data-array.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2018 The University of Chicago
 * All rights reserved.
 * 
 * This file exists for stupid reasons; the ir.spec needs arrays of basis, but the gen for the ir script doesn't like that. This is easier than fixing that.
*)

structure BasisDataArray : sig
	   type meta
	   val vars : meta -> int list
	   val dim : meta -> int
	   val maxDegree : meta -> int
	   val minDegree : meta -> int
				     
	   type t
	   (* basis and var map i.e first var here is first integer in the list*)
	   (*polys, uniform dim, var acceses, *)
	   val same : t * t -> bool
	   val hash : t -> word
	   val toString : t -> string
	   val isAffine : t -> bool
	   val explode : t -> BasisData.t ArrayNd.ArrayNd * meta
	   val makeUniform : BasisData.t ArrayNd.ArrayNd * int -> t
	   val makeVarSubset : BasisData.t ArrayNd.ArrayNd * int list -> t
	   val shape : t -> int list
	   val domainDim : t -> int

	   val d : t -> t
	   val D : t -> t

	   datatype basisInfo =  Unknown
	   val analyzeBasis : t -> basisInfo
			
				 
	  end = struct
type meta = {vars : int list, dim : int, maxDegree : int, minDegree : int}
fun vars({vars,...} : meta) = vars
fun dim({dim,...} : meta) = dim
fun maxDegree({maxDegree,...} : meta) = maxDegree
fun minDegree({minDegree,...} : meta) = minDegree
	      
fun sameMeta(m1, m2) = (vars(m1) = vars(m2)) andalso (dim(m1) = dim(m2))
		       andalso (maxDegree(m1) = maxDegree(m2))
		       andalso (minDegree(m1) = minDegree(m2))
fun hashMeta(m1) = (List.foldr (fn (d, s) => 0w5 * Word.fromInt d + s) 0w3 (vars(m1)))
		  + 0w7 * (Word.fromInt (maxDegree(m1)))
		  + 0w11 * (Word.fromInt (minDegree(m1)))
datatype t = Array of BasisData.t ArrayNd.ArrayNd * meta
(*basis, vars used, maxDegree*)
							
fun same(Array(a1, m1),Array(a2,m2)) = sameMeta(m1,m2) andalso ArrayNd.all (BasisData.same) (ArrayNd.zip(a1,a2))
fun hash(Array(a1, m1)) = (ArrayNd.foldr (fn (d, s) => 0w5 * BasisData.hash d + s) 0w3 a1)
			  + 0w7 * hashMeta(m1)
fun intStrings(xs) = String.concatWith "," (List.map (Int.toString) xs)
fun toString(Array(a1,_)) = "BasisArray(" ^ intStrings(ArrayNd.shape a1) ^")"

fun isAffine(Array(a1, m1)) = (maxDegree m1 <= 1) orelse (ArrayNd.all BasisData.isAffine a1)
						   

fun explode(Array(a1,m1)) = (a1, m1)
fun makeUniform(a1, n) =
    let
     val degrees = List.map BasisData.degree (ArrayNd.toList a1) (* dumb *)
     val max = List.foldr (fn (x,y) => Int.max(x,y)) 0 degrees
     val min = List.foldr (fn (x,y) => Int.min(x,y)) 0 degrees
     val vars = List.tabulate(n, fn x => x)
     val meta : meta = {vars=vars,minDegree=min,maxDegree=max, dim =n}
    in
     Array(a1, meta)
    end

fun makeVarSubset(a1, vars) =
    let
     val n = List.length vars
     val degrees = List.map BasisData.degree (ArrayNd.toList a1) (* dumb *)
     val max = List.foldr (fn (x,y) => Int.max(x,y)) 0 degrees
     val min = List.foldr (fn (x,y) => Int.min(x,y)) 0 degrees
     val meta : meta = {vars=vars,minDegree=min,maxDegree=max, dim =n}
    in
     Array(a1, meta)
    end      
fun shape(Array(a1,_)) = ArrayNd.shape a1
fun domainDim(Array(a1, meta)) = dim meta

fun dervMeta({vars, dim, maxDegree, minDegree}) = {vars=vars, dim=dim, maxDegree=Int.max(0, maxDegree - 1), minDegree = Int.max(0, minDegree - 1)}

fun d(Array(basisArray, meta)) =
    let
     val first = ArrayNd.sub(basisArray,0)
     val dim = List.length (vars(meta))
     val funcs = (List.tabulate(dim, fn z => (fn y => BasisData.dx(y, z))))
     fun expandMap a = Array.fromList (List.map (fn y => y a) funcs)
     val new = ArrayNd.expandMap expandMap dim basisArray
    in
     Array(new, dervMeta(meta))
    end
fun D(Array(basisArray, meta)) =
    let
     (*val first = ArrayNd.sub(basisArray,0) *)
     val dim = List.length (vars(meta))
     val preConcat = List.tabulate(dim, fn _ => basisArray)
     val new = ArrayNd.concat(preConcat)
     fun modFunc(idx, a) = BasisData.dx(a, List.nth(idx, 0))
     val _   = ArrayNd.modifyi' modFunc new

    in
     Array(new, dervMeta(meta))
    end

datatype basisInfo =  Unknown
fun analyzeBasis(a) = Unknown

end
