(* api-types.sml
 *
 * A representation of the types of values that can be communicated to and from a
 * Diderot program.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *)

structure APITypes =
  struct

    datatype t
      = IntTy
      | BoolTy
      | TensorTy of int list
      | StringTy
      | ImageTy of int * int list
      | FemData of FemData.femType
      | SeqTy of t * int option
      | TupleTy of t list

    val realTy = TensorTy[]

  (* does a type have a non-static size? *)
    fun hasDynamicSize StringTy = true
      | hasDynamicSize (ImageTy _) = true
      | hasDynamicSize (FemData _) = false (*TODO: Clarify the use of this function for this entry.*)
      | hasDynamicSize (SeqTy(_, NONE)) = true
      | hasDynamicSize (TupleTy(tys)) = (List.exists hasDynamicSize tys)
      | hasDynamicSize _ = false
			 
    fun hasFem (SeqTy(ty,_)) = hasFem ty
      | hasFem (FemData(_)) = true
      | hasFem (TupleTy(tys)) = not (List.all (not o hasFem) tys)
      | hasFem _ = false

    fun hasTupleEsq (SeqTy(ty, SOME(n))) = hasDynamicSize ty orelse hasTupleEsq ty (*has tuple or [][n]*)
      | hasTupleEsq (SeqTy(ty, NONE)) = hasTupleEsq ty
      | hasTupleEsq (TupleTy(tys)) = true
      | hasTupleEsq _ =  false


    fun depth ty =
	let
	 fun depth'(SeqTy(ty, SOME(s)), ds) = depth'(ty, s::ds)
	   | depth'(_, ds) = List.rev ds
	in
	 depth'(ty, [])
	end
    fun toString IntTy = "int"
      | toString BoolTy = "bool"
      | toString (TensorTy[]) = "real"
      | toString (TensorTy[2]) = "vec2"
      | toString (TensorTy[3]) = "vec3"
      | toString (TensorTy[4]) = "vec4"
      | toString (TensorTy dd) = concat["tensor[", String.concatWithMap "," Int.toString dd, "]"]
      | toString (TupleTy(tys)) = concat ["(", concat (List.map toString tys), ")"]
      | toString StringTy = "string"
      | toString (ImageTy(d, dd)) = concat[
         "image(", Int.toString d, ")[", String.concatWithMap "," Int.toString dd, "]"
        ]
      | toString (SeqTy(ty, NONE)) = toString ty ^ "[]"
      | toString (SeqTy(ty, SOME d)) = concat[toString ty, "[", Int.toString d, "]"]
      | toString (FemData data) = "FemData:" ^ (FemData.toString data)

    fun toOutputAbleTypes(ty) =
	(*Converts types to the form Tuple(base[]) so that base = base Ty ([n][n]...) ad so that base is in pre-order traversal order *)
	let
	 fun isTuple (n, TupleTy _) = true
	   | isTuple _ = false
	 (*inefficient but simple:
	  *First, flatten a nested tupple i.e Tuple(..., Tuple, ...)
	  *Second, traverse tuple to analyze next level
	  *If a seq type is ever found with a tuple directly inside it, switch them.
	  *If a finite seq type is ever found with an infinite one, make an equivelant tuple of infinite ones.
	  *Otherwise, continue down
	  *preverse other types, except fem types wher we fail
	  **By applying thes rules, we put the type in the suitable form and we maintain pre-order rules on base types as nodes (i.e if you wrote out the pre-order, took out all non-base types, the orders between the original type and toOutpuTableType would be the same.)
	  **Further, the [n] inside of [] nodes stay in the same order so we can itterate through them correctly.
	  *)
	 fun toat(ty) =
	     (case ty
	       of TupleTy(tys) => (case List.find isTuple (List.tabulate(List.length tys, fn x => (x, List.nth(tys, x))))
				    of SOME((n, TupleTy(tys'))) => let
				    in
				     TupleTy(List.take(tys, n) @ tys' @ List.drop(tys, n + 1))
				    end
				     | NONE => TupleTy(List.map toat tys)
				  (*end case*))
		| SeqTy(ty', r) => (case (ty', r) (*Fix me*)
				     of (TupleTy(tys), _) => TupleTy(List.map (fn t => SeqTy(t, r)) tys)
				       (*flip T[][5] to (T[], T[], T[], T[], T[])*)
				       |(SeqTy(ty'', NONE), SOME(k)) => TupleTy(List.tabulate(k, fn _ => SeqTy(ty'', NONE)))
				       | _ => SeqTy(toat ty', r)
				   (*end case*))
		| FemData _ => raise Fail "unexpected Fem Data"
		| _ => ty 
	     (*end case*))
	 fun loop ty = let
	  val ty' = toat ty
	 in
	  if (toString ty) = (toString ty')
          (*TODO: Fix stupid lazy hack -> will fail on FEM data,
	    	  but this will crash earlier when making access pattern *)
	  then ty'
	  else loop ty'
	 end
	in
	 (case loop ty
	   of TupleTy(tys) => tys
	    | t => [t]
	 (*end case*))
	end	  
		     

    fun isBase(IntTy) = true
      | isBase (BoolTy) = true
      | isBase (TensorTy _) = true
      | isBase (StringTy) = true
      | isBase (SeqTy(ty', SOME(n))) = isBase ty'
      | isBase _ = false
    fun baseSpec ty =
	let
	 fun bs (IntTy) = (1, IntTy)
	   | bs (BoolTy) = (1, BoolTy)
	   | bs (TensorTy(s)) = (List.foldr op* 1 s, TensorTy[1])
	   | bs (StringTy) = (1, StringTy)
	   | bs (SeqTy(ty', SOME(n))) = let val (size, ty'') = bs ty' in (n * size, ty'') end
	   | bs _ = raise Fail "nonbase ty has no base spec"
	 val ret = bs ty
	in
	 (ret)
	end

		 
    fun isOutputAble(ty) =
	(case ty
	  of IntTy => true
	   | BoolTy => true
	   | TensorTy _ => true
	   | StringTy => true
	   | ImageTy(_, _) => raise Fail "ImageTy should not be considered in output!"
	   | FemData _ => raise Fail "FemData should not be considered in tuple involved output"
	   | SeqTy(ty', SOME(n)) => isOutputAble ty' andalso not(hasDynamicSize ty') (*[][6] needs to be [6][]*)
	   | SeqTy(ty', NONE) => isOutputAble ty'
	   | TupleTy _ => false 
	(*end case*))
    (*is of form T[n][n]...[]*)
    fun isSingleOutputWithFem(ty) = (case ty
				      of SeqTy(ty', SOME(n)) => not(hasDynamicSize ty') andalso isSingleOutputWithFem ty'
				       | SeqTy(ty', NONE) => isSingleOutputWithFem ty'
				       | TupleTy _ => false
				       | ImageTy _ => raise Fail "ImageTy shoud not exist in output!"
				       | _ => true)
    (*for both inputs and outputs*)
    datatype path = TuplePath of int (*Tuple*)
		  | SeqArray of int option (*Seq or Array containing things other than a seq*)
		  | ArraySeq of int (*Array containing a seq*)
		  | BasePath of t * int * t (*copyable base modulo*)
    type iterate = (int * path) list (*itteration depth x access type*)
    datatype acc = TupleAcc of int | FixedArrayAcc of int | VarArrayAcc of int * int option | BaseCopy of t * int * t
    datatype loops = Fixed of int * int | From of int * acc list (*for x_int in (0..int-1) vs for x_int in something(acc list ...)*)
    (*for outputs:*)
    type copyOut = { loop : loops list,
		     accs : acc list,
		     outputTy : t,
		     dataItteration : iterate,
		     nrrdNum : int}

    (*for inputs:*)
    datatype inputType = BaseInput of acc list * int * (t * int * t)(*might want type here...*)
		       | NrrdSeqInput of acc list * int * acc list * (t * int * t)
    type inputJoin = (acc list * int list * inputType list)
    (* datatype inputNrrds = BaseInput of acc list * int * (t * int * t) | Join of inputJoin *)
    type copyIn = {
     diderotInputTy : t,
     inputTys : t list,
     inputs : inputType list,
     joins : (inputJoin * t) list
     

    }
    (*NrrdSeqInput - path how to go through the original type:
     The first part before the seq describes how index into the dest
     the seq part then tells us we want to load
     the rest of the part tells you where to set info from the load*)

    (*We create an iterate, which describes how to go through the original data*)
    (*Normal form: Tuple(Seq(Seq[n]...T))*)
    (*We can build an access pattern: Translate is fairly obvious.*)
    (*We can build the loop (var, range) where range is either fixed int or an access list*)
    (*We can use base to do the copy*)
    (*Picture:Base [n]Tuples[n]Tuple [] Tuples [n] Tuples  *)
    (*So the part outside the [], we have to make these fixed nrrds; the tuples parts inside also have to be listed. The Base[n] comes a copy!*)

    local
     fun convt(ty, itterDepth) =
	 if isBase ty
	 then let val (baseSize, valueTy) = baseSpec(ty)
	      in
	       [[(~1, BasePath(ty, baseSize, valueTy))]] : iterate list
	      end
	 else
	  (case (ty)
	    of SeqTy(ty', SOME(n)) => if hasDynamicSize ty'
				      then let
				       val extraMap =  (fn i =>  ArraySeq i)
				       val rest : iterate list = convt(ty', itterDepth + 1)
				       val conversions : (int * iterate list) list = (List.tabulate(n, fn i => (i, rest)))
				       val mapFunc = (fn (i,x) =>  List.map (fn j => (itterDepth, extraMap(i)) :: j) x)
				       val ret : iterate list = List.concatMap mapFunc conversions
				      in ret end
				      else let val next = convt(ty', itterDepth + 1)

		 			   in List.map (fn x => (itterDepth, SeqArray(SOME(n))) :: x) next end
	     | SeqTy(ty', NONE) => let val next = convt(ty', itterDepth + 1)
		 		   in List.map (fn x => (itterDepth, SeqArray(NONE)) :: x) next end
				     
	     | TupleTy(tys) => let val nexts : iterate list list = List.map (fn t => convt(t, itterDepth)) tys
				   val nexts' : (int * iterate list) list =
				       List.tabulate(List.length(tys), fn x => (x, List.nth(nexts, x)))
		 	       in List.concatMap (fn (i,x) => List.map (fn j => (~1, TuplePath(i)) :: j)  x) nexts' end
	     | _ => raise Fail "impossible ill formed output type"
	  (*end case*))
     fun buildAcceses(iter : iterate) =
	 let
	  fun bas (i, TuplePath(j)) = (TupleAcc(j))
	    | bas (i, ArraySeq(j)) = (FixedArrayAcc(j)) 
	    | bas (i, SeqArray(r)) = (VarArrayAcc(i, r))
	    | bas (_, BasePath(a,b,c)) = BaseCopy(a,b,c)
	 in
	  List.map bas iter
	 end

     (* fun normalizeAcc(accs : acc list) = *)
     (* 	 let *)
     (* 	  (*just bubble sort with a particular order? *)
     (* 	  feels like it. Tuples at the top (die) *)
     (* 	  SeqArray(NONE) *)
     (* 	  ArraySeq *)
     (* 	  Base *)
     (* 	   *) *)
     (* 	 in *)
     (* 	 end *)
		     

     fun buildLoop(sourceIter : iterate, sourceAccs : acc list) =
	 let
	  (*for each itter, drop tuple, collect acc*)
	  fun foldToLoop((iter, acc), (loop, accs)) =
	      (case (iter, acc)
		of ((_, TuplePath(_)), a) => (loop, a :: accs)
		 | ((_, ArraySeq(j)), a) => (loop, a :: accs)
		 | ((j, SeqArray(SOME(k))), a) => ((Fixed (j,k)) :: loop, a :: accs)
		 | ((j, SeqArray(NONE)), a) => ((From(j, (List.rev accs))) :: loop, a :: accs)
		 | ((_, BasePath _), _) => (loop, (List.rev accs))
	      (*end case*))
	  val (loo, lacc) = List.foldr foldToLoop ([],[]) (List.rev(ListPair.zip(sourceIter, sourceAccs)))
				       
	 in
	  (List.rev loo, sourceAccs, sourceIter)
	 end
     fun scanPath (_, TuplePath(_)) = false
       | scanPath (_, SeqArray(NONE)) = true
       | scanPath (_, SeqArray(SOME(_))) = false
       | scanPath (_, ArraySeq(_)) = false
       | scanPath (_, BasePath(_)) = false

     fun categorizeInput(idx : int, iterate : iterate) =
	 (case (List.find (fn (a,b) => scanPath b) (List.tabulate(List.length iterate, fn x => (x, List.nth(iterate, x)))))
	   of NONE => let val BaseCopy(r) = (List.last(buildAcceses iterate)) in BaseInput(buildAcceses iterate, idx, r) end
	    | SOME((k, (i, p))) =>
	      let
	       val toSeq = buildAcceses (List.take(iterate, k))
	       val seqSave = List.rev(List.drop(List.rev(buildAcceses (List.drop(iterate, k+1))), 1))
	       val BaseCopy(r) = (List.last(buildAcceses iterate))
	      in
	       NrrdSeqInput(toSeq, idx, seqSave, r)
	      end
	 (*end case*))
     fun accPrefixEq(TupleAcc(j), TupleAcc(k)) = j=k
       | accPrefixEq (FixedArrayAcc(j), FixedArrayAcc(k)) = j=k
       | accPrefixEq (VarArrayAcc(a,SOME(b)), VarArrayAcc(c,SOME(d))) = a=c andalso b=d
       | accPrefixEq (VarArrayAcc(a,NONE), _) = raise Fail "accPrefix can't contain []"
       | accPrefixEq (_, VarArrayAcc(a,NONE)) = raise Fail "accPrefix can't contain []"
       | accPrefixEq (BaseCopy _, _) = raise Fail "accPrefix can't contain base"
       | accPrefixEq (_, BaseCopy _) = raise Fail "accPrefix can't contain base"
       | accPrefixEq _ = false

     fun equalLoad (BaseInput(_), _) = false
       | equalLoad (_, BaseInput(_)) = false
       | equalLoad (NrrdSeqInput(a, _, _ , _), NrrdSeqInput(b, _, _, _)) =
	 (List.length(a) = List.length(b))
	 andalso (ListPair.all (accPrefixEq) (a,b))
     fun buildJoins(ins : inputType list, tyList) =
	 let
	  (*with current list, *)
	  fun joinRecurse(BaseInput(_)::ilist, joins, t::tl) = joinRecurse(ilist, joins, tl)
	    | joinRecurse ((i as NrrdSeqInput(prefix, id, _, _))::ilist, joins, t::tl) =
	      let
	       val detect = fn x => equalLoad(i,x)
	       val same = List.filter detect ilist
	       val diff = List.filter (not o detect) ilist
	       val newJoin = ((prefix, id :: (List.map (fn NrrdSeqInput(_, id', _, _ ) => id') same), i::same), t)

	      in
	       joinRecurse(diff, newJoin :: joins, tl)
	      end
	    | joinRecurse ([], joins, []) = List.rev joins
	 in
	  joinRecurse(ins, [], tyList)
	 end
		   


	   
    in
    fun buildOuputputConversionRoutine(ty) =
	let
	 val typeItteration : iterate list = convt(ty, 0)
	 val accItteration = List.map buildAcceses typeItteration
	 val loopAndAccAndItter = ListPair.map buildLoop (typeItteration, accItteration)
	 val outputTypes = ListPair.zip (toOutputAbleTypes ty, List.tabulate(List.length typeItteration, fn x => x))
	 val conversions = ListPair.map (fn (((ot, n), (lo, sa, it))) =>
					    {loop = lo, accs = sa,
					     outputTy = ot, dataItteration=it,
					     nrrdNum = n})
					(outputTypes, loopAndAccAndItter)
	in
	 conversions : copyOut list(*loop to iterate over the API form instance, acces on the variable, how to itterate over the source*)
	end
    fun buildInputConversionRoutine(ty :t) : copyIn =
	let
	 val typeItteration : iterate list = convt(ty, 0)
	 val tabs = List.tabulate(List.length typeItteration, fn x => x)
	 val inputs = ListPair.map categorizeInput (tabs, typeItteration)
	 val inputTys = toOutputAbleTypes ty
	 val joins = buildJoins(inputs, inputTys)

	 val _ = if List.length inputTys <> List.length typeItteration
		 then raise Fail "error in conversion for inputs"
		 else ()
	in
	 {diderotInputTy = ty, inputTys=inputTys, inputs = inputs, joins = joins}
	end
    end

  end
