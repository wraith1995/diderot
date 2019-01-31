(* array-nd.sml 
 *
 * This is a numpy N-d array. Think of it as a generalization of Array2.sml.
 * The compiler already uses the ideas in here in several places so I'm generalizing it. 
 *)


structure ArrayNd : sig
	   datatype 'a ArrayNd = AND of {
		     dims : int,
		     shape : int list,
		     elems : 'a Array.array,
		     index : (int list -> int option) option ref,
		     inverseIndex : ((int list) Array.array) option ref
		    }


	   (* Key indexing features*)
	   val computeIndex : 'a ArrayNd -> 'a ArrayNd
	   val getIndex : 'a ArrayNd -> (int list -> int option)
	   val computeInverseIndex : 'a ArrayNd -> 'a ArrayNd
	   val getInverseIndex : 'a ArrayNd -> (int list) Array.array
	   val buildIndexInfo : int list -> (int list -> int option) * (int list) Array.array

	   (*a constant*)
	   val maxLen : int
	   (* ways to build these guys*)
	   val fromList : 'a list -> 'a ArrayNd
	   val fromList' : 'a list * int list -> 'a ArrayNd
	   val fromArray : 'a Array.array -> 'a ArrayNd
	   val fromArray' : 'a Array.array * int list -> 'a ArrayNd
	   val duplicate : 'a ArrayNd -> 'a ArrayNd
	   val zip : 'a ArrayNd * 'b ArrayNd -> ('a * 'b) ArrayNd
	   val concat : 'a ArrayNd list -> 'a ArrayNd
	   val array' : 'a * int list -> 'a ArrayNd
	   val subregion : 'a ArrayNd * int list * int list -> 'a ArrayNd
	   (* get meta data*)
	   val length : 'a ArrayNd -> int
	   val compatibleMeta : 'a ArrayNd * 'b ArrayNd -> bool
	   val indexInside : 'a ArrayNd * int list -> bool
	   val sameData : 'a ArrayNd * 'a ArrayNd -> bool
	   val same : 'a ArrayNd * 'a ArrayNd -> bool
	   val shape : 'a ArrayNd -> int list
	   (* the usual*)
	   val sub : ('a ArrayNd * int) -> 'a
	   val sub' : ('a ArrayNd * int list) -> 'a
	   val app : ('a -> unit) -> 'a ArrayNd -> unit
	   val appi : (int * 'a -> unit) -> 'a ArrayNd -> unit
	   val appi' : (int list * 'a -> unit) -> 'a ArrayNd -> unit
	   val foldr : ('a * 'b -> 'b) -> 'b -> 'a ArrayNd -> 'b
	   val foldri : (int * 'a * 'b -> 'b) -> 'b -> 'a ArrayNd -> 'b
	   val foldri' : (int list * 'a * 'b -> 'b) -> 'b -> 'a ArrayNd -> 'b
	   val modify : ('a -> 'a) -> 'a ArrayNd -> unit
	   val modifyi : (int * 'a -> 'a) -> 'a ArrayNd -> unit
	   val modifyi' : (int list * 'a -> 'a) -> 'a ArrayNd -> unit

	   val expandMap : ('a -> 'a Array.array) -> int -> 'a ArrayNd -> 'a ArrayNd
	   val reshape: 'a ArrayNd * int list -> 'a ArrayNd

								  
	   val all : ('a -> bool) -> 'a ArrayNd -> bool
	   (*discouraged misc*)
	   val toList : 'a ArrayNd -> 'a list
	   val convertToTree : 'b ArrayNd * ('b -> 'a) * ('a list * int -> 'a) -> 'a

	  end = struct
datatype 'a ArrayNd = AND of {
	  dims : int,
	  shape : int list,
	  elems : 'a Array.array,
	  index : ((int list -> int option) option) ref,
	  inverseIndex : (((int list) Array.array) option) ref
	 }
(* datatype 'a TreeRep = Val of 'a 1| Tree of 'a TreeRep list *)
(*use "a.sml"; structure A = ArrayNd; val aaaa = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; val z = A.fromList'(aaaa, [2,3,4]); val i = A.getIndex z; val ii = A.getInverseIndex z; val i' = Option.valOf o i; fun ii' j = Array.sub(ii,j);
use "a.sml"; structure A = ArrayNd; val aaaa = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; val z = A.fromList'(aaaa, [2,3,4]); val i = A.getIndex z; val ii = A.getInverseIndex z; val i' = Option.valOf o i; fun ii' j = Array.sub(ii,j); val test = i' o ii'; val test' = ii' o i';*)
fun getShape(AND{shape,...}) = shape
fun getDims(AND{dims,...}) = dims
fun getArray(AND{elems,...}) = elems
fun getIndex(AND{index,...}) = index
fun getInverseIndex(AND{inverseIndex,...}) = inverseIndex			       

fun valid(AND({dims, shape, ...})) = let val len = (List.length shape) in len = dims andalso len > 0 andalso (List.all (fn x => x >= 0) shape) end
fun computeIndex(t as AND{dims, shape, elems, index, inverseIndex}) = (
 case !index
  of SOME(_) => t
   | NONE => 
    let
     fun product xs = List.foldr (fn (x,y) => x*y) 1 xs
     fun sum xs = List.foldr (fn (x,y) => x+y) 0 xs
     val subArraySizes = List.foldr (fn (x,y) => (case y
						   of y'::ys => (y'*x)::y)) [1] ((shape))
     val (subArraySizes', size) = (case subArraySizes
				    of [] => raise Fail "invalid array"
				     | x::xs => (xs, x))
     fun index'(idx) = SOME(ListPair.foldr (fn (x,y,z) => z + x*y) 0 (subArraySizes', idx)) handle exn => NONE
    in
     AND{dims=dims, shape=shape, elems=elems, index=ref (SOME(index')), inverseIndex=inverseIndex}
    end
(* end case*))


fun getIndex(t) = (case computeIndex t
		    of AND{index,...} => Option.valOf (!index))


fun computeInverseIndex(t as AND{dims, shape, elems, index, inverseIndex}) = (
 case !inverseIndex
  of SOME(_) => t
   | NONE => 
     let
      fun product xs = List.foldr (fn (x,y) => x*y) 1 xs
      val subArraySizes = List.foldr (fn (x,y) => (case y
						    of y'::ys => (y'*x)::y)) [1] ((shape))
      val (subArraySizes', size) = (case subArraySizes
				     of [] => raise Fail "invalid array"
				      | x::xs => (xs, x))

      fun tabulator(num, [1], collected) = List.rev (num :: collected)
	| tabulator(num, modder::modding, collected) = let
	 val (mod', div') = (Int.mod(num, modder), Int.div(num, modder))
	in
	 tabulator(mod', modding, div' :: collected)
	end
      fun tab(num) = tabulator(num, subArraySizes', [])
      val inverseIndex' = Array.tabulate(size, tab)
     in
      AND{dims=dims, shape=shape, elems=elems, index=index, inverseIndex=ref (SOME(inverseIndex'))}
     end
(* end case*))

fun getInverseIndex t = (case computeInverseIndex t
			  of AND{inverseIndex,...} => Option.valOf (!inverseIndex)
			)

fun buildIndexInfo shape =
    let
     val temp = AND{dims = List.length shape, shape = shape, elems = Array.fromList [], index = ref NONE, inverseIndex = ref NONE}
		   
    in
     (getIndex temp, getInverseIndex temp)
    end

val maxLen = Array.maxLen
fun fromList(list) =
    let
     val len = List.length list
    in
     AND{
      dims = 1,
      shape = [len],
      elems = Array.fromList list,
      index = ref NONE, inverseIndex= ref NONE
     }
    end
fun fromArray(array) =
    let
     val len = Array.length array
    in
     AND{
      dims = 1,
      shape = [len],
      elems = array,
      index = ref NONE, inverseIndex= ref NONE
     }
    end
fun fromList'(list, shape) = AND{dims = List.length shape, shape=shape, elems = Array.fromList list, index = ref NONE, inverseIndex = ref NONE}
fun fromArray'(array, shape) = AND{dims = List.length shape, shape=shape, elems = array, index = ref NONE, inverseIndex = ref NONE}
fun length(AND{shape,...})  = List.foldr op* 1 shape
fun update'(t as AND{elems,...}, idx, elem) =
    let
     val index = getIndex t
     val intIndex = index idx		  
    in
     (case intIndex
       of NONE => ()
	| SOME(i) => Array.update(elems, i, elem))
    end
fun sub (t as AND{elems, ...}, idx) = Array.sub(elems, idx)
fun sub'(t as AND{elems,...}, idx) =
    let
     val index = getIndex t
     val intIndex = index idx
    in
     (case intIndex
       of NONE => raise Fail "impossible"
	| SOME(i) => Array.sub(elems,i))
    end
fun arrayToList arr = Array.foldr (op ::) [] arr
fun toList (AND{elems,...}) = arrayToList elems

(* get indexing scheme, mapSame, hash?*)

fun duplicate (t as AND{dims, shape, elems, index, inverseIndex}) =
    AND{dims = dims, shape = shape, index = ref (!index), inverseIndex = ref (!inverseIndex),
	elems =Array.tabulate(length t , fn idx => Array.sub(elems,idx))}




fun compatibleMeta((AND({dims=dims1, shape=shape1, ...})),
		   (AND({dims=dims2, shape=shape2, ...}))) =
    dims1=dims2 andalso (ListPair.all op= (shape1,shape2))
fun sameData (AND{elems=elems1,...}, AND{elems=elems2, ...}) = elems1 = elems2
fun same(t1,t2) = compatibleMeta(t1, t2) andalso sameData(t1, t2)
fun shape(AND{shape,...}) = shape
fun indexInside(AND{shape,...}, idx) = List.all (fn x => x > 0) (ListPair.map op- (shape, idx))
				       handle exn => false
fun zip(t1 as AND{dims,shape, elems=elems1, index=index1, inverseIndex=inverseIndex1},
	t2 as AND{elems=elems2,index=index2, inverseIndex=inverseIndex2,...}) =
    let
     val _ = if compatibleMeta(t1,t2)
	     then ()
	     else raise Fail "arrays to zip not of same type"
     fun tabFunc idx = (Array.sub(elems1,idx), Array.sub(elems2, idx))
     fun copy a = (ref (!a))
     val newArray = Array.tabulate(length t1, tabFunc)
     val index = (case (!index1, !index2)
		   of (SOME(a), _) => copy index1
		    | (_, SOME(a)) => copy index2
		    | _ => ref NONE)
     val inverseIndex = (case (!inverseIndex1,!inverseIndex2 )
			  of (SOME(a),_) => copy inverseIndex1
			   | (_, SOME(a)) => copy inverseIndex2
			   | _ => ref NONE)
    in
     AND{dims=dims, shape=shape, elems=newArray, inverseIndex=inverseIndex, index=index}
    end
fun applyToElems f (AND{dims, shape, elems, index, inverseIndex}) = f elems 
fun all f (AND({elems,...})) = Array.all f elems

fun app f (AND{elems,...}) = Array.app f elems
fun appg f (AND{elems,...}) = Array.app f elems
fun appi f (AND{elems, ...}) = Array.appi f elems
(* function takes index list*)
fun appi' f (t as AND{elems,...}) = let val inverseIndex = getInverseIndex t
				    in Array.appi (fn (x,y) => f(Array.sub(inverseIndex, x), y)) elems end
fun modify f t = applyToElems (Array.modify f) t
fun modifyi f t = applyToElems (Array.modifyi f) t
fun modifyi' f (t as AND{elems,...}) = let val inverseIndex = getInverseIndex t
		   in Array.modifyi (fn (x,y) => f(Array.sub(inverseIndex, x), y)) elems end

fun foldr f s (AND{elems, ...}) = Array.foldr f s elems
fun foldri f s (AND{elems, ...}) = Array.foldri f s elems
fun foldri' f s (t as AND{elems, ...}) = let val inverseIndex = getInverseIndex t
					 in Array.foldri (fn (x,y,z) => f(Array.sub(inverseIndex, x), y, z)) s elems end
				      
(* fun subArray idx (AND{dims, shape, elems, index, inverseIndex}) = *)
(*     let *)
(*      val idxLen = List.length idx *)
(*      val subIdx = idxLen < dims *)
(*      val subDim =  *)
(*     in *)
(*     end *)

(* ListOfArrays*)

fun array'(a, shape) =
    let
     val size = List.foldr op* 1 shape
     val new = Array.array(size, a)
    in
     AND{shape=shape, dims = List.length shape, elems = new, index = ref NONE, inverseIndex = ref NONE}
    end
fun subregion(t as AND{elems, shape, ...}, offset, limit) =
    let
     val newShape = ListPair.map (op-) (limit, offset)
     val first = Array.sub(elems, 0)
     val new = array'(first, newShape)
     fun set(idx, a) = sub'(t, ListPair.map (op+)(offset, idx))
    in
     (modifyi' set new; new)
    end


fun concat(arrays) =
    let
     val len = List.length arrays
     val _ = if len <= 1 then raise Fail "Can't concat one array!" else ()
     val dimss = List.map getDims arrays
     val startDim = List.nth(dimss,0)
     val _ = if (List.all (fn x => startDim = x ) dimss) then () else
	     raise Fail "can't concat arrays of different sizes"
     val shapes = List.map getShape arrays
     val startShape = List.nth(shapes, 0)
     fun sameShape x = ListPair.allEq op= (x, startShape)
     val _ = if (List.all sameShape shapes) then ()
	     else raise Fail "can't concat arrays of different shapes"
     val size = List.foldr op* 1 startShape
     val newSize = size * len
     val start = Array.sub(getArray(List.nth(arrays, 0)),0)
     val newArray = Array.array(newSize, start)
     val units = List.tabulate(len, fn x =>
				       let
					val a = getArray(List.nth(arrays, x))
				       in
					Array.copy({src=a,dst=newArray,di=x*size})
				       end
			      )

     (* val index' = ref (List.find (Option.isSome) (List.map (fn x => !getIndex(x)) arrays)) *)
     (* val inverseIndex' = ref (List.find (Option.isSome) (List.map (fn x => !getInverseIndex(x)) arrays)) *)
     val result = AND{dims=startDim+1, shape=len::startShape, elems=newArray, index = ref NONE, inverseIndex = ref NONE}
    in
     (computeIndex o computeInverseIndex) result
    end
(*do a list version, then reshape*)
fun expandMap f size (t as AND{dims, shape, elems, ...}) =
    let
     val newShape = List.@(shape, [size])
     val newSize = (length t) * size
     val newDims = dims + 1
     val newArray = Array.array(newSize, Array.sub(elems, 0))
     val newT = AND{dims = newDims, shape = newShape, elems = newArray, index = ref NONE, inverseIndex = ref NONE}
			    
     (* val indexArray = getInverseIndex t *)
     (* val indexFunc = (fn i => Array.sub(indexArray, i)) *)
     val inverseIndex = Option.valOf o (getIndex newT)


     fun appMap(idx, a) =
	 let
	  val newIndex = List.@(idx, [0])
	  val offset = inverseIndex newIndex
	  val newArray = f a
	 in
	  Array.copy({src=elems, dst = newArray, di = offset})
	 end
			    
	       (* expand the array *)
	       (* copy each one over...*)
	       (* send result*)
    in
     (appi' appMap t; newT)
    end
fun reshape(AND{dims, shape, elems, ...}, shape') = (computeIndex o computeInverseIndex) (AND{dims=List.length shape', shape=shape', elems = elems, index = ref NONE, inverseIndex = ref NONE})
fun toArrayMap expandMap insertShape (t as AND{dims, shape, elems, ...}) =
    let
     val a = ()
    in
     ()
    end
(* n array grab at a time*)
fun convertToTree(t as AND{dims, shape, elems, ...}, convert, convertFunc) =
    let
     val listRep = List.map convert (toList t)
     val parseOrder = List.rev shape
     fun groupLists(n, [], lists) = List.rev lists
       | groupLists (n, list, lists) =
	 groupLists(n, List.drop(list, n), convertFunc(List.take(list,n), n)::lists)
	 handle exn => raise exn

     fun doDims([d], list) =
	 let
	  val int = groupLists(d, list, [])
	 in
	  (case int
	    of [x] => x
	     | _ => raise Fail "impossible")
	 end
       | doDims(d::ds, list ) =
	 let
	  val int = groupLists(d, list, [])
	 in
	  doDims(ds, int)
	 end

    in
     doDims(parseOrder, listRep)
    end

					  
end
