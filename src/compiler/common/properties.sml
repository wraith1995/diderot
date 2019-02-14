(* properties.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *
 * Various properties of Diderot programs that we want to track.
 *)

structure Properties =
  struct

  (* program properties *)
    datatype t
      = HasConsts               (* present if the program has const variables *)
      | HasInputs               (* present if the program has input variables *)
      | HasGlobals              (* present if the program has global variables *)
      | StrandArray             (* present if the strands are organized in a array *)
      | DynamicSeq              (* present if strands compute dynamic sequences *)
      | StrandsMayDie           (* present if strands may die *)
      | HasStartMethod          (* present if there is a user-defined start method *)
      | HasStabilizeMethod      (* present if there is a user-defined stabilize method *)
      | NewStrands              (* present if new strands may be created dynamically *)
      | StrandCommunication     (* present if strands read the state of other strands *)
      | GlobalInit              (* present if there is a global initialization block *)
      | GlobalStart             (* present if there is a global start block *)
      | GlobalUpdate            (* present if there is a global update block *)
      | GlobalReduce            (* present if there is a global reduction (implies GlobalUpdate) *)
      | KillAll                 (* present if there is a global die statement *)
      | StabilizeAll            (* present if there is a global stabilize statement *)
      | StrandSets              (* present if there are strand sets in the global update *)
      | NeedsBSP                (* present if there is a feature (other than those listed above)
                                 * that requires BSP execution.
                                 *)
      | NeedsExtraLibs of string list * string list * string list * string list

    fun toString HasConsts = "HasConsts"
      | toString HasInputs = "HasInputs"
      | toString HasGlobals = "HasGlobals"
      | toString StrandArray = "StrandArray"
      | toString DynamicSeq = "DynamicSeq"
      | toString StrandsMayDie = "StrandsMayDie"
      | toString HasStartMethod = "HasStartMethod"
      | toString HasStabilizeMethod = "HasStabilizeMethod"
      | toString NewStrands = "NewStrands"
      | toString StrandCommunication = "StrandCommunication"
      | toString GlobalInit = "GlobalInit"
      | toString GlobalStart = "GlobalStart"
      | toString GlobalUpdate = "GlobalUpdate"
      | toString GlobalReduce = "GlobalReduce"
      | toString StabilizeAll = "StabilizeAll"
      | toString KillAll = "KillAll"
      | toString StrandSets = "StrandSets"
      | toString NeedsBSP = "NeedsBSP"
      | toString (NeedsExtraLibs(includes, includeDirs, links, linkDirs)) =
	"ExtraFlags(-#[" ^ (String.concatWith "," includes) ^"],-I["
	^ (String.concatWith "," includeDirs) ^ "],-l["
	^ (String.concatWith "," links) ^ "],-L["
	^ (String.concatWith "," linkDirs) ^ "])"


    fun propsToString [] = "none"
      | propsToString [p] = toString p
      | propsToString props = String.concatWithMap "," toString props

    fun hasProp (prop : t) = List.exists (fn p => prop = p)

    fun getProps (prop :t) = List.filter (fn p => prop = p)

    fun clearProp (prop : t) = List.filter (fn p => (p <> prop))

    fun collectUserCompileInserts(props : t list) =
	let
	 fun f (NeedsExtraLibs(a1,b1,c1,d1), (a,b,c,d)) = (List.@(a1,a),
							   List.@(b1,b),
							   List.@(c1,c),
							   List.@(d1,d))
	   | f (_, s) = s
	 val (r1, r2, r3, r4) = List.foldr f ([],[],[],[]) props
	 fun nub([]) = []
	   | nub (x::xs) = if List.exists (fn y => y=x) xs
			   then nub(xs)
			   else x::nub(xs)
	in
	 (nub r1, nub r2, nub r3, nub r4)
	end

  (* NOTE: this function is less precise than the TargetSpec.dualState function *)
    val dualState = hasProp StrandCommunication

  end
