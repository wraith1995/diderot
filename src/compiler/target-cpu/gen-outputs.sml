(* gen-outputs.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *
 * Generate strand output functions.  The output formats always have a single axis for the
 * data elements followed by one, or more, axes for the output structure.  There are four
 * cases that we handle:
 *
 *      grid, fixed-size elements:
 *              nrrd has object axis followed by grid axes
 *
 *      collection, fixed-size elements
 *              nrrd has object axis followed by a single axis
 *
 *      grid, dynamic-size elements
 *              nLengths nrrd has size 2 for objects (offset, length) followed by grid axes
 *              nData nrrd has object axis followed by a single axis
 *
 *      collection, dynamic-size elements
 *              nLengths nrrd has size 2 for objects (offset, length) followed by a single axis
 *              nData nrrd has object axis followed by a single axis
 *
 * The object axis kind depends on the output type, but it will either be one of the tensor types
 * that Teem knows about or else nrrdKindList.  In any case, the data elements are written as a
 * flat vector following the in-memory layout.  The other axes in the file will have nrrdKindSpace
 * as their kind.
 *
 * TODO: some of this code will be common across all targets (e.g., writing outputs to files), so
 * we will want to refactor it.
 *
 * TODO: for sequences of tensors (e.g., tensor[3][2]), we should use a separate axis for the
 * sequence dimension with kind nrrdKindList.
 *
 * TODO: since the runtime tracks numbers of strands in various states, we should be
 * able to use that information directly from the world without having to recompute it!
 *)

structure GenOutputs : sig

  (* gen (props, nAxes, outputs)
   *    returns a list of function declarations for getting the output/snapshot nrrds from
   *    the program state.  The arguments are:
   *        props       - the target information
   *        nAxes       - the number of axes in the grid of strands (NONE for a collection)
   *        outputs     - the list of output state variables paired with their API types
   *)
    val gen : CodeGenEnv.t * int option * OutputUtil.output_info list -> CLang.decl list

  end = struct

    structure IR = TreeIR
    structure V = TreeVar
    structure Ty = APITypes
    structure CL = CLang
    structure Nrrd = NrrdEnums
    structure U = GenOutputsUtil
    structure RN = CxxNames
    structure Env = CodeGenEnv
    structure GenAPI = GenLibraryInterface

    fun mapi f l = let
          fun mapf (i, [], l) = List.rev l
            | mapf (i, x::xs, l) = mapf (i+1, xs, f(i, x)::l)
          in
            mapf (0, l, [])
          end

    val nrrdPtrTy = CL.T_Ptr(CL.T_Named "Nrrd")
    val sizeTy = CL.T_Named "size_t"
    val wrldPtr = RN.worldPtrTy
    fun mkInt i = CL.mkInt(IntInf.fromInt i)

  (* variables in the generated code *)
    val wrldV = CL.mkVar "wrld"
    val sizesV = fn x => CL.mkVar ("sizes_" ^ Int.toString(x))
    val iV = CL.mkVar "ix"
    val nV = CL.mkVar "n"
    val cpV = CL.mkVar "cp"
    val ipV = CL.mkVar "ip"
    val msgV = CL.mkVar "msg"
    val offsetV = fn i => CL.mkVar ("offset_" ^ (Int.toString i))
    val nDataV = fn i => CL.mkVar ("nData_"  ^ (Int.toString(i)))
    val nLengthsV = fn i => CL.mkVar ("nLengths_" ^ (Int.toString i))
    val numElemsV = fn i => CL.mkVar ("numElems_" ^ (Int.toString i))
    val DIDEROT_DEAD = CL.mkVar "diderot::kDead"
    val DIDEROT_STABLE = CL.mkVar "diderot::kStable"
    val varName = (fn (CL.E_Var(x)) => x)

    fun strandMeth (f, args) = CL.mkDispatch(CL.mkIndirect(wrldV, "_strands"), f, args)

  (* dymanic sequence operations *)
    fun seqLength seq = CL.mkDispatch(seq, "length", [])
    fun seqCopy (seq, dst) = CL.mkDispatch(seq, "copy_to", [dst])

    (*analyze a type to determine the nature of the loop needed to copy it over*)
    fun tyAnalysis(ty : APITypes.t, flag) =
	let
	 fun tyAnalysis'(APITypes.SeqTy(ty', a), xs, flag) = tyAnalysis'(ty', a::xs, flag)
	   | tyAnalysis'(APITypes.FemData(_), xs, flag) = if flag
							  then (true, NONE :: (List.rev xs))
							  else (true, List.rev xs)
	   | tyAnalysis'(_, xs, _) = (false, List.rev xs)
	in
	 tyAnalysis'(ty, [], flag)
	end
    fun copy(dynSeq, ty, copySource, nElems, elemCTy, SOME(loopinfo), dynSeqFunc) = 
	let
	 val (fem, seqDims) = tyAnalysis( ty, dynSeq)  (*TODO: Clean this crap up*)
	 fun baseTy(accs) = (case List.last accs
			      of Ty.BaseCopy(ty, num, elemTy) => (ty, num, elemTy)
			       | _ => raise Fail "missing base ty")
	 val (_, nElems, _) = baseTy (#accs loopinfo) (*replace nElemns with base info to account for outTy over esimtate whereas outTy is corresct correct for estimating sizes.*)
	 val _ = if fem
		 then raise Fail "Fem and complex output types not compatible currently"
		 else ()
	 fun buildAcc(a::accs : APITypes.acc list, base) =
	     (case a
	       of APITypes.TupleAcc(j) => buildAcc(accs, CL.mkSelect(base, "t_" ^ Int.toString j))
		| APITypes.FixedArrayAcc(j) => buildAcc(accs, CL.E_Subscript(base, CL.E_Int(IntLit.fromInt j, CL.intTy)))
		| APITypes.VarArrayAcc(j, SOME _) => buildAcc(accs, CL.E_Subscript(base, CL.mkVar ("idx_" ^ Int.toString j)))
		| APITypes.VarArrayAcc(j, NONE) => (case dynSeqFunc
						     of NONE => buildAcc(accs, CL.E_Subscript(base, CL.mkVar ("idx_" ^ Int.toString j)))
						      | SOME(f) => base)
		| APITypes.BaseCopy _ => base
	     (*end case*))
	   | buildAcc([], base) = base
	 fun loopMap(copySource, loop : APITypes.loops) =
	     let
	      val loopVar = fn x => CL.mkVar ("idx_" ^ Int.toString x)
	      val getStr = fn (CL.E_Var(s)) => s
	      in
	     (case loop
	       of APITypes.Fixed(i, j) => (fn x => CL.S_For(CL.T_Named("auto"),
						   [((getStr o loopVar) i, mkInt 0)],
						   CL.E_BinOp(loopVar i, CL.#<, mkInt j),
						   [CL.E_PostOp(loopVar i, CL.^++)],
						   CL.S_Block([x])))
		| APITypes.From(i, accs) =>
		  let val acc = buildAcc(accs, copySource)
		  in (fn x => CL.S_For(CL.T_Named("auto"),
				       [((getStr o loopVar) i, mkInt 0)],
				       CL.E_BinOp(loopVar i, CL.#<, CL.mkDispatch(acc, "length", [])),
				       [CL.E_PostOp(loopVar i, CL.^++)],
				       CL.S_Block([x])))
		  end
		     
	     (*end case*))
	     end
	 fun buildLoopFunction(loops : APITypes.loops list, copySource) =
	     let
	      val loopFuncs = List.map (fn x => loopMap(copySource, x)) loops
	      val loopFuncs = (case dynSeqFunc
				of NONE => loopFuncs
				 | SOME _ =>
				   let
				    val tabed = List.tabulate(List.length loops, fn x => (x, List.nth(loops, x)))
				    val (id, _) = Option.valOf (List.find (fn (j, Ty.From _) => true | _ => false) tabed)
				   in
				    List.take( loopFuncs, id) 
				   end
			      (*end case*))
	      fun clId(x : CL.stm ) : CL.stm = x
	      val loopFunc = List.foldr (fn (a,b) => a o b) clId loopFuncs
	     in
	      loopFunc
	     end
	 fun buildLoop(loopInfo : APITypes.copyOut, copyTarget, copySource) =
	     let
	      val accDesc = buildAcc(#accs loopInfo, copySource)
	      val loopFunc = buildLoopFunction(#loop loopInfo, copySource)
	      val transferStm = if nElems <> 1 (*Warning: might not be the best test...*)
				then
				 CL.mkCall("memcpy",
					    [
					      copyTarget,
					      CL.mkUnOp(CL.%&, CL.E_Grp(accDesc)),
					      CL.mkBinOp(mkInt nElems, CL.#*, CL.mkSizeof elemCTy)
					  ])
				else
				 CL.S_Exp(
				  CL.E_AssignOp(CL.E_UnOp(CL.%*, CL.E_Grp(CL.E_Cast(CL.T_Ptr(elemCTy), copyTarget))),
						CL.$=, accDesc)
				 )
	      val copyStm = (case List.last (#accs loopInfo)
			      of APITypes.BaseCopy(baseTy, baseSize, valueTy) =>
				 CL.S_Block[
				  transferStm,
				  CL.mkExpStm(CL.mkAssignOp(copyTarget,
							    CL.+=,
							       CL.mkBinOp(mkInt nElems, CL.#*, CL.mkSizeof elemCTy)))
				 ]
			      | _ => raise Fail "invalid final accesor! Base should be present"
			    (*end class*))
	     in
	      (case dynSeqFunc
		of NONE => loopFunc copyStm
		 | SOME f => loopFunc (f accDesc)
	      (*end case*))
	     end
	in
	 buildLoop(loopinfo, cpV, copySource)
        end
		(*Copy for simple output types or copy for simple fem types*)
      | copy(dynSeq, ty, copyTarget, nElems, elemCTy, NONE, dynSeqFunc) = 
	let
	 val (fem, seqDims) = tyAnalysis( ty, dynSeq)
	 fun mkLoop([], vars, copyTarget) =
	     let
	      val vars' = List.rev vars
	      val acc = List.foldr (fn (x,y) => CL.E_Subscript(y, CL.mkVar x)) copyTarget vars'
	     in
	      
	     CL.mkBlock([
				     CL.mkExpStm(CL.mkApply("copy_to", [acc, cpV])),
				     CL.mkExpStm(CL.mkAssignOp(cpV,
							       CL.+=,
								  CL.mkBinOp(mkInt nElems, CL.#*, CL.mkSizeof elemCTy)))
		       ])
	     end
	   | mkLoop(SOME(x)::xs, vars, copyTarget) =
	     let
	      val n = List.length (SOME(x)::xs)
	      val var = "ix"^(Int.toString n)
	      val rest = mkLoop(xs, var::vars, copyTarget)
	     in
	      CL.mkFor(CL.intTy, [(var, CL.E_Int(IntLit.fromInt 0, CL.intTy))],
		       CL.E_BinOp(CL.E_Var(var), CL.#<, CL.E_Int(IntLit.fromInt x, CL.intTy)),
		       [CL.E_PostOp(CL.E_Var(var), CL.^++)],
		       rest
		      )
	     end
	   | mkLoop(NONE::xs, vars, copyTarget) =
	     let
	      val n = List.length (NONE::xs)
	      val var = "ix"^(Int.toString n)
	      val copyTarget' = "seqVar"^(Int.toString n)
	      val rest = mkLoop(xs, [], CL.E_Subscript(copyTarget, CL.mkVar var))
	     in
	      CL.mkFor(CL.intTy, [(var, CL.E_Int(IntLit.fromInt 0, CL.intTy))],
		       CL.E_BinOp(CL.E_Var(var), CL.#<, CL.mkDispatch(copyTarget, "length", [])),
		       [CL.E_PostOp(CL.E_Var(var), CL.^++)],
		       CL.S_Block(
		       [rest]
		      ))
	     end
	 (*need a case for NONE*)
				    
	 fun copy(true, false, NONE) = CL.mkBlock[
	      CL.mkAssign(cpV, seqCopy(copyTarget, cpV))]
	   | copy(false, false, NONE) =  CL.mkBlock[
              CL.mkCall("memcpy", [
                        cpV,
                        CL.mkUnOp(CL.%&, CL.E_Grp(copyTarget)),
                        CL.mkBinOp(mkInt nElems, CL.#*, CL.mkSizeof elemCTy)
                       ]),
              CL.mkExpStm(CL.mkAssignOp(cpV,
					CL.+=,
					   CL.mkBinOp(mkInt nElems, CL.#*, CL.mkSizeof elemCTy)))
             ]
	   | copy (false, true, NONE) = mkLoop(seqDims, [], copyTarget)
	   | copy (true, true, NONE) = mkLoop(seqDims, [], copyTarget)
	   | copy (_, _, SOME(f)) = f copyTarget (*dumb one becomes id if copy use is needed.*)
	in
	 copy(dynSeq, fem, dynSeqFunc)
	end

    (* utility functions for initializing the sizes array *)
    fun sizes j i = CL.mkSubscript(sizesV j, mkInt i)
    fun setSizes j (i, v) = CL.mkAssign(sizes j i, v)

  (* get the number of alive or stable strands strands *)
    fun numStrandsExp snapshot = strandMeth (if snapshot then "num_alive" else "num_stable", [])

  (* code to access state variable
        wrld->outState[i].name
   * or
        wrld->state[i].name
   *)
    fun stateVar spec = let
          fun singleV _ name = CL.mkIndirect(strandMeth("strand", [iV]), "sv_"^name)
          fun dualV U.Shared name = CL.mkIndirect(strandMeth("in_state", [iV]), "sv_"^name)
            | dualV _ name = CL.mkIndirect(strandMeth("local_state", [iV]), "sv_"^name)
          in
            if TargetSpec.dualState spec then dualV else singleV
          end

  (* code fragment to loop over strands
        for (auto ix = wrld->_strands.begin_MODE(), ix != wrld->_strands.end_MODE(), ...) ...
   *
   * where "mode" is either "alive" or "stable".
   *)
    fun forStrands mode stm = CL.mkFor(
          CL.T_Named "auto", [("ix", strandMeth("begin_"^mode, []))],
          CL.mkBinOp(iV, CL.#!=, strandMeth("end_"^mode, [])),
          [CL.mkAssignOp(iV, CL.$=, strandMeth("next_"^mode, [iV]))],
          stm)

  (* code fragment to initialize the axes kinds; the data axis (axis[0]) is given, but we skip it
   * (by convention) if it is scalar. The other axes are the specified domAxisKind.
   *)
    fun initAxisKinds (nrrd, dataAxisKind, nAxes, domAxisKind) = let
        (* nData->axis[0].kind *)
          fun axisKind i = CL.mkSelect(CL.mkSubscript(CL.mkIndirect(nrrd, "axis"), mkInt i), "kind")
          fun init (i, k) = CL.mkAssign (axisKind i, CL.mkVar(Nrrd.kindToEnum k))
          val (firstSpace, dataAxis) = (case dataAxisKind
                 of Nrrd.KindScalar => (0, [])
                  | _ => (1, [init(0, dataAxisKind)])
                (* end case *))
          in
            dataAxis @ List.tabulate(nAxes, fn i => init(i+firstSpace, domAxisKind))
          end

  (* create the body of an output function for dynamic-sequence outputs.  The parameter 'ty'
   * is the element type of the dynamic sequence.  The structure of the
   * function body is:
   *
   *    declarations
   *    compute sizes array for nLengths
   *    allocate nrrd for nLengths
   *    compute sizes array for nData
   *    allocate nrrd for nData
   *    copy data from strands to nrrd
   *)
    (*TODO: needs the path and the copy, but also needs the original seq path! If we have  -1 -2 -> what do we do?*)
    fun genDynOutput (env, snapshot, nAxes, ty, name, kind, num,
		      (elemCTy, nrrdType, axisKind, nElems),
		      loopInfo) = let
     val numString = Int.toString num
     val numElemsV = numElemsV num
     val offsetV = offsetV num
     val nLengthsV = nLengthsV num
     val nDataV = nDataV num
     val spec = Env.target env
     val setSizes = setSizes num
     val stateVar = stateVar spec kind
     val sizesV = "sizes_" ^ numString
          val (nAxes, domAxisKind) = (case nAxes
                 of NONE => (1, Nrrd.KindList)
                  | SOME n => (n, Nrrd.KindSpace)
                (* end case *))
        (* declarations *)
          val sizesDecl = CL.mkDecl(CL.T_Array(sizeTy, SOME(nAxes+1)), sizesV, NONE)
          (* count number of elements (and stable strands) *)
	  (* loop at each seq node*)
	  val name = stateVar name
          val countElems = let
                val nElemsInit = CL.mkDeclInit(CL.uint32, varName numElemsV, CL.mkInt 0)
                val cntElems = fn x => CL.S_Exp(CL.mkAssignOp(numElemsV, CL.+=, seqLength(x)))
                val forLoop = forStrands (if snapshot then "alive" else "stable")
                in [
                  CL.mkComment["count number of elements"],
                  nElemsInit,
                  forLoop (copy(true, ty, name, nElems, elemCTy, loopInfo, SOME(cntElems)))
                ] end
        (* generate code to allocate the nLengths nrrd *)
          val lengthsNrrd = let
                val dimSizes = setSizes (0, CL.mkInt 2)  (* nLengths is 2-element vector *)
                in
                  CL.mkComment["allocate nLengths nrrd"] ::
                  (if #isGrid spec
                    then dimSizes ::
                      List.tabulate (nAxes, fn i =>
                        setSizes (i+1, CL.mkSubscript(CL.mkIndirect(wrldV, "_size"), mkInt(nAxes-i-1)))) @
                      [U.maybeAlloc (env, nLengthsV, Nrrd.tyToEnum Nrrd.TypeInt, nAxes+1, num)]
                    else [
                        dimSizes, setSizes (1, numStrandsExp snapshot),
                        U.maybeAlloc (env, nLengthsV, Nrrd.tyToEnum Nrrd.TypeInt, 2, num)
                      ])
                end
        (* code to check for no data to output (i.e., all of the output sequences are empty) *)
          val checkForEmpty = [
                  CL.mkComment["check for empty output"],
                  CL.mkIfThen(
                    CL.mkBinOp(mkInt 0, CL.#==, numElemsV),
                    CL.mkBlock[
                        CL.mkCall("nrrdEmpty", [nDataV]) (*TODO: Figure out if non-return for emtpy seqs allowed?*)
                      ])
                ]
        (* generate code to allocate the data nrrd *)
          val dataNrrd = if (axisKind = Nrrd.KindScalar)
                then [ (* drop data axis for scalar data by convention *)
                    CL.mkComment["allocate nData nrrd"],
                    setSizes(0, numElemsV),
                    U.maybeAlloc (env, nDataV, Nrrd.tyToEnum nrrdType, 1, num)
                  ]
                else [
                    CL.mkComment["allocate nData nrrd"],
                    setSizes(0, mkInt nElems),
                    setSizes(1, numElemsV),
                    U.maybeAlloc (env, nDataV, Nrrd.tyToEnum nrrdType, 2, num)
                  ]
        (* generate the nLengths copy code *)
          val copyLengths = let
                val pInit = CL.mkDeclInit(CL.T_Ptr CL.uint32, "ip",
                      CL.mkReinterpretCast(CL.T_Ptr(CL.uint32), CL.mkIndirect(nLengthsV, "data")))
                val offsetDecl = CL.mkDeclInit(CL.uint32, varName offsetV, CL.mkInt 0)
                val copyBlk = fn x => CL.mkBlock[
                        CL.mkDeclInit(CL.uint32, "n", seqLength(x)),
                        CL.mkAssign(CL.mkUnOp(CL.%*, CL.mkPostOp(ipV, CL.^++)), offsetV),
                        CL.mkAssign(CL.mkUnOp(CL.%*, CL.mkPostOp(ipV, CL.^++)), nV),
                        CL.S_Exp(CL.mkAssignOp(offsetV, CL.+=, nV))
                      ]
                val mode = if #isGrid spec orelse snapshot
                      then "alive"
                      else "stable"
                in
                  CL.mkComment["initialize nLengths nrrd"] ::
                  pInit ::
                  offsetDecl ::
                  forStrands mode ((copy(true, ty, name, nElems, elemCTy, loopInfo, SOME(copyBlk)))) ::
                  initAxisKinds (nLengthsV, Nrrd.Kind2Vector, nAxes, domAxisKind)
                end
          (* generate the nData copy code *)
	  val targetVar = name
          val copyData = let
                val pInit = CL.mkDeclInit(CL.charPtr, "cp",
					  CL.mkReinterpretCast(CL.charPtr, CL.mkIndirect(nDataV, "data")))
		val copyFunc = fn x => CL.mkAssign(cpV, seqCopy(x, cpV))
                val copyStm = copy(true, ty, targetVar, nElems, elemCTy, loopInfo, NONE)
                val mode = if #isGrid spec
                      then "alive"
                      else if #snapshot spec
                        then "alive"
                        else "stable"
                in
                  CL.mkComment["initialize nLengths nrrd"] ::
                  pInit ::
                  forStrands mode copyStm ::
                  initAxisKinds (nDataV, axisKind, 1, Nrrd.KindList)
                end
        (* the function body *)
          val stms =
                sizesDecl ::
                countElems @
                lengthsNrrd @
                checkForEmpty @
                dataNrrd @
                copyLengths @
                copyData
                
    in
     (*modify name*)
            ([CL.PARAM([], nrrdPtrTy, varName nLengthsV), CL.PARAM([], nrrdPtrTy, varName nDataV)], CL.mkBlock stms)
          end

  (* create the body of an output function for fixed-size outputs.  The structure of the
   * function body is:
   *
   *    declare and compute sizes array
   *    allocate nrrd nData
   *    copy data from strands to nrrd
   *)
    fun genFixedOutput (env, snapshot, nAxes, ty, name, kind, num,
			(elemCTy, nrrdType, axisKind, nElems),
			loopinfo) = let
     val numString = Int.toString num
     val numElemsV = numElemsV num
     val offsetV = offsetV num
     val nLengthsV = nLengthsV num
     val nDataV = nDataV num
     val sizesV = "sizes_" ^ numString
     val setSizes = setSizes num
          val spec = Env.target env
          val stateVar = stateVar spec kind
          val (nAxes, domAxisKind) = (case nAxes
                 of NONE => (1, Nrrd.KindList)
                  | SOME n => (n, Nrrd.KindSpace)
                (* end case *))
          val nDataAxes = if (axisKind = Nrrd.KindScalar) then 0 else 1
        (* generate the sizes initialization code *)
          val initSizes = let
                val dimSizes = let
                 val dcl = CL.mkDecl(CL.T_Array(sizeTy, SOME(nAxes+nDataAxes)), sizesV, NONE)
                      in
                        if (axisKind = Nrrd.KindScalar)
                          then [dcl]
                          else [dcl, setSizes(0, mkInt nElems)]
                      end
                in
                  if #isGrid spec
                    then dimSizes @
                      List.tabulate (nAxes, fn i =>
                        setSizes(i+nDataAxes, CL.mkSubscript(CL.mkIndirect(wrldV, "_size"), mkInt(nAxes-i-1))))
                    else dimSizes @ [setSizes(nDataAxes, numStrandsExp snapshot)]
                end
        (* generate the copy code *)
          val copyCode = let
                val pDecl = CL.mkDeclInit(CL.charPtr, "cp",
					  CL.mkReinterpretCast(CL.charPtr, CL.mkIndirect(nDataV, "data")))
		val targetVar = stateVar name
                val copyBlk = copy(false, ty, targetVar, nElems, elemCTy, loopinfo, NONE) (*fix me....accPath, accTy*)

                val mode = if #isGrid spec orelse snapshot
                      then "alive"
                      else "stable"
                in
                  pDecl :: forStrands mode copyBlk :: initAxisKinds (nDataV, axisKind, nAxes, domAxisKind)
                end
        (* the function body *)
          val stms =
                CL.mkComment["Compute sizes of nrrd file"] ::
                initSizes @
                CL.mkComment["Allocate nData nrrd"] ::
                U.maybeAlloc (env, nDataV, Nrrd.tyToEnum  nrrdType, nAxes+nDataAxes, num) ::
                CL.mkComment["copy data to output nrrd"] ::
                copyCode
          in
            ([CL.PARAM([], nrrdPtrTy, varName nDataV)], CL.mkBlock stms)
          end

    fun gen (env, nAxes, outputs) = let
          val spec = Env.target env
          val mkFunc = if (#exec spec)
                then (fn (funcName, params, body) => CL.D_Func(
                        [], CL.boolTy, [], funcName,
                        CL.PARAM([], RN.worldPtrTy, "wrld")::params,
                        body))
                else let
                  val wrldParam = CL.PARAM([], GenAPI.worldTy spec, "cWrld")
                  val wrldCastStm = CL.mkDeclInit(RN.worldPtrTy, "wrld",
                        CL.mkReinterpretCast(RN.worldPtrTy, CL.mkVar "cWrld"))
                  in
                    fn (funcName, params, body) => CL.D_Func(
                        ["extern \"C\""], CL.boolTy, [], funcName,
                        wrldParam :: params,
                        CL.prependStm(wrldCastStm, body))
                       end
	  (*TODO: This is needlessly complicated*)
	  fun preProcOutput(ty) =
	      let
	       val loopsInfo = Ty.buildOuputputConversionRoutine(ty)
	       val numOutputs = List.length loopsInfo
	      in
	       (loopsInfo, numOutputs)
	      end
	  fun getFnDispatch(env, snapshot, nAxes, ty, name, kind, num, tyinfo, NONE) =
	      let
	       val (ps, CL.S_Block(stms)) =
		   (case ty
		     of Ty.SeqTy(ty', NONE) => genDynOutput(env, snapshot, nAxes, ty', name, kind, num, tyinfo, NONE)
		      | _ => genFixedOutput(env, snapshot, nAxes, ty, name, kind, num, tyinfo, NONE)
		   (*end case*))
	      in
	       (ps, CL.S_Block(stms @ [CL.mkReturn(SOME(CL.mkVar "false"))]))
	      end
	    | getFnDispatch (env, snapshot, nAxes, ty, name, kind, num, tyinfo, SOME(info)) = 
			     (*, path, SOME(pathTy), SOME(fscheck), tyinfo) =*)
			    if Ty.hasDynamicSize ty
	      then genDynOutput(env, snapshot, nAxes, ty, name, kind, num, tyinfo, SOME(info)) (*FIX ME*)
	      else genFixedOutput(env, snapshot, nAxes, ty, name, kind, num, tyinfo, SOME(info))
	  fun oldSystem(env, snapshot, nAxes, ty, name, kind) = let
	   val elemOutTy = (case ty
			     of (Ty.SeqTy(s, NONE)) => s
			      | s => s)
	   val tyinfo = OutputUtil.infoOf (env, elemOutTy)
	  in
	   getFnDispatch(env, snapshot, nAxes, ty, name, kind, 0, tyinfo, NONE)
	  end

	  fun genFnHelper(env, snapshot, nAxes, ty, name, kind) =
	      if Ty.isSingleOutputWithFem(ty)
	      then oldSystem(env, snapshot, nAxes, ty, name, kind)
	      else
	       let
		val (loopInfos, numOutputs) : APITypes.copyOut list * int = preProcOutput ty
	       in
		if numOutputs  = 1
		then oldSystem(env, snapshot, nAxes, #outputTy (List.nth(loopInfos, 0)), name, kind)
		else
		 let
		  
		  (* val itterList = List.tabulate(List.length loopInfo, fn x => x) *)
		  (* val foldList = ListPair.zip(ListPair.zip(ListPair.zip (outTys, loopInfo), fixedSize), itterList) *)
		  fun ouptutFolder (loopinfo : APITypes.copyOut, (params, block)) =
		      let
		       val outTy = #outputTy loopinfo

		       val num = #nrrdNum loopinfo
		       (*correction: This should be the base ty!*)
		       val tyCopyInfo = OutputUtil.infoOf (env, (case outTy
								     of Ty.SeqTy(s, NONE) => s
								      | _ => outTy))
		       val (params', block') = getFnDispatch(env, snapshot, nAxes, outTy, name, kind, num, tyCopyInfo, SOME(loopinfo))
								
							       (* val (elemCTy, nrrdType, axisKind, nElems) =  *)
		      in (params@params', CL.mkBlock([block, block'])) end
			
		 in
		  List.foldr ouptutFolder ([], CL.mkBlock([])) loopInfos
		 end
	       end
          fun getFn snapshot {name, ty, kind} = let
                val funcName = if snapshot
                      then GenAPI.snapshotGet(spec, name)
			       else GenAPI.outputGet(spec, name)
                val (params, body) = genFnHelper(env, snapshot, nAxes, ty, name, kind)
                in
                  mkFunc (funcName, params, body)
                end
          val getFns = List.map (getFn false) outputs
          val getFns = if (#snapshot spec)
                then List.map (getFn true) outputs @ getFns
                else getFns
          in
            if (#exec spec)
              then getFns @ U.genOutput(env, outputs)
              else getFns
          end

  end
