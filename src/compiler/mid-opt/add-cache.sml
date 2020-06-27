(* add-cache.sml
 *
 * Add cache for dofs
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2020 The University of Chicago
 * All rights reserved.
 *)

structure AddCache : sig
	   val translate : MidIR.program -> MidIR.program
	  end = struct

structure SrcIR = MidIR
structure SrcTy = MidTypes
structure DstIR = MidIR
structure DstTy = MidTypes
structure SrcOp = MidOps
structure DstOp = MidOps	    


val cvtTy =  fn x => x
structure Env = TranslateEnvFn (
 struct
 structure SrcIR = SrcIR
 structure DstIR = DstIR
 val cvtTy = cvtTy
 end)		    

structure IR = MidIR
structure Ty = MidTypes
structure Op = MidOps
structure Var = IR.Var


(*FIXME: allow more control of these so we can preserve caching for repeated evaluations downstream.

Direction: F and x has a cache assocaited to it right now. We want to group F and several x.

X
x+h
x+h'

So the question is when should you do that:
Use the position:

When to use a previous check:
When it is exactly the same.
When there is a dependence that doesn't have anything to do
When we can trace it back to a previous pos, we use the previous check.
Get bindings: []... Deps

Function:
map from vars to check info.
if none, we trace it back to potential check info (later...)

Check: if there are dependencies on it. If no dependencies, natural thing. Otherwise, if there are, scann for previous ints involved (but no findPos). Then try.

Strats:
-One cache per. 
-One per unique instruction
-determine locality
---Phi means check both if they are used.
---meshPos -- same


Plan:
--make sure we don't create uneeded cache checks and saves and such

if same as last one, just use previous load and save the load...




*)
val currentSize = ref 0
val numberChecks = ref 0


fun sameTriplet((x1, x2, x3), (y1, y2, y3)) = Var.same(x1, y1) andalso Var.same(x2, y2) andalso Var.same(x3, y3)
fun hashTriplet(x, y, z) = Var.hash x + 0w997 * Var.hash y + 0w1999 * Var.hash z

val table : ((MidIR.var * MidIR.var * MidIR.var), (int * int * SrcTy.ty * int)) HashTable.hash_table = HashTable.mkTable (hashTriplet, sameTriplet)  (256, Fail "triplet not found")
datatype lookupInfo = CHECK of int * int * SrcTy.ty * int
fun findCheckInfo((x,y,z), (dofTy, dofSize)) =
    (case HashTable.find table (x, y, z)
      of SOME(a,b,c,d) => CHECK(a,b,c,d)
       | NONE =>
	 let
	  val addr = !numberChecks
	  val _ = numberChecks := (!numberChecks) + 1
	  val sizeAddr = !currentSize
	  val _ = currentSize := (!currentSize) + dofSize
						    
	  val _ = HashTable.insert table ((x,y,z), (addr, sizeAddr, dofTy, dofSize))
	 in
	  CHECK((addr, sizeAddr, dofTy, dofSize))
	 end
    )




		       
fun loadOpToIndex ((Op.LoadVoxels(info, n), [img, idx])) =
    let val v = IR.Var.new ("voxelBase", Ty.IntTy) in
     ([(v, (Op.LoadVoxelsBase(info, n), [img, idx]))], v, (img, img)) end
  | loadOpToIndex ((Op.LoadVoxelsWithCtl(info, n, ctl), [img, idx])) =
    let val v = IR.Var.new ("voxelBase", Ty.IntTy) in
     ([(v, (Op.LoadVoxelsBaseWithCtl(info, n, ctl), [img, idx]))], v, (img, img)) end
  (*FIX ME: this needs to handle mulitple possible src1, src2 -- best option is to retrieve them, get index and multiply by this times*numCells before - which you can via access to an array and the structure of the IR before... eep**)
  | loadOpToIndex ((Op.ExtractFemItemN(tys, ty, (FemOpt.ExtractDofs, data), _, _, _, _), [src1, src2, idx])) =  ([], idx, (src1, src2))
  | loadOpToIndex _ = raise Fail "impossible: dofs from non-dofs"


fun loadsDofs (IR.OP(Op.LoadVoxels(info, n), [img, idx])) = true
  | loadsDofs (IR.OP(Op.LoadVoxelsWithCtl(info, n, ctl), [img, idx])) = true
  | loadsDofs (IR.OP(Op.ExtractFemItemN(tys, ty, (FemOpt.ExtractDofs, data), _, _, _, _), [src1, src2, idx])) =  true
  | loadsDofs _ = false
fun power (x, 0) = 1  
  | power (x, n) = x * power(x, n-1);
fun size lst = List.foldr (op*) 1 lst
fun dofsInfo ((Op.LoadVoxels(info, n), [img, idx])) =
    let
     val shape = ImageInfo.voxelShape info
     val dim = ImageInfo.dim info
     val lastSize = power(n, dim)
     val shape' = shape@[lastSize]
     val ten = Ty.TensorTy shape'
    in
     (size shape', ten)
    end
  | dofsInfo ((Op.LoadVoxelsWithCtl(info, n, ctl), [img, idx])) =
    let
     val shape = ImageInfo.voxelShape info
     val dim = ImageInfo.dim info
     val lastSize = power(n, dim)
     val shape' = shape@[lastSize]
     val ten = Ty.TensorTy shape'
    in
     (size shape', ten)
    end    
  | dofsInfo ((Op.ExtractFemItemN(tys, ty, opt as (FemOpt.ExtractDofs, data), _, _, _, _), [src1, src2, idx])) =
    let
     val shape = FemOpt.findTargetShape opt
     val ten = Ty.TensorTy shape
    in
     (size shape, ten)
    end
  | dofsInfo _ = raise Fail "impossible"		    
		 
fun buildCondition(calculateInt : IR.assignment list, intVar : IR.var,
		   dofLoad : IR.rhs, dofType : Ty.ty, dofSize : int, dstVar : IR.var,
		   address :int, sizeAddress : int) =
    let
     (* val address = !numberChecks *)
     (* val _ = numberChecks := (!numberChecks) + 1 *)
     (* val sizeAddress = !currentSize *)
     (* val _ = currentSize := (!currentSize) + dofSize *)
     

     val checkBool = IR.Var.new ("check", Ty.BoolTy)
     val checkBoolAssign = IR.ASSGN(checkBool, IR.OP(Op.Check(address), [intVar]))
     val loadCfg = IR.CFG.mkBlock (calculateInt@[checkBoolAssign])

     
     val cond = IR.Node.mkCOND checkBool
     val _ = IR.Node.addEdge(IR.CFG.exit loadCfg, cond)

     val cachedDofs = IR.Var.new("cachedDofs", dofType)
     val cachedDofsLoad = IR.ASSGN(cachedDofs, IR.OP(Op.Load(address, sizeAddress, dofType, dofSize), []))
     val cachedDofBranch = IR.CFG.mkBlock [cachedDofsLoad]
     val _ = IR.Node.setTrueBranch (cond, IR.CFG.entry cachedDofBranch)
     val _ = IR.Node.setPred(IR.CFG.entry cachedDofBranch, cond)

     val intermediateDof = IR.Var.new ("loadedDofs", dofType)
     val savedDofs = IR.Var.new ("savedDofs", dofType)
     val loadDofs = IR.ASSGN((intermediateDof, dofLoad))
     val saveDofs = IR.MASSGN(([], (IR.OP(Op.Save(address, sizeAddress, dofType, dofSize), [intVar, intermediateDof]))))
     val saveDofsBranch = IR.CFG.mkBlock [loadDofs, saveDofs]
     val _ = IR.Node.setFalseBranch (cond, IR.CFG.entry saveDofsBranch)
     val _ = IR.Node.setPred(IR.CFG.entry saveDofsBranch, cond)

     (* val phiDofs = IR.Var.new("phiDofs", dofType) *)
     val phi = [(dstVar, [SOME(cachedDofs), SOME(intermediateDof)])]
     val joinNd = IR.Node.mkJOIN phi
     val _ = IR.Node.addEdge(IR.CFG.exit cachedDofBranch, joinNd)
     val _ = IR.Node.addEdge(IR.CFG.exit saveDofsBranch, joinNd)
     val _ = IR.Node.setEdgeMask(joinNd, [false, false])
     val cfg = IR.CFG{entry=IR.CFG.entry loadCfg, exit=joinNd}
    in
     cfg
    end

fun expandDofs(env, (y, IR.OP(rator, ys)) : IR.var * IR.rhs) =
    let
     val y' = Env.rename (env, y)
     val ys' = Env.renameList(env, ys)
     val rhs' = IR.OP(rator, ys')
     val (assigns, intVar, (v1, v2)) = loadOpToIndex(rator, ys')
     val assigns' = List.map (IR.ASSGN o (fn (x,y) => (x, IR.OP(y)))) assigns
     val (dofSize, dofType) = dofsInfo (rator, ys')
     val cfgRet = (case  findCheckInfo((intVar, v1, v2), (dofType, dofSize))
		    of CHECK(address, sizeAddress, _, _) =>
		       buildCondition(assigns', intVar, rhs', dofType, dofSize, y', address, sizeAddress)
		  (* end case*))
    in
     cfgRet
    end
  | expandDofs _ = raise Fail "impossible"

fun expand(env, (y, rhs) : (IR.var * IR.rhs), b) =
    let
     fun expand' (env, (y, rhs)) = let
      fun assign rhs = [DstIR.ASSGN(Env.rename (env, y), rhs)]
     in
      case rhs
       of SrcIR.GLOBAL x => assign (DstIR.GLOBAL(Env.renameGV(env, x)))
        | SrcIR.STATE(NONE, fld) => assign (DstIR.STATE(NONE, Env.renameSV(env, fld)))
        | SrcIR.STATE(SOME x, fld) =>
          assign (DstIR.STATE(SOME(Env.rename(env, x)), Env.renameSV(env, fld)))
        | SrcIR.VAR x => assign (DstIR.VAR(Env.rename(env, x)))
        | SrcIR.LIT lit => assign (DstIR.LIT lit)
        | SrcIR.OP(rator, args) =>
	  assign (DstIR.OP(rator, Env.renameList(env, args)))
        | SrcIR.CONS(args, ty) => assign (DstIR.CONS(Env.renameList(env, args), cvtTy ty))
        | SrcIR.SEQ(args, ty) => assign (DstIR.SEQ(Env.renameList(env, args), cvtTy ty))
        | SrcIR.EINAPP(rator, args) => assign (DstIR.EINAPP(rator, Env.renameList(env, args)))
        | SrcIR.APPLY(f, args) =>
          assign (DstIR.APPLY(Env.renameFV(env, f), Env.renameList(env, args)))
        | _ => raise Fail("bogus rhs for ASSIGN: " ^ SrcIR.RHS.toString rhs)
		     (* end case *)
     end
    in
     if b andalso loadsDofs(rhs)
     then expandDofs(env, (y, rhs))
     else IR.CFG.mkBlock (expand' (env, (y, rhs)))
    end

fun mexpand (env, (ys, rhs), b) = let
 fun massign rhs = let
  val nd = DstIR.Node.mkMASSIGN(Env.renameList(env, ys), rhs)
 in
  DstIR.CFG{entry=nd, exit=nd}
 end
 fun mkOP (rator, xs) = massign(DstIR.OP(rator, Env.renameList(env, xs)))
in
 case rhs
  of SrcIR.OP(SrcOp.EigenVecs2x2, xs) => mkOP (DstOp.EigenVecs2x2, xs)
   | SrcIR.OP(SrcOp.EigenVecs3x3, xs) => mkOP (DstOp.EigenVecs3x3, xs)
   | SrcIR.OP(SrcOp.KillAll, []) => mkOP (DstOp.KillAll, [])
   | SrcIR.OP(SrcOp.StabilizeAll, []) => mkOP (DstOp.StabilizeAll, [])
   | SrcIR.OP(SrcOp.Print tys, xs) => mkOP (DstOp.Print(List.map cvtTy tys), xs)
   | SrcIR.MAPREDUCE mrs => let
    val mrs = List.map
                (fn (r, f, xs) => (r, Env.renameFV(env, f), Env.renameList(env, xs)))
                mrs
   in
    massign (DstIR.MAPREDUCE mrs)
   end
   | _ => raise Fail("bogus rhs for MASSIGN: " ^ SrcIR.RHS.toString rhs)
		(* end case *)
end

(*rewrite just touches the key ones*)
structure Trans = TranslateScopedFn(
 struct
 open Env
 type scope = bool
 val expand = expand
 val mexpand = mexpand

 val constInitScope = false
 val funcScope = false
 val globInitScope = false
 val globStart = false
 val globUpdate = false
 val createScope = false


 val strandInit = false
 val strandStart = false
 val strandUpScope = true
 val strandStabScope = false
 end
)

structure Promote = PromoteFn (DstIR)
fun translate prog =
    let
     val prog = Trans.translate prog
     val IR.Program{props, consts, inputs, constInit, globals, funcs, globInit, strand, create, start, update} = prog
     val props = if (!currentSize) <> 0
		 then Properties.UsesCache(!numberChecks, !currentSize)::props
		 else props
     val prog = IR.Program{props=props, consts=consts, inputs=inputs, constInit=constInit, globals=globals, funcs=funcs,
			  globInit=globInit, strand=strand, create=create, start=start, update=update}
    in
     (MidCensus.init prog;
      Promote.transform prog)
    end
end
