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
(* add cache opts: check(int)(v), load(size, ty), dofs = save(size, ty, dofs) *)
(*expand and mexapnd and the ref (current entry, )*)
(*build the vars and do a rewrite*)
(*add the new statevar and the init - rewrite program*)

(*check = check... (assign node)
if check (cond node)
then dofsTemp=load (true branch)
else loadDofs = ...; dofsTemp = save(loadDofs) (false branch)


dofsTemp = 

function a(checkInt, dofLoad, dofSize, dofType, dstVar) = build the cfg
function checkInt
function buildDofLoad, etcc
 *)
structure IR = MidIR
structure Ty = MidTypes
structure Op = MidOps
val currentSize = ref 0
val numberChecks = ref 0


(*build assigments from op*)
(*calculate dofSize from op*)
(*call things and add cfg*)
(*Borning cfg otherwise*)

(*
NEW IDEA:
0. Continue as I was - need to do image glue and line everything up ^ + MAKE SAVE AN MASSIGN
1. Make a prepass on the program:
*Make a giant tensor and giant array
-Continue translation as is: NO MASSIGN
-Build giant tensor of zero
-build giant array of -1
-init at top of stateInit - manual sergery

2. Do a rewrite where we rewrite as follows:
check -> access statevar and compare for equality
save ->  state = copyPart(a,b) -> save this....
read -> read from the state var and copy specific places

UGG^the save part irks me -> I don't like that 

PROPERTY PLAN WITH MASSIGN AND 

YEAH: SV PLAN MUST GO DUE TO COPY SEMANTICS^^^^^^^^^^^

NEW PLAN:
0. Continue as I was - need to do image glue and line everything up ^ + MAKE SAVE AN MASSIGN
1. Build property and no rewrite pass
2. In cxx:
--header files includes cache.hxx with cachemarker and cache --allocated as neg1 -- http://www.cplusplus.com/reference/cstring/memset/
--push these things down to tree:
---use property, in tree, add the state vars and new tree-types--at least leads to an error
---low to tree we do some translation of the ops

--in sv, based on property, add these to the local -- init handled 
--in the multiassign, save grabs the tensor and just memcopies to the write place
--in the read, grab the pointer
--in the cache -- just array access to int array




OTHER stuff:
Plan:
Okay maybe later:
*maybe worry about F \circ InvT at pos
*try again as thing bugs me
---
*inside fix (comp, at position)
*world pos  (caching in simple)
*cache (new strand var)
*re-enable basis vectorization eval
*maybe remove cell list if unneeded
--ONWAYD!
---
affine
Pos (facet operation so you could build the other one)/ffi cooperation infustructure
---
Firedrake load curvey mesh
--roll something dumb
Patrick opencascade



06107
InsideFix:
-we n
-Need to add new datatype: Image(a, b) | Compose of insides*insides | Mesh(findPos, mesh,..) | Cell(inverse, args) | RefCell
inside(pos) - newimageType = Data {imagesupport, meshvars+func, CellsVars+func == check if inside}  | Comp of newimagetype * newimageType
*
Plan:
--AST: add inside and fix associated things (pp, etc)
--Keep the inside basis vars to be 1
--places where I make prim thing - I am going to check for inside--do the specifics and generate that
---make a func for this
--propogate this and replace in simple
--make key replace in fields
--run tests (okay...)
--Do the damb thing:
---rewrite image type
---rewrite ops on it
---do the analysis
---write generator


Cache:
Plan:
0. Make translate that passes stage...
1.
translate with refs to record number of reals and ref for number 
at the dof assigns:
-compute index
-dof
-send index to tester 
-if true
load
else do dofs and set
2.
With the refs involved, make the varying non-shared statevar
Init the var to zeros a

Make cfg functor to do update or whatever else only or maybe give scope info (sig with scope data) - copy mexpand
Go through update function - add the cache - keep track of total side and number
Allocate non-varying strand of array of int for cache hitting and
--strand vars
--strand init
Go through the update function again to do this again




*)


		       
fun loadOpToIndex ((Op.LoadVoxels(info, n), [img, idx])) =
    let val v = IR.Var.new ("voxelBase", Ty.IntTy) in
     ([(v, (Op.LoadVoxelsBase(info, n), [img, idx]))], v) end
  | loadOpToIndex ((Op.LoadVoxelsWithCtl(info, n, ctl), [img, idx])) =
    let val v = IR.Var.new ("voxelBase", Ty.IntTy) in
     ([(v, (Op.LoadVoxelsBaseWithCtl(info, n, ctl), [img, idx]))], v) end
  (*FIX ME: this needs to handle mulitple possible src1, src2 -- best option is to retrieve them, get index and multiply by this times*numCells before - which you can via access to an array and the structure of the IR before... eep**)
  | loadOpToIndex ((Op.ExtractFemItemN(tys, ty, (FemOpt.ExtractDofs, data), _, _, _, _), [src1, src2, idx])) =  ([], idx)
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
     val dofsCount = FemOpt.findIndexLength opt
     val shape' = dofsCount::shape
     val ten = Ty.TensorTy shape'
    in
     (size shape', ten)
    end
  | dofsInfo _ = raise Fail "impossible"		    
		 
fun buildCondition(calculateInt : IR.assignment list, intVar : IR.var,
		   dofLoad : IR.rhs, dofType : Ty.ty, dofSize : int, dstVar : IR.var) =
    let
     val address = !numberChecks
     val _ = numberChecks := (!numberChecks) + 1
     val sizeAddress = !currentSize
     val _ = currentSize := (!currentSize) + dofSize
     

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
     val (assigns, intVar) = loadOpToIndex(rator, ys')
     val assigns' = List.map (IR.ASSGN o (fn (x,y) => (x, IR.OP(y)))) assigns
     val (dofSize, dofType) = dofsInfo (rator, ys')
     val cfgRet = buildCondition(assigns', intVar, rhs', dofType, dofSize, y')
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


fun translate prog =
    let
     val prog = Trans.translate prog
    in
     (MidCensus.init prog;
     prog)
    end
end
