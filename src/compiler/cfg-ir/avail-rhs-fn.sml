(* avail-rhs-fn.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2016 The University of Chicago
 * All rights reserved.
 *
 * Generic infrastructure for tracking available RHS expressions.  We use this in
 * high-to-mid and mid-to-low phases to do on-the-fly redundant-computation
 * elimination.
 *)

functor AvailRHSFn (

    val phase : string

    structure IR : SSA



  ) : sig

    type rhs = IR.rhs

  (* a table for tracking available applications *)
    type t

  (* create a new table *)
    val new : unit -> t

  (* add a LowIR assignment to the table and return the lhs variable.  We specify a name and
   * type for the lhs.  If the assignment is redundant, then we return the lhs of the previous
   * assignment.  Otherwise, we use the name and type to create a new local variable.
   *)
    val addAssign : t * string * IR.Ty.ty * rhs -> IR.var

  (* adds an assignment to the list of assignments without checking for redundancy *)
    val addAssignToList : t * IR.assign -> unit

  (* get the assignments from the table.  Note that the assignments will be in _reverse_
   * order from how they were added to the table.
   *)
    val getAssignments : t -> IR.assign list

    val assignOp : t * string * IR.Ty.ty * IR.Op.rator * IR.var list -> IR.var
    val assignEin : t * string * IR.Ty.ty * Ein.ein * IR.var list -> IR.var
    val assignCons : t * string   * IR.var list * IR.Ty.ty-> IR.var

    val makeRealLit : t * int list * RealLit.t -> IR.var
    val makeNan : t * int list -> IR.var
    val makeIntLit : t * int -> IR.var

  end = struct

    structure ST = Stats

    val cntNewAssign    = ST.newCounter (phase ^ ":new-assignment")
    val cntReuseAssign  = ST.newCounter (phase ^ ":reuse-assignment")

    datatype rhs = datatype IR.rhs

    structure Tbl = HashTableFn (
      struct
        type hash_key = rhs
        fun addHashVar (x, h) = IR.Var.hash x + h
        fun hashVal rhs = (case rhs
               of IR.GLOBAL x => 0w9941 + IR.GlobalVar.hash x
                | IR.STATE(NONE, x) => 0w7477 + IR.StateVar.hash x
                | IR.STATE(SOME x, y) => 0w7477 + IR.Var.hash x + IR.StateVar.hash y
                | IR.VAR x => 0w7919 + IR.Var.hash x
                | IR.LIT lit => 0w6997 + Literal.hash lit
                | IR.OP(rator, args) => List.foldl addHashVar (IR.Op.hash rator) args
                | IR.CONS(args, _) => List.foldl addHashVar 0w5987 args
                | IR.SEQ(args, _) => List.foldl addHashVar 0w6011 args
                | IR.EINAPP(ein, args) => List.foldl addHashVar (EinUtil.hash ein) args
                | IR.APPLY(f, args) => List.foldl addHashVar (IR.Func.hash f) args
                | IR.MAPREDUCE mrs => List.foldl
                    (fn ((r, f, args), h) =>
                      (List.foldl addHashVar (Reductions.hash r + IR.Func.hash f) args))
                        0w0 mrs
              (* end case *))
        fun sameKey (rhs1, rhs2) = (case (rhs1, rhs2)
               of (IR.GLOBAL x, IR.GLOBAL y) => IR.GlobalVar.same(x, y)
                | (IR.STATE(NONE, x1), IR.STATE(NONE, x2)) => IR.StateVar.same(x1, x2)
                | (IR.STATE(SOME x1, y1), IR.STATE(SOME x2, y2)) =>
                    IR.Var.same(x1, x2) andalso IR.StateVar.same(y1, y2)
                | (IR.VAR x, IR.VAR y) => IR.Var.same(x, y)
                | (IR.LIT a, IR.LIT b) => Literal.same(a, b)
                | (IR.OP(op1, xs), IR.OP(op2, ys)) =>
                    IR.Op.same(op1, op2) andalso ListPair.allEq IR.Var.same (xs, ys)
                | (IR.CONS(xs, _), IR.CONS(ys, _)) => ListPair.allEq IR.Var.same (xs, ys)
                | (IR.SEQ(xs, ty1), IR.SEQ(ys, ty2)) => IR.Ty.same(ty1, ty2) andalso ListPair.allEq IR.Var.same (xs, ys)
                | (IR.EINAPP(ein1, xs), IR.EINAPP(ein2, ys)) =>
                    EinUtil.same(ein1, ein2) andalso ListPair.allEq IR.Var.same (xs, ys)
                | (IR.APPLY(f, xs), IR.APPLY(g, ys)) =>
                    IR.Func.same(f, g) andalso ListPair.allEq IR.Var.same (xs, ys)
                | (IR.MAPREDUCE mrs1, IR.MAPREDUCE mrs2) =>
                    ListPair.allEq
                      (fn ((r1, f1, xs1), (r2, f2, xs2)) =>
                        IR.Func.same(f1, f2) andalso Reductions.same (r1, r2)
                        andalso ListPair.allEq IR.Var.same (xs1, xs2))
                        (mrs1, mrs2)
                | _ => false
              (* end case *))
      end)

    datatype t = TBL of {
        assigns : IR.assign list ref,
        avail : IR.var Tbl.hash_table
     }

    fun new () = TBL{
            assigns = ref[],
            avail = Tbl.mkTable (32, Fail(phase ^ ": AvailRHS"))
          }

    fun addAssign (TBL{assigns, avail}, lhsName, lhsTy, rhs) = (
          case Tbl.find avail rhs
           of SOME y => (ST.tick cntReuseAssign; y)
            | NONE => let
                val lhs = IR.Var.new (lhsName, lhsTy)
                in
                  ST.tick cntNewAssign;
                  Tbl.insert avail (rhs, lhs);
                  assigns := (lhs, rhs) :: !assigns;
                  lhs
                end
    (* end case *))


								 

    fun addAssignToList (TBL{assigns, avail}, assgn) =
        assigns := assgn :: !assigns


    fun assignOp (avail, pre, ty, opss, args) =
        addAssign(avail, pre, ty, IR.OP(opss, args))
    fun assignEin (avail, pre, ty, ein, args) =
        addAssign(avail, pre, ty, IR.EINAPP(ein, args))
    fun assignCons (_, _, [], ty) = raise Fail "empty cons"
      | assignCons (avail, pre, args as x::_, ty) = let
      in
       addAssign(avail, pre, ty, IR.CONS(args, ty))
      end

    fun makeRealLit(t, intlist, rlit) =
	let
	 val nanVar = addAssign(t, "lit", IR.Ty.realTy, IR.LIT(Literal.Real (rlit)))
	 val reved = (List.rev intlist)
	 fun conses([], s, x) = x
	   | conses ((d)::ds, s, x) = conses(ds, d::s, assignCons(t, "lit", List.tabulate(d, fn y => x ), IR.Ty.tensorTy' (d::s)))
	in
	 conses(reved, [], nanVar)
	end						      

    fun makeNan(t, intlist) =
	let
	 val nanVar = addAssign(t, "nan", IR.Ty.realTy, IR.LIT(Literal.Real (RealLit.nan)))
	 val reved = (List.rev intlist)
	 fun conses([], s, x) = x
	   | conses ((d)::ds, s, x) = conses(ds, d::s, assignCons(t, "nans", List.tabulate(d, fn y => x ), IR.Ty.tensorTy' (d::s)))
	in
	 conses(reved, [], nanVar)
	end

    fun makeIntLit(t, int) =
	let
	 val intLit = Literal.intLit int
	in
	 addAssign(t, "intlit", IR.Ty.intTy, IR.LIT(intLit))
	end
			      

    fun getAssignments (TBL{assigns, ...}) = !assigns

  end

