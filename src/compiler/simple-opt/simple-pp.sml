(* simple-pp.sml
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *
 * Pretty printing for the Simple-AST representation.
 *)

structure SimplePP : sig

    val output : (SimpleVar.t -> string ) option -> TextIO.outstream * string * Simple.program  -> unit

    val outputFunc : TextIO.outstream * string * Simple.func_def * (SimpleVar.t -> string ) option  -> unit

  end = struct

    structure PP = TextIOPP
    structure Ty = SimpleTypes
    structure S = Simple
    structure V = SimpleVar

    val indent = PP.Abs 2

    fun ppList ppFn (left, sep, right) (ppStrm, list) = let
          fun sp () = PP.space ppStrm 1
          val string = PP.string ppStrm
          fun pp [] = string right
            | pp [x] = (ppFn(ppStrm, x); string right)
            | pp (x::xs) = (ppFn(ppStrm, x); string sep; sp(); pp xs)
          in
            string left; pp list
          end

  (* print type arguments; we use "#" to denote differentiation arguments, "$" to denote
   * shape arguments, and "%" to denote dimension arguments.
   *)
    fun ppTyArgs (ppStrm, mvs) = let
          val string = PP.string ppStrm
          fun ppTyArg (_, mv) = (case mv
                 of Ty.TY ty => string(Ty.toString ty)
                  | Ty.DIFF NONE => string "#∞"
                  | Ty.DIFF(SOME k) => string("#"^Int.toString k)
                  | Ty.SHAPE shp => string(concat[
                        "$[", String.concatWith "," (List.map Int.toString shp), "]"
                      ])
                  | Ty.DIM d => string("%"^Int.toString d)
                (* end case *))
          in
            ppList ppTyArg ("<", ";", ">") (ppStrm, mvs)
          end

    fun ppVar (ppStrm, x, f) = PP.string ppStrm (case f
						  of NONE => (V.uniqueNameOf x)
						   | SOME(f') => (V.uniqueNameOf x) ^ " (" ^ (f' x) ^ ")"
						(* end case*))

    fun ppVarBind (ppStrm, x, f) = (
          PP.string ppStrm (Ty.toString(V.typeOf x)); PP.space ppStrm 1; ppVar (ppStrm, x, f))

    fun ppFunc (ppStrm, f) = PP.string ppStrm (SimpleFunc.uniqueNameOf f)

    fun ppASTVar (ppStrm, x) = PP.string ppStrm (Var.uniqueNameOf x)

    fun ppVarDecl ppStrm = let
          fun sp () = PP.space ppStrm 1
          val string = PP.string ppStrm
          in
            fn x => (
              PP.openHBox ppStrm;
                case V.kindOf x
                 of V.ConstVar => (string "const"; sp()) (* do we need this? *)
                  | V.InputVar => (string "input"; sp())
                  | V.StrandOutputVar => (string "output"; sp())
                  | _ => ()
                (* end case *);
                string(Ty.toString(V.typeOf x)); sp(); string(V.uniqueNameOf x); string ";";
              PP.closeBox ppStrm)
          end

    fun ppExp (ppStrm, e, f) = let
          fun sp () = PP.space ppStrm 1
          val string = PP.string ppStrm
          fun var x = ppVar (ppStrm, x, f)
          fun ppIndex (ppStrm, NONE) = PP.string ppStrm ":"
            | ppIndex (ppStrm, SOME i) = string(Int.toString i)
	  fun ppVar' (a,b) = ppVar(a, b, f)
          fun pp e = (case e
                 of S.E_Var x => var x
                  | S.E_Lit lit => string (Literal.toString lit)
                  | S.E_Kernel k => string (Kernel.name k)
                  | S.E_Select(x, fld) => (var x; string "."; var fld)
                  | S.E_Apply(f, args) => (string "func ";ppFunc (ppStrm, f); sp(); ppArgs (ppStrm, args))
                  | S.E_Prim(f, [], args, ty) => (string "prim "; ppASTVar(ppStrm, f); sp(); ppArgs (ppStrm, args); string(" ty=");string(Ty.toString(ty)))
                  | S.E_Prim(f, mvs, args, ty) => (
                   string "prim ";ppASTVar(ppStrm, f); ppTyArgs (ppStrm, mvs); sp(); ppArgs (ppStrm, args); string(" ty=");string(Ty.toString(ty))
		  )
                  | S.E_Tensor(es, _) => ppList ppVar' ("[", ",", "]") (ppStrm, es)
                  | S.E_Seq(es, _) => ppList ppVar' ("{", ",", "}") (ppStrm, es)
                  | S.E_Tuple es => ppList ppVar' ("(", ",", ")") (ppStrm, es)
                  | S.E_Project(x, i) => (var x; string "."; string(Int.toString i))
                  | S.E_Slice(x, indices, _) => (
                      var x;
                      ppList ppIndex ("[", ",", "]") (ppStrm, indices))
                  | S.E_Coerce{srcTy, dstTy, x} => (
                      string "("; string(Ty.toString dstTy); string ")"; var x)
                  | S.E_BorderCtl(ctl, x) => (
                      string(BorderCtl.fmt V.uniqueNameOf ctl); string "("; var x; string ")")
                  | S.E_LoadSeq(ty, nrrd) => (
                      string "loadSeq<"; string(Ty.toString ty); string ">(\"";
                      string(String.toString nrrd); string "\")")
                  | S.E_LoadImage(ty, nrrd, _) => (
                      string "loadImage<"; string(Ty.toString ty); string ">(\"";
                      string(String.toString nrrd); string "\")")
		  | S.E_LoadFem(data, x,y) => (
		   string "loadFem("; string ((FemData.toString data) ^ ", "); var x; string ", "; var y;string ")"
		  )
		  | S.E_ExtractFem(v, data) => (
		   string "extractFem("; var v; string ", "; string (FemData.toString data); string ")"
		  )
		  | S.E_ExtractFemItem(v, ty, opt) => (
		   string "extractFemItem("; string (FemOpt.toString opt); string ", "; string (Ty.toString ty); string ", ";
		   var v; string ")"
		  )
		  | S.E_ExtractFemItem2(v1, v2, ty, outTy, opt) => (
		   string "extractFemItem2("; string (FemOpt.toString opt); string ", "; string (Ty.toString ty); string ", ";
		   var v1; string", "; var v2; string " : ";string (Ty.toString ty);  string ")"
		  )
		  | S.E_FemField(v1,v1', v2, ty, field, func) => (
		   string "femField("; var v1; string ","; Option.app var v2; string ","; string (Ty.toString ty); string ","; string (FemOpt.fieldString field);string ",";
		   (Option.map (fn x => (ppFunc (ppStrm, x); sp())) func); string ")"
		  )
		  | S.E_ExtractFemItemN(vars, tys, ty, opt, func) => (
		   string "extractFemItemN("; ppArgs (ppStrm, vars); string ",";string (FemOpt.toString opt); string ",";
		   (Option.map (fn x => (ppFunc (ppStrm, x); sp())) func); string ")"
		  )

                  | S.E_InsideImage(pos, img, s) => (
                      string "insideImage("; var pos; string ","; var img;
                      string ","; string(Int.toString s); string ")")
                  | S.E_FieldFn f => (string "field_fn("; ppFunc (ppStrm, f); string ")")
                (* end case *))
          in
            pp e
          end

    and ppArgs (ppStrm, args) = ppList (fn (a,b) => ppVar(a,b,NONE)) ("(", ",", ")") (ppStrm, args)

    fun ppBlock (ppStrm, [], S.Block{code=[], ...}, f) = PP.string ppStrm "{ }"
      | ppBlock (ppStrm, vars, S.Block{code, ...}, f) = let
          fun sp () = PP.space ppStrm 1
          fun nl () = PP.newline ppStrm
          val string = PP.string ppStrm
          fun var x = ppVar (ppStrm, x, f)
	  fun ppBlock' (a,b, c) = ppBlock(a,b,c,f)
          fun ppStmt stmt = (
                nl();
                case stmt
                 of S.S_Var(x, optE) => (
                      PP.openHBox ppStrm;
                        ppVarBind (ppStrm, x, f);
                        case optE
                         of SOME e => (sp(); string "="; sp(); ppExp(ppStrm, e, f))
                          | NONE => ()
                        (* end case *);
                        string ";";
                      PP.closeBox ppStrm)
                  | S.S_Assign(x, e) => (
                      PP.openHBox ppStrm;
                        var x; sp(); string "="; sp(); ppExp(ppStrm, e, f); string ";";
                      PP.closeBox ppStrm)
                  | S.S_IfThenElse(x, blk, S.Block{code=[], ...}) => (
                      PP.openHBox ppStrm;
                        string "if"; sp(); ppVar(ppStrm, x, f);
                        sp(); ppBlock' (ppStrm, [], blk);
                      PP.closeBox ppStrm)
                  | S.S_IfThenElse(x, blk1, blk2) => (
                      PP.openHBox ppStrm;
                        string "if"; sp(); ppVar(ppStrm, x, f);
                        sp(); ppBlock' (ppStrm, [], blk1);
                      PP.closeBox ppStrm;
                      PP.openHBox ppStrm;
                        sp(); string "else"; sp(); ppBlock' (ppStrm, [], blk2);
                      PP.closeBox ppStrm)
                  | S.S_Foreach(x, xs, blk) => (
                      PP.openHBox ppStrm;
                        string "foreach"; sp(); ppVar(ppStrm, x, f);
                        sp(); string "in"; sp(); ppVar(ppStrm, xs, f);
                        sp(); ppBlock' (ppStrm, [], blk);
                      PP.closeBox ppStrm)
                  | S.S_New(strand, args) => (
                      PP.openHBox ppStrm;
                        string "new"; sp(); string(Atom.toString strand); sp();
                        ppArgs (ppStrm, args); string ";";
                      PP.closeBox ppStrm)
                  | S.S_KillAll => string "kill_all;"
                  | S.S_StabilizeAll => string "stabilize_all;"
                  | S.S_Continue => string "continue;"
                  | S.S_Die => string "die;"
                  | S.S_Stabilize => string "stabilize;"
                  | S.S_Return x => (
                      PP.openHBox ppStrm;
                        string "return"; sp(); ppVar(ppStrm, x, f); string ";";
                      PP.closeBox ppStrm)
                  | S.S_Print args => (
                      PP.openHBox ppStrm;
                        string "print"; sp(); ppArgs (ppStrm, args); string ";";
                      PP.closeBox ppStrm)
                  | S.S_MapReduce[S.MapReduce{result, reduction, mapf, args, source, domain}] => let
                      val S.Func{params, body, ...} = mapf
                      in
                        PP.openHBox ppStrm;
                          ppVarBind(ppStrm, result, f); string "="; sp();
                          string(Reductions.toString reduction);
                          sp(); string "{"; sp();
                          PP.openVBox ppStrm indent;
                            PP.openHBox ppStrm;
                              ppList (fn (a,b) => ppVar(a,b, f)) ("(", ",", ")") (ppStrm, params);
                              sp(); string "=>"; sp();
                            PP.closeBox ppStrm;
                            ppBlock' (ppStrm, [], body);
                          PP.closeBox ppStrm;
                          sp(); ppList (fn (a,b) => ppVar(a,b, f)) ("(", ",", ")") (ppStrm, source::args); nl();
                          string "|";
                          sp(); var source; sp(); string "in"; sp();
                          string(StrandSets.toString domain);
                          sp(); string "};";
                        PP.closeBox ppStrm
                      end
                  | S.S_MapReduce mrs => raise Fail "FIXME: fused map-reduce"
(*
                      fun ppRes (ppStrm, x) = (
                            string(Ty.toString(V.typeOf x)); sp(); string(V.uniqueNameOf x))
                      val S.Func{params, body, ...} = body
                      fun ppRed (ppStrm, r) = string(Reductions.toString r)
                      in
                        PP.openHBox ppStrm;
                          case (results, reductions)
                           of ([x], [r]) => (
                                ppRes (ppStrm, x); sp(); string "="; sp(); ppRed(ppStrm, r))
                            | _ => (
                                ppList ppRes ("(", ",", ")") (ppStrm, results);
                                sp(); string "="; sp();
                                ppList ppRed ("(", ",", ")") (ppStrm, reductions))
                          (* end case *);
                          sp(); string "{"; sp();
                          PP.openVBox ppStrm indent;
                            PP.openHBox ppStrm;
                              ppList ppVar ("(", ",", ")") (ppStrm, params);
                              sp(); string "=>"; sp();
                            PP.closeBox ppStrm;
                            ppBlock (ppStrm, [], body);
                          PP.closeBox ppStrm;
                          sp(); ppList ppVar ("(", ",", ")") (ppStrm, args); nl();
                          string "|";
                          sp(); var source; sp(); string "in"; sp(); string "strands";
                          sp(); string "};";
                        PP.closeBox ppStrm
                      end
*)
                (* end case *))
          in
            PP.openVBox ppStrm (PP.Abs 0);
              string "{";
              PP.openVBox ppStrm indent;
                List.app (fn vdcl => (nl(); ppVarDecl ppStrm vdcl)) vars;
                List.app ppStmt code;
              PP.closeBox ppStrm;
              nl(); string "}";
            PP.closeBox ppStrm
          end

    fun ppParams (ppStrm, params) = let
          fun sp () = PP.space ppStrm 1
          val string = PP.string ppStrm
          in
            ppList
              (fn (_, x) => (string(Ty.toString(V.typeOf x)); sp(); ppVar (ppStrm, x, NONE)))
              ("(", ",", ")")
              (ppStrm, params)
          end

    fun ppFuncDef ppStrm (S.Func{f, params, body}) f' = let
          fun sp () = PP.space ppStrm 1
          fun nl () = PP.newline ppStrm
          val string = PP.string ppStrm
          fun var x = ppVar (ppStrm, x, NONE)
          in
            PP.openHBox ppStrm;
              if SimpleFunc.isDifferentiable f then (string "field"; sp()) else ();
              string "function"; sp();
              string(Ty.toString(SimpleFunc.resultTypeOf f));
              string("#" ^ Int.toString(SimpleFunc.useCount f));
              sp(); ppFunc(ppStrm, f); sp(); ppParams (ppStrm, params);
            PP.closeBox ppStrm;
            nl();
            ppBlock (ppStrm, [], body, f');
            nl()
          end

    fun ppStrand (ppStrm, f, strand) = let
          val S.Strand{
                  name, params, spatialDim, state, stateInit, startM, updateM, stabilizeM
                } = strand
          fun sp () = PP.space ppStrm 1
          fun nl () = PP.newline ppStrm
          val string = PP.string ppStrm
          fun var x = ppVar (ppStrm, x, f)
          fun ppMethod name body = (
                nl(); string name; sp(); ppBlock (ppStrm, [], body, f))
          in
            PP.openHBox ppStrm;
              string "strand"; sp(); string(Atom.toString name); sp();
              ppParams (ppStrm, params);
              case spatialDim
               of SOME d => (sp(); string "in"; sp(); string(Int.toString d ^ "D"))
                | NONE => ()
              (* end case *);
              sp();
            PP.closeBox ppStrm;
            PP.openVBox ppStrm indent;
              string "{";
              List.app (fn vdcl => (nl(); ppVarDecl ppStrm vdcl)) state;
              nl();
              ppBlock (ppStrm, [], stateInit, f);
              Option.app (ppMethod "start") startM;
              ppMethod "update" updateM;
              Option.app (ppMethod "stabilize") stabilizeM;
            PP.closeBox ppStrm;
            nl();
            string "}";  nl()
          end

    fun ppInput f ppStrm = let
          fun sp () = PP.space ppStrm 1
          fun nl () = PP.newline ppStrm
          val string = PP.string ppStrm
          fun var x = ppVar (ppStrm, x, f)
          fun pp (S.INP{var, name, ty, desc, init}) = (
                PP.openHBox ppStrm;
                  string "input"; sp();
                  string (APITypes.toString ty); sp();
                  ppVar (ppStrm, var, f);
                  case desc
                   of SOME desc => string(concat["(\"", desc, "\")"])
                    | NONE => ()
                  (* end case *);
                  case init
                   of S.NoDefault => string ";"
                    | _ => (sp(); string "="; sp(); string(Inputs.initToString init); string ";")
                  (* end case *);
                  nl();
                PP.closeBox ppStrm)
          in
            pp
          end

    fun output f (outS, message, prog) = let
          val S.Program{
                  props, consts, inputs, constInit, globals, funcs, globInit,
                  strand, create, start, update
                } = prog
          val ppStrm = PP.openOut {dst = outS, wid = 120}
          fun sp () = PP.space ppStrm 1
          fun nl () = PP.newline ppStrm
          val string = PP.string ppStrm
          fun ppTopBlock (prefix, SOME blk) = (
                PP.openHBox ppStrm;
                  string prefix; sp();
                  ppBlock (ppStrm, [], blk,f );
                PP.closeBox ppStrm;
                nl())
            | ppTopBlock _ = ()
          fun ppVarDecl prefix x = (
                PP.openHBox ppStrm;
                  string prefix; sp(); string(Ty.toString(V.typeOf x)); sp();
                  ppVar (ppStrm, x, f); string ";";
                PP.closeBox ppStrm;
                PP.newline ppStrm)
          in
            PP.openVBox ppStrm (PP.Abs 0);
              PP.string ppStrm (concat[
                  "/* Simplified Program (after ", message, ") start */"
                ]); nl();
              PP.openHBox ppStrm;
                string "properties:";
                sp();
                string (Properties.propsToString props);
                PP.newline ppStrm;
              PP.closeBox ppStrm;
              List.app (ppVarDecl "const") consts;
              List.app (ppInput f ppStrm) inputs;
              ppTopBlock ("constants", SOME constInit);
              List.app (ppVarDecl "global") globals;
              List.app (fn a => ppFuncDef ppStrm a f) funcs;
              ppTopBlock ("globalInit", SOME globInit);
              ppStrand (ppStrm, f, strand);
              case create
               of Create.Create{dim=SOME d, code} =>
                    ppTopBlock (concat["grid(", Int.toString d, ")"], SOME code)
                | Create.Create{code, ...} => ppTopBlock ("collection", SOME code)
              (* end case *);
              ppTopBlock ("start", start);
              ppTopBlock ("update", update);
              string "/* Program end */"; PP.newline ppStrm;
            PP.closeBox ppStrm;
            PP.closeStream ppStrm
          end

    fun outputFunc (outS, msg, func, f) = let
          val ppStrm = PP.openOut {dst = outS, wid = 120}
          in
            PP.openVBox ppStrm (PP.Abs 0);
              PP.string ppStrm (concat[
                  "/* SimpleAST: ", msg, " */"
                ]); PP.newline ppStrm;
              ppFuncDef ppStrm func;
            PP.closeBox ppStrm;
            PP.closeStream ppStrm
          end

  end
