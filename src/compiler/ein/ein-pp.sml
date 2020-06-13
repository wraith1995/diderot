(* ein-pp.sml
 *
 * NOTE: this code probably should be merged into the EinUtil module.
 *
 * This code is part of the Diderot Project (http://diderot-language.cs.uchicago.edu)
 *
 * COPYRIGHT (c) 2015 The University of Chicago
 * All rights reserved.
 *)

structure EinPP : sig

    val toString : Ein.ein -> string

    val expToString : Ein.ein_exp -> string

  end = struct

    structure E = Ein

    val i2s = Int.toString
    val shp2s = String.concatWithMap " " i2s

    fun ti2s (id, E.T) = "t_"^Int.toString id
      | ti2s (id, E.F) = "f_"^Int.toString id

    fun index2s (E.C cx) = concat["'", i2s cx, "'"]
      | index2s (E.V ix) = "i" ^ i2s ix

    fun multiIndex2s [] = ""
      | multiIndex2s alpha = concat ["_{", String.concatWithMap "," index2s alpha, "}"]

    fun delta (a, b) = concat["δ_{", index2s a, ",", index2s b,"}"]
    fun deltaKrn (a, b) = concat["δ_{", index2s a, ",", index2s b,"}"]

    fun deriv [] = ""
      | deriv alpha = concat["∇",multiIndex2s alpha]

    fun expToString e = (case e
           of E.Const r => i2s r
            | E.ConstR r => Rational.toString r
            | E.Tensor(id, []) => "T" ^ i2s id
            | E.Tensor(id, alpha) => concat["T", i2s id, multiIndex2s alpha]
            | E.Zero(alpha) => concat["Z", multiIndex2s alpha]
            | E.Delta ix => delta ix
            | E.Epsilon(ix, jx, kx) => concat["ϵ_{i", index2s ix, ",i", index2s jx, ",i", index2s kx, "}"]
            | E.Eps2(ix, jx) => concat["ϵ_{i", index2s ix, ",i", index2s jx, "}"]
            | E.Field(id, []) => "F" ^ i2s id
            | E.Field(id, alpha) => concat["F", i2s id, multiIndex2s alpha]
            | E.Lift e1 => concat["«", expToString e1, "»"]
            | E.Conv(img, alpha, kern, beta) => let
                val alpha = if null alpha then "" else multiIndex2s alpha
                val beta = if null beta then "" else "dx" ^ multiIndex2s beta
                in
                  concat ["V", i2s img, alpha, "⊛", beta, "H", i2s kern]
            end
	    | E.Identity(i, mu, opt)	=>
	      let
	       val marker = (case opt
		 of NONE => ""
		  | SOME(s, i) => if s
				  then "world_" ^ (i2s i)
				  else "ref_" ^ (i2s i)
	       (* end case*))
		 
	      in
	       concat ["Id_", i2s i, "_{", index2s mu, "}_", marker]
	      end
	    | E.Fem(E.Plain(basis,len, stamp), cell, index, dofs, coeffShape, dxes) =>
	      let
	       val basisLength = Int.toString len
	       val basisString = "Basis(" ^ ")"
	       val alpha = if null coeffShape then "" else multiIndex2s coeffShape
               val beta = if null dxes then "" else "dx" ^ multiIndex2s dxes
	       val s = Option.getOpt(Option.map Stamp.toString stamp, "noFunc")
				  
	      in
	       concat ["femV(",s,", ", i2s cell,", ", i2s index,", ", i2s dofs,")", alpha,  "⊛", basisString, beta]
	      end
	    | E.Fem(E.Invert(basis,len, stamp), cell, index, dofs, coeffShape, dxes) =>
	      let
	       val basisLength = Int.toString len
	       val basisString = "Basis(" ^ ")"
	       val alpha = if null coeffShape then "" else multiIndex2s coeffShape
               val beta = if null dxes then "" else "dx" ^ multiIndex2s dxes
	       val s = Option.getOpt(Option.map Stamp.toString stamp, "NoStamp")
				  
	      in
	       concat ["femInvV(",s,", ", i2s cell,", ", i2s index,", ", i2s dofs,")", alpha,  "⊛", basisString, beta]
	      end
            | E.Partial alpha => "∂/∂x" ^ multiIndex2s alpha
            | E.Apply(e1, e2) => concat [ expToString e1, "@(", expToString e2, ")"]
            | E.Probe(e1, e2) => concat ["Probe(", expToString e1, ",", expToString e2, ")"]
            | E.Comp(e1, es) => let
                fun f(e2, n1) = concat ["[", expToString e2, "{", shp2s n1, "}", "]"]
                in 
                  concat ["Cmp(", expToString e1,")", String.concatWithMap ", " f es] 
                end
            | E.Value ix => "i" ^ i2s ix
            | E.Img(fid, alpha, pos, s) => concat [
                  "V", i2s fid, multiIndex2s alpha, "(", i2s s, ")[",
                  String.concatWithMap "," expToString pos, "]#"]
            | E.Krn(tid, [], dim) => concat["H", i2s tid, "(", Int.toString dim, ")"]
            | E.Krn(tid, betas, dim) => concat[
                  "H", i2s tid, "^{", String.concatWithMap "" deltaKrn betas, "}(", Int.toString dim, ")"
                ]
            | E.OField(E.CFExp es, e1, E.Partial []) => concat [
                  "CFExp[Tids:", String.concatWithMap " ," ti2s es, "](exp:", expToString e1, ")"
                ]
            | E.OField(E.CFExp(es), e1, E.Partial alpha) => concat [
                  expToString (E.OField(E.CFExp(es), e1, E.Partial [])),
                  "dx", multiIndex2s alpha, ")"
                ]
            | E.Poly(E.Tensor(tid, cx), 1, dx) => concat [deriv dx,"(P", i2s tid, multiIndex2s  cx, ")"]
            | E.Poly(E.Tensor(tid, cx), n, dx) => concat [deriv dx,"(P", i2s tid, multiIndex2s  cx, ")^",  i2s n]
            | E.If (E.Var id, e3, e4) =>    concat[ "if(", Int.toString(id), ") then ", expToString e3," else ", expToString e4]
            | E.If (E.Compare(op1, e1, e2), e3, e4) => let
                val c = (case op1
                    of E.GT => ">"
                    | E.LT => "<"
                    | E.GTE => "=>"
                    | E.LTE => "<="
                    | E.EQ => "="
                    (*end case*))
               in
                concat[ "if(", expToString e1, c, expToString e2, ") then ", expToString e3," else ", expToString e4]
               end
            | E.Sum(sx, e) => let
                val sx = List.map
                      (fn (v, lb, ub) => concat ["(i", i2s v, "=", i2s lb, "..", i2s ub, ")"])
                        sx
                in
                  concat ("Σ" :: sx @ ["<(", expToString e, ")>"]@sx)
                end
            | E.Op1(E.PowInt n, e) => concat["(", expToString e , ")^", i2s n]
            | E.Op1(f, e) => let
                val f = (case f
                       of E.Neg => "Neg"
                        | E.Exp => "Exp"
                        | E.Sqrt => "Sqrt"
                        | E.Cosine => "Cosine"
                        | E.ArcCosine => "ArcCosine"
                        | E.Sine => "Sine"
                        | E.ArcSine => "ArcSine"
                        | E.Tangent => "Tangent"
                        | E.ArcTangent => "ArcTangent"
                        | E.Abs => "Abs"
                        | E.Sgn => "Sign"
                        | _ => raise Fail "impossible"
                      (* end case *))
                in
                  concat[f, "(", expToString e, ")"]
                end
           | E.Op2(E.Sub, e1, e2) => concat ["(", expToString e1, ") - (", expToString e2, ")"]
           | E.Op2(E.Div, e1, e2) => concat ["(", expToString e1, ") / ( ", expToString e2, ")"]
           | E.Op2(E.Max, e1, e2) => concat ["Max(", expToString e1, ", ", expToString e2, ")"]
           | E.Op2(E.Min, e1, e2) => concat ["Min(", expToString e1, ", ", expToString e2, ")"]
           | E.Op3(E.Clamp, e1, e2, e3) => concat["Clamp <", expToString e1, ",", expToString e2, ",", expToString e3,">"]
           | E.Opn(E.Add, el) => concat["(", String.concatWithMap " + " expToString el,")"]
           | E.Opn(E.Prod, el) => concat["(", String.concatWithMap " * " expToString el, ")"]
          (* end case *))

    fun toString (Ein.EIN{params, index, body}) = let
          fun paramToString (i, E.TEN(t, shp)) = concat["T", i2s i, "[", shp2s shp, "]"]
            | paramToString (i, E.FLD d) = concat["F", i2s i, "[", i2s d, "]"]
            | paramToString (i, E.KRN) = "H" ^ i2s i
            | paramToString (i, E.IMG(d, shp)) = concat["V", i2s i, "(", i2s d, ")[", shp2s shp, "]"]
	    | paramToString (i, E.INT) = "INT"^ i2s i
	    | paramToString (i, E.FEM(d)) = "FEM(" ^ (Atom.toString( FemData.nameOf d)) ^ ")" ^ i2s i
          val params = String.concatWith "," (List.mapi paramToString params)
          val index = if null index then "" else concat["_{", shp2s index, "}"]
          in
            concat["λ(", params, ")<", expToString body, ">", index]
          end

  end
