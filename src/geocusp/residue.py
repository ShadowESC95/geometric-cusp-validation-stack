from __future__ import annotations
import sympy as sp
from dataclasses import dataclass

@dataclass
class ResidueAudit:
    form: sp.Expr
    pullback: sp.Expr
    residue: sp.Expr
    cross_ratio_exceptional: sp.Expr


def source_one_loop_form():
    x,y,z,w = sp.symbols("x y z w")
    # Coefficient of dx∧dy∧dz∧dw in dlog x∧dlog y∧dlog z∧dlog w.
    return x,y,z,w,sp.Integer(1)/(x*y*z*w)


def source_chart_brackets():
    x,y,z,w = sp.symbols("x y z w")
    return {
        "12": y, "23": w, "34": z, "14": x,
        "13": -(w*x+y*z), "24": -1
    }


def blowup_intersection_yw():
    """Blow up the two adjacent massless corners y=w=0 using y=rho*eta,w=rho."""
    x,z0,rho,eta = sp.symbols("x z_0 rho eta", nonzero=True)
    y = rho*eta
    w = rho
    # dlog form. Use symbolic wedge bookkeeping via coordinates (x,rho,eta,z0).
    # dlog x ∧ dlog y ∧ dlog z0 ∧ dlog w
    # = dlog x ∧ (dlog rho + dlog eta) ∧ dlog z0 ∧ dlog rho
    # = - dlog x ∧ dlog eta ∧ dlog z0 ∧ dlog rho.
    # Reorder to dlog rho first: 3 swaps -> sign flips, yielding + dlog rho∧dlog x∧dlog eta∧dlog z0.
    residue = 1/(x*eta*z0)
    chi = sp.simplify((y*z0)/(w*x))
    return sp.factor(residue), sp.factor(chi)


def strict_cusp_blowup_zero_residue():
    """Use U_12 gauge near strict line b=c=d=0,a=alpha.

    The source one-loop form in this chart is da∧db∧dc∧dd/[a*d*(a*d-b*c)].
    Blowing up b,c,d simultaneously yields a vanishing coefficient of d rho/rho:
    the strict collinear line is not itself a simple log divisor of this chart.
    """
    a,b,c,d,rho,B,C,D = sp.symbols("a b c d rho B C D")
    form_coeff = 1/(a*d*(a*d-b*c))
    sub = {b:rho*B, c:rho*C, d:rho*D}
    # Jacobian for (b,c,d) wrt (rho,B,C,D) wedge, choosing d rho ∧ dB ∧ dC ∧ dD.
    # Determinant of 3x3 derivatives after treating rho as radial coordinate is rho^2.
    pull = sp.simplify(form_coeff.subs(sub) * rho**2)
    residue = sp.simplify(sp.limit(rho * pull, rho, 0))
    return sp.factor(pull), sp.factor(residue)


def chart_transition_source_to_U12():
    x,y,z,w,a,b,c,d = sp.symbols("x y z_0 w a b c d")
    mapping = {
        a: -x/z,
        b: -y-x*w/z,
        c: 1/z,
        d: w/z,
    }
    inverse = {
        z: 1/c,
        x: -a/c,
        w: d/c,
        y: a*d/c-b
    }
    # Jacobian of log coordinates to confirm the canonical form transition.
    xs=-a/c; ys=a*d/c-b; zs=1/c; ws=d/c
    J=sp.Matrix([[sp.diff(sp.log(expr), var) for var in (a,b,c,d)] for expr in (xs,ys,zs,ws)]).det()
    return mapping, inverse, sp.factor(J)
