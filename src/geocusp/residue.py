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
    """Return the residue with the radial differential last, and the cross-ratio."""
    x,z0,rho,eta = sp.symbols("x z_0 rho eta", nonzero=True)
    y = rho*eta
    w = rho
    jacobian = sp.Matrix([x, y, z0, w]).jacobian([x, eta, z0, rho]).det()
    residue = rho * jacobian / (x * y * z0 * w)
    chi = sp.simplify((y*z0)/(w*x))
    return sp.factor(residue), sp.factor(chi)


def strict_cusp_blowup_zero_residue():
    """Use U_12 gauge near strict line b=c=d=0,a=alpha.

    The source one-loop form in this chart is da∧db∧dc∧dd/[a*d*(a*d-b*c)].
    The homogeneous scaling of b,c,d yields a vanishing coefficient of
    d rho/rho at a generic point of this exceptional divisor.
    """
    a,b,c,d,rho,B,C,D = sp.symbols("a b c d rho B C D")
    form_coeff = 1/(a*d*(a*d-b*c))
    sub = {b:rho*B, c:rho*C, d:rho*D}
    # The radial Jacobian has degree two on an affine chart of the blow-up.
    pull = sp.simplify(form_coeff.subs(sub) * rho**2)
    residue = sp.simplify(sp.limit(rho * pull, rho, 0))
    return sp.factor(pull), sp.factor(residue)


def strict_cusp_weighted_blowup():
    """Use b=-rho, c=rho*C, d=rho^2*D in the U_12 chart.

    Express the residue in the radial-last order (a, z, D, rho), where
    z is the cross-ratio. This is a local form on the weighted divisor.
    """
    a, rho, C, D, z = sp.symbols("a rho C D z", nonzero=True)
    b, c, d = -rho, rho * C, rho**2 * D
    cross_ratio = sp.factor((b * c - a * d) / (a * d))
    chart = sp.Matrix([a, b, -rho * a * D * (z + 1), d])
    jacobian = chart.jacobian([a, z, D, rho]).det()
    a_chart, b_chart, c_chart, d_chart = chart
    coefficient = jacobian / (a_chart * d_chart * (a_chart * d_chart - b_chart * c_chart))
    return sp.factor(sp.residue(coefficient, rho, 0)), cross_ratio


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
