from __future__ import annotations
import sympy as sp
import numpy as np
from dataclasses import dataclass
from .twistor import line_from_xyzw


def symbolic_two_loop_integrand():
    x1,y1,z1,w1,x2,y2,z2,w2 = sp.symbols("x1 y1 z1 w1 x2 y2 z2 w2")
    n12 = w1*x1 + w2*x2 + y1*z1 + y2*z2
    D1 = w1*x1*y1*z1
    D2 = w2*x2*y2*z2
    X12 = -(x1-x2)*(w1-w2) - (y1-y2)*(z1-z2)
    integrand = -n12/(D1*D2*X12)
    return (x1,y1,z1,w1,x2,y2,z2,w2), sp.factor(integrand), sp.factor(X12), sp.factor(n12)


def collective_scaling(L: int, graph: list[tuple[int,int]], xvals=None, zvals=None, dirs=None):
    """Power counting under y_i=rho*a_i, w_i=rho*b_i for a tree graph.

    This is bookkeeping, not a computation on an integrand: it assumes
      D_ext ~ rho^(2L), each mutual bracket ~ rho, numerator ~ rho^(edges),
      measure ~ rho^(2L-1) d rho,
    so the scalar power is -2L and the form power is -1 for every input.
    The assumptions are verified against the explicit integrand only at
    L=2, by `two_loop_blowup_leading_coefficient`.
    """
    edges = len(graph)
    scalar_power = -(2*L + edges) + edges
    radial_jacobian_power = 2*L - 1
    form_power = scalar_power + radial_jacobian_power
    return {
        "L": L,
        "edges": edges,
        "scalar_integrand_power": scalar_power,
        "radial_jacobian_power": radial_jacobian_power,
        "radial_form_power": form_power,
        "logarithmic": form_power == -1,
    }


def two_loop_blowup_leading_coefficient():
    rho = sp.symbols("rho")
    (x1,y1,z1,w1,x2,y2,z2,w2), I, X12, n12 = symbolic_two_loop_integrand()
    aa,bb,cc,dd = sp.symbols("a b c d", nonzero=True)
    subs = {y1:rho*aa, w1:rho*bb, y2:rho*cc, w2:rho*dd}
    scaled = sp.factor(I.subs(subs))
    lead = sp.factor(sp.limit(rho**4*scaled, rho, 0))
    return scaled, lead


def iterated_residue_order_scan() -> dict:
    """Scan all 4! residue orders on the four small cusp variables.

    Spectator variables are fixed to generic rational values so that the result is
    exact. If the boundary is a normal-crossing divisor with an order-independent
    iterated residue, all 24 orders should agree. Differences are therefore a sharp
    diagnostic that a blow-up/resolution or a different contour construction is needed.
    """
    import itertools
    y1,w1,y2,w2 = sp.symbols('y1 w1 y2 w2')
    x1,z1,x2,z2 = sp.symbols('x1 z1 x2 z2')
    allvars, I, X, n = symbolic_two_loop_integrand()
    aa,yy,zz,ww,bb,qq,rr,ss = allvars
    expr = I.subs({aa: sp.Rational(11,10), zz: sp.Rational(4,5), bb: sp.Rational(17,10), rr: sp.Rational(13,10)})
    small = [yy, ww, qq, ss]
    values = {}
    for perm in itertools.permutations(small):
        r = expr
        for v in perm:
            r = sp.residue(r, v, 0)
        values[tuple(str(v) for v in perm)] = sp.simplify(r)
    uniq = {str(v) for v in values.values()}
    return {
        "n_orders": len(values),
        "n_unique_results": len(uniq),
        "unique_results": sorted(uniq),
        "spectator_values": {"x1": "11/10", "z1": "4/5", "x2": "17/10", "z2": "13/10"},
        "order_results": {" -> ".join(k): str(v) for k, v in values.items()},
    }


# --------------------------------------------------------------------------
# The two-loop form derived from the inequalities
# --------------------------------------------------------------------------
def two_loop_form_from_inequalities():
    """Canonical form of the two-loop negative geometry, from its inequalities alone.

    The region is
        x_i, y_i, z_i, w_i > 0   and   Q = (x1-x2)(w1-w2) + (y1-y2)(z1-z2) > 0,
    that is <AB_1 AB_2> < 0 in the chart of `symbolic_two_loop_integrand`.
    With A = x1 - x2 and c = (y1-y2)(z1-z2), Q = A (w1 - w2) + c is linear in
    w1 and in w2, so the region is cut into pieces on each of which w1 and
    then w2 run over an interval with endpoints rational in the remaining
    variables, over a base that is a product of simplices. The form of each
    piece is the product of the d log forms of its intervals. Nothing about
    the integrand is used.

    Returns (variables, list of (label, density), their sum).
    """
    x1, y1, z1, w1, x2, y2, z2, w2 = v = sp.symbols("x1 y1 z1 w1 x2 y2 z2 w2")
    A, c = x1 - x2, (y1 - y2) * (z1 - z2)
    Q = A * (w1 - w2) + c
    # bases: the orderings of (x1, x2) and of (y1, y2), (z1, z2)
    x_gt, x_lt = 1 / (x2 * (x1 - x2)), 1 / (x1 * (x2 - x1))
    y_gt, y_lt = 1 / (y2 * (y1 - y2)), 1 / (y1 * (y2 - y1))
    z_gt, z_lt = 1 / (z2 * (z1 - z2)), 1 / (z1 * (z2 - z1))
    c_pos = y_gt * z_gt + y_lt * z_lt            # c > 0
    c_neg = y_gt * z_lt + y_lt * z_gt            # c < 0
    pieces = [
        # A > 0: w1 > w2 - c/A
        ("A>0, c>0, 0<w2<c/A, w1>0", x_gt * c_pos * c / (w2 * (c - A * w2)) / w1),
        ("A>0, c>0, w2>c/A, w1>w2-c/A", x_gt * c_pos * A / (A * w2 - c) * A / Q),
        ("A>0, c<0, w2>0, w1>w2-c/A", x_gt * c_neg / w2 * A / Q),
        # A < 0: 0 < w1 < w2 - c/A
        ("A<0, c>0, w2>0, 0<w1<w2-c/A", x_lt * c_pos / w2 * (c - A * w2) / (w1 * Q)),
        ("A<0, c<0, w2>c/A, 0<w1<w2-c/A", x_lt * c_neg * A / (A * w2 - c) * (c - A * w2) / (w1 * Q)),
    ]
    return v, pieces, sp.factor(sp.together(sum(d for _, d in pieces)))


def region_piece(label_index: int, point: dict) -> bool:
    """Is `point` (a dict of the eight variables) in piece number `label_index`?"""
    A = point["x1"] - point["x2"]
    c = (point["y1"] - point["y2"]) * (point["z1"] - point["z2"])
    w1, w2 = point["w1"], point["w2"]
    if A == 0 or c == 0:
        return False
    star = w2 - c / A
    return [A > 0 and c > 0 and w2 < c / A,
            A > 0 and c > 0 and w2 > c / A and w1 > star,
            A > 0 and c < 0 and w1 > star,
            A < 0 and c > 0 and w1 < star,
            A < 0 and c < 0 and w2 > c / A and w1 < star][label_index]


def two_loop_boundary_check():
    """One level of the residue recursion for the two-loop negative geometry.

    The residue of the form at the facet w1 = 0 must be the canonical form
    of that facet, the region
        x_i, y_i, z_i > 0, w2 > 0,  Q0 = -(x1-x2) w2 + (y1-y2)(z1-z2) > 0.
    That form is computed independently, by the same interval rule in w2:
        A > 0, c > 0: 0 < w2 < c/A;   A < 0, c > 0: w2 > 0;   A < 0, c < 0: w2 > c/A.
    Returns (residue, form of the facet, their difference).
    """
    (x1, y1, z1, w1, x2, y2, z2, w2), integrand, _, _ = symbolic_two_loop_integrand()
    residue = sp.factor(sp.limit(w1 * integrand, w1, 0))
    A, c = x1 - x2, (y1 - y2) * (z1 - z2)
    x_gt, x_lt = 1 / (x2 * (x1 - x2)), 1 / (x1 * (x2 - x1))
    y_gt, y_lt = 1 / (y2 * (y1 - y2)), 1 / (y1 * (y2 - y1))
    z_gt, z_lt = 1 / (z2 * (z1 - z2)), 1 / (z1 * (z2 - z1))
    c_pos, c_neg = y_gt * z_gt + y_lt * z_lt, y_gt * z_lt + y_lt * z_gt
    facet = (x_gt * c_pos * c / (w2 * (c - A * w2))
             + x_lt * c_pos / w2
             + x_lt * c_neg * A / (A * w2 - c))
    facet = sp.factor(sp.together(facet))
    return residue, facet, sp.simplify(residue - facet)
