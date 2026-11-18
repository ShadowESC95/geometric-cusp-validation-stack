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


def two_loop_mutual_facet_check():
    """The residue of the two-loop form on the facet Q = 0 against the
    canonical form of that facet.

    On Q = 0, w1 = w2 - c/A with A = x1 - x2, c = (y1-y2)(z1-z2). The facet is
    the set of the other seven positive variables for which this w1 is
    positive:
        A > 0, c > 0: w2 > c/A;   A > 0, c < 0: w2 > 0;
        A < 0, c > 0: w2 > 0;     A < 0, c < 0: w2 > c/A.
    Its form is computed by the interval rule in w2. With Q = A (w1 - w2) + c,
    dw1 = dQ/A, so the residue is f/A at w1 = w2 - c/A, and the orientation
    induced on the facet agrees with that of the chart for A > 0 and is
    opposite for A < 0: the pieces with A < 0 enter with a minus sign.
    Returns (residue, form of the facet, ratio).
    """
    x1, y1, z1, w1, x2, y2, z2, w2 = sp.symbols("x1 y1 z1 w1 x2 y2 z2 w2")
    A, c = x1 - x2, (y1 - y2) * (z1 - z2)
    n12 = w1 * x1 + w2 * x2 + y1 * z1 + y2 * z2
    D1, D2 = w1 * x1 * y1 * z1, w2 * x2 * y2 * z2
    residue = sp.factor(sp.simplify((n12 / (D1 * D2 * A)).subs(w1, w2 - c / A)))
    x_gt, x_lt = 1 / (x2 * (x1 - x2)), 1 / (x1 * (x2 - x1))
    y_gt, y_lt = 1 / (y2 * (y1 - y2)), 1 / (y1 * (y2 - y1))
    z_gt, z_lt = 1 / (z2 * (z1 - z2)), 1 / (z1 * (z2 - z1))
    c_pos, c_neg = y_gt * z_gt + y_lt * z_lt, y_gt * z_lt + y_lt * z_gt
    shifted = A / (A * w2 - c)                      # d w2 / (w2 - c/A)
    facet = sp.factor(sp.together(x_gt * c_pos * shifted + x_gt * c_neg / w2
                                  - x_lt * c_pos / w2 - x_lt * c_neg * shifted))
    return residue, facet, sp.simplify(residue / facet)


def two_loop_second_level_check():
    """Second level of the residue recursion at every face of codimension two.

    Starting from the facet w1 = 0 (whose form is already matched), take the
    further residue at each of the other seven coordinates and at Q = 0, and
    compare with the form of the face computed from its inequalities. On
    w1 = 0 the region is  -A w2 + c > 0  with A = x1 - x2, c = (y1-y2)(z1-z2):
        w2 = 0:  x1, x2 > 0 free and c > 0;
        x1 = 0:  x2 w2 + c > 0;      x2 = 0:  -x1 w2 + c > 0;
        y1 = 0:  c = -y2 (z1 - z2);  y2 = 0:  c = y1 (z1 - z2);
        z1 = 0:  c = -z2 (y1 - y2);  z2 = 0:  c = z1 (y1 - y2);
        Q = 0:   w2 = c/A, which needs sign(A) = sign(c); since dw2 = -dQ/A
                 the piece with A > 0 enters with a minus sign.
    The symmetries of the form carry w1 = 0 to every coordinate facet, so
    these eight faces represent all faces of codimension two. Returns the
    eight differences, which should vanish.
    """
    x1, y1, z1, w1, x2, y2, z2, w2 = sp.symbols("x1 y1 z1 w1 x2 y2 z2 w2")
    residue, _, _ = two_loop_boundary_check()
    A, c = x1 - x2, (y1 - y2) * (z1 - z2)
    x_gt, x_lt = 1 / (x2 * (x1 - x2)), 1 / (x1 * (x2 - x1))
    y_gt, y_lt = 1 / (y2 * (y1 - y2)), 1 / (y1 * (y2 - y1))
    z_gt, z_lt = 1 / (z2 * (z1 - z2)), 1 / (z1 * (z2 - z1))
    c_pos, c_neg = y_gt * z_gt + y_lt * z_lt, y_gt * z_lt + y_lt * z_gt

    def face(cc, positive, negative, free):
        """Interval rule in w2 for -A w2 + cc > 0; `positive` is the form of cc > 0."""
        return free * (x_gt * positive * cc / (w2 * (cc - A * w2)) + x_lt * positive / w2
                       + x_lt * negative * A / (A * w2 - cc))

    faces = {
        "w1=w2=0": (sp.limit(w2 * residue, w2, 0), c_pos / (x1 * x2)),
        "w1=x1=0": (sp.limit(x1 * residue, x1, 0), (c_pos / w2 + c_neg * x2 / (x2 * w2 + c)) / x2),
        "w1=x2=0": (sp.limit(x2 * residue, x2, 0), (c_pos * c / (w2 * (c - x1 * w2))) / x1),
        "w1=y1=0": (sp.limit(y1 * residue, y1, 0), face(-y2 * (z1 - z2), z_lt, z_gt, 1 / y2)),
        "w1=y2=0": (sp.limit(y2 * residue, y2, 0), face(y1 * (z1 - z2), z_gt, z_lt, 1 / y1)),
        "w1=z1=0": (sp.limit(z1 * residue, z1, 0), face(-z2 * (y1 - y2), y_lt, y_gt, 1 / z2)),
        "w1=z2=0": (sp.limit(z2 * residue, z2, 0), face(z1 * (y1 - y2), y_gt, y_lt, 1 / z1)),
        "w1=Q=0": (sp.limit((w2 - c / A) * residue, w2, c / A), -x_gt * c_pos + x_lt * c_neg),
    }
    return {name: sp.simplify(a - b) for name, (a, b) in faces.items()}


def two_loop_infinity_check():
    """The poles of the two-loop form at infinity in this chart.

    The form also has poles where a coordinate tends to infinity. Each is
    compared with the form of the corresponding face:
        w1 = infinity:        Q > 0 needs x1 > x2, the rest free;
        w1 = 0, w2 = infinity: -A w2 + c > 0 needs x1 < x2, the rest free;
        Q = 0, x1 = x2:       w1 = w2 - c/A runs to infinity; x2 and the rest free.
    Residues are taken as lim v f (v the coordinate, or x1 - x2), which is
    the residue up to the orientation of the upper end of an interval.
    Returns the three ratios, which should be 1.
    """
    (x1, y1, z1, w1, x2, y2, z2, w2), form, _, _ = symbolic_two_loop_integrand()
    x_gt, x_lt = 1 / (x2 * (x1 - x2)), 1 / (x1 * (x2 - x1))
    facet, _, _ = two_loop_boundary_check()
    mutual, _, _ = two_loop_mutual_facet_check()
    cases = {
        "w1=inf": (sp.limit(w1 * form, w1, sp.oo), x_gt / (y1 * z1 * y2 * z2 * w2)),
        "w1=0,w2=inf": (sp.limit(w2 * facet, w2, sp.oo), x_lt / (y1 * y2 * z1 * z2)),
        "Q=0,x1=x2": (sp.limit((x1 - x2) * mutual, x1, x2), 1 / (w2 * x2 * y1 * y2 * z1 * z2)),
    }
    return {name: sp.simplify(a / b) for name, (a, b) in cases.items()}


# --------------------------------------------------------------------------
# The residue recursion at every depth
# --------------------------------------------------------------------------
_REGION_MEMO: dict = {}


def _definite_sign(f):
    """+1 or -1 if the polynomial has that sign on the positive orthant by
    inspection of its coefficients, 0 otherwise."""
    coefficients = sp.Poly(f, *sorted(f.free_symbols, key=str)).coeffs() if f.free_symbols else [f]
    if all(k > 0 for k in coefficients):
        return 1
    if all(k < 0 for k in coefficients):
        return -1
    return 0


def region_form(constraints, variables):
    """The form of {variables > 0, every constraint > 0}, as a density.

    Computed from the inequalities alone. A constraint that factors is
    split into the sign patterns of its factors (a disjoint union). Then a
    variable v occurring in exactly one constraint, alpha v + beta > 0, is
    integrated out by the interval rule:
        alpha > 0, beta > 0:  v > 0              1/v
        alpha > 0, beta < 0:  v > -beta/alpha    alpha/(alpha v + beta)
        alpha < 0, beta > 0:  0 < v < -beta/alpha   beta/(v (alpha v + beta))
    each times the form of the region of the remaining variables on which
    alpha and beta have those signs. Every piece is a positive density on
    its region. Raises NotImplementedError if no variable occurs in exactly
    one constraint.
    """
    key = (frozenset(sp.expand(c) for c in constraints), tuple(variables))
    if key not in _REGION_MEMO:
        _REGION_MEMO[key] = _region_form(list(constraints), tuple(variables))
    return _REGION_MEMO[key]


def _region_form(constraints, variables):
    import itertools

    clean = set()
    for index, c in enumerate(constraints):
        c = sp.expand(c)
        if c == 0:
            return sp.Integer(0)
        coefficient, factors = sp.factor_list(c)
        sign = 1 if coefficient > 0 else -1
        odd = []
        for f, multiplicity in factors:
            d = _definite_sign(f)
            if d:
                sign *= d ** multiplicity
            elif multiplicity % 2:
                odd.append(sp.expand(f))
        if not odd:
            if sign < 0:
                return sp.Integer(0)
            continue
        if len(odd) > 1:
            rest = [k for j, k in enumerate(constraints) if j != index]
            total = sp.Integer(0)
            for eps in itertools.product((1, -1), repeat=len(odd)):
                if sp.prod(eps) == sign:
                    total += region_form(tuple(rest + [e * f for e, f in zip(eps, odd)]), variables)
            return total
        clean.add(sp.expand(sign * odd[0]))
    if any(sp.expand(-c) in clean for c in clean):
        return sp.Integer(0)
    clean = sorted(clean, key=str)
    if not clean:
        return sp.prod([1 / v for v in variables]) if variables else sp.Integer(1)
    single = [v for v in variables if sum(1 for c in clean if c.has(v)) == 1]
    if not single:
        raise NotImplementedError("no variable occurs in exactly one constraint: " + str(clean))
    v = single[0]
    P = next(c for c in clean if c.has(v))
    if sp.degree(P, v) != 1:
        raise NotImplementedError("constraint not linear in " + str(v))
    alpha, beta = sp.expand(P.coeff(v, 1)), sp.expand(P.coeff(v, 0))
    others = [c for c in clean if c is not P]
    rest = tuple(u for u in variables if u != v)
    if beta == 0:
        return region_form(tuple(others + [alpha]), rest) / v
    total = (region_form(tuple(others + [alpha, beta]), rest) / v
             + region_form(tuple(others + [alpha, -beta]), rest) * alpha / (alpha * v + beta)
             + region_form(tuple(others + [-alpha, beta]), rest) * beta / (v * (alpha * v + beta)))
    return sp.cancel(sp.together(total))


def _two_loop_quadric():
    variables = sp.symbols("x1 y1 z1 w1 x2 y2 z2 w2")
    x1, y1, z1, w1, x2, y2, z2, w2 = variables
    return variables, sp.expand((x1 - x2) * (w1 - w2) + (y1 - y2) * (z1 - z2))


def face_constraint(constraint, v):
    """The inequality on the face v = 0 of {constraint > 0, v > 0}: the
    constraint at v = 0, or, if that vanishes identically, the coefficient
    of v. In the second case the face is the limit of the region and not
    the restriction of the inequality."""
    constraint = sp.expand(constraint)
    at_zero = constraint.subs(v, 0)
    return at_zero if at_zero != 0 else sp.expand(constraint.coeff(v, 1))


def _residue(form, v, at=0):
    return sp.cancel(sp.cancel((v - at) * form).subs(v, at)) if form != 0 else sp.Integer(0)


def coordinate_recursion_check():
    """The residue recursion on every coordinate face, to every depth.

    Starting from the region {all eight variables > 0, Q > 0}, whose form
    from the inequalities is the two-loop integrand, every face reached by
    setting coordinates to zero in any order is visited. On each, the
    residue of its form at each remaining coordinate is compared with the
    form of the next face computed from its own inequality. Returns the
    number of faces, the number of residues compared, the number of steps
    at which the inequality vanished identically, and the mismatches.
    """
    variables, Q = _two_loop_quadric()
    _, integrand, _, _ = symbolic_two_loop_integrand()
    states, frontier = {(Q, variables)}, [(Q, variables)]
    compared, degenerate, mismatches = 0, 0, []
    top = sp.cancel(region_form((Q,), variables) - integrand) == 0
    while frontier:
        c, remaining = frontier.pop()
        form = region_form((c,) if c != 1 else (), remaining)
        for v in remaining:
            rest = tuple(u for u in remaining if u != v)
            if c == 1:
                c2 = sp.Integer(1)
            else:
                degenerate += sp.expand(c).subs(v, 0) == 0
                c2 = face_constraint(c, v)
                if c2.is_number:
                    c2 = sp.Integer(1 if c2 > 0 else -1)
            face = sp.Integer(0) if c2 == -1 else region_form((c2,) if c2 != 1 else (), rest)
            compared += 1
            if sp.cancel(_residue(form, v) - face) != 0:
                mismatches.append((c, v))
            if rest and (c2, rest) not in states:
                states.add((c2, rest))
                frontier.append((c2, rest))
    return {"top_form_is_the_integrand": bool(top), "faces": len(states), "residues": compared,
            "degenerate_steps": int(degenerate), "mismatches": mismatches}


def _quadric_face_form(c, remaining, u):
    """The form of the face c = 0 in the chart that eliminates u, or None.

    c = alpha u + beta, so u = -beta/alpha, positive where alpha and beta
    have opposite signs; du = dc/alpha, so the piece with alpha < 0 enters
    with a minus sign."""
    c = sp.expand(c)
    if not c.has(u) or sp.degree(c, u) != 1:
        return None
    alpha, beta = sp.expand(c.coeff(u, 1)), sp.expand(c.coeff(u, 0))
    if beta == 0:
        return None
    rest = tuple(w for w in remaining if w != u)
    return sp.cancel(region_form((alpha, -beta), rest) - region_form((-alpha, beta), rest)), alpha, beta, rest


def quadric_recursion_check():
    """The residue at the quadric on every coordinate face, and one level below.

    On each face with a non-trivial inequality c > 0 and in each chart u in
    which c is linear: (1) the residue of the face's form at c = 0 against
    the form of {c = 0} from its inequalities; (2) the residue of that form
    at each other coordinate v against the form of the quadric face of the
    face v = 0, in the same chart. Steps where the chart degenerates on the
    deeper face (u no longer usable) are skipped and counted.
    """
    variables, Q = _two_loop_quadric()
    states, frontier = {(Q, variables)}, [(Q, variables)]
    first = second = skipped = 0
    mismatches = []
    while frontier:
        c, remaining = frontier.pop()
        form = region_form((c,), remaining)
        for u in remaining:
            chart = _quadric_face_form(c, remaining, u)
            if chart is None:
                continue
            on_quadric, alpha, beta, rest = chart
            first += 1
            if sp.cancel(_residue(form, u, -beta / alpha) - on_quadric) != 0:
                mismatches.append(("quadric", c, u))
            for v in rest:
                c2 = face_constraint(c, v)
                deeper = None if c2.is_number else _quadric_face_form(c2, tuple(w for w in remaining if w != v), u)
                if deeper is None:
                    skipped += 1
                    continue
                second += 1
                if sp.cancel(_residue(on_quadric, v) - deeper[0]) != 0:
                    mismatches.append(("inside", c, u, v))
        for v in remaining:
            c2, rest = face_constraint(c, v), tuple(w for w in remaining if w != v)
            if c2.is_number or not rest or (c2, rest) in states:
                continue
            states.add((c2, rest))
            frontier.append((c2, rest))
    return {"faces": len(states), "quadric_residues": first, "residues_inside": second,
            "skipped": skipped, "mismatches": mismatches}


# ---- zeros and infinities, in every order ------------------------------------
_T_INF = sp.Symbol("t_inf")


def _leading(p, v, at_infinity=False):
    """First non-vanishing coefficient of p in v at v = 0, or the leading
    coefficient at v = infinity. All polynomials here have degree at most one
    in each variable."""
    p = sp.expand(p)
    if p == 0 or not p.has(v):
        return p
    if at_infinity:
        return sp.expand(p.coeff(v, 1))
    at_zero = p.subs(v, 0)
    return at_zero if at_zero != 0 else sp.expand(p.coeff(v, 1))


def _end_residue(form, v, at_infinity=False):
    """lim v f at v = 0 or at v = infinity (the residue up to the orientation
    of the upper end of an interval)."""
    if form == 0:
        return sp.Integer(0)
    if at_infinity:
        return sp.cancel(sp.cancel(form.subs(v, 1 / _T_INF) / _T_INF).subs(_T_INF, 0))
    return sp.cancel(sp.cancel(v * form).subs(v, 0))


def _normalised(c):
    return sp.Integer(1 if c > 0 else -1) if c.is_number else c


def _all_coordinate_faces():
    """Every face reached from the full region by sending coordinates to zero
    or to infinity in any order, as (inequality, remaining variables)."""
    variables, Q = _two_loop_quadric()
    states, frontier = {(Q, variables)}, [(Q, variables)]
    while frontier:
        c, remaining = frontier.pop()
        for v in remaining:
            rest = tuple(u for u in remaining if u != v)
            for at_infinity in (False, True):
                c2 = c if c in (1, -1) else _normalised(_leading(c, v, at_infinity))
                if rest and (c2, rest) not in states:
                    states.add((c2, rest))
                    frontier.append((c2, rest))
    return states


def _face_form(c, remaining):
    return sp.Integer(0) if c == -1 else region_form((c,) if c != 1 else (), remaining)


def full_recursion_check():
    """The residue recursion with every coordinate sent to zero or to
    infinity, in every order, to every depth.

    The chart is compactified as eight intervals [0, infinity]. On each face
    the residue of its form at each end of each remaining coordinate is
    compared with the form of the next face from its own inequality (the
    leading coefficient at that end). Returns counts and mismatches;
    `nonzero` counts the comparisons in which the residue is not zero.
    """
    states = _all_coordinate_faces()
    compared = nonzero = 0
    mismatches = []
    for c, remaining in states:
        form = _face_form(c, remaining)
        for v in remaining:
            rest = tuple(u for u in remaining if u != v)
            for at_infinity in (False, True):
                c2 = c if c in (1, -1) else _normalised(_leading(c, v, at_infinity))
                residue = _end_residue(form, v, at_infinity)
                compared += 1
                nonzero += residue != 0
                if sp.cancel(residue - _face_form(c2, rest)) != 0:
                    mismatches.append((c, v, at_infinity))
    return {"faces": len(states), "residues": compared, "nonzero": int(nonzero), "mismatches": mismatches}


def _pair_form(alpha, beta, rest):
    """Form of the face alpha u + beta = 0 in the chart eliminating u."""
    return sp.cancel(region_form((alpha, -beta), rest) - region_form((-alpha, beta), rest))


def quadric_full_recursion_check():
    """The quadric faces, in every chart, followed to every depth.

    A quadric face in the chart eliminating u is given by (alpha, beta) with
    u = -beta/alpha; its form is that of {alpha > 0, beta < 0} minus that of
    {alpha < 0, beta > 0}. First, on every face of `full_recursion_check`
    and in every chart, the residue of the face's form at the quadric is
    compared with this. Then every quadric face is followed down: at each end
    of each remaining coordinate v the residue is compared with the form
    given by the leading coefficients of alpha and beta at that end, which
    describe where u stays positive in the limit. Returns counts and
    mismatches.
    """
    pairs, first, mismatches = set(), 0, []
    for c, remaining in _all_coordinate_faces():
        if c in (1, -1):
            continue
        form = region_form((c,), remaining)
        for u in remaining:
            chart = _quadric_face_form(c, remaining, u)
            if chart is None:
                continue
            on_quadric, alpha, beta, rest = chart
            first += 1
            if sp.cancel(_residue(form, u, -beta / alpha) - on_quadric) != 0:
                mismatches.append(("quadric", c, u))
            pairs.add((alpha, beta, rest))
    initial, frontier = len(pairs), list(pairs)
    compared = nonzero = 0
    while frontier:
        alpha, beta, rest = frontier.pop()
        form = _pair_form(alpha, beta, rest)
        for v in rest:
            rest2 = tuple(w for w in rest if w != v)
            for at_infinity in (False, True):
                a2, b2 = _leading(alpha, v, at_infinity), _leading(beta, v, at_infinity)
                residue = _end_residue(form, v, at_infinity)
                compared += 1
                nonzero += residue != 0
                if sp.cancel(residue - _pair_form(a2, b2, rest2)) != 0:
                    mismatches.append(("inside", alpha, beta, v, at_infinity))
                if rest2 and (a2, b2, rest2) not in pairs:
                    pairs.add((a2, b2, rest2))
                    frontier.append((a2, b2, rest2))
    return {"quadric_residues": first, "quadric_faces_in_charts": initial, "all_depths": len(pairs),
            "residues": compared, "nonzero": int(nonzero), "mismatches": mismatches}


def _unexpected_poles(form, allowed):
    """Denominator factors of a form that are not simple or not in `allowed`
    (compared up to sign)."""
    if form == 0:
        return []
    _, denominator = sp.fraction(sp.cancel(form))
    _, factors = sp.factor_list(denominator)
    permitted = set()
    for a in allowed:
        a = sp.expand(a)
        if a.is_number:
            continue
        for f, _m in sp.factor_list(a)[1]:
            permitted.add(sp.expand(f))
            permitted.add(sp.expand(-f))
    return [(f, m) for f, m in factors if m != 1 or sp.expand(f) not in permitted]


def pole_check():
    """No face form has a pole other than a simple one on its own boundary.

    For every face of `full_recursion_check` the denominator of its form
    must consist of remaining coordinates and factors of its inequality,
    each to the first power; for every quadric face of
    `quadric_full_recursion_check`, of remaining coordinates and factors of
    alpha and beta. Returns the numbers of forms examined and the offenders.
    """
    offenders, faces, quadric = [], 0, 0
    pairs = set()
    for c, remaining in _all_coordinate_faces():
        form = _face_form(c, remaining)
        faces += 1
        if _unexpected_poles(form, list(remaining) + ([c] if c not in (1, -1) else [])):
            offenders.append(("face", c, remaining))
        if c in (1, -1):
            continue
        for u in remaining:
            chart = _quadric_face_form(c, remaining, u)
            if chart is not None:
                pairs.add((chart[1], chart[2], chart[3]))
    frontier = list(pairs)
    while frontier:
        alpha, beta, rest = frontier.pop()
        quadric += 1
        if _unexpected_poles(_pair_form(alpha, beta, rest), list(rest) + [alpha, beta]):
            offenders.append(("quadric", alpha, beta, rest))
        for v in rest:
            rest2 = tuple(w for w in rest if w != v)
            for at_infinity in (False, True):
                nxt = (_leading(alpha, v, at_infinity), _leading(beta, v, at_infinity), rest2)
                if rest2 and nxt not in pairs:
                    pairs.add(nxt)
                    frontier.append(nxt)
    return {"face_forms": faces, "quadric_forms": quadric, "offenders": offenders}


def leading_singularity_two():
    """An iterated residue of the two-loop form equal to 2.

    Take x1 = 0, x2 = 0: the face is {(y1-y2)(z1-z2) > 0}, two pieces that
    meet along y1 = y2, z1 = z2. Then the residues at y1 = y2 and z1 = z2,
    and the remaining coordinates to zero. Each piece contributes 1. For
    comparison the chain through the eight coordinates gives 1. Returns
    (the form on y1 = y2, the residue through the crossing, the residue
    through coordinates only).
    """
    (x1, y1, z1, w1, x2, y2, z2, w2), form, _, _ = symbolic_two_loop_integrand()
    r = _end_residue(_end_residue(form, x1), x2)
    on_diagonal = sp.cancel(sp.cancel((y1 - y2) * r).subs(y1, y2))
    r = sp.cancel(sp.cancel((z1 - z2) * on_diagonal).subs(z1, z2))
    for v in (y2, z2, w1, w2):
        r = _end_residue(r, v)
    plain = form
    for v in (w1, w2, x1, x2, y1, z1, y2, z2):
        plain = _end_residue(plain, v)
    return sp.factor(on_diagonal), r, plain


def chart_brackets():
    """The chart coordinates as twistor brackets, with Z the identity.

    For the line A = (1, x, 0, -y), B = (0, z, 1, w):
        <AB14> = x,  <AB12> = y,  <AB34> = z,  <AB23> = w,
    and for two lines <AB1 AB2> = -Q. So the hyperplanes x_i = 0 etc. are
    the divisors where a loop line meets an edge of the square, and Q = 0
    is where the two loop lines meet. Returns the four brackets and the
    ratio <AB1 AB2>/Q.
    """
    x, y, z, w, x2, y2, z2, w2 = sp.symbols("x y z w x2 y2 z2 w2")
    A, B = sp.Matrix([1, x, 0, -y]), sp.Matrix([0, z, 1, w])
    A2, B2 = sp.Matrix([1, x2, 0, -y2]), sp.Matrix([0, z2, 1, w2])
    unit = sp.eye(4)
    bracket = lambda i, j: sp.factor(sp.Matrix.hstack(A, B, unit[:, i - 1], unit[:, j - 1]).det())
    Q = (x - x2) * (w - w2) + (y - y2) * (z - z2)
    mutual = sp.Matrix.hstack(A, B, A2, B2).det()
    return ({"14": bracket(1, 4), "12": bracket(1, 2), "34": bracket(3, 4), "23": bracket(2, 3)},
            sp.cancel(mutual / Q), (x, y, z, w))


def orientation_flip_residues():
    """Residues of the quadric-face forms where alpha and beta share a factor.

    There the face is a whole orthant as a set and its orientation flips
    across the common factor. Returns a tally of residue / (orthant form):
    every entry should be +2 or -2.
    """
    import collections

    pairs = set()
    for c, remaining in _all_coordinate_faces():
        if c in (1, -1):
            continue
        for u in remaining:
            chart = _quadric_face_form(c, remaining, u)
            if chart is not None:
                pairs.add((chart[1], chart[2], chart[3]))
    frontier, tally = list(pairs), collections.Counter()
    while frontier:
        alpha, beta, rest = frontier.pop()
        common = sp.gcd(alpha, beta)
        if not common.is_number:
            form = _pair_form(alpha, beta, rest)
            for h, _m in sp.factor_list(common)[1]:
                if h in rest:
                    continue
                t = sorted(h.free_symbols, key=str)[0]
                root = sp.solve(h, t)[0]
                residue = _residue(form, t, root)
                orthant = sp.prod([1 / v for v in rest if v != t])
                tally[sp.cancel(residue / orthant)] += 1
        for v in rest:
            rest2 = tuple(q for q in rest if q != v)
            for at_infinity in (False, True):
                nxt = (_leading(alpha, v, at_infinity), _leading(beta, v, at_infinity), rest2)
                if rest2 and nxt not in pairs:
                    pairs.add(nxt)
                    frontier.append(nxt)
    return dict(tally)
