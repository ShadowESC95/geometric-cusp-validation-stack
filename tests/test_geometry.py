import sympy as sp

from geocusp.residue import (
    blowup_intersection_yw,
    chart_transition_source_to_U12,
    source_chart_brackets,
    strict_cusp_blowup_zero_residue,
    strict_cusp_weighted_blowup,
)
from geocusp.negative_geometry import (
    collective_scaling,
    iterated_residue_order_scan,
    symbolic_two_loop_integrand,
    two_loop_blowup_leading_coefficient,
)


def _line(x, y, z, w):
    return sp.Matrix([[1, x, 0, -y], [0, z, 1, w]])


def _bracket(line, i, j):
    unit = sp.eye(4)
    return sp.Matrix.vstack(line, unit[i - 1, :], unit[j - 1, :]).det()


def test_source_chart_brackets_are_determinants():
    x, y, z, w = sp.symbols("x y z w")
    line = _line(x, y, z, w)
    for label, value in source_chart_brackets().items():
        i, j = (int(c) for c in label)
        assert sp.expand(_bracket(line, i, j) - value) == 0


def test_exceptional_cross_ratio():
    res, chi = blowup_intersection_yw()
    x, z0, rho, eta = sp.symbols('x z_0 rho eta', nonzero=True)
    assert sp.simplify(chi - eta * z0 / x) == 0
    assert res == 1 / (x * eta * z0)


def test_strict_cusp_not_simple_radial_divisor():
    pull, residue = strict_cusp_blowup_zero_residue()
    assert sp.simplify(residue) == 0


def test_strict_cusp_weighted_blowup():
    a, C, D, z = sp.symbols("a C D z", nonzero=True)
    residue, cross_ratio = strict_cusp_weighted_blowup()
    assert sp.simplify(residue + 1 / (a * z * D)) == 0
    assert sp.simplify(cross_ratio + 1 + C / (a * D)) == 0


def test_chart_transition_gives_U12_gauge_and_form():
    """The stated map is the row reduction to (1 0 a b; 0 1 c d), and it
    carries dlog x dlog y dlog z0 dlog w to da db dc dd / [a d (ad - bc)]."""
    mapping, inverse, _ = chart_transition_source_to_U12()
    a, b, c, d = sp.symbols("a b c d")
    x, y, z, w = sp.symbols("x y z_0 w")
    source = [inverse[s] for s in (x, y, z, w)]
    line = _line(*source)
    gauge = (line[:, :2].inv() * line).applyfunc(sp.simplify)
    assert gauge == sp.Matrix([[1, 0, a, b], [0, 1, c, d]])
    for symbol, expr in mapping.items():
        assert sp.simplify(expr.subs(inverse, simultaneous=True) - symbol) == 0
    jac = sp.Matrix(source).jacobian([a, b, c, d]).det()
    coefficient = jac / sp.prod(source)
    assert sp.cancel(coefficient - 1 / (a * d * (a * d - b * c))) == 0


def test_two_loop_integrand_matches_brackets():
    """Rebuild the log integrand from 4x4 determinants of the two lines."""
    variables, integrand, mutual, numerator = symbolic_two_loop_integrand()
    x1, y1, z1, w1, x2, y2, z2, w2 = variables
    ab, cd = _line(x1, y1, z1, w1), _line(x2, y2, z2, w2)
    assert sp.expand(sp.Matrix.vstack(ab, cd).det() - mutual) == 0
    cross = (_bracket(ab, 1, 3) * _bracket(cd, 2, 4)
             + _bracket(ab, 2, 4) * _bracket(cd, 1, 3))
    assert sp.expand(cross - numerator) == 0
    poles = sp.Integer(1)
    for line in (ab, cd):
        for i, j in ((1, 2), (2, 3), (3, 4), (1, 4)):
            poles *= _bracket(line, i, j)
    assert sp.cancel(integrand + numerator / (poles * mutual)) == 0


def test_two_loop_collective_scaling_is_logarithmic():
    """The scalar scales exactly as rho^-4, so with the rho^3 Jacobian the
    form has a simple radial pole with a finite, nonzero coefficient."""
    rho = sp.symbols("rho")
    scaled, lead = two_loop_blowup_leading_coefficient()
    assert lead != 0 and not lead.has(rho)
    assert sp.cancel(rho**4 * scaled - lead) == 0
    count = collective_scaling(2, [(0, 1)])
    assert count["scalar_integrand_power"] == -4
    assert count["radial_form_power"] == -1


def test_iterated_residues_depend_on_order():
    """The four cusp divisors are not normal crossing for the mutual factor:
    the 24 residue orders give four distinct values, each on six orders."""
    scan = iterated_residue_order_scan()
    assert scan["n_orders"] == 24
    values = [sp.Rational(v) for v in scan["order_results"].values()]
    expected = {sp.Rational(-2000, 2431), sp.Rational(250, 187),
                sp.Rational(625, 429), sp.Rational(-625, 663)}
    assert set(values) == expected
    assert all(values.count(v) == 6 for v in expected)
