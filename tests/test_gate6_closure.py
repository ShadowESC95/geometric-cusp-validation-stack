"""Checks for docs/GATE_6_CLOSURE.md: the Lorentzian chain at a cusp."""
import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from geocusp import lorentzian_cusp as lc

ONE = lambda L: 1
WEIGHTS = {
    "z^0.3": lambda L: mp.exp(mp.mpf("0.3") * L),
    "z^-0.6": lambda L: mp.exp(mp.mpf("-0.6") * L),
    "log^2 + pi^2": lambda L: L**2 + mp.pi**2,
    "log^3": lambda L: L**3,
}


@pytest.fixture(autouse=True)
def precision():
    with mp.workdps(20):
        yield


# ---- step 1: kinematics --------------------------------------------------
@pytest.mark.parametrize("alpha,beta,yperp", [
    (0.03, -0.02, (0.011, 0.007)), (-0.4, 0.25, (0.3, -0.2)), (1e-3, 2e-3, (0.0, 0.0)),
])
def test_light_cone_decomposition_is_exact(alpha, beta, yperp):
    inv = lc.cusp_invariants(alpha, beta, yperp)
    s, w = inv["s"], inv["x02"]
    assert s < 0 and inv["t"] < 0
    assert w == pytest.approx(alpha * beta * s - inv["rho"], abs=1e-12)
    assert inv["x01"] - w == pytest.approx(beta * s, abs=1e-12)
    assert inv["x03"] - w == pytest.approx(-alpha * s, abs=1e-12)


def test_fourth_propagator_tends_to_t_at_the_cusp():
    for scale in (1e-2, 1e-4):
        inv = lc.cusp_invariants(scale, -scale, (scale, 0.0))
        assert inv["x04"] == pytest.approx(inv["t"], rel=10 * scale)


def test_local_identities_are_exact():
    assert all(r == 0 for r in lc.local_identities().values())


# ---- step 2: dictionary with the weighted divisor -------------------------
def test_weighted_chart_brackets():
    a, b, c, d = sp.symbols("a b c d")
    br = lc.weighted_chart_dictionary()
    assert br["12"] == a * d - b * c       # x02^2: weight 2 under (b, c, d) -> (r b, r c, r^2 d)
    assert br["23"] == d                   # x03^2: weight 2
    assert br["14"] == -a                  # x01^2: stays finite along the edge
    assert br["34"] == 1                   # x04^2: frozen
    assert sp.simplify(br["one_plus_z"] - b * c / (a * d)) == 0


def test_twistor_and_position_cross_ratios_agree():
    rng = np.random.default_rng(3)
    frame = rng.normal(size=(4, 4))
    for _ in range(5):
        a, d = rng.normal(size=2)
        b, c = rng.normal(size=2) + 1j * rng.normal(size=2)
        x0, xs = lc.chart_point(a, b, c, d, frame)
        np.testing.assert_allclose(lc.position_cross_ratio(x0, xs),
                                   (b * c - a * d) / (a * d), rtol=1e-9)


def test_strict_cusp_line_is_the_null_edge_through_x2_and_x3():
    rng = np.random.default_rng(5)
    frame = rng.normal(size=(4, 4))
    det = lambda p, q: abs(np.linalg.det(p - q))
    _, xs = lc.chart_point(0.3, 0, 0, 0, frame)
    assert det(xs[2], xs[1]) < 1e-12                       # x2, x3 null separated
    for a in (0.0, 0.5, -1.7):
        x0, _ = lc.chart_point(a, 0, 0, 0, frame)
        assert det(x0, xs[1]) < 1e-10 and det(x0, xs[2]) < 1e-10
        span = np.column_stack([(x0 - xs[1]).ravel(), (xs[2] - xs[1]).ravel()])
        assert np.linalg.svd(span, compute_uv=False)[1] < 1e-10   # on the line x2 x3
    x0, _ = lc.chart_point(0.0, 0, 0, 0, frame)
    np.testing.assert_allclose(x0, xs[1], atol=1e-12)            # a = 0 is the cusp x2


# ---- step 3: sheets visited by the Lorentzian contour ---------------------
@pytest.mark.parametrize("q1,q3,expected", [
    (1e-3, 1e-7, [-1, -2, -1]),     # x01^2 > 0, x03^2 > 0
    (-1e-3, -1e-7, [1, 0]),         # both spacelike
    (-1e-3, 1e-7, [-1, 0]),
    (1e-3, -1e-7, [-1]),
])
def test_argument_of_z_along_the_contour(q1, q3, expected):
    ell, phi = lc.sheet_path(q1, q3)
    assert lc.plateaus(phi) == expected
    assert ell[0] == pytest.approx(0.0, abs=2e-3)          # every path starts at z = -1


# ---- step 4: the two differences ------------------------------------------
RATIO = mp.mpf("1e-12")          # x03^2 / x01^2: how far into the collinear limit


def test_constant_weight_reproduces_the_box_double_logarithm():
    np.testing.assert_allclose(complex(lc.quadrant_sum(ONE, 1e-3, 1e-15)), 2j * np.pi, atol=1e-8)
    np.testing.assert_allclose(complex(lc.quadrant_sum(ONE, 1e-15, 1e-3)), 2j * np.pi, atol=1e-8)


@pytest.mark.parametrize("name", list(WEIGHTS))
@pytest.mark.parametrize("beta", [0.1, 0.5, 0.9])
def test_on_the_edge_segment_the_contour_is_the_unit_circle(name, beta):
    """At x2 + beta p2 with 0 < beta < 1 the collinear logarithm has density
    -2 pi i times the unit-circle mean of F, at every point of the edge."""
    f = WEIGHTS[name]
    mean = complex(lc.circle_mean(f))
    q1 = -mp.mpf(beta)                                   # q1 = beta s with s = -1
    jump = complex(lc.edge_jump(f, q1, abs(q1) * RATIO))
    assert jump == pytest.approx(-2j * np.pi * mean, abs=1e-7 * max(1.0, abs(mean)))


@pytest.mark.parametrize("name", list(WEIGHTS))
@pytest.mark.parametrize("beta", [-0.3, 1.4])
def test_off_the_edge_segment_the_contour_is_contractible(name, beta):
    """On the same null line outside the segment the loop closes up on one
    sheet and the collinear logarithm is absent for every weight."""
    f = WEIGHTS[name]
    q1 = -mp.mpf(beta)
    jump = complex(lc.edge_jump(f, q1, abs(q1) * RATIO))
    assert abs(jump) < 1e-7 * max(1.0, abs(complex(lc.circle_mean(f))))


@pytest.mark.parametrize("name", list(WEIGHTS))
@pytest.mark.parametrize("q1,q3", [(1e-3, 1e-15), (1e-15, 1e-3)])
def test_quadrant_sum_is_the_circle_mean_in_both_collinear_sectors(name, q1, q3):
    """Double-logarithm density at the cusp, summed over the four sign
    quadrants, in the sector collinear to either edge."""
    f = WEIGHTS[name]
    ratio = complex(lc.quadrant_sum(f, q1, q3) / lc.quadrant_sum(ONE, q1, q3))
    mean = complex(lc.circle_mean(f))
    assert ratio == pytest.approx(mean, abs=1e-7 * max(1.0, abs(mean)))


def test_corrections_are_linear_in_the_ratio_of_scales():
    """The error is proportional to x03^2 / x01^2, so the identity is the
    leading term of the collinear limit, not a coincidence at one point."""
    f = WEIGHTS["log^2 + pi^2"]
    mean = complex(lc.circle_mean(f))
    q1 = 1e-4
    errors = [abs(complex(lc.quadrant_sum(f, q1, q1 * r) / lc.quadrant_sum(ONE, q1, q1 * r)) - mean)
              for r in (1e-3, 1e-6, 1e-9)]
    assert errors[1] < errors[0] / 100 and errors[2] < errors[1] / 100
    assert errors[2] < 1e-5


def test_freezing_the_fourth_propagator_only_rescales_the_circle():
    """With x04^2 held at t instead of its edge value, the contour is the
    circle |z| = 1/(1 - beta): the same class, a different representative."""
    f = WEIGHTS["z^0.3"]
    beta = mp.mpf("0.2")
    q1, q3 = -beta, beta * RATIO
    frozen = (q1 * q3 * lc.kernel(f, q1, q3, x04=-1.0)
              - q1 * (-q3) * lc.kernel(f, q1, -q3, x04=-1.0))
    expected = -2j * mp.pi * lc.circle_mean(f, 1 / (1 - beta))
    assert complex(frozen) == pytest.approx(complex(expected), abs=1e-7)
