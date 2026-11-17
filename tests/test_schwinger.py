"""Checks for the parametric representation, docs/GATE_6_CONTOUR.md section 11."""
import numpy as np
import pytest

from geocusp import feynman_contour as fc
from geocusp import schwinger as sw


@pytest.fixture(scope="module")
def square():
    return fc.Contour().x


def test_rotated_invariants_lie_in_a_proper_cone(square):
    """|Y_i| >= sin(phi) (mu + |x0 - x_i|_E^2) and 0 <= arg Y_i <= pi - 2 phi."""
    rng = np.random.default_rng(0)
    mu = 0.3
    for phi in (0.05, 0.6, 1.3):
        for _ in range(2000):
            x = rng.normal(size=4) * 10 ** rng.uniform(-1, 1.5)
            for xi in square:
                Y = mu + np.sum((x[1:] - xi[1:]) ** 2) - np.exp(-2j * phi) * (x[0] - xi[0]) ** 2
                euclid = mu + np.sum((x - xi) ** 2)
                assert abs(Y) >= np.sin(phi) * euclid * (1 - 1e-12)
                assert -1e-12 <= np.angle(Y) <= np.pi - 2 * phi + 1e-12


@pytest.mark.parametrize("powers", [(1, 1, 1, 1), (1.2, 0.8, 1.2, 0.8)])
@pytest.mark.parametrize("phi", [np.pi / 2, 1.0])
def test_parametric_representation_at_positive_rotation(square, powers, phi):
    """Direct four-dimensional integration against the Feynman-parameter formula."""
    left, left_err = sw.lhs_position_space(square, powers, 0.5, phi, log2_samples=19)
    right, right_err = sw.rhs_parametric(square, powers, 0.5, phi, log2_samples=17)
    assert abs(left - right) < 5 * (left_err + right_err) + 2e-3 * abs(right)
    assert left_err < 5e-3 * abs(right)                 # the comparison has teeth


def test_right_side_is_continuous_down_to_the_minkowski_limit(square):
    powers = (1.2, 0.8, 1.2, 0.8)
    limit, _ = sw.minkowski_limit(square, powers, 0.5)
    gaps = [abs(sw.rhs_parametric(square, powers, 0.5, phi)[0] - limit) for phi in (0.1, 0.01, 0.001)]
    assert gaps[0] > 5 * gaps[1] > 25 * gaps[2]
    assert abs(limit.real) < 1e-9 and limit.imag < 0     # -i pi^2/Gamma(4) times a positive mean


@pytest.mark.parametrize("region", ["--", "+-"])
def test_no_invariant_vanishes_between_the_rotated_integral_and_the_contour(region):
    C = fc.tuned_contour(region, 1.1)
    bound = fc.certificate(C)["excursion_bound"]
    for mu in (0.5, 0.001):
        result = sw.family_check(C, 2500, mu=mu, seed=2)
        assert result["min_abs_over_scale"] > 0.5
        assert result["max_excursion"] < bound


# ---- the contour in extra dimensions, and step 6 --------------------------
def test_real_shift_is_not_the_feynman_prescription_but_the_extended_contour_is():
    """With the four-dimensional field and a real shift mu, the surface
    x0i^2 = mu is passed on the wrong side at some points once mu > delta0.
    Displacing the extra coordinate as well removes every such point."""
    C = fc.tuned_contour("--", 1.1)
    small = fc.shifted_surface_sides(C, 0.5 * C.delta0, 8000, seed=0)
    assert small["real_shift"] == 0 and small["extended"] == 0
    large = fc.shifted_surface_sides(C, 5 * C.delta0, 8000, seed=0)
    assert large["real_shift"] > 100                    # the earlier claim fails here
    assert large["extended"] == 0
    for mu in (1.0, 4.0, 50.0):
        assert fc.shifted_surface_sides(C, mu, 4000, seed=1)["extended"] == 0


@pytest.mark.parametrize("region", ["--", "+-", "-+"])
def test_contour_in_five_dimensions(region):
    C = fc.tuned_contour(region, 1.1)
    result = fc.survey_extended(C, 8000, seed=4)
    assert result["near_count"] > 3000 and result["largest_mu_near"] > 100
    assert result["max_im_near"] <= 0.0
    assert result["max_excursion"] <= fc.certificate(C)["excursion_bound"]
    assert result["min_ratio"] > 0


@pytest.mark.parametrize("region", ["--", "+-", "-+"])
def test_step_6_identity_nonvanishing_and_uniform_bound(region):
    """Lemma 10: Im(e^{i phi} Y) = sin(phi)(|n|^2 - lam^2|V|^2) - 2 lam cos(phi) V.n,
    |V| <= J |n| with J < 1, no Y vanishes for phi > 0, and at lam = 1 the
    lower bound of Lemma 6 survives uniformly in phi."""
    C = fc.tuned_contour(region, 1.1)
    result = sw.step6_check(C, 4000, seed=3)
    assert result["identity"] < 1e-9
    assert result["v_over_n"] < result["jacobian_bound"] < 1
    assert result["nonvanishing"] > 0.5
    assert result["uniform_count"] > 1000 and result["uniform"] > 0.15


def test_step_6_also_holds_at_larger_rotation_angles():
    C = fc.tuned_contour("--", 1.1)
    result = sw.step6_check(C, 2500, seed=5, phi_max=0.3)
    assert result["identity"] < 1e-9 and result["nonvanishing"] > 0.5


# ---- Feynman parameters to Mellin-Barnes, the timelike diagonal, non-integer D
def test_dirichlet_exponents_of_the_mellin_barnes_derivation():
    """After the Mellin-Barnes split of (a1 a3 S + a2 a4 T)^(eta-2) the
    parameter integral is a Dirichlet integral whose exponents sum to 2 eta,
    and the contour Re w = -1 - Re p + eta/2 keeps every exponent positive."""
    import sympy as sp

    p, w, eta = sp.symbols("p w eta")
    b13 = 1 + p + w                       # exponents of a1, a3 (power a - 1 = p, plus w)
    b24 = (1 - p) + (eta - 2 - w)         # exponents of a2, a4 (power -p, plus eta - 2 - w)
    assert sp.simplify(2 * b13 + 2 * b24 - 2 * eta) == 0
    on_contour = {w: -1 - p + eta / 2}
    assert sp.simplify(b13.subs(on_contour) - eta / 2) == 0
    assert sp.simplify(b24.subs(on_contour) - eta / 2) == 0
    # the split itself needs -(2 - eta) < Re w < 0: true for |Re p| < 1 and small eta
    for pv in (-0.9, 0.0, 0.9):
        centre = -1 - pv + 0.05 / 2
        assert -(2 - 0.05) < centre < 0


def test_beta_integral_behind_the_mellin_barnes_split():
    """int_0^oo x^(c-1) (1+x)^(-lam) dx = Gamma(c) Gamma(lam-c)/Gamma(lam), 0 < c < lam."""
    import mpmath as mp

    for c, lam in ((0.3, 1.6), (1.1, 1.95)):
        c, lam = mp.mpf(c), mp.mpf(lam)
        lhs = mp.quad(lambda v: mp.exp(c * v) * (1 + mp.exp(v)) ** (-lam), [-mp.inf, 0, mp.inf])
        assert abs(lhs - mp.gamma(c) * mp.gamma(lam - c) / mp.gamma(lam)) < mp.mpf("1e-12")


def test_mellin_barnes_form_is_continuous_up_to_the_timelike_point():
    import mpmath as mp

    from geocusp import lorentzian_edges as le

    with mp.workdps(15):
        real = le.regulated_double_pole(0.2, 0.4, 2.5)
        assert abs(le.regulated_double_pole_complex(0.2, 0.4, 2.5, 0) - real) < 1e-10
        at_pi = le.regulated_double_pole_complex(0.2, 0.4, 1.0, mp.pi)
        gaps = [abs(le.regulated_double_pole_complex(0.2, 0.4, 1.0, mp.pi - d) - at_pi)
                for d in (0.1, 0.01, 0.001)]
    assert gaps[0] > 5 * gaps[1] > 25 * gaps[2] and gaps[2] < 1e-3
    assert abs(at_pi.imag) > 1                         # the absorptive part appears


def test_timelike_diagonal_rotated_parameter_integral_equals_mellin_barnes():
    """With S = |s| e^{i pi} the Feynman denominator vanishes inside the
    simplex. On the rotated contour it does not, and the integral equals the
    Mellin-Barnes form at the same point."""
    import mpmath as mp

    from geocusp import lorentzian_edges as le

    with mp.workdps(15):
        barnes = complex(le.regulated_double_pole_complex(-0.3, 0.3, 2.5, mp.pi))
    rotated = le.parametric_double_pole_rotated(-0.3, 0.3, 2.5, np.pi, 0.5)
    assert abs(rotated - barnes) < 1e-6 * abs(barnes)


def test_extended_field_is_radial_in_the_extra_dimensions():
    """V4 depends on e only through |e| and Ve is a multiple of e: the
    D-dimensional integral reduces to one over x0 and m = |e|."""
    C = fc.tuned_contour("--", 1.1)
    rng = np.random.default_rng(7)
    largest_sigma = 0.0
    for _ in range(300):
        x0 = C.centre + rng.normal(size=4) * 10 ** rng.uniform(-1, 1.3)
        e = rng.normal(size=3) * 10 ** rng.uniform(-3, 1)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        v4, ve = fc.field_extended(C, x0, e)
        v4r, ver = fc.field_extended(C, x0, q @ e)
        np.testing.assert_allclose(v4r, v4, rtol=0, atol=1e-12 * (1 + np.abs(v4).max()))
        np.testing.assert_allclose(ver, q @ ve, rtol=0, atol=1e-12 * (1 + np.abs(ve).max()))
        sigma = ve @ e / (e @ e)
        np.testing.assert_allclose(ve, sigma * e, rtol=0, atol=1e-14)
        largest_sigma = max(largest_sigma, abs(sigma))
    assert largest_sigma < 0.01                         # m(1 + i sigma) stays in Re > 0
