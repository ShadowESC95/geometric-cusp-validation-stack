"""Checks for docs/GATE_6_BOUNDS.md: all cusps, all regions, bounds, regulator."""
import mpmath as mp
import numpy as np
import pytest

from geocusp import lorentzian_edges as le

ONE = lambda L: 1
POWER = lambda L: mp.exp(mp.mpf("0.3") * L)          # not symmetric under z -> 1/z
LOG2 = lambda L: L**2 + mp.pi**2
REGIONS = {"--": (-1.0, -1.7), "+-": (1.3, -1.0), "-+": (-1.0, 0.8)}


@pytest.fixture(autouse=True)
def precision():
    with mp.workdps(18):
        yield


# ---- item 3: every cusp, every real kinematic region ----------------------
@pytest.mark.parametrize("region", list(REGIONS))
@pytest.mark.parametrize("cusp", [1, 2, 3, 4])
def test_local_formulas_are_exact_at_every_cusp(region, cusp):
    rng = np.random.default_rng(10 * cusp + len(region))
    for _ in range(4):
        alpha, beta = rng.normal(size=2) * 0.3
        rho, theta = abs(rng.normal()) * 0.2, rng.uniform(0, 2 * np.pi)
        exact, local, info = le.exact_invariants(cusp, alpha, beta, rho, theta, region)
        for j in exact:
            assert exact[j] == pytest.approx(local[j], abs=1e-12)
        assert info["Dp"] * info["u"] / info["D"] > 0        # |p_perp|^2 of the far edge


def test_exactly_one_mandelstam_invariant_is_positive():
    for region in REGIONS:
        _, s, t = le.light_like_square(region)
        assert sum(x > 0 for x in (s, t, -s - t)) == 1


@pytest.mark.parametrize("region", list(REGIONS))
@pytest.mark.parametrize("cusp", [1, 2, 3, 4])
def test_lorentzian_phases_never_exceed_two_pi(region, cusp):
    c = le.Cusp(cusp, *REGIONS[region])
    visited = set()
    for a in (1, -1):
        for b in (1, -1):
            visited |= set(c.sheet_sequence(a * 1e-3 * abs(c.D), b * 1e-7 * abs(c.D)))
            visited |= set(c.sheet_sequence(a * 1e-7 * abs(c.D), b * 1e-3 * abs(c.D)))
    assert max(abs(k) for k in visited) <= 2


@pytest.mark.parametrize("region", list(REGIONS))
@pytest.mark.parametrize("cusp", [1, 2, 3, 4])
@pytest.mark.parametrize("edge", ["next", "prev"])
def test_edge_statement_at_every_cusp_and_region(region, cusp, edge):
    c = le.Cusp(cusp, *REGIONS[region])
    mean = complex(le.circle_mean(POWER))
    on = complex(c.edge_jump(POWER, edge, 0.3))
    unit = complex(c.edge_jump(ONE, edge, 0.3))
    assert unit == pytest.approx(-2j * np.pi, abs=1e-8)
    assert on == pytest.approx(-2j * np.pi * mean, abs=1e-8)
    for outside in (-0.4, 1.5):
        assert abs(complex(c.edge_jump(POWER, edge, outside))) < 1e-8


# ---- item 3: the two-loop weight on its continued sheets ------------------
@pytest.mark.parametrize("L", [(0.3, -4.0), (-1.2, -5.5), (2.0, -3.5), (0.05, -3.3), (-0.02, -3.2)])
def test_two_loop_weight_is_regular_at_z_minus_one(L):
    """The continuations across -1 < z < 0 and across z < -1 agree, so there
    is no monodromy around z = -1 on the physical sheet."""
    with mp.workdps(25):
        point = mp.mpc(*L)
        inner, outer = le.f2_cover(point, "inner"), le.f2_cover(point, "outer")
        assert abs(inner - outer) < mp.mpf("1e-18") * abs(inner)


def test_two_loop_weight_is_continuous_across_the_lips_and_real_analytic():
    for ell in (-0.7, 0.9):
        up = le.f2_cover(mp.mpc(ell, -mp.pi + mp.mpf("1e-9")))
        down = le.f2_cover(mp.mpc(ell, -mp.pi - mp.mpf("1e-9")))
        assert abs(up - down) < 1e-5
    L = mp.mpc("0.4", "4.4")
    assert abs(le.f2_cover(L) - mp.conj(le.f2_cover(mp.conj(L)))) < 1e-14
    value = le.f2_cover(mp.mpc(0, mp.pi) - mp.mpf("1e-7"))
    assert abs(value - 11 * mp.pi**4 / 45) < 1e-4           # F2(-1) = 11 pi^4 / 45


def test_three_loop_coefficient_from_the_lorentzian_edge_chain():
    """With the two-loop weight, the edge chain gives -11 pi^4/15 on the
    segment, hence Gamma^(3) = 44 pi^4/45, and nothing off the segment."""
    with mp.workdps(15):
        c = le.Cusp(2, -1.0, -1.7)
        ratio = mp.mpf("1e-9")
        on = complex(c.edge_jump(le.f2_cover, "next", 0.4, ratio) / (-2j * mp.pi))
        off = complex(c.edge_jump(le.f2_cover, "next", -0.4, ratio) / (-2j * mp.pi))
    assert on == pytest.approx(-11 * np.pi**4 / 15, abs=1e-5)
    assert abs(off) < 1e-5
    assert -4 * on.real / 3 == pytest.approx(44 * np.pi**4 / 45, abs=1e-5)


# ---- item 1: bounds -------------------------------------------------------
def test_lemma_1_error_is_bounded_by_r_times_a_power_of_log():
    """|J(r) - J(0)| <= C r (1 + log(1/r))^(k+1), uniformly along the edge."""
    c = le.Cusp(2, -1.0, -1.7)
    limit = -2j * mp.pi * le.circle_mean(LOG2)
    for position in (0.05, 0.4, 0.8):
        scaled = []
        for r in (1e-3, 1e-5, 1e-7):
            err = abs(c.edge_jump(LOG2, "next", position, mp.mpf(r)) - limit)
            scaled.append(float(err / (r * (1 + mp.log(1 / r)) ** 3)))
        assert max(scaled) < 5.0
        assert scaled[2] < 2 * scaled[0]                 # no growth relative to the bound


def test_lemma_2_exact_fourth_propagator_changes_the_jump_by_c_log_c():
    """Replacing x04^2 by its edge value costs O(c log(1/c)) in the collinear
    invariant c, including the term linear in the transverse direction, which
    averages out over the angle."""
    c = le.Cusp(2, -1.0, -1.7)
    scaled = []
    for small in (1e-4, 1e-6):
        exact = le.exact_edge_jump(c, POWER, 0.4, small, 0.05, True, n_theta=12)
        model = le.exact_edge_jump(c, POWER, 0.4, small, 0.05, False)
        scaled.append(float(abs(exact - model) / (small * mp.log(1 / small))))
    assert 0.3 < scaled[1] / scaled[0] < 3.0
    mean = complex(le.circle_mean(POWER))
    assert complex(exact) == pytest.approx(-2j * np.pi * mean, abs=2e-3)


def test_hard_cutoff_is_not_an_admissible_regulator_for_growing_weights():
    """Where the two collinear invariants are comparable, the quadrant sum
    departs from its collinear limit by an amount that grows with the cutoff
    when the weight grows at z -> infinity. The edge statement is unaffected;
    the double pole must be taken with an analytic regulator."""
    c = le.Cusp(2, -1.0, -1.0)
    limit = 2j * mp.pi * le.circle_mean(LOG2)
    gaps = [abs(c.quadrant_sum(LOG2, q, q) - limit) for q in (mp.mpf("1e-8"), mp.mpf("1e-16"))]
    assert gaps[1] > 3 * gaps[0] > 30
    flat = [abs(c.quadrant_sum(ONE, q, q) - 2j * mp.pi) for q in (mp.mpf("1e-8"), mp.mpf("1e-16"))]
    # the constant weight departs by a fixed amount (exactly 2), which does not grow
    assert flat[0] == pytest.approx(2.0, abs=1e-5) and flat[1] == pytest.approx(2.0, abs=1e-5)


# ---- item 2: dimensional regularisation -----------------------------------
@pytest.mark.parametrize("p,eps", [(-0.2, 0.5), (-0.4, 0.3)])
def test_dimreg_density_matches_the_closed_form(p, eps):
    c = le.Cusp(2, -1.0, -1.7)
    numeric = complex(le.dimreg_quadrant_sum(c, p, eps, 1e-9, 3e-9))
    closed = complex(le.dimreg_density_closed_form(p, eps))
    assert numeric == pytest.approx(closed, rel=2e-4)


@pytest.mark.parametrize("p", [0.3, -0.6, 0.85])
def test_dimreg_density_tends_to_the_unregulated_value(p):
    target = 2j * mp.pi * mp.sin(mp.pi * p) / (mp.pi * p)
    gaps = [abs(le.dimreg_density_closed_form(p, mp.mpf(10) ** -k) - target) for k in (2, 4, 6)]
    assert gaps[0] > 50 * gaps[1] > 2500 * gaps[2]
    assert gaps[2] < 1e-4


# ---- analytic regularity of the two-loop weight at z = -1 -----------------
def test_two_loop_weight_has_no_logarithm_at_z_minus_one():
    """After the inversion relations the log^3 and log^4 terms vanish
    identically and the log and log^2 terms are multiples of Landen's
    combination, which is zero. This is exact, not numerical."""
    import sympy as sp

    coefficients, landen, _ = le.f2_expansion_at_minus_one()
    assert coefficients[4] == 0 and coefficients[3] == 0
    for k in (1, 2):
        quotient, remainder = sp.div(coefficients[k], landen, *landen.free_symbols)
        assert remainder == 0 and quotient != 0
    assert le.landen_derivative() == 0                   # and the combination is 0 at e = 0


def test_regular_part_reproduces_the_weight_on_both_sides_of_z_minus_one():
    with mp.workdps(25):
        assert abs(le.f2_regular_part(0) - 11 * mp.pi**4 / 45) < mp.mpf("1e-20")
        # upper lip: Im(1+z) > 0 on the principal sheet
        for e in (mp.mpc("0.05", "0.02"), mp.mpc("-0.08", "0.03"), mp.mpc("0.001", "0.2")):
            L = mp.log(-1 + e)
            assert abs(le.f2_regular_part(e) - le.f2_cover(L)) < mp.mpf("1e-18")
        # continued through z = -1 to Im(1+z) < 0: the sheet with Im L > pi
        for e in (mp.mpc("0.05", "-0.02"), mp.mpc("-0.08", "-0.03")):
            L = mp.log(-1 + e) + 2j * mp.pi
            assert abs(le.f2_regular_part(e) - le.f2_cover(L)) < mp.mpf("1e-18")


# ---- uniform bound on the regulated integral (Lemma 4) --------------------
def test_regulated_integral_reduces_to_the_existing_real_case():
    from geocusp import solver

    for p in (0.3, -0.6):
        mine = complex(le.regulated_double_pole(p, 0.02))
        assert mine.real == pytest.approx(solver.double_pole(p, 0.02), abs=2e-6)
        assert abs(mine.imag) < 1e-12


@pytest.mark.parametrize("sigma", [0.3, -0.5])
def test_regulated_integral_grows_no_faster_than_exp_pi_im_p(sigma):
    """|eps^2 I_eps(sigma + i tau)| (1 + tau) exp(-pi tau) stays bounded,
    uniformly as the regulator is removed."""
    for eta in (0.1, 0.02, 0.004):
        scaled = [float(abs(le.regulated_double_pole(mp.mpc(sigma, tau), eta))
                        * (1 + tau) * mp.exp(-mp.pi * tau)) for tau in (0, 2, 4, 8)]
        assert max(scaled) < 4.0 and scaled[-1] < 1.0


def test_regulated_integral_converges_for_complex_exponent():
    p = mp.mpc("0.3", "2.0")
    target = 4 * mp.sin(mp.pi * p) / (mp.pi * p)
    gaps = [abs(le.regulated_double_pole(p, eta) - target) for eta in (0.02, 0.002)]
    assert gaps[1] < gaps[0] / 5 and gaps[1] < 0.5


# ---- pinch surfaces --------------------------------------------------------
def test_only_the_edges_can_pinch_the_contour():
    import sympy as sp

    result = le.landau_pinch_analysis()
    a, c, d, s, t = result["symbols"]
    assert result["adjacent_pair_on_edge"] == [0, 0]          # the whole edge x2 -> x3
    first, second = result["diagonal_pair"]                   # the diagonal x1 -> x3
    assert sp.simplify(first - s * (a - 1) ** 2) == 0 and sp.simplify(second - s * a**2) == 0
    assert result["triple_123"] == [{a: 0, c: 0}]             # only the cusp x2
    assert result["quadruple"] == []


# ---- the regulated integral: Feynman parameters against Mellin-Barnes -----
@pytest.mark.parametrize("p,eta,ratio", [(0.2, 0.4, 1.0), (-0.3, 0.3, 2.5), (0.0, 0.25, 1.0)])
def test_feynman_parameter_form_equals_the_mellin_barnes_form(p, eta, ratio):
    """Two representations of eps^2 s t times the regulated integral of z^p:
    the Feynman-parameter integral in D = 4 + 2 eta, and the one-fold
    Mellin-Barnes integral used elsewhere. They agree, normalisation included."""
    parametric = le.parametric_double_pole(p, eta, ratio)
    barnes = complex(le.regulated_double_pole(p, eta, ratio))
    assert abs(barnes.imag) < 1e-12
    assert parametric == pytest.approx(barnes.real, rel=1e-7)


def test_extra_dimensions_integrate_to_the_dimensional_power():
    """int_0^oo d mu mu^(eta-1)/Gamma(eta) * Gamma(2)/(A + mu)^2 = Gamma(2-eta) A^(eta-2):
    integrating out the 2 eta extra dimensions is a shift x0i^2 -> x0i^2 - mu."""
    for eta, A in ((0.3, 0.7), (0.05, 2.5)):
        eta, A = mp.mpf(eta), mp.mpf(A)
        # mu = v^(1/eta) removes the endpoint singularity: mu^(eta-1) d mu = dv / eta
        lhs = mp.quad(lambda v: 1 / (A + v ** (1 / eta)) ** 2,
                      [0, A**eta / 2, A**eta, 2 * A**eta, mp.inf]) / (eta * mp.gamma(eta))
        assert abs(lhs - mp.gamma(2 - eta) * A ** (eta - 2)) < mp.mpf("1e-10")


# ---- repairs after review -------------------------------------------------
@pytest.mark.parametrize("region", list(REGIONS))
@pytest.mark.parametrize("cusp", [1, 2, 3, 4])
def test_far_invariant_on_the_plane_is_exact(region, cusp):
    """x0,i+2^2 at y_perp = 0 is D'(1 + alpha - beta) + alpha beta D exactly,
    not D'(1 - beta)(1 + alpha)."""
    rng = np.random.default_rng(3 * cusp + len(region))
    for _ in range(5):
        alpha, beta = rng.normal(size=2) * 0.4
        exact, _, info = le.exact_invariants(cusp, alpha, beta, 0.0, 0.0, region)
        c = le.Cusp(cusp, info["D"], info["Dp"]) if cusp % 2 == 0 else le.Cusp(cusp, info["Dp"], info["D"])
        far = (cusp + 1) % 4 + 1
        assert c.far_edge_value(beta * info["D"], -alpha * info["D"]) == pytest.approx(exact[far], abs=1e-12)


@pytest.mark.parametrize("sigma", [0.3, -0.5])
def test_uniform_bound_at_a_timelike_diagonal(sigma):
    """Lemma 4 at arg(S/T) = pi: |N_eta| (1+tau) exp(-pi tau) stays bounded,
    uniformly in eta."""
    mp.mp.dps = 15
    values = [abs(le.regulated_double_pole_complex(mp.mpc(sigma, tau), eta, 1.0, mp.pi))
              * (1 + tau) * mp.exp(-mp.pi * tau)
              for eta in (0.1, 0.004) for tau in (0, 2, 4, 8)]
    assert 0.4 < min(values) and max(values) < 5


def test_two_loop_weight_is_polynomially_bounded_on_the_cover():
    """|F2(e^L)| <= C (1 + |L|)^4 for |Im L| <= 2.9 pi (Lemma 3'). The ratio
    is 128.7 at L = 0."""
    mp.mp.dps = 15
    rng = np.random.default_rng(0)
    points = [mp.mpc(0, 0)] + [mp.mpc(rng.uniform(-40, 40), rng.uniform(-2.9, 2.9) * np.pi)
                               for _ in range(300)]
    points += [mp.mpc(x, y * np.pi) for x in (-1000, -3, 3, 1000) for y in (-2.9, -1.5, 0.5, 2.9)]
    worst = max(float(abs(le.f2_cover(L)) / (1 + abs(L)) ** 4) for L in points)
    assert 128 < worst < 135


# ---- meromorphy in eta (Theorem M) -------------------------------------------
def test_continued_box_agrees_with_the_mellin_barnes_form():
    import mpmath as mp
    from geocusp import lorentzian_edges as le
    for p, eta, ratio in ((0.2, 0.1, 2.5), (0.1 + 0.7j, 0.2, 0.4)):
        a, b = le.regulated_double_pole(p, eta, ratio), le.box_continued(p, eta, ratio, 0)
        assert abs(a - b) < 1e-12 * abs(a)
    a = le.regulated_double_pole_complex(0.2, 0.1, 1.7, mp.pi)
    assert abs(a - le.box_continued(0.2, 0.1, 1.7, mp.pi)) < 1e-12 * abs(a)


def test_continued_box_is_analytic_through_eta_zero():
    import mpmath as mp
    from geocusp import lorentzian_edges as le
    p = 0.2
    limit = 4 * mp.sin(mp.pi * p) / (mp.pi * p)
    above, below = le.box_continued(p, 1e-4, 2.5), le.box_continued(p, -1e-4, 2.5)
    assert abs(above - limit) < 1e-3 and abs(below - limit) < 1e-3
    assert abs((above + below) / 2 - limit) < 1e-6          # no pole, no jump
    assert abs(le.box_continued(p, 0, 2.5) - limit) < 1e-12
    for eta in (-0.2, 0.1j, -0.3 + 0.2j):
        assert mp.isfinite(le.box_continued(p, eta, 2.5))
    with pytest.raises(ValueError):
        le.box_continued(0.2, 0.4, 2.5)


def test_continued_box_grows_no_faster_than_the_bound():
    import mpmath as mp
    from geocusp import lorentzian_edges as le
    for eta in (-0.2, 0.15j):
        scaled = [abs(le.box_continued(0.1 + 1j * t, eta, 2.5)) * (1 + t) * mp.exp(-mp.pi * t) for t in (2, 4, 8)]
        assert max(scaled) < 5 and min(scaled) > 0.1


@pytest.mark.parametrize("p,eta,timelike", [
    (0.2, 0.25, False), (0.2, 0.25, True), (0.2, 0.03, False), (-0.3, 1e-3, False),
    (-0.3, 1e-3, True), (0.2, 1e-6, False), (0.2, 1e-6, True)])
def test_parameter_quadratures_are_valid_down_to_eta_zero(p, eta, timelike):
    """The two direct parameter quadratures against the continued form, which is
    analytic through eta = 0. With fixed ranges they used to fail silently
    below eta ~ 0.01 (spacelike) and ~ 0.1 (timelike): the double pole comes
    from tails of length 1/eta, now integrated in closed form."""
    if timelike:
        value = le.parametric_double_pole_rotated(p, eta, 2.5, np.pi, 0.5)
        continued = complex(le.box_continued(p, eta, 2.5, np.pi))
    else:
        value = le.parametric_double_pole(p, eta, 2.5)
        continued = complex(le.box_continued(p, eta, 2.5, 0))
    assert abs(value - continued) < 1e-8 * abs(continued)
    if eta <= 1e-6:
        assert abs(value - 4 * np.sin(p * np.pi) / (p * np.pi)) < 1e-4


@pytest.mark.parametrize("p,eta,modulus,timelike", [
    (-0.96, 0.05, 2.5, False),      # 1 + p small: a range of length 1/(1+p), closed form
    (0.9, 0.05, 2.5, False),        # 1 - p - eta small: a range of length 1/(1-p-eta)
    (0.9, 0.05, 2.5, True),
    (0.0, 0.9, 1.0, False),         # eta outside the range of the continued form
])
def test_parameter_quadratures_at_the_ends_of_their_range(p, eta, modulus, timelike):
    """Against the Mellin-Barnes form (valid here: 1 + p > eta/2). A first rewrite
    returned a wrong value near p = -1 and did not finish near p + eta = 1."""
    theta = np.pi if timelike else 0.0
    barnes = complex(le.regulated_double_pole_complex(p, eta, modulus, theta))
    if timelike:
        value = le.parametric_double_pole_rotated(p, eta, modulus, np.pi, 0.5)
    else:
        value = le.parametric_double_pole(p, eta, modulus)
    assert abs(value - barnes) < 5e-9 * abs(barnes)


def test_parameter_quadratures_refuse_what_they_cannot_do():
    timelike = 2.5 * np.exp(1j * np.pi)
    for gamma in (-0.5, 0.0, 0.05, 3.2):                 # wrong way, on the zero, too close, too far
        with pytest.raises(ValueError):
            le._parametric_box(0.2, 0.1, timelike, gamma)
    for eta in (1e-160, 1e-13, 0.0, 0.9):                # overflow, below the floor, zero, p + eta >= 1
        with pytest.raises(ValueError):
            le._parametric_box(0.2, eta, 2.5)
    with pytest.raises(ValueError):
        le._parametric_box(0.2, 0.1, 2.5 * np.exp(-0.3j))  # arg S < 0
    # a real positive S needs no rotation and ignores one
    assert abs(le._parametric_box(0.2, 0.1, 2.5, 0.5) - le._parametric_box(0.2, 0.1, 2.5, 0.0)) < 1e-12


def test_parameter_quadratures_do_not_depend_on_the_rotation_or_the_window():
    a = le.parametric_double_pole_rotated(0.2, 0.01, 2.5, np.pi, 0.5)
    b = le.parametric_double_pole_rotated(0.2, 0.01, 2.5, np.pi, 0.9)
    assert abs(a - b) < 1e-8
    c = le._parametric_box(0.2, 0.01, 2.5, 0.0)
    d = le._parametric_box(0.2, 0.01, 2.5, 0.0, window=60.0)
    assert abs(c - d) < 1e-8
    for eta in (0.0, -0.1):
        with pytest.raises(ValueError):
            le.parametric_double_pole(0.2, eta, 2.5)
