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
