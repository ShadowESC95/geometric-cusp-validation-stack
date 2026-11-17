"""Checks for docs/CUSP_LOCAL_MODEL.md."""
import mpmath as mp
import pytest

from geocusp import cusp_local as cl


@pytest.fixture(autouse=True)
def precision():
    mp.mp.dps = 20


def test_measure_identity_is_exact():
    ratio, residual = cl.measure_identity()
    assert ratio in (1, -1) and residual == 0


@pytest.mark.parametrize("p", [0.3, -0.5, mp.mpc(0.2, 1.5)])
def test_triangle_double_pole_is_the_circle_mean(p):
    limit = mp.sin(mp.pi * p) / (mp.pi * p)
    errors = [abs(cl.triangle_closed_form(p, eta) - limit) for eta in (1e-2, 1e-3, 1e-4)]
    assert errors[2] < 1e-3 * abs(limit) * 5
    assert errors[1] / errors[2] == pytest.approx(10, rel=0.1)      # linear in eta


def test_triangle_is_symmetric_in_the_limit_under_p_to_minus_p():
    assert abs(cl.triangle_closed_form(0.4, 1e-6) - cl.triangle_closed_form(-0.4, 1e-6)) < 1e-4


def test_dirichlet_integral_of_the_triangle():
    numeric, closed = cl.dirichlet_check(0.3, 0.25)
    assert numeric == pytest.approx(closed, rel=1e-10)


@pytest.mark.parametrize("p", [0.3, -0.3, mp.mpc(0.2, 1.0)])
def test_vertex_integral_identity(p):
    numeric, closed, main = cl.vertex_main_term(p, 0.25, 0.1)
    assert abs(numeric - closed) < 1e-9 * abs(closed)
    assert abs(main) > 0


@pytest.mark.parametrize("p,ratio", [(0.3, 1), (-0.5, 2.5), (mp.mpc(0.2, 1.5), 1)])
def test_remainder_is_a_single_pole(p, ratio):
    """(N_eta - 4 sin(pi p)/(pi p)) / eta stays bounded as eta -> 0."""
    values = [abs(cl.remainder_ratio(p, eta, ratio)) for eta in (0.2, 0.1, 0.05, 0.025)]
    assert max(values) < 60
    assert abs(values[3] - values[2]) < abs(values[1] - values[0])     # settling


def test_four_cusps_add_to_the_box():
    """2 x triangle(p) + 2 x triangle(-p) and N_eta have the same limit."""
    p, eta = mp.mpf("0.3"), mp.mpf("1e-4")
    cusps = 2 * cl.triangle_closed_form(p, eta) + 2 * cl.triangle_closed_form(-p, eta)
    assert cusps == pytest.approx(float(4 * mp.sin(mp.pi * p) / (mp.pi * p)), rel=2e-3)


# ---- timelike diagonal and the exact formula for a general weight ----------
@pytest.mark.parametrize("p", [0.3, mp.mpc(0.2, 1.5), mp.mpc(0.2, -1.5)])
def test_timelike_remainder_is_a_single_pole(p):
    values = [abs(cl.timelike_remainder_ratio(p, eta)) for eta in (0.2, 0.1, 0.05)]
    assert max(values) < 120
    assert abs(values[2] - values[1]) < abs(values[1] - values[0]) + 1


@pytest.mark.parametrize("p,eta", [(0.3, 0.2), (0.2 + 1.0j, 0.15), (-0.4, 0.05)])
def test_regulated_circle_formula_reproduces_the_triangle(p, eta):
    import numpy as np
    value = cl.regulated_circle_formula(lambda L: np.exp(p * L), eta)
    assert value == pytest.approx(complex(cl.triangle_closed_form(p, eta)), rel=2e-6)


def test_regulated_circle_formula_for_a_logarithmic_weight():
    """F = log^2 z + pi^2: the formula agrees with the second p-derivative of
    the closed form, and tends to M[F] = 2 pi^2/3 linearly in eta."""
    import numpy as np
    mean = 2 * np.pi**2 / 3
    ratios = []
    for eta in (0.2, 0.05, 0.0125):
        value = cl.regulated_circle_formula(lambda L: L**2 + np.pi**2, eta)
        reference = complex(mp.diff(lambda q: cl.triangle_closed_form(q, eta), 0, 2)
                            + mp.pi**2 * cl.triangle_closed_form(0, eta))
        assert value == pytest.approx(reference, rel=1e-6)
        ratios.append(abs(value - mean) / eta)
    assert max(ratios) < 12 and min(ratios) > 9


def test_smearing_density_is_a_probability_density_of_width_eta():
    """Mass 1 (with x = y^(1/eta) to remove the endpoint singularity) and
    mean psi(1) - psi(1 - eta), which is about eta pi^2/6."""
    for eta in (0.2, 0.05):
        e = mp.mpf(eta)
        near = mp.quad(lambda y: cl.smearing_density(y ** (1 / e), e) * y ** (1 / e - 1) / e, [mp.mpf("1e-12"), 1])
        far = mp.quad(lambda x: cl.smearing_density(x, e), [1, 10, mp.inf])
        assert near + far == pytest.approx(1, abs=1e-8)
        mean = mp.quad(lambda x: x * cl.smearing_density(x, e), [0, 1, 10, mp.inf])
        assert mean == pytest.approx(float(mp.digamma(1) - mp.digamma(1 - e)), rel=1e-3)
        assert 0.5 * eta < mean < 3 * eta


def test_circle_formula_prefactor_tends_to_one():
    assert cl.circle_formula_prefactor(1e-6) == pytest.approx(1, abs=1e-5)


# ---- derivation from the Lorentzian domain (section 8) ---------------------
import numpy as np

RHO = 0.3 * np.exp(0.9j)        # a phase on rho is the Feynman prescription


@pytest.mark.parametrize("beta", [-0.4, 1.3, 1.6])
def test_light_cone_integral_vanishes_off_the_edge_segment(beta):
    for f in (lambda L: np.exp((-0.35 + 0.4j) * L), lambda L: (L**2 + np.pi**2) * np.exp(-0.2 * L)):
        assert abs(cl.light_cone_integral(f, beta, RHO)) < 1e-9


@pytest.mark.parametrize("beta", [0.3, 0.7])
def test_light_cone_integral_is_a_discontinuity_on_the_segment(beta):
    f = lambda L: np.exp((-0.35 + 0.4j) * L)
    direct = cl.light_cone_integral(f, beta, RHO)
    closed = cl.discontinuity_integral(f, beta, RHO)
    assert abs(direct) > 1 and direct == pytest.approx(closed, rel=1e-5)


def test_discontinuity_form_for_a_logarithmic_weight():
    f = lambda L: (L**2 + np.pi**2) * np.exp(-0.2 * L)
    assert cl.light_cone_integral(f, 0.4, RHO) == pytest.approx(
        cl.discontinuity_integral(f, 0.4, RHO, reach=150), rel=1e-7)


@pytest.mark.parametrize("p,eta", [(0.3, 0.2), (0.2 + 1.0j, 0.15), (-0.4, 0.05)])
def test_lorentzian_product_formula_gives_the_closed_form(p, eta):
    value = cl.lorentzian_product_formula(lambda L: np.exp(p * L), eta)
    assert value == pytest.approx(complex(cl.triangle_closed_form(p, eta)), rel=5e-4)


def test_lorentzian_and_circle_formulas_agree_for_a_logarithmic_weight():
    f = lambda L: L**2 + np.pi**2
    for eta in (0.2, 0.05):
        assert cl.lorentzian_product_formula(f, eta) == pytest.approx(
            cl.regulated_circle_formula(f, eta), rel=5e-4)


def test_three_beta_integrals_give_the_closed_form():
    """B(-p, 1+eta) B(eta, 1+p-eta) B(eta, 1+eta-p) times -2i sin(pi p)."""
    p, eta = mp.mpc(-0.35, 0.4), mp.mpf("0.2")
    m = mp.beta(-p, 1 + eta) * mp.beta(eta, 1 + p - eta) * mp.beta(eta, 1 + eta - p)
    lorentzian = -2j * mp.sin(mp.pi * p) * m
    closed = 2j * mp.pi * mp.gamma(1 + eta) / eta**2 * cl.triangle_closed_form(p, eta)
    assert abs(lorentzian - closed) < 1e-12 * abs(closed)


# ---- the full four-point integrand (section 9) -----------------------------
FULL = dict(rho=0.35, kappa=0.3, T=1.7, eps=2e-2)
WEIGHTS = [lambda L: np.exp((-0.35 + 0.4j) * L), lambda L: L**2 + np.pi**2, lambda L: 1 + 0 * L]


@pytest.mark.parametrize("beta", [-2.5, -1.9, 1.2, 2.0])
def test_full_integrand_vanishes_outside_the_null_slab(beta):
    for f in WEIGHTS:
        assert abs(cl.full_light_cone_integral(f, beta, **FULL)) < 1e-9


@pytest.mark.parametrize("beta,apart", [(-2.5, None), (-1.0, "u4"), (-0.3, "u4"), (0.4, "u3"), (0.8, "u3"), (1.2, None)])
def test_which_zero_lies_apart(beta, apart):
    zeros = cl.full_alpha_zeros(beta, FULL["rho"], FULL["kappa"], FULL["T"], FULL["eps"])
    upper = [k for k, v in zeros.items() if v.imag > 0]
    lower = [k for k, v in zeros.items() if v.imag < 0]
    if apart is None:
        assert not upper or not lower
    else:
        assert [apart] in (upper, lower)


@pytest.mark.parametrize("beta", [-1.0, -0.3, 0.4, 0.8])
def test_full_integrand_is_a_keyhole_inside_the_slab(beta):
    for f in WEIGHTS:
        direct = cl.full_light_cone_integral(f, beta, **FULL)
        keyhole = cl.full_keyhole_integral(f, beta, **FULL)
        assert abs(direct) > 1e-3 and direct == pytest.approx(keyhole, rel=1e-6)


# ---- the far end of an edge piece (section 10) -----------------------------
@pytest.mark.parametrize("v,p,eta", [(0.3, -0.3, 0.5), (0.02, -0.6, 0.3)])
def test_far_end_transverse_integral(v, p, eta):
    assert cl.far_end_transverse(v, p, eta) == pytest.approx(
        cl.far_end_transverse_feynman(v, p, eta), rel=1e-8)


@pytest.mark.parametrize("p,eta", [(-0.3, 0.5), (-0.6, 0.3), (-0.2, 0.1)])
def test_far_end_integral_has_the_closed_form(p, eta):
    assert cl.far_end_numerical(p, eta) == pytest.approx(float(cl.far_end_closed_form(p, eta)), rel=1e-7)


@pytest.mark.parametrize("p", [0.3, -0.5, mp.mpc(0.2, 1.0)])
def test_far_end_carries_one_unit_like_the_near_end(p):
    limit = mp.sin(mp.pi * p) / (mp.pi * p)
    far = cl.far_end_coefficient(p, mp.mpf("1e-5"), T=1.7)
    near = cl.triangle_closed_form(p, mp.mpf("1e-5"))
    assert abs(far - limit) < 1e-3 * abs(limit) and abs(near - limit) < 1e-3 * abs(limit)


def test_far_end_coefficient_is_the_scaled_closed_form():
    """(1/(2 eta)) 2 i sin(pi p) J, in the units of the near end."""
    p, eta, T = mp.mpf("-0.3"), mp.mpf("0.2"), mp.mpf("1.7")
    scaled = (T / (1 + T)) ** eta / (2 * eta) * 2j * mp.sin(mp.pi * p) * cl.far_end_closed_form(p, eta)
    near_units = -2j * mp.pi ** (2 + eta) / (mp.gamma(1 + eta) * eta**2) * cl.far_end_coefficient(p, eta, T)
    assert abs(scaled - near_units) < 1e-12 * abs(scaled)
