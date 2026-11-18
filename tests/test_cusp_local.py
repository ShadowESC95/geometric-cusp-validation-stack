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


# ---- the edge piece: reduction and remainder (section 11) ------------------
@pytest.mark.parametrize("p,eta,T", [(0.3, 0.25, 1.0), (0.15, 0.1, 1.0), (0.3, 0.2, 3.0), (-0.25, 0.1, 0.4)])
def test_two_edge_pieces_add_up_to_the_box(p, eta, T):
    """A(p; T) + T^eta A(-p; 1/T) equals the box computed by Mellin-Barnes."""
    from geocusp.lorentzian_edges import regulated_double_pole
    total = cl.box_from_two_edges(p, eta, T)
    box = -2j * np.pi ** (2 + eta) * T**eta * complex(regulated_double_pole(p, eta, 1 / T)) / eta**2
    assert total == pytest.approx(box, rel=1e-8)


@pytest.mark.parametrize("p,T", [(-0.3, 1.0), (0.3, 1.0), (-0.3, 3.0), (0.5, 0.4)])
def test_edge_piece_minus_both_end_terms_is_a_single_pole(p, T):
    values = []
    for eta in (0.1, 0.05, 0.025):
        r = cl.edge_piece(p, eta, T) - cl.edge_near_term(p, eta) - cl.edge_far_term(p, eta, T)
        values.append((eta * r / (-2j * np.pi**2)).real)
    assert max(abs(v) for v in values) < 1.5
    assert abs(values[2] - values[1]) < abs(values[1] - values[0])        # settling


def test_far_end_term_from_the_cube_agrees_with_section_ten():
    """(T+1)^(p-1) tau^(eta-1+p) h0/(2 eta) in the units of A is the far-end closed form."""
    from scipy.special import gamma as G
    p, eta, T = -0.3, 0.15, 1.7
    tau = T / (T + 1)
    h0 = G(1 + p) * G(1 - p - eta) / G(2 - eta) * G(eta) ** 2 / G(2 * eta)
    j_far = (T + 1) ** (p - 1) * tau ** (eta - 1 + p) * h0 / (2 * eta)
    a_far = -2j * np.pi ** (2 + eta) * G(2 - eta) * T ** (1 - p) / (G(1 - p) ** 2 * G(1 + p) ** 2) * j_far
    assert a_far == pytest.approx(cl.edge_far_term(p, eta, T), rel=1e-12)


def test_near_end_term_from_the_cube_has_the_same_double_pole():
    from scipy.special import gamma as G
    p, T = 0.3, 2.0
    for eta in (1e-3, 1e-4):
        j_near = T ** (p - 1) * G(1 + p) * G(1 - p) / eta**2
        a_near = -2j * np.pi ** (2 + eta) * G(2 - eta) * T ** (1 - p) / (G(1 - p) ** 2 * G(1 + p) ** 2) * j_near
        assert abs(a_near / cl.edge_near_term(p, eta) - 1) < 20 * eta


# ---- general weights (section 12) ------------------------------------------
def test_edge_piece_with_a_logarithmic_weight_has_double_pole_two_means():
    """F = log^2 z + pi^2: A[F] = A''(0) + pi^2 A(0); eta^2 A[F]/(-2 pi^2 i) -> 2 M[F]."""
    target = 2 * (2 * np.pi**2 / 3)
    slopes = []
    for eta in (0.1, 0.05):
        A = lambda q: cl.edge_piece(q, eta, 1.0)
        h = 0.08
        a0 = A(1e-9)
        d2 = (-(A(2 * h) + A(-2 * h)) + 16 * (A(h) + A(-h)) - 30 * a0) / (12 * h * h)
        value = (eta**2 * (d2 + np.pi**2 * a0) / (-2j * np.pi**2)).real
        slopes.append((value - target) / eta)
    assert 20 < slopes[0] < 30 and 20 < slopes[1] < 30          # O(eta) correction with a stable coefficient
    assert abs(slopes[1] - slopes[0]) < 3


def test_far_end_factor_tends_to_one_uniformly_on_a_line():
    """eta^2 Gamma(eta)^2/Gamma(1+2 eta) * Gamma(1-p-eta)/Gamma(1-p) - 1 = O(eta log|tau|)."""
    eta = mp.mpf("0.01")
    for tau in (0, 3, 30, 300):
        p = mp.mpc(0.1, tau)
        factor = eta**2 * mp.gamma(eta) ** 2 / mp.gamma(1 + 2 * eta) * mp.gamma(1 - p - eta) / mp.gamma(1 - p)
        assert abs(factor - 1) < 10 * eta * (1 + mp.log(2 + tau))


# ---- any region: support and keyholes (section 14) -------------------------
REGIONS14 = {"--": (-1, -1.7), "-+": (-1, 2.2), "+-": (1, -0.45)}
LOGW = lambda L: L**2 + np.pi**2


@pytest.mark.parametrize("region,beta", [("--", -2.5), ("--", 1.5), ("-+", -0.3), ("-+", 2.6), ("+-", -0.3), ("+-", 1.5)])
def test_alpha_integral_vanishes_outside_the_slab_in_every_region(region, beta):
    s, t = REGIONS14[region]
    assert abs(cl.general_light_cone_integral(LOGW, beta, 0.35, 0.3, s, t)) < 1e-9
    assert cl.general_keyhole_sum(LOGW, beta, 0.35, 0.3, s, t)[1] == []


@pytest.mark.parametrize("region,beta,apart", [
    ("--", -1.0, {"u4"}), ("--", 0.4, {"u3"}),
    ("-+", 0.4, {"u3", "u4"}), ("-+", 1.5, {"u4"}),
    ("+-", 0.2, {"w", "u1"}), ("+-", 0.8, {"u3"})])
def test_alpha_integral_is_a_sum_of_keyholes_in_every_region(region, beta, apart):
    s, t = REGIONS14[region]
    direct = cl.general_light_cone_integral(LOGW, beta, 0.35, 0.3, s, t)
    value, group, _ = cl.general_keyhole_sum(LOGW, beta, 0.35, 0.3, s, t)
    assert set(group) == apart
    assert abs(direct) > 1e-3 and direct == pytest.approx(value, rel=1e-5)


def test_far_invariant_is_positive_on_the_cut_with_a_timelike_diagonal():
    """s = -1, t = T > 1: x04^2 = [(T - beta) u3 + (T - 1)|y + y*|^2]/b on the cut."""
    T, beta, rho, kappa, u3 = 2.2, 0.4, 0.35, 0.3, 0.7
    b = 1 - beta
    _, coeff, const = cl.general_alpha_zeros(beta, rho, kappa, -1, T, 0.0)
    alpha = (u3 + rho) / b
    ystar = b * np.sqrt(T / (T - 1))
    formula = ((T - beta) * u3 + (T - 1) * (rho + 2 * kappa * ystar + ystar**2)) / b
    assert (coeff["u4"] * alpha + const["u4"]).real == pytest.approx(formula, rel=1e-12)
    assert formula > 0 and (coeff["w"] * alpha + const["w"]).real < 0 and (coeff["u1"] * alpha + const["u1"]).real < 0


def test_first_invariant_changes_sign_on_the_cut_when_s_is_positive():
    """s = +1: on the cut around the zero of x03^2, x01^2 takes both signs."""
    T, beta, rho, kappa = 0.45, 0.6, 0.05, 0.02
    _, coeff, const = cl.general_alpha_zeros(beta, rho, kappa, 1, -T, 0.0)
    values = []
    for u3 in (0.01, 2.0):
        alpha = (u3 + rho) / (coeff["u3"])                        # u3 = coeff * alpha - rho
        values.append((coeff["u1"] * alpha + const["u1"]).real)
    assert values[0] * values[1] < 0


@pytest.mark.parametrize("p,T", [(-0.3, 2.2), (0.3, 2.2)])
def test_timelike_edge_piece_minus_end_terms_is_a_single_pole(p, T):
    values = []
    for eta in (0.1, 0.05, 0.025):
        r = cl.edge_piece(p, eta, T, sign=-1) - cl.edge_near_term(p, eta) - cl.edge_far_term(p, eta, T, sign=-1)
        values.append((eta * r / (-2j * np.pi**2)).real)
    assert max(abs(v) for v in values) < 2
    assert abs(values[2] - values[1]) < abs(values[1] - values[0])


# ---- the edge jump as the residue of the collinear pole --------------------
@pytest.mark.parametrize("beta", [0.3, 0.7])
@pytest.mark.parametrize("p", [0.25, -0.4, 0.2 + 0.6j])
def test_collinear_residue_for_powers(beta, p):
    value = cl.collinear_residue(lambda L: np.exp(p * L), beta)
    reference = 2j * np.pi * np.sin(np.pi * p) / (np.pi * p) * (1 - beta) ** (-p) / beta
    assert value == pytest.approx(reference, rel=1e-7)


@pytest.mark.parametrize("beta", [0.3, 0.7])
def test_collinear_residue_for_a_logarithmic_weight(beta):
    value = cl.collinear_residue(lambda L: L**2 + np.pi**2, beta)
    reference = 2j * np.pi / beta * (2 * np.pi**2 / 3 + np.log(1 - beta) ** 2)
    assert value == pytest.approx(reference, rel=1e-7)


# ---- the edge piece as a one-dimensional integral -----------------------------
@pytest.mark.parametrize("p,eta,T", [(0.3, 0.2, 2.5), (-0.4, 0.1, 0.3)])
def test_one_dimensional_form_equals_the_cube_integral(p, eta, T):
    cube = complex(cl.edge_piece(p, eta, T))
    one = complex(cl.edge_piece_one_dim(p, eta, T / (1 + T)))
    assert abs(cube - one) < 1e-10 * abs(cube)


def test_one_dimensional_form_equals_the_piece_with_a_timelike_diagonal_that_works():
    p, eta, T = -0.3, 0.2, 2.2
    cube = complex(cl.edge_piece(p, eta, T, sign=-1))
    one = complex(cl.edge_piece_one_dim(p, eta, T / (T - 1)))
    assert abs(cube - one) < 1e-10 * abs(cube)


@pytest.mark.parametrize("p,eta,modulus", [(0.3, 0.2, 0.45), (0.3, 0.2, 2.2), (-0.25, 0.1, 0.6)])
def test_two_continued_pieces_add_up_to_the_box_with_a_timelike_diagonal(p, eta, modulus):
    import mpmath as mp
    from geocusp import lorentzian_edges as le
    first, second = cl.timelike_two_edges(p, eta, modulus, sigma=-1)
    box = (-2j * mp.pi ** (2 + mp.mpf(eta)) * mp.mpf(modulus) ** mp.mpf(eta) * mp.expjpi(-mp.mpf(eta))
           * le.regulated_double_pole_complex(p, eta, 1 / mp.mpf(modulus), mp.pi) / mp.mpf(eta) ** 2)
    assert abs(first + second - box) < 1e-12 * abs(box)
    # the other side is the mirror image
    f2, s2 = cl.timelike_two_edges(p, eta, modulus, sigma=1)
    assert abs(f2 + s2 + mp.conj(box)) < 1e-12 * abs(box)


@pytest.mark.parametrize("p,rho,side", [(0.3, -0.8, -1), (-0.3, -2.5, 1)])
def test_continued_piece_carries_two_units_and_a_bounded_remainder(p, rho, side):
    import mpmath as mp
    unit = -2j * mp.pi ** 2
    sinc = mp.sin(mp.pi * p) / (mp.pi * p)
    scaled = []
    for eta in (0.1, 0.05, 0.025):
        e = mp.mpf(eta)
        a = cl.edge_piece_one_dim(p, eta, rho, side=side)
        scaled.append(e * (a - 2 * unit * sinc / e ** 2) / unit)
    assert all(abs(v) < 12 for v in scaled)
    assert abs(scaled[-1] - scaled[-2]) < 0.5        # eta (A - 2 sinc/eta^2) converges
    with pytest.raises(ValueError):
        cl.edge_piece_one_dim(0.3, 0.1, -0.5)


@pytest.mark.parametrize("p,eta,T", [(0.3, 0.31, 2.2), (-0.4, 0.2, 3.0)])
def test_direct_timelike_computation_equals_the_two_continued_pieces(p, eta, T):
    """No contour, no Mellin-Barnes, no continuation in T on the left side."""
    one, two = cl.timelike_direct(p, eta, T)
    first, second = cl.timelike_two_edges(p, eta, T, sigma=1)
    assert abs(one + two - complex(first + second)) < 1e-5 * abs(first + second)
    # the other side of the cut is not what the direct computation gives
    f2, s2 = cl.timelike_two_edges(p, eta, T, sigma=-1)
    assert abs(one + two - complex(f2 + s2)) > 0.1 * abs(f2 + s2)
    with pytest.raises(ValueError):
        cl.timelike_direct(p, eta, 0.5)


@pytest.mark.parametrize("p,eta,T", [(0.3, 0.31, 2.2), (0.15, 0.2, 3.0)])
def test_timelike_keyhole_piece_is_the_second_continued_piece(p, eta, T):
    piece = cl.timelike_keyhole_piece(p, eta, T)
    first, second = cl.timelike_two_edges(p, eta, T, sigma=1)
    assert abs(piece - complex(second)) < 1e-3 * abs(second)
    # and the rest of the direct computation is the first piece
    one, two = cl.timelike_direct(p, eta, T)
    assert abs(one + two - piece - complex(first)) < 2e-3 * abs(first)


@pytest.mark.parametrize("beta,y,eps", [(0.4, (0.3, 0.5), 0.05), (1.8, (-0.9, 0.3), 0.03)])
@pytest.mark.parametrize("p", [0.3, -0.25])
def test_timelike_alpha_integral_splits_into_two_keyholes_at_positive_displacement(beta, y, eps, p):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        line, k3, k4, formula = cl.timelike_alpha_split(p, 2.2, beta, y, eps)
    assert abs(line - k3 - k4) < 1e-7 * abs(line)
    assert abs(k4 - formula) < 1e-7 * abs(k4)          # the tau formula, both signs of p
    if beta > 1:
        assert k3 == 0


@pytest.mark.parametrize("p,tol", [(0.3, 5e-4), (-0.3, 5e-3)])
def test_transverse_integral_of_the_keyhole_is_the_closed_form_at_positive_displacement(p, tol):
    """Lemma E of section 14: the exchange of integrals, with eps > 0."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid = cl.timelike_keyhole_transverse(p, 0.3, 2.2, 0.4, 0.3, n=80)
    closed = complex(cl.timelike_second_piece_at_beta(p, 0.3, 2.2, 0.4, 0.3))
    assert abs(grid - closed) < tol * abs(closed)


@pytest.mark.parametrize("p", [0.3, -0.3])
def test_closed_form_converges_as_the_displacement_is_removed(p):
    limit = cl.timelike_second_piece_at_beta(p, 0.3, 2.2, 0.4, 0.0)
    gaps = [abs(cl.timelike_second_piece_at_beta(p, 0.3, 2.2, 0.4, e) - limit) for e in (1e-2, 1e-4, 1e-6)]
    assert gaps[0] > gaps[1] > gaps[2]
    assert gaps[2] < 0.02 * abs(limit)
    # the rate is eps^eta: two decades in eps^(eta) between successive entries
    assert 2.0 < gaps[1] / gaps[2] < 6.0
