"""Checks for docs/GATE_6_CONTOUR.md: an explicit global integration contour."""
import mpmath as mp
import numpy as np
import pytest

from geocusp import feynman_contour as fc


@pytest.fixture(scope="module")
def contour():
    return fc.Contour()


REGIONS = ("--", "+-", "-+")


def test_potential_is_convex_smooth_and_vanishes_exactly_on_the_triangle():
    """d = grad Psi satisfies d.(P - V) >= Psi for every vertex V and
    |d|^2 <= 16 Psi; Psi = 0 exactly on T; and d is continuously differentiable."""
    rng = np.random.default_rng(0)
    vertices = np.array([[0.0, 0.0], [0.0, -1.0], [1.0, 0.0]])
    for _ in range(20000):
        P = rng.normal(size=2) * 10 ** rng.uniform(-3, 1.2)
        d, psi = fc.triangle_gradient(P), fc.triangle_potential(P)
        assert min(d @ (P - v) for v in vertices) >= psi - 1e-12
        assert d @ d <= 16 * psi + 1e-12
        inside = P[0] >= 0 and P[1] <= 0 and P[0] - P[1] <= 1
        assert (psi == 0) == inside
    jac = lambda P, h=1e-7: np.array([(fc.triangle_gradient(P + h * e) - fc.triangle_gradient(P - h * e))
                                      / (2 * h) for e in np.eye(2)])
    worst = 0.0
    for _ in range(3000):
        P = rng.normal(size=2) * 2
        worst = max(worst, np.abs(jac(P) - jac(P + rng.normal(size=2) * 1e-5)).max())
    assert worst < 1e-3                                # no jump in the derivative
    assert np.linalg.norm(jac(np.array([3.0, 4.0])), 2) <= fc.HESSIAN_BOUND + 1e-6


def test_frame_field_controls_its_three_invariants_everywhere(contour):
    """V_j.(x0 - x_i) <= -(kappa|D|/2) Psi - nu rho for i = j-1, j, j+1."""
    rng = np.random.default_rng(1)
    worst = -np.inf
    for _ in range(5000):
        x0 = contour.centre + rng.normal(size=4) * 10 ** rng.uniform(-1, 1.5)
        for j in range(4):
            v = contour.frame_field(j, x0)
            bound = contour.frame_bound(j, x0)
            for i in ((j - 1) % 4, j, (j + 1) % 4):
                worst = max(worst, fc.mdot(v, x0 - contour.x[i]) - bound)
    assert worst < 1e-9


def test_no_real_point_is_close_to_all_four_cones(contour):
    """So some frame field is always fully switched on near a cone."""
    from scipy.optimize import minimize

    rng = np.random.default_rng(2)
    best = np.inf
    for _ in range(80):
        res = minimize(lambda x: np.max(np.abs(contour.invariants(x))),
                       contour.centre + rng.normal(size=4) * 4, method="Nelder-Mead",
                       options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        best = min(best, res.fun)
    assert best > 2 * contour.delta0
    # the analytic bound used in the proof is weaker than the true minimum but suffices
    assert 2 * contour.delta0 < fc.four_cone_gap(contour) < best


def test_global_properties_of_the_contour(contour):
    result = fc.survey(contour, 20000, seed=3)
    assert result["near_count"] > 10000
    assert result["max_im_near"] <= 0.0                 # correct side of every cone
    assert result["max_excursion"] < 0.01               # |arg z| <= 2 pi + 0.04
    assert result["min_ratio"] > 1e-9                   # zeros only on the edges


def test_displacement_is_small_and_the_map_is_a_local_diffeomorphism(contour):
    """The Jacobian of x0 -> x0 + i V is 1 + i dV with dV real. Its
    determinant can vanish only if dV has an eigenvalue of modulus 1, so a
    spectral norm below 1 is enough. The largest values come from the
    derivative of the cutoffs that switch the frame fields."""
    worst = 0.0
    for x0 in fc.sample_points(contour, 1500, seed=4):
        if contour.distance_to_edges(x0) < 1e-3:
            continue
        h = 1e-6 * max(1.0, np.linalg.norm(x0))
        jac = np.array([(contour.field(x0 + h * e) - contour.field(x0 - h * e)) / (2 * h)
                        for e in np.eye(4)])
        worst = max(worst, np.linalg.norm(jac, 2))
    assert worst < 0.5


@pytest.mark.parametrize("cusp", [0, 1, 2, 3])
def test_lower_bounds_that_give_absolute_integrability(contour, cusp):
    """Near each cusp every nearby invariant is at least c nu times the scale
    max(rho, ...) used in the convergence estimate."""
    assert fc.edge_lower_bounds(contour, 6000, seed=5 + cusp, cusp=cusp) > 0.5


def test_elementary_gamma_bounds():
    for a in (0.05, 0.3, 0.7, 1.0):
        for v in (0.0, 0.2, 1.0, 3.0, 9.0):
            value, bound = fc.gamma_product_bounds(a, v)
            assert value <= bound * (1 + mp.mpf("1e-12"))
    for a in (1.0, 1.4, 1.95, 2.0):
        for v in (0.0, 0.2, 1.0, 3.0, 9.0):
            value, bound = fc.gamma_product_bounds(a, v)
            assert value <= bound * (1 + mp.mpf("1e-12"))
    # the shift used for the weight: |Gamma(a + i y)|^2 = |Gamma(1 + a + i y)|^2 / (a^2 + y^2)
    a, y = mp.mpf("0.01"), mp.mpf("0.37")
    lhs = abs(mp.gamma(mp.mpc(a, y))) ** 2
    rhs = abs(mp.gamma(mp.mpc(1 + a, y))) ** 2 / (a**2 + y**2)
    assert abs(lhs - rhs) < mp.mpf("1e-12") * lhs


@pytest.mark.parametrize("region", REGIONS)
def test_every_explicit_condition_of_the_proof_holds(region):
    """Section 9 of GATE_6_CONTOUR reduces Proposition 6 to a list of
    inequalities between named constants. They hold for the contour used."""
    result = fc.certificate(fc.Contour(region=region))
    failed = [k for k, v in result.items() if isinstance(v, bool) and not v]
    assert not failed, failed
    assert result["jacobian_bound"] < 1 and result["excursion_bound"] < 0.2


@pytest.mark.parametrize("region", REGIONS)
def test_contour_in_every_kinematic_region(region):
    """One sign changes where a diagonal is timelike; the properties are the same."""
    C = fc.Contour(region=region)
    rng = np.random.default_rng(6)
    worst = -np.inf
    for _ in range(2000):
        x0 = C.centre + rng.normal(size=4) * 10 ** rng.uniform(-1, 1.5)
        for j in range(4):
            v, bound = C.frame_field(j, x0), C.frame_bound(j, x0)
            for i in ((j - 1) % 4, j, (j + 1) % 4):
                worst = max(worst, fc.mdot(v, x0 - C.x[i]) - bound)
    assert worst < 1e-12
    result = fc.survey(C, 8000, seed=8)
    assert result["max_im_near"] <= 0.0 and result["near_count"] > 4000
    assert result["max_excursion"] < fc.certificate(C)["excursion_bound"]
    assert result["min_ratio"] > 1e-9


def test_measured_values_respect_the_proved_bounds(contour):
    cert = fc.certificate(contour)
    worst = 0.0
    for x0 in fc.sample_points(contour, 1500, seed=9):
        if contour.distance_to_edges(x0) < 1e-3:
            continue
        h = 1e-6 * max(1.0, np.linalg.norm(x0))
        jac = np.array([(contour.field(x0 + h * e) - contour.field(x0 - h * e)) / (2 * h)
                        for e in np.eye(4)])
        worst = max(worst, np.linalg.norm(jac, 2))
    assert worst <= cert["jacobian_bound"]
    assert fc.survey(contour, 6000, seed=10)["max_excursion"] <= cert["excursion_bound"]


# ---- every square: constants chosen from the geometry ---------------------
@pytest.mark.parametrize("region", REGIONS)
def test_certificate_holds_across_the_family_of_squares(region):
    """Up to Poincare transformations and scale a light-like square is fixed
    by its region and one angle. With constants derived from the geometry
    every condition of the proof holds across the range of angles."""
    ratios = []
    for angle in np.linspace(0.15, 2.99, 30):
        C = fc.tuned_contour(region, angle)
        result = fc.certificate(C)
        failed = [k for k, v in result.items() if isinstance(v, bool) and not v]
        assert not failed, (angle, failed)
        assert result["jacobian_bound"] < 0.5 and result["excursion_bound"] < 0.1
        ratios.append(C.s / C.t)
    assert max(np.abs(ratios)) / min(np.abs(ratios)) > 150      # two decades of cross-ratio


def test_certificate_is_scale_invariant():
    for scale in (0.01, 37.0):
        assert fc.certificate(fc.tuned_contour("--", 0.7, scale=scale))["all"]


@pytest.mark.parametrize("region,angle", [("--", 0.2), ("--", 2.9), ("+-", 0.3), ("-+", 2.8)])
def test_extreme_squares_have_the_same_properties(region, angle):
    C = fc.tuned_contour(region, angle)
    result = fc.survey(C, 5000, seed=5)
    assert result["max_im_near"] <= 0 and result["min_ratio"] > 0
    assert result["max_excursion"] <= fc.certificate(C)["excursion_bound"]


# ---- absolute integrability, line by line ---------------------------------
@pytest.mark.parametrize("region,angle", [("--", 1.1), ("+-", 1.1), ("-+", 1.1), ("--", 0.3)])
def test_lemma_6_displaced_invariant_dominates_its_real_part_and_its_imaginary_part(region, angle):
    """Near a cone, |V^2| <= (theta/2)|Im x0i^2| and so
    |x0i^2 displaced| >= (1/2) max(|x0i^2|, |Im x0i^2|)."""
    C = fc.tuned_contour(region, angle)
    result = fc.lemma6_check(C, 8000, seed=1)
    assert result["count"] > 2000 and result["theta"] < 1e-3
    assert result["max_v2_over_im"] <= result["theta"] / 2
    assert result["min_ratio"] >= 0.5


@pytest.mark.parametrize("region,angle", [("--", 1.1), ("+-", 1.1), ("-+", 1.1), ("+-", 2.6)])
def test_bounds_on_the_strips_around_the_edges(region, angle):
    C = fc.tuned_contour(region, angle)
    result = fc.strip_check(C, 8000, seed=2)
    assert result["frame_on"] == 1.0                    # the cusp's frame is fully on
    for key in ("j", "next", "prev", "far"):
        assert result[key] >= 1.0, key
    assert 0 < result["b0"] < 0.5                       # the two strips of an edge overlap


def test_the_strip_integral_converges():
    """int d alpha d beta d rho rho^(eta-1) / (max(rho,|alpha|) max(rho,|beta|))
    <= C int d rho rho^(eta-1) (1 + log(1/rho))^2 < infinity, for eta > 0."""
    import mpmath as mp

    eta = mp.mpf("0.05")
    inner = lambda rho: (2 * (1 + mp.log(1 / rho))) ** 2          # exact value of the two inner integrals
    # rho = exp(-x) turns the endpoint singularity into exponential decay
    total = mp.quad(lambda x: mp.exp(-eta * x) * inner(mp.exp(-x)), [0, 10, 100, 1000, mp.inf])
    closed = 4 * (1 / eta + 2 / eta**2 + 2 / eta**3)             # int rho^(eta-1) 4 (1+L)^2
    assert abs(total - closed) < 1e-8 * closed


@pytest.mark.parametrize("region", REGIONS)
def test_lower_bound_at_infinity(region):
    """For r >= 2R, |x0i^2 displaced| >= c kappa' max(r, |(x0 - c)^2|)."""
    assert fc.infinity_check(fc.tuned_contour(region, 1.1), 6000, seed=3) > 1.0


@pytest.mark.parametrize("region", REGIONS)
def test_certificate_covers_the_strips_and_the_bound_at_infinity(region):
    """The widths of the strips of Lemma 7 and the constant of Lemma 9 are
    conditions of the certificate, and the sampled bounds agree with them."""
    C = fc.tuned_contour(region, 1.1)
    cert = fc.certificate(C)
    for key in ("F_theta_small", "S0_strips_overlap", "S1_frame_fully_on", "S2_near_invariants",
                "S3_strips_inside", "I0_c_u_positive", "I1_real_part", "I2_theta_is_one"):
        assert cert[key], key
    strips = fc.strip_check(C, 3000, seed=2)
    assert strips["b0"] == pytest.approx(cert["b0"]) and strips["a1"] == pytest.approx(cert["a1"])
    assert strips["far"] >= 2.0                         # |x0,j+2^2| >= 2 delta0, as (S) asserts
    assert fc.lemma6_check(C, 3000, seed=1)["max_v2_over_im"] <= cert["theta"] / 2


def test_field_is_continuously_differentiable(contour):
    """The Jacobian of V has no jumps: its change over a short step is small
    everywhere on a stress sample, including near the triangles and cutoffs."""
    rng = np.random.default_rng(12)
    worst = 0.0
    for x0 in fc.sample_points(contour, 600, seed=12):
        if contour.distance_to_edges(x0) < 1e-3:
            continue
        scale = max(1.0, np.linalg.norm(x0))
        h = 1e-6 * scale
        jac = lambda x: np.array([(contour.field(x + h * e) - contour.field(x - h * e)) / (2 * h)
                                  for e in np.eye(4)])
        step = rng.normal(size=4) * 1e-4 * scale
        worst = max(worst, np.abs(jac(x0) - jac(x0 + step)).max())
    assert worst < 5e-3


# ---- repairs after review -------------------------------------------------
@pytest.mark.parametrize("region,low,high", [("--", -2, 2), ("+-", -1, 2), ("-+", -2, 1)])
def test_phase_of_z_at_real_points_and_on_the_contour(region, low, high):
    """Lemma P: at real points Im log z / pi is an integer in the stated
    range. On the contour it stays within the same range plus a margin."""
    out = fc.phase_survey(fc.tuned_contour(region, 1.1), 6000, seed=6)
    assert low - 1e-9 <= out["real_min"] and out["real_max"] <= high + 1e-9
    assert low - 0.06 <= out["contour_min"] and out["contour_max"] <= high + 0.06


def test_causal_fact_behind_lemma_p():
    """If x13 is timelike and x0 is timelike from x2, then x0 is timelike
    from x1 or from x3."""
    from geocusp import lorentzian_edges as le
    p, s, t = le.light_like_square("+-")
    assert s > 0
    x = [np.zeros(4), p[0], p[0] + p[1], p[0] + p[1] + p[2]]
    rng = np.random.default_rng(1)
    hits = 0
    for _ in range(20000):
        x0 = x[1] + rng.normal(size=4) * rng.choice([0.3, 1.0, 3.0])
        if le.msq(x0 - x[1]) > 0:
            hits += 1
            assert le.msq(x0 - x[0]) > 0 or le.msq(x0 - x[2]) > 0
    assert hits > 1000


def test_strip_integral_with_powers_converges_only_below_the_regulator():
    """Lemma 8': with propagator powers 1 -+ sigma the strip bound is
    integrable if and only if |sigma| < eta."""
    eta = 0.2
    inside = [fc.strip_power_integral(eta, 0.1, c) for c in (1e-12, 1e-24, 1e-48)]
    outside = [fc.strip_power_integral(eta, 0.3, c) for c in (1e-12, 1e-24, 1e-48)]
    assert inside[2] - inside[1] < 0.2 * (inside[1] - inside[0])      # increments shrink
    assert outside[2] - outside[1] > 5 * (outside[1] - outside[0])    # increments grow
