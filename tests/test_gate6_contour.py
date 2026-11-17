"""Checks for docs/GATE_6_CONTOUR.md: an explicit global integration contour."""
import mpmath as mp
import numpy as np
import pytest

from geocusp import feynman_contour as fc


@pytest.fixture(scope="module")
def contour():
    return fc.Contour()


def test_projection_onto_the_triangle_has_the_obtuse_angle_property():
    """d = P - proj(P) satisfies d.(P - V) >= |d|^2 for every vertex V."""
    rng = np.random.default_rng(0)
    vertices = np.array([[0.0, 0.0], [0.0, -1.0], [1.0, 0.0]])
    for _ in range(20000):
        P = rng.normal(size=2) * 2
        d = P - fc.project_to_triangle(P)
        assert min(d @ (P - v) - d @ d for v in vertices) > -1e-12


def test_frame_field_controls_its_three_invariants_everywhere(contour):
    """V_j.(x0 - x_i) <= -(kappa|D|/2)|d|^2 - nu rho for i = j-1, j, j+1."""
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


def test_global_properties_of_the_contour(contour):
    result = fc.survey(contour, 20000, seed=3)
    assert result["near_count"] > 10000
    assert result["max_im_near"] <= 0.0                 # correct side of every cone
    assert result["max_excursion"] < 0.01               # |arg z| <= 2 pi + 0.04
    assert result["min_ratio"] > 1e-7                   # zeros only on the edges


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
