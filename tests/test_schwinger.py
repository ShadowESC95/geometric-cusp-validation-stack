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
