import numpy as np

import pytest

from geocusp.benchmark import (
    Fbar0,
    Fbar1,
    Fbar2,
    circle_mean_mp,
    contour_I_Fbar,
    contour_I_insertion,
    cusp_coefficient_from_mean,
    insertion_F,
)
from geocusp.contours import benchmark_F0, benchmark_F1, contour_average


def test_source_contour_reproduces_first_two_cusp_coefficients():
    mean0, _, _ = contour_average(benchmark_F0)
    mean1, _, _ = contour_average(benchmark_F1)
    np.testing.assert_allclose(
        [cusp_coefficient_from_mean(1, mean0), cusp_coefficient_from_mean(2, mean1)],
        [4, -4 * np.pi**2 / 3],
        rtol=1e-7,
        atol=1e-12,
    )


def test_insertion_and_normalized_contour_agree():
    g = 0.2
    for order, normalized in enumerate((Fbar0, Fbar1)):
        norm_action, _, _ = contour_I_Fbar(g, normalized)
        insertion_action, _, _ = contour_I_insertion(
            lambda z: insertion_F(g, z, order)
        )
        np.testing.assert_allclose(insertion_action, norm_action, atol=1e-12)


def test_three_loop_coefficient_from_published_two_loop_function():
    """Gamma^(3) = 44 pi^4/45 from Fbar2 through the rule fixed at one loop.

    Nothing in the rule was adjusted for this order, so the test can fail.
    """
    mean = circle_mean_mp(Fbar2)
    np.testing.assert_allclose(mean, -11 * np.pi**4 / 15, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(
        cusp_coefficient_from_mean(3, mean), 44 * np.pi**4 / 45, rtol=1e-12
    )


def test_mp_circle_mean_reproduces_lower_orders():
    np.testing.assert_allclose(circle_mean_mp(lambda z: -1 + 0*z), -1, atol=1e-14)
    np.testing.assert_allclose(
        circle_mean_mp(lambda z: Fbar1(complex(z))), 2 * np.pi**2 / 3, rtol=1e-12
    )


@pytest.mark.parametrize("z", [0.2, 0.37, 1.0, 3.5])
def test_finite_functions_are_symmetric_under_inversion(z):
    """Cyclic symmetry of the Wilson loop requires F(z) = F(1/z)."""
    np.testing.assert_allclose(Fbar1(z), Fbar1(1 / z), rtol=1e-13)
    np.testing.assert_allclose(
        complex(Fbar2(z)), complex(Fbar2(1 / z)), rtol=1e-13, atol=1e-13
    )
