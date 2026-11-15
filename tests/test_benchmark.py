import numpy as np

from geocusp.benchmark import (
    Fbar0,
    Fbar1,
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
