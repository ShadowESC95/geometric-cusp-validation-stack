"""Resolve the complex regulator transform at the interior spectral cusp."""

from cmath import exp, sqrt as complex_sqrt
from math import pi, sqrt

import numpy as np

from .complex_regulator import cut_boundary, regulator_transform
from .regulator_spectrum import spectral_density
from .regulator_upper_edge import cusp_coefficient, edges


def check_regulator_interior_cusp():
    for z in (0.5, 1.5, 5.):
        cusp, upper = edges(z)
        plus, minus = cut_boundary(z, cusp)
        sigma = spectral_density(cusp, z)
        np.testing.assert_allclose(plus-minus, 2j*pi*sigma, rtol=0, atol=1e-9)

        coefficient = cusp_coefficient(z)
        target = -pi*coefficient
        # Approach the cusp from both half-planes and on both sides of
        # Re(w)=H_minus. The displacement uses the principal square root.
        worst = 0.
        for angle in (pi/4, 3*pi/4, -pi/4, -3*pi/4):
            baseline = plus if angle > 0 else minus
            errors = []
            for radius in (1e-4, 1e-5, 1e-6, 1e-7):
                delta = radius*exp(1j*angle)
                value = regulator_transform(z, cusp+delta)
                scaled = (value-baseline)/complex_sqrt(delta)
                errors.append(abs(scaled/target-1))
            assert all(new < old for old, new in zip(errors, errors[1:]))
            worst = max(worst, errors[-1])

        # On the right lip the leading correction is real. On the left
        # it is purely imaginary and changes sign between the two lips.
        distance = 1e-6
        for direction in (-1, 1):
            at = cusp+direction*distance
            values = cut_boundary(z, at)
            for side, (value, baseline) in enumerate(zip(values, (plus, minus))):
                if direction > 0:
                    root = sqrt(distance)
                else:
                    root = (1-2*side)*1j*sqrt(distance)
                scaled = (value-baseline)/root
                np.testing.assert_allclose(scaled, target, rtol=0.025, atol=0.)
        assert worst < 0.01
        print(f"z={z:g}, H_minus={cusp:.6f}, H_plus={upper:.6f}: "
              f"coefficient={target:.9f}, max ray error={worst:.3e}")

    print("PASS: interior spectral cusp on both regulator boundary values")


if __name__ == "__main__":
    check_regulator_interior_cusp()
