"""Continue the critical kernel in the regulator using its positive spectrum."""

from cmath import exp, sqrt as complex_sqrt
from math import log, pi

import numpy as np
from scipy.integrate import quad

from .critical_regulator_kernel import critical_kernel
from .regulator_spectrum import (
    density_from_edge, spectral_density, spectral_integral, spectral_support,
)
from .regulator_upper_edge import upper_coefficient


def regulator_transform(z, w):
    """Cauchy transform with w=i*eta, away from the real spectral cut."""
    upper = spectral_support(z)
    if w.imag == 0 and 0 <= w.real <= upper:
        raise ValueError("Use a boundary value on the spectral cut")
    offsets = (-64, -8, -1, 0, 1, 8, 64)
    breakpoints = (w.real+offset*abs(w.imag) for offset in offsets)
    return spectral_integral(z, lambda h: 1/(h-w), breakpoints=breakpoints)


def cut_boundary(z, w):
    """Return the upper and lower boundary values at an interior cut point."""
    upper = spectral_support(z)
    if not 0 < w < upper:
        raise ValueError("The cut point must lie strictly inside the support")
    sigma = spectral_density(w, z)
    cusp = min(1, z)/(2*(1+z))
    # Subtract the value at the pole, and add back the constant PV analytically.
    # The square-root substitution also resolves the divergent upper edge at z=1.
    pole_distance = upper-w

    def subtracted(s):
        distance = upper*s*s
        return 2*upper*s*(density_from_edge(distance, z)-sigma)/(pole_distance-distance)

    s_pole = np.sqrt(1-w/upper)
    knots = sorted({0., s_pole, np.sqrt(1-cusp/upper), 1.})
    pv = sigma*log((upper-w)/w)
    pv += sum(quad(subtracted, lo, hi, epsabs=1e-9, epsrel=1e-9)[0]
              for lo, hi in zip(knots, knots[1:]))
    return pv+1j*pi*sigma, pv-1j*pi*sigma


def check_complex_regulator():
    for z in (1., 1.5, 5.):
        upper = spectral_support(z)
        np.testing.assert_allclose(regulator_transform(z, 0.02j),
                                   critical_kernel(z, 0.02), rtol=0, atol=1e-9)
        for fraction in (0.37, 0.88):
            w = upper*fraction
            plus, minus = cut_boundary(z, w)
            jump = plus-minus
            np.testing.assert_allclose(jump, 2j*pi*spectral_density(w, z),
                                       rtol=0, atol=1e-9)
            errors = []
            for height in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
                numerically_above = regulator_transform(z, w+1j*height)
                # Conjugation follows directly from the real spectrum.
                numerically_below = numerically_above.conjugate()
                errors.append(max(abs(numerically_above-plus), abs(numerically_below-minus)))
            assert all(new < old for old, new in zip(errors, errors[1:]))
            assert errors[-1] < 5e-4
            np.testing.assert_allclose(
                regulator_transform(z, w-1j*height), numerically_below, rtol=0, atol=1e-9,
            )
            print(f"z={z:g}, w/H={fraction:g}: cut jump={jump:.9f}, "
                  f"boundary error={errors[-1]:.3e}")

        C = upper_coefficient(z)
        target = -pi*C if z == 1 else pi*C
        for angle in (0., pi/3, -pi/3):
            errors = []
            for radius in (1e-4, 1e-5, 1e-6, 1e-7):
                delta = radius*exp(1j*angle)
                root = complex_sqrt(delta)

                def integrand(s):
                    d = upper*s*s
                    weight = 2*upper*s*density_from_edge(d, z)
                    if z == 1:
                        return -root*weight/(d+delta)
                    return root*weight/(d*(d+delta))

                cusp = min(1, z)/(2*(1+z))
                knots = sorted({0., 1., np.sqrt(1-cusp/upper),
                                *(np.sqrt(m*radius/upper) for m in (1, 8, 64))})
                result = sum(quad(integrand, lo, hi, complex_func=True,
                                  epsabs=1e-8, epsrel=1e-8)[0]
                             for lo, hi in zip(knots, knots[1:]))
                errors.append(abs(result-target)/abs(target))
            assert all(new < old for old, new in zip(errors, errors[1:]))
            assert errors[-1] < 0.004
        print(f"z={z:g}: endpoint coefficient={target:.9f}, "
              f"scaled transform={result:.9f}, relative error={errors[-1]:.3e}")

    print("PASS: complex regulator cut and upper-endpoint singularity")


if __name__ == "__main__":
    check_complex_regulator()
