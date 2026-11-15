"""Positive spectral density for the critical regulator kernel."""

from math import log, log1p, pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .critical_regulator_kernel import critical_kernel
from .large_regulator_ratio import critical_coefficients


def spectral_support(z):
    """Upper edge of the pushforward of the folded critical kernel."""
    return max(1, z)/(2*(1+z))


def spectral_density(h, z):
    """Density of the positive measure in F(z,eta)=int dh*density/(h-i*eta)."""
    z = max(z, 1/z)
    upper_edge = spectral_support(z)
    if not 0 < h < upper_edge:
        return 0.

    def roots(g):
        discriminant = sqrt(max(0., 1-2*(1+z)*h/g))
        lo = (1+z)*h/(2*g*(1+discriminant))
        hi = (1+discriminant)/4
        return lo, hi, discriminant

    if z == 1:
        a, b, discriminant = roots(1)
        # h=4*t*(1-2*t)/(1+z) is independent of q at z=1.
        return 8*(a*(log1p(-a)-log(a))+b*(log1p(-b)-log(b)))/discriminant

    a, b, _ = roots(z)
    intervals = ((a, b),)
    if h < 1/(2*(1+z)):
        c, d, _ = roots(1)
        intervals = ((a, c), (d, b))

    def angular_factor(t):
        q = t+((1+z)*h/(4*t)-(1-2*t))/(z-1)
        return 1/q+1/(1-q)

    # For z>1, integrating against delta(h-h(t,q)) cancels the t
    # in the weight; the remaining integral is along the allowed t intervals.
    return 2*(1+z)/(z-1)*sum(
        quad(angular_factor, lo, hi, epsabs=1e-11, epsrel=1e-10)[0]
        for lo, hi in intervals
    )


def spectral_integral(z, kernel):
    """Integrate a function of h against the positive spectral density."""
    upper = spectral_support(z)
    inner_edge = min(1, z)/(2*(1+z))
    # h=upper*(1-s^2) handles the square-root edge at z=1.
    def density(s):
        h = upper*(1-s)*(1+s)
        return 2*upper*s*spectral_density(h, z)*kernel(h)

    knots = sorted({0., sqrt(1-inner_edge/upper), 1.})
    return sum(quad(density, lo, hi, epsabs=1e-9, epsrel=1e-9,
                    complex_func=True)[0]
               for lo, hi in zip(knots, knots[1:]))


def certify_transform():
    z, t, q, eta = sp.symbols("z t q eta", positive=True)
    g = 1-(1+z)*t+(z-1)*q
    h = 4*t*g/(1+z)
    weight = 8*t/(q*(1-q))
    assert sp.cancel(
        weight/(h-sp.I*eta)-2*(1+z)/(q*(1-q)*(g-sp.I*eta*(1+z)/(4*t)))
    ) == 0
    assert sp.cancel(h.subs(q, t)-4*t*(1-2*t)/(1+z)) == 0
    assert sp.cancel(h.subs(q, 1-t)-4*z*t*(1-2*t)/(1+z)) == 0
    assert sp.cancel(sp.diff(h, q)-4*t*(z-1)/(1+z)) == 0
    assert sp.cancel(h-h.subs({q: 1-q, z: 1/z}, simultaneous=True)) == 0
    primitive = 8*t*sp.log(q/(1-q))
    assert sp.cancel(sp.diff(primitive, q)-weight) == 0
    assert (primitive.subs({q: sp.Rational(3, 4), t: sp.Rational(1, 4)})
            -primitive.subs({q: sp.Rational(1, 4), t: sp.Rational(1, 4)})) == 4*sp.log(3)


def check_regulator_spectrum():
    certify_transform()
    mass, first, second = critical_coefficients()
    regulators = (0.001, 0.02, 0.2, 1., 64.)
    for z in (0.5, 1., 1.5, 5.):
        upper = spectral_support(z)
        inverse = pi*pi+log(z)**2
        target_moments = (inverse, mass, first, second(z))
        moments = [spectral_integral(z, lambda h: h**n) for n in (-1, 0, 1, 2)]
        np.testing.assert_allclose(moments, target_moments, rtol=0, atol=1e-8)

        previous_real = inverse
        previous_scaled_imag = previous_scaled_real = 0.
        mismatches = []
        for eta in regulators:
            F = spectral_integral(z, lambda h: 1/(h-1j*eta))
            mismatches.append(abs(F-critical_kernel(z, eta)))
            assert 0 < F.real < previous_real
            assert mass/(1+(upper/eta)**2) <= eta*F.imag <= mass
            assert first/(1+(upper/eta)**2) <= eta**2*F.real <= first
            assert eta*F.imag > previous_scaled_imag
            assert eta**2*F.real > previous_scaled_real
            previous_real, previous_scaled_imag, previous_scaled_real = (
                F.real, eta*F.imag, eta**2*F.real,
            )
        assert max(mismatches) < 1e-8
        print(f"z={z:g}: support=[0,{upper:.9f}], mass={moments[1].real:.9f}, "
              f"max transform mismatch={max(mismatches):.3e}")

    print("PASS: positive regulator spectrum, moment sum rules and uniform bounds")


if __name__ == "__main__":
    check_regulator_spectrum()
