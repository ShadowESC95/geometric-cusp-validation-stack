"""Upper-edge asymptotics and the crossover as z approaches one."""

from math import log, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import beta

from .regulator_spectrum import spectral_density, spectral_support


def edges(z):
    """Return the two maxima obtained at the angular endpoints."""
    return min(1, z)/(2*(1+z)), spectral_support(z)


def upper_coefficient(z):
    lower, upper = edges(z)
    if z == 1:
        return 2*log(3)
    return 8/(3*(upper-lower)*sqrt(upper))


def crossover(xi):
    """Limit of sqrt(H_plus-H_minus)*sigma at h=H_plus-xi*(H_plus-H_minus)."""
    right = sqrt(xi)
    left = sqrt(max(xi-1, 0))

    def integrand(v):
        q = 0.75-(xi-v*v)/2
        return 1/(q*(1-q))

    # The two allowed intervals in v are related by reflection.
    return quad(integrand, left, right, epsabs=1e-12, epsrel=1e-12)[0]


def upper_moment(z, order):
    """Integrate the upper part of the (h/H)^n moment for order >= 40."""
    lower, upper = edges(z)
    # Use u=n*(1-h/H), then u=s^2. The far region u>40 contributes
    # at most exp(-40) times the zeroth moment.
    def integrand(s):
        h = upper*(1-s*s/order)
        return (2*upper*s/order)*np.exp(order*np.log1p(-s*s/order))*spectral_density(h, z)

    inner = sqrt(order*(1-lower/upper))
    knots = sorted({0., sqrt(40.), min(inner, sqrt(40.))})
    return sum(quad(integrand, a, b, epsabs=1e-12, epsrel=1e-9)[0]
               for a, b in zip(knots, knots[1:]))


def certify_rescaling():
    z, t, q = sp.symbols("z t q", positive=True)
    upper, lower = z/(2*(1+z)), 1/(2*(1+z))
    gap = upper-lower
    h = 4*t*(1-(1+z)*t+(z-1)*q)/(1+z)
    # Exact distance to the upper edge, and to the angular endpoint at fixed t.
    assert sp.cancel(upper-h-upper*(4*t-1)**2-8*t*gap*(1-t-q)) == 0
    assert sp.cancel(lower-h.subs(q, t)-lower*(4*t-1)**2) == 0
    assert sp.cancel(8*t/(q*(1-q))/sp.diff(h, q)-1/(gap*q*(1-q))) == 0


def check_regulator_upper_edge():
    certify_rescaling()
    for z in (1., 1.5, 5.):
        lower, upper = edges(z)
        coefficient = upper_coefficient(z)
        exponent = -0.5 if z == 1 else 0.5
        errors = []
        for distance in (1e-4, 1e-5, 1e-6, 1e-7):
            measured = spectral_density(upper-distance, z)/distance**exponent
            errors.append(abs(measured/coefficient-1))
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 3e-6
        print(f"z={z:g}: edge power={exponent:+.1f}, coefficient={coefficient:.9f}, "
              f"estimate={measured:.9f}")

        ratios = []
        for order in (64, 256, 1024):
            denominator = coefficient*upper**(exponent+1)*beta(order+1, exponent+1)
            ratios.append(upper_moment(z, order)/denominator)
        assert abs(ratios[-1]-1) < 0.007
        assert abs(ratios[-1]-1) < abs(ratios[0]-1)
        print(f"z={z:g}: n=1024 moment/edge prediction={ratios[-1]:.9f}")

        if z != 1:
            baseline = spectral_density(lower, z)
            predicted = 8/(3*(upper-lower)*sqrt(lower))
            distance = 1e-7
            below = (baseline-spectral_density(lower-distance, z))/sqrt(distance)
            above = (spectral_density(lower+distance, z)-baseline)/sqrt(distance)
            np.testing.assert_allclose(below, predicted, rtol=0.004, atol=0)
            assert abs(above/predicted) < 0.004
            print(f"z={z:g}: interior cusp coefficient={predicted:.9f}, estimate={below:.9f}")

    for xi in (0.1, 0.5, 1., 1.5, 4.):
        target = crossover(xi)
        errors = []
        for z in (1.01, 1.001, 1.0001):
            lower, upper = edges(z)
            gap = upper-lower
            h = lower if xi == 1 else upper-gap*xi
            value = sqrt(gap)*spectral_density(h, z)
            errors.append(abs(value-target))
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 0.002
        print(f"xi={xi:g}: crossover={target:.9f}, z=1.0001 estimate={value:.9f}")

    np.testing.assert_allclose(crossover(1e-6)/sqrt(1e-6), 16/3, rtol=1e-5)
    np.testing.assert_allclose(sqrt(1e4)*crossover(1e4), 2*log(3), rtol=1e-4)
    print("PASS: upper spectral edge, interior cusp and the merger at z=1")


if __name__ == "__main__":
    check_regulator_upper_edge()
