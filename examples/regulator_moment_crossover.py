"""High-order spectral moments near the merger of the two edge locations."""

from math import log, pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import beta

from .regulator_upper_edge import upper_coefficient, upper_moment
from .regulator_spectrum import spectral_support


def moment_profile(gap_parameter):
    """Return G and J in sqrt(n)*M_n = G(lambda) + J(lambda)/n + o(1/n)."""
    # Change variable to u=v/(1+lambda) to resolve the large-lambda endpoint.
    scale = 1+gap_parameter

    def density(v, correction):
        u = v/scale
        q = 0.75-u/2
        f = 1/q+1/(1-q)
        value = f
        if correction:
            c = (2*u-1)/4
            f_second = 2/q**3+2/(1-q)**3
            value = c*c*f_second/4-(7/8+(gap_parameter*u)**2/2)*f
        return np.exp(-gap_parameter*u)*value

    # At v>80 the omitted exponential tail is below quadrature tolerance.
    return tuple(sqrt(pi)/(4*scale)*quad(
        density, 0, min(scale, 80), args=(correction,), epsabs=1e-11, epsrel=1e-11,
    )[0] for correction in (False, True))


def direct_moment(order, gap_parameter):
    """Compute sqrt(n)*M_n from the original folded pushforward coordinates."""
    delta = gap_parameter/order
    # s=sqrt(n)*(4*t-1); the omitted Gaussian tail is negligible at this scale.
    limit = min(8., sqrt(order))
    v, weights = np.polynomial.legendre.leggauss(64)
    u = (1+v)/2
    angular_power = np.exp(order*np.log1p(-delta*u))

    def integrand(s):
        radial = 1-s*s/order
        t = (1+s/sqrt(order))/4
        q = 1-t-u*(1-2*t)
        return radial**(order+1)*(weights @ (angular_power/(q*(1-q))))/8

    return quad(integrand, -limit, limit, epsabs=1e-12, epsrel=1e-11)[0]


def certify_factorization():
    t, u, z = sp.symbols("t u z", positive=True)
    q = 1-t-u*(1-2*t)
    h = 4*t*(1-(1+z)*t+(z-1)*q)/(1+z)
    H = z/(2*(1+z))
    delta = (z-1)/z
    assert sp.cancel(h/H-(1-(4*t-1)**2)*(1-delta*u)) == 0
    assert sp.expand(8*t*(-sp.diff(q, u))-(1-(4*t-1)**2)) == 0

    s, a, x = sp.symbols("s a x", real=True)
    exponent = sp.series((1/x+1)*sp.log(1-x*s*s)+sp.log(1-x*a)/x, x, 0, 2)
    assert sp.expand(exponent.removeO()+s*s+a+x*(s*s+s**4/2+a*a/2)) == 0
    gaussian = sp.exp(-s*s)
    assert sp.integrate(gaussian*s*s, (s, -sp.oo, sp.oo)) == sp.sqrt(sp.pi)/2
    assert sp.integrate(gaussian*s**4, (s, -sp.oo, sp.oo)) == 3*sp.sqrt(sp.pi)/4


def check_regulator_moment_crossover():
    certify_factorization()
    for gap_parameter in (0., 0.5, 2., 8., 24.):
        profile, correction = moment_profile(gap_parameter)
        errors, corrected_errors = [], []
        for order in (64, 256, 1024):
            value = direct_moment(order, gap_parameter)
            errors.append(abs(value-profile))
            corrected_errors.append(abs(value-profile-correction/order))
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert all(new < old for old, new in zip(corrected_errors, corrected_errors[1:]))
        assert corrected_errors[-1] < 1e-5
        assert corrected_errors[-1] < errors[-1]/100
        z = 1/(1-gap_parameter/order)
        spectral = sqrt(order)*upper_moment(z, order)
        np.testing.assert_allclose(spectral, value, rtol=0, atol=2e-8)
        print(f"lambda={gap_parameter:g}: G={profile:.9f}, J={correction:.9f}, "
              f"n={order} moment={value:.9f}, corrected error={corrected_errors[-1]:.3e}")

    at_zero, _ = moment_profile(0)
    np.testing.assert_allclose(at_zero, sqrt(pi)*log(3), rtol=0, atol=1e-12)
    far = 8192
    large, _ = moment_profile(far)
    np.testing.assert_allclose(far*large, 4*sqrt(pi)/3, rtol=3e-4, atol=0)

    # At fixed z>1 the large-lambda prediction reproduces the n^-3/2
    # coefficient obtained directly from the upper spectral edge.
    z, order = 1.5, 1024
    H = spectral_support(z)
    delta = 1-1/z
    edge = upper_coefficient(z)*H**1.5*beta(order+1, 1.5)
    matched = (4*sqrt(pi)/3)/(delta*order**1.5)
    np.testing.assert_allclose(edge, matched, rtol=0.002, atol=0)
    print("PASS: high-order moments across the merging spectral edges")


if __name__ == "__main__":
    check_regulator_moment_crossover()
