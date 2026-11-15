"""Integrate the excised projective CD family before reducing its scalar weight."""

from math import atan2, cos, log, pi, sin

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .cd_projective_infinity import rotated_integral


def certify_circle_reduction():
    z, u, x = sp.symbols("z u x", positive=True)
    kernel = 2*(1+z)/(z*(1+u*z)*(1+x)*(x+u))
    scalar = kernel*(1+x)*(x+u)/2
    assert sp.residue(scalar, z, 0) == 1
    assert sp.cancel((z+1/u)*scalar).subs(z, -1/u) == 1/u-1

    phase = sp.symbols("phi", real=True)
    numerator_norm_sq = 2+2*sp.cos(phase)
    denominator_norm_sq = 1+u*u+2*u*sp.cos(phase)
    assert sp.expand(denominator_norm_sq-u*numerator_norm_sq-(1-u)**2) == 0


def x_integral(u):
    """Integrate 1/((1+x)*(x+u)) on the positive x axis."""
    if u == 1.:
        return 1.
    d = u-1
    return (np.log1p(d) if abs(d) < .25 else log(u))/d


def retained_arc_weight(u, half_angle):
    """Integral of (1+z)/(1+u*z) with dphi/(2*pi) on the retained arc."""
    phase = atan2(u*sin(half_angle), 1+u*cos(half_angle))
    return (half_angle+(1-u)*phase/u)/pi


def compactified_period(half_angle):
    """Perform the base arc integral before the positive CD integrations."""
    # Fold u>1 back onto (0,1), using x_integral(1/u)=u*x_integral(u).
    def folded(u):
        return 2*x_integral(u)*(retained_arc_weight(u, half_angle)
                               +retained_arc_weight(1/u, half_angle)/u)

    return quad(folded, 0, 1, epsabs=1e-11, epsrel=1e-11)[0]


def check_cd_projective_period():
    certify_circle_reduction()
    for u in (.2, .9, 1.1, 5.):
        half_angle = pi-.11
        measured = quad(lambda phi: (1+np.exp(1j*phi))/(1+u*np.exp(1j*phi)),
                        -half_angle, half_angle, complex_func=True,
                        epsabs=1e-11, epsrel=1e-11)[0]/(2*pi)
        np.testing.assert_allclose(measured, retained_arc_weight(u, half_angle),
                                   rtol=0, atol=1e-11)
        np.testing.assert_allclose(
            quad(lambda x: 1/((1+x)*(x+u)), 0, np.inf)[0],
            x_integral(u), rtol=0, atol=1e-10,
        )

    # The second z residue is inside precisely when u>1. Folding that
    # region gives the same integral as for u<1.
    period = 4*quad(x_integral, 0, 1, epsabs=1e-11, epsrel=1e-11)[0]
    np.testing.assert_allclose(period, 2*pi*pi/3, rtol=0, atol=1e-10)

    # Absolute domination uses |(1+z)/(1+u*z)| <= 1/sqrt(u).
    # The integral of x_integral(u)/sqrt(u) over the half-line gives pi^2.
    # Folding about one and substituting u=v^2 removes the square root.
    mass = 4*quad(lambda v: x_integral(v*v), 0, 1,
                  epsabs=1e-11, epsrel=1e-11)[0]
    np.testing.assert_allclose(mass, pi*pi, rtol=0, atol=1e-10)
    print(f"Circle-first CD period={period:.12f}, majorant mass={mass:.12f}")

    previous = np.inf
    for width in (.2, .08, .025, .008):
        retained = compactified_period(pi-width)
        omitted = period-retained
        assert 0 < omitted < previous
        assert omitted <= 2*pi*width
        # Integrating the independently verified scalar benchmark across
        # the two omitted arcs gives this exact finite-width expression.
        expected = width**2-width**3/(3*pi)
        np.testing.assert_allclose(omitted, expected, rtol=0, atol=1e-9)
        previous = omitted
        print(f"Seam half-width={width:g}: omitted period={omitted:.9e}, "
              f"omitted/width^2={omitted/width**2:.9f}")

    # Check the truncated result with the CD-first integration order.
    nodes, weights = np.polynomial.legendre.leggauss(48)
    half_angle = pi-.08
    cd_first = half_angle/(2*pi)*(weights @ np.array([
        rotated_integral(np.exp(1j*half_angle*t)) for t in nodes
    ]))
    np.testing.assert_allclose(cd_first, compactified_period(half_angle),
                               rtol=0, atol=1e-8)
    print("PASS: circle-first period and the shrinking seam contribution")


if __name__ == "__main__":
    check_cd_projective_period()
