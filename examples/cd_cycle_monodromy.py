"""Measure the obstruction to a periodic CD contour on the fixed positive-x chart."""

from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from geocusp.contours import winding_number

from .complex_cross_ratio import integrated_projective


def certify_winding_obstruction():
    x, t, z = sp.symbols("x t z", nonzero=True)
    quotient = (t+z)/(t+1/x)
    assert quotient.subs(t, 0) == x*z
    assert sp.limit(quotient, t, sp.oo) == 1
    # This map extends at t=infinity. On a periodic family of paths it
    # would give a homotopy in C* between these two endpoint loops.
    circle = np.exp(2j*pi*np.arange(256)/256)
    for positive_x in (0.25, 1., 4.):
        assert winding_number(positive_x*circle) == 1
        assert winding_number(np.ones_like(circle)) == 0

    coefficient = 2*(1+z)/((1+x)*(1+t*x)*(t+z))
    assert sp.cancel(sp.residue(coefficient, t, -z)
                     -2*(1+z)/((1+x)*(1-z*x))) == 0
    v = sp.symbols("v")
    compactified = (1/((1+x)*(1-z*x))).subs(x, v/(1-v))/(1-v)**2
    assert sp.cancel(compactified-1/(1-(1+z)*v)) == 0


def meridian_period(z):
    """Period of the clockwise t meridian integrated over positive x."""
    pole = 1/(1+z)
    knots = [0., 1.]
    if 0 < pole.real < 1:
        width = abs(pole.imag)
        knots += [pole.real+k*width for k in (-8, -1, 0, 1, 8)
                  if 0 < pole.real+k*width < 1]
    knots.sort()
    residue = sum(quad(lambda v: 1/(1-(1+z)*v), a, b,
                       complex_func=True, epsabs=1e-10, epsrel=1e-10)[0]
                  for a, b in zip(knots, knots[1:]))
    return -4j*pi*(1+z)*residue


def check_cd_cycle_monodromy():
    certify_winding_obstruction()
    print("Endpoint winding of (t+z)/(t+1/x): t=0 gives +1, t=infinity gives 0")

    for angle in (-2.9, -pi/2, -0.15):
        z = np.exp(1j*angle)
        for x in (0.25, 1., 4.):
            # Position and radius of a representative meridian. Both
            # the stationary pole and the endpoint at t=0 stay outside.
            radius = min(abs(z), abs(z-1/x))/4
            loop = -z+radius*np.exp(-2j*pi*np.arange(256)/256)
            np.testing.assert_allclose(
                [winding_number(loop, point=pole) for pole in (-z, 0, -1/x)],
                [-1, 0, 0], atol=1e-12, rtol=0,
            )
        period = meridian_period(z)
        np.testing.assert_allclose(period, 4j*pi*np.log(-z), rtol=0, atol=1e-9)
        continuous = integrated_projective(z)+period
        np.testing.assert_allclose(continuous, pi*pi+(np.log(z)+2j*pi)**2,
                                   rtol=0, atol=1e-9)
        print(f"angle={angle:.3f}: clockwise meridian period={period:.9f}")

    # Approach z=1 on the returned sheet. The fixed-x poles collide at
    # x=1 exactly at the endpoint, so take the lateral limit.
    errors = []
    for angle in (-0.01, -0.001, -0.0001):
        period = meridian_period(np.exp(1j*angle))
        errors.append(abs(period+4*pi*pi))
    assert errors[2] < errors[1] < errors[0] and errors[-1] < 0.002
    print(f"One-turn scalar shift at z=1: numerical={period:.9f}, limit={-4*pi*pi:.9f}")

    # Compare using the principal branch over the whole circle with
    # continuing from phi=0 through the cut and keeping the same sheet.
    nodes, weights = np.polynomial.legendre.leggauss(40)
    positive_angles = (1+nodes)*pi/2
    upper = np.exp(1j*positive_angles)
    lower = upper.conjugate()
    principal = weights @ np.array([
        integrated_projective(z)+integrated_projective(z.conjugate()) for z in upper
    ])/4
    carried = principal + weights @ np.array([meridian_period(z) for z in lower])/4
    np.testing.assert_allclose(principal, 2*pi*pi/3, rtol=0, atol=1e-8)
    np.testing.assert_allclose(carried-principal, -pi*pi, rtol=0, atol=1e-8)
    print(f"Principal contour mean={principal.real:.9f}; "
          f"mean with continued branch={carried.real:.9f}")
    print("PASS: fixed-chart winding obstruction and CD monodromy")


if __name__ == "__main__":
    check_cd_cycle_monodromy()
