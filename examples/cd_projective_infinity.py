"""Follow the fixed-corner winding through a projective chart crossing."""

from math import log, pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from geocusp.contours import winding_number

from .complex_cross_ratio import integrated_projective


def certify_rotation():
    x, y, u, z, v = sp.symbols("x y u z v", nonzero=True)
    kernel = 2*(1+z)/((1+x)*(1+y)*(z*x+y))
    rotated = sp.cancel(z*kernel.subs(y, z*u))
    assert sp.cancel(rotated-2*(1+z)/((1+x)*(1+z*u)*(x+u))) == 0
    # The right-hand side is the residue of the rotated form in the u variable.
    assert sp.cancel((u+1/z)*rotated
                     -2*(1+z)/(z*(1+x)*(x+u))) == 0

    # Follow a diagonal connecting the two nonsingular corners of the
    # compactified positive square, with y/z=u and x=1/u.
    u_of_v = v/(1-v)
    X = 1-v
    Y = z*u_of_v/(1+z*u_of_v)
    D = sp.cancel(z*X*(1-Y)+Y*(1-X))
    infinity_factor = 1-v+z*v
    numerator = z*((1-v)**2+v**2)
    assert sp.cancel(D*infinity_factor-numerator) == 0
    assert sp.cancel(D.subs(v, 0)-z) == sp.cancel(D.subs(v, 1)-1) == 0
    # At z=-1 the denominator vanishes only at v=1/2. The derivatives
    # with respect to (v, arg z) fix the sign of the intersection.
    assert sp.diff(infinity_factor, v).subs(z, -1) == -2
    phase_derivative = sp.I*z*sp.diff(infinity_factor, z)
    assert phase_derivative.subs({z: -1, v: sp.Rational(1, 2)}) == -sp.I/2


def rotated_integral(z):
    """Integrate in x analytically, then along positive u in the rotated family."""
    # Compactify the semi-infinite u integral and resolve the near-cut pole.
    def compact(v):
        d = 2*v-1
        log_odds = (np.log1p(2*d/(1-d)) if abs(d) < .25
                    else np.log(v)-np.log1p(-v))
        return log_odds/d if d else 2.

    if z == 1:
        return 4*quad(compact, 0, 1, epsabs=1e-10, epsrel=1e-10)[0]
    pole = 1/(1-z)
    # Cancel the nearby complex pole before quadrature and integrate the
    # subtracted constant analytically on the principal sheet.
    subtract = compact(pole) if .05 < pole.real < .95 and abs(pole.imag)<.1 else 0.
    knots = (0., pole.real, 1.) if subtract else (0., 1.)
    # Absolute tolerances for the two pieces avoid an unnecessary
    # refinement when their real parts nearly cancel after recombination.
    residual = sum(quad(lambda v: (compact(v)-subtract)/(1+(z-1)*v), lo, hi,
                        complex_func=True, epsabs=3e-10, epsrel=1e-10)[0]
                   for lo, hi in zip(knots, knots[1:]))
    return 2*(1+z)*(residual+subtract*np.log(complex(z))/(z-1))


def check_cd_projective_infinity():
    certify_rotation()
    theta = 2*pi*np.arange(256)/256
    for radius in (0.12, 0.025):
        v = .5+radius*np.cos(theta)
        z = np.exp(1j*(pi+radius*np.sin(theta)))
        crossing = 1-v+z*v
        D = z*((1-v)**2+v*v)/crossing
        np.testing.assert_allclose(
            [winding_number(crossing), winding_number(D)],
            [1, -1], atol=1e-12, rtol=0,
        )
    print("One signed infinity crossing: +1 for the chart factor, -1 for D")

    for z in (1., 2., np.exp(0.7j), np.exp(-2.9j), -0.4+0.03j, -2.5-0.03j):
        np.testing.assert_allclose(rotated_integral(z), integrated_projective(z),
                                   atol=1e-9, rtol=0)

    # The rotated u pole approaches the real axis from above when Im(z)>0.
    # Its two bypasses give +2*pi*i times the residue in the normalization
    # of the u integration. On the seam the scalar jump vanishes.
    for a in (0.4, 1., 2.5):
        residue = quad(lambda x: -2*(1-a)/((1+x)*(1+a*x)), 0, np.inf,
                       epsabs=1e-11)[0]
        jump = 2j*pi*residue
        np.testing.assert_allclose(jump, 4j*pi*log(a), atol=1e-10, rtol=0)
        errors = []
        for height in (1e-3, 1e-4, 1e-5, 1e-6):
            above = rotated_integral(complex(-a, height))
            below = rotated_integral(complex(-a, -height))
            errors.append(abs(above-below-jump))
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 4e-5
        print(f"z=-{a:g}: infinity residue jump={jump:.9f}, lateral error={errors[-1]:.3e}")

    print("PASS: periodic projective family and its infinity crossing")


if __name__ == "__main__":
    check_cd_projective_infinity()
