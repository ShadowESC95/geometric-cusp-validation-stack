"""Locate finite-regulator poles on the diagonal complex CD cap."""

from math import pi

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .finite_cross_ratio import frame_parameters
from .full_finite_frame import ROTATION


def certify_diagonal_denominator():
    z = sp.symbols("z", positive=True)
    x, A, s, Q, eps, k, radius, chi = sp.symbols("x A s Q eps k radius chi")
    _, _, (p1, p2, p3, p4, S, T) = frame_parameters(z)
    a1, a2 = p3*x/S, -p4*x/T
    C = p1*a1+p2*a2
    a3 = C*(1-A-s)/(A*(p3+S*a1))
    a4 = C*Q*s/(A*(-p4+T*a2))
    delta = C+(p3+S*a1)*a3-(-p4+T*a2)*a4
    regulator = -sp.I*eps*(1+a1+a2+a3+a4)**2

    H = 1-(1+Q)*s
    G = 1+(Q-1)*s
    U = A*((5*z+3)/(1+z)+4*x)+x*(G-A)/(1+x)
    normalized = H-sp.I*eps*U**2/(4*x*A)
    assert sp.cancel((delta+regulator)*A/C-normalized) == 0

    # The punctured cap has z=-1+radius/k. The normalized denominator
    # is a quadratic in k, with a finite limit on x=eps*chi.
    on_cap = sp.cancel(normalized.subs(z, -1+radius/k))
    U_cap = A*(5-2*k/radius+4*x)+x*(G-A)/(1+x)
    assert sp.cancel(on_cap-(H-sp.I*eps*U_cap**2/(4*x*A))) == 0
    lower = sp.cancel(on_cap.subs(x, eps*chi))
    assert sp.cancel(lower.subs(eps, 0)-(H-sp.I*A*(5-2*k/radius)**2/(4*chi))) == 0


def denominator_on_cap(k, radius, x, A, s, eps):
    H = 1-(1+ROTATION)*s
    G = 1+(ROTATION-1)*s
    U = A*(5-2*k/radius+4*x)+x*(G-A)/(1+x)
    return H-1j*eps*U**2/(4*x*A)


def original_denominator(z, x, A, s, eps):
    """Evaluate the regulated polynomial using the original affine parameters."""
    p1, p2 = 2*z/(5*z+3), 2/(5*z+3)
    p3 = 4*(1+z)/(5*z+3)
    a1 = a2 = p3*x/2
    C = p1*a1+p2*a2
    a3 = C*(1-A-s)/(A*(p3+2*a1))
    a4 = C*ROTATION*s/(A*(p3+2*a2))
    delta = C+(p3+2*a1)*a3-(p3+2*a2)*a4
    regulator = -1j*eps*(1+a1+a2+a3+a4)**2
    return A*(delta+regulator)/C


def check_cap_roots():
    radius, A, s, chi = .06, .5, .25, .25
    H = 1-(1+ROTATION)*s
    G = 1+(ROTATION-1)*s
    offset = radius/2*np.sqrt(-4j*chi*H/A)
    limiting_roots = 2.5*radius+np.array([-offset, offset])
    print(f"Joint limit poles in k: {limiting_roots[0]:.9f}, {limiting_roots[1]:.9f}")
    unit = np.exp(2j*pi*np.arange(256)/256)
    previous = np.inf
    for eps in (.3, .08, .01, 1e-4, 1e-8):
        x = eps*chi
        center = radius/2*(5+4*x+x*(G-A)/(A*(1+x)))
        roots = center+np.array([-offset, offset])

        # Exclude the two base-pole holes and the frame singularity d=0.
        holes = np.array([0., radius, 2.5*radius])
        assert np.max(abs(roots)) < .25
        assert np.min(abs(roots[:, None]-holes)) > .022
        z = -1+radius/roots
        np.testing.assert_allclose(denominator_on_cap(roots, radius, x, A, s, eps),
                                   0., atol=3e-14, rtol=0)
        np.testing.assert_allclose(original_denominator(z, x, A, s, eps),
                                   0., atol=3e-14, rtol=0)

        # The derivative at either root is nonzero. Check the local
        # degree and the degree on the rim, including the punctures.
        derivative = -2j*eps*A/(x*radius**2)*(roots-center)
        step = .001
        divided_difference = (denominator_on_cap(roots+step, radius, x, A, s, eps)
                              -denominator_on_cap(roots-step, radius, x, A, s, eps))/(2*step)
        np.testing.assert_allclose(divided_difference, derivative, rtol=1e-12, atol=1e-12)
        assert np.min(abs(derivative)) > 1
        for root in roots:
            around_root = denominator_on_cap(root+.005*unit, radius, x, A, s, eps)
            np.testing.assert_allclose(winding_number(around_root), 1, atol=1e-12)
        for hole in holes:
            around_hole = denominator_on_cap(hole+.02*unit, radius, x, A, s, eps)
            np.testing.assert_allclose(winding_number(around_hole), 0, atol=1e-12)
        rim = denominator_on_cap(unit, radius, x, A, s, eps)
        np.testing.assert_allclose(winding_number(rim), 2, atol=1e-12)
        error = np.max(abs(roots-limiting_roots))
        assert error < previous
        previous = error
        print(f"eps={eps:g}: pole-limit error={error:.3e}; rim winding=2")
    assert previous < 1e-8


def check_cd_regulated_base_cap():
    certify_diagonal_denominator()
    check_cap_roots()
    print("PASS: two simple regulated poles obstruct the undeformed complex base-cap lift")


if __name__ == "__main__":
    check_cd_regulated_base_cap()
