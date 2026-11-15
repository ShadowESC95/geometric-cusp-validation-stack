"""Check the regulated auxiliary infinity edge on the z=1 diagonal."""

from math import pi, sin

import numpy as np
import sympy as sp
from scipy.integrate import quad_vec

from .finite_cross_ratio import frame_parameters
from .full_finite_frame import ROTATION


def certify_regulated_edge():
    x, A, s, Q, b = sp.symbols("x A s Q b", nonzero=True)
    _, _, (p1, p2, p3, p4, S, T) = frame_parameters(sp.Integer(1))
    a1, a2 = p3*x/S, -p4*x/T
    C = p1*a1+p2*a2
    assert sp.cancel(C-x/4) == 0
    a3 = C*(1-A-s)/(A*(p3+S*a1))
    a4 = C*Q*s/(A*(-p4+T*a2))
    projective_sum = 1+a1+a2+a3+a4
    v = x/(1+x)**2
    G = 1+(Q-1)*s
    assert sp.cancel(A*projective_sum-(1+x)*(A+v*(G-A)/4)) == 0

    # The normalized regulator survives in the joint corner A=v*b.
    multiplier = sp.cancel((A*projective_sum**2/C).subs(A, v*b))
    corner = 4/b*(b+G/4)**2
    assert sp.cancel(multiplier-4/b*(b+(G-v*b)/4)**2) == 0
    assert sp.cancel(multiplier.subs(x, 0)-corner) == 0
    assert sp.cancel(sp.limit(multiplier, x, sp.oo)-corner) == 0
    assert sp.cancel(multiplier.subs(x, 1/x)-multiplier) == 0

    # The edge coefficient is integrated before taking the regulator limit.
    H = 1-(1+Q)*s
    primitive = sp.integrate(G**2/H**4, s)
    exact = sp.cancel(primitive.subs(s, 1)-primitive.subs(s, 0))
    assert sp.cancel(exact+1/(3*Q)) == 0


def normalized_denominator(s, A, x, eps):
    v = x/(1+x)**2
    H = 1-(1+ROTATION)*s
    G = 1+(ROTATION-1)*s
    multiplier = 4*(A+v*(G-A)/4)**2/(v*A)
    return H-1j*eps*multiplier


def check_gap_and_corner_limits():
    s = np.linspace(0, 1, 129)[1:-1]
    A = (1-s)*np.geomspace(1e-10, 1., 81)[:, None]
    gap = sin(pi/8)
    for x in (.03, 1., 14.):
        for eps in (1e-6, .02, .3):
            denominator = normalized_denominator(s, A, x, eps)
            projected = np.imag(np.exp(-1j*pi/8)*denominator)
            assert np.max(projected) <= -gap+1e-14
    print(f"Uniform regulated denominator gap >= {gap:.9f}")

    eps = .08
    G = 1+(ROTATION-1)*s
    H = 1-(1+ROTATION)*s
    for b in (.25, 2.):
        target = H-4j*eps/b*(b+G/4)**2
        previous = np.inf
        for x in (.01, .0001, .000001):
            A = b*x/(1+x)**2
            lower = normalized_denominator(s, A, x, eps)
            upper = normalized_denominator(s, A, 1/x, eps)
            np.testing.assert_allclose(lower, upper, atol=1e-14, rtol=0)
            error = np.max(abs(lower-target))
            assert error < previous
            previous = error
        assert previous < 2e-6
        assert np.max(abs(target-H)) > .05
        print(f"joint ratio={b:g}: two-endpoint regulated limit error={previous:.3e}")


def fiber_correction(x, eps, order=64):
    """Integrate the regulated minus unregulated fiber, resolving A near zero."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    # Set A=(1-s)*t and reuse the s nodes for vector-valued t quadrature.
    s = (1+nodes)/2
    H = 1-(1+ROTATION)*s

    def correction(t):
        D = normalized_denominator(s, (1-s)*t, x, eps)
        # Difference of inverse cubes, with the regulator shift kept explicit.
        diff = ((H-D)/H)*(1+H/D+(H/D)**2)/(H**2*D)
        return -2*ROTATION*(1-s)*diff

    v = x/(1+x)**2
    split = eps*v/4
    integral, _ = quad_vec(correction, 0, 1, points=[split, np.sqrt(split)],
                           epsabs=2e-13, epsrel=2e-10, limit=300)
    return weights @ integral/2


def check_edge_logarithm():
    for x in (.35, 1.):
        target = .5j*x/(1+x)**2
        scaled = []
        for eps in (1e-4, 1e-6, 1e-8, 1e-10):
            diff = fiber_correction(x, eps)
            scaled.append(diff/eps)
        # Divided differences in log(1/eps) remove the O(eps) term.
        slopes = np.diff(scaled)/np.log(100.)
        errors = abs(slopes-target)
        assert np.all(np.diff(errors) < 0)
        assert errors[-1] < 5e-6
        print(f"x={x:g}: edge log coefficient={slopes[-1]:.9f}, target={target:.9f}")


def check_cd_regulated_auxiliary():
    certify_regulated_edge()
    check_gap_and_corner_limits()
    check_edge_logarithm()
    print("PASS: regulated auxiliary gap, joint corners and edge logarithm at z=1")


if __name__ == "__main__":
    check_cd_regulated_auxiliary()
