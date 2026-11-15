"""Compactify the zero-regulator auxiliary fiber and resolve two map corners."""

from math import pi, sin

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .complex_cross_ratio import check_projective_pullback
from .finite_cross_ratio import frame_parameters
from .full_finite_frame import ROTATION


def certify_fiber():
    """Use (r,s)=(w*(1-t),t) on the normalized projective triangle."""
    r, s, w, t, k, Q = sp.symbols("r s w t k Q", nonzero=True)
    ray_weight = -2*Q/(1-Q*k)**3
    compactified = sp.cancel(ray_weight.subs(k, t/(1-t))/(1-t)**2)
    density = -2*Q/(1-(1+Q)*s)**3
    assert sp.cancel(compactified-(1-t)*density.subs(s, t)) == 0
    assert sp.Matrix([w*(1-t), t]).jacobian([w, t]).det() == 1-t
    A, B, C = 1-r-s, r, Q*s
    assert sp.expand(A+B-C-(1-(1+Q)*s)) == 0

    primitive = -(1-s)**2/(1-(1+Q)*s)**2
    assert sp.cancel(sp.diff(primitive, s)-(1-s)*density) == 0
    assert primitive.subs(s, 1)-primitive.subs(s, 0) == 1

    # Counterclockwise boundary in dr^ds orientation: s=0, then A=0,
    # then r=0. Each edge pulls back the fiber two-form to zero.
    for edge, start, end in (([k, 0], [0, 0], [1, 0]),
                             ([1-k, k], [1, 0], [0, 1]),
                             ([0, 1-k], [0, 1], [0, 0])):
        curve = sp.Matrix(edge)
        assert curve.subs(k, 0) == sp.Matrix(start)
        assert curve.subs(k, 1) == sp.Matrix(end)


def certify_diagonal_map():
    """Return symbolic limits in the charts A=x*b and A=X*b, X=1/x."""
    z = sp.symbols("z", positive=True)
    x, X, A, s, b, Q = sp.symbols("x X A s b Q")
    _, _, (p1, p2, p3, p4, S, T) = frame_parameters(z)
    a1, a2 = p3*x/S, -p4*x/T
    C = p1*a1+p2*a2
    B = 1-s-A
    a3 = C*B/(A*(p3+S*a1))
    a4 = C*Q*s/(A*(-p4+T*a2))
    d, h, F = 5*z+3, z+1, 1+x
    original = sp.Matrix([a1, a2, a3, a4, 1])
    homogeneous = sp.Matrix([2*h*x*F*A, 2*h*x*F*A, h*x*B, h*x*Q*s, d*F*A])
    assert all(sp.cancel(value) == 0 for value in original*d*F*A-homogeneous)

    lower = homogeneous.subs(A, x*b).applyfunc(lambda v: sp.cancel(v/x))
    upper = homogeneous.subs({A: X*b, x: 1/X}).applyfunc(lambda v: sp.cancel(X*v))
    expected_lower = sp.Matrix([0, 0, h*(1-s), h*Q*s, d*b])
    expected_upper = sp.Matrix([2*h*b, 2*h*b, h*(1-s), h*Q*s, 0])
    for actual, expected in ((lower.subs(x, 0), expected_lower),
                             (upper.subs(X, 0), expected_upper)):
        assert all(sp.cancel(value) == 0 for value in actual-expected)

    # Opposite affine charts x=A*c and X=A*c cover b=infinity.
    c = sp.symbols("c")
    lower_other = homogeneous.subs(x, A*c).applyfunc(lambda v: sp.cancel(v/A))
    upper_other = homogeneous.subs(x, 1/(A*c)).applyfunc(lambda v: sp.cancel(A*c*c*v))
    assert all(sp.cancel(value) == 0 for value in
               lower_other.subs(A, 0)-sp.Matrix([0, 0, c*h*(1-s), c*h*Q*s, d]))
    assert all(sp.cancel(value) == 0 for value in
               upper_other.subs(A, 0)-sp.Matrix([2*h, 2*h, c*h*(1-s), c*h*Q*s, 0]))


def check_auxiliary_triangle():
    Q = ROTATION

    def density(s):
        return -2*Q/(1-(1+Q)*s)**3

    values = np.linspace(0, 1, 1025)
    gap = sin(pi/8)
    np.testing.assert_allclose(np.min(abs(1-(1+Q)*values)), gap, atol=1e-14, rtol=0)
    full = quad(lambda s: (1-s)*density(s), 0, 1, complex_func=True, epsabs=1e-11)[0]
    np.testing.assert_allclose(full, 1, atol=1e-12, rtol=0)
    print(f"Auxiliary triangle period={full:.12f}; denominator gap={gap:.9f}")

    for width in (.04, .004, .0004):
        retained = quad(lambda s: (1-s-3*width)*density(s), width, 1-2*width,
                        complex_func=True, epsabs=1e-11)[0]
        # Retain r,s,A >= width; the discarded area is O(width).
        area = (1-(1-3*width)**2)/2
        assert abs(full-retained) <= 2*area/gap**3
        print(f"edge trim={width:g}: omitted auxiliary period={abs(full-retained):.3e}")


def check_corner_limits():
    """Compare finite points of the diagonal cap with its new ratio faces."""
    s = .25
    for ratio in (.25, 2.):
        k = np.array([.5j, .5, -.2])
        z = -1+.06/k
        h, d = 1+z, 5*z+3
        zeros = np.zeros_like(z)
        expected = np.array([zeros, zeros, h*(1-s), h*ROTATION*s, d*ratio])
        expected_upper = np.array([2*h*ratio, 2*h*ratio, h*(1-s), h*ROTATION*s, zeros])
        previous = np.inf
        for distance in (1e-2, 1e-4, 1e-6):
            A = distance*ratio
            B = 1-s-A

            def coordinates(x):
                leading = 2*h*x*(1+x)*A
                return np.array([leading, leading, h*x*B, h*x*ROTATION*s, d*(1+x)*A])

            low_error = np.max(abs(coordinates(distance)/distance-expected))
            high_error = np.max(abs(distance*coordinates(1/distance)-expected_upper))
            error = max(low_error, high_error)
            assert error < previous
            previous = error
        assert previous < 3e-5
        print(f"joint corner A/distance={ratio:g}: chart-limit error={previous:.3e}")


def check_cd_auxiliary_faces():
    check_projective_pullback()
    certify_fiber()
    certify_diagonal_map()
    check_auxiliary_triangle()
    check_corner_limits()
    print("PASS: normalized auxiliary fiber and the two diagonal-cap corner charts")


if __name__ == "__main__":
    check_cd_auxiliary_faces()
