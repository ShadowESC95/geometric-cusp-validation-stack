"""Compare the resolved CD frame with the rotated Feynman parameter chart."""

from math import pi, sin

import numpy as np
import sympy as sp

from .cd_regulated_base_cap import denominator_on_cap
from .cd_regulated_cap_meridians import RADIUS
from .full_finite_frame import ROTATION


GAP = sin(pi/8)


def certify_transverse_chart():
    x, t, w, A, s, Q, eps = sp.symbols("x t w A s Q eps")
    l = x*w/A
    z = -1+1/(sp.Rational(5, 2)+l)
    p1, p2, p3 = 1+3/(2*l), -1-5/(2*l), -2/l
    a1, a2 = p3*x/2, p3*t*x/2
    C = p1*a1+p2*a2
    B = 1-A-s
    a3 = C*B/(A*(p3+2*a1))
    a4 = C*Q*s/(A*(p3+2*a2))
    J = (3-5*t)/2+(1-t)*l
    assert sp.cancel(J+(sp.Rational(5, 2)+l)*(z+t)) == 0
    assert sp.cancel(C+x*J/l**2) == 0
    image = w*sp.Matrix([a1, a2, a3, a4, 1])
    expected = sp.Matrix([-A, -A*t, B*J/(2*(1+x)), Q*s*J/(2*(1+t*x)), w])
    assert all(sp.cancel(entry) == 0 for entry in image-expected)
    H = 1-(1+Q)*s
    denominator = H-sp.I*eps*A*(1+a1+a2+a3+a4)**2/C
    assert sp.cancel(denominator-(H+sp.I*eps*x*sum(expected)**2/(A*J))) == 0

    # At fixed z, the ray parameters have positive orientation in
    # (x,t,r,s), and replacing r by 1-A-s reverses the fiber orientation.
    fiber_coordinates = sp.Matrix([(1-A-s)/A, s/A])
    assert sp.cancel(fiber_coordinates.jacobian([A, s]).det()+1/A**3) == 0

    # r=1-A-s. With (x,t,w,A,s) ordered as written, the two
    # negative Jacobian signs from z and r cancel.
    jacobian = sp.Matrix([x, t*x, z, B, s]).jacobian([x, t, w, A, s]).det()
    assert sp.cancel(jacobian-x*x/(A*(sp.Rational(5, 2)+l)**2)) == 0
    cd_density = 2*(1+z)/(z*(1+x)*(1+t*x)*x*(z+t))
    K = sp.Rational(5, 2)+l
    reduced_factor = -2/((1+x)*(1+t*x)*K*(1-K)*J)
    assert sp.cancel(cd_density*jacobian-x*reduced_factor/A) == 0


def frame_chart(x, t, w, A, s, eps):
    l = x*w/A
    z = -1+1/(2.5+l)
    J = (3-5*t)/2+(1-t)*l
    B = 1-A-s
    params = (-A, -A*t, B*J/(2*(1+x)), ROTATION*s*J/(2*(1+t*x)), w)
    H = 1-(1+ROTATION)*s
    denominator = H+1j*eps*x*sum(params)**2/(A*J)
    return z, J, params, denominator


def check_real_overlap():
    # Continue the same real parameter ray across the base point z=0.
    # The unreduced affine polynomial is regular there on this slice.
    for z in (2., 0., -.2, -.59, -.59999):
        l = 1/(1+z)-2.5
        x, s, A = .02, .25, .4
        w = A*l/x
        for t in (.8, 1., 1.2):
            actual_z, J, params, normalized = frame_chart(x, t, w, A, s, .08)
            np.testing.assert_allclose(actual_z, z, atol=3e-15, rtol=0)
            assert w < 0 and J < 0
            a1, a2, a3, a4 = (value/w for value in params[:4])
            assert min(a1, a2, a3, (a4/ROTATION).real) > 0
            np.testing.assert_allclose((a4/ROTATION).imag, 0., atol=1e-12, rtol=0)
            p1, p2, p3 = 1+3/(2*l), -1-5/(2*l), -2/l
            C = p1*a1+p2*a2
            delta = 2*(a1*a3-a2*a4)+C+p3*(a3-a4)
            np.testing.assert_allclose(normalized, A*(delta-.08j*(1+a1+a2+a3+a4)**2)/C,
                                       atol=2e-12, rtol=2e-14)
            projected = (np.exp(-1j*pi/8)*normalized).imag
            assert projected <= -GAP+1e-13
    print("Negative real frame chart agrees with the rotated Feynman ray for 0.8 <= y/x <= 1.2")


def check_transverse_density():
    # |l| <= .03, .494 <= |J| <= 1.506, and |sum(params)| <= 6.4.
    # Thus |D-H| <= .1*.01*6.4^2/.494 < GAP/2. The density
    # bound uses |2.5+l| >= 2.47 and |1-(2.5+l)| >= 1.47.
    assert .1*.01*6.4**2/.494 < GAP/2
    assert 4/(2.47*1.47*.494*(GAP/2)**3) < 320
    theta = np.linspace(0, 2*pi, 129)
    w = 3*np.exp(1j*theta)
    for s in (0., .25, .5857864376, .95):
        for fraction in (.05, .8):
            A = fraction*(1-s)
            for t in (.8, 1., 1.2):
                previous = np.inf
                for ratio in (.01, .001, .0001):
                    x = A*ratio
                    z, J, params, D = frame_chart(x, t, w, A, s, .1)
                    assert np.max(abs(sum(params))) <= 6.4
                    assert np.min(abs(J)) >= .494-1e-15
                    assert np.min(abs(D)) >= GAP/2
                    if t == 1.:
                        k = RADIUS*(2.5+ratio*w)
                        np.testing.assert_allclose(D, denominator_on_cap(k, RADIUS, x, A, s, .1),
                                                   atol=1e-14, rtol=0)
                    K = 1/(1+z)
                    density = 4*ROTATION*ratio/((1+x)*(1+t*x)*K*(1-K)*J*D**3)
                    magnitude = np.max(abs(density))
                    assert magnitude <= 320*ratio
                    assert magnitude < previous
                    previous = magnitude
    print("Transverse five-form density is bounded by 320*x/A across the bounded frame chart")


def check_cd_frame_feynman_overlap():
    certify_transverse_chart()
    check_real_overlap()
    check_transverse_density()
    print("PASS: local rotated Feynman ray overlap and transverse frame bound")


if __name__ == "__main__":
    check_cd_frame_feynman_overlap()
