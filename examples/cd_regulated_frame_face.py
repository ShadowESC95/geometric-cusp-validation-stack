"""Resolve the regulated cap roots and the affine pole at the frame puncture."""

from math import pi

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_regulated_base_cap import denominator_on_cap, original_denominator
from .cd_regulated_cap_meridians import RADIUS
from .full_finite_frame import ROTATION


def certify_frame_pullback():
    eps, A, s, tau, u, Q, w, v = sp.symbols("eps A s tau u Q w v")
    x = eps*A*tau**2
    B = 1-A-s
    G = 1+(Q-1)*s
    H = 1-(1+Q)*s
    p1 = 1+3/(2*tau*u)
    p2 = -1-5/(2*tau*u)
    p3 = -2/(tau*u)
    p4 = -p3
    a1, a2 = -eps*A*tau, -eps*A*tau
    a3 = -eps*tau*B/(2*(1+x))
    a4 = -eps*tau*Q*s/(2*(1+x))
    a5 = u
    homogeneous_sum = a1+a2+a3+a4+a5
    shift = eps*tau*(4*A+(G-A)/(1+x))/2
    assert sp.cancel(homogeneous_sum-(u-shift)) == 0

    delta = 2*a1*a3-2*a2*a4+a5*(p1*a1+p2*a2+p3*a3+p4*a4)
    assert sp.cancel(delta-eps*H) == 0
    assert sp.cancel(delta-sp.I*eps*homogeneous_sum**2
                     -eps*(H-sp.I*(u-shift)**2)) == 0

    # The homogeneous map has another ratio face where tau=u=0.
    lifted = sp.Matrix([a1, a2, a3, a4, a5]).subs(u, tau*w)
    lifted = lifted.applyfunc(lambda entry: sp.cancel(entry/tau).subs(tau, 0))
    expected = sp.Matrix([-eps*A, -eps*A, -eps*B/2, -eps*Q*s/2, w])
    assert all(sp.cancel(entry) == 0 for entry in lifted-expected)
    reciprocal = sp.Matrix([a1, a2, a3, a4, a5]).subs(tau, u*v)
    reciprocal = reciprocal.applyfunc(lambda entry: sp.cancel(entry/u).subs(u, 0))
    assert all(sp.cancel(entry) == 0 for entry in reciprocal-v*expected.subs(w, 1/v))
    # In this second chart the normalized denominator tends to H.
    assert (H-sp.I*(tau*w-shift)**2).subs(tau, 0) == H


def check_winding_at_frame():
    unit = np.exp(2j*pi*np.arange(256)/256)
    eps = .08
    for s in (0., .25, 1/(1+ROTATION.real), .9, .99):
        A = .6*(1-s)
        G = 1+(ROTATION-1)*s
        H = 1-(1+ROTATION)*s
        offsets = -1j*np.sqrt(1j*H)*np.array([1, -1])
        for tau in (.1, .01, .001):
            x = eps*A*tau**2
            shift = eps*tau*(4*A+(G-A)/(1+x))/2
            roots = shift+offsets
            assert np.min(abs(roots)) > .5 and np.max(abs(roots)) < 1.1
            loops = np.array([1.5*unit, .25*unit, roots[0]+.15*unit, roots[1]+.15*unit])
            k = RADIUS*(2.5+tau*loops)
            np.testing.assert_allclose(H-1j*(roots-shift)**2, 0., atol=6e-16, rtol=0)
            normalized = H-1j*(loops-shift)**2
            np.testing.assert_allclose(denominator_on_cap(k, RADIUS, x, A, s, eps),
                                       normalized, atol=8e-12, rtol=0)
            np.testing.assert_allclose(original_denominator(-1+RADIUS/k, x, A, s, eps),
                                       normalized, atol=8e-12, rtol=0)

            affine = eps*normalized/loops**2
            normalized_indices = [winding_number(loop) for loop in normalized]
            affine_indices = [winding_number(loop) for loop in affine]
            np.testing.assert_allclose(normalized_indices, [2, 0, 1, 1], atol=1e-12, rtol=0)
            np.testing.assert_allclose(affine_indices, [0, -2, 1, 1], atol=1e-12, rtol=0)
            for indices in (normalized_indices, affine_indices):
                np.testing.assert_allclose(indices[0]-sum(indices[1:]), 0., atol=1e-12, rtol=0)
    print("Normalized winding (outer, frame, zeros): (2, 0, 1, 1)")
    print("Affine winding (outer, frame, zeros): (0, -2, 1, 1)")


def check_cd_regulated_frame_face():
    certify_frame_pullback()
    check_winding_at_frame()
    print("PASS: resolved frame face and regulated meridian balance in both gauges")


if __name__ == "__main__":
    check_cd_regulated_frame_face()
