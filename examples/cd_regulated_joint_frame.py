"""Resolve the simultaneous regulator and frame limit of the diagonal cap."""

from math import pi, sin

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_regulated_cap_meridians import RADIUS
from .cd_regulated_meridian_sweep import root_direction
from .full_finite_frame import ROTATION


def certify_joint_frame_charts():
    x, A, s, Q, R, k, w, v, scale = sp.symbols("x A s Q R k w v scale")
    B = 1-A-s
    H = 1-(1+Q)*s
    G = 1+(Q-1)*s
    # The zero-regulator affine formulas for the coordinate map are also
    # the coordinates used before inserting the finite regulator.
    z = -1+R/k
    p3 = sp.cancel(4*(1+z)/(5*z+3))
    a1 = p3*x/2
    C = p3*a1/2
    a3 = C*B/(A*(p3+2*a1))
    a4 = C*Q*s/(A*(p3+2*a1))
    image = sp.Matrix([a1, a1, a3, a4, 1])
    # w = u/(eps*tau) = A*(k/R-5/2)/x. This expression removes
    # the regulator from the chart of the projective map.
    resolved = image.subs(k, R*(sp.Rational(5, 2)+x*w/A))
    resolved = resolved.applyfunc(lambda entry: sp.cancel(w*entry))
    expected = sp.Matrix([-A, -A, -B/(2*(1+x)), -Q*s/(2*(1+x)), w])
    assert all(sp.cancel(entry) == 0 for entry in resolved-expected)

    # The reciprocal chart v=1/w has last coordinate one and extends
    # across v=0. The limiting image recovers A, s and w on L=1.
    reciprocal = resolved.subs(w, 1/v).applyfunc(lambda entry: sp.cancel(v*entry))
    assert reciprocal.subs(v, 0) == sp.Matrix([0, 0, 0, 0, 1])
    limit = resolved.subs(x, 0)
    section = -limit[0]-2*limit[2]-2*limit[3]/Q
    assert sp.cancel(section-1) == 0
    assert sp.cancel(-2*limit[3]/Q-s) == 0 and limit[0] == -A and limit[4] == w

    U = A*(5-2*k/R+4*x)+x*(G-A)/(1+x)
    denominator = H-sp.I*scale*U**2/(4*x**2)
    # scale = (eps*tau)^2 = eps*x/A. Keeping w bounded gives a
    # uniform perturbation of H, including at w=0.
    pulled_back = denominator.subs(k, R*(sp.Rational(5, 2)+x*w/A))
    shift = 2*A+(G-A)/(2*(1+x))
    assert sp.cancel(pulled_back-(H-sp.I*scale*(w-shift)**2)) == 0


def check_split_scales():
    radius = 3.
    unit = np.exp(2j*pi*np.arange(256)/256)
    w = radius*unit
    gap = sin(pi/8)
    for s in (0., .25, 1/(1+ROTATION.real), .95, .999):
        H = 1-(1+ROTATION)*s
        G = 1+(ROTATION-1)*s
        for A in (.02*(1-s), .8*(1-s)):
            B = 1-s-A
            previous = np.inf
            for eps, tau in ((.08, .1), (.008, .08), (.0001, .04), (1e-6, .02)):
                x = eps*A*tau**2
                scale = eps*tau
                shift = 2*A+(G-A)/(2*(1+x))
                resolved = H-1j*scale**2*(w-shift)**2
                assert np.max(abs(resolved-H)) <= scale**2*(radius+2)**2
                assert np.min(abs(resolved)) >= gap/2
                np.testing.assert_allclose(winding_number(resolved), 0., atol=1e-12, rtol=0)
                # Restore the affine denominator only on the punctured loop.
                affine = resolved/(eps*tau**2*w**2)
                np.testing.assert_allclose(winding_number(affine), -2., atol=1e-12, rtol=0)

                u_roots = scale*shift+np.array([1, -1])*root_direction(s)/RADIUS
                assert np.min(abs(u_roots/scale)) > radius
                v_roots = scale/u_roots
                np.testing.assert_allclose(H-1j*(scale/v_roots-scale*shift)**2,
                                           0., atol=8e-16, rtol=0)
                # Resolve the approach to v=0 at the pole scale.
                reciprocal_targets = RADIUS/root_direction(s)*np.array([1, -1])
                assert np.max(abs(v_roots/scale-reciprocal_targets)) <= 6*scale
                error = np.max(abs(v_roots))
                assert error < previous
                previous = error
                # In reciprocal projective coordinates, each root maps
                # toward the same point [0:0:0:0:1]. Retain its u label.
                direction = np.array([-A, -A, -B/(2*(1+x)), -ROTATION*s/(2*(1+x))])
                coordinates = direction[:, None]*v_roots
                assert np.max(abs(coordinates)) <= error
            assert previous < 4e-8
    print("Bounded frame chart: normalized winding 0, affine winding -2")
    print("Both regulated roots approach the reciprocal projective point, with distinct rescaled positions")


def check_cd_regulated_joint_frame():
    certify_joint_frame_charts()
    check_split_scales()
    print("PASS: simultaneous regulator and frame limit of the homogeneous map")


if __name__ == "__main__":
    check_cd_regulated_joint_frame()
