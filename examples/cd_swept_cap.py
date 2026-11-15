"""Sweep the compactified CD caps and account for their oriented boundaries."""

from math import pi

import numpy as np
import sympy as sp

from .cd_resolved_seam import disk_contribution, normal_coordinates


def certify_cap_pullbacks():
    """Check the boundary orientations and the zero-period cap pieces."""
    x, y, z, q, s, t, f, base, delta, k = sp.symbols("x y z q s t f base delta k")
    kernel = 2*(1+z)/(z*(1+x)*(1+y)*(z*x+y))

    # For the first collar, columns give derivatives in (theta, s).
    # y=f-1, df/dtheta=i*f, and z=(1-s)*base-s.
    jacobian = -sp.I*f*(1+base)
    first = sp.cancel(kernel.subs({y: f-1, z: (1-s)*base-s})*jacobian/(2*sp.pi*sp.I))
    compact = sp.cancel(first.subs(x, t/(1-t))/(1-t)**2)
    expected = -(1-s)*(1+base)**2/(sp.pi*z*(z*t+(f-1)*(1-t)))
    assert sp.cancel(compact-expected.subs(z, (1-s)*base-s)) == 0

    # The second collar is on z=-1, hence dz pulls back to zero.
    assert sp.Matrix([x, y, -1]).jacobian([x, y, s]).det() == 0

    # The third collar uses (x,q,s). Its pulled-back coefficient is
    # holomorphic near q=0 even though the moving section vanishes there.
    moved_y, moved_z = x+(1-s)*q, -1+s*q/x
    jacobian = sp.Matrix([x, moved_y, moved_z]).jacobian([x, q, s]).det()
    assert sp.cancel(jacobian-q/x) == 0
    third = sp.cancel(kernel.subs({y: moved_y, z: moved_z})*jacobian)
    expected = 2*s*q/(x**2*(1+x)*moved_z*(1+moved_y))
    assert sp.cancel(third-expected) == 0
    assert sp.denom(third).subs(q, 0) != 0

    # For the last disk and its puncture tubes, y duplicates x. The same
    # holds in the reciprocal z chart, since that alters only the last row.
    assert sp.Matrix([x, x, z]).jacobian([x, z, s]).det() == 0

    # Fiber incidence in the basis seam, three joins, and two pole circles.
    # The transverse disk has two positively oriented meridians removed.
    seam, join1, join2, join3, zero, infinity = sp.eye(6).columnspace()
    boundaries = (seam-join1, join1-join2, join2-join3, join3-zero-infinity)
    fiber_boundary = sum(boundaries, sp.zeros(6, 1))
    assert fiber_boundary == seam-zero-infinity
    # Product boundary: d(I x cap) = dI x cap - I x d(cap).
    assert -fiber_boundary == -seam+zero+infinity

    # In the final disk x=y=t/(1-t). The moving section vanishes on both
    # compactified endpoint caps, independently of k.
    moving_section = (delta-k)*t*(1-t)+k*(1-t)*t
    assert sp.cancel(moving_section-delta*t*(1-t)) == 0
    assert moving_section.subs(t, 0) == moving_section.subs(t, 1) == 0


def cap_period(radius, lower, upper, order=24):
    """Integrate the first collar in the coordinate order (t,theta,s)."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = lower+(upper-lower)*(nodes[:, None, None]+1)/2
    theta = pi*(nodes[None, :, None]+1)
    s = (nodes[None, None, :]+1)/2
    normal, _, base = normal_coordinates(radius, theta)
    z = (1-s)*base-s
    integrand = -(1-s)*(1+base)**2/(pi*z*(z*t+(normal-1)*(1-t)))
    weighted = np.einsum("i,j,k,ijk->", weights, weights, weights, integrand)
    return (upper-lower)*pi/4*weighted


def excised_period(radius, lower, upper, order=28):
    """Independently integrate the excised disk in (t,r,theta)."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = lower+(upper-lower)*(nodes[:, None, None]+1)/2
    r = radius*(nodes[None, :, None]+1)/2
    theta = pi*(nodes[None, None, :]+1)
    normal, u, _ = normal_coordinates(r, theta)
    u_minus_one = (r*r-2*normal.real)/(u+1)
    one_plus_z = (u_minus_one+normal)/u
    integrand = one_plus_z/normal*r/(pi*u*(t+u*(1-t)))
    weighted = np.einsum("i,j,k,ijk->", weights, weights, weights, integrand)
    return (upper-lower)*radius*pi/4*weighted


def check_cd_swept_cap():
    certify_cap_pullbacks()
    intervals = ((.2, .7), (.05, .95), (.002, .998), (0., 1.))
    previous = np.inf
    for radius in (.2, .06, .015, .004):
        for lower, upper in intervals:
            cap = cap_period(radius, lower, upper)
            disk = excised_period(radius, lower, upper)
            np.testing.assert_allclose(cap, disk, rtol=1e-10, atol=2e-13)
            # For R<=.2 the collar pointwise bound is
            # 25*(1-s)*R^2/pi in t. Integrate over either omitted strip.
            np.testing.assert_allclose(cap_period(radius, 0., 1.),
                                       cap+cap_period(radius, 0., lower)
                                       +cap_period(radius, upper, 1.), atol=2e-13, rtol=0)
            omitted = abs(cap_period(radius, 0., 1.)-cap)
            assert omitted <= 25*radius**2*(lower+1-upper)+2e-13

        np.testing.assert_allclose(cap, disk_contribution(radius, 48), atol=2e-13, rtol=0)
        error = abs(cap/radius**2-.5)
        assert error < previous
        previous = error
        print(f"normal radius={radius:g}: swept cap/R^2={cap.real/radius**2:.12f}")

    assert previous < 1e-5
    print("PASS: oriented swept cap, zero-period base tubes and scalar endpoint limit")


if __name__ == "__main__":
    check_cd_swept_cap()
