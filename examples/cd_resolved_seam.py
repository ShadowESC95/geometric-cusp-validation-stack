"""Resolve the normal circle of the polar seam on the projective CD family."""

from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from geocusp.contours import winding_number


def normal_coordinates(radius, angle):
    """Write f=1+u*z=r*exp(i*angle), with u>0 and |z|=1."""
    f = radius*np.exp(1j*angle)
    y = f-1
    u = abs(y)
    z = y/u
    return f, u, z


def disk_contribution(radius, order):
    """Normalized CD period over |1+u*z|<radius in the real sweep."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    distance = radius*(nodes[:, None]+1)/2
    angle = pi*(nodes[None, :]+1)
    f, u, _ = normal_coordinates(distance, angle)
    # u-1=(r^2-2*Re(f))/(u+1) avoids cancellation when r is small.
    difference = (distance**2-2*f.real)/(u+1)
    one_plus_z = (difference+f)/u
    x_weight = np.divide(np.log1p(difference), difference, out=np.ones_like(u),
                         where=difference != 0)
    # du dphi=(r/u)dr dangle. The z and CD factors use dphi/(2*pi).
    integrand = one_plus_z/f*x_weight*distance/(pi*u)
    return radius*pi/2*(weights @ integrand @ weights)


def certify_normal_form():
    x, u, z, f = sp.symbols("x u z f", nonzero=True)
    # At fixed z, du=df/z. Keep the dx,df orientation inherited from dx,du.
    rotated = 2*(1+z)/((1+x)*(x+u)*(1+u*z))
    coefficient = sp.cancel(rotated.subs(u, (f-1)/z)/z)
    residue = sp.cancel(f*coefficient).subs(f, 0)
    assert sp.cancel(residue-2*(1+z)/((1+x)*(z*x-1))) == 0
    assert sp.cancel(residue.subs(z, -1)) == 0
    stripped = sp.cancel(residue/(1+z))
    assert sp.cancel(stripped.subs(z, -1)+2/(1+x)**2) == 0

    r, theta = sp.symbols("r theta", real=True)
    modulus = sp.sqrt(1-2*r*sp.cos(theta)+r*r)
    one_plus_z = 1+(r*sp.exp(sp.I*theta)-1)/modulus
    linear = sp.simplify(sp.diff(one_plus_z, r).subs(r, 0).expand(complex=True))
    assert linear == sp.I*sp.sin(theta)
    angular_mean = sp.integrate(linear*sp.exp(-sp.I*theta), (theta, 0, 2*sp.pi))
    assert sp.simplify(angular_mean-sp.pi) == 0


def certify_cap_obstruction(radius, angles):
    """Track the boundary in the complement of both y poles."""
    f, u, z = normal_coordinates(radius, angles)
    for x in (0.3, 1., 4.):
        # The ratio (1+y)/(zx+y) extends to 1 at projective y=infinity
        # when x and z remain finite. On the normal circle y=z*u.
        ratio = f/(z*(x+u))
        np.testing.assert_allclose(winding_number(ratio), 1, rtol=0, atol=1e-12)

    # At the limiting seam z=-1, fill from the exterior using k=r/f.
    # The point k=0 is y=infinity; the other divisor y=x meets the disk.
    k, r, x = sp.symbols("k r x", positive=True)
    y = r/k-1
    chart_ratio = sp.cancel((1+y)/(y-x))
    assert chart_ratio.subs(k, 0) == 1
    assert sp.cancel((chart_ratio*(r-(1+x)*k))-r) == 0
    assert sp.solve(sp.cancel(k*(y-x)), k) == [r/(1+x)]


def check_cd_resolved_seam():
    certify_normal_form()
    angles = 2*pi*np.arange(256)/256
    previous = np.inf
    for radius in (0.2, 0.08, 0.025, 0.008):
        f, u, z = normal_coordinates(radius, angles)
        np.testing.assert_allclose(1+u*z, f, rtol=0, atol=5e-16)
        np.testing.assert_allclose(
            [winding_number(f), winding_number(z)], [1, 0], atol=1e-12, rtol=0,
        )
        value = disk_contribution(radius, 32)
        np.testing.assert_allclose(disk_contribution(radius, 48), value, rtol=0, atol=1e-10)
        error = abs(value/radius**2-0.5)
        assert error < previous
        previous = error
        print(f"normal radius={radius:g}: winding=+1, "
              f"disk contribution/radius^2={value.real/radius**2:.12f}")
    assert previous < 1e-4
    certify_cap_obstruction(radius, angles)
    print("Cap obstruction at finite x,z: winding +1 even through y=infinity")

    # A counterclockwise normal meridian at the complexified seam.
    # Removing (1+z) leaves the coefficient of dx df/f at z=-1.
    stripped_period = 2j*pi*quad(lambda x: -2/(1+x)**2, 0, np.inf)[0]
    np.testing.assert_allclose(stripped_period, -4j*pi, rtol=0, atol=1e-12)
    print(f"Seam normal residue period divided by (1+z): {stripped_period:.9f}")
    print("PASS: resolved CD seam, normal winding and disk contribution")


if __name__ == "__main__":
    check_cd_resolved_seam()
