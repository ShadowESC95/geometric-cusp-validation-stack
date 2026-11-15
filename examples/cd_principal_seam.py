"""Match the reduced principal contour's seam to the compactified CD cap."""

from math import pi

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_resolved_seam import disk_contribution, normal_coordinates


def principal_normal_coordinates(radius, angle):
    """Resolve h=y/x+z; the positive disk orientation has clockwise h."""
    h = radius*np.exp(-1j*angle)
    cosine = np.sqrt(1-h.imag**2)
    z = -cosine+1j*h.imag
    u = h.real+cosine
    return h, z, u, cosine


def certify_principal_pullbacks():
    x, y, z, u, h, s, base, b, a = sp.symbols("x y z u h s base b a")
    kernel = 2*(1+z)/(z*(1+x)*(1+y)*(z*x+y))
    in_ratio = sp.cancel(x*kernel.subs(y, x*u))
    assert sp.cancel(in_ratio-2*(1+z)/(z*(1+x)*(1+x*u)*(z+u))) == 0

    # The rotated pullback uses y=z*u. Invert both x and z, so the
    # two orientation reversals cancel and the forms agree exactly.
    rotated = sp.cancel(z*kernel.subs(y, z*u))
    corrected = rotated.subs({x: 1/x, z: 1/z}, simultaneous=True)/(x*z)**2
    assert sp.cancel(corrected-in_ratio) == 0
    assert sp.cancel(1+u/z-(u+z)/z) == 0

    # With h=a+i*b, u=a+sqrt(1-b^2) and phi=pi-arcsin(b).
    coordinates = sp.Matrix([a+sp.sqrt(1-b*b), sp.pi-sp.asin(b)])
    assert sp.simplify(coordinates.jacobian([a, b]).det()+1/sp.sqrt(1-b*b)) == 0

    # On the first collar, y=x*(h-z_s) keeps the moving factor equal to x*h.
    moved = (1-s)*base-s
    jacobian = sp.Matrix([x, x*(h-moved), moved]).jacobian([x, h, s]).det()
    assert sp.expand(jacobian+x*(1+base)) == 0
    angular = -sp.I*h
    actual = sp.cancel(kernel.subs({y: x*(h-moved), z: moved})*jacobian*angular/(2*sp.pi*sp.I))
    expected = (1+moved)*(1+base)/(sp.pi*moved*(1+x)*(1+x*(h-moved)))
    assert sp.cancel(actual-expected) == 0


def ratio_weight(offset):
    """Integrate dx/((1+x)*(1+(1+offset)*x)) over the positive half-line."""
    return np.divide(np.log1p(offset), offset, out=np.ones_like(offset), where=offset != 0)


def principal_disk_period(radius, order=40):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    r = radius*(nodes[:, None]+1)/2
    angle = pi*(nodes[None, :]+1)
    h, _, _, cosine = principal_normal_coordinates(r, angle)
    one_plus_z = h.imag**2/(1+cosine)+1j*h.imag
    offset = h.real-h.imag**2/(1+cosine)
    integrand = one_plus_z/h*ratio_weight(offset)*r/(pi*cosine)
    return radius*pi/2*(weights @ integrand @ weights)


def principal_collar_period(radius, order=36):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    angle = pi*(nodes[:, None]+1)
    s = (nodes[None, :]+1)/2
    h, _, _, cosine = principal_normal_coordinates(radius, angle)
    one_plus_base = h.imag**2/(1+cosine)+1j*h.imag
    one_plus_z = (1-s)*one_plus_base
    z = one_plus_z-1
    integrand = one_plus_z*one_plus_base/z*ratio_weight(h-one_plus_z)/pi
    return pi/2*(weights @ integrand @ weights)


def check_common_meridian(radius, angle):
    h, z, u, _ = principal_normal_coordinates(radius, angle)
    np.testing.assert_allclose(u+z, h, atol=1e-15, rtol=0)
    np.testing.assert_allclose(winding_number(h), -1, atol=1e-12, rtol=0)
    f = h/z
    np.testing.assert_allclose(abs(f), radius, atol=1e-15, rtol=0)
    _, rotated_u, rotated_base = normal_coordinates(radius, np.angle(f))
    np.testing.assert_allclose(rotated_u, u, atol=1e-15, rtol=0)
    np.testing.assert_allclose(rotated_base, 1/z, atol=1e-15, rtol=0)
    for x in (0.3, 1., 4.):
        for s in (0., .25, .5, .75, 1.):
            base = (1-s)*z-s
            y = x*(h-base)
            assert np.all(base.real < 0) and np.all(y.real > 0)
            np.testing.assert_allclose(base*x+y, x*h, atol=1e-15, rtol=0)
        # At z=-1, y=x*(1+h). This is precisely the clockwise
        # meridian y=x+x*q in the fixed-radius variant of example 46.
        w = (1+y)/(y-x)
        np.testing.assert_allclose(winding_number(w), 1, atol=1e-12, rtol=0)


def check_principal_arc_homotopy():
    """Exercise the pole-free homotopy on retained base arcs."""
    for width in (.3, .05, .002):
        phase = np.linspace(-pi+width, pi-width, 129)[None, :]
        s = np.linspace(0, 1, 21)[:, None]
        rotation = np.exp(1j*s*phase)
        relative = np.exp(1j*(1-s)*phase)
        # For |arg(v)|<=pi-width, |a+b*v| >= sin(width/2)*(a+b)
        # when a,b>0 and |v|=1. These two distances control both poles
        # at every fixed x and positive real y radius in the homotopy.
        for a, b in ((.2, 3.), (1., 1.), (4., .25)):
            bound = np.sin(width/2)*(a+b)
            assert np.min(abs(a+b*rotation)) >= bound-1e-15
            assert np.min(abs(a+b*relative)) >= bound-1e-15


def check_cd_principal_seam():
    certify_principal_pullbacks()
    check_principal_arc_homotopy()
    angle = 2*pi*np.arange(128)/128
    previous = np.inf
    for radius in (.2, .06, .015, .004):
        check_common_meridian(radius, angle)
        disk = principal_disk_period(radius)
        collar = principal_collar_period(radius)
        rotated = disk_contribution(radius, 48)
        np.testing.assert_allclose(collar, disk, atol=2e-13, rtol=0)
        difference = disk-rotated
        np.testing.assert_allclose(rotated, disk, atol=2e-13, rtol=0)
        error = abs(disk/radius**2-.5)
        assert error < previous
        previous = error
        print(f"normal radius={radius:g}: principal disk/R^2={disk.real/radius**2:.12f}, "
              f"difference from rotated disk={difference.real:+.3e}")
    assert previous < 1e-5
    print("PASS: principal arc transport and the shared clockwise CD meridian")


if __name__ == "__main__":
    check_cd_principal_seam()
