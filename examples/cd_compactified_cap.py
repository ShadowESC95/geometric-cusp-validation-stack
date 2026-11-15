"""Track the CD seam through projective infinity and the base poles."""

from math import pi

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_resolved_seam import normal_coordinates


def certify_compactification():
    """Check the global winding function, pole residues and intersection counts."""
    x0, x1, y0, y1, z0, z1, a, b, c = sp.symbols("x0 x1 y0 y1 z0 z1 a b c")
    section = z1*x1*y0+z0*x0*y1
    scaled_section = section.subs({x0: a*x0, x1: a*x1, y0: b*y0, y1: b*y1,
                                  z0: c*z0, z1: c*z1}, simultaneous=True)
    assert sp.expand(scaled_section-a*b*c*section) == 0
    divisor_x = x0+x1
    divisor_y = y0+y1
    for base_divisor in (z0, z1):
        denominator = base_divisor*divisor_x*divisor_y
        assert sp.cancel(scaled_section/(a*b*c*denominator)-section/denominator) == 0

    x, y, z, v = sp.symbols("x y z v")
    form = 2*(1+z)/(z*(1+x)*(1+y)*(z*x+y))
    assert sp.cancel(sp.residue(form, z, 0)-2/((1+x)*(1+y)*y)) == 0
    near_infinity = sp.cancel(-form.subs(z, 1/v)/v**2)
    assert sp.cancel(sp.residue(near_infinity, v, 0)+2/((1+x)*(1+y)*x)) == 0

    # Rows are the moving divisor, x=-1, y=-1, z=0 and z=infinity.
    # Columns are the three projective coordinate spheres.
    intersections = sp.Matrix([[1, 1, 1], [1, 0, 0], [0, 1, 0],
                               [0, 0, 1], [0, 0, 1]])
    seam_disk = sp.Matrix([0, 0, 1, 0, 0])
    degrees = intersections[:3, :].inv()*(-seam_disk[:3, :])
    assert degrees == sp.Matrix([0, -1, 1])
    cap = seam_disk+intersections*degrees
    assert cap == sp.Matrix([0, 0, 0, 1, 1])


def check_cap_for_x(x, phase):
    """Follow a meridian into an explicit disk through the two base poles."""
    radius = .06
    normal, _, base = normal_coordinates(radius, phase)
    y = normal-1
    for fraction in (0., .25, .5, .75, 1.):
        z = (1-fraction)*base-fraction
        # |z+1| < 1/2 on this straight contraction. Thus z is nonzero
        # and Re(z*x+y) < 0 for x>0 and |normal|=.06.
        assert np.all(abs(z+1) < .5)
        moving = z*x+y
        assert np.all(moving.real < 0)
        D = moving/((1+x)*normal)
        np.testing.assert_allclose([winding_number(z), winding_number(D), winding_number(D/z)],
                                   [0, -1, -1], rtol=0, atol=1e-12)

    # At z=-1, w=(y+1)/(y-x) identifies the twice-punctured y sphere
    # with C*. Express both endpoints as exp(i*phase) times a factor
    # with a continuous logarithm and interpolate those logarithms.
    q_radius = radius*x/(1+x)
    q = q_radius*np.exp(-1j*phase)
    start = np.log(radius/(1+x))+1j*pi-np.log1p(-normal/(1+x))
    end = np.log((1+x)/q_radius)+np.log1p(q/(1+x))
    np.testing.assert_allclose(np.exp(1j*phase+start), normal/(y-x), atol=1e-12)
    np.testing.assert_allclose(np.exp(1j*phase+end), (x+1+q)/q, atol=1e-12)
    for fraction in (0., .25, .5, .75, 1.):
        w = np.exp(1j*phase+(1-fraction)*start+fraction*end)
        np.testing.assert_allclose(winding_number(w), 1, rtol=0, atol=1e-12)
        # w=1 is y=infinity, and D extends there.
        np.testing.assert_allclose(winding_number(1/((1+x)*w)), -1,
                                   rtol=0, atol=1e-12)

    # Keep the clockwise moving-divisor meridian q fixed and slide it
    # from the y plane into the z plane. The endpoint has y=x.
    delta = q_radius/x
    for fraction in (0., .25, .5, .75, 1.):
        z = -1+fraction*q/x
        y = x+(1-fraction)*q
        np.testing.assert_allclose(z*x+y, q, rtol=0, atol=3e-15)
        assert np.all(y.real > 0)
        assert abs(z).min() >= 1-delta-1e-12

    # Disk |k|<=1, z=-1+delta/k. Use homogeneous coordinates, so the
    # moving section is x*delta even at k=0. The base poles are the
    # two simple zeros of z1 and z0, respectively.
    for k in (0., delta, .5j, np.exp(1j*phase)):
        z0, z1 = k, delta-k
        np.testing.assert_allclose(x*(z1+z0), q_radius, rtol=0, atol=1e-12)
    for center, local_coordinate in ((delta, lambda k: (delta-k)/k),
                                     (0., lambda k: k/(delta-k))):
        local_loop = center+delta/4*np.exp(1j*phase)
        np.testing.assert_allclose(winding_number(local_coordinate(local_loop)),
                                   1, rtol=0, atol=1e-12)


def check_cd_compactified_cap():
    certify_compactification()
    phase = 2*pi*np.arange(256)/256
    for x in (.3, 1., 4.):
        check_cap_for_x(x, phase)
        print(f"x={x:g}: projective quotient winding=-1, base-pole indices=(+1, +1)")
    print("PASS: compactified CD seam and the cap through both base poles")


if __name__ == "__main__":
    check_cd_compactified_cap()
