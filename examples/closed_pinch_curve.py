"""Integrate the intersecting-pinch contribution around its closed curve."""

from math import log1p, pi, sqrt

import numpy as np
import sympy as sp

from .energy_contour import finite_frame, finite_numerator, interval

HALF_WIDTH = 1/320


def tube_coordinates(u, v, theta):
    """Parametrize D4=0 by (D2, D_insertion, angle) in the tube."""
    z_center = (79*u/4 - 18*v - 42)/41
    radius_sq = (
        -377*u**2/576 + 7*u*v/4 + 17*u/12 - v**2 - 2*v
        + 41*z_center**2/36 - 1
    )
    radius = np.sqrt(radius_sq)
    z = z_center + 6*radius/sqrt(41)*np.cos(theta)
    y = radius*np.sin(theta)
    energy = -(u + z)/3
    x = 7*u/8 - v - z/2 - 1
    return energy, np.stack((x, y, z), axis=-1)


def regulated_d4_residue(u, v, theta, centers, eps):
    """Lower D4 energy residue including the spatial Jacobian of the tube."""
    energy, point = tube_coordinates(u, v, theta)
    distances_sq = np.sum((point[..., None, :] - centers[:, 1:])**2, axis=-1)
    omega = np.sqrt(distances_sq[..., 3] - 1j*eps)
    pole = centers[3, 0] + omega
    denominators = (pole[..., None] - centers[:, 0])**2 - distances_sq + 1j*eps
    spectators = np.prod(denominators[..., (0, 1, 2, 4)], axis=-1)
    spatial_jacobian = 2*(energy - centers[3, 0])/sqrt(41)
    return (-4*point[..., 0])*spatial_jacobian/(2*omega*spectators)


def integrate_tube(centers, eps, n, angles):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    extent = np.arcsinh(HALF_WIDTH/eps)
    u = eps*np.sinh(extent*nodes)
    w = weights*extent*eps*np.cosh(extent*nodes)
    theta = 2*pi*np.arange(angles)/angles
    values = regulated_d4_residue(u[:, None, None], u[:, None], theta, centers, eps)
    return 2*pi*(w @ values.mean(axis=-1) @ w)


def d4_tube_boundary_value():
    """Analytic integral of the lower D4 residue as eps tends to zero."""
    pv_u = pi*log1p(2*HALF_WIDTH/(2-HALF_WIDTH))
    return -pi**3 - 1j*pi*pv_u


def check_closed_pinch_curve():
    R = sp.Rational
    centers = finite_frame()
    E, x, y, z = sp.symbols("E x y z", real=True)
    u, v, w = sp.symbols("u v w", real=True)
    loop = (E, x, y, z)
    assert sp.expand(finite_numerator(loop, centers)) == -4*x
    d = [sp.expand(interval(loop, center)) for center in centers]
    affine = sp.solve((d[1]-d[3]-u, d[4]-d[3]-v), (E, x))
    assert affine == {E: -(u+z)/3, x: 7*u/8-v-z/2-1}
    d4 = sp.expand(d[3].subs(affine))
    z_center = (79*u/4-18*v-42)/41
    radius_sq = sp.expand(d4.subs({z: z_center, y: 0}))
    assert sp.expand(d4 - radius_sq + y**2 + R(41, 36)*(z-z_center)**2) == 0
    assert radius_sq.subs({u: 0, v: 0}) == R(8, 41)
    assert sp.Matrix([affine[E], affine[x]]).jacobian([u, v]).det() == R(1, 3)
    # In elliptical polar coordinates, dz dy = (3/sqrt(41))*|dw dtheta|,
    # where w=D4. The energy derivative dD4/dE is 2*(E+1/2).
    angle = sp.symbols("theta", real=True)
    radius = sp.sqrt(radius_sq-w)
    ellipse = sp.Matrix([z_center+6*radius*sp.cos(angle)/sp.sqrt(41), radius*sp.sin(angle)])
    assert sp.simplify(ellipse.jacobian([w, angle]).det()) == -3/sp.sqrt(41)

    # Certify a positive ellipse radius and the required energy branches
    # over the whole box using a coefficient bound for the quadratic.
    width = R(1, 320)
    radius_variation = sum(
        abs(coefficient)*width**sum(powers)
        for powers, coefficient in sp.Poly(radius_sq, u, v).terms()
        if sum(powers) > 0
    )
    assert 0 < R(8, 41)-radius_variation < R(8, 41)+radius_variation < R(1, 4)
    z_span = (R(79, 4)+18)*width/41 + 3/sp.sqrt(41)
    e_low = (R(42, 41)-z_span-width)/3
    e_high = (R(42, 41)+z_span+width)/3
    assert -R(1, 2) < e_low < e_high < R(13, 16)

    # On the ellipse both spectators are affine in z. Their angular
    # discriminants are positive squares, and they are negative at its center.
    d1, d3 = (sp.expand((d[i]-d[3]).subs(affine)) for i in (0, 2))
    angular_discriminants = [
        sp.factor(f.subs(z, z_center)**2 - sp.diff(f, z)**2*36*radius_sq/41)
        for f in (d1, d3)
    ]
    assert angular_discriminants == [(5*u-8*v-2)**2/41, 64*(u-v-1)**2/41]
    assert 13*width < 2 and 2*width < 1
    assert all(f.subs(z, z_center).subs({u: 0, v: 0}) < 0 for f in (d1, d3))

    # Angular integration of N/(sqrt(41)*D1*D3) by partial fractions.
    a, b = sp.symbols("a b")
    numerator = -4*affine[x]
    coefficients = sp.solve(sp.Poly(numerator-a*d1-b*d3, z).all_coeffs(), (a, b))
    angle_integral = -2*sp.pi*(
        coefficients[a]/(8*(1-u+v)) + coefficients[b]/(2-5*u+8*v)
    )
    assert sp.simplify(angle_integral - 2*sp.pi/(2-u)) == 0

    # The angular integral depends only on u. On a symmetric (u,v) box the
    # double principal value and the PV in v vanish. The PV in u survives.
    exact_contact = sp.simplify(-sp.pi**2*angle_integral.subs(u, 0))
    assert exact_contact == -sp.pi**3
    contact = float(exact_contact)
    pv_primitive = sp.pi*sp.log((2+u)/(2-u))
    assert sp.simplify(sp.diff(pv_primitive, u)
                       - (angle_integral-angle_integral.subs(u, -u))/u) == 0
    pv_u = float(pv_primitive.subs(u, width))
    target = d4_tube_boundary_value()
    np.testing.assert_allclose(target, contact - 1j*pi*pv_u, rtol=0, atol=1e-12)
    numeric_centers = np.array(centers, dtype=float)
    theta = 2*pi*np.arange(128)/128
    for uu, vv in ((0, 0), (HALF_WIDTH/2, -HALF_WIDTH/2), (-HALF_WIDTH, HALF_WIDTH)):
        energy, spatial = tube_coordinates(uu, vv, theta)
        distances_sq = np.sum((spatial[:, None, :]-numeric_centers[:, 1:])**2, axis=-1)
        denominators = (energy[:, None]-numeric_centers[:, 0])**2 - distances_sq
        shell = np.broadcast_to([uu, 0, vv], (len(theta), 3))
        np.testing.assert_allclose(denominators[:, [1, 3, 4]], shell, atol=1e-12)
        density = -4*spatial[:, 0]/(sqrt(41)*denominators[:, 0]*denominators[:, 2])
        np.testing.assert_allclose(2*pi*density.mean(), 2*pi/(2-uu), rtol=0, atol=1e-11)

    print(f"Closed-curve double-delta term={contact:.9f}, tube target={target:.9f}")
    errors = []
    for eps in (1e-5, 1e-6, 1e-7):
        value = integrate_tube(numeric_centers, eps, n=72, angles=96)
        errors.append(abs(value-target))
        print(f"eps={eps:.0e}: D4 residue over the full tube={value:.9f}")
    np.testing.assert_allclose(
        integrate_tube(numeric_centers, 1e-7, n=88, angles=128), value,
        rtol=0, atol=1e-6,
    )
    assert errors[2] < errors[1] < errors[0] and errors[2] < 0.003
    print("PASS: regulated D4 residue and the closed intersection curve")


if __name__ == "__main__":
    check_closed_pinch_curve()
