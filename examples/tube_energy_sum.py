"""Integrate the complete energy contour over the closed-pinch tube."""

from math import pi, sqrt

import numpy as np
import sympy as sp

from .closed_pinch_curve import (
    HALF_WIDTH,
    d4_tube_boundary_value,
    integrate_tube,
    regulated_d4_residue,
    tube_coordinates,
)
from .energy_contour import direct_energy_integral, finite_frame, interval

CONTOUR_LEFT = 49/80
CONTOUR_RIGHT = 21/10
CHECK_REGULATOR = 1e-4


def grouped_residue(energy, point, centers, eps, nodes):
    """Counterclockwise contour enclosing the four lower poles other than D4."""
    distances_sq = np.sum((point[..., None, :] - centers[:, 1:])**2, axis=-1)
    jacobian = 2*(energy - centers[3, 0])/sqrt(41)
    midpoint = (CONTOUR_LEFT + CONTOUR_RIGHT)/2
    radius = (CONTOUR_RIGHT - CONTOUR_LEFT)/2
    total = np.zeros(np.shape(energy), dtype=complex)
    for offset in radius*np.exp(2j*pi*np.arange(nodes)/nodes):
        denominators = (midpoint+offset-centers[:, 0])**2 - distances_sq + 1j*eps
        total += offset/np.prod(denominators, axis=-1)
    return (-4*point[..., 0])*jacobian*total/nodes


def integrate_group(centers, eps, order, angles, nodes):
    u, weights = np.polynomial.legendre.leggauss(order)
    u, weights = HALF_WIDTH*u, HALF_WIDTH*weights
    theta = 2*pi*np.arange(angles)/angles
    energy, point = tube_coordinates(u[:, None, None], u[:, None], theta)
    values = grouped_residue(energy, point, centers, eps, nodes)
    return 2*pi*(weights @ values.mean(axis=-1) @ weights)


def certify_contour(centers):
    """Bound all real pole positions on the unregulated tube."""
    R = sp.Rational
    F, u, v, a = sp.symbols("F u v a", real=True)
    width = R(1, 320)
    z_origin = -R(42, 41)
    z = z_origin + a
    z_center_offset = (79*u/4-18*v)/41
    E = -(u+z)/3
    x = 7*u/8-v-z/2-1
    d4_at_y_zero = interval((E, x, 0, z), centers[3])
    radius_sq = d4_at_y_zero.subs(a, z_center_offset)
    constant = R(8, 41)
    radius_variation = sum(
        abs(coefficient)*width**sum(powers)
        for powers, coefficient in sp.Poly(radius_sq-constant, u, v).terms()
    )
    a_bound = R(17, 40)
    assert (R(79, 4)+18)*width/41 + 6*sp.sqrt((constant+radius_variation)/41) < a_bound
    y_sq = d4_at_y_zero
    e_low = (R(42, 41)-a_bound-width)/3
    e_high = (R(42, 41)+a_bound+width)/3
    assert -R(1, 2) < e_low < e_high < R(1, 2)

    left_outside, left_inside = R(3, 5), R(5, 8)
    right_test = R(83, 40)
    for index in (0, 1, 2, 4):
        d = sp.expand(interval((F, x, sp.sqrt(y_sq), z), centers[index]))
        bounds = []
        for test in (left_outside, left_inside, right_test):
            poly = sp.Poly(d.subs(F, test), u, v, a)
            value = poly.coeff_monomial(1)
            variation = sum(
                abs(coefficient)*width**(j+k)*a_bound**l
                for (j, k, l), coefficient in poly.terms() if j+k+l > 0
            )
            bounds.append((value-variation, value+variation))
        assert bounds[0][1] < 0  # upper-half-plane pole lies left of 3/5
        assert bounds[1][1] < 0  # grouped lower pole lies right of 5/8
        assert bounds[2][0] > 0 and right_test > centers[index][0]

    # Even at the largest sampled eps, the square-root pole shift has
    # magnitude at most sqrt(eps). The four included poles have clearance
    # on both sides, and the six excluded poles stay to the left.
    max_shift = sqrt(CHECK_REGULATOR)
    assert float(left_outside)+max_shift < CONTOUR_LEFT < float(left_inside)-max_shift
    assert float(right_test)+max_shift < CONTOUR_RIGHT
    assert float(e_high)+max_shift < CONTOUR_LEFT


def check_tube_energy_sum():
    exact_centers = finite_frame()
    certify_contour(exact_centers)
    centers = np.array(exact_centers, dtype=float)

    # D1 and D3 have coincident lower poles here at eps=0.
    z_crossing = (-15-3*sqrt(5))/20
    theta_crossing = np.arccos((z_crossing+42/41)*41/(12*sqrt(2)))
    energy, point = tube_coordinates(0, 0, theta_crossing)
    poles = centers[:, 0] + np.linalg.norm(point-centers[:, 1:], axis=-1)
    np.testing.assert_allclose(poles[0], poles[2], rtol=0, atol=1e-12)
    assert np.isfinite(grouped_residue(energy, point, centers, 0, nodes=240))

    for u, v, theta in ((0, 0, theta_crossing), (HALF_WIDTH/2, -HALF_WIDTH/2, 1.3)):
        energy, point = tube_coordinates(u, v, theta)
        # Use finite eps for the independent real-axis energy quadrature.
        eps = CHECK_REGULATOR
        group = grouped_residue(energy, point, centers, eps, nodes=240)
        singular = regulated_d4_residue(u, v, theta, centers, eps)
        jacobian = 2*(energy-centers[3, 0])/sqrt(41)
        direct = direct_energy_integral(point, centers, eps)*jacobian
        np.testing.assert_allclose(-2j*pi*(group+singular), direct, rtol=1e-10, atol=1e-8)

    regular_limit = integrate_group(centers, 0, order=7, angles=96, nodes=192)
    np.testing.assert_allclose(
        regular_limit, integrate_group(centers, 0, order=5, angles=64, nodes=160),
        rtol=0, atol=1e-9,
    )
    target = -2j*pi*(d4_tube_boundary_value() + regular_limit)
    print(f"Other four lower poles at eps=0: {regular_limit.real:.12f}")
    print(f"Full tube energy integral limit: {target:.9f}")
    errors = []
    for eps in (1e-5, 1e-6, 1e-7):
        group = integrate_group(centers, eps, order=5, angles=64, nodes=160)
        singular = integrate_tube(centers, eps, n=72, angles=96)
        value = -2j*pi*(group + singular)
        errors.append(abs(value-target))
        print(f"eps={eps:.0e}: energy integral={value:.9f}, error={errors[-1]:.3e}")
    assert errors[2] < errors[1] < errors[0] and errors[2] < 0.014
    print("PASS: full energy contour over the bounded tube")


if __name__ == "__main__":
    check_tube_energy_sum()
