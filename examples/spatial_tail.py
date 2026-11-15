"""Bound and integrate the large-spatial-momentum tail in the finite frame."""

from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .energy_contour import direct_energy_integral, finite_frame, finite_numerator

MIN_RADIUS = 6


def reflected_centers(centers):
    """Append the insertion reflected across p_x=0."""
    e, x, y, z = centers[-1]
    return [*centers, (e, -x, y, z)]


def wick_energy(spatial, centers, eps):
    radius = np.linalg.norm(spatial)
    scaled_centers = np.asarray(centers, dtype=float)/radius
    time_centers = scaled_centers[:, 0]
    distances_sq = np.sum((spatial/radius-scaled_centers[:, 1:])**2, axis=-1)

    def scaled_integrand(t):
        factors = distances_sq - (1j*t-time_centers)**2 - 1j*eps/radius**2
        return 1/np.prod(factors)

    value = quad(scaled_integrand, -np.inf, np.inf, complex_func=True,
                 epsabs=1e-11, epsrel=1e-11)[0]
    return 4j*spatial[0]/radius**9*value


def tail_quadrature(radius, centers, n, eps=0, x_coefficient=-4):
    """Integrate the normalized tail using the reflection-symmetrized kernel."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    # r = radius/s. Write Euclidean energy as r*tan(alpha).
    s = ((1+nodes)/2).reshape(n, 1, 1, 1)
    mu = nodes.reshape(1, n, 1, 1)
    alpha = (pi/2*nodes).reshape(1, 1, n, 1)
    phi = (2*pi*np.arange(n)/n).reshape(1, 1, 1, n)
    nx = np.sqrt(1-mu**2)*np.cos(phi)
    inverse_radius = s/radius
    times_cos = 2j*np.sin(alpha)*np.cos(alpha)
    cos_sq = np.cos(alpha)**2
    factors = np.ones((n,)*4, dtype=complex)
    for e, x, y, z in np.asarray(reflected_centers(centers), dtype=float):
        assert y == 0
        linear = times_cos*e - 2*(nx*x+mu*z)*cos_sq
        constant = x*x+y*y+z*z-e*e-1j*eps
        factors *= 1 + inverse_radius*linear + inverse_radius**2*cos_sq*constant

    density = s**5*nx**2*cos_sq**5/factors
    integral = np.einsum("i,j,k,ijkl->", weights, weights, weights, density)/n
    # The radial, energy-angle and azimuth weights supply pi^2/(2*n)
    # before including the reflected numerator coefficient
    # -2*x_coefficient*p_insertion_x/pi^2.
    return -x_coefficient*float(centers[-1][1])*integral/radius**6


def asymptotic_coefficients(centers):
    """Exact coefficients at orders R^-6 and R^-8 for eps=0."""
    x, z, t = sp.symbols("x z t", real=True)
    all_centers = reflected_centers(centers)
    linear = [2*(sp.I*t*e-x*a-z*c) for e, a, _, c in all_centers]
    constant = [a*a+b*b+c*c-e*e for e, a, b, c in all_centers]
    quadratic = sp.expand((sum(linear)**2 + sum(l*l for l in linear))/2)
    # Spherical moments weighted by n_x^2. Terms odd in t or n_z vanish.
    weighted = (sp.Rational(3, 5)*quadratic.coeff(x, 2)
                + quadratic.coeff(z, 2)/5 + quadratic.coeff(t, 2)*t**2)

    def moment(power):
        # Integral of (1+t^2)^(-power) over the real line.
        return sp.pi*sp.binomial(2*power-2, power-1)/4**(power-1)

    second = (-sum(constant)*moment(7) + weighted.subs(t, 0)*moment(8)
              + weighted.coeff(t, 2)*(moment(7)-moment(8)))
    measure = -sp.Rational(16, 3)/sp.pi
    return sp.simplify(measure*moment(6)/6), sp.simplify(measure*second/8)


def check_spatial_tail():
    centers = finite_frame()
    E, x, y, z = sp.symbols("E x y z", real=True)
    assert sp.expand(finite_numerator((E, x, y, z), centers)) == -4*x
    for e, a, b, c in reflected_centers(centers):
        assert a*a+b*b+c*c <= sp.Rational(9, 4)
        assert abs(e) <= 1 and a*a+b*b+c*c-e*e >= 0
    assert MIN_RADIUS-sp.Rational(3, 2) > 1
    assert MIN_RADIUS/2 >= 3
    positive_x = sp.symbols("p_x", real=True)
    common, insertion = sp.symbols("A_common A_insertion", nonzero=True)
    # A_insertion(-p_x) = A_insertion(p_x) - 2*p_x.
    reflected = (
        4*positive_x/(common*insertion)
        - 4*positive_x/(common*(insertion-2*positive_x))
    )/2
    assert sp.cancel(
        reflected + 4*positive_x**2/(common*insertion*(insertion-2*positive_x))
    ) == 0

    # For r >= 6 each rotated factor has Re A >= (r^2+t^2)/2.
    # After reflection the normalized kernel is -4*p_x^2/(pi^2*prod(A)).
    # Integrating this majorant yields 14/R^6, uniformly in eps >= 0.
    leading, second = asymptotic_coefficients(centers)
    assert leading == -sp.Rational(7, 32)
    assert second == -sp.Rational(50809, 40960)
    bound = float(-2**6*leading)
    assert bound == 14

    for point in (np.array([6., 1., -2.]), np.array([-5., 4., 1.])):
        eps = 0.2
        direct = direct_energy_integral(point, np.asarray(centers, dtype=float), eps)
        rotated = wick_energy(point, centers, eps)
        np.testing.assert_allclose(direct, rotated, rtol=0, atol=1e-12)

    errors = []
    for radius in (MIN_RADIUS, 12, 24, 48):
        value = tail_quadrature(radius, centers, n=24)
        refined = tail_quadrature(radius, centers, n=32)
        np.testing.assert_allclose(
            value*radius**6, refined*radius**6, rtol=0, atol=3e-9,
        )
        assert abs(refined) <= bound/radius**6
        scaled = refined.real*radius**6
        correction = radius**2*(scaled - float(leading))
        errors.append(abs(correction-float(second)))
        print(f"R={radius:2d}: R^6 tail={scaled:.10f}, "
              f"scaled subleading term={correction:.9f}")
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 0.003
    for eps in (0.01, 0.001):
        value = tail_quadrature(MIN_RADIUS, centers, n=24, eps=eps)
        assert abs(value) <= bound/MIN_RADIUS**6
    print("PASS: Wick rotation, uniform tail bound and both asymptotic coefficients")


if __name__ == "__main__":
    check_spatial_tail()
