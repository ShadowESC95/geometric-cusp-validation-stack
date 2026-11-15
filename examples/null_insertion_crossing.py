"""Resolve a crossing between an external null pinch and the insertion pole."""

from math import pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .energy_contour import finite_frame, interval, upper_residues

WIDTH = 1/32
RADIUS = 1/64
CONTOUR = (1.3, 0.14)


def lower_pinch_residues(u, r, x, eps, centers):
    """Two isolated lower energy residues, including reflection and normalization."""
    lam = (1+4*u)/5
    small = r*r-1j*eps
    dz = (lam-1.5)[..., None] - centers[:, 3]
    total = 0j
    for length, energy, slope, spectators in (
        (1-lam, 0.5-lam, -2, (1, 2)),
        (lam, 0.5+lam, -1, (2, 3)),
    ):
        omega = np.sqrt(length*length+small)
        displacement = small/(omega+length)
        pole = energy+displacement
        d = (pole[..., None]-centers[:, 0])**2-dz**2-small[..., None]
        insertion = d[..., 4]-0.25
        # The other propagator of the null pair equals slope*(E-E_axis).
        total += 1/(2*omega*slope*displacement*d[..., spectators[0]]
                    *d[..., spectators[1]]*(insertion**2-x*x))
    return 8/pi*x*x*total


def regular_lower_residues(u, r, x, eps, centers):
    """Group the remaining lower poles of both reflected integrands."""
    dz = ((1+4*u)/5-1.5)[..., None] - centers[:, 3]
    center, radius = CONTOUR
    total = np.zeros(np.broadcast_shapes(u.shape, x.shape), dtype=complex)
    for offset in radius*np.exp(2j*pi*np.arange(32)/32):
        d = (center+offset-centers[:, 0])**2-dz**2-r[..., None]**2+1j*eps
        insertion = d[..., 4]-0.25
        total += offset/(np.prod(d[..., :4], axis=-1)*(insertion**2-x*x))
    return 8/pi*x*x*total/32


def integrate_cylinder(radius, eps, centers, n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    phi = pi*(1+nodes)/4
    angles, axial_weights = np.polynomial.legendre.leggauss(32)
    # Upper semicircle oriented from -WIDTH to WIDTH.
    u = WIDTH*np.exp(1j*pi*(1-angles)/2)
    axial_weights = axial_weights*(-1j*pi/2)*u*(4/5)
    cutoff = sqrt(eps)
    knots = sorted({0, min(cutoff, radius), min(4*cutoff, radius), radius})
    total = 0j
    for a, b in zip(knots, knots[1:]):
        r = ((a+b)/2+(b-a)*nodes/2)[:, None, None]
        x = r*np.cos(phi)[None, :, None]
        values = lower_pinch_residues(u, r, x, eps, centers)
        values += regular_lower_residues(u, r, x, eps, centers)
        total += np.einsum("i,j,k,ijk->", weights*r[:, 0, 0], weights,
                           axial_weights, values)*(b-a)/2
    # Four quadrants after averaging across p_x=0.
    return pi*total


def certify_deformation(centers, eps):
    """Verify clearance for the energy contours and upper axial deformation."""
    lam_center = 1/5
    lam_shift = 4*WIDTH/5
    dz = lam_center-1.5-centers[:, 3]
    lengths = np.sqrt(dz**2+centers[:, 1]**2)
    squared_shift = 2*abs(dz)*lam_shift+lam_shift**2+RADIUS**2+eps
    squared_shift[4] += RADIUS
    assert np.all(squared_shift < lengths**2)
    shifts = squared_shift/lengths
    upper, lower = centers[:, 0]-lengths, centers[:, 0]+lengths
    center, radius = CONTOUR
    assert np.all(abs(lower[[1, 2, 4]]-center)+shifts[[1, 2, 4]] < radius)
    assert np.all(abs(lower[[0, 3]]-center)-shifts[[0, 3]] > radius)
    assert np.all(abs(upper-center)-shifts > radius)
    for lower_index, omitted_upper in ((0, 1), (3, 0)):
        for branch, energies in enumerate((lower, upper)):
            for index in range(len(centers)):
                if (branch, index) in ((0, lower_index), (1, omitted_upper), (1, 4)):
                    continue
                assert abs(energies[index]-lower[lower_index]) > shifts[index]+shifts[lower_index]
    assert abs(upper[4]-lower[0]) > shifts[4]+shifts[0]

    # At the lower D4 pole, the reflected insertion factors are H(u) +/- x,
    # where H=u-(21/8)*delta, delta=sqrt((1-lambda)^2+r^2-i*eps)-(1-lambda).
    # Im H on the real axis is nonnegative. The derivative bound gives
    # Re H' > 0 throughout the half-disk, so neither factor vanishes above it.
    min_length = 1-lam_center-lam_shift
    max_small = RADIUS**2+eps
    min_omega = sqrt(min_length**2-lam_shift**2-max_small)
    assert (21/10)*max_small/(min_omega*(min_length+min_omega)) < 1
    assert WIDTH-(21/8)*max_small/(min_length+min_omega) > RADIUS


def shrinking_cylinder_coefficient(exact_centers):
    R = sp.Rational
    lam, E, r = sp.symbols("lambda E r", real=True)
    z = lam-R(3, 2)
    d = [sp.expand(interval((E, r, 0, z), center)) for center in exact_centers]
    assert sp.expand(d[0]-d[3]) == -2*(E-(R(1, 2)-lam))
    assert sp.expand(d[1]-d[0]) == -(E-(R(1, 2)+lam))
    assert sp.expand(d[4]-d[3]+r) == -21*E/8-11*z/8-1

    low = [sp.factor(p.subs({E: R(1, 2)-lam, r: 0})) for p in d]
    high = [sp.factor(p.subs({E: R(1, 2)+lam, r: 0})) for p in d]
    u = sp.symbols("u", real=True)
    assert list(map(sp.expand, low)) == [0, 2*lam, 2*lam-2, 0, (5*lam-1)/4]
    assert list(map(sp.expand, high)) == [0, 0, 4*lam-2, 4*lam, -R(1, 4)]
    f = sp.cancel((1/(low[1]*low[2])).subs(lam, (1+4*u)/5)*R(4, 5))
    # For the D4 pole the coefficient is -2 times the boundary integral
    # of f(u)/(u+i0)^2. Subtract its constant and linear Taylor terms.
    derivative = sp.diff(f, u).subs(u, 0)
    remainder = sp.cancel((f-f.subs(u, 0)-u*derivative)/u**2)
    finite_part = (quad(sp.lambdify(u, remainder, "numpy"), -WIDTH, WIDTH)[0]
                   - 2*float(f.subs(u, 0))/WIDTH)
    delta_derivative = 2*sp.pi*derivative
    assert delta_derivative == 15*sp.pi/2
    # The second null pair, isolated at the lower D1 pole, is smooth in u.
    second = -4/(high[2]*high[3]*high[4]**2)
    second_integral = quad(sp.lambdify(lam, second, "numpy"),
                           (1-4*WIDTH)/5, (1+4*WIDTH)/5)[0]
    return second_integral - 2*finite_part + 1j*float(delta_derivative)


def check_null_insertion_crossing():
    exact_centers = finite_frame()
    centers = np.array(exact_centers, dtype=float)
    certify_deformation(centers, eps=1e-4)
    u, r, eps = np.array(0.0), np.array(RADIUS), 1e-4
    point = np.array([r, 0.0, -1.3])
    reflected = point*np.array([-1, 1, 1])
    upper = sum(upper_residues(p, centers, eps).sum() for p in (point, reflected))/pi
    lower = (lower_pinch_residues(u, r, r, eps, centers)
             + regular_lower_residues(u, r, r, eps, centers))
    np.testing.assert_allclose(lower, upper, rtol=0, atol=1e-7)

    target = shrinking_cylinder_coefficient(exact_centers)
    print(f"Predicted shrinking-cylinder coefficient: {target:.9f}")
    errors = []
    for radius in (RADIUS, 1/128, 1/256, 1/512, 1/1024):
        value = integrate_cylinder(radius, 0, centers, n=12)/radius**2
        errors.append(abs(value-target))
        print(f"radius={radius:.6f}: normalized integral/radius^2={value:.9f}")
    assert all(new < old for old, new in zip(errors, errors[1:])) and errors[-1] < 0.022

    boundary = integrate_cylinder(RADIUS, 0, centers, n=16)
    errors = []
    for eps in (1e-7, 1e-8, 1e-9):
        regulated = integrate_cylinder(RADIUS, eps, centers, n=16)
        errors.append(abs(regulated-boundary)/RADIUS**2)
        print(f"eps={eps:.0e}: scaled difference from boundary={errors[-1]:.3e}")
    refined = integrate_cylinder(RADIUS, 1e-9, centers, n=26)
    assert abs(refined-regulated)/RADIUS**2 < 2e-4
    assert errors[2] < errors[1] < errors[0] and errors[2] < 0.01
    print("PASS: insertion crossing and the finite-part cylinder limit")


if __name__ == "__main__":
    check_null_insertion_crossing()
