"""Five-propagator energy integral and a transverse Feynman pinch."""

from itertools import combinations
from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .geometry import EXTERNAL, affine, lorentzian_line


def finite_frame(cross_ratio=1):
    """Return the rho=1/2 frame with the requested cross-ratio."""
    z = sp.sympify(cross_ratio)
    ab = lorentzian_line(-1, 1/(2*(1+z)), -sp.Rational(1, 2))
    matrices = [affine(line)[1] for line in (*EXTERNAL, ab)]
    return [
        ((X[0, 0]+X[1, 1])/2, sp.re(X[0, 1]), -sp.im(X[0, 1]),
         (X[0, 0]-X[1, 1])/2)
        for X in matrices
    ]


def interval(first, second):
    """Mostly-minus squared interval, usable with exact or numeric inputs."""
    return (first[0]-second[0])**2 - sum(
        (a-b)**2 for a, b in zip(first[1:], second[1:])
    )


def finite_numerator(loop, centers):
    """Combined insertion numerator of arXiv:2112.06956, Eq. (2.30)."""
    d1, d2, d3, d4, d0 = (interval(loop, center) for center in centers)
    a1, a2, a3, a4 = (interval(centers[-1], center) for center in centers[:-1])
    s = interval(centers[0], centers[2])
    t = interval(centers[1], centers[3])
    return t*(a3*d1+a1*d3) + s*(a4*d2+a2*d4) - s*t*d0


def energy_poles(spatial, centers, eps):
    squared_momenta = np.sum((np.asarray(spatial)[..., None, :] - centers[:, 1:])**2, axis=-1)
    omega = np.sqrt(squared_momenta - 1j*eps)
    return centers[:, 0]-omega, centers[:, 0]+omega


def upper_residues(spatial, centers, eps, x_coefficient=-4):
    spatial = np.asarray(spatial)
    centers = np.asarray(centers)
    upper, lower = energy_poles(spatial, centers, eps)
    differences = upper[..., :, None]-upper[..., None, :]
    indices = np.arange(len(centers))
    differences[..., indices, indices] = 1
    factors = differences*(upper[..., :, None]-lower[..., None, :])
    return x_coefficient*spatial[..., :1]/np.prod(factors, axis=-1)


def direct_energy_integral(spatial, centers, eps):
    upper, lower = energy_poles(spatial, centers, eps)
    bounds = np.r_[-np.inf, np.sort(np.r_[upper.real, lower.real]), np.inf]

    def integrand(energy):
        loop = (energy, *spatial)
        denominators = [interval(loop, point) + 1j*eps for point in centers]
        return finite_numerator(loop, centers) / np.prod(denominators)

    return sum(
        quad(integrand, lo, hi, complex_func=True,
             epsabs=1e-10, epsrel=1e-10, limit=200)[0]
        for lo, hi in zip(bounds, bounds[1:])
    )


def check_energy_contour():
    exact_centers = finite_frame()
    energy, x, y, z = sp.symbols("E x y z", real=True)
    assert sp.expand(finite_numerator((energy, x, y, z), exact_centers)) == -4*x
    centers = np.array(exact_centers, dtype=float)

    for point in ([0.4, 0.35, -0.65], [-0.3, -0.25, 0.4], [0.2, 0.2, -1.2]):
        for eps in (0.08, 0.02):
            residues = upper_residues(point, centers, eps)
            direct = direct_energy_integral(point, centers, eps)
            result = 2j*pi*residues.sum()
            np.testing.assert_allclose(direct, result, rtol=0, atol=1e-9)
            print(f"spatial={point}, eps={eps:g}: energy mismatch={abs(direct-result):.2e}")

    # An upper pole of D2 crosses a lower pole of D4 when p_z=-3/20.
    # Keep the transverse coordinates fixed and integrate across the crossing.
    R = sp.Rational
    crossing = (R(1, 20), R(3, 10), R(3, 10), -R(3, 20))
    width = R(1, 25)
    slice_point = (energy, *crossing[1:3], z)
    denominators = [sp.expand(interval(slice_point, point)) for point in exact_centers]
    pinches = []
    for i, j in combinations(range(len(denominators)), 2):
        meeting_energy, = sp.solve(denominators[i] - denominators[j], energy)
        intersection = sp.Poly(denominators[i].subs(energy, meeting_energy), z)
        for pz in intersection.real_roots():
            if abs(pz - crossing[3]) > width:
                continue
            at = meeting_energy.subs(z, pz)
            if (at - exact_centers[i][0]) * (at - exact_centers[j][0]) < 0:
                pinches.append((i, j, at, pz))
    assert pinches == [(1, 3, crossing[0], crossing[3])]

    d = [interval(crossing, point) for point in exact_centers]
    assert d == [-R(9, 5), 0, -R(1, 5), 0, -R(49, 40)]
    omega2 = exact_centers[1][0]-crossing[0]
    omega4 = crossing[0]-exact_centers[3][0]
    slope = (-(crossing[3]-exact_centers[1][3])/omega2
             -(crossing[3]-exact_centers[3][3])/omega4)
    assert slope == -R(320, 209)
    # The energy residue gives 2*pi*i. The spatial boundary value
    # of 1/(gap+i0) then contributes -i*pi*delta(gap).
    spike_weight = sp.factor(2*sp.pi**2*finite_numerator(crossing, exact_centers)
                             /(-4*omega2*omega4*abs(slope)*d[0]*d[2]*d[4]))
    assert spike_weight == -250*sp.pi**2/147

    center = float(crossing[3])
    radius = float(width)
    errors = []
    for eps in (1e-3, 1e-4, 1e-5):
        def real_energy_integral(pz):
            residues = upper_residues([0.3, 0.3, pz], centers, eps)
            return (2j*pi*residues.sum()).real

        nodes = np.array([-radius, -16*eps, -4*eps, -eps, 0,
                          eps, 4*eps, 16*eps, radius])
        bounds = center + np.unique(np.clip(nodes, -radius, radius))
        value = sum(
            quad(real_energy_integral, lo, hi, epsabs=1e-10, epsrel=1e-10)[0]
            for lo, hi in zip(bounds, bounds[1:])
        )
        errors.append(abs(value-float(spike_weight)))
        print(f"eps={eps:.0e}: integrated peak={value:.9f}, "
              f"pinch prediction={float(spike_weight):.9f}")
    assert errors[2] < errors[1] < errors[0] and errors[2] < 0.003
    print("PASS: energy residues and the integrated spatial pinch")


if __name__ == "__main__":
    check_energy_contour()
