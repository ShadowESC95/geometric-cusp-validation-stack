"""Resolve a merging pair of Feynman pinches in a local spatial patch."""

from itertools import combinations
from math import pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .energy_contour import finite_frame, finite_numerator, interval, upper_residues


Q_LIMIT = 0.055
LAMBDA_LIMIT = 0.004
RATIO_LIMITS = (0.8, 1.25)


def spatial_point(lam, q, ratio):
    """Use lambda=p_perp^2-1/2, q=p_z+3/4, and ratio=p_y/p_x."""
    x = sqrt((0.5 + lam) / (1 + ratio**2))
    return x, ratio*x, q - 0.75


def shell_lambda(q):
    return -8*q**2/9


def energy_real(lam, q, ratio, centers, eps):
    residues = upper_residues(spatial_point(lam, q, ratio), centers, eps)
    return (2j*pi*residues.sum()).real


def limiting_surface_density(q, ratio, centers):
    """Integrate out energy and lambda on the D2=D4=0 surface."""
    point = (0.25 - q/3, *spatial_point(shell_lambda(q), q, ratio))
    spectators = np.prod([interval(point, centers[i]) for i in (0, 2, 4)])
    transverse_jacobian = 1/(2*(1 + ratio**2))
    numerator = -4*point[1]
    # The energy residue and spatial boundary value give 2*pi^2.
    # D2 has negative energy derivative and |d(D2,D4)/d(E,lambda)| = 3.
    return 2*pi**2*numerator/(-3*spectators) * transverse_jacobian


def integrate_patch(centers, n, eps=None):
    """Use Gauss-Legendre in q and p_y/p_x and adaptive quadrature in lambda."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    lo, hi = RATIO_LIMITS
    ratios = (lo + hi)/2 + (hi - lo)/2*nodes
    total = 0.0
    for q, q_weight in zip(Q_LIMIT*nodes, weights):
        for ratio, ratio_weight in zip(ratios, weights):
            if eps is None:
                density = limiting_surface_density(q, ratio, centers)
            else:
                jacobian = 1/(2*(1 + ratio**2))
                shell = shell_lambda(q)
                density = jacobian * quad(
                    energy_real, -LAMBDA_LIMIT, LAMBDA_LIMIT,
                    args=(q, ratio, centers, eps), points=[shell],
                    epsabs=1e-9, epsrel=1e-9, limit=160,
                )[0]
            total += q_weight*ratio_weight*density
    return total * Q_LIMIT*(hi - lo)/2


def check_spatial_pinch():
    exact_centers = finite_frame()
    # Only timelike pairs can have opposite-energy branches meet away from a
    # null segment. The other timelike pair (D4, insertion) cannot meet here:
    # |p-P4| + |p-P0| exceeds the energy separation 21/16 throughout the patch.
    assert interval(exact_centers[1], exact_centers[3]) == 2
    assert interval(exact_centers[3], exact_centers[4]) == 1
    null_pairs = set()
    timelike_pairs = set()
    for i, j in combinations(range(len(exact_centers)), 2):
        separation = interval(exact_centers[i], exact_centers[j])
        if separation == 0:
            null_pairs.add((i, j))
        elif separation > 0:
            timelike_pairs.add((i, j))
    assert timelike_pairs == {(1, 3), (3, 4)}
    assert all(exact_centers[k][1:3] == (0, 0) for pair in null_pairs for k in pair)
    assert Q_LIMIT < 0.25 and 8*Q_LIMIT**2/9 < LAMBDA_LIMIT
    x_min = sqrt((0.5 - LAMBDA_LIMIT)/(1 + RATIO_LIMITS[1]**2))
    omega4_min = sqrt(0.5 - LAMBDA_LIMIT + (0.25 - Q_LIMIT)**2)
    assert omega4_min + x_min + 0.5 > 21/16

    R = sp.Rational
    E, q, lam = sp.symbols("E q lambda", real=True)
    x, y = sp.symbols("x y", positive=True)
    to_patch = sp.Matrix([x**2 + y**2 - R(1, 2), y/x])
    assert sp.simplify(to_patch.jacobian([x, y]).det()) == 2 + 2*y**2/x**2
    d2 = (E - 1)**2 - (q + R(1, 4))**2 - (R(1, 2) + lam)
    d4 = (E + R(1, 2))**2 - (q - R(1, 4))**2 - (R(1, 2) + lam)
    meeting, = sp.solve(d2 - d4, E)
    assert meeting == R(1, 4) - q/3
    assert sp.expand(d2.subs(E, meeting)) == -lam - 8*q**2/9
    assert sp.Matrix([d2, d4]).jacobian([E, lam]).det() == 3
    assert (meeting - 1).subs(q, -Q_LIMIT) < 0
    assert (meeting + R(1, 2)).subs(q, Q_LIMIT) > 0

    omega2 = sp.sqrt((q + R(1, 4))**2 + R(1, 2) + lam)
    omega4 = sp.sqrt((q - R(1, 4))**2 + R(1, 2) + lam)
    gap = R(3, 2) - omega2 - omega4
    fold = {q: 0, lam: 0}
    assert gap.subs(fold) == sp.diff(gap, q).subs(fold) == 0
    kappa = -sp.diff(gap, q, 2).subs(fold)/2
    regulator_scale = (1/(2*omega2) + 1/(2*omega4)).subs(fold)
    loop = (R(1, 4), R(1, 2), R(1, 2), -R(3, 4))
    spectators = sp.prod(interval(loop, exact_centers[i]) for i in (0, 2, 4))
    amplitude = finite_numerator(loop, exact_centers) / (-R(9, 4)*spectators)
    # At lambda=0 the gap is -kappa*q^2 + i*regulator_scale*eps.
    # With q proportional to sqrt(eps), integral(dx/(1+x^4)) = pi/sqrt(2).
    fold_limit = sp.simplify(2*sp.pi**2*amplitude/sp.sqrt(2*kappa*regulator_scale))
    assert fold_limit == -8*sp.pi**2/9

    centers = np.array(exact_centers, dtype=float)
    errors = []
    for eps in (1e-4, 1e-5, 1e-6):
        integral = quad(
            lambda q: energy_real(0, q, 1, centers, eps), -Q_LIMIT, Q_LIMIT,
            points=[0], epsabs=1e-8, epsrel=1e-9,
        )[0]
        scaled = sqrt(eps)*integral
        errors.append(abs(scaled - float(fold_limit)))
        print(f"eps={eps:.0e}: sqrt(eps)*fold integral={scaled:.9f}, "
              f"target={float(fold_limit):.9f}")
    assert errors[2] < errors[1] < errors[0] and errors[2] < 1e-5

    target = integrate_patch(centers, n=9)
    np.testing.assert_allclose(integrate_patch(centers, n=7), target, rtol=0, atol=1e-9)
    errors = []
    for eps in (1e-4, 1e-5, 1e-6):
        value = integrate_patch(centers, n=5, eps=eps)
        errors.append(abs(value-target))
        print(f"eps={eps:.0e}: spatial patch={value:.9f}, surface target={target:.9f}")
    np.testing.assert_allclose(
        integrate_patch(centers, n=7, eps=1e-6), value, rtol=0, atol=1e-7,
    )
    assert errors[2] < errors[1] < errors[0] and errors[2] < 4e-4
    print("PASS: merging pinch roots and the integrated spatial patch")


if __name__ == "__main__":
    check_spatial_pinch()
