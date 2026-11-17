"""The parametric representation of the shifted four-point integral.

See docs/GATE_6_CONTOUR.md, section 10. For external points x_i = (t_i, xi_i),
powers a_i with sum 4, mu > 0 and a rotation angle 0 < phi <= pi/2, put

    Y_i(phi) = mu + |xi - xi_i|^2 - exp(-2 i phi) (tau - t_i)^2 .

All Y_i lie in the closed cone spanned by 1 and exp(i(pi - 2 phi)), which has
opening less than pi, and |Y_i| >= sin(phi) (mu + |x0 - x_i|_E^2). The
integral over real (tau, xi) is therefore absolutely convergent, Feynman's
formula applies pointwise, and the Gaussian-type integral is done by rotating
tau inside the cone. The result is

    int d tau d^3 xi  prod_i Y_i^(-a_i)
      = -i exp(i phi) pi^2 / prod Gamma(a_i)
        * int_simplex prod alpha_i^(a_i - 1) / c_alpha(phi)^2,
    c_alpha(phi) = mu - sum_{i<k} alpha_i alpha_k x_ik^2(phi),
    x_ik^2(phi) = exp(-2 i phi) (t_i - t_k)^2 - |xi_i - xi_k|^2 .

As phi -> 0+ the left side is the Feynman (-i0) integral and the right side
is the parametric formula with the Minkowski invariants.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.stats import qmc


def rotated_square(xi, xk, phi):
    """x_ik^2 with the time difference squared multiplied by exp(-2 i phi)."""
    d = np.asarray(xi) - np.asarray(xk)
    return np.exp(-2j * phi) * d[0] ** 2 - d[1:] @ d[1:]


def lhs_position_space(points, powers, mu, phi, log2_samples: int = 20, scale: float = 2.0,
                       seed: int = 0):
    """Quasi-Monte Carlo value of int d tau d^3 xi prod Y_i^(-a_i), with an error estimate.

    Each coordinate is mapped by x = centre + scale tan(pi (u - 1/2)).
    """
    points = np.asarray(points, dtype=float)
    centre = points.mean(axis=0)
    estimates = []
    for rep in range(8):
        sampler = qmc.Sobol(d=4, scramble=True, seed=seed + rep)
        u = sampler.random_base2(log2_samples - 3)
        ang = np.pi * (u - 0.5)
        x = centre + scale * np.tan(ang)
        jac = np.prod(scale * np.pi / np.cos(ang) ** 2, axis=1)
        value = np.ones(len(x), dtype=complex)
        for xi, a in zip(points, powers):
            Y = mu + np.sum((x[:, 1:] - xi[1:]) ** 2, axis=1) - np.exp(-2j * phi) * (x[:, 0] - xi[0]) ** 2
            value *= Y ** (-a)
        estimates.append(np.mean(value * jac))
    estimates = np.array(estimates)
    return estimates.mean(), np.std(estimates) / np.sqrt(len(estimates))


def rhs_parametric(points, powers, mu, phi, log2_samples: int = 18, seed: int = 0):
    """-i exp(i phi) pi^2 / Gamma(4) times the Dirichlet(a) mean of 1 / c_alpha^2."""
    from scipy.stats import gamma as gamma_dist

    points = np.asarray(points, dtype=float)
    powers = np.asarray(powers, dtype=float)
    estimates = []
    for rep in range(8):
        sampler = qmc.Sobol(d=4, scramble=True, seed=seed + rep)
        u = sampler.random_base2(log2_samples - 3)
        g = np.column_stack([gamma_dist.ppf(u[:, k], powers[k]) for k in range(4)])
        alpha = g / g.sum(axis=1, keepdims=True)
        c = np.full(len(alpha), mu, dtype=complex)
        for i in range(4):
            for k in range(i + 1, 4):
                c -= alpha[:, i] * alpha[:, k] * rotated_square(points[i], points[k], phi)
        estimates.append(np.mean(1 / c**2))
    estimates = np.array(estimates)
    factor = -1j * np.exp(1j * phi) * np.pi**2 / gamma_fn(powers.sum())
    return factor * estimates.mean(), abs(factor) * np.std(estimates) / np.sqrt(len(estimates))


def minkowski_limit(points, powers, mu, log2_samples: int = 18, seed: int = 0):
    """The phi -> 0+ limit of the right side: the Minkowski parametric formula."""
    return rhs_parametric(points, powers, mu, 0.0, log2_samples, seed)


def family_check(contour, n: int, mu: float = 0.05, seed: int = 0, phi_floor: float = 1e-3):
    """The two-parameter family joining the rotated-metric integral to the contour.

    For phi >= 0 and 0 <= lam <= 1 put
        Y_i = mu - [exp(-2 i phi) (z0 - t_i)^2 - |zvec - xi_i|^2],  z = x0 + i lam V(x0).
    phi > 0, lam = 0 is the rotated-metric integral on the real domain;
    phi = 0, lam = 1 is the contour of this document.

    For phi > 0 no Y_i vanishes: Im Y_i >= sin(2 phi) (z0 - t_i)^2 up to terms
    of order V. At phi = 0 and mu > 0 a single Y_i vanishes on curves inside
    the flat triangles, where V = 0; the singularity is integrable.

    Returns, over the sample with phi >= phi_floor: the smallest
    |Y_i| / (sin(2 phi) (t - t_i)^2 + mu sin(phi)), which is bounded below if
    the family has no zeros, and the largest excursion of arg Y_i outside
    [0, pi] on the branch continued from the Feynman prescription.
    """
    from .feynman_contour import sample_points

    rng = np.random.default_rng(seed)
    smallest, excursion = np.inf, 0.0
    for x0 in sample_points(contour, n, seed):
        if contour.distance_to_edges(x0) < 1e-6:
            continue
        v = contour.field(x0)
        for _ in range(6):
            phi = 10 ** rng.uniform(np.log10(phi_floor), np.log10(0.6))
            lam = rng.choice([0.0, 1.0]) if rng.random() < 0.3 else rng.uniform(0, 1)
            z = x0 + 1j * lam * v
            for xi in contour.x:
                d = z - xi
                Y = mu - (np.exp(-2j * phi) * d[0] ** 2 - d[1:] @ d[1:])
                scale = np.sin(2 * phi) * (x0[0] - xi[0]) ** 2 + mu * np.sin(phi)
                smallest = min(smallest, abs(Y) / scale)
                a = np.angle(Y)
                if a < -np.pi / 2:
                    a += 2 * np.pi
                excursion = max(excursion, -a, a - np.pi, 0.0)
    return {"min_abs_over_scale": float(smallest), "max_excursion": float(excursion)}
