"""Check the full finite-frame loop integral at non-unit cross-ratios."""

from itertools import combinations
from math import log, pi

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.stats import qmc

from .cd_integral import reduced_cd
from .energy_contour import finite_frame, finite_numerator, interval
from .full_finite_frame import ROTATION, SPATIAL_CUTOFF, estimate, spatial_integral
from .spatial_tail import tail_quadrature


def frame_parameters(z):
    """Exact affine numerator and projective coefficients for the rho=1/2 family."""
    centers = finite_frame(z)
    p1, p2, p3, p4 = [sp.factor(-interval(centers[4], p)) for p in centers[:4]]
    assert all(p.is_positive for p in (p1, p2, p3, -p4))
    S = -interval(centers[0], centers[2])
    T = interval(centers[1], centers[3])
    assert S == T == 2
    assert sp.simplify(-p1*p3/(p2*p4)-z) == 0
    loop = sp.symbols("E x y p_z", real=True)
    x = loop[1]
    numerator = sp.Poly(sp.cancel(finite_numerator(loop, centers)), x).as_expr()
    coefficient = numerator.coeff(x)
    assert numerator == coefficient*x

    a1, a2, a3, a4 = a = sp.symbols("a1 a2 a3 a4", positive=True)
    projective = (*a, 1)
    delta = -sum(projective[i]*projective[j]*interval(centers[i], centers[j])
                 for i, j in combinations(range(5), 2))
    expected = S*a1*a3-T*a2*a4+sum(ai*pj for ai, pj in zip(a, (p1, p2, p3, p4)))
    assert sp.cancel(delta-expected) == 0
    assert sp.simplify(2*coefficient*centers[4][1]-4*(T*p1*p3-S*p2*p4)) == 0
    return centers, coefficient, (p1, p2, p3, p4, S, T)


def parameter_integral(z, parameters, eps, power, seed):
    """Regulated projective integral with a4 rotated into the upper quadrant."""
    p1, p2, p3, p4, S, T = map(float, parameters)
    u, v, w, h = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power).T
    s = -2*np.log(u)
    k = 1/np.sqrt(1-h)-1
    total = 0j
    for log_y in (-s, s):
        y = np.exp(log_y)
        log_ratio = log_y-log(z)
        conditional = np.exp((1-v)*log_ratio)
        near = abs(log_ratio) < 1e-12
        safe_log = np.where(near, 1., log_ratio)
        x = np.expm1(v*safe_log)*conditional/np.expm1((1-v)*safe_log)
        x = np.where(near, v/(1-v), x)
        # Jacobian-folded weight after mapping dy on either side of y=1
        # and dx conditional on log(y/z).
        quotient = np.divide(log_ratio, y-z, out=np.full_like(y, 1/z), where=~near)
        weight = 4*(1+z)*y*quotient/(u*(1+y))
        a1, a2 = p3*x/S, -p4*y/T
        C = p1*a1+p2*a2
        A = C/(1-w)
        a3 = w*A/(p3+S*a1)
        a4 = ROTATION*A*k/(-p4+T*a2)
        P = 1+a1+a2+a3+a4
        denominator = 1-ROTATION*k-1j*eps*P**2/A
        total += np.mean(-ROTATION*weight*((1+k)/denominator)**3)
    return total


def boundary_integral(z):
    """Integrate the two remaining projective parameters at zero regulator."""
    def density(y):
        relative = (y-z)/z
        quotient = np.log1p(relative)/(y-z) if abs(relative) < 0.5 and y != z else (
            1/z if y == z else log(y/z)/(y-z)
        )
        return 2*(1+z)*quotient/(1+y)

    knots = (0, *sorted({1, z}), np.inf)
    return sum(quad(density, lo, hi, epsabs=1e-10, epsrel=1e-10)[0]
               for lo, hi in zip(knots, knots[1:]))


def check_finite_cross_ratio():
    symbol = sp.symbols("z", positive=True)
    frame, affine_coefficient, projective = frame_parameters(symbol)
    for z in (sp.Rational(1, 2), sp.Rational(3, 2), sp.Integer(2)):
        exact = [tuple(sp.sympify(value).subs(symbol, z) for value in center) for center in frame]
        coefficient = float(affine_coefficient.subs(symbol, z))
        parameters = tuple(value.subs(symbol, z) for value in projective)
        centers = np.asarray(exact, dtype=float)
        for energy, *spatial in exact:
            assert SPATIAL_CUTOFF-sp.sqrt(sum(a*a for a in spatial)) > abs(energy)
        expected = boundary_integral(float(z))
        np.testing.assert_allclose(expected, reduced_cd(float(z)), rtol=0, atol=1e-8)
        np.testing.assert_allclose(expected, pi*pi+log(float(z))**2, rtol=0, atol=1e-8)
        print(f"z={z}: reduced boundary quadrature={expected:.12f}")

        eps = 0.3
        tail = tail_quadrature(SPATIAL_CUTOFF, exact, 32, eps, coefficient)
        spatial, spread = estimate(
            lambda seed: spatial_integral(eps, 19, seed, centers, coefficient),
            (23, 37, 61),
        )
        parameter, parameter_spread = estimate(
            lambda seed: parameter_integral(float(z), parameters, eps, 18, seed),
            (13, 21, 36, 47),
        )
        mismatch = abs(spatial+tail-parameter)
        print(f"z={z}, eps={eps}: spatial={spatial+tail:.9f}, parameter={parameter:.9f}, "
              f"mismatch={mismatch:.3e}, scramble spreads={spread:.3e}/{parameter_spread:.3e}")
        assert max(mismatch, spread, parameter_spread) < 0.002

        errors = []
        for eps in (1e-3, 1e-5, 1e-6, 0):
            value, spread = estimate(
                lambda seed: parameter_integral(float(z), parameters, eps, 18, seed),
                (13, 21, 36, 47),
            )
            errors.append(abs(value-expected))
            print(f"z={z}, eps={eps:.0e}: parameter={value:.9f}, "
                  f"boundary difference={errors[-1]:.3e}, scramble spread={spread:.3e}")
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-2] < 0.005 and max(errors[-1], spread) < 1e-4
    print("PASS: full finite-frame integral away from z=1")


if __name__ == "__main__":
    check_finite_cross_ratio()
