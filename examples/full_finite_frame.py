"""Compute the full finite-frame integral by two regulated representations."""

from itertools import combinations
from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.stats import qmc

from .cd_integral import reduced_cd
from .energy_contour import finite_frame, finite_numerator, interval, upper_residues
from .spatial_tail import tail_quadrature

SPATIAL_CUTOFF = 3
ROTATION = np.exp(1j*pi/4)


def certify_parameter_reduction(centers):
    """Verify the projective denominator and shifted numerator symbolically."""
    a1, a2, a3, a4 = weights = sp.symbols("a1 a2 a3 a4", positive=True)
    projective = (*weights, 1)
    P = sum(projective)
    centroid = tuple(sum(a*p[k] for a, p in zip(projective, centers))/P for k in range(4))
    loop = sp.symbols("E x y z", real=True)
    x = loop[1]
    assert sp.expand(finite_numerator(loop, centers)) == -4*x
    assert sp.simplify(finite_numerator(centroid, centers)-2/P) == 0
    delta = -sum(projective[i]*projective[j]*interval(centers[i], centers[j])
                 for i, j in combinations(range(5), 2))
    expected = (a1+a2)/4+(1+2*a1)*a3-(1+2*a2)*a4
    assert sp.expand(delta-expected) == 0

    t, eps, c, b, p, B = sp.symbols("t eps c b p B", positive=True)
    zz = c+sp.I*b
    imaginary_part = sp.im(t-B*zz-sp.I*eps*(p+zz)**2)
    assert sp.expand(imaginary_part+B*b+eps*((p+c)**2-b*b)) == 0
    real_part = sp.re(t-B*zz-sp.I*eps*(p+zz)**2)
    assert sp.expand(real_part-t+B*c-2*eps*b*(p+c)) == 0
    # In the sector 0 <= Im(a4) <= Re(a4) there is no zero.
    # At eps=0, the analytic a4 integral on the rotated ray is -1/(2*B*A^2).
    A, B = sp.symbols("A B", positive=True)
    a4_primitive = 1/(2*B*(A-B*a4)**2)
    assert sp.simplify(sp.diff(a4_primitive, a4)-(A-B*a4)**-3) == 0
    assert sp.limit(a4_primitive, a4, sp.oo) - a4_primitive.subs(a4, 0) == -1/(2*B*A**2)


def boundary_integral():
    """Zero-regulator limit after the rotated a4 and positive a3 integrals."""
    # The two remaining parameters give 4*dx*dy/((1+x)*(1+y)*(x+y)).
    # Integrate x and fold y>1 back by reciprocal substitution.
    return 8*quad(lambda y: -np.log(y)/(1-y*y), 0, 1, epsabs=1e-11)[0]


def parameter_integral(eps, power, seed):
    """Four-dimensional quadrature with a4 rotated through +pi/4."""
    u, v, w, h = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power).T
    s = -2*np.log(u)
    prefactor = 8*u*s/(-np.expm1(-2*s))
    total = 0j
    for log_y in (-s, s):
        y = np.exp(log_y)
        # Conditional map for dx/((1+x)*(x+y)) after integrating at eps=0.
        x = (np.expm1(v*log_y)*np.exp((1-v)*log_y)
             / np.expm1((1-v)*log_y))
        a1, a2 = x/2, y/2
        C = (a1+a2)/4
        A = C/(1-w)
        a3 = w*A/(1+2*a1)
        B = 1+2*a2
        k = 1/np.sqrt(1-h)-1
        a4 = ROTATION*A*k/B
        P = 1+a1+a2+a3+a4
        denominator = 1-ROTATION*k-1j*eps*P**2/A
        density = -ROTATION*prefactor*((1+k)/denominator)**3
        total += density.mean()
    return total


def spatial_integral(eps, power, seed, centers, x_coefficient=-4):
    """Spatial ball integral of the normalized upper energy residue sum."""
    radius = SPATIAL_CUTOFF
    u, v, w = qmc.Sobol(3, scramble=True, seed=seed).random_base2(power).T
    r = radius*np.cbrt(u)
    mu = 2*v-1
    phi = pi*w
    transverse = r*np.sqrt(1-mu**2)
    points = np.stack((transverse*np.cos(phi), transverse*np.sin(phi), r*mu), axis=-1)
    reflected = points*np.array([-1, 1, 1])
    total = 0j
    chunk_size = 2**15
    for start in range(0, len(points), chunk_size):
        chunk = slice(start, start+chunk_size)
        total += upper_residues(points[chunk], centers, eps, x_coefficient).sum()
        total += upper_residues(reflected[chunk], centers, eps, x_coefficient).sum()
    # 2/pi from upper closure and the loop normalization, times the
    # half-ball volume (the other half is included by reflection).
    return 4*radius**3/3*total/len(points)


def estimate(method, seeds):
    samples = np.array([method(seed) for seed in seeds])
    mean = samples.mean()
    return mean, np.max(abs(samples-mean))


def check_full_finite_frame():
    exact = finite_frame()
    certify_parameter_reduction(exact)
    centers = np.asarray(exact, dtype=float)
    for energy, *spatial in exact:
        assert SPATIAL_CUTOFF-sp.sqrt(sum(a*a for a in spatial)) > abs(energy)

    # At these regulators, the ball residue integral is sufficiently smooth
    # for a direct numerical comparison with the independent parameter map.
    for eps, tolerance in ((0.3, 0.001), (0.16, 0.004)):
        tail = tail_quadrature(SPATIAL_CUTOFF, exact, 32, eps)
        np.testing.assert_allclose(tail_quadrature(SPATIAL_CUTOFF, exact, 24, eps),
                                   tail, rtol=0, atol=5e-7)
        spatial, spread = estimate(lambda seed: spatial_integral(eps, 19, seed, centers),
                                   (23, 37, 61))
        parameter, parameter_spread = estimate(lambda seed: parameter_integral(eps, 18, seed),
                                               (13, 21, 36, 47))
        mismatch = abs(spatial+tail-parameter)
        print(f"eps={eps:g}: spatial={spatial+tail:.9f}, parameter={parameter:.9f}, "
              f"mismatch={mismatch:.3e}, scramble spreads={spread:.3e}/{parameter_spread:.3e}")
        assert max(mismatch, spread, parameter_spread) < tolerance

    expected = boundary_integral()
    print(f"Analytically reduced boundary quadrature: {expected:.12f}")
    np.testing.assert_allclose(expected, reduced_cd(1), rtol=0, atol=1e-9)
    errors = []
    for eps in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0):
        value, spread = estimate(lambda seed: parameter_integral(eps, 18, seed),
                                 (13, 21, 36, 47))
        errors.append(abs(value-expected))
        print(f"eps={eps:.0e}: parameter={value:.9f}, "
              f"benchmark difference={errors[-1]:.3e}, scramble spread={spread:.3e}")
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-2] < 0.003 and max(errors[-1], spread) < 5e-5
    print("PASS: full finite-frame comparison and the z=1 benchmark")


if __name__ == "__main__":
    check_full_finite_frame()
