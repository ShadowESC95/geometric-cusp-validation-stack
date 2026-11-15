"""Finite CD parameter integral, its contour average, and regulator checks."""

from cmath import exp
from math import asinh, inf, log, pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad


def density(r, z):
    """Density after combining the finite five-propagator integral."""
    if r == 0:
        return 2 * (1 + 1 / z)
    logs = np.log1p(r / z) + np.log1p(r)
    return 2 * (1 + z) * logs / (r * (1 + z + r))


def reduced_cd(z):
    return sum(
        quad(density, lo, hi, args=(z,), complex_func=True,
             epsabs=1e-10, epsrel=1e-10, limit=250)[0]
        for lo, hi in ((0, 1), (1, inf))
    )


def paired_parameter_cd(z, tol):
    """Evaluate the earlier paired-parameter quadrature on the positive axis."""
    a = sqrt(z)

    def outer(t):
        def inner(u):
            b = 4 * u * (1 - u)
            return asinh(t / sqrt(b)) / sqrt(b + t * t)

        uv = quad(inner, 0, 0.5, epsabs=tol / 10,
                  epsrel=tol / 10, limit=250)[0]
        return uv * (1 / (1 + a * t)**2 + 1 / (a + t)**2)

    value = quad(outer, 0, 1, epsabs=tol / 10,
                 epsrel=tol / 10, limit=250)[0]
    return 16 * (1 + z) * value


def cd_tail(z, cutoff):
    return quad(lambda t: cutoff * density(cutoff * t, z), 1, inf,
                complex_func=True, epsabs=1e-12, epsrel=1e-10, limit=250)[0]


def check_cd_benchmark():
    for z in (0.5, 1., 2.):
        coarse = paired_parameter_cd(z, 1e-7)
        fine = paired_parameter_cd(z, 1e-10)
        np.testing.assert_allclose(coarse, fine, rtol=0, atol=1e-6)
        np.testing.assert_allclose(fine, pi**2 + log(z)**2, rtol=0, atol=1e-6)
        np.testing.assert_allclose(reduced_cd(z), fine, rtol=0, atol=1e-7)
        print(f"z={z:g}: F1={fine:.12f}")

    for rho in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        near = {}
        tails = []
        cutoff = 0.1 / rho
        for z in (0.5, 1., 2.):
            near[z] = quad(density, 0, cutoff, args=(z,),
                           points=[1] if cutoff > 1 else None,
                           epsabs=1e-10, epsrel=1e-10)[0]
            tail = cd_tail(z, cutoff).real
            tails.append(tail)
            np.testing.assert_allclose(near[z] + tail, pi**2 + log(z)**2,
                                       rtol=0, atol=1e-7)
        increment = near[2.] - near[1.]
        print(f"rho={rho:.0e}: near increment={increment:.9f}, max tail={max(tails):.3e}")
    assert abs(increment - log(2)**2) < 5e-4
    assert max(tails) < 2e-3


def check_contour_and_tail():
    errors = []
    for n in (32, 128, 512):
        z = np.exp(1j * (np.arange(n) + 0.5) * pi / n)
        weights = np.array([reduced_cd(point) for point in z])
        np.testing.assert_allclose(weights.imag, 0, atol=1e-9)
        coefficient = -2 * weights.real.mean()
        errors.append(abs(coefficient + 4*pi**2/3))
        print(f"n={n}: Gamma_cusp[2]={coefficient:.12f}, error={errors[-1]:.3e}")
    assert errors[2] < errors[1] < errors[0] and errors[2] < 1e-5

    z = np.r_[1 + 0j, np.exp(1j * (np.arange(64) + 0.5) * pi / 64)]
    for cutoff in (100, 10_000, 1_000_000):
        bound = 4 * (2*log(cutoff) + 2 + pi + log(3)) / cutoff
        largest = max(abs(cd_tail(point, cutoff)) for point in z)
        assert largest <= bound
        print(f"R={cutoff:,}: max parameter tail={largest:.6e}, bound={bound:.6e}")


def check_numerator_reduction():
    a = sp.Matrix(sp.symbols("a1:5", positive=True))
    q, z = sp.symbols("q z", positive=True)
    S, T = 2, 2*q
    p = sp.Matrix([z/4, sp.Rational(1, 4), 1, q])
    delta = S*a[0]*a[2] + T*a[1]*a[3] + p.dot(a)
    c = sp.Matrix([T*p[2], S*p[3], T*p[0], S*p[1]])
    prefactor = 4*(T*p[0]*p[2] + S*p[1]*p[3])
    flow = (c + S*T*a) / delta**2
    divergence = sum(sp.diff(flow[i], a[i]) for i in range(4))
    assert sp.cancel(divergence + prefactor/delta**3) == 0
    assert sp.cancel(prefactor/(2*q) - (1 + z)) == 0


def alpha4_tail(z, endpoint):
    """Omitted parameter integral at alpha4=endpoint/q after analytic reductions.

    Integrating alpha3 and then alpha4 from the endpoint to infinity leaves
    a two-parameter integral. Integrating x=2*alpha1 gives
    log(c/z)/(c-z), with c=y+8*endpoint*(1+y) and y=2*alpha2.
    """
    def outer(y):
        c = y + 8 * endpoint * (1 + y)
        relative = (c - z) / z
        if abs(relative) < 0.5:
            quotient = np.log1p(relative) / (z * relative) if relative else 1 / z
        else:
            quotient = np.log(c / z) / (c - z)
        return 2 * (1 + z) * quotient / (1 + y)

    return quad(outer, 0, inf, complex_func=True,
                epsabs=1e-11, epsrel=1e-10, limit=250)[0]


def finite_alpha4_interval(z, endpoint):
    """Integrate the finite interval after eliminating alpha3 and alpha4."""
    def outer(y):
        def inner(x):
            base = z*x + y
            shift = 8*endpoint*(1+y)
            return shift / ((1+x)*base*(base+shift))

        return quad(inner, 0, inf, complex_func=True, epsabs=1e-10,
                    epsrel=1e-10, limit=200)[0] / (1+y)

    return 2*(1+z)*quad(outer, 0, inf, complex_func=True, epsabs=1e-9,
                        epsrel=1e-9, limit=200)[0]


def check_continuation():
    check_numerator_reduction()
    for z in (0.5, 1., 2.):
        full = reduced_cd(z)
        np.testing.assert_allclose(full, pi**2 + log(z)**2, rtol=0, atol=1e-9)
        for phi in (0, pi/2, pi - 1e-3, pi - 1e-6):
            q = exp(-1j * phi)
            # alpha4=t/q rotates the cutoff endpoint along with the ray.
            errors = []
            for cutoff in (1., 100., 10_000.):
                transported = alpha4_tail(z, q * (cutoff/q))
                fixed = alpha4_tail(z, q * cutoff)
                np.testing.assert_allclose(transported, alpha4_tail(z, cutoff),
                                           rtol=0, atol=1e-10)
                # With c=y+8*q*cutoff*(1+y), use |c| >= (8*cutoff-1)*(1+y).
                # The upper bound on |log(c/z)| adds pi for its argument.
                bound = (2*(1+z)*(log((8*cutoff+1)/z) + pi + 1)
                         / (8*cutoff-1-z))
                error = max(abs(transported), abs(fixed))
                assert error <= bound
                errors.append(error)
            assert errors[2] < errors[1] < errors[0] and errors[2] < 0.002
            if z == 1:
                print(f"phi={phi:.9f}: full={full.real:.12f}, "
                      f"fixed-cutoff tail={abs(fixed):.3e}")
    # Check the subtraction against direct finite-interval quadrature
    # at two phases away from the propagator singularities.
    for endpoint in (1., -1j):
        direct = finite_alpha4_interval(1., endpoint)
        np.testing.assert_allclose(direct + alpha4_tail(1., endpoint),
                                   reduced_cd(1.), rtol=0, atol=1e-8)
    print("PASS: both parameter cutoff paths approach the continued integral")


if __name__ == "__main__":
    check_cd_benchmark()
    check_contour_and_tail()
    check_continuation()
