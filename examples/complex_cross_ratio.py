"""Continue the finite-frame projective integral across the principal sheet."""

from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.stats import qmc

from .cd_integral import reduced_cd
from .finite_cross_ratio import frame_parameters
from .full_finite_frame import ROTATION, estimate


def check_projective_pullback():
    z = sp.symbols("z", positive=True)
    _, _, (p1, p2, p3, p4, S, T) = frame_parameters(z)
    x, y, w, k = sp.symbols("x y w k", positive=True)
    Q = sp.symbols("Q", nonzero=True)
    a1, a2 = p3*x/S, -p4*y/T
    C = sp.factor(p1*a1+p2*a2)
    A = C/(1-w)
    B = -p4+T*a2
    a3 = w*A/(p3+S*a1)
    a4 = Q*A*k/B
    # The coordinate map is triangular in (x,y,w,k).
    variables = (x, y, w, k)
    coordinates = (a1, a2, a3, a4)
    jacobian = sp.prod(sp.diff(ai, v) for ai, v in zip(coordinates, variables))
    delta = C+(p3+S*a1)*a3-B*a4
    numerator = -4*(T*p1*p3-S*p2*p4)
    assert sp.cancel(delta-A*(1-Q*k)) == 0
    result = (2*(1+z)/((1+x)*(1+y)*(z*x+y)))*(-2*Q)/(1-Q*k)**3
    assert sp.cancel(numerator*jacobian/delta**3-result) == 0


def projective_integral(z, power, seed):
    """Unreduced projective quadrature on the transported four-parameter domain."""
    z_abs = abs(z)
    u, v, _, t = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power).T
    s = -2*np.log(u)
    # t=k/(1+k) resolves the infinite end of the a4 ray.
    # The w Jacobian cancels the dependence on w in the transformed density.
    ray_weight = -2*ROTATION*(1-t)/(1-t-ROTATION*t)**3
    result = 0j
    for log_y in (-s, s):
        y = np.exp(log_y)
        log_ratio = log_y-np.log(z_abs)
        near = abs(log_ratio) < 1e-12
        safe = np.where(near, 1., log_ratio)
        x = np.expm1(v*safe)*np.exp((1-v)*log_ratio)/np.expm1((1-v)*safe)
        x = np.where(near, v/(1-v), x)
        quotient = np.divide(log_ratio, y-z_abs, out=np.full_like(y, 1/z_abs), where=~near)
        jacobian_weight = 4*(1+z)*y*quotient/(u*(1+y))
        result += np.mean(jacobian_weight*(z_abs*x+y)/(z*x+y)*ray_weight)
    return result


def integrated_projective(z, log_z=None):
    """One-dimensional quadrature after integrating x and the two auxiliary parameters."""
    if log_z is None:
        log_z = np.log(complex(z))

    def density(y):
        offset = (y-z)/z
        if offset == 0:
            quotient = 1/z
        elif abs(offset) < 0.5:
            quotient = np.log1p(offset)/(y-z)
        else:
            quotient = (np.log(y)-log_z)/(y-z)
        return 2*(1+z)*quotient/(1+y)

    scale = abs(z)
    knots = (0, *sorted({1, scale}), np.inf)
    return sum(quad(density, lo, hi, complex_func=True, epsabs=1e-10,
                    epsrel=1e-10)[0] for lo, hi in zip(knots, knots[1:]))


def check_complex_cross_ratio():
    check_projective_pullback()
    points = (0.35+0.8j, -0.65+0.4j, -1.5-0.35j, 1j)
    for z in points:
        expected = pi*pi+np.log(z)**2
        reduced = integrated_projective(z)
        np.testing.assert_allclose(reduced, expected, rtol=0, atol=1e-9)
        np.testing.assert_allclose(reduced, reduced_cd(z), rtol=0, atol=1e-9)
        value, spread = estimate(lambda seed: projective_integral(z, 19, seed), (13, 21, 36, 47))
        mismatch = abs(value-reduced)
        assert max(mismatch, spread) < 0.002
        print(f"z={z}: F={reduced:.9f}, transported integral mismatch={mismatch:.3e}")

    for a in (0.4, 1., 2.5):
        above, below = (integrated_projective(-a, np.log(a)+side*1j*pi)
                        for side in (1, -1))
        np.testing.assert_allclose(above, np.log(a)**2+2j*pi*np.log(a), rtol=0, atol=1e-9)
        np.testing.assert_allclose(below, above.conjugate(), rtol=0, atol=1e-9)
        jump = above-below
        np.testing.assert_allclose(jump, 4j*pi*np.log(a), rtol=0, atol=1e-9)
        # In the two-parameter integral the jump comes from y=a*x.
        # The imaginary part of z*x+y has the sign of Im(z) for x>0.
        on_pole = quad(lambda x: 1/((1+x)*(1+a*x)), 0, np.inf, epsabs=1e-11)[0]
        np.testing.assert_allclose(jump, -4j*pi*(1-a)*on_pole, rtol=0, atol=1e-9)

        errors = []
        for eta in (1e-2, 1e-4, 1e-6):
            values = [integrated_projective(complex(-a, side*eta)) for side in (1, -1)]
            errors.append(max(abs(values[0]-above), abs(values[1]-below)))
        assert errors[2] < errors[1] < errors[0] and errors[-1] < 2e-5
        print(f"z=-{a:g}: cut jump={jump:.9f}, last approach error={errors[-1]:.3e}")
    print("PASS: transported projective integral and both sides of the principal cut")


if __name__ == "__main__":
    check_complex_cross_ratio()
