"""Large-eta asymptotics and the distinct fixed-regulator cusp path."""

import numpy as np
import sympy as sp
from scipy.stats import qmc

from .critical_regulator_kernel import critical_kernel
from .cusp_regulator_scaling import SEEDS, frame_scaling
from .finite_cross_ratio import scaled_parameter_integral
from .full_finite_frame import ROTATION, estimate


def critical_coefficients():
    """Coefficients of i/eta, 1/eta^2 and -i/eta^3 in the reduced kernel."""
    z, t, q = sp.symbols("z t q", positive=True)
    g = 1-(1+z)*t+(z-1)*q
    angular = (g.subs(q, 0)**2/q+g.subs(q, 1)**2/(1-q)-(z-1)**2)
    assert sp.cancel(g*g/(q*(1-q))-angular) == 0
    L = sp.log((1-t)/t)
    interval = (t, 0, sp.Rational(1, 2))
    first = 16*sp.integrate(t*L, interval)
    second = 32*sp.integrate(t*t*(1-2*t)*L, interval)
    third = 128/(1+z)**2*sp.integrate(
        t**3*((g.subs(q, 0)**2+g.subs(q, 1)**2)*L-(z-1)**2*(1-2*t)),
        interval,
    )
    assert first == 8*sp.log(2)-4
    assert second == 4-sp.Rational(16, 3)*sp.log(2)
    expected = 8*(-91*z*z+67*z-91+(132*z*z-96*z+132)*sp.log(2))/(45*(1+z)**2)
    assert sp.simplify(third-expected) == 0
    return float(first), float(second), sp.lambdify(z, expected, "numpy")


def compactified_integral(z, rho, eps, power, seed):
    """Compute the finite-frame integral divided by rho^2, including at rho=0."""
    u, v, w, k = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power).T
    cuts = np.sort(np.stack((u, v, w), axis=1), axis=1)
    boundaries = np.column_stack((np.zeros(len(cuts)), cuts, np.ones(len(cuts))))
    b0, b1, b2, b3 = np.diff(boundaries, axis=1).T
    D = 1+z+rho*rho*(z-1)
    p1, p2 = 2*rho*rho*z/D, 2*rho*rho/D
    p3 = (1+z)/D
    A = (p1*b1+p2*b2)*b0+(p3*b0+2*b1)*b3
    B = p3*b0+2*b2
    tail = ROTATION*k/(1-k)
    denominator = A-B*tail-1j*eps*(1+tail)**2
    density = -16*p3*p3*ROTATION*b0/((1-k)**2*denominator**3)
    # The three sorted uniforms give a uniform point of the unit simplex.
    return density.mean()/6


def certify_compactification():
    b1, b2, b3, s = sp.symbols("b1 b2 b3 s", positive=True)
    b0 = 1-b1-b2-b3
    coordinates = sp.Matrix([b1/b0, b2/b0, b3/b0, s/b0])
    assert sp.cancel(coordinates.jacobian([b1, b2, b3, s]).det()-b0**-5) == 0
    assert sp.cancel(1+sum(coordinates)-(1+s)/b0) == 0


def check_large_regulator_ratio():
    first, second, third = critical_coefficients()
    print(f"Universal coefficients: eta*Im -> {first:.12f}, eta^2*Re -> {second:.12f}")
    for z in (0.5, 1., 1.5, 2.):
        C = third(z)
        previous = np.inf
        for eta in (1., 4., 16., 64.):
            value = critical_kernel(z, eta)
            # Geometric-series bounds with |g*t/a| <= 1/(2*eta).
            assert 0 <= second/eta**2-value.real <= second/(4*eta**4)
            imag_correction = eta**3*(first/eta-value.imag)
            assert -1e-9 <= C-imag_correction <= C/(4*eta**2)+1e-9
            error = abs(imag_correction-C)
            assert error < previous
            previous = error
        print(f"z={z:g}: cubic coefficient={C:.9f}, eta=64 estimate={imag_correction:.9f}")

    certify_compactification()
    rho_symbol, z_symbol, projective = frame_scaling()
    eps, z = 0.3, 1.5

    def quotient(rho, z):
        return estimate(lambda seed: compactified_integral(z, rho, eps, 18, seed), SEEDS)

    # Cross-check this coordinate map against the earlier finite-rho evaluator.
    rho = 0.5
    parameters = [sp.sympify(p).subs({rho_symbol: rho, z_symbol: z}) for p in projective]
    direct, _ = estimate(lambda seed: scaled_parameter_integral(
        z, parameters, eps/rho**2, rho**2, 18, seed), SEEDS)
    mapped, spread = quotient(rho, z)
    np.testing.assert_allclose(mapped*rho**2, direct, rtol=0, atol=4e-4)
    assert spread < 0.002

    limit, spread = quotient(0, z)
    assert spread < 0.001
    np.testing.assert_allclose(limit, quotient(0, 0.5)[0], rtol=0, atol=1e-12)
    print(f"Fixed eps={eps:g}: rho^-2 scaled limit={limit:.9f}, scramble spread={spread:.3e}")
    errors = []
    for rho in (1/4, 1/16, 1/64, 1/256):
        value, _ = quotient(rho, z)
        errors.append(abs(value-limit))
        print(f"rho={rho:.8f}: fixed-eps integral/rho^2={value:.9f}")
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 5e-5
    formal_critical_coefficient = 1j*first/eps
    assert abs(limit-formal_critical_coefficient) > 2
    print(f"Critical-kernel coefficient at the same eps: {formal_critical_coefficient:.9f}")
    print("PASS: large-ratio expansion and fixed-regulator cusp scaling")


if __name__ == "__main__":
    check_large_regulator_ratio()
