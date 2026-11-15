"""Compare fixed-frame regulator paths in the rho -> 0 cusp limit."""

from math import log, pi

import numpy as np
import sympy as sp

from .energy_contour import finite_frame, finite_numerator, interval
from .finite_cross_ratio import boundary_integral, scaled_parameter_integral
from .full_finite_frame import estimate

SEEDS = (13, 21, 36, 47)


def frame_scaling():
    rho, z = sp.symbols("rho z", positive=True)
    D = 1+z+rho**2*(z-1)
    frame = finite_frame(z, rho)
    p = [sp.factor(-interval(frame[4], center)) for center in frame[:4]]
    scaled_p = (2*z/D, 2/D, (1+z)/D, -(1+z)/D, 2, 2)
    assert all(sp.cancel(a-b) == 0 for a, b in
               zip(p, (rho**2*scaled_p[0], rho**2*scaled_p[1], *scaled_p[2:4])))
    assert sp.simplify(-p[0]*p[2]/(p[1]*p[3])-z) == 0
    loop = sp.symbols("E x y p_z", real=True)
    numerator = sp.cancel(finite_numerator(loop, frame))
    assert sp.cancel(numerator+8*rho*(1+z)/D*loop[1]) == 0
    assert sp.cancel(2*sp.diff(numerator, loop[1])*frame[4][1]
                    -8*(p[0]*p[2]-p[1]*p[3])) == 0

    # In this projective chart the prescribed +i*eps shift becomes
    # -i*eps*(1+sum(a))^2. Substituting a3=rho^2*b3, a4=rho^2*b4
    # leaves the regulator ratio eps/rho^2 after factoring out rho^2.
    a1, a2, b3, b4, eps = sp.symbols("a1 a2 b3 b4 eps", positive=True)
    A = scaled_p[0]*a1+scaled_p[1]*a2+(p[2]+2*a1)*b3
    B = -p[3]+2*a2
    projective_sum = 1+a1+a2+rho**2*(b3+b4)
    delta = p[0]*a1+p[1]*a2+(p[2]+2*a1)*rho**2*b3-B*rho**2*b4
    regulated = (delta-sp.I*eps*projective_sum**2)/rho**2
    assert sp.cancel(regulated-(A-B*b4-sp.I*(eps/rho**2)*projective_sum**2)) == 0
    return rho, z, scaled_p


def check_cusp_regulator_scaling():
    rho_symbol, z_symbol, symbolic_parameters = frame_scaling()
    z = 1.5

    def integral(rho, eta, power=18):
        parameters = [sp.sympify(p).subs({rho_symbol: rho, z_symbol: z})
                      for p in symbolic_parameters]
        return estimate(lambda seed: scaled_parameter_integral(
            z, parameters, eta, rho*rho, power, seed), SEEDS)

    boundary = boundary_integral(z)
    np.testing.assert_allclose(boundary, pi*pi+log(z)**2, rtol=0, atol=1e-9)
    critical_eta = 0.02
    critical, critical_spread = integral(0, critical_eta)
    print(f"z={z:g}: zero-ratio limit={boundary:.9f}, "
          f"eta={critical_eta:g} limiting kernel={critical:.9f}, spread={critical_spread:.3e}")
    assert abs(critical-boundary) > 0.5
    critical_refined, _ = integral(0, critical_eta, power=19)
    np.testing.assert_allclose(critical, critical_refined, rtol=0, atol=0.002)
    zero_ratio, zero_spread = integral(0, 0)
    assert max(abs(zero_ratio-boundary), zero_spread) < 5e-5

    fixed_values, critical_errors, vanishing_errors = [], [], []
    for rho in (0.25, 0.0625, 0.015625, 0.00390625, 1/2048):
        fixed, fixed_spread = integral(rho, 0.001/rho**2)
        along_critical, critical_spread = integral(rho, critical_eta)
        vanishing, vanishing_spread = integral(rho, rho**2)
        fixed_values.append(abs(fixed))
        critical_errors.append(abs(along_critical-critical))
        vanishing_errors.append(abs(vanishing-boundary))
        print(f"rho={rho:.8f}: fixed eps={fixed:.9f}, eps=eta*rho^2={along_critical:.9f}, "
              f"eps=rho^4={vanishing:.9f}")
        assert max(fixed_spread, critical_spread, vanishing_spread) < 0.004

    assert all(next_value < value for value, next_value in zip(fixed_values, fixed_values[1:]))
    assert fixed_values[-1] < 0.001
    assert all(next_value < value for value, next_value in zip(critical_errors, critical_errors[1:]))
    assert critical_errors[-1] < 0.001
    assert all(next_value < value for value, next_value in zip(vanishing_errors, vanishing_errors[1:]))
    assert vanishing_errors[-1] < 0.03
    print("PASS: three fixed-frame regulator paths in the cusp limit")


if __name__ == "__main__":
    check_cusp_regulator_scaling()
