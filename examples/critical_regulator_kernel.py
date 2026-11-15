"""One-dimensional reduction and small-regulator asymptotics of the cusp kernel."""

from math import log, pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .cusp_regulator_scaling import SEEDS, frame_scaling
from .finite_cross_ratio import boundary_integral, scaled_parameter_integral
from .full_finite_frame import estimate


def folded_density(t, z, eta):
    """Integrate along one coordinate in each half of the unit square."""
    b = z-1
    c = 1-(z+1)*t-1j*eta*(z+1)/(4*t)
    low, high = c+b*t, c+b*(1-t)
    span = 1-2*t
    left = np.log1p(span*c/(t*high))/c
    right = np.log1p(span*(c+b)/(t*low))/(c+b)
    return 2*(1+z)*(left+right)


def critical_kernel(z, eta):
    """Return the fixed-ratio limit, with eta=epsilon/rho^2."""
    if eta == 0:
        return boundary_integral(z)
    # Resolve the endpoint layers at t=0 and t=1/2.
    scales = (0.1, 1, 10, 100, 1000, 10000)
    knots = sorted({0., .25, .5,
                    *(min(.25, s*eta) for s in scales),
                    *(max(.25, .5-s*eta) for s in scales)})
    return sum(quad(folded_density, lo, hi, args=(z, eta), complex_func=True,
                    epsabs=2e-11, epsrel=2e-11)[0]
               for lo, hi in zip(knots, knots[1:]))


def check_reduction():
    z, eta, x, y, b3, b4 = sp.symbols("z eta x y b3 b4", positive=True)
    Q = sp.symbols("Q", nonzero=True)
    a1, a2 = x/2, y/2
    C = (2*z*a1+2*a2)/(1+z)-sp.I*eta*(1+a1+a2)**2
    A = C+(1+x)*b3
    B = 1+y
    # First use the antiderivative on the rotated b4 ray, then on b3>=0.
    first = 1/(2*B*(A-B*Q*b4)**2)
    second = -1/((1+x)*(C+(1+x)*b3))
    assert sp.cancel(sp.diff(first, b4)-Q/(A-B*Q*b4)**3) == 0
    assert sp.cancel(sp.diff(second, b3)-1/A**2) == 0
    # The scaled numerator is -16. The x,y Jacobian contributes 1/4.
    factor = -4
    result = factor*(-1/(2*B))*1/((1+x)*C)
    expected = 2*(1+z)/((1+x)*(1+y)*(z*x+y-sp.I*eta*(1+z)*(1+(x+y)/2)**2))
    assert sp.cancel(result-expected) == 0

    # x+y=s, x/s=v, t=1/(s+2). The radial integration weight becomes
    # 2*(1+z)*(1-2*t)/((t+v*(1-2*t))*(1-t-v*(1-2*t))).
    t, v, q, c, b = sp.symbols("t v q c b", positive=True)
    s = (1-2*t)/t
    transform = {x: s*v, y: s*(1-v)}
    kernel = expected.subs(transform)*s/t**2
    weight = t+v*(1-2*t)
    mapped = 2*(1+z)*(1-2*t)/(weight*(1-weight)*(1-(z+1)*t
                                                +(z-1)*weight-sp.I*eta*(1+z)/(4*t)))
    assert sp.cancel(kernel-mapped) == 0
    # Partial fractions in q=t+v*(1-2*t) give the two logarithms.
    assert sp.cancel(1/(q*(1-q)*(c+b*q))
                     -1/(q*(c+b*q))-1/((1-q)*(c+b*q))) == 0
    assert sp.simplify(sp.diff(-sp.log((c+b*q)/q)/c, q)-1/(q*(c+b*q))) == 0
    assert sp.simplify(sp.diff(sp.log((c+b*q)/(1-q))/(c+b), q)
                       -1/((1-q)*(c+b*q))) == 0


def check_critical_regulator_kernel():
    check_reduction()
    rho_symbol, z_symbol, symbolic_parameters = frame_scaling()
    for z in (0.5, 1., 1.5, 2.):
        parameters = [sp.sympify(p).subs({rho_symbol: 0, z_symbol: z})
                      for p in symbolic_parameters]
        eta = 0.02
        value = critical_kernel(z, eta)
        np.testing.assert_allclose(value, critical_kernel(1/z, eta), rtol=0, atol=1e-10)
        direct, spread = estimate(lambda seed: scaled_parameter_integral(
            z, parameters, eta, 0., 18, seed), SEEDS)
        mismatch = abs(direct-value)
        assert max(mismatch, spread) < 0.002
        print(f"z={z:g}, eta={eta:g}: reduced kernel={value:.12f}, "
              f"four-parameter mismatch={mismatch:.3e}")

    z = 1.5
    boundary = critical_kernel(z, 0)
    coefficient = (1+z)**2*(1+1/z**2)/4
    errors = []
    for eta in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        value = critical_kernel(z, eta)
        logarithm = log(1/eta)
        imag_ratio = value.imag/(eta*logarithm**2)
        real_ratio = (boundary-value.real)/(eta*logarithm)
        errors.append(abs(imag_ratio-coefficient)+abs(real_ratio-pi*coefficient))
        print(f"eta={eta:.0e}: Im/(eta*log^2)={imag_ratio:.9f}, "
              f"Re deficit/(eta*log)={real_ratio:.9f}")
    print(f"Predicted ratios: {coefficient:.9f}, {pi*coefficient:.9f}")
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 1
    print("PASS: reduced critical kernel and its logarithmic regulator asymptotics")


if __name__ == "__main__":
    check_critical_regulator_kernel()
