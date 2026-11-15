"""Match the fixed-regulator cusp coefficient to the critical kernel."""

from math import log

import numpy as np
import sympy as sp
from scipy.special import expit
from scipy.stats import qmc

from .critical_regulator_kernel import critical_kernel
from .cusp_regulator_scaling import SEEDS
from .full_finite_frame import ROTATION, estimate
from .large_regulator_ratio import compactified_integral


def matched_integral(z, eps, inverse_eta, power, seed):
    """Evaluate eta*F on rho^2=eps*inverse_eta, including boundary kernels."""
    u, v, w, h = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power).T
    t = np.sqrt(u)/2
    log_odds = np.log1p(-t)-np.log(t)
    q = expit((2*v-1)*log_odds)
    x, y = q/t-1, (1-q)/t-1
    P = 1/(2*t)
    D = 1+z+eps*inverse_eta*(z-1)
    p3 = (1+z)/D
    U, V = p3+x, p3+y
    s = w/(1-w)
    alpha = ROTATION.conjugate()*s*h
    beta = ROTATION*s*(1-h)
    source = inverse_eta*(z*x+y)/(D*P**2)
    denominator = source+alpha-beta-1j*(1+eps*P*(alpha/U+beta/V))**2
    # The map for (x,y) samples t and the logit of q. The sum s of the
    # two ray radii gives the extra s in their Jacobian.
    weight = -4*p3**2*log_odds*q*(1-q)/(U*V*t*t)
    return np.mean(weight*s/((1-w)**2*denominator**3))


def certify_matching():
    v = sp.symbols("v", positive=True)
    a, b = sp.symbols("a b", nonnegative=True)
    # When the two auxiliary parameters rotate in opposite directions
    # through at most pi/4, the real part of the squared regulator factor
    # is bounded below by 1 in scaled coordinates.
    assert sp.expand((1+a+b)**2-(a-b)**2) == 1+2*a+2*b+4*a*b

    # At eps=0 the two ray integrations give -1/(2*(v-i)).
    # Here v is the positive source carried into the boundary kernel.
    x, y = sp.symbols("x y", real=True)
    Q1, Q2 = sp.symbols("Q1 Q2", nonzero=True)
    denominator = v+Q1*x-Q2*y-sp.I
    primitive_y = Q1/(2*denominator**2)
    assert sp.cancel(sp.diff(primitive_y, y)-Q1*Q2/denominator**3) == 0
    primitive_x = 1/(2*(v+Q1*x-sp.I))
    assert sp.cancel(sp.diff(primitive_x, x)+Q1/(2*(v+Q1*x-sp.I)**2)) == 0

    # At inverse_eta=0 the remaining weight is 2/((1+x)*(1+y)*P^2).
    t, q = sp.symbols("t q", positive=True)
    x, y = q/t-1, (1-q)/t-1
    jacobian = sp.Matrix([x, y]).jacobian([t, q]).det()
    assert jacobian == t**-3
    assert sp.cancel(2*jacobian/((1+x)*(1+y)*(1+x/2+y/2)**2)
                     -8*t/(q*(1-q))) == 0
    coefficient = 16*sp.integrate(t*sp.log((1-t)/t), (t, 0, sp.Rational(1, 2)))
    assert coefficient == 8*sp.log(2)-4
    return coefficient


def check_matched_regulator_limit():
    first = float(certify_matching())
    analytic = 1j*(8*log(2)-4)
    assert abs(1j*first-analytic) < 1e-14

    def integral(eps, inverse_eta, z=1.5):
        return estimate(lambda seed: matched_integral(z, eps, inverse_eta, 18, seed), SEEDS)

    eps = 0.3
    fixed, spread = integral(eps, 0)
    ordinary, _ = estimate(
        lambda seed: eps*compactified_integral(1.5, 0, eps, 18, seed), SEEDS,
    )
    np.testing.assert_allclose(fixed, ordinary, rtol=0, atol=4e-4)
    assert spread < 0.001

    errors = []
    for eps in (0.3, 0.16, 0.05, 0.02, 0.001):
        value, spread = integral(eps, 0)
        errors.append(abs(value-analytic))
        print(f"eps={eps:g}: eps*fixed-regulator coefficient={value:.9f}, "
              f"difference={errors[-1]:.3e}, spread={spread:.3e}")
        assert spread < 0.001
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 2e-4

    inverse_eta = 0.1
    target = critical_kernel(1.5, 1/inverse_eta)/inverse_eta
    for eps in (0.05, 0.001, 0.0001, 0):
        value, spread = integral(eps, inverse_eta)
        assert spread < 0.001
        print(f"eta=10, eps={eps:g}: eta*F={value:.9f}, difference={abs(value-target):.3e}")
    np.testing.assert_allclose(value, target, rtol=0, atol=2e-4)
    limit, spread = integral(0, 0)
    assert max(abs(limit-analytic), spread) < 2e-4
    for delta in (0.02, 0.001):
        value, spread = integral(delta, delta)
        print(f"rho=eps={delta:g}: eta*F={value:.9f}, spread={spread:.3e}")
        assert spread < 0.001 and abs(value-analytic) < delta+2e-4
    print(f"Common matched limit: {analytic:.12f}")
    print("PASS: matched small-regulator and large-ratio limits")


if __name__ == "__main__":
    check_matched_regulator_limit()
