"""Confine regulated diagonal-cap poles to the lower endpoint layer."""

from math import pi, sin

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_regulated_base_cap import cap_root_coefficients, denominator_on_cap
from .cd_regulated_cap_meridians import RADIUS
from .cd_regulated_upper_cap import UPPER_CUTOFF
from .full_finite_frame import ROTATION


MAX_REGULATOR = 1e-4
LAYER_WIDTH = 2500.
GAP = sin(pi/8)
PHASE = np.exp(-1j*pi/8)
LINEAR_BOUND = 5+2/RADIUS
ERROR_SLOPE = 5*LINEAR_BOUND/2+4*UPPER_CUTOFF


def certify_regulator_expansion():
    x, A, eps, B, G, H, a = sp.symbols("x A eps B G H a", nonzero=True)
    original = H-sp.I*eps*(A*B+x*G/(1+x))**2/(4*x*A)
    dominant = H-sp.I*eps*x*G**2/(4*A*(1+x)**2)
    correction = -sp.I*eps*(B*G/(2*(1+x))+A*B**2/(4*x))
    assert sp.cancel(original-dominant-correction) == 0
    # A proportional to eps leaves a finite regulator term, independent of k.
    scaled = sp.cancel(original.subs(A, eps*a))
    limit = H-sp.I*x*G**2/(4*a*(1+x)**2)
    assert sp.cancel(scaled.subs(eps, 0)-limit) == 0


def check_uniform_confinement():
    error_bound = LINEAR_BOUND**2/(4*LAYER_WIDTH)+MAX_REGULATOR*ERROR_SLOPE
    assert error_bound < GAP/2
    unit = np.exp(2j*pi*np.arange(256)/256)
    k = np.concatenate((unit, .45*unit, [0.]))
    fractions = np.array([1e-10, 1e-6, .02, .5, 1.])[:, None]
    for eps in (MAX_REGULATOR, 1e-6, 1e-9):
        lower = LAYER_WIDTH*eps
        for x in (lower, 2*lower, .5, 2., UPPER_CUTOFF):
            bound = GAP-eps*(LINEAR_BOUND**2/(4*x)+ERROR_SLOPE)
            for s in (0., .25, 1/(1+ROTATION.real), .9, .999):
                A = fractions*(1-s)
                values = denominator_on_cap(k, RADIUS, x, A, s, eps)
                # The relative slack is only for floating-point projection
                # when the edge term makes the denominator very large.
                slack = 5*np.finfo(float).eps*abs(values)
                assert np.all((PHASE*values).imag <= -bound+slack)
                center, slope = cap_root_coefficients(RADIUS, x, A, s)
                offset = slope*np.sqrt(x/eps)
                assert np.all(abs(center+offset) > 1)
                assert np.all(abs(center-offset) > 1)
                for loop in values[:, :len(unit)]:
                    np.testing.assert_allclose(winding_number(loop), 0., atol=1e-12, rtol=0)
    print(f"Outside x < {LAYER_WIDTH:g}*eps: uniform gap >= {GAP-error_bound:.9f}")


def check_auxiliary_edge_balance():
    k = np.array([.4j, -.7, .9])
    for x, s, a in ((.4, .25, .3), (2., .75, 1.2), (12., .9, .25)):
        G = 1+(ROTATION-1)*s
        H = 1-(1+ROTATION)*s
        target = H-1j*x*G**2/(4*a*(1+x)**2)
        assert (PHASE*target).imag < -GAP
        errors = []
        for eps in (1e-4, 1e-6, 1e-8):
            A = eps*a
            assert A < 1-s
            result = denominator_on_cap(k, RADIUS, x, A, s, eps)
            errors.append(np.max(abs(result-target)))
        assert errors[0] > errors[1] > errors[2]
        assert errors[-1] < 1e-6
        print(f"x={x:g}, A/eps={a:g}: finite edge-limit error={errors[-1]:.3e}")

    # Keep one lower-layer control where the two previously found poles persist.
    eps, s, A, chi = 1e-5, .25, .5, .25
    loop = np.exp(2j*pi*np.arange(256)/256)
    values = denominator_on_cap(loop, RADIUS, eps*chi, A, s, eps)
    np.testing.assert_allclose(winding_number(values), 2., atol=1e-12, rtol=0)


def check_cd_regulated_finite_cap():
    certify_regulator_expansion()
    check_uniform_confinement()
    check_auxiliary_edge_balance()
    print("PASS: uniform finite-x confinement and the surviving auxiliary edge limit")


if __name__ == "__main__":
    check_cd_regulated_finite_cap()
