"""Bound the regulated diagonal CD cap uniformly at the upper x endpoint."""

from math import pi, sin

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_regulated_base_cap import cap_root_coefficients, denominator_on_cap
from .cd_regulated_cap_meridians import RADIUS
from .full_finite_frame import ROTATION


UPPER_CUTOFF = 14.


def certify_upper_scaling():
    """Check the reciprocal chart and the two regulator balances."""
    X, b, k, R, G, H, eps, beta, ell = sp.symbols("X b k R G H eps beta ell")
    # For X=1/x and b=x*A, the regulator numerator on the cap is U.
    U = 4*b+b*X*(5-2*k/R)+(G-b*X)/(1+X)
    multiplier = U**2/(4*b)
    assert sp.cancel(multiplier.subs(X, 0)-(4*b+G)**2/(4*b)) == 0
    # The two balances of eps*(4*b+G)**2/(4*b) at the ends of the face.
    low = sp.cancel((eps*multiplier).subs(b, eps*beta))
    high = sp.cancel((eps*multiplier).subs(b, ell/eps))
    assert sp.cancel(low.subs({X: 0, eps: 0})-G**2/(4*beta)) == 0
    assert sp.cancel(high.subs({X: 0, eps: 0})-4*ell) == 0

    # For real x, decompose U into a positive real term perturbed by k,
    # plus a positive multiple of G. The real coefficient is increasing.
    x, A = sp.symbols("x A", positive=True)
    linear = 5+4*x-x/(1+x)
    direct = A*(5-2*k/R+4*x)+x*(G-A)/(1+x)
    assert sp.cancel(direct-A*(linear-2*k/R)-x*G/(1+x)) == 0
    assert sp.cancel(direct.subs({x: 1/X, A: b*X}, simultaneous=True)-U) == 0
    assert sp.cancel(sp.diff(linear, x)-(4-1/(1+x)**2)) == 0
    assert sp.factor(sp.diff(linear, x)).is_positive


def check_uniform_gap():
    """Check the sector estimate at the cutoff and under parameter refinement."""
    cutoff = UPPER_CUTOFF
    # The allowed sector for U is [-3*pi/16,5*pi/16]. The distance
    # from the real center to its nearer bounding ray exceeds 2/R.
    sector_margin = (5+4*cutoff-cutoff/(1+cutoff))*sin(3*pi/16)-2/RADIUS
    assert sector_margin > 0
    theta = 2*pi*np.arange(256)/256
    rim = np.exp(1j*theta)
    k = np.concatenate((rim, .5*rim, [0.]))
    phase = np.exp(-1j*pi/8)
    gap = sin(pi/8)
    for x in (cutoff, 30., 1e4, 1e8):
        for eps in (1e-8, 1e-4, .3):
            for s in (0., .25, 1/(1+ROTATION.real), .95, .999):
                As = (1-s)*np.array([1e-9, .01, .5, 1.])[:, None]
                values = denominator_on_cap(k, RADIUS, x, As, s, eps)
                assert np.max((phase*values).imag) <= -gap+2e-14
                for loop in values[:, :len(rim)]:
                    np.testing.assert_allclose(winding_number(loop), 0., atol=1e-12, rtol=0)
                centers, slopes = cap_root_coefficients(RADIUS, x, As, s)
                offsets = slopes*np.sqrt(x/eps)
                assert np.min(abs(centers+offsets)) > 1
                assert np.min(abs(centers-offsets)) > 1
    print(f"For x >= {cutoff:g}: denominator gap >= {gap:.9f}, rim winding=0")


def check_upper_corner_limits():
    k = np.array([.5j, -.4, .9])
    s = .25
    G = 1+(ROTATION-1)*s
    H = 1-(1+ROTATION)*s
    eps = .08
    for b in (.25, 2.):
        target = H-1j*eps*(4*b+G)**2/(4*b)
        previous = np.inf
        for x in (1e2, 1e4, 1e6):
            values = denominator_on_cap(k, RADIUS, x, b/x, s, eps)
            error = np.max(abs(values-target))
            assert error < previous
            previous = error
        assert previous < 1e-5
        print(f"Upper ratio x*A={b:g}: fixed-regulator limit error={previous:.3e}")

    beta, ell = .5, .25
    targets = np.array([H-1j*G**2/(4*beta), H-4j*ell])[:, None]
    previous = np.inf
    for eps in (1e-2, 1e-3, 1e-4, 1e-5):
        x = 1/eps**2
        A = np.array([eps*beta, ell/eps])[:, None]/x
        values = denominator_on_cap(k, RADIUS, x, A, s, eps)
        error = np.max(abs(values-targets))
        assert error < previous
        previous = error
    assert previous < 1e-4
    print(f"Two joint regulator balances: final limit error={previous:.3e}")


def check_cd_regulated_upper_cap():
    certify_upper_scaling()
    check_uniform_gap()
    check_upper_corner_limits()
    print("PASS: pole-free regulated diagonal cap at the upper endpoint")


if __name__ == "__main__":
    check_cd_regulated_upper_cap()
