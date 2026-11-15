"""Integrate the complete lower one-meridian band on diagonal slices."""

from math import sqrt

import numpy as np
from scipy.optimize import brentq

from .cd_regulated_cap_meridians import RADIUS, limiting_crossing
from .cd_regulated_meridian_sweep import CENTER, root_direction, swept_roots
from .cd_remote_meridian_period import circle_period, period_on_slice
from .full_finite_frame import ROTATION


SWITCH = 1/(1+ROTATION.real)


def rim_thresholds(A, s, eps):
    slope = root_direction(s)
    limiting = limiting_crossing(CENTER, np.array([slope, -slope]))
    if eps == 0:
        return limiting
    thresholds = []
    for branch, reference in enumerate(limiting):
        def rim_error(lam):
            return abs(swept_roots(lam, A, s, eps)[branch])-1

        crossing = brentq(rim_error, .98*reference, 1.02*reference, xtol=1e-11)
        k = swept_roots(crossing, A, s, eps)[branch]
        x = eps*A*crossing
        G = 1+(ROTATION-1)*s
        velocity = RADIUS/2*eps*(4*A+(G-A)/(1+x)**2)
        velocity += (1-2*branch)*slope/(2*sqrt(crossing))
        assert np.real(k.conjugate()*velocity) > 0
        thresholds.append(crossing)
    return np.array(thresholds)


def band_slice(A, s, eps, nodes, weights):
    """Integrate the surviving residue between the two outward crossings."""
    thresholds = rim_thresholds(A, s, eps)
    lo, hi = sorted(thresholds)
    survivor = np.argmax(thresholds)
    result = 0j
    for node, weight in zip(nodes, weights):
        lam = lo+(hi-lo)*node
        tau = sqrt(lam)
        density, roots, _ = period_on_slice(1., A, s, eps, tau)
        cap_roots = RADIUS*(2.5+tau*roots)
        tracked = swept_roots(lam, A, s, eps)
        assert abs(tracked[survivor]) < 1 < abs(tracked[1-survivor])
        inside = np.argmin(abs(cap_roots-tracked[survivor]))
        np.testing.assert_allclose(cap_roots[inside], tracked[survivor], atol=1e-13, rtol=0)
        root = roots[inside]
        circle = min(abs(root-roots[1-inside]), abs(root+2.5/tau), abs(root+1.5/tau))/4
        result += weight*circle_period(density, root, circle)
    return (hi-lo)*result


def band_coefficient(eps, order):
    """Return the period per unit t, divided by eps, at t=y/x=1."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    result = 0j
    if eps == 0:
        # The limiting slice has no A dependence, so integrate A exactly.
        for lo, hi in ((0., SWITCH), (SWITCH, 1.)):
            for node, weight in zip(nodes, weights):
                s = lo+(hi-lo)*node
                result += (hi-lo)*weight*(1-s)**2/2*band_slice(0., s, eps, nodes, weights)
        return result

    for fraction, weight_a in zip(nodes, weights):
        def threshold_difference(s):
            first, second = rim_thresholds(fraction*(1-s), s, eps)
            return first-second

        switch = brentq(threshold_difference, SWITCH-.02, SWITCH+.02, xtol=1e-12)
        for lo, hi in ((0., switch), (switch, 1.)):
            for node, weight_s in zip(nodes, weights):
                s = lo+(hi-lo)*node
                A = fraction*(1-s)
                # dx=eps*A*d(lambda), and dA=(1-s)*d(fraction) at fixed s.
                result += ((hi-lo)*weight_a*weight_s*A*(1-s)
                           *band_slice(A, s, eps, nodes, weights))
    return result


def check_cd_diagonal_meridian_band():
    limit = band_coefficient(0., 40)
    np.testing.assert_allclose(band_coefficient(0., 28), limit, atol=2e-8, rtol=0)
    np.testing.assert_allclose(limit, -11.242017568-7.416123969j, atol=2e-8, rtol=0)
    print(f"Full diagonal one-meridian band per eps per unit t: {limit:.10f}")
    previous = np.inf
    for eps in (1e-5, 1e-6, 1e-7):
        value = band_coefficient(eps, 14)
        error = abs(value-limit)
        assert error < previous
        previous = error
        print(f"eps={eps:g}: scaled band={value:.10f}, error={error:.3e}")
    assert previous < .005
    np.testing.assert_allclose(band_coefficient(1e-7, 20), value, atol=6e-5, rtol=0)
    print("PASS: full diagonal one-meridian band and its nonzero linear coefficient")


if __name__ == "__main__":
    check_cd_diagonal_meridian_band()
