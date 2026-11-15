"""Evaluate paired regulated meridians in the transverse pole-scale chart."""

from math import pi

import numpy as np
import sympy as sp

from .cd_frame_feynman_overlap import frame_chart
from .full_finite_frame import ROTATION


def certify_rational_form():
    u, tau, Q, H, J, K, P, x, t, h = sp.symbols("u tau Q H J K P x t h")
    D = H+sp.I*(u+h*P)**2/J
    numerator = H*J+sp.I*(u+h*P)**2
    previous = 4*Q*tau/((1+x)*(1+t*x)*K*(1-K)*J*D**3)
    reduced = 4*Q*tau*J**2/((1+x)*(1+t*x)*K*(1-K)*numerator**3)
    assert sp.cancel(previous-reduced) == 0

    # The leading pair has equal prefactors and opposite residues.
    q, root, factor = sp.symbols("q root factor", nonzero=True)
    residue = sp.diff(factor/(q+2*root)**3, q, 2).subs(q, 0)/2
    assert sp.cancel(residue-6*factor/(2*root)**5) == 0
    assert sp.cancel(residue+residue.subs(root, -root)) == 0


def period_on_slice(t, A, s, eps, tau):
    """Return a slice density, its two local poles and its paired period."""
    Q = ROTATION
    x = eps*A*tau**2
    h = eps*tau
    B = 1-A-s
    H = 1-(1+Q)*s
    j0 = (3-5*t)/2
    j1 = (1-t)*tau
    g = B/(2*(1+x))+Q*s/(2*(1+t*x))
    p0 = -A*(1+t)+g*j0
    p1 = g*j1
    m0, m1 = h*p0, 1+h*p1
    polynomial = np.array([1j*m1**2, H*j1+2j*m0*m1, H*j0+1j*m0**2])
    roots = np.roots(polynomial)
    prefactor = 4*Q*tau/((1+x)*(1+t*x))

    def density(u):
        K = 2.5+tau*u
        J = j0+j1*u
        N = np.polyval(polynomial, u)
        return prefactor*J**2/(K*(1-K)*N**3)

    # The only other poles are at K=0,1. The differential has no
    # residue at infinity since its coefficient decays as u**-6.
    base = np.array([0., 1.])
    remote = (base-2.5)/tau
    N = np.polyval(polynomial, remote)
    residues = prefactor*np.array([1., t*t])/(tau*(1-2*base)*N**3)
    return density, roots, -residues.sum()


def circle_period(density, center, radius):
    """Integrate density du/(2*pi*i) counterclockwise."""
    normal = radius*np.exp(2j*pi*np.arange(512)/512)
    return np.mean(normal*density(center+normal))


def check_paired_period():
    for t, A, s in ((.8, .4, .25), (1., .1, .6), (1.2, .04, .9)):
        eps = .08
        Q = ROTATION
        H = 1-(1+Q)*s
        j0 = (3-5*t)/2
        target_root = np.sqrt(1j*H*j0)
        leading = 6*4*Q/((2.5)*(1-2.5)*j0*(1j/j0)**3*(2*target_root)**5)
        pair_target = 4j*Q*(t*t/1.5**6-1/2.5**6)
        errors = []
        for tau in (.18, .09, .045, .0225, .01125):
            density, roots, pair = period_on_slice(t, A, s, eps, tau)
            roots = roots[np.argsort(abs(roots-target_root))]
            assert np.max(abs(roots)) < 1.5 and np.min(abs(roots)) > .3
            assert abs(roots[0]-roots[1]) > .36
            single = [circle_period(density, root, .18) for root in roots]
            outer = circle_period(density, 0., 2.)
            np.testing.assert_allclose(outer, pair, atol=3e-14, rtol=2e-5)
            np.testing.assert_allclose(sum(single), pair, atol=3e-14, rtol=2e-5)
            errors.append(abs(pair/tau**6-pair_target))

            # Check the complete Jacobian against the bounded-frame formula
            # continued algebraically to w=u/(eps*tau), at fixed x,A,s.
            u = np.array([.3j, -.2, 1.+.7j])
            x, h = eps*A*tau**2, eps*tau
            z, J, _, D = frame_chart(x, t, u/h, A, s, eps)
            K = 1/(1+z)
            benchmark = 4*Q*(x/A)/((1+x)*(1+t*x)*K*(1-K)*J*D**3*h)
            np.testing.assert_allclose(density(u), benchmark, atol=2e-14, rtol=1e-13)

        assert np.all(np.diff(errors) < 0)
        np.testing.assert_allclose(np.array(single)/tau, [leading, -leading], rtol=.09)
        np.testing.assert_allclose(pair/tau**6, pair_target, atol=2e-4, rtol=0)
        print(f"y/x={t:g}: paired period/tau^6={pair/tau**6:.9f}, "
              f"limit={pair_target:.9f}")


def check_cd_remote_meridian_period():
    certify_rational_form()
    check_paired_period()
    print("PASS: opposite leading meridians and the sixth-order paired period")


if __name__ == "__main__":
    check_cd_remote_meridian_period()
