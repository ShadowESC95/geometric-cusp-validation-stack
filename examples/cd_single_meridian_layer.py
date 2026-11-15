"""Integrate a compact one-meridian patch in the lower regulator layer."""

from math import prod, sqrt

import numpy as np

from .cd_regulated_cap_meridians import RADIUS
from .cd_remote_meridian_period import circle_period, period_on_slice
from .full_finite_frame import ROTATION


# The patch keeps strictly between the two rim crossings. Its order in
# the five-form is (lambda,t,u,A,s) with x=eps*A*lambda.
PATCH = ((455., 485.), (.999, 1.001), (.25, .35), (.268, .272))
MAX_REGULATOR = 1e-5


def certify_patch_clearance():
    """Bound the pole displacement from the reference slice by Rouché's theorem."""
    (lo, hi), (t_lo, t_hi), (_, A_hi), (s_lo, s_hi) = PATCH
    tau_max = sqrt(hi)
    H0 = 1-(1+ROTATION)*.27
    root = np.sqrt(-1j*H0)
    assert root.real > 0 and root.imag < 0
    radius = .025
    u_bound = abs(root)+radius
    shift_k = RADIUS*tau_max*radius
    for sign in (1, -1):
        direction = sign*root
        # |5/2+tau*direction| is increasing on the interval.
        assert 2.5*direction.real+sqrt(lo)*abs(direction)**2 > 0
    assert abs(RADIUS*(2.5+sqrt(lo)*root))-shift_k > 1
    assert abs(RADIUS*(2.5-sqrt(hi)*root))+shift_k < 1
    assert abs(root.imag) > radius  # avoids all real base punctures

    # On the two reference circles compare N=H*J+i*(u+h*P0)^2
    # with N0=-H0+i*u^2. The perturbation bound holds on the full box.
    t_error = max(1-t_lo, t_hi-1)
    s_error = max(.27-s_lo, s_hi-.27)
    J_error = t_error*(2.5+tau_max*u_bound)
    H_error = abs(1+ROTATION)*s_error
    P_bound = A_hi*(1+t_hi)+(1+J_error)/2
    shift_u = MAX_REGULATOR*tau_max*P_bound
    N_error = H_error+(abs(H0)+H_error)*J_error+shift_u*(2*u_bound+shift_u)
    assert N_error < radius*(2*abs(root)-radius)
    print("Uniform pole clearance on the compact one-meridian patch")
    return -root, radius


def patch_coefficient(eps, order, reference_circle=None):
    """Compute the counterclockwise period divided by eps and patch volume."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    rules = [zip(lo+(hi-lo)*(1+nodes)/2, weights/2) for lo, hi in PATCH]
    lambdas, ratios, auxiliaries, phases = [tuple(rule) for rule in rules]
    integral = 0j
    for lam, wl in lambdas:
        tau = sqrt(lam)
        for t, wt in ratios:
            for A, wa in auxiliaries:
                for s, ws in phases:
                    density, roots, _ = period_on_slice(t, A, s, eps, tau)
                    inside = abs(RADIUS*(2.5+tau*roots)) < 1
                    assert np.count_nonzero(inside) == 1
                    root, other = roots[inside][0], roots[~inside][0]
                    circle = min(abs(root-other), abs(root+2.5/tau), abs(root+1.5/tau))/4
                    center, radius = reference_circle or (root, circle)
                    # dx=eps*A*d(lambda) modulo a term proportional to dA.
                    integral += wl*wt*wa*ws*A*circle_period(density, center, radius)
    return integral


def check_cd_single_meridian_layer():
    reference_circle = certify_patch_clearance()
    limit = patch_coefficient(0., 7, reference_circle)
    np.testing.assert_allclose(patch_coefficient(0., 5), limit, rtol=1e-9, atol=0)
    assert limit.real < -.3
    volume = prod(hi-lo for lo, hi in PATCH)
    print(f"Limiting patch period/eps (counterclockwise) = {volume*limit:.10e}")
    errors = []
    for eps in (1e-5, 1e-6, 1e-7):
        estimate = patch_coefficient(eps, 5, reference_circle)
        errors.append(abs(estimate-limit))
        print(f"eps={eps:g}: scaled patch error={errors[-1]:.3e}")
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < 2e-5
    print("PASS: compact single-meridian period has a nonzero linear regulator coefficient")


if __name__ == "__main__":
    check_cd_single_meridian_layer()
