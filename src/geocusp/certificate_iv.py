"""The certificate of GATE_6_CONTOUR section 9 in interval arithmetic.

`certificate` in feynman_contour.py evaluates the conditions in floating
point at one square. Here the same conditions are evaluated with mpmath
intervals, for every square whose shape parameter lies in an interval
[angle_lo, angle_hi], with the contour constants fixed at those of the
midpoint. Every comparison is made between the upper end of one interval
and the lower end of the other, so a True is a statement about all squares
in the range, with rounding accounted for.

Two bounds are replaced by weaker, simpler ones that imply the originals:
arctan(x) <= x, and B^2 <= trace of the Gram matrix of the transverse basis.
"""
from __future__ import annotations

import numpy as np
from mpmath import iv, mp

from . import feynman_contour as fc

iv.dps = 30


def _I(x):
    # a binary float is converted exactly, so the interval contains the number actually used
    return x if isinstance(x, iv.mpf) else iv.mpf(float(x) if isinstance(x, (float, np.floating)) else x)


def _max(xs):
    xs = [_I(x) for x in xs]
    return iv.mpf([max(x.a for x in xs), max(x.b for x in xs)])


def _min(xs):
    xs = [_I(x) for x in xs]
    return iv.mpf([min(x.a for x in xs), min(x.b for x in xs)])


def _abs(x):
    x = _I(x)
    if x.a >= 0:
        return x
    if x.b <= 0:
        return -x
    return iv.mpf([0, max(-x.a, x.b)])


def lt(a, b) -> bool:
    """Certainly a < b."""
    return bool(_I(a).b < _I(b).a)


def le(a, b) -> bool:
    return bool(_I(a).b <= _I(b).a)


def _mdot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def _enorm(a):
    return iv.sqrt(sum((x * x for x in a), iv.mpf(0)))


def _square(region, angle):
    c, sn = iv.cos(angle), iv.sin(angle)
    one, zero = iv.mpf(1), iv.mpf(0)
    a = [one, zero, zero, one]
    b = [one, zero, zero, -one]
    m = [-one, -sn, zero, -c]
    q = [-one, sn, zero, c]
    if region == "--":
        return [a, m, b, q]
    if region == "+-":
        return [a, b, m, q]
    if region == "-+":
        return [q, a, b, m]
    raise ValueError("region must be '--', '+-' or '-+'")


def _transverse_bound(pm, pp, seeds):
    """Upper bound for B: sqrt of the trace of the Gram matrix of a transverse
    basis with yperp^2 = -1, built from the given seed vectors."""
    D = 2 * _mdot(pm, pp)
    basis = []
    for k in seeds:
        seed = [iv.mpf(1 if i == k else 0) for i in range(4)]
        cp, cm = 2 * _mdot(seed, pp) / D, 2 * _mdot(seed, pm) / D
        v = [seed[i] - cp * pm[i] - cm * pp[i] for i in range(4)]
        for e in basis:
            coeff = _mdot(v, e)
            v = [v[i] + coeff * e[i] for i in range(4)]
        norm = -_mdot(v, v)
        if not norm.a > 0:
            raise ArithmeticError("transverse seed not certainly spacelike")
        root = iv.sqrt(norm)
        basis.append([x / root for x in v])
    return iv.sqrt(sum((x * x for e in basis for x in e), iv.mpf(0)))


def _float_seeds(C, j):
    """Which seed vectors the floating-point Gram-Schmidt accepts."""
    pm, pp = C.p[(j - 1) % 4], C.p[j]
    D = 2 * fc.mdot(pm, pp)
    basis, used = [], []
    for k, seed in enumerate(np.eye(4)):
        v = seed - 2 * fc.mdot(seed, pp) / D * pm - 2 * fc.mdot(seed, pm) / D * pp
        for e in basis:
            v = v + fc.mdot(v, e) * e
        norm = -fc.mdot(v, v)
        if norm > 1e-3:
            basis.append(v / np.sqrt(norm))
            used.append(k)
        if len(basis) == 2:
            break
    return used


def interval_certificate(region: str, angle_lo: float, angle_hi: float) -> dict:
    """Conditions of section 9 for all squares with angle in [angle_lo, angle_hi]."""
    mid = 0.5 * (angle_lo + angle_hi)
    C = fc.tuned_contour(region, mid)
    ref = fc.certificate(C)
    angle = iv.mpf([float(angle_lo), float(angle_hi)])
    p = _square(region, angle)
    x = [[iv.mpf(0)] * 4]
    for k in range(3):
        x.append([x[-1][i] + p[k][i] for i in range(4)])
    centre = [sum((xi[i] for xi in x), iv.mpf(0)) / 4 for i in range(4)]
    Delta = _max([_enorm([xi[i] - centre[i] for i in range(4)]) for xi in x])
    Pn = _max([_enorm(pj) for pj in p])
    Dsigned = [2 * _mdot(p[(j - 1) % 4], p[j]) for j in range(4)]
    Ds = [_abs(d) for d in Dsigned]
    Dmin, Dmax = _min(Ds), _max(Ds)
    s, t = Dsigned[1], Dsigned[2]                    # s = 2 p1.p2, t = 2 p2.p3
    u_abs = _abs(s + t)
    kap, nu, kf = _I(C.kappa), _I(C.nu), _I(C.kappa_far)
    d0, R, tau0 = _I(C.delta0), _I(C.radius), _I(C.tau0)
    r17, r5 = iv.sqrt(17), iv.sqrt(5)
    out = {"region": region, "angle": [angle_lo, angle_hi]}

    A1 = 2 * Pn / Dmin
    A2 = 1 + 2 * A1 * Pn
    lam0 = 2 * fc.HESSIAN_BOUND * kap * A1 * Pn + nu * A2

    # H1: a value eps, slightly below the floating-point gap, for which the contradiction is certain
    eps = _I(0.99 * ref["four_cone_gap"])
    holds = False
    for j in range(4):
        D, Dp = Ds[j], Ds[(j + 1) % 4]
        rho = eps + 4 * eps**2 / D
        shift = 4 * eps * Dp / D + eps + 2 * iv.sqrt(Dp * u_abs / D * rho)
        holds = holds or lt(shift, Dp - eps)
    out["H1_gap_certified"] = holds
    out["H1_frame_on"] = lt(2 * d0, eps)

    out["C0_centres"] = le(3 * Delta, tau0)
    margin = (tau0 - 2 * Delta) * (R / r17 - Delta) - (d0 + Delta * iv.sqrt(d0))
    out["C1_far_admissible"] = lt(Delta, R / r17) and lt(0, margin)
    xi0 = 2 * (2 * Delta + iv.sqrt(d0))
    out["C2_theta_is_one"] = le(r5 / 2 * xi0, R)
    xi_min = R * 4 / r17
    worst = (xi_min / 4 + Delta) ** 2 - (xi_min - Delta) ** 2
    out["C3_wedge_free"] = lt(Delta, xi_min) and lt(worst, -d0)

    reach = 2 * R + Delta
    frame_term = 2 * lam0 * reach**2
    v_inner = 4 * lam0 * reach + kf * (2 * R + tau0 + Delta)
    xi = R * 2 / r5
    out["C4_wedge_margin"] = lt(4 * Delta, xi)
    case_b = 1 + 2 * Delta / (3 * (R / r17 - Delta))
    case_c = 1 + (r5 / 2 + Delta / xi) * (tau0 + Delta) / (0.75 * (xi - 4 * Delta))
    # arctan(x) <= x
    e_bound = frame_term / (d0 - v_inner**2) + 4 * kf * _max([iv.mpf(1), case_b, case_c])
    out["D_excursion_small"] = lt(v_inner**2, d0 / 2) and lt(e_bound, 0.2)

    j_inner = 4 * (3 * reach**2 * lam0 / d0 + lam0)
    j_outer = kf * (1 + 7.5 * (1 + (tau0 + Delta) / R))
    j_blend = 1.5 / R * (4 * lam0 * reach + kf * (2 * R + tau0 + Delta))
    out["E_local_diffeomorphism"] = lt(j_inner + j_outer + j_blend, 1)

    B = _max([_transverse_bound(p[(j - 1) % 4], p[j], _float_seeds(C, j)) for j in range(4)])
    theta = _max([512 * kap * Pn**2 / Dmin, 8 * nu * B**2])
    out["F_theta_small"] = le(theta, 1)

    b0 = 4 * d0 / Dmin
    a1 = d0 / (8 * Dmax)
    rho1 = (d0 / 8) ** 2 / _max([iv.mpf(1), Dmax * u_abs / Dmin]) / 4
    out["S0_strips_overlap"] = lt(b0, 0.5)
    out["S2_near_invariants"] = le(a1 * Dmax * (1 + a1) + rho1, d0)
    far_lower = Dmin * (b0 - a1) - a1 * Dmax - rho1 - 2 * iv.sqrt(rho1 * Dmax * u_abs / Dmin)
    out["S1_frame_fully_on"] = le(2 * d0, far_lower)
    out["S3_strips_inside"] = le(Delta + a1 * Pn + B * iv.sqrt(rho1), R)

    omega = _I(0.98 * ref["omega"])      # any omega with 2 omega <= c_u(omega) will do
    c_u = (tau0 - 2 * Delta) * (1 / r17 - Delta / (2 * R)) - Delta * iv.sqrt(omega / (2 * R))
    out["I0_c_u_positive"] = lt(0, c_u) and le(2 * omega, c_u)
    r0 = 2 * R
    lhs = kf**2 * (2 * (r0 + Delta) * (tau0 + Delta) + (tau0 + Delta) ** 2)
    out["I1_real_part"] = le(lhs, omega * r0 / 4) and le(2 * kf**2 * (tau0 + Delta), omega / 4)
    out["I2_theta_is_one"] = (le(2 * (2 * Delta + iv.sqrt(omega * r0)), 2 / r5 * r0)
                              and lt(5 * omega / 4, r0))

    mu = _max([64 * kap * Pn**2 / Dmin, nu * B**2])
    K = kf * (2 * R + tau0 + Delta)
    if margin.a > 0:
        mu_far = kf * (1 + (1.5 * d0 + (tau0 + Delta) ** 2) / (2 * margin))
        mu_q = _max([3 * mu + iv.sqrt(3 * mu), mu_far + iv.sqrt(3 * mu) * K**2 / (2 * kf * margin)])
        gamma = fc.PHASE_CORNER_MARGIN + 3 * e_bound
        out["Q1_corner_excluded"] = lt(gamma, iv.pi / 2) and lt(iv.sin(gamma) * mu_q, iv.cos(gamma))
        out["mu_Q_upper"] = float(mu_q.b)
    else:
        out["Q1_corner_excluded"] = False
    out["Q3_corner_excluded_far"] = lt(2 * kf**2 * Delta**2 / (1 - kf**2),
                                       (2 * R / r17 - Delta) * (tau0 - 2 * Delta))
    n0_min = R / r17 - Delta
    out["Q4_outer_zone"] = (lt(kf * (tau0 + Delta), (1 - kf) * n0_min)
                            and lt(8 * kf**2 * Delta / (3 * (1 - kf**2)), 2 * n0_min))

    out["excursion_upper"] = float(e_bound.b)
    out["ratio_s_over_t"] = [float((s / t).a), float((s / t).b)]
    out["all"] = all(v for v in out.values() if isinstance(v, bool))
    return out


def cover(region: str, angle_lo: float, angle_hi: float, boxes: int) -> dict:
    """Run the interval certificate on `boxes` adjacent angle intervals."""
    edges = np.linspace(angle_lo, angle_hi, boxes + 1)
    failed, ratios, worst_e, worst_mu = [], [], 0.0, 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        try:
            r = interval_certificate(region, float(lo), float(hi))
        except (ArithmeticError, ZeroDivisionError) as error:
            failed.append((float(lo), float(hi), [str(error)]))
            continue
        ratios += r["ratio_s_over_t"]
        worst_e = max(worst_e, r["excursion_upper"])
        worst_mu = max(worst_mu, r.get("mu_Q_upper", np.inf))
        if not r["all"]:
            failed.append((float(lo), float(hi), [k for k, v in r.items() if v is False]))
    return {"region": region, "angle": [angle_lo, angle_hi], "boxes": boxes,
            "all": not failed, "failed": failed[:5], "failed_count": len(failed),
            "ratio_s_over_t": [min(ratios), max(ratios)],
            "excursion_upper": worst_e, "mu_Q_upper": worst_mu}


# --------------------------------------------------------------------------
# The inner field near the edges, in high precision
# --------------------------------------------------------------------------
def edge_sign_check(region: str = "--", angle: float = 1.1, dps: int = 100,
                    samples: int = 60, seed: int = 0) -> dict:
    """Sign of Im x0i^2 on the contour at points within 10^-6 ... 10^-30 of
    the four edges, in `dps`-digit arithmetic with exactly null edges.

    In float64 the edges are null only to 1e-16 and the sign is lost within
    about 1e-9 of an edge. Here the inner field is re-evaluated with mpmath.
    Returns the largest Im x0i^2 over invariants within delta0 of their
    cone (must be negative) and the number of such invariants.
    """
    contour = fc.tuned_contour(region, angle)
    rng = np.random.default_rng(seed)
    old = mp.dps
    mp.dps = dps
    try:
        c, sn = mp.cos(mp.mpf(angle)), mp.sin(mp.mpf(angle))
        one, zero = mp.mpf(1), mp.mpf(0)
        a, b = [one, zero, zero, one], [one, zero, zero, -one]
        m, q = [-one, -sn, zero, -c], [-one, sn, zero, c]
        p = {"--": [a, m, b, q], "+-": [a, b, m, q], "-+": [q, a, b, m]}[region]
        dot = lambda u, v: u[0] * v[0] - u[1] * v[1] - u[2] * v[2] - u[3] * v[3]
        add = lambda u, v, k=1: [u[i] + k * v[i] for i in range(4)]
        x = [[zero] * 4]
        for k in range(3):
            x.append(add(x[-1], p[k]))
        kap, nu, d0 = mp.mpf(contour.kappa), mp.mpf(contour.nu), mp.mpf(contour.delta0)
        dphi = lambda t: zero if t <= 0 else (t * t if t <= 1 else 2 * t - 1)
        step = lambda u: zero if u <= 0 else (one if u >= 1 else u * u * (3 - 2 * u))

        def field(x0):
            X = [dot(add(x0, xi, -1), add(x0, xi, -1)) for xi in x]
            V = [zero] * 4
            for j in range(4):
                weight = step(abs(X[(j + 2) % 4]) / d0 - 1)
                if weight == 0:
                    continue
                pm, pp = p[(j - 1) % 4], p[j]
                D = 2 * dot(pm, pp)
                y = add(x0, x[j], -1)
                al, be = 2 * dot(y, pp) / D, 2 * dot(y, pm) / D
                yperp = add(add(y, pm, -al), pp, -be)
                l1, l2, l3 = dphi(-be), dphi(al), dphi(be - al - 1)
                d_beta, d_alpha = -l1 + l3, l2 - l3
                Vj = add([-mp.sign(D) * kap * (d_beta * pm[i] + d_alpha * pp[i]) for i in range(4)],
                         yperp, nu)
                V = add(V, Vj, weight)
            return X, V

        worst, count, smallest = -mp.inf, 0, mp.inf
        for _ in range(samples):
            j = int(rng.integers(0, 4))
            along = mp.mpf(float(rng.uniform(0.05, 0.95)))
            size = mp.mpf(10) ** (-int(rng.integers(6, 31)))
            offset = [mp.mpf(float(v)) * size for v in rng.normal(size=4)]
            x0 = add(add(x[j], p[j], along), offset)
            X, V = field(x0)
            for i in range(4):
                if abs(X[i]) <= d0:
                    im = 2 * dot(V, add(x0, x[i], -1))
                    worst, count = max(worst, im), count + 1
                    smallest = min(smallest, abs(im))
        return {"region": region, "dps": dps, "pairs": count,
                "largest_im": float(worst), "smallest_abs_im": float(smallest),
                "all_negative": bool(worst < 0)}
    finally:
        mp.dps = old
