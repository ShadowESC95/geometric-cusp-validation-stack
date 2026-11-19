"""The conditions and the constant of GATE_6_PHASE in interval arithmetic.

`phase_bound.phase_bound` evaluates them in floating point at one square.
Here they are evaluated with mpmath intervals for every square whose angle
lies in [angle_lo, angle_hi], with the contour constants fixed at those of
the midpoint, exactly as `certificate_iv.interval_certificate` does for the
conditions of GATE_6_CONTOUR section 9, and with its helpers.

What is enclosed, and how:

* Delta, P, D, D', u, A1, lambda0, the margin M of (C1): interval evaluation of
  their closed formulas on the interval square (as in certificate_iv).
* e: the upper bound used by certificate_iv (arctan x <= x).
* The transverse norms: `certificate_iv._transverse_bound_planar` (the square
  root of the trace of the Gram matrix), per frame; B_c is the largest over
  the frames of the cusps with the timelike diagonal, B_n over the other two,
  B over all four.  theta and mu are formed from B.  These are upper bounds,
  and every constant below is monotone in them in the safe direction.
* b_0, eps, delta_1, Lambda_c, Lambda_n, a_min, mu_-, mu_+, omega_0, s_0..s_3 and
  the rest: interval evaluation of their formulas.
* arctan x is replaced by x, and arcsin x by t - t^3/3 + t^5/5 with
  t = x / sqrt(1 - x^2) (by t itself if t > 1); both are upper bounds, so
  `bound_upper` is an upper bound of the float `bound` (slightly larger).
* `bound_upper` is the zone-by-zone bound of the theorem: the largest of b_in
  (Proposition R8), b_bl (Proposition R6) and, for r >= 2R, 0 if (R7)
  holds on the box and e + max(e, 3 e_sp) otherwise.

Every comparison is between the upper end of one interval and the lower end
of the other.  A True says: for every square with angle in the range, the
conditions hold for the contour with the constants of the midpoint square.
"""
from __future__ import annotations

import numpy as np
from mpmath import iv

from . import certificate_iv as ci
from . import feynman_contour as fc

_I, _max, _min, _abs, lt, le = ci._I, ci._max, ci._min, ci._abs, ci.lt, ci.le

CONDITIONS = ["certificate", "R4_centres", "R4_trivial", "G1", "G2", "R5_trivial",
              "B1", "B2", "B3", "B4", "B5", "B6", "R6_trivial"]

from .phase_bound import PAIRING


def _asin_upper(x):
    """An interval whose upper end bounds arcsin(x) for 0 <= x < 1.

    arcsin x = arctan t with t = x / sqrt(1 - x^2), and for 0 <= t <= 1 the alternating
    series gives arctan t <= t - t^3/3 + t^5/5 (an increasing function of t on [0, 1]);
    for t > 1 the bound arctan t <= t is used."""
    if not (x.a >= 0 and x.b < 1):
        raise ArithmeticError("arcsin argument not certainly in [0, 1)")
    t = iv.mpf((x / iv.sqrt(1 - x * x)).b)
    return t - t**3 / 3 + t**5 / 5 if t.b <= 1 else t


def interval_phase_bound(region: str, angle_lo: float, angle_hi: float, with_certificate: bool = True) -> dict:
    """Conditions of GATE_6_PHASE and an upper bound for |arg z| - 2 pi, for all squares
    with angle in [angle_lo, angle_hi] and the constants of the midpoint square."""
    if region == "--":
        raise ValueError("no timelike diagonal in the region (-,-)")
    mid = 0.5 * (angle_lo + angle_hi)
    C = fc.tuned_contour(region, mid)
    angle = iv.mpf([float(angle_lo), float(angle_hi)])
    p = ci._square(region, angle)
    x = [[iv.mpf(0)] * 4]
    for k in range(3):
        x.append([x[-1][i] + p[k][i] for i in range(4)])
    centre = [sum((xi[i] for xi in x), iv.mpf(0)) / 4 for i in range(4)]
    Delta = _max([ci._enorm([xi[i] - centre[i] for i in range(4)]) for xi in x])
    Pn = _max([ci._enorm(pj) for pj in p])
    Dsigned = [2 * ci._mdot(p[(j - 1) % 4], p[j]) for j in range(4)]
    tl = [j for j in range(4) if Dsigned[j].a > 0]
    sl = [j for j in range(4) if Dsigned[j].b < 0]
    if len(tl) != 2 or len(sl) != 2:
        raise ArithmeticError("signs of the diagonals not certain")
    D = iv.mpf([min(Dsigned[j].a for j in tl), max(Dsigned[j].b for j in tl)])           # > 0
    Dp = iv.mpf([min(-Dsigned[j].b for j in sl), max(-Dsigned[j].a for j in sl)])        # |D'| > 0
    Dmin, Dmax = _min([D, Dp]), _max([D, Dp])
    u_abs = D - Dp                                                                       # |u| = D - |D'|
    kap, nu, kf = _I(C.kappa), _I(C.nu), _I(C.kappa_far)
    d0, R, tau0 = _I(C.delta0), _I(C.radius), _I(C.tau0)
    r2, r5, r17 = iv.sqrt(2), iv.sqrt(5), iv.sqrt(17)
    out = {"region": region, "angle": [angle_lo, angle_hi]}
    out["certificate"] = bool(ci.interval_certificate(region, angle_lo, angle_hi)["all"]) if with_certificate else True

    # ---- quantities of section 9, as in certificate_iv
    A1 = 2 * Pn / Dmin
    A2 = 1 + 2 * A1 * Pn
    lam0 = 2 * fc.HESSIAN_BOUND * kap * A1 * Pn + nu * A2
    M = (tau0 - 2 * Delta) * (R / r17 - Delta) - (d0 + Delta * iv.sqrt(d0))
    reach = 2 * R + Delta
    K = kf * (2 * R + tau0 + Delta)
    v_in = 4 * lam0 * reach
    v2 = (v_in + K) ** 2
    xi = R * 2 / r5
    if not (lt(Delta, R / r17) and lt(4 * Delta, xi) and lt(v2, d0) and lt(0, M) and lt(2 * Delta, tau0)
            and lt(Delta, R) and lt(0, u_abs)):
        raise ArithmeticError("a denominator of section 9 is not certainly positive")
    case_b = 1 + 2 * Delta / (3 * (R / r17 - Delta))
    case_c = 1 + (r5 / 2 + Delta / xi) * (tau0 + Delta) / (0.75 * (xi - 4 * Delta))
    e = 2 * lam0 * reach**2 / (d0 - v2) + 4 * kf * _max([iv.mpf(1), case_b, case_c])    # upper bound of e
    Bj = [ci._transverse_bound_planar(p[(j - 1) % 4], p[j]) for j in range(4)]
    B, B_c, B_n = _max(Bj), _max([Bj[j] for j in tl]), _max([Bj[j] for j in sl])
    theta = _max([512 * kap * Pn**2 / Dmin, 8 * nu * B**2])
    mu = _max([64 * kap * Pn**2 / Dmin, nu * B**2])

    # ---- outer zone (R4)
    q = (tau0 + Delta) / (tau0 - 2 * Delta)
    eps_k = kf**2 * (1 + 2 * r2 * q + q * q)
    out["R4_centres"] = le(3 * Delta, tau0) and lt(r2 * tau0, 2 * R)
    out["R4_trivial"] = lt(eps_k, 1) and le(e, iv.pi / 4)
    if not out["R4_trivial"]:
        raise ArithmeticError("eps_kappa or e too large")
    e_sp = 2 * kf * (1 + r2 * q) / (1 - eps_k)                                           # tan e_sp >= e_sp
    K_out = _max([e, 3 * e_sp])
    out["R4_trivial"] = out["R4_trivial"] and lt(4 * e_sp, iv.pi)                      # tan e_sp >= e_sp
    # Proposition R7 (|arg z| <= 2 pi exactly for r >= 2R); reported, not used in bound_upper
    out["R7_outer_exact"] = out["R4_centres"] and le(2 * Delta * tau0, (R - tau0 / 2) * (tau0 - 2 * Delta))

    # ---- inner zone (R5)
    out["G1"] = lt(4 * d0 * u_abs, D * Dp)
    out["G2"] = lt(4 * d0, D)
    iota = 2 * lam0 * (R + Delta) ** 2
    v_R2 = (4 * lam0 * (R + Delta)) ** 2
    out["R5_trivial"] = lt(v_R2, d0) and le(theta, 2) and lt(nu, 1) and lt(4 * kap**2, 1)
    if not out["R5_trivial"]:
        raise ArithmeticError("a trivial condition of R5 fails")
    e4 = 4 * iota / (d0 - v_R2)
    beta_c = 2 * Pn * (R + Delta) / D
    G = 2 * D / d0 * _max([iv.mpf(1), beta_c])
    w = u_abs / Dp
    c_3d = G * (2 * nu * (1 + beta_c * w) + 2 * kap * (beta_c * D / Dp + u_abs / D + u_abs**2 * beta_c / (Dp * D)))
    eps_V = _max([_asin_upper(theta / 2), _asin_upper(v_R2 / d0)])
    K_in = _max([2 * e + e4, e + 2 * eps_V, (2 * nu + 2 * c_3d + 2 * iota / d0) + 2 * iota / d0 + 3 * eps_V])
    e_fr = iota / (d0 - v_R2)                                                            # arctan x <= x
    K_3d = (2 * nu + 2 * c_3d + 2 * iota / d0) + 2 * iota / d0 + 3 * eps_V
    b_in = e_fr + _max([2 * e_fr + e4, e_fr + 2 * eps_V, K_3d])                         # Proposition R8

    # ---- blend zone (R6)
    iota2 = 2 * lam0 * reach**2
    pi0 = 2 * lam0 * reach * Pn
    mu_m = 3 * mu + K * iv.sqrt(3 * mu / (2 * kf * M))
    mu_p = mu_m + K**2 / (2 * kf * M)
    omega0 = v_in * (2 * K + v_in)
    dv = d0 - v2
    eps = r2 * Pn / (R - Delta)
    b0 = eps * d0 / Dmin
    Lam_c = (R - Delta - B_c * iv.sqrt(d0 + b0**2 * D)) / Pn
    Lam_n = (R - Delta - B_n * iv.sqrt(d0 * (1 + eps) + 2 * b0**2 * Dp)) / Pn
    a_min = Lam_c - b0
    d1 = d0 * (1 - eps)
    out["B1"] = le(tau0, R) and lt(r5 * Delta, R)
    out["B2"] = le(kap, nu)
    out["B3"] = le(b0, 0.5) and le(2, Lam_c) and le(2, Lam_n) and le(4 * d0, Dmin)
    out["B4"] = lt(v2, d1)
    gap = (tau0 - 2 * Delta) ** 2
    out["R6_trivial"] = lt(v2, d0 / 2) and lt(v2, gap)
    # (B6)
    tau1 = 2 * Delta * tau0 / (tau0 - 2 * Delta)
    W1 = Delta * (tau0 + Delta) + Delta * (tau1 + tau0 + Delta)
    out["B6"] = False
    rho1 = None
    if lt(tau1, R):
        root = iv.sqrt(R * R - tau1 * tau1) - Delta
        if lt(tau1 + Delta, root):
            rho1 = root * root - (tau1 + Delta) ** 2
            out["B6"] = lt(0, rho1)
    out["B5"] = False
    if out["B3"] and out["B4"] and out["B6"] and out["R6_trivial"] and a_min.a > 0:
        unit = iota2 / dv
        half = _I(PAIRING) * unit
        a_c = 2 * Pn * reach / D
        a_n = 2 * Pn * reach / Dp
        vs_c = 2 * nu + 2 * D / d0 * (kap * a_c**2 + nu * a_c**2 + kap * a_c * _max([a_c + 1, iv.mpf(2)]))
        vs_n = 2 * Dp / d0 * (kap * a_n**2 + nu * (a_n**2 + a_n) + kap * a_n * (a_n + 1))
        far0 = 2 * K * Pn * omega0 / dv**2
        vterm = 2 * v2 * pi0 / dv**2
        s0 = _max([3 * half + vterm + far0, 2 * half + far0,
                   unit + 2 * half + (vs_c + vs_n) * (d0 / dv) ** 2 + vterm + far0])
        vs = (2 * nu + kap * b0**2 / (a_min * (1 + a_min))) + (2 * nu + kap * b0**2) / a_min
        F = 2 * K * Pn * mu_m / dv
        chi_eps = eps * eps * (3 - 2 * eps)
        s1 = _max([half + 2 * mu_p * pi0 / dv + F,
                   (1 + mu_p) * d0 / dv * vs + 2 * mu_m * pi0 / dv + chi_eps * iota2 * (1 / dv + 1 / (d1 - v2)) + F])
        s2 = half + 2 * mu_p * pi0 / dv + F
        s3 = mu_m + (1 + mu_p) ** 2 * vs
        s = _max([s0, s1, s2, s3])
        if s.b < 1:
            asin_s = _asin_upper(s)
            e_sp2 = (4 * iota2 / gap + 2 * kf * (1 + r2 * q)) / (1 - v2 / gap)
            out["B5"] = lt(e + asin_s + mu_m, iv.pi / 2) and lt(4 * e_sp2, iv.pi)
            e_bl_hi = half                                                               # arctan x <= x
            e_bl_lo = half + 2 * kf * (1 + W1 / rho1) / (1 - v2 / d0)
            b_bl = e_bl_hi + e_bl_lo + asin_s
            b_out = iv.mpf(0) if out["R7_outer_exact"] else e + K_out
            bound = _max([b_in, b_bl, b_out])
            out.update(s_upper=float(s.b), b_bl_upper=float(b_bl.b), bound_upper=float(bound.b))
    out.update(e_upper=float(e.b), K_in_upper=float(K_in.b), b_in_upper=float(b_in.b), K_out_upper=float(K_out.b),
               c_3d_upper=float(c_3d.b),
               Lambda_lower=float(min(Lam_c.a, Lam_n.a)), B_n_upper=float(B_n.b),
               ratio_s_over_t=[float((-D / Dp).a), float((-D / Dp).b)] if region == "+-"
               else [float((-Dp / D).a), float((-Dp / D).b)])
    out["all_conditions"] = all(out[k] for k in CONDITIONS)
    if "bound_upper" not in out:
        out["bound_upper"] = float("inf")
    return out


def cover(region: str, angle_lo: float, angle_hi: float, boxes: int, with_certificate: bool = True) -> dict:
    """`interval_phase_bound` on `boxes` adjacent angle intervals (the subdivision of
    `certificate_iv.cover`).  Reports where all conditions hold and the largest bound there."""
    edges = np.linspace(angle_lo, angle_hi, boxes + 1)
    failed, ratios, worst, where, r7_failed = [], [], 0.0, None, 0
    parts = {"e_upper": 0.0, "K_in_upper": 0.0, "b_in_upper": 0.0, "b_bl_upper": 0.0, "K_out_upper": 0.0,
             "c_3d_upper": 0.0, "s_upper": 0.0}
    lam = np.inf
    for lo, hi in zip(edges[:-1], edges[1:]):
        try:
            r = interval_phase_bound(region, float(lo), float(hi), with_certificate)
        except (ArithmeticError, ZeroDivisionError) as error:
            failed.append((float(lo), float(hi), [str(error)]))
            continue
        if not r["all_conditions"]:
            failed.append((float(lo), float(hi), [k for k in CONDITIONS if not r[k]]))
            continue
        ratios += r["ratio_s_over_t"]
        r7_failed += not r["R7_outer_exact"]
        lam = min(lam, r["Lambda_lower"])
        for k in parts:
            parts[k] = max(parts[k], r.get(k, 0.0))
        if r["bound_upper"] > worst:
            worst, where = r["bound_upper"], (float(lo), float(hi))
    out = {"region": region, "angle": [angle_lo, angle_hi], "boxes": boxes, "all": not failed,
           "failed": failed[:5], "failed_count": len(failed),
           "ratio_s_over_t": [min(ratios), max(ratios)] if ratios else None,
           "bound_upper": worst, "bound_upper_at": where, "Lambda_lower": float(lam),
           "R7_outer_exact": not failed and r7_failed == 0, "R7_failed_count": int(r7_failed)}
    out.update(parts)
    return out
