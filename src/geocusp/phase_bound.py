"""The phase bound with a timelike diagonal: constants, conditions, and a point-wise survey helper.

See docs/GATE_6_PHASE.md.  With one diagonal timelike, on the contour of
docs/GATE_6_CONTOUR.md,

    |arg z|  <=  2 pi + e + max(K_in, K_bl, e, 3 e_sp),

provided the conditions of GATE_6_CONTOUR section 9 and the conditions (G1), (G2), (R4),
(B1)-(B5) below hold.  K_c = u_c - u_{c-1} - u_{c+1} with u_i = pi + arg x0i^2 on the
Feynman branch, at a cusp c whose diagonal is timelike; K_in, K_bl and max(e, 3 e_sp)
bound it for r <= R, R < r < 2R and r >= 2R.

`phase_bound(contour)` evaluates every constant and condition in floating point.
`k_values(contour, x0)` evaluates K_c and arg z at a real point x0, in float64.
"""
from __future__ import annotations

import numpy as np

from . import feynman_contour as fc

__all__ = ["phase_bound", "k_values", "timelike_cusps"]


def timelike_cusps(contour: fc.Contour) -> list[int]:
    """Indices j (0-based) of the cusps with D_j = 2 p_{j-1}.p_j > 0."""
    return [j for j in range(4) if 2 * fc.mdot(contour.p[(j - 1) % 4], contour.p[j]) > 0]


def phase_bound(contour: fc.Contour) -> dict:
    """Constants and conditions of GATE_6_PHASE for a square with one timelike diagonal.

    Returns a dict with the geometric constants, the constants of Propositions R4, R5, R6,
    one boolean per condition, `all_conditions` (their conjunction, including the
    certificate of GATE_6_CONTOUR), and `bound`, the proved bound on |arg z| minus 2 pi.
    Raises ValueError in the region (-,-), where 2 pi + 4e is already proved.
    """
    C = contour
    if not (C.s > 0 or C.t > 0):
        raise ValueError("phase_bound needs a timelike diagonal; in the region (-,-) "
                         "|arg z| <= 2 pi + 4e holds by Proposition 6")
    cert = fc.certificate(C)
    eu = np.linalg.norm
    Ds = [2 * fc.mdot(C.p[(j - 1) % 4], C.p[j]) for j in range(4)]
    Dabs = [abs(d) for d in Ds]
    Dmin, Dmax = min(Dabs), max(Dabs)
    D, Dp = max(Ds), min(Ds)                       # timelike diagonal, the other one
    u = -D - Dp
    Delta, P, A1, lam0 = cert["Delta"], cert["P"], cert["A1"], cert["lambda0"]
    B = max(fc.transverse_norm(C, j) for j in range(4))
    R, d0, kf, tau0, kap, nu = C.radius, C.delta0, C.kappa_far, C.tau0, C.kappa, C.nu
    e, theta, M = cert["excursion_bound"], cert["theta"], cert["C1_margin"]
    mu = max(64 * kap * P**2 / Dmin, nu * B**2)
    out = {"Delta": Delta, "P": P, "D": D, "D_other": Dp, "u": u, "B": B, "lambda0": lam0,
           "e": e, "theta": theta, "mu": mu, "M": M, "certificate": bool(cert["all"])}
    sqrt, atan, asin = np.sqrt, np.arctan, lambda x: np.arcsin(min(x, 1.0))

    # ---- outer zone (Proposition R4)
    q = (tau0 + Delta) / (tau0 - 2 * Delta) if tau0 > 2 * Delta else np.inf
    eps_k = kf**2 * (1 + 2 * sqrt(2) * q + q * q)
    e_sp = atan(2 * kf * (1 + sqrt(2) * q) / (1 - eps_k)) if eps_k < 1 else np.inf
    out.update(eps_kappa=eps_k, e_sp=e_sp, K_out=max(e, 3 * e_sp))
    out["R4_centres"] = bool(tau0 >= 3 * Delta and 2 * R > sqrt(2) * tau0)
    out["R4_trivial"] = bool(eps_k < 1 and e < np.pi / 2)

    # ---- inner zone (Proposition R5)
    out["G1"] = bool(4 * d0 * abs(u) < D * abs(Dp))
    out["G2"] = bool(4 * d0 < D)
    gamma1 = 1 - 2 * sqrt(d0 * abs(u) / (D * abs(Dp)))
    iota = 2 * lam0 * (R + Delta) ** 2
    v_R = 4 * lam0 * (R + Delta)
    out["R5_trivial"] = bool(v_R**2 < d0 and theta / 2 <= 1 and nu < 1 and 4 * kap**2 < 1)
    e4 = atan(4 * iota / (d0 - v_R**2)) if v_R**2 < d0 else np.inf
    beta_max = A1 * (R + Delta)
    alpha_max = beta_max * D / abs(Dp)
    beta_star = abs(Dp) / abs(u) * gamma1
    c1 = (2 * nu + 2 * kap * alpha_max) * 2 * D / d0 * max(1.0, beta_max)
    if gamma1 > 0:
        c1p = (2 * abs(Dp) * (kap * alpha_max**2 + (1 + beta_max) * (nu * alpha_max + kap * beta_max))
               / (0.5 * d0 * min(1.0, beta_star)))
    else:
        c1p = np.inf
    c_star = max(c1, c1p)
    eps_V = max(asin(theta / 2), asin(v_R**2 / d0))
    K_in_12 = 2 * e + e4
    K_in_3ab = e + 2 * eps_V
    K_in_3d = atan(2 * nu + 2 * c_star + 2 * iota / d0) + 2 * atan(iota / d0) + 3 * eps_V
    K_in = max(K_in_12, K_in_3ab, K_in_3d)
    out.update(gamma1=gamma1, iota=iota, v_R=v_R, e4=e4, beta_max=beta_max, alpha_max=alpha_max,
               beta_star=beta_star, c1=c1, c1_prime=c1p, c_star=c_star, eps_V=eps_V,
               K_in_cases_1_2=K_in_12, K_in_cases_3ab=K_in_3ab, K_in_case_3d=K_in_3d, K_in=K_in,
               K_near_cusp=2 * nu / (1 - nu**2) + 4 * kap**2 / (1 - 4 * kap**2))

    # ---- blend zone (Proposition R6)
    reach = 2 * R + Delta
    K = kf * (2 * R + tau0 + Delta)
    v_in = 4 * lam0 * reach
    v2 = (v_in + K) ** 2
    iota2 = 2 * lam0 * reach**2
    pi0 = 2 * lam0 * reach * P
    mu_m = 3 * mu + K * sqrt(3 * mu / (2 * kf * M)) if M > 0 else np.inf
    mu_p = mu_m + K**2 / (2 * kf * M) if M > 0 else np.inf
    omega0 = v_in * (2 * K + v_in)
    dv = d0 - v2
    eps = sqrt(2) * P / (R - Delta)
    b0 = eps * d0 / Dmin
    rad = R - Delta - B * sqrt(d0 + 2 * b0 * Dmax)
    Lam = rad / P
    a_min = Lam - b0
    d1 = d0 * (1 - eps)
    out.update(K_far=K, v_in=v_in, v=v_in + K, iota2=iota2, pi0=pi0, mu_minus=mu_m, mu_plus=mu_p,
               omega0=omega0, delta_v=dv, eps=eps, b0=b0, Lambda_b=Lam, a_min=a_min, delta1=d1)
    out["B1"] = bool(R >= tau0 and R > sqrt(5) * Delta)
    out["B2"] = bool(kap <= nu)
    out["B3"] = bool(b0 <= 0.5 and Lam >= 2 and 4 * d0 <= Dmin)
    out["B4"] = bool(d1 > v2)
    ok = out["B3"] and out["B4"] and dv > 0 and a_min > 0
    if ok:
        vs = (2 * nu + kap * b0**2 / (a_min * (1 + a_min))) + (2 * nu + kap * b0**2) / a_min
        F = 2 * K * P * mu_m / dv
        s0 = 8 * iota2 / dv + 2 * K * P * omega0 / dv**2
        s1 = max((iota2 + 2 * mu_p * pi0) / dv + F,
                 (1 + mu_p) * d0 / dv * vs + 2 * mu_m * pi0 / dv + iota2 * (1 / dv + 1 / (d1 - v2)) + F)
        s2 = (iota2 + 2 * mu_p * pi0) / dv + F
        s3 = mu_m + (1 + mu_p) ** 2 * vs
    else:
        vs = F = s0 = s1 = s2 = s3 = np.inf
    s = max(s0, s1, s2, s3)
    gap = (tau0 - 2 * Delta) ** 2
    out["R6_trivial"] = bool(v2 < d0 / 2 and v2 < gap)
    e_sp2 = (atan((4 * iota2 / gap + 2 * kf * (1 + sqrt(2) * q)) / (1 - v2 / gap)) if v2 < gap else np.inf)
    out["B5"] = bool(s < 1 and e + asin(s) + atan(mu_m) < np.pi / 2)
    K_bl = max(e + asin(s), 3 * e_sp2) if s < 1 else np.inf
    out.update(varsigma=vs, F=F, s0=s0, s1=s1, s2=s2, s3=s3, s=s, e_sp_blend=e_sp2, K_bl=K_bl)

    conditions = ["certificate", "R4_centres", "R4_trivial", "G1", "G2", "R5_trivial",
                  "B1", "B2", "B3", "B4", "B5", "R6_trivial"]
    out["all_conditions"] = all(out[k] for k in conditions)
    out["bound"] = float(e + max(K_in, K_bl, out["K_out"]))
    return {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in out.items()}


def k_values(contour: fc.Contour, x0) -> dict:
    """K_c at the cusps with a timelike diagonal, and arg z, at the real point x0 (float64).

    Returns {"K": {c: K_c}, "arg_z": Im log z, "u": the four u_i, "zone": 0, 1 or 2,
    "min_modulus": min_i |x0i^2| on the contour}.  Float64 loses the real parts within
    about 1e-15 of a tangency line of two adjacent cones; callers should skip points
    where `min_modulus` is not well above the rounding error of the invariants.
    """
    C = contour
    x0 = np.asarray(x0, dtype=float)
    X, _ = C.displaced_invariants(x0)
    u = np.array([C.branch_argument(val) for val in X]) + np.pi
    ks = {c: float(u[c] - u[(c - 1) % 4] - u[(c + 1) % 4]) for c in timelike_cusps(C)}
    r = np.linalg.norm(x0 - C.centre)
    zone = 0 if r <= C.radius else (1 if r < 2 * C.radius else 2)
    return {"K": ks, "arg_z": float(fc.log_z_imaginary_part(C, X)), "u": u, "zone": zone,
            "min_modulus": float(np.min(np.abs(X)))}
