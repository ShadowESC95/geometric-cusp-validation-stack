"""The phase bound with a timelike diagonal: constants, conditions, and a point-wise survey helper.

See docs/GATE_6_PHASE.md.  With one diagonal timelike, on the contour of
docs/GATE_6_CONTOUR.md,

    |arg z| <= 2 pi + b_in         for r <= R        (Proposition R8),
    |arg z| <= 2 pi + b_bl         for R < r < 2R    (Proposition R6)        ,
    |arg z| <= 2 pi                for r >= 2R       (Proposition R7),

provided the conditions of GATE_6_CONTOUR section 9 and the conditions (R4), (R7), (G1),
(G2), (B1)-(B5) and the trivial ones hold.  K_c = u_c - u_{c-1} - u_{c+1} with
u_i = pi + arg x0i^2 on the Feynman branch, at a cusp c whose diagonal is timelike.
`bound` is the largest of the three zone bounds (`bound_inner`, `bound_blend`,
`bound_outer`); without (R7) the outer zone has e + max(e, 3 e_sp) instead of 0.

`phase_bound(contour)` evaluates every constant and condition in floating point.
`k_values(contour, x0)` evaluates K_c and arg z at a real point x0, in float64.
"""
from __future__ import annotations

import numpy as np

from . import feynman_contour as fc

__all__ = ["phase_bound", "k_values", "timelike_cusps", "PAIRING"]

# chi(t) <= PAIRING * (1 + t) for the cubic smoothstep chi: (F5) of GATE_6_PHASE section 5
PAIRING = 0.52


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
    Bj = [fc.transverse_norm(C, j) for j in range(4)]
    B = max(Bj)
    B_c = max(Bj[j] for j in range(4) if Ds[j] > 0)      # frames of the cusps with the timelike diagonal
    B_n = max(Bj[j] for j in range(4) if Ds[j] < 0)      # frames of their neighbours
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
    out["R4_trivial"] = bool(eps_k < 1 and e <= np.pi / 4 and 4 * e_sp < np.pi)
    # Proposition R7: |arg z| <= 2 pi exactly for r >= 2R
    out["R7_outer_exact"] = bool(out["R4_centres"] and (R - tau0 / 2) * (tau0 - 2 * Delta) >= 2 * Delta * tau0)

    # ---- inner zone (Proposition R5)
    out["G1"] = bool(4 * d0 * abs(u) < D * abs(Dp))
    out["G2"] = bool(4 * d0 < D)
    gamma1 = 1 - 2 * sqrt(d0 * abs(u) / (D * abs(Dp)))
    iota = 2 * lam0 * (R + Delta) ** 2
    v_R = 4 * lam0 * (R + Delta)
    out["R5_trivial"] = bool(v_R**2 < d0 and theta / 2 <= 1 and nu < 1 and 4 * kap**2 < 1)
    e4 = atan(4 * iota / (d0 - v_R**2)) if v_R**2 < d0 else np.inf
    beta_c = 2 * P * (R + Delta) / D               # |alpha|, |beta| <= beta_c in the frame of c, r <= R
    G = 2 * D / d0 * max(1.0, beta_c)              # |A|/|r_c|, |B|/|r_c| <= G where the neighbour frame is on
    w = abs(u) / abs(Dp)
    c_3d = G * (2 * nu * (1 + beta_c * w)
                + 2 * kap * (beta_c * D / abs(Dp) + abs(u) / D + u * u * beta_c / (abs(Dp) * D)))
    eps_V = max(asin(theta / 2), asin(v_R**2 / d0))
    K_in_12 = 2 * e + e4
    K_in_3ab = e + 2 * eps_V
    # for r <= R only the frame i+2 can push arg x0i^2 out of [-pi, 0]: excursion e_fr <= e
    e_fr = atan(iota / (d0 - v_R**2)) if v_R**2 < d0 else np.inf
    K_in_3d = atan(2 * nu + 2 * c_3d + 2 * iota / d0) + 2 * atan(iota / d0) + 3 * eps_V
    K_in = max(K_in_12, K_in_3ab, K_in_3d)
    out.update(gamma1=gamma1, iota=iota, v_R=v_R, e4=e4,
               beta_c=beta_c, G=G, c_3d=c_3d, eps_V=eps_V,
               K_in_cases_1_2=K_in_12, K_in_cases_3ab=K_in_3ab, K_in_case_3d=K_in_3d, K_in=K_in,
               K_near_cusp=2 * nu / (1 - nu**2) + 4 * kap**2 / (1 - 4 * kap**2), e_fr=e_fr,
               K_in_sharp=max(2 * e_fr + e4, e_fr + 2 * eps_V, K_in_3d))

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
    Lam_c = (R - Delta - B_c * sqrt(d0 + b0**2 * D)) / P
    Lam_n = (R - Delta - B_n * sqrt(d0 * (1 + eps) + 2 * b0**2 * abs(Dp))) / P
    a_min = Lam_c - b0
    d1 = d0 * (1 - eps)
    out.update(K_far=K, v_in=v_in, v=v_in + K, iota2=iota2, pi0=pi0, mu_minus=mu_m, mu_plus=mu_p,
               omega0=omega0, delta_v=dv, eps=eps, b0=b0, Lambda_c=Lam_c, Lambda_n=Lam_n, a_min=a_min, delta1=d1,
               B_c=B_c, B_n=B_n)
    out["B1"] = bool(R >= tau0 and R > sqrt(5) * Delta)
    out["B2"] = bool(kap <= nu)
    out["B3"] = bool(b0 <= 0.5 and Lam_c >= 2 and Lam_n >= 2 and 4 * d0 <= Dmin)
    out["B4"] = bool(d1 > v2)
    gap = (tau0 - 2 * Delta) ** 2
    out["R6_trivial"] = bool(v2 < d0 / 2 and v2 < gap)
    # (B6): where n_i.w_i < 0 in the blend zone the invariants are far on the spacelike side
    tau1 = 2 * Delta * tau0 / (tau0 - 2 * Delta) if tau0 > 2 * Delta else np.inf
    W1 = Delta * (tau0 + Delta) + Delta * (tau1 + tau0 + Delta)
    root = sqrt(R * R - tau1 * tau1) - Delta if tau1 < R else -np.inf
    rho1 = root * root - (tau1 + Delta) ** 2 if root > tau1 + Delta else -np.inf
    out["B6"] = bool(rho1 > 0)
    ok = out["B3"] and out["B4"] and out["B6"] and out["R6_trivial"] and dv > 0 and a_min > 0
    if ok:
        unit = iota2 / dv                              # |J_i^(j)| / |Re X_i| for an invariant beyond delta_0
        half = PAIRING * unit                          # the same for the frame i+2 ((F5) of GATE_6_PHASE section 5)
        a_c = 2 * P * reach / D                        # coordinates in the frame of c, r < 2R
        a_n = 2 * P * reach / abs(Dp)                  # coordinates in the frames of its neighbours
        vs_c = 2 * nu + 2 * D / d0 * (kap * a_c**2 + nu * a_c**2 + kap * a_c * max(a_c + 1, 2.0))
        vs_n = 2 * abs(Dp) / d0 * (kap * a_n**2 + nu * (a_n**2 + a_n) + kap * a_n * (a_n + 1))
        far0 = 2 * K * P * omega0 / dv**2
        vterm = 2 * v2 * pi0 / dv**2
        s0_pp = 3 * half + vterm + far0
        s0_mix = 2 * half + far0
        s0_nn = unit + 2 * half + (vs_c + vs_n) * (d0 / dv) ** 2 + vterm + far0
        s0 = max(s0_pp, s0_mix, s0_nn)
        vs = (2 * nu + kap * b0**2 / (a_min * (1 + a_min))) + (2 * nu + kap * b0**2) / a_min
        F = 2 * K * P * mu_m / dv
        chi_eps = eps * eps * (3 - 2 * eps)            # weight of the frame c+1 in case N1c
        s1 = max(half + 2 * mu_p * pi0 / dv + F,
                 (1 + mu_p) * d0 / dv * vs + 2 * mu_m * pi0 / dv + chi_eps * iota2 * (1 / dv + 1 / (d1 - v2)) + F)
        s2 = half + 2 * mu_p * pi0 / dv + F
        s3 = mu_m + (1 + mu_p) ** 2 * vs
        e_bl_hi = atan(half)                                                  # u_i <= pi + e_bl_hi
        e_bl_lo = atan(half + 2 * kf * (1 + W1 / rho1) / (1 - v2 / d0))       # u_i >= -e_bl_lo
    else:
        unit = half = a_c = a_n = vs_c = vs_n = s0_pp = s0_mix = s0_nn = np.inf
        vs = F = s0 = s1 = s2 = s3 = e_bl_hi = e_bl_lo = np.inf
    s = max(s0, s1, s2, s3)
    e_sp2 = (atan((4 * iota2 / gap + 2 * kf * (1 + sqrt(2) * q)) / (1 - v2 / gap)) if v2 < gap else np.inf)
    out["B5"] = bool(s < 1 and e + asin(s) + atan(mu_m) < np.pi / 2 and 4 * e_sp2 < np.pi)
    K_bl = max(e_bl_lo + asin(s), 3 * e_sp2) if s < 1 else np.inf
    b_bl = e_bl_hi + e_bl_lo + asin(s) if s < 1 else np.inf
    out.update(varsigma=vs, F=F, s0=s0, s0_pp=s0_pp, s0_mix=s0_mix, s0_nn=s0_nn, s1=s1, s2=s2, s3=s3, s=s,
               unit=unit, a_c=a_c, a_n=a_n, varsigma_c=vs_c, varsigma_n=vs_n, tau1=tau1, W1=W1, rho1=rho1,
               e_bl_hi=e_bl_hi, e_bl_lo=e_bl_lo, e_sp_blend=e_sp2, K_bl=K_bl, b_bl=b_bl)

    conditions = ["certificate", "R4_centres", "R4_trivial", "G1", "G2", "R5_trivial",
                  "B1", "B2", "B3", "B4", "B5", "B6", "R6_trivial"]
    out["all_conditions"] = all(out[k] for k in conditions)
    out["bound_inner"] = float(min(e, e_fr) + min(K_in, out["K_in_sharp"]))
    out["bound_blend"] = float(b_bl)
    out["bound_outer"] = 0.0 if out["R7_outer_exact"] else float(e + out["K_out"])
    out["bound"] = max(out["bound_inner"], out["bound_blend"], out["bound_outer"])
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
