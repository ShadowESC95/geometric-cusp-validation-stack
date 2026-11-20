"""The conditions of certificate_iv and phase_bound_iv in scaled variables, uniformly
up to the degenerate squares (docs/GATE_6_PHASE.md section 7.1).

    PYTHONPATH=src python -m geocusp.scaled_limit --extend 0.1 --boxes 100

Setting.  With d the distance of the angle from the degenerate end (d = pi - angle at the
end angle -> pi in the regions (+,-), (-,+), (-,-); d = angle at the end angle -> 0 in
(-,-)) put sigma = sin(d/2).  Every quantity Q used by the two interval modules is written

    Q(sigma) = m(sigma) * sigma**k,      k rational (multiples of 1/4 occur),

where m is given by an explicit expression in sigma that is continuous on [0, sigma0]
(no division by sigma anywhere).  The class `S` carries an interval enclosure of m over a
box of sigma and the exponent k.  Arithmetic factors out the smaller exponent, e.g.

    m1 s^k1 + m2 s^k2 = (m1 + m2 s^(k2-k1)) s^k1     (k1 <= k2),

and a comparison Q1 < Q2 is decided at the common exponent k = min(k1, k2), i.e. after
dividing both sides by s^k > 0.  So a True is a statement about every sigma in the box with
sigma > 0, and the box may contain sigma = 0: the powers of d have been divided out.

Geometry.  Only two closed forms are not computed from the edge vectors in intervals,
because evaluating them that way loses the order in sigma by cancellation:
  * 2 x.y = -2 f for the four pairs of edges of the square that are not {a,b}, {m,q},
    with f = 1 + cos(angle) or 1 - cos(angle), i.e. f = 2 sigma^2 (the small diagonal
    |D'| = 4 sigma^2) or f = 2 cos^2(d/2);
  * the bound 1 + |v|^2/(x.y)^2 of certificate_iv._transverse_bound_planar for those
    pairs, which equals 4/f exactly (|v|^2 = 2 sin^2(angle) + f^2 = 4 f - f^2).
Both are algebraic identities in the edge vectors (checked with sympy and against the
interval code in test_scaled_limit.py).  Delta, P, the timelike D = 4 and the transverse
bound for the pairs {a,b}, {m,q} are computed from the edge vectors in intervals.

Family of constants (FAMILY), with |D'| = 4 sigma^2 the small diagonal:
    delta0 = c1 |D'|,  kappa = nu = c2 delta0 |D'|,  kappa' = c3 |D'|^(1/2),  R, tau0 fixed;
(H1) is verified at eps = e_h1 |D'| and (I) at the fixed omega.

The transcriptions `scaled_certificate` and `scaled_phase_bound` follow
certificate_iv.interval_certificate and phase_bound_iv.interval_phase_bound line by line,
with the same variable names (test_scaled_limit.py compares every named intermediate with
the repository code at several d > 0).
"""
from __future__ import annotations

import sys
from fractions import Fraction

from mpmath import iv, mp

from . import certificate_iv as ci
from . import feynman_contour as fc
from . import phase_bound_iv as pbi
from .phase_bound import PAIRING

iv.dps = 30

# --------------------------------------------------------------------------------------
# The scaled family.  All numbers are dyadic, so that kappa = nu is a point interval and
# (B2) kappa <= nu is decided exactly.
# --------------------------------------------------------------------------------------
FAMILY = {
    "c1": 0.046875,          # delta0 / |D'|                       (3/64)
    "c2": 2.0 ** -20,        # kappa / (delta0 |D'|) = nu / (...)   (9.54e-7)
    "c3": 2.0 ** -9,         # kappa' / |D'|^(1/2)  (~ kappa'/d)    (1.95e-3)
    "R": 9.0,
    "tau0": 9.0,
    "e_h1": 0.125,           # (H1) verified at eps = e_h1 |D'|
    "omega": 0.25,           # (I) any omega with 2 omega <= c_u(omega)
}


# --------------------------------------------------------------------------------------
# Scaled interval numbers
# --------------------------------------------------------------------------------------
def _rpow(x, q):
    """x**q for an interval x >= 0 and rational q >= 0 with a power-of-two denominator."""
    q = Fraction(q)
    if q < 0:
        raise ValueError("negative power")
    if q == 0:
        return iv.mpf(1)
    den = q.denominator
    while den > 1:
        if den % 2:
            raise ValueError("denominator must be a power of two")
        x = iv.sqrt(x)
        den //= 2
    return x ** q.numerator


class Box:
    """An interval [lo, hi] of sigma, 0 <= lo <= hi < 1."""

    def __init__(self, lo, hi=None):
        hi = lo if hi is None else hi
        self.sig = iv.mpf([lo, hi])
        if not (self.sig.a >= 0 and self.sig.b < 1):
            raise ValueError("need 0 <= sigma < 1")
        self.lo, self.hi = self.sig.a, self.sig.b
        self._cache = {}

    def pw(self, q):
        """An enclosure of sigma**q over the box (q >= 0)."""
        q = Fraction(q)
        if q not in self._cache:
            self._cache[q] = _rpow(self.sig, q)
        return self._cache[q]

    def pw_hi(self, q):
        return _rpow(iv.mpf(self.hi), q)

    def pw_lo(self, q):
        return _rpow(iv.mpf(self.lo), q)


def _dec(x):
    """An interval containing the binary float x and the decimal number it prints as
    (so that d_hi = 1e-6 covers the real number 10^-6, which is above the float)."""
    if isinstance(x, str):
        return iv.mpf(x)
    a, b = iv.mpf(float(x)), iv.mpf(repr(float(x)))
    return iv.mpf([min(a.a, b.a), max(a.b, b.b)])


def box_from_d(d_lo, d_hi):
    """The sigma-box containing sin(d/2) for every d in [d_lo, d_hi] (0 <= d_lo <= d_hi < pi)."""
    s = iv.sin(iv.mpf([_dec(d_lo).a, _dec(d_hi).b]) / 2)
    return Box(max(s.a, 0), s.b)


class S:
    """m * sigma**k for every sigma in `box`, with m an interval."""

    __slots__ = ("m", "k", "box")

    def __init__(self, m, k=0, box=None):
        self.m = ci._I(m)
        self.k = Fraction(k)
        self.box = box

    def __repr__(self):
        return f"S({self.m}, k={self.k})"

    def is_zero(self):
        return self.m.a == 0 and self.m.b == 0

    # arithmetic -----------------------------------------------------------------
    def __add__(self, other):
        mx, my, k, box = _pair(self, other)
        return S(mx + my, k, box)

    __radd__ = __add__

    def __sub__(self, other):
        mx, my, k, box = _pair(self, other)
        return S(mx - my, k, box)

    def __rsub__(self, other):
        return _co(other) - self

    def __neg__(self):
        return S(-self.m, self.k, self.box)

    def __mul__(self, other):
        o = _co(other)
        return S(self.m * o.m, self.k + o.k, self.box or o.box)

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = _co(other)
        if o.m.a <= 0 <= o.m.b:
            raise ZeroDivisionError("scaled denominator not certainly non-zero")
        return S(self.m / o.m, self.k - o.k, self.box or o.box)

    def __rtruediv__(self, other):
        return _co(other) / self

    def __pow__(self, n):
        if not isinstance(n, int):
            raise TypeError("integer powers only")
        return S(self.m ** n, self.k * n, self.box)


def _co(x):
    return x if isinstance(x, S) else S(x)


def _pair(x, y):
    """Both numbers at the common exponent min(kx, ky); an exact zero adopts the other's."""
    x, y = _co(x), _co(y)
    box = x.box or y.box
    if x.is_zero():
        return iv.mpf(0), y.m, y.k, box
    if y.is_zero():
        return x.m, iv.mpf(0), x.k, box
    k = min(x.k, y.k)
    mx = x.m if x.k == k else x.m * box.pw(x.k - k)
    my = y.m if y.k == k else y.m * box.pw(y.k - k)
    return mx, my, k, box


RECORD = None      # a list to log every labelled comparison, or None


def _log(label, op, x, y, mx, my, k, ok):
    if RECORD is not None and label:
        RECORD.append({"label": label, "op": op, "k_lhs": _co(x).k, "k_rhs": _co(y).k, "k": k,
                       "lhs": mx, "rhs": my, "ok": ok})


def lt(x, y, label=None) -> bool:
    """Certainly x < y for every sigma > 0 in the box."""
    mx, my, k, _ = _pair(x, y)
    ok = bool(mx.b < my.a)
    _log(label, "<", x, y, mx, my, k, ok)
    return ok


def le(x, y, label=None) -> bool:
    mx, my, k, _ = _pair(x, y)
    ok = bool(mx.b <= my.a)
    _log(label, "<=", x, y, mx, my, k, ok)
    return ok


def _nonzero(xs):
    return [x for x in xs if not x.is_zero()]


def smax(xs):
    """Enclosure of the maximum, at the smallest exponent (as certificate_iv._max)."""
    xs = [_co(x) for x in xs]
    nz = _nonzero(xs) or xs
    k = min(x.k for x in nz)
    box = next((x.box for x in xs if x.box is not None), None)
    ms = [x.m if (x.is_zero() or x.k == k) else x.m * box.pw(x.k - k) for x in xs]
    return S(iv.mpf([max(m.a for m in ms), max(m.b for m in ms)]), k, box)


def smin(xs):
    """Enclosure of the minimum.  With different exponents all terms must be positive; the
    minimum is then taken at the largest exponent K: x sigma^k = (x sigma^(k-K)) sigma^K with
    sigma^(k-K) in [sigma_hi^(k-K), sigma_lo^(k-K)] (infinite upper end if sigma_lo = 0)."""
    xs = [_co(x) for x in xs]
    ks = {x.k for x in xs}
    box = next((x.box for x in xs if x.box is not None), None)
    if len(ks) == 1:
        return S(iv.mpf([min(x.m.a for x in xs), min(x.m.b for x in xs)]), xs[0].k, box)
    if not all(x.m.a > 0 for x in xs):
        raise ArithmeticError("smin with different exponents needs positive terms")
    K = max(ks)
    lows, highs = [], []
    for x in xs:
        n = K - x.k
        if n == 0:
            lows.append(x.m.a)
            highs.append(x.m.b)
        else:
            if box.hi > 0:                      # at sigma = 0 alone this term is +infinite
                lows.append((iv.mpf(x.m.a) / box.pw_hi(n)).a)
            if box.lo > 0:
                highs.append((iv.mpf(x.m.b) / box.pw_lo(n)).b)
    return S(iv.mpf([min(lows), min(highs)]), K, box)


def shull(xs):
    """Interval hull, as phase_bound_iv builds D and D' from the frames."""
    xs = [_co(x) for x in xs]
    k = min(x.k for x in xs)
    box = next((x.box for x in xs if x.box is not None), None)
    ms = [x.m if x.k == k else x.m * box.pw(x.k - k) for x in xs]
    return S(iv.mpf([min(m.a for m in ms), max(m.b for m in ms)]), k, box)


def sabs(x):
    x = _co(x)
    return S(ci._abs(x.m), x.k, x.box)


def ssqrt(x):
    x = _co(x)
    if not x.m.a >= 0:
        raise ArithmeticError("square root of a number not certainly non-negative")
    return S(iv.sqrt(x.m), x.k / 2, x.box)


def collapse(x):
    """The ordinary interval m * sigma^k (needs k >= 0)."""
    x = _co(x)
    if x.k == 0 or x.is_zero():
        return x.m
    if x.k < 0:
        raise ArithmeticError("negative exponent: unbounded at sigma = 0")
    return x.m * x.box.pw(x.k)


def value(x):
    """m * sigma^k as an interval, also for k < 0 when the box does not contain 0."""
    x = _co(x)
    if x.k >= 0:
        return collapse(x)
    if not x.box.lo > 0:
        raise ArithmeticError("negative exponent and sigma = 0 in the box")
    return x.m / x.box.pw(-x.k)


def asin_upper(x):
    """phase_bound_iv._asin_upper applied to the ordinary interval of x (k >= 0)."""
    return S(pbi._asin_upper(collapse(x)))


# --------------------------------------------------------------------------------------
# Geometry of the square
# --------------------------------------------------------------------------------------
EDGE_ORDER = {"+-": "abmq", "-+": "qabm", "--": "ambq"}        # as certificate_iv._square
ENDS = {"pi": "angle = pi - d", "0": "angle = d"}
CASES = [("+-", "pi"), ("-+", "pi"), ("--", "pi"), ("--", "0")]


def geometry(region, end, box):
    """Edges, vertices and the scaled primitives of the square at distance d from the end."""
    if region not in EDGE_ORDER or end not in ENDS:
        raise ValueError("unknown region or end")
    if region != "--" and end != "pi":
        raise ValueError("the end angle -> 0 is degenerate only in the region (-,-)")
    sig = box.sig
    gam2 = 1 - sig * sig                       # cos^2(d/2)
    gam = iv.sqrt(gam2)
    sin_d, cos_d = 2 * sig * gam, 1 - 2 * sig * sig
    small, large = S(2, 2, box), S(2 * gam2, 0, box)            # 2 sigma^2, 2 cos^2(d/2)
    if end == "pi":
        sn, c = sin_d, -cos_d
        f_plus, f_minus = small, large          # 1 + cos(angle), 1 - cos(angle)
    else:
        sn, c = sin_d, cos_d
        f_plus, f_minus = large, small
    one, zero = iv.mpf(1), iv.mpf(0)
    vec = {"a": [one, zero, zero, one], "b": [one, zero, zero, -one],
           "m": [-one, -sn, zero, -c], "q": [-one, sn, zero, c]}
    names = EDGE_ORDER[region]
    p = [vec[ch] for ch in names]
    x = [[iv.mpf(0)] * 4]
    for k in range(3):
        x.append([x[-1][i] + p[k][i] for i in range(4)])
    centre = [sum((xi[i] for xi in x), iv.mpf(0)) / 4 for i in range(4)]
    Delta = S(ci._max([ci._enorm([xi[i] - centre[i] for i in range(4)]) for xi in x]), 0, box)
    Pn = S(ci._max([ci._enorm(pj) for pj in p]), 0, box)
    Dsigned, Bj, f_of = [], [], []
    for j in range(4):
        pair = {names[(j - 1) % 4], names[j]}
        pm, pp = p[(j - 1) % 4], p[j]
        if pair in ({"a", "b"}, {"m", "q"}):                    # 2 x.y = 4, no cancellation
            Dsigned.append(S(2 * ci._mdot(pm, pp), 0, box))
            Bj.append(S(ci._transverse_bound_planar(pm, pp), 0, box))
            f_of.append(None)
        else:
            f = f_plus if pair in ({"q", "a"}, {"b", "m"}) else f_minus
            Dsigned.append(-2 * f)                              # 2 x.y = -2 f
            Bj.append(ssqrt(4 / f))                             # 1 + |v|^2/(x.y)^2 = 4/f
            f_of.append(f)
    return {"p": p, "x": x, "centre": centre, "Delta": Delta, "Pn": Pn, "Dsigned": Dsigned,
            "Bj": Bj, "Dsmall": 2 * small, "sn": sn, "c": c}


def constants(g, fam=FAMILY):
    Dsmall = g["Dsmall"]                                        # |D'| = 4 sigma^2
    d0 = S(fam["c1"]) * Dsmall
    nu = S(fam["c2"]) * d0 * Dsmall
    return {"delta0": d0, "kappa": nu, "nu": nu, "kappa_far": S(fam["c3"]) * ssqrt(Dsmall),
            "R": S(fam["R"]), "tau0": S(fam["tau0"])}


# --------------------------------------------------------------------------------------
# certificate_iv.interval_certificate, scaled
# --------------------------------------------------------------------------------------
def scaled_certificate(region, end, box, fam=FAMILY, vals=None):
    """Conditions of GATE_6_CONTOUR section 9 for every square with sin(d/2) in the box,
    with the constants of the family."""
    g = geometry(region, end, box)
    Delta, Pn, Dsigned = g["Delta"], g["Pn"], g["Dsigned"]
    Ds = [sabs(dd) for dd in Dsigned]
    Dmin, Dmax = smin(Ds), smax(Ds)
    s, t = Dsigned[1], Dsigned[2]
    u_abs = sabs(s + t)
    cst = constants(g, fam)
    kap, nu, kf = cst["kappa"], cst["nu"], cst["kappa_far"]
    d0, R, tau0 = cst["delta0"], cst["R"], cst["tau0"]
    r17, r5 = S(iv.sqrt(17)), S(iv.sqrt(5))
    out = {"region": region, "end": end, "sigma": [float(box.lo), float(box.hi)]}

    A1 = 2 * Pn / Dmin
    A2 = 1 + 2 * A1 * Pn
    lam0 = 2 * fc.HESSIAN_BOUND * kap * A1 * Pn + nu * A2

    # H1 at eps = e_h1 |D'|
    eps = S(fam["e_h1"]) * g["Dsmall"]
    holds = False
    for j in range(4):
        D, Dp = Ds[j], Ds[(j + 1) % 4]
        rho = eps + 4 * eps**2 / D
        shift = 4 * eps * Dp / D + eps + 2 * ssqrt(Dp * u_abs / D * rho)
        holds = lt(shift, Dp - eps, f"H1 contradiction, frame {j}") or holds
    out["H1_gap_certified"] = holds
    out["H1_frame_on"] = lt(2 * d0, eps, "H1 2 delta0 < eps")

    out["C0_centres"] = le(3 * Delta, tau0, "C0 3 Delta <= tau0")
    margin = (tau0 - 2 * Delta) * (R / r17 - Delta) - (d0 + Delta * ssqrt(d0))
    out["C1_far_admissible"] = (lt(Delta, R / r17, "C1 Delta < R/sqrt17")
                                and lt(0, margin, "C1 margin > 0"))
    xi0 = 2 * (2 * Delta + ssqrt(d0))
    out["C2_theta_is_one"] = le(r5 / 2 * xi0, R, "C2")
    xi_min = R * 4 / r17
    worst = (xi_min / 4 + Delta) ** 2 - (xi_min - Delta) ** 2
    out["C3_wedge_free"] = lt(Delta, xi_min, "C3 Delta < xi_min") and lt(worst, -d0, "C3 worst < -delta0")

    reach = 2 * R + Delta
    frame_term = 2 * lam0 * reach**2
    v_inner = 4 * lam0 * reach + kf * (2 * R + tau0 + Delta)
    xi = R * 2 / r5
    out["C4_wedge_margin"] = lt(4 * Delta, xi, "C4")
    case_b = 1 + 2 * Delta / (3 * (R / r17 - Delta))
    case_c = 1 + (r5 / 2 + Delta / xi) * (tau0 + Delta) / (0.75 * (xi - 4 * Delta))
    e_bound = frame_term / (d0 - v_inner**2) + 4 * kf * smax([S(1), case_b, case_c])
    out["D_excursion_small"] = (lt(v_inner**2, d0 / 2, "D v^2 < delta0/2")
                                and lt(e_bound, 0.2, "D e < 0.2"))
    out["D_far_real_part"] = (lt(kf**2, 0.5, "D kf^2 < 1/2")
                              and lt(kf**2 * (2 * case_b - 1), 0.5, "D kf^2 (2 C_b - 1) < 1/2")
                              and lt(kf**2 * (r5 / 2 * xi + 2 * Delta + tau0) ** 2,
                                     0.375 * xi * (xi - 4 * Delta), "D third case"))

    j_inner = 4 * (3 * reach**2 * lam0 / d0 + lam0)
    j_outer = kf * (1 + 7.5 * (1 + (tau0 + Delta) / R))
    j_blend = 1.5 / R * (4 * lam0 * reach + kf * (2 * R + tau0 + Delta))
    out["E_local_diffeomorphism"] = lt(j_inner + j_outer + j_blend, 1, "E ||dV|| < 1")

    B = smax(g["Bj"])
    theta = smax([512 * kap * Pn**2 / Dmin, 8 * nu * B**2])
    out["F_theta_small"] = le(theta, 1, "F theta <= 1")

    b0 = 4 * d0 / Dmin
    a1 = d0 / (8 * Dmax)
    rho1 = (d0 / 8) ** 2 / smax([S(1), Dmax * u_abs / Dmin]) / 4
    out["S0_strips_overlap"] = lt(b0, 0.5, "S0 b0 < 1/2")
    out["S2_near_invariants"] = le(a1 * Dmax * (1 + a1) + rho1, d0, "S2")
    far_lower = Dmin * (b0 - a1) - a1 * Dmax - rho1 - 2 * ssqrt(rho1 * Dmax * u_abs / Dmin)
    out["S1_frame_fully_on"] = le(2 * d0, far_lower, "S1 far invariant >= 2 delta0")
    out["S3_strips_inside"] = le(Delta + a1 * Pn + B * ssqrt(rho1), R, "S3")

    omega = S(fam["omega"])
    c_u = (tau0 - 2 * Delta) * (1 / r17 - Delta / (2 * R)) - Delta * ssqrt(omega / (2 * R))
    out["I0_c_u_positive"] = lt(0, c_u, "I0 c_u > 0") and le(2 * omega, c_u, "I0 2 omega <= c_u")
    r0 = 2 * R
    lhs = kf**2 * (2 * (r0 + Delta) * (tau0 + Delta) + (tau0 + Delta) ** 2)
    out["I1_real_part"] = (le(lhs, omega * r0 / 4, "I1 a")
                           and le(2 * kf**2 * (tau0 + Delta), omega / 4, "I1 b"))
    out["I2_theta_is_one"] = (le(2 * (2 * Delta + ssqrt(omega * r0)), 2 / r5 * r0, "I2 a")
                              and lt(5 * omega / 4, r0, "I2 b"))

    mu = smax([64 * kap * Pn**2 / Dmin, nu * B**2])
    K = kf * (2 * R + tau0 + Delta)
    if margin.m.a > 0:
        mu_far = kf * (1 + (1.5 * d0 + (tau0 + Delta) ** 2) / (2 * margin))
        mu_q = smax([3 * mu + ssqrt(3 * mu), mu_far + ssqrt(3 * mu) * K**2 / (2 * kf * margin)])
        gamma = fc.PHASE_CORNER_MARGIN + 3 * e_bound
        g_iv = collapse(gamma)
        out["Q1_corner_excluded"] = (lt(gamma, iv.pi / 2, "Q1 gamma < pi/2")
                                     and lt(S(iv.sin(g_iv)) * mu_q, S(iv.cos(g_iv)), "Q1 tan(gamma) mu_Q < 1"))
        out["mu_Q_upper"] = float(collapse(mu_q).b)
    else:
        out["Q1_corner_excluded"] = False
    out["Q3_corner_excluded_far"] = lt(2 * kf**2 * Delta**2 / (1 - kf**2),
                                       (2 * R / r17 - Delta) * (tau0 - 2 * Delta), "Q3")
    n0_min = R / r17 - Delta
    out["Q4_outer_zone"] = (lt(kf * (tau0 + Delta), (1 - kf) * n0_min, "Q4 a")
                            and lt(8 * kf**2 * Delta / (3 * (1 - kf**2)), 2 * n0_min, "Q4 b"))

    out["excursion_upper"] = float(collapse(e_bound).b)
    out["jacobian_upper"] = float(collapse(j_inner + j_outer + j_blend).b)
    out["all"] = all(v for v in out.values() if isinstance(v, bool))
    if vals is not None:
        vals.update({k: v for k, v in locals().items() if isinstance(v, S)})
        vals.update(Dsigned=Dsigned, Ds=Ds)
    return out


# --------------------------------------------------------------------------------------
# phase_bound_iv.interval_phase_bound, scaled
# --------------------------------------------------------------------------------------
def scaled_phase_bound(region, end, box, fam=FAMILY, with_certificate=True, vals=None):
    """Conditions of GATE_6_PHASE section 6 and an upper bound for |arg z| - 2 pi, for every
    square with sin(d/2) in the box, with the constants of the family."""
    if region == "--":
        raise ValueError("no timelike diagonal in the region (-,-)")
    g = geometry(region, end, box)
    Delta, Pn, Dsigned = g["Delta"], g["Pn"], g["Dsigned"]
    tl = [j for j in range(4) if Dsigned[j].m.a > 0]
    sl = [j for j in range(4) if Dsigned[j].m.b < 0]
    if len(tl) != 2 or len(sl) != 2:
        raise ArithmeticError("signs of the diagonals not certain")
    D = shull([Dsigned[j] for j in tl])
    Dp = shull([-Dsigned[j] for j in sl])
    Dmin, Dmax = smin([D, Dp]), smax([D, Dp])
    u_abs = sabs(D - Dp)
    cst = constants(g, fam)
    kap, nu, kf = cst["kappa"], cst["nu"], cst["kappa_far"]
    d0, R, tau0 = cst["delta0"], cst["R"], cst["tau0"]
    r2, r5, r17 = S(iv.sqrt(2)), S(iv.sqrt(5)), S(iv.sqrt(17))
    out = {"region": region, "end": end, "sigma": [float(box.lo), float(box.hi)]}
    out["certificate"] = bool(scaled_certificate(region, end, box, fam)["all"]) if with_certificate else True

    A1 = 2 * Pn / Dmin
    A2 = 1 + 2 * A1 * Pn
    lam0 = 2 * fc.HESSIAN_BOUND * kap * A1 * Pn + nu * A2
    M = (tau0 - 2 * Delta) * (R / r17 - Delta) - (d0 + Delta * ssqrt(d0))
    reach = 2 * R + Delta
    K = kf * (2 * R + tau0 + Delta)
    v_in = 4 * lam0 * reach
    v2 = (v_in + K) ** 2
    xi = R * 2 / r5
    if not (lt(Delta, R / r17) and lt(4 * Delta, xi) and lt(v2, d0, "pre v^2 < delta0") and lt(0, M)
            and lt(2 * Delta, tau0) and lt(Delta, R)):
        raise ArithmeticError("a denominator of section 9 is not certainly positive")
    case_b = 1 + 2 * Delta / (3 * (R / r17 - Delta))
    case_c = 1 + (r5 / 2 + Delta / xi) * (tau0 + Delta) / (0.75 * (xi - 4 * Delta))
    e = 2 * lam0 * reach**2 / (d0 - v2) + 4 * kf * smax([S(1), case_b, case_c])
    Bj = g["Bj"]
    B, B_c, B_n = smax(Bj), smax([Bj[j] for j in tl]), smax([Bj[j] for j in sl])
    theta = smax([512 * kap * Pn**2 / Dmin, 8 * nu * B**2])
    mu = smax([64 * kap * Pn**2 / Dmin, nu * B**2])

    # outer zone (R4), (R7)
    q = (tau0 + Delta) / (tau0 - 2 * Delta)
    eps_k = kf**2 * (1 + 2 * r2 * q + q * q)
    out["R4_centres"] = le(3 * Delta, tau0, "R4 3 Delta <= tau0") and lt(r2 * tau0, 2 * R, "R4 sqrt2 tau0 < 2R")
    out["R4_trivial"] = lt(eps_k, 1, "R4 eps_kappa < 1") and le(e, iv.pi / 4, "R4 e <= pi/4")
    if not out["R4_trivial"]:
        raise ArithmeticError("eps_kappa or e too large")
    e_sp = 2 * kf * (1 + r2 * q) / (1 - eps_k)
    K_out = smax([e, 3 * e_sp])
    out["R4_trivial"] = out["R4_trivial"] and lt(4 * e_sp, iv.pi, "R4 4 e_sp < pi")
    out["R7_outer_exact"] = out["R4_centres"] and le(2 * Delta * tau0, (R - tau0 / 2) * (tau0 - 2 * Delta), "R7")

    # inner zone (R5)
    out["G1"] = lt(4 * d0 * u_abs, D * Dp, "G1 4 delta0 |u| < D |D'|")
    out["G2"] = lt(4 * d0, D, "G2 4 delta0 < D")
    iota = 2 * lam0 * (R + Delta) ** 2
    v_R2 = (4 * lam0 * (R + Delta)) ** 2
    out["R5_trivial"] = (lt(v_R2, d0, "R5 v_R^2 < delta0") and le(theta, 2, "R5 theta <= 2")
                         and lt(nu, 1, "R5 nu < 1") and lt(4 * kap**2, 1, "R5 4 kappa^2 < 1"))
    if not out["R5_trivial"]:
        raise ArithmeticError("a trivial condition of R5 fails")
    e4 = 4 * iota / (d0 - v_R2)
    beta_c = 2 * Pn * (R + Delta) / D
    G = 2 * D / d0 * smax([S(1), beta_c])
    w = u_abs / Dp
    c_3d = G * (2 * nu * (1 + beta_c * w) + 2 * kap * (beta_c * D / Dp + u_abs / D + u_abs**2 * beta_c / (Dp * D)))
    eps_V = smax([asin_upper(theta / 2), asin_upper(v_R2 / d0)])
    K_in = smax([2 * e + e4, e + 2 * eps_V, (2 * nu + 2 * c_3d + 2 * iota / d0) + 2 * iota / d0 + 3 * eps_V])
    e_fr = iota / (d0 - v_R2)
    K_3d = (2 * nu + 2 * c_3d + 2 * iota / d0) + 2 * iota / d0 + 3 * eps_V
    b_in = e_fr + smax([2 * e_fr + e4, e_fr + 2 * eps_V, K_3d])

    # blend zone (R6)
    iota2 = 2 * lam0 * reach**2
    pi0 = 2 * lam0 * reach * Pn
    mu_m = 3 * mu + K * ssqrt(3 * mu / (2 * kf * M))
    mu_p = mu_m + K**2 / (2 * kf * M)
    omega0 = v_in * (2 * K + v_in)
    dv = d0 - v2
    eps = r2 * Pn / (R - Delta)
    b0 = eps * d0 / Dmin
    Lam_c = (R - Delta - B_c * ssqrt(d0 + b0**2 * D)) / Pn
    Lam_n = (R - Delta - B_n * ssqrt(d0 * (1 + eps) + 2 * b0**2 * Dp)) / Pn
    a_min = Lam_c - b0
    d1 = d0 * (1 - eps)
    out["B1"] = le(tau0, R, "B1 tau0 <= R") and lt(r5 * Delta, R, "B1 sqrt5 Delta < R")
    out["B2"] = le(kap, nu, "B2 kappa <= nu")
    out["B3"] = (le(b0, 0.5, "B3 b0 <= 1/2") and le(2, Lam_c, "B3 Lambda_c >= 2")
                 and le(2, Lam_n, "B3 Lambda_n >= 2") and le(4 * d0, Dmin, "B3 4 delta0 <= min|D|"))
    out["B4"] = lt(v2, d1, "B4 v^2 < delta1")
    gap = (tau0 - 2 * Delta) ** 2
    out["R6_trivial"] = lt(v2, d0 / 2, "R6 v^2 < delta0/2") and lt(v2, gap, "R6 v^2 < (tau0-2Delta)^2")
    tau1 = 2 * Delta * tau0 / (tau0 - 2 * Delta)
    W1 = Delta * (tau0 + Delta) + Delta * (tau1 + tau0 + Delta)
    out["B6"] = False
    rho1 = None
    if lt(tau1, R, "B6 tau1 < R"):
        root = ssqrt(R * R - tau1 * tau1) - Delta
        if lt(tau1 + Delta, root, "B6 second"):
            rho1 = root * root - (tau1 + Delta) ** 2
            out["B6"] = lt(0, rho1, "B6 rho1 > 0")
    out["B5"] = False
    if out["B3"] and out["B4"] and out["B6"] and out["R6_trivial"] and a_min.m.a > 0:
        unit = iota2 / dv
        half = S(PAIRING) * unit
        a_c = 2 * Pn * reach / D
        a_n = 2 * Pn * reach / Dp
        vs_c = 2 * nu + 2 * D / d0 * (kap * a_c**2 + nu * a_c**2 + kap * a_c * smax([a_c + 1, S(2)]))
        vs_n = 2 * Dp / d0 * (kap * a_n**2 + nu * (a_n**2 + a_n) + kap * a_n * (a_n + 1))
        far0 = 2 * K * Pn * omega0 / dv**2
        vterm = 2 * v2 * pi0 / dv**2
        s0 = smax([3 * half + vterm + far0, 2 * half + far0,
                   unit + 2 * half + (vs_c + vs_n) * (d0 / dv) ** 2 + vterm + far0])
        vs = (2 * nu + kap * b0**2 / (a_min * (1 + a_min))) + (2 * nu + kap * b0**2) / a_min
        F = 2 * K * Pn * mu_m / dv
        chi_eps = eps * eps * (3 - 2 * eps)
        s1 = smax([half + 2 * mu_p * pi0 / dv + F,
                   (1 + mu_p) * d0 / dv * vs + 2 * mu_m * pi0 / dv + chi_eps * iota2 * (1 / dv + 1 / (d1 - v2)) + F])
        s2 = half + 2 * mu_p * pi0 / dv + F
        s3 = mu_m + (1 + mu_p) ** 2 * vs
        s = smax([s0, s1, s2, s3])
        if lt(s, 1, "B5 s < 1"):
            asin_s = asin_upper(s)
            e_sp2 = (4 * iota2 / gap + 2 * kf * (1 + r2 * q)) / (1 - v2 / gap)
            out["B5"] = (lt(e + asin_s + mu_m, iv.pi / 2, "B5 e + arcsin s + mu_- < pi/2")
                         and lt(4 * e_sp2, iv.pi, "B5 4 e'_sp < pi"))
            e_bl_hi = half
            e_bl_lo = half + 2 * kf * (1 + W1 / rho1) / (1 - v2 / d0)
            b_bl = e_bl_hi + e_bl_lo + asin_s
            b_out = S(0) if out["R7_outer_exact"] else e + K_out
            bound = smax([b_in, b_bl, b_out])
            out.update(s_upper=float(collapse(s).b), b_bl_upper=float(collapse(b_bl).b),
                       bound_upper=float(collapse(bound).b))
    out.update(e_upper=float(collapse(e).b), K_in_upper=float(collapse(K_in).b),
               b_in_upper=float(collapse(b_in).b), K_out_upper=float(collapse(K_out).b),
               c_3d_upper=float(collapse(c_3d).b),
               Lambda_lower=float(min(collapse(Lam_c).a, collapse(Lam_n).a)),
               Bn2_Dp_upper=float(collapse(B_n**2 * Dp).b))
    out["all_conditions"] = all(out[k] for k in pbi.CONDITIONS)
    if "bound_upper" not in out:
        out["bound_upper"] = float("inf")
    if vals is not None:
        vals.update({k: v for k, v in locals().items() if isinstance(v, S)})
        vals.update(Dsigned=Dsigned, Bj=Bj)
    return out


def check_box(region, end, box, fam=FAMILY):
    """Everything the graded cover requires on one box: the certificate (all regions) and,
    with a timelike diagonal, all conditions of the phase bound and (R7)."""
    cert = scaled_certificate(region, end, box, fam)
    res = {"region": region, "end": end, "sigma": cert["sigma"], "certificate": cert["all"],
           "cert": cert}
    if region == "--":
        res["ok"] = cert["all"]
        return res
    try:
        ph = scaled_phase_bound(region, end, box, fam, with_certificate=False)
    except (ArithmeticError, ZeroDivisionError) as error:
        res.update(ok=False, phase=None, error=str(error))
        return res
    res["phase"] = ph
    res["ok"] = bool(cert["all"] and ph["all_conditions"] and ph["R7_outer_exact"])
    return res


def scaled_cover(region, end, d_hi, boxes=1, d_lo=0.0, fam=FAMILY, geometric=False):
    """check_box on `boxes` adjacent sigma-boxes covering sin(d/2) for d in [d_lo, d_hi].
    With `geometric` the boxes after the first grow geometrically (d_lo must be 0)."""
    full = box_from_d(d_lo, d_hi)
    lo, hi = full.lo, full.hi                 # exact endpoints of the enclosure of sin(d/2)
    flo, fhi = float(lo), float(hi)
    if geometric and boxes > 1:
        inner = [fhi * 10.0 ** (i - boxes + 1) for i in range(boxes - 1)]
    else:
        inner = [flo + (fhi - flo) * i / boxes for i in range(1, boxes)]
    # adjacent boxes share their endpoints exactly, the first starts at lo and the last ends at hi
    edges = [lo] + [x for x in inner if flo < x < fhi] + [hi]
    failed, worst = [], {"bound_upper": 0.0, "excursion_upper": 0.0, "mu_Q_upper": 0.0,
                         "b_in_upper": 0.0, "b_bl_upper": 0.0, "jacobian_upper": 0.0}
    lam = float("inf")
    for a, b in zip(edges[:-1], edges[1:]):
        try:
            r = check_box(region, end, Box(a, b), fam)
        except (ArithmeticError, ZeroDivisionError) as error:
            failed.append((float(a), float(b), str(error)))
            continue
        if not r["ok"]:
            bad = [k for k, v in r["cert"].items() if v is False]
            if r.get("phase"):
                bad += [k for k in pbi.CONDITIONS + ["R7_outer_exact"] if not r["phase"][k]]
            failed.append((float(a), float(b), bad or r.get("error")))
            continue
        for k in ("excursion_upper", "mu_Q_upper", "jacobian_upper"):
            worst[k] = max(worst[k], r["cert"].get(k, 0.0))
        if r.get("phase"):
            for k in ("bound_upper", "b_in_upper", "b_bl_upper"):
                worst[k] = max(worst[k], r["phase"][k])
            lam = min(lam, r["phase"]["Lambda_lower"])
    return {"region": region, "end": end, "d": [float(d_lo), float(d_hi)],
            "sigma": [float(lo), float(hi)], "boxes": boxes, "all": not failed,
            "failed": failed[:5], "failed_count": len(failed), "Lambda_lower": lam, **worst}


# --------------------------------------------------------------------------------------
# The scaled four-cone gap
# --------------------------------------------------------------------------------------
def scaled_gap(region, end, box, tol=1e-9):
    """Largest e (to tol) for which the contradiction of (H1) at eps = e |D'| is certain on the
    whole box, maximised over the four frames: the scaled analogue of
    feynman_contour.four_cone_gap (which bisects eps in [0, |D'|] in float64)."""
    g = geometry(region, end, box)
    Ds = [sabs(dd) for dd in g["Dsigned"]]
    u_abs = sabs(g["Dsigned"][1] + g["Dsigned"][2])
    best = 0.0
    for j in range(4):
        D, Dp = Ds[j], Ds[(j + 1) % 4]

        def holds(e):
            eps = S(e) * g["Dsmall"]
            rho = eps + 4 * eps**2 / D
            shift = 4 * eps * Dp / D + eps + 2 * ssqrt(Dp * u_abs / D * rho)
            return lt(shift, Dp - eps)

        lo, hi = 0.0, 1.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if holds(mid) else (lo, mid)
        best = max(best, lo)
    return best


# --------------------------------------------------------------------------------------
# The repository code at one square, with the constants of the family
# --------------------------------------------------------------------------------------
def family_contour(region, end, angle, fam=FAMILY):
    """A feynman_contour.Contour at the float `angle`, with the family's constants rounded
    to float64."""
    d = (iv.pi - iv.mpf(angle)) if end == "pi" else iv.mpf(angle)
    sig = iv.sin(d / 2)
    box = Box(sig.a, sig.b)
    cst = constants(geometry(region, end, box), fam)
    mid = lambda x: float(mp.mpf(collapse(x).mid))
    return fc.Contour(kappa=mid(cst["kappa"]), nu=mid(cst["nu"]), kappa_far=mid(cst["kappa_far"]),
                      delta0=mid(cst["delta0"]), radius=fam["R"], tau0=fam["tau0"],
                      angle=float(angle), region=region), box


def repository_reference(region, end, angle, fam=FAMILY, exact=False):
    """Run certificate_iv.interval_certificate and phase_bound_iv.interval_phase_bound on the
    one-point box [angle, angle] with feynman_contour.tuned_contour temporarily replaced, and
    capture their local variables.  Nothing on disk is changed.

    exact=False: a real feynman_contour.Contour with the family's constants rounded to float64
    (the repository code then also runs its float certificate for eps and omega).
    exact=True: the constants are passed as the family's intervals (certificate_iv._I passes
    intervals through), and feynman_contour.certificate is replaced so that the repository's
    eps = 0.99 * gap and omega = 0.98 * omega enclose the family's e_h1 |D'| and omega; then
    every intermediate must agree with the scaled evaluation to the working precision."""
    C, box = family_contour(region, end, angle, fam)
    patches = {"tuned_contour": lambda region_, angle_, *a, **k: C}
    if exact:
        g = geometry(region, end, box)
        cst = constants(g, fam)
        C = type("ExactContour", (), {})()
        C.kappa, C.nu, C.kappa_far = (collapse(cst[k]) for k in ("kappa", "nu", "kappa_far"))
        C.delta0, C.radius, C.tau0 = collapse(cst["delta0"]), iv.mpf(fam["R"]), iv.mpf(fam["tau0"])
        eps = collapse(S(fam["e_h1"]) * g["Dsmall"])
        patches = {"tuned_contour": lambda region_, angle_, *a, **k: C,
                   "certificate": lambda C_: {"four_cone_gap": eps / 0.99, "omega": iv.mpf(fam["omega"]) / 0.98}}
    codes = {ci.interval_certificate.__code__: "certificate", pbi.interval_phase_bound.__code__: "phase"}
    captured = {}

    def prof(frame, event, arg):
        if event == "return" and frame.f_code in codes:
            captured[codes[frame.f_code]] = dict(frame.f_locals)

    saved = {name: getattr(fc, name) for name in patches}
    for name, f in patches.items():
        setattr(fc, name, f)
    sys.setprofile(prof)
    try:
        cert = ci.interval_certificate(region, float(angle), float(angle))
        phase = pbi.interval_phase_bound(region, float(angle), float(angle)) if region != "--" else None
    finally:
        sys.setprofile(None)
        for name, f in saved.items():
            setattr(fc, name, f)
    return {"certificate": cert, "phase": phase, "locals": captured, "box": box, "contour": C}


def compare_with_repository(region, end, d, fam=FAMILY, exact=False):
    """Every named intermediate of the scaled evaluation at the one-point sigma-box of d against
    the repository's local variable of the same name.  Returns the largest relative difference
    of the midpoints, the number compared, the names not found, boolean mismatches, and the
    repository locals that have no counterpart here (should be only 'angle')."""
    angle = float(mp.pi - mp.mpf(d)) if end == "pi" else float(d)
    ref = repository_reference(region, end, angle, fam, exact=exact)
    box = ref["box"]
    mine = {"certificate": {}, "phase": {}}
    outs = {"certificate": scaled_certificate(region, end, box, fam, vals=mine["certificate"])}
    if region != "--":
        outs["phase"] = scaled_phase_bound(region, end, box, fam, vals=mine["phase"])
    g = geometry(region, end, box)
    skip = set() if exact else {"eps", "rho", "shift", "omega", "c_u"}       # depend on eps, omega
    worst, where, count, missing, bool_mismatch, overlap_fail = 0.0, None, 0, [], [], []
    for key, out in outs.items():
        L = ref["locals"][key]
        ours = dict(mine[key])
        if key == "certificate":
            ours = {k: v for k, v in ours.items() if k not in skip}
        ours["centre"] = [S(c, 0, box) for c in g["centre"]]
        for name, v in ours.items():
            if name not in L:
                missing.append((key, name))
                continue
            mv = v if isinstance(v, list) else [v]
            tv = L[name] if isinstance(L[name], list) else [L[name]]
            for a, b in zip(mv, tv):
                a, b = value(a), ci._I(b)
                with mp.workdps(iv.dps + 10):
                    ma, mb = mp.mpf(a.mid), mp.mpf(b.mid)
                    rel = float(abs(ma - mb) / max(abs(mb), mp.mpf(10) ** -300))
                count += 1
                if rel > worst:
                    worst, where = rel, (key, name)
                if exact and not (a.a <= b.b and b.a <= a.b):
                    overlap_fail.append((key, name))
        R = ref[key]
        for k2, v2 in out.items():
            if isinstance(v2, bool) and k2 in R and R[k2] != v2 and (exact or k2 not in
                                                                    ("H1_gap_certified", "H1_frame_on")):
                bool_mismatch.append((key, k2, v2, R[k2]))
        unmatched = sorted(k for k, v in L.items()
                           if (isinstance(v, iv.mpf) or (isinstance(v, list) and v and isinstance(v[0], iv.mpf)))
                           and k not in ours)
    return {"worst_rel": worst, "where": where, "count": count, "missing": missing,
            "bool_mismatch": bool_mismatch, "overlap_fail": overlap_fail, "unmatched": unmatched,
            "repository_all": [ref["certificate"]["all"]] + ([ref["phase"]["all_conditions"]] if ref["phase"] else []),
            "scaled_all": [outs["certificate"]["all"]] + ([outs["phase"]["all_conditions"]] if "phase" in outs else [])}


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------
def rescaling_table(region, end, box, fam=FAMILY):
    """Every labelled comparison of the two transcriptions on the box: the exponents of the
    two sides, the exponent divided out, and the two scaled sides."""
    global RECORD
    RECORD = []
    try:
        scaled_certificate(region, end, box, fam)
        if region != "--":
            scaled_phase_bound(region, end, box, fam, with_certificate=False)
        rows = RECORD
    finally:
        RECORD = None
    return rows


def _fmt(x):
    return f"[{float(x.a):.6g}, {float(x.b):.6g}]"


def _cover_line(r):
    line = (f"{r['region']} end {r['end']}, d in {r['d']}, {r['boxes']} box(es): all={r['all']}  "
            f"e<={r['excursion_upper']:.5f}  mu_Q<={r['mu_Q_upper']:.3g}  ||dV||<={r['jacobian_upper']:.4f}")
    if r["region"] != "--":
        line += (f"  bound-2pi<={r['bound_upper']:.5f}  b_in<={r['b_in_upper']:.5f}  "
                 f"b_bl<={r['b_bl_upper']:.5f}  Lambda>={r['Lambda_lower']:.4f}")
    return line + (f"  failed({r['failed_count']})={r['failed'][:2]}" if r["failed"] else "")


def main(d0=1e-6):
    print(f"family: {FAMILY}")
    print(f"\n1. One box d in [0, {d0}] (sigma in [0, sin(d0/2)]) per case")
    for region, end in CASES:
        print("  " + _cover_line(scaled_cover(region, end, d0, boxes=1)))
    print("\n2. The limit d = 0 (sigma = [0, 0]), scaled values")
    for region, end in CASES:
        box = Box(0, 0)
        vals = {}
        c = scaled_certificate(region, end, box, vals=vals)
        line = (f"  {region} end {end}: cert all={c['all']}, scaled gap e_* = "
                f"{scaled_gap(region, end, box):.9f} (1 - sqrt3/2 = {1 - 3 ** 0.5 / 2:.9f}), "
                f"Delta={_fmt(vals['Delta'].m)}, e={c['excursion_upper']:.6f}, ||dV||={c['jacobian_upper']:.6f}")
        if region != "--":
            pv = {}
            ph = scaled_phase_bound(region, end, box, vals=pv)
            line += (f"\n      phase: all={ph['all_conditions']} R7={ph['R7_outer_exact']} "
                     f"b_in={ph['b_in_upper']:.6f} b_bl={ph['b_bl_upper']:.6f} s={ph['s_upper']:.6f} "
                     f"c_3d={ph['c_3d_upper']:.6f} Lambda={ph['Lambda_lower']:.4f} B_n^2|D'|<={ph['Bn2_Dp_upper']:.6f}")
        print(line)
    print("\n3. Rescaling table, (+,-) at the box d in [0, d0]: label, k_lhs, k_rhs, k divided out, scaled sides")
    for row in rescaling_table("+-", "pi", box_from_d(0, d0)):
        print(f"  {row['label']:<36} {str(row['k_lhs']):>5} {str(row['k_rhs']):>5} {str(row['k']):>5}  "
              f"{_fmt(row['lhs'])} {row['op']} {_fmt(row['rhs'])}  {row['ok']}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--d0", type=float, default=1e-6)
    ap.add_argument("--extend", type=float, default=None,
                    help="also cover d in [0, EXTEND] with geometric boxes")
    ap.add_argument("--boxes", type=int, default=400)
    args = ap.parse_args()
    main(args.d0)
    if args.extend:
        print(f"\n4. Extended: d in [0, {args.extend}], {args.boxes} equal sigma-boxes")
        for region, end in CASES:
            print("  " + _cover_line(scaled_cover(region, end, args.extend, boxes=args.boxes)))
