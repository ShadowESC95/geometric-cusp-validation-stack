"""All four cusps, general kinematics, and weights on the cover of the z plane.

Extends `geocusp.lorentzian_cusp` (which treats cusp 2 with spacelike
diagonals) so that the statements of docs/GATE_6_CLOSURE.md can be computed,
not inferred, at every cusp and in every real kinematic region of the
light-like square. See docs/GATE_6_BOUNDS.md.

Conventions. Metric (+,-,-,-). Edges p_i = x_{i+1} - x_i, all null, summing
to zero. s = x13^2, t = x24^2. Near the cusp x_i write

    x0 = x_i + alpha p_{i-1} + beta p_i + yperp,   rho = |yperp|^2,
    D  = 2 p_{i-1}.p_i   (= s at cusps 2, 4 and t at cusps 1, 3),
    D' = the other diagonal.

Then exactly
    x0,i^2   = w,                 w = alpha beta D - rho,
    x0,i-1^2 = w + c_prev,        c_prev = beta D,
    x0,i+1^2 = w + c_next,        c_next = -alpha D,
    x0,i+2^2 = D'(1 + alpha - beta) + w + 2 sqrt(rho D' u / D) cos(theta),

with u = -s - t. Every invariant, including s and t, carries -i0, so

    arg z = pi (tau_2 + tau_4 + tau_s - tau_1 - tau_3 - tau_t),

tau = 1 for a timelike (positive) invariant and 0 for a spacelike one.
"""
from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp
import numpy as np

METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def msq(v):
    v = np.asarray(v, dtype=float)
    return float(v @ METRIC @ v)


def mdot(a, b):
    return float(np.asarray(a, dtype=float) @ METRIC @ np.asarray(b, dtype=float))


# --------------------------------------------------------------------------
# Explicit light-like squares in the three real kinematic regions
# --------------------------------------------------------------------------
def light_like_square(region: str = "--", angle: float = 1.1):
    """Null edges p1..p4 (rows) with s = (p1+p2)^2, t = (p2+p3)^2.

    region is the pair of signs of (s, t): "--", "+-" or "-+". Real null
    vectors summing to zero have exactly one positive Mandelstam invariant,
    so "++" does not exist.
    """
    c, sn = np.cos(angle), np.sin(angle)
    if region == "--":
        p = np.array([[1.0, 0, 0, 1.0], [-1.0, -sn, 0, -c], [1.0, 0, 0, -1.0], [-1.0, sn, 0, c]])
    elif region == "+-":
        p = np.array([[1.0, 0, 0, 1.0], [1.0, 0, 0, -1.0], [-1.0, -sn, 0, -c], [-1.0, sn, 0, c]])
    elif region == "-+":
        p = np.array([[-1.0, sn, 0, c], [1.0, 0, 0, 1.0], [1.0, 0, 0, -1.0], [-1.0, -sn, 0, -c]])
    else:
        raise ValueError("region must be '--', '+-' or '-+'")
    return p, msq(p[0] + p[1]), msq(p[1] + p[2])


def _transverse_basis(a, b):
    basis = []
    dot = mdot(a, b)
    for seed in np.eye(4):
        v = seed - mdot(seed, b) / dot * a - mdot(seed, a) / dot * b
        for e in basis:
            v = v + mdot(v, e) * e
        n = -msq(v)
        if n > 1e-12:
            basis.append(v / np.sqrt(n))
        if len(basis) == 2:
            return basis
    raise RuntimeError("no transverse basis")


def exact_invariants(cusp: int, alpha, beta, rho, theta, region="--", angle=1.1):
    """The four x0j^2 from explicit vectors, and the local formulas for them."""
    p, s, t = light_like_square(region, angle)
    i = cusp - 1
    pm, pi_, pn = p[(i - 1) % 4], p[i], p[(i + 1) % 4]
    e1, e2 = _transverse_basis(pm, pi_)
    # orient the transverse frame so that theta = 0 points along -p_{i+1} transverse part
    far_dir = np.array([-mdot(pn, e1), -mdot(pn, e2)])     # components of p_{i+1} on (e1, e2)
    phi0 = np.arctan2(far_dir[1], far_dir[0])
    y = (alpha * pm + beta * pi_
         + np.sqrt(rho) * (np.cos(theta + phi0) * e1 + np.sin(theta + phi0) * e2))
    pos = {i: np.zeros(4), (i + 1) % 4: pi_, (i + 2) % 4: pi_ + pn, (i - 1) % 4: -pm}
    exact = {j + 1: msq(y - pos[j]) for j in range(4)}
    D = 2 * mdot(pm, pi_)
    Dp = t if abs(D - s) < 1e-9 else s
    u = -s - t
    w = alpha * beta * D - rho
    local = {
        cusp: w,
        (i - 1) % 4 + 1: w + beta * D,
        (i + 1) % 4 + 1: w - alpha * D,
        (i + 2) % 4 + 1: Dp * (1 + alpha - beta) + w + 2 * np.sqrt(rho * Dp * u / D) * np.cos(theta),
    }
    return exact, local, {"s": s, "t": t, "D": D, "Dp": Dp, "u": u}


# --------------------------------------------------------------------------
# Sheet bookkeeping
# --------------------------------------------------------------------------
def lg(x):
    """log(x - i0): argument -pi for negative real x, 0 for positive real x.

    A complex x must lie in the closed lower half plane.
    """
    if isinstance(x, (mp.mpc, complex)) and mp.im(x) != 0:
        if mp.im(x) > 0:
            raise ValueError("invariant left the lower half plane")
        return mp.log(x)
    x = mp.re(x)
    return mp.log(abs(x)) - (1j * mp.pi if x < 0 else 0)


@dataclass(frozen=True)
class Cusp:
    index: int
    s: float = -1.0
    t: float = -1.0

    @property
    def D(self):
        return self.s if self.index % 2 == 0 else self.t

    @property
    def Dfar(self):
        return self.t if self.index % 2 == 0 else self.s

    @property
    def sigma(self):
        """+1 if x0,cusp^2 is in the numerator of z (cusps 2, 4), else -1."""
        return 1 if self.index % 2 == 0 else -1

    def far_edge_value(self, c_prev, c_next):
        """x0,i+2^2 on the plane yperp = 0, exactly: D'(1 + alpha - beta) + alpha beta D.

        With beta = c_prev/D and alpha = -c_next/D this is
        D'(1 - c_prev/D - c_next/D) - c_prev c_next/D. On the edge (alpha = 0)
        it is D'(1 - beta).
        """
        return self.Dfar * (1 - c_prev / self.D - c_next / self.D) - c_prev * c_next / self.D

    def log_z(self, w, c_prev, c_next, far):
        near = lg(w) + lg(far) - lg(w + c_prev) - lg(w + c_next)
        return self.sigma * near + lg(self.s) - lg(self.t)

    def w_start(self, c_prev, c_next):
        return -c_prev * c_next / self.D

    def kernel(self, f, c_prev, c_next, far=None, tilt=0.6):
        """int_0^oo d rho f(L) / [w (w + c_prev)(w + c_next)], w = w_start - rho - i0."""
        c_prev, c_next = mp.mpf(c_prev), mp.mpf(c_next)
        if far is None:
            far = self.far_edge_value(c_prev, c_next)
        w0 = self.w_start(c_prev, c_next)
        direction = -mp.exp(1j * mp.mpf(tilt))

        def integrand(r):
            w = w0 + r * direction
            return (f(self.log_z(w, c_prev, c_next, far))
                    / (w * (w + c_prev) * (w + c_next)) * (-direction))

        marks = sorted({mp.mpf(0), abs(w0), abs(c_prev), abs(c_next), mp.mpf(1), mp.inf})
        return mp.quad(integrand, marks)

    def edge_jump(self, f, edge: str, position, ratio=mp.mpf("1e-12")):
        """Density of the collinear logarithm at a point of a null line.

        edge = "next": the line through x_i and x_{i+1}, at x_i + position*p_i.
        edge = "prev": the line through x_{i-1} and x_i, at x_i - position*p_{i-1}.
        In both cases the segment is 0 < position < 1. Returns the value of
        c_prev c_next K at (collinear invariant > 0) minus (collinear invariant < 0),
        with the far propagator at its exact edge value D'(1 - position).
        """
        position = mp.mpf(position)
        large = position * self.D                      # beta D, or -alpha D with alpha = -position
        small = abs(large) * ratio
        far = self.Dfar * (1 - position)
        if edge == "next":
            plus = large * small * self.kernel(f, large, small, far)
            minus = large * (-small) * self.kernel(f, large, -small, far)
        elif edge == "prev":
            plus = small * large * self.kernel(f, small, large, far)
            minus = (-small) * large * self.kernel(f, -small, large, far)
        else:
            raise ValueError("edge must be 'next' or 'prev'")
        return plus - minus

    def quadrant_sum(self, f, c_prev, c_next):
        c_prev, c_next = abs(mp.mpf(c_prev)), abs(mp.mpf(c_next))
        return sum(c_prev * c_next * self.kernel(f, a * c_prev, b * c_next)
                   for a in (1, -1) for b in (1, -1))

    def sheet_sequence(self, c_prev, c_next, eta=1e-9, n=300, reach=0.02):
        """Distinct values of arg z / pi along rho in [0, reach * max(|c|)]."""
        c_prev, c_next = mp.mpf(c_prev), mp.mpf(c_next)
        far = self.far_edge_value(c_prev, c_next)
        w0 = self.w_start(c_prev, c_next)
        small = min(abs(c_prev), abs(c_next))
        scale = min(abs(w0), small)
        top = reach * max(abs(c_prev), abs(c_next))
        rhos = [mp.mpf(0)] + [mp.mpf(10) ** x for x in
                              np.linspace(float(mp.log10(scale)) - 3, float(mp.log10(top)), n)]
        out = []
        for r in rhos:
            w = w0 - r
            near = min(abs(w), abs(w + c_prev), abs(w + c_next))
            L = self.log_z(w - 1j * eta * max(near, scale * mp.mpf("1e-3")), c_prev, c_next, far)
            value = float(mp.im(L) / mp.pi)
            k = round(value)
            if abs(value - k) < 0.02 and (not out or out[-1] != k):
                out.append(int(k))
        return out


def circle_mean(f, radius=1):
    ell = mp.log(radius)
    return mp.quad(lambda p: f(ell + 1j * p), mp.linspace(-mp.pi, mp.pi, 9)) / (2 * mp.pi)


# --------------------------------------------------------------------------
# The two-loop weight on the cover of the z plane
# --------------------------------------------------------------------------
def _f2_from_letters(L, li_a, li_b):
    """Eq. (2.8) of arXiv:2112.06956 with log z = L and given Li_n(1/(1+z)), Li_n(z/(1+z))."""
    pi = mp.pi
    s2 = li_a[2] + li_b[2]
    return (-L**4 / 2
            + L**2 * (mp.mpf(2) / 3 * s2 - 19 * pi**2 / 9)
            + L * (4 * li_a[3] - 4 * li_b[3])
            + mp.mpf(2) / 3 * (s2 - pi**2 / 6) ** 2
            + mp.mpf(8) / 3 * pi**2 * (s2 - pi**2 / 6)
            + 8 * li_a[4] + 8 * li_b[4]
            - 23 * pi**4 / 18)


def f2_cover(L, route: str = "inner"):
    """The two-loop weight as a function of L = log z for |Im L| < 3 pi.

    On |Im L| <= pi this is the published function on the principal sheet.
    For -3pi < Im L < -pi it is the analytic continuation across the lower
    lip of the cut. Two routes are available: across -1 < z < 0 ("inner",
    where Li_n(1/(1+z)) crosses its cut) or across z < -1 ("outer", where
    Li_n(z/(1+z)) does). They agree if and only if the function is regular at
    z = -1 on the physical sheet. The upper sheets follow by conjugation.

    Exactly on the lines Im L = +-pi the polylogarithms are evaluated on
    their cuts and the value returned is that of one lip only; evaluate
    slightly off the line.
    """
    L = mp.mpmathify(L)
    im = mp.im(L)
    if im > mp.pi:
        return mp.conj(f2_cover(mp.conj(L), route))
    z = mp.exp(L)
    a, b = 1 / (z + 1), z / (z + 1)
    li_a = {n: mp.polylog(n, a) for n in (2, 3, 4)}
    li_b = {n: mp.polylog(n, b) for n in (2, 3, 4)}
    if im < -mp.pi:
        if im <= -3 * mp.pi:
            raise ValueError("continuation implemented for |Im L| < 3 pi only")
        for n in (2, 3, 4):
            if route == "inner":
                li_a[n] += 2j * mp.pi * mp.log(a) ** (n - 1) / mp.factorial(n - 1)
            else:
                li_b[n] -= 2j * mp.pi * mp.log(b) ** (n - 1) / mp.factorial(n - 1)
    return _f2_from_letters(L, li_a, li_b)


# --------------------------------------------------------------------------
# Item 1: the exact fourth propagator, with its angular dependence
# --------------------------------------------------------------------------
def exact_edge_kernel(cusp: Cusp, f, c_prev, c_next, rho_max, exact=True,
                      n_theta: int = 16, tilt=0.6):
    """Angular average of the transverse integral up to rho_max.

    With exact=True the far propagator is
        D'(1 + alpha - beta) + w + 2 sqrt(rho D' u / D) cos(theta),
    both in log z and in the measure; with exact=False it is held at its edge
    value D'(1 - beta)(1 + alpha). The contour is rho = r exp(i tilt),
    0 <= r <= rho_max, the same for both, so that their difference isolates
    the effect of the approximation. Returns the average of
        (far_edge / far) f(L) / [w (w + c_prev)(w + c_next)].
    """
    c_prev, c_next = mp.mpf(c_prev), mp.mpf(c_next)
    D, Dp = mp.mpf(cusp.D), mp.mpf(cusp.Dfar)
    u = -mp.mpf(cusp.s) - mp.mpf(cusp.t)
    alpha, beta = -c_next / D, c_prev / D
    far_edge = cusp.far_edge_value(c_prev, c_next)
    w0 = cusp.w_start(c_prev, c_next)
    phase = mp.exp(1j * mp.mpf(tilt))
    amp = 2 * mp.sqrt(Dp * u / D)
    marks = sorted({mp.mpf(0), abs(w0), abs(c_prev), abs(c_next), mp.mpf(rho_max)})
    marks = [m for m in marks if m <= rho_max]

    def one_angle(theta):
        cos_t = mp.cos(theta)

        def integrand(r):
            rho = r * phase
            w = w0 - rho
            base = cusp.log_z(w, c_prev, c_next, far_edge)
            if exact:
                far = Dp * (1 + alpha - beta) + w + amp * mp.sqrt(rho) * cos_t
                ratio = far / far_edge                     # close to 1: principal log is safe
                L = base + cusp.sigma * mp.log(ratio)
                weight = 1 / ratio
            else:
                L, weight = base, 1
            return weight * f(L) / (w * (w + c_prev) * (w + c_next)) * phase

        return mp.quad(integrand, marks)

    if not exact:
        return one_angle(mp.mpf(0))
    thetas = [2 * mp.pi * (k + mp.mpf("0.5")) / n_theta for k in range(n_theta)]
    return sum(one_angle(th) for th in thetas) / n_theta


def exact_edge_jump(cusp: Cusp, f, position, small, rho_max, exact=True, **kw):
    """Jump across the collinear invariant = 0 on the 'next' edge at `position`."""
    large = mp.mpf(position) * cusp.D
    small = abs(mp.mpf(small))
    plus = large * small * exact_edge_kernel(cusp, f, large, small, rho_max, exact, **kw)
    minus = large * (-small) * exact_edge_kernel(cusp, f, large, -small, rho_max, exact, **kw)
    return plus - minus


# --------------------------------------------------------------------------
# Item 2: the same density with the dimensional-regularisation weight
# --------------------------------------------------------------------------
def dimreg_quadrant_sum(cusp: Cusp, p, eps, c_prev, c_next):
    """Quadrant sum for F = z^p with the transverse measure rho^(-eps) d rho,
    rescaled by |c_prev c_next / D|^eps so that it has a limit at the cusp."""
    p, eps = mp.mpf(p), mp.mpf(eps)
    c_prev, c_next = abs(mp.mpf(c_prev)), abs(mp.mpf(c_next))
    direction = -mp.exp(1j * mp.mpf("0.6"))
    total = mp.mpc(0)
    for a in (1, -1):
        for b in (1, -1):
            cp, cn = a * c_prev, b * c_next
            far = cusp.far_edge_value(cp, cn)
            w0 = cusp.w_start(cp, cn)

            def integrand(r, cp=cp, cn=cn, far=far, w0=w0):
                w = w0 + r * direction
                rho = -r * direction                       # rho = w0 - w
                return (rho ** (-eps) * mp.exp(p * cusp.log_z(w, cp, cn, far))
                        / (w * (w + cp) * (w + cn)) * (-direction))

            marks = sorted({mp.mpf(0), abs(w0), c_prev, c_next, mp.mpf(1), mp.inf})
            total += c_prev * c_next * mp.quad(integrand, marks)   # sign(cp cn) * cp * cn
    return total * (c_prev * c_next / abs(cusp.D)) ** eps


def dimreg_density_closed_form(p, eps):
    """Cusp limit of `dimreg_quadrant_sum` for eps > |Re p|, in closed form.

    In the region rho ~ |c_prev c_next / D| the invariants reduce to
    z = -1 + rho |D| / (c_prev c_next), so each sign quadrant contributes a
    ray integral of z^(p-1) |1+z|^(-eps) on the sheet fixed by the -i0 rule.
    The rays are Euler Beta integrals. As eps -> 0 the expression tends to
    2 pi i sin(pi p)/(pi p), the value obtained without the regulator.
    """
    p, eps = mp.mpmathify(p), mp.mpmathify(eps)
    g = mp.gamma
    phase = mp.exp(-1j * mp.pi * p)
    rays = (2 * mp.cos(mp.pi * p) * g(p)
            * (phase * g(eps - p) / g(eps) - g(1 - eps) / g(1 + p - eps))
            - 2 * phase * g(1 - eps) * g(eps - p) / g(1 - p))
    return -rays


# --------------------------------------------------------------------------
# Analytic regularity of the two-loop weight at z = -1
# --------------------------------------------------------------------------
def f2_expansion_at_minus_one():
    """Expand the two-loop weight in lambda = log(1 + z) near z = -1.

    Write e = 1 + z with Im e > 0 (the upper lip), x = 1/e. The inversion
    relations

        Li_2(y) + Li_2(1/y) = -pi^2/6 - log^2(-y)/2
        Li_3(y) - Li_3(1/y) = -pi^2 log(-y)/6 - log^3(-y)/6
        Li_4(y) + Li_4(1/y) = -7 pi^4/360 - pi^2 log^2(-y)/12 - log^4(-y)/24

    applied to y = 1/(1+z) and y = z/(1+z) express the weight through
        a_n = Li_n(e),  b_n = Li_n(-e/(1-e)),  c = log(1-e)
    which are analytic at e = 0, and through lambda = log e, with
    log(-1/e) = i pi - lambda, log((1-e)/e) = c - lambda, log z = i pi + c.

    Returns (coefficients, landen), where coefficients[k] multiplies
    lambda^k and landen = 2 a_2 + 2 b_2 + c^2, which vanishes identically
    (Landen's identity). The coefficients of lambda^3 and lambda^4 are zero
    as polynomials; those of lambda and lambda^2 are multiples of `landen`.
    """
    import sympy as sp

    lam, c = sp.symbols("lambda c")
    a2, a3, a4, b2, b3, b4 = sp.symbols("a2 a3 a4 b2 b3 b4")
    pi, I = sp.pi, sp.I

    def inverted(n, li, log_minus):
        if n == 2:
            return -li - pi**2 / 6 - log_minus**2 / 2
        if n == 3:
            return li - pi**2 / 6 * log_minus - log_minus**3 / 6
        return -li - 7 * pi**4 / 360 - pi**2 / 12 * log_minus**2 - log_minus**4 / 24

    la, lb = I * pi - lam, c - lam
    A = {2: inverted(2, a2, la), 3: inverted(3, a3, la), 4: inverted(4, a4, la)}
    B = {2: inverted(2, b2, lb), 3: inverted(3, b3, lb), 4: inverted(4, b4, lb)}
    L = I * pi + c
    s2 = A[2] + B[2]
    F = (-L**4 / 2 + L**2 * (sp.Rational(2, 3) * s2 - 19 * pi**2 / 9)
         + L * (4 * A[3] - 4 * B[3])
         + sp.Rational(2, 3) * (s2 - pi**2 / 6) ** 2
         + sp.Rational(8, 3) * pi**2 * (s2 - pi**2 / 6)
         + 8 * A[4] + 8 * B[4] - 23 * pi**4 / 18)
    poly = sp.Poly(sp.expand(F), lam)
    coefficients = {k: sp.expand(poly.coeff_monomial(lam**k)) for k in range(5)}
    landen = 2 * a2 + 2 * b2 + c**2
    symbols = {"c": c, "a": (a2, a3, a4), "b": (b2, b3, b4)}
    return coefficients, landen, symbols


def landen_derivative():
    """d/de of Li_2(e) + Li_2(-e/(1-e)) + log^2(1-e)/2, which must vanish."""
    import sympy as sp

    e = sp.symbols("e")
    c = sp.log(1 - e)
    # d Li_2(u)/de = -log(1-u) u'/u ;  u = -e/(1-e) has 1 - u = 1/(1-e)
    u = -e / (1 - e)
    d_a2 = -sp.log(1 - e) / e
    d_b2 = c * sp.diff(u, e) / u
    return sp.simplify(d_a2 + d_b2 + c * sp.diff(c, e))


def f2_regular_part(e):
    """The lambda^0 coefficient evaluated at a complex e = 1 + z."""
    import sympy as sp

    coefficients, _, sym = f2_expansion_at_minus_one()
    e = mp.mpmathify(e)
    u = -e / (1 - e)
    values = {sym["c"]: mp.log(1 - e)}
    for n, (an, bn) in zip((2, 3, 4), zip(sym["a"], sym["b"])):
        values[an] = mp.polylog(n, e)
        values[bn] = mp.polylog(n, u)
    names = list(values)
    fn = sp.lambdify(names, coefficients[0], "mpmath")
    return fn(*[values[n] for n in names])


# --------------------------------------------------------------------------
# Uniform bound on the regulated integral of z^p (Lemma 4)
# --------------------------------------------------------------------------
def regulated_double_pole(p, eta, ratio=1):
    """eps^2 s t times the regulated integral of z^p at eps = -eta, complex p.

    The Mellin-Barnes contour of GATE_6_FUNCTIONAL, w = -1 - p + eta/2 + i y,
    makes the p-dependence explicit:

        eta^2 / Gamma(2 eta) / (2 pi)
          * int dy |Gamma(eta/2 + i y)|^4 ratio^(eta/2 + i y)
            Gamma(1 + p - eta/2 - i y) Gamma(1 - p - eta/2 + i y)
            / (Gamma(1+p)^2 Gamma(1-p)^2).

    Tends to 4 sin(pi p)/(pi p) as eta -> 0.
    """
    p, eta = mp.mpmathify(p), mp.mpf(eta)

    def integrand(y):
        weight = abs(mp.gamma(eta / 2 + 1j * y)) ** 4 * mp.mpf(ratio) ** (eta / 2 + 1j * y)
        return weight * mp.gamma(1 + p - eta / 2 - 1j * y) * mp.gamma(1 - p - eta / 2 + 1j * y)

    knots = [-mp.inf, -5, -1, -eta, 0, eta, 1, 5, mp.inf]
    integral = mp.quad(integrand, knots) / (2 * mp.pi)
    return eta**2 * mp.rgamma(2 * eta) * mp.rgamma(1 + p) ** 2 * mp.rgamma(1 - p) ** 2 * integral


# --------------------------------------------------------------------------
# Pinch surfaces (Landau analysis)
# --------------------------------------------------------------------------
def landau_pinch_analysis():
    """Where can the integration contour be pinched?

    All four propagators carry the same -i0, so a pinch at a real point needs
    positive weights lambda_i on the vanishing invariants with
    sum_i lambda_i (x0 - x_i) = 0: x0 is a convex combination of the x_i
    involved and is null-separated from each of them. With
    x0 = x2 - a p1 + (c + d) p2 + d p3 the four invariants are polynomials in
    (a, c, d); the function returns them and the solution sets.
    """
    import sympy as sp

    a, c, d, s, t = sp.symbols("a c d s t", real=True)

    def sq(k1, k2, k3):          # (k1 p1 + k2 p2 + k3 p3)^2 with 2p1p2=s, 2p2p3=t, 2p1p3=-s-t
        return sp.expand(k1 * k2 * s + k2 * k3 * t + k1 * k3 * (-s - t))

    inv = {1: sq(1 - a, c + d, d), 2: sq(-a, c + d, d),
           3: sq(-a, c + d - 1, d), 4: sq(-a, c + d - 1, d - 1)}
    return {
        "invariants": inv,
        "adjacent_pair_on_edge": [sp.simplify(inv[2].subs({a: 0, d: 0})),
                                  sp.simplify(inv[3].subs({a: 0, d: 0}))],
        "diagonal_pair": [sp.factor(inv[1].subs({d: 0, c: 1 - a})),
                          sp.factor(inv[3].subs({d: 0, c: 1 - a}))],
        "triple_123": sp.solve([inv[1].subs(d, 0), inv[2].subs(d, 0), inv[3].subs(d, 0)],
                               [a, c], dict=True),
        "quadruple": sp.solve([inv[1], inv[2], inv[3], inv[4]], [a, c, d], dict=True),
        "symbols": (a, c, d, s, t),
    }


# --------------------------------------------------------------------------
# Feynman-parameter form of the regulated integral
# --------------------------------------------------------------------------
def parametric_double_pole(p, eta, ratio=1):
    """eps^2 s t times the regulated integral of z^p from Feynman parameters.

    With propagator powers 1+p on legs 1, 3 and 1-p on legs 2, 4, in
    D = 4 + 2 eta dimensions,

        Gamma(2 - eta) / (Gamma(1+p)^2 Gamma(1-p)^2)
          * int d^4a delta(1 - a2 - a4) (a1 a3)^p (a2 a4)^(-p)
            (a1 + a2 + a3 + a4)^(-2 eta) / (a1 a3 s + a2 a4 t)^(2 - eta),

    using the Cheng-Wu freedom to put the delta function on a2 + a4. The
    extra-dimensional integral has already been done: it is
        int_0^oo d mu mu^(eta-1)/Gamma(eta) * Gamma(2)/(A + mu)^2
          = Gamma(2 - eta) A^(eta - 2).
    Substituting a1 = sqrt(q) e^y, a3 = sqrt(q) e^(-y) leaves three nested
    one-dimensional integrals. Agrees with `regulated_double_pole`, which is
    the Mellin-Barnes form of the same quantity.
    """
    from scipy.integrate import quad
    from scipy.interpolate import CubicSpline
    from scipy.special import gamma as gamma_fn

    p, eta, ratio = float(p), float(eta), float(ratio)

    def h(q):
        root = 2 * np.sqrt(q)
        with np.errstate(over="ignore"):
            value = quad(lambda y: (1 + root * np.cosh(y)) ** (-2 * eta), 0, np.inf,
                         epsabs=0, epsrel=1e-11, limit=400)[0]
        return 2 * value

    grid = np.linspace(-34, 34, 681)                       # log q
    spline = CubicSpline(grid, np.log([h(np.exp(g)) for g in grid]))

    def over_q(X):
        # q = X exp(v): the integrand in v decays like exp(-(1-p-eta)|v|) or faster
        def integrand(v):
            q = X * np.exp(v)
            return q ** (1 + p) * np.exp(spline(np.log(q))) / (q * ratio + X) ** (2 - eta)
        return quad(integrand, -30 - np.log(X), 30 - np.log(X), epsabs=0, epsrel=1e-10,
                    limit=400, points=[0.0])[0]

    total = 2 * quad(lambda x: (x * (1 - x)) ** (-p) * over_q(x * (1 - x)), 0, 0.5,
                     epsabs=0, epsrel=1e-9, limit=200)[0]
    # s t from the prefactor and (s/t)^p from the definition of z, with t = 1
    return (eta**2 * ratio ** (1 + p) * gamma_fn(2 - eta)
            / (gamma_fn(1 + p) ** 2 * gamma_fn(1 - p) ** 2) * total)


# --------------------------------------------------------------------------
# Complex cross-ratio: the timelike diagonal
# --------------------------------------------------------------------------
def regulated_double_pole_complex(p, eta, modulus=1, theta=0):
    """The Mellin-Barnes form at S/T = modulus * exp(i theta), 0 <= theta <= pi.

    The integrand is |Gamma(eta/2 + i y)|^4 (S/T)^(eta/2 + i y) times two
    Gamma functions. The Gamma factors decay like exp(-3 pi |y|) in total,
    and |(S/T)^(i y)| = exp(-theta y), so the integral converges absolutely
    for |theta| < 3 pi and is continuous in theta up to theta = pi, which is
    the timelike diagonal s -> s - i0. Returns eps^2 S T (S/T)^p times the
    regulated integral with T = 1.
    """
    p, eta = mp.mpmathify(p), mp.mpf(eta)
    log_ratio = mp.log(mp.mpf(modulus)) + 1j * mp.mpf(theta)

    def integrand(y):
        weight = abs(mp.gamma(eta / 2 + 1j * y)) ** 4 * mp.exp((eta / 2 + 1j * y) * log_ratio)
        return weight * mp.gamma(1 + p - eta / 2 - 1j * y) * mp.gamma(1 - p - eta / 2 + 1j * y)

    knots = [-mp.inf, -5, -1, -eta, 0, eta, 1, 5, mp.inf]
    integral = mp.quad(integrand, knots) / (2 * mp.pi)
    return eta**2 * mp.rgamma(2 * eta) * mp.rgamma(1 + p) ** 2 * mp.rgamma(1 - p) ** 2 * integral


def parametric_double_pole_rotated(p, eta, modulus=1, theta=np.pi, gamma=0.5):
    """The Feynman-parameter form with S = modulus exp(i theta), T = 1, on a rotated contour.

    In the Cheng-Wu gauge a2 + a4 = 1 with q = a1 a3,
        F = q S + X,  X = a2 a4.
    For theta = pi (timelike diagonal) F vanishes at the real point q = X/|S|.
    Rotating q -> r exp(-i gamma) puts every term of F in the cone
    0 <= arg <= pi - gamma, so |F| >= sin(gamma/2)(X + r |S|) and the integral
    converges absolutely; by Cauchy's theorem it equals the boundary value
    from theta < pi. The y integral defines h(q) for complex q.
    """
    from scipy.integrate import quad
    from scipy.interpolate import CubicSpline
    from scipy.special import gamma as gamma_fn

    p, eta = float(p), float(eta)
    S = modulus * np.exp(1j * theta)
    rot = np.exp(-1j * gamma)

    def h(r):
        root = 2 * np.sqrt(r) * np.exp(-0.5j * gamma)
        with np.errstate(over="ignore", invalid="ignore"):
            f = lambda y: (1 + root * np.cosh(y)) ** (-2 * eta)
            re = quad(lambda y: f(y).real, 0, 60, epsabs=0, epsrel=1e-11, limit=600)[0]
            im = quad(lambda y: f(y).imag, 0, 60, epsabs=0, epsrel=1e-11, limit=600)[0]
        return 2 * (re + 1j * im)

    grid = np.linspace(-34, 34, 681)                         # log r
    values = np.array([h(np.exp(g)) for g in grid])
    spline_re, spline_im = CubicSpline(grid, values.real), CubicSpline(grid, values.imag)

    def over_q(X):
        def integrand(v):
            r = X * np.exp(v)
            q = r * rot
            hq = spline_re(np.log(r)) + 1j * spline_im(np.log(r))
            return r * rot * q**p * hq * (q * S + X) ** (eta - 2)     # dq = r rot dv
        lo, hi = -30 - np.log(X), 30 - np.log(X)
        re = quad(lambda v: integrand(v).real, lo, hi, epsabs=0, epsrel=1e-10, limit=600, points=[0.0])[0]
        im = quad(lambda v: integrand(v).imag, lo, hi, epsabs=0, epsrel=1e-10, limit=600, points=[0.0])[0]
        return re + 1j * im

    def outer(x):
        return (x * (1 - x)) ** (-p) * over_q(x * (1 - x))

    re = quad(lambda x: outer(x).real, 0, 0.5, epsabs=0, epsrel=1e-9, limit=200)[0]
    im = quad(lambda x: outer(x).imag, 0, 0.5, epsabs=0, epsrel=1e-9, limit=200)[0]
    total = 2 * (re + 1j * im)
    return (eta**2 * S ** (1 + p) * gamma_fn(2 - eta)
            / (gamma_fn(1 + p) ** 2 * gamma_fn(1 - p) ** 2) * total)


def box_continued(p, eta, modulus=1, theta=0):
    """N_eta(p) continued in eta to -1/2 < Re eta < 1/2 - |Re p|.

    The Mellin-Barnes line Re zeta = eta/2 of `regulated_double_pole` is moved
    to Re zeta = -1/2 across the double pole of Gamma(zeta)^2 at zeta = 0:

        N = eta^2 / (Gamma(2 eta) Gamma(1+p)^2 Gamma(1-p)^2) * (R + J),
        R = Gamma(eta)^2 Gamma(1+p) Gamma(1-p-eta)
            * (log(S/T) - 2 psi(eta) - 2 gamma - psi(1+p) + psi(1-p-eta)),
        J = (1/2 pi) int dy  G(-1/2 + i y),
        G(zeta) = Gamma(zeta)^2 Gamma(eta-zeta)^2 (S/T)^zeta
                  Gamma(1+p-zeta) Gamma(1-p-eta+zeta).

    J is analytic for Re eta > -1/2, and the zero of eta^2/Gamma(2 eta) cancels
    the triple pole of R at eta = 0, so N is analytic in the whole region,
    including eta = 0, where it equals 4 sin(pi p)/(pi p). S/T = modulus *
    exp(i theta), 0 <= theta <= pi; theta = pi is a timelike diagonal.
    """
    p, eta = mp.mpmathify(p), mp.mpmathify(eta)
    log_ratio = mp.log(mp.mpf(modulus)) + 1j * mp.mpf(theta)
    if not -0.5 < mp.re(eta) < 0.5 - abs(mp.re(p)):
        raise ValueError("need -1/2 < Re eta < 1/2 - |Re p|")

    def G(y):
        zeta = mp.mpf(-0.5) + 1j * y
        return (mp.gamma(zeta) ** 2 * mp.gamma(eta - zeta) ** 2 * mp.exp(zeta * log_ratio)
                * mp.gamma(1 + p - zeta) * mp.gamma(1 - p - eta + zeta))

    centre = mp.im(p)
    knots = sorted({-mp.inf, -6, -1, 0, 1, 6, mp.inf, centre})
    J = mp.quad(G, knots) / (2 * mp.pi)
    front = mp.rgamma(1 + p) ** 2 * mp.rgamma(1 - p) ** 2
    if abs(eta) < mp.mpf(10) ** (-mp.mp.dps // 2):
        return 4 * mp.gamma(1 + p) * mp.gamma(1 - p) * front
    R = (mp.gamma(eta) ** 2 * mp.gamma(1 + p) * mp.gamma(1 - p - eta)
         * (log_ratio - 2 * mp.psi(0, eta) - 2 * mp.euler - mp.psi(0, 1 + p) + mp.psi(0, 1 - p - eta)))
    return eta ** 2 * mp.rgamma(2 * eta) * front * (R + J)
