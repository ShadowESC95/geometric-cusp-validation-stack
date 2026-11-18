"""Lorentzian integration chain at a cusp, in position space.

This module supports docs/GATE_6_CLOSURE.md. It works with the frozen-loop
integral

    x13^2 x24^2  int d^4x0  F(z) / (x01^2 x02^2 x03^2 x04^2),
    z = x02^2 x04^2 x13^2 / (x01^2 x03^2 x24^2),

every x0i^2 carrying the Feynman -i0, near the cusp x0 -> x2 of a light-like
square with x13^2 = s < 0 and x24^2 = t < 0. Metric (+,-,-,-).

Write x0 = x2 + alpha p1 + beta p2 + yperp, with p1 = x2 - x1 and
p2 = x3 - x2 the two null edges at the cusp. Then, exactly,

    x02^2 = w,   x01^2 = w + q1,   x03^2 = w + q3,
    w = q1 q3 / |s| - rho,   q1 = beta s,   q3 = -alpha s,   rho = |yperp|^2,

and x04^2 -> t at the cusp. With the fourth propagator frozen at t, the
transverse integral at fixed (q1, q3) is the kernel

    K_f(q1, q3) = int_0^oo d rho  f(L(w)) / [ w (w + q1) (w + q3) ],

where L = log z on the sheet fixed by the four -i0 prescriptions and
F(z) = f(log z). `quadrant_sum` forms the combination of the four sign
quadrants of (q1, q3) that multiplies the double logarithm.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import sympy as sp

# --------------------------------------------------------------------------
# Kinematics: explicit vectors, to check the decomposition above
# --------------------------------------------------------------------------
METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def msq(v):
    """Minkowski square, signature (+,-,-,-)."""
    v = np.asarray(v, dtype=float)
    return float(v @ METRIC @ v)


def light_like_square(angle: float = 1.1):
    """Null edges p1..p4 with p1+p2+p3+p4 = 0 and both diagonals spacelike.

    Returns (p, s, t) with s = (p1+p2)^2 < 0 and t = (p2+p3)^2 < 0.
    """
    c, sn = np.cos(angle), np.sin(angle)
    p = np.array([[1.0, 0.0, 0.0, 1.0],
                  [-1.0, -sn, 0.0, -c],
                  [1.0, 0.0, 0.0, -1.0],
                  [-1.0, sn, 0.0, c]])
    return p, msq(p[0] + p[1]), msq(p[1] + p[2])


def transverse_basis(p1, p2):
    """Two unit spacelike vectors orthogonal to both null edges."""
    basis = []
    for seed in np.eye(4):
        v = seed.copy()
        # remove the components along p1 and p2 (they are null, p1.p2 != 0)
        dot12 = p1 @ METRIC @ p2
        v = v - (v @ METRIC @ p2) / dot12 * p1 - (v @ METRIC @ p1) / dot12 * p2
        for b in basis:
            v = v + (v @ METRIC @ b) * b          # b.b = -1
        norm = -msq(v)
        if norm > 1e-12:
            basis.append(v / np.sqrt(norm))
        if len(basis) == 2:
            break
    return basis


def cusp_invariants(alpha, beta, yperp, angle: float = 1.1):
    """The four x0i^2 for x0 = x2 + alpha p1 + beta p2 + yperp, and (s, t)."""
    p, s, t = light_like_square(angle)
    e1, e2 = transverse_basis(p[0], p[1])
    y = alpha * p[0] + beta * p[1] + yperp[0] * e1 + yperp[1] * e2
    return {"x01": msq(y + p[0]), "x02": msq(y), "x03": msq(y - p[1]),
            "x04": msq(y - p[1] - p[2]), "s": s, "t": t,
            "rho": yperp[0] ** 2 + yperp[1] ** 2}


# --------------------------------------------------------------------------
# Exact local identities (sympy)
# --------------------------------------------------------------------------
def local_identities():
    """Residuals of the algebraic identities used in the proof; all are zero.

    measure_split : 1/[w (w+q1)(w+q3)] dw
                    = (1/(q1 q3)) dlog z_edge + (1/(q1 (q1-q3))) [dw/(w+q1) - dw/(w+q3)]
    z_factor      : z = z_edge / (1 + w/q1)
    start         : z_edge at rho = 0 equals -1 / (1 + q1/|s|)
    cap           : 1 + z_edge = -rho |s| ... / (q1 (w + q3)), so z_edge = -1
                    only on the plane yperp = 0 up to the stated O(q1/|s|) shift
    """
    w, q1, q3, S, rho = sp.symbols("w q1 q3 S rho", real=True, nonzero=True)
    z = -S * w / ((w + q1) * (w + q3))
    z_edge = -(S / q1) * w / (w + q3)
    dlog_edge = sp.diff(sp.log(z_edge), w)
    measure = 1 / (w * (w + q1) * (w + q3))
    split = dlog_edge / (q1 * q3) + (1 / (w + q1) - 1 / (w + q3)) / (q1 * (q1 - q3))
    w0 = q1 * q3 / S
    return {
        "measure_split": sp.simplify(measure - split),
        "z_factor": sp.simplify(z - z_edge / (1 + w / q1)),
        "start": sp.simplify(z_edge.subs(w, w0) + 1 / (1 + q1 / S)),
        "cap": sp.simplify(1 + z_edge.subs(w, w0 - rho)
                           - (q1 * (w0 - rho + q3) - S * (w0 - rho)) / (q1 * (w0 - rho + q3))),
    }


def weighted_chart_dictionary():
    """Brackets in the U_12 chart A = e1 + a e3 + b e4, B = e2 + c e3 + d e4.

    Returns the four propagator brackets and 1 + z. The strict cusp line
    b = c = d = 0 is the null edge through x2 and x3; `a` moves along it.
    """
    a, b, c, d = sp.symbols("a b c d")
    line = sp.Matrix([[1, 0, a, b], [0, 1, c, d]])
    unit = sp.eye(4)
    br = lambda i, j: sp.Matrix.vstack(line, unit[i - 1, :], unit[j - 1, :]).det()
    out = {"12": sp.expand(br(1, 2)), "23": sp.expand(br(2, 3)),
           "34": sp.expand(br(3, 4)), "14": sp.expand(br(1, 4))}
    z = out["12"] * out["34"] / (out["23"] * out["14"])
    out["one_plus_z"] = sp.factor(sp.simplify(1 + z))
    return out


# --------------------------------------------------------------------------
# Sheet bookkeeping and the transverse kernel
# --------------------------------------------------------------------------
def edge_value_x04(q1, q3, s=-1.0, t=-1.0):
    """x04^2 on the plane yperp = 0, exactly: t (1 + alpha - beta) + alpha beta s.

    With beta = -q1/|s| and alpha = q3/|s| (s < 0). On the null edge
    (alpha = 0) it is t (1 - beta); it changes sign at the far cusp beta = 1.
    """
    return t * (1 + q1 / abs(s) + q3 / abs(s)) + q1 * q3 / abs(s)


def log_z(w, q1, q3, s=-1.0, t=-1.0, x04=None):
    """log z with the sheet fixed by the Feynman prescription.

    Each x0i^2 - i0 is continued from the spacelike region through the lower
    half plane, so its logarithm has argument -pi when spacelike and 0 when
    timelike. `w` must have a (small) negative imaginary part. x04^2 is a
    real constant: `x04`, or its edge value when omitted.
    """
    if not (s < 0 and t < 0):
        raise ValueError("this module assumes spacelike diagonals: s < 0 and t < 0")
    if x04 is None:
        x04 = edge_value_x04(q1, q3, s, t)
    log4 = mp.log(abs(x04)) - (1j * mp.pi if x04 < 0 else 0)
    # z = (s/t) x02^2 x04^2 / (x01^2 x03^2), with s/t > 0
    return mp.log(w) + log4 - mp.log(w + q1) - mp.log(w + q3) + mp.log(s / t)


def sheet_path(q1, q3, s=-1.0, n: int = 400, eta: float = 1e-9, reach: float = 0.02):
    # uses t = s, so that s/t = 1
    """Sample (log|z|, arg z) along the Lorentzian contour in the edge regime.

    rho runs from 0 until |w| = reach * |q1|, so that x01^2 keeps its sign.
    eta is the size of the -i0 displacement relative to the local scale.
    """
    q1, q3 = mp.mpf(q1), mp.mpf(q3)
    w0 = q1 * q3 / abs(s)
    scale = min(abs(w0), abs(q3))
    rhos = [mp.mpf(0)] + [mp.mpf(10) ** x for x in
                          np.linspace(float(mp.log10(scale)) - 3,
                                      float(mp.log10(reach * abs(q1))), n)]
    ell, phi = [], []
    for r in rhos:
        w = w0 - r
        L = log_z(w - 1j * eta * max(abs(w), abs(w + q3), scale * mp.mpf("1e-3")), q1, q3, s, s)
        ell.append(float(mp.re(L)))
        phi.append(float(mp.im(L)))
    return np.array(ell), np.array(phi)


def plateaus(phi, tol: float = 0.02):
    """The sequence of distinct values of arg z / pi visited along a path."""
    out = []
    for value in np.asarray(phi) / np.pi:
        k = round(value)
        if abs(value - k) < tol and (not out or out[-1] != k):
            out.append(int(k))
    return out


def kernel(f, q1, q3, s=-1.0, t=-1.0, tilt: float = 0.6, x04=None):
    """K_f(q1, q3) = int_0^oo d rho f(L) / [w (w+q1)(w+q3)], w = q1 q3/|s| - rho - i0.

    The constant factor 1/x04^2 of the measure is left out; it cancels in
    every ratio to the constant weight.

    The integrand is holomorphic for Im w < 0, so the contour is rotated to
    the ray w = w0 - r exp(i tilt), which keeps clear of the three poles.
    """
    q1, q3 = mp.mpf(q1), mp.mpf(q3)
    w0 = q1 * q3 / abs(s)
    direction = -mp.exp(1j * mp.mpf(tilt))

    def integrand(r):
        w = w0 + r * direction
        return f(log_z(w, q1, q3, s, t, x04)) / (w * (w + q1) * (w + q3)) * (-direction)

    marks = sorted({mp.mpf(0), abs(w0), abs(q3), abs(q1), mp.mpf(abs(s)), mp.inf})
    return mp.quad(integrand, marks)


def quadrant_sum(f, q1, q3, s=-1.0, t=-1.0):
    """Sum over the sign quadrants of sign(q1 q3) * q1 q3 * K_f.

    This is the coefficient of d|q1|/|q1| d|q3|/|q3| after the principal-value
    pairing of the two collinear variables, i.e. the density of the double
    logarithm. q1 and q3 are taken positive here.
    """
    q1, q3 = abs(mp.mpf(q1)), abs(mp.mpf(q3))
    total = mp.mpc(0)
    for s1 in (1, -1):
        for s3 in (1, -1):
            total += q1 * q3 * kernel(f, s1 * q1, s3 * q3, s, t)   # sign(q1 q3)*(s1 q1)(s3 q3)
    return total


def edge_jump(f, q1, q3, s=-1.0, t=-1.0):
    """q1 q3 K at q3 > 0 minus the same at q3 < 0, for fixed q1 (either sign).

    q1 = beta s need not be small: this is the density of the collinear
    logarithm at the point x2 + beta p2 of the null line through x2 and x3.
    The fourth propagator is held at its exact edge value t (1 - beta).
    """
    q3 = abs(mp.mpf(q3))
    q1 = mp.mpf(q1)
    x04 = t * (1 + q1 / abs(s))
    return (q1 * q3 * kernel(f, q1, q3, s, t, x04=x04)
            - q1 * (-q3) * kernel(f, q1, -q3, s, t, x04=x04))


def circle_mean(f, radius=1):
    """(1/2pi) int_{-pi}^{pi} f(log r + i phi) d phi: mean of F = f(log z) on |z| = r."""
    ell = mp.log(radius)
    return mp.quad(lambda p: f(ell + 1j * p), mp.linspace(-mp.pi, mp.pi, 5)) / (2 * mp.pi)


# --------------------------------------------------------------------------
# Momentum twistors <-> position space (complexified, signature-free)
# --------------------------------------------------------------------------
def point_of_line(za, zb):
    """2x2 matrix X with mu = X lambda for both twistors (rows: lambda, mu)."""
    lam = np.array([za[:2], zb[:2]], dtype=complex).T
    mu = np.array([za[2:], zb[2:]], dtype=complex).T
    return mu @ np.linalg.inv(lam)


def chart_point(a, b, c, d, frame):
    """x0 and the four cusp points for the U_12 chart, in a generic frame.

    `frame` is an invertible 4x4 matrix applied to all twistors so that every
    point is finite. Returns (X0, [X1, X2, X3, X4]) with x_i = (Z_{i-1} Z_i).
    """
    ext = np.eye(4, dtype=complex) @ frame
    line = np.array([[1, 0, a, b], [0, 1, c, d]], dtype=complex) @ frame
    x0 = point_of_line(line[0], line[1])
    xs = [point_of_line(ext[(i - 1) % 4], ext[i]) for i in range(4)]
    return x0, xs


def position_cross_ratio(x0, xs):
    d = lambda p, q: np.linalg.det(p - q)
    x1, x2, x3, x4 = xs
    return d(x0, x2) * d(x0, x4) * d(x1, x3) / (d(x0, x1) * d(x0, x3) * d(x2, x4))


# --------------------------------------------------------------------------
# Reality conditions of the chart for a real Lorentzian square
# --------------------------------------------------------------------------
_SIGMA = (np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]]))


def hermitian_point(x):
    """The 2x2 Hermitian matrix of a real point; its determinant is x^2 in the metric (+,-,-,-)."""
    return sum(x[m] * _SIGMA[m] for m in range(4))


def real_square_twistors(region: str = "--", angle: float = 1.1):
    """Momentum twistors Z_1..Z_4 (rows) of a real light-like square, with
    x_i = (Z_{i-1} Z_i), each normalised to a unit spinor lambda_i."""
    from .lorentzian_edges import light_like_square as square

    p, s, t = square(region, angle)
    x = [np.zeros(4)]
    for k in range(3):
        x.append(x[-1] + p[k])
    Z = []
    for i in range(4):                               # (x_{i+1} - x_i) lambda_i = 0
        w, v = np.linalg.eig(hermitian_point(p[i]))
        lam = v[:, np.argmin(abs(w))]
        Z.append(np.concatenate([lam, hermitian_point(x[i]) @ lam]))
    return np.array(Z), np.array(x)


def twistor_form(Z):
    """H_ij = h(Z_i, Z_j) with h(Z, Z') = lambda^dagger mu' - mu^dagger lambda'.
    A line is a real point exactly when h vanishes on it."""
    lam, mu = Z[:, :2], Z[:, 2:]
    return lam.conj() @ mu.T - mu.conj() @ lam.T


def chart_of_real_point(x0, Z):
    """(a, b, c, d) with A = Z_1 + a Z_3 + b Z_4, B = Z_2 + c Z_3 + d Z_4 spanning the line of x0."""
    X0 = hermitian_point(x0)
    on_line = [np.concatenate([e, X0 @ e]) for e in np.eye(2)]
    coeff = np.array([np.linalg.solve(Z.T, v) for v in on_line])
    red = np.linalg.solve(coeff[:, :2], coeff)
    return red[0, 2], red[0, 3], red[1, 2], red[1, 3]
