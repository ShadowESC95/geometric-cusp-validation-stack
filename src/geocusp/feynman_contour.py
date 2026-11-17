"""An explicit integration contour for the frozen-loop integral.

See docs/GATE_6_CONTOUR.md. The real integration domain is displaced by a
purely imaginary vector field,

    x0  ->  x0 + i V(x0),        V real,

so that  x0i^2 -> x0i^2 - V^2 + 2 i V.(x0 - x_i)  exactly. The imaginary part
is linear in V, so any convex combination of admissible fields is admissible.

Building blocks (spacelike diagonals, s < 0 and t < 0):

* Frame fields. At the cusp x_j write x0 = x_j + alpha p_{j-1} + beta p_j
  + yperp and let T be the triangle {beta >= 0, alpha <= 0, beta - alpha <= 1}
  in the (beta, alpha) plane, the projection of conv(x_{j-1}, x_j, x_{j+1}).
  With d = P - proj_T(P), P = (beta, alpha),
      V_j = kappa (d_beta p_{j-1} + d_alpha p_j) + nu yperp
  satisfies, for i in {j-1, j, j+1},
      V_j.(x0 - x_i) = (kappa D/2) d.(P - vertex_i) - nu rho
                     <= -(kappa |D|/2) |d|^2 - nu rho <= 0,
  by the obtuse-angle property of the projection onto a convex set.

* Far fields. V = -kappa' (x0 - c) with c in the common past of the four
  cusps is admissible near every future nappe, and with c in the common
  future near every past nappe.

The frame field V_j is switched off near the cone of the opposite cusp
x_{j+2}, the only invariant it does not control.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .lorentzian_edges import METRIC, light_like_square


def mdot(a, b):
    return a @ METRIC @ b


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


def project_to_triangle(P):
    """Nearest point of T = {beta >= 0, alpha <= 0, beta - alpha <= 1} to P = (beta, alpha)."""
    verts = np.array([[0.0, 0.0], [0.0, -1.0], [1.0, 0.0]])
    b, a = P
    if b >= 0 and a <= 0 and b - a <= 1:
        return np.array(P, dtype=float)
    best, best_d = None, np.inf
    for k in range(3):
        u, v = verts[k], verts[(k + 1) % 3]
        tpar = np.clip(((P - u) @ (v - u)) / ((v - u) @ (v - u)), 0.0, 1.0)
        q = u + tpar * (v - u)
        dist = np.sum((P - q) ** 2)
        if dist < best_d:
            best, best_d = q, dist
    return best


@dataclass
class Contour:
    kappa: float = 3e-6          # longitudinal displacement of the frame fields
    nu: float = 3e-6             # transverse displacement of the frame fields
    kappa_far: float = 2e-3      # far field
    delta0: float = 0.05         # width of the neighbourhood of a cone, inside the ball
    radius: float = 8.0          # frame fields inside, far fields outside 2*radius
    tau0: float = 8.0            # time offset of the far-field centres
    angle: float = 1.1
    p: np.ndarray = field(init=False)
    x: np.ndarray = field(init=False)

    def __post_init__(self):
        self.p, self.s, self.t = light_like_square("--", self.angle)
        pts = [np.zeros(4)]
        for k in range(3):
            pts.append(pts[-1] + self.p[k])
        self.x = np.array(pts)                         # x[0..3] = x1..x4
        self.centre = self.x.mean(axis=0)
        e_t = np.array([1.0, 0, 0, 0])
        self.c_past = self.centre - self.tau0 * e_t
        self.c_future = self.centre + self.tau0 * e_t

    # ---- invariants -------------------------------------------------------
    def invariants(self, x0):
        """x0i^2 for (possibly complex) x0; bilinear, no conjugation."""
        return np.array([mdot(x0 - xi, x0 - xi) for xi in self.x])

    # ---- frame data -------------------------------------------------------
    def frame(self, j, x0):
        pm, pp = self.p[(j - 1) % 4], self.p[j]
        D = 2 * mdot(pm, pp)
        y = x0 - self.x[j]
        alpha, beta = 2 * mdot(y, pp) / D, 2 * mdot(y, pm) / D
        yperp = y - alpha * pm - beta * pp
        return pm, pp, D, alpha, beta, yperp, -mdot(yperp, yperp)

    def frame_field(self, j, x0):
        pm, pp, D, alpha, beta, yperp, rho = self.frame(j, x0)
        P = np.array([beta, alpha])
        d = P - project_to_triangle(P)
        return self.kappa * (d[0] * pm + d[1] * pp) + self.nu * yperp

    def frame_bound(self, j, x0):
        """-(kappa |D|/2) |d|^2 - nu rho: upper bound on V_j.(x0 - x_i), i = j-1, j, j+1."""
        pm, pp, D, alpha, beta, yperp, rho = self.frame(j, x0)
        P = np.array([beta, alpha])
        d = P - project_to_triangle(P)
        return -(self.kappa * abs(D) / 2) * (d @ d) - self.nu * rho

    # ---- the field --------------------------------------------------------
    def weights(self, x0):
        X = self.invariants(x0)
        return np.array([smoothstep(abs(X[(j + 2) % 4]) / self.delta0 - 1.0) for j in range(4)])

    def inner_field(self, x0):
        m = self.weights(x0)
        return sum(m[j] * self.frame_field(j, x0) for j in range(4))

    def outer_field(self, x0):
        r = x0 - self.centre
        space = np.linalg.norm(r[1:])
        theta = r[0] / max(space, 1e-300)
        up = smoothstep((theta - 0.25) / 0.25)
        down = smoothstep((-theta - 0.25) / 0.25)
        return (-self.kappa_far * up * (x0 - self.c_past)
                - self.kappa_far * down * (x0 - self.c_future))

    def blend(self, x0):
        return smoothstep(np.linalg.norm(x0 - self.centre) / self.radius - 1.0)

    def field(self, x0):
        g = self.blend(x0)
        v = np.zeros(4)
        if g < 1:
            v = v + (1 - g) * self.inner_field(x0)
        if g > 0:
            v = v + g * self.outer_field(x0)
        return v

    def displaced_invariants(self, x0):
        """The four invariants on the contour, and the field used."""
        v = self.field(x0)
        return self.invariants(x0 + 1j * v), v

    # ---- diagnostics ------------------------------------------------------
    @staticmethod
    def branch_argument(value):
        """Argument on the branch continued from the Feynman prescription.

        The physical values are 0 (timelike) and -pi (spacelike); the
        returned argument lies in (-3pi/2, pi/2].
        """
        a = np.angle(value)
        return a - 2 * np.pi if a > np.pi / 2 else a

    def excursion(self, value):
        """How far the argument leaves the closed interval [-pi, 0]."""
        a = self.branch_argument(value)
        return max(a, -np.pi - a, 0.0)

    def distance_to_edges(self, x0):
        """Euclidean distance from x0 to the union of the four closed edge segments."""
        best = np.inf
        for j in range(4):
            u, v = self.x[j], self.x[(j + 1) % 4]
            tpar = np.clip((x0 - u) @ (v - u) / ((v - u) @ (v - u)), 0.0, 1.0)
            best = min(best, float(np.linalg.norm(x0 - (u + tpar * (v - u)))))
        return best


# --------------------------------------------------------------------------
# Sampling diagnostics
# --------------------------------------------------------------------------
_G = np.diag([1.0, -1.0, -1.0, -1.0])


def sample_points(contour: Contour, n: int, seed: int = 0):
    """Real points chosen to stress the contour: generic points at three
    scales, points close to one cone, points close to two cones, and points
    close to the four null lines on and off the edge segments."""
    rng = np.random.default_rng(seed)
    C = contour

    def null_direction():
        v = rng.normal(size=3)
        return np.array([rng.choice([-1.0, 1.0]), *(v / np.linalg.norm(v))])

    for _ in range(n):
        kind = rng.integers(0, 6)
        if kind == 0:
            yield C.centre + rng.normal(size=4) * 2
        elif kind == 1:
            yield C.centre + rng.normal(size=4) * 10
        elif kind == 2:
            yield C.centre + rng.normal(size=4) * 60
        elif kind == 3:
            i = rng.integers(0, 4)
            yield (C.x[i] + 10 ** rng.uniform(-2, 2) * null_direction()
                   + rng.normal(size=4) * 10 ** rng.uniform(-4, -1))
        elif kind == 4:
            i = rng.integers(0, 4)
            k = (i + 2) % 4 if rng.random() < 0.7 else rng.integers(0, 4)
            x = C.x[i] + 10 ** rng.uniform(-1, 1.5) * null_direction()
            for _ in range(40):                      # Newton step onto both cones
                F = np.array([mdot(x - C.x[k], x - C.x[k]), mdot(x - C.x[i], x - C.x[i])])
                J = np.array([2 * _G @ (x - C.x[k]), 2 * _G @ (x - C.x[i])])
                x = x - np.linalg.lstsq(J, F, rcond=None)[0]
            yield x + rng.normal(size=4) * 10 ** rng.uniform(-6, -2)
        else:
            j = rng.integers(0, 4)
            yield (C.x[j] + rng.uniform(-1.5, 2.5) * C.p[j]
                   + rng.normal(size=4) * 10 ** rng.uniform(-4, -0.5))


def survey(contour: Contour, n: int, seed: int = 0, edge_exclusion: float = 1e-5):
    """Global properties of the contour on a sample.

    Returns
      max_excursion   largest amount by which arg(x0i^2) leaves [-pi, 0]
      max_im_near     largest Im(x0i^2) at points close to the cone of x_i
      near_count      number of (point, cone) pairs entering max_im_near
      min_ratio       smallest min_i |x0i^2| / min(1, dist_to_edges^2)
    """
    C = contour
    out = {"max_excursion": 0.0, "max_im_near": -np.inf, "near_count": 0, "min_ratio": np.inf}
    for x0 in sample_points(C, n, seed):
        dist = C.distance_to_edges(x0)
        if dist < edge_exclusion:
            continue
        real = C.invariants(x0)
        displaced, _ = C.displaced_invariants(x0)
        r = np.linalg.norm(x0 - C.centre)
        tol = C.delta0 / 2 if r < C.radius else 0.1 * r
        for i in range(4):
            out["max_excursion"] = max(out["max_excursion"], C.excursion(displaced[i]))
            if abs(real[i]) <= tol:
                out["near_count"] += 1
                out["max_im_near"] = max(out["max_im_near"], displaced[i].imag)
        out["min_ratio"] = min(out["min_ratio"], np.min(np.abs(displaced)) / min(1.0, dist**2))
    return out


def edge_lower_bounds(contour: Contour, n: int, seed: int = 0, cusp: int = 1):
    """Near a cusp, compare the displaced invariants with the scales that make
    the regulated integrand absolutely integrable.

    With x0 = x_j + alpha p_{j-1} + beta p_j + yperp the test quantities are
        |x0,j^2|   / max(rho, |alpha beta D|)
        |x0,j-1^2| / max(rho, |beta D| |1 + alpha|)
        |x0,j+1^2| / max(rho, |alpha D| |1 - beta|)
    Returns their minimum over the sample divided by nu.
    """
    rng = np.random.default_rng(seed)
    C = contour
    j = cusp
    pm, pp = C.p[(j - 1) % 4], C.p[j]
    D = 2 * mdot(pm, pp)
    # transverse unit vectors
    basis = []
    for seed_vec in np.eye(4):
        v = seed_vec - 2 * mdot(seed_vec, pp) / D * pm - 2 * mdot(seed_vec, pm) / D * pp
        for e in basis:
            v = v + mdot(v, e) * e
        norm = -mdot(v, v)
        if norm > 1e-12:
            basis.append(v / np.sqrt(norm))
        if len(basis) == 2:
            break
    worst = np.inf
    for _ in range(n):
        alpha = rng.choice([-1, 1]) * 10 ** rng.uniform(-5, np.log10(0.4))
        beta = rng.choice([-1, 1]) * 10 ** rng.uniform(-5, np.log10(0.4))
        rho = 10 ** rng.uniform(-10, -1)
        th = rng.uniform(0, 2 * np.pi)
        x0 = (C.x[j] + alpha * pm + beta * pp
              + np.sqrt(rho) * (np.cos(th) * basis[0] + np.sin(th) * basis[1]))
        X, _ = C.displaced_invariants(x0)
        scales = {j: max(rho, abs(alpha * beta * D)),
                  (j - 1) % 4: max(rho, abs(beta * D) * abs(1 + alpha)),
                  (j + 1) % 4: max(rho, abs(alpha * D) * abs(1 - beta))}
        worst = min(worst, min(abs(X[i]) / sc for i, sc in scales.items()))
    return worst / C.nu


# --------------------------------------------------------------------------
# Elementary bounds on the Gamma function (used in place of Stirling)
# --------------------------------------------------------------------------
def gamma_product_bounds(a, v):
    """(|Gamma(a + i v)|^2, bound) with the bound from the Weierstrass product.

    |Gamma(a + i v)|^2 = Gamma(a)^2 / prod_{n>=0} (1 + v^2/(a+n)^2). Comparing
    factor by factor with prod_{n>=1} (1 + v^2/n^2) = sinh(pi v)/(pi v):
        0 < a <= 1:  |Gamma(a+iv)|^2 <= Gamma(a)^2 pi v / sinh(pi v)
        1 <= a <= 2: |Gamma(a+iv)|^2 <= Gamma(a)^2 (1 + v^2) pi v / sinh(pi v)
    """
    import mpmath as mp

    a, v = mp.mpf(a), mp.mpf(v)
    value = abs(mp.gamma(mp.mpc(a, v))) ** 2
    ratio = mp.pi * v / mp.sinh(mp.pi * v) if v != 0 else mp.mpf(1)
    bound = mp.gamma(a) ** 2 * ratio * (1 if a <= 1 else 1 + v**2)
    return value, bound
