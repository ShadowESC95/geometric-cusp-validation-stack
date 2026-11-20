"""An explicit integration contour for the frozen-loop integral.

See docs/GATE_6_CONTOUR.md. The real integration domain is displaced by a
purely imaginary vector field,

    x0  ->  x0 + i V(x0),        V real,

so that  x0i^2 -> x0i^2 - V^2 + 2 i V.(x0 - x_i)  exactly. The imaginary part
is linear in V, so any convex combination of admissible fields is admissible.

Building blocks (written for D = 2 p_{j-1}.p_j < 0; where a diagonal is
timelike the longitudinal part of the frame field changes sign):

* Frame fields. At the cusp x_j write x0 = x_j + alpha p_{j-1} + beta p_j
  + yperp and let T be the triangle {beta >= 0, alpha <= 0, beta - alpha <= 1}
  in the (beta, alpha) plane, the projection of conv(x_{j-1}, x_j, x_{j+1}).
  With d = grad Psi(P), P = (beta, alpha), Psi a C^2 convex function that
  vanishes exactly on T,
      V_j = kappa (d_beta p_{j-1} + d_alpha p_j) + nu yperp
  satisfies, for i in {j-1, j, j+1},
      V_j.(x0 - x_i) = (kappa D/2) d.(P - vertex_i) - nu rho
                     <= -(kappa |D|/2) Psi(P) - nu rho <= 0,
  by convexity: Psi(v) >= Psi(P) + d.(v - P) and Psi(v) = 0. The field is
  continuously differentiable, so the ordinary Stokes theorem applies.

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


def _phi(tv):
    """C^2 convex function: t^3/3 on [0, 1], t^2 - t + 1/3 beyond, 0 for t <= 0."""
    tv = max(tv, 0.0)
    return tv**3 / 3 if tv <= 1 else tv * tv - tv + 1.0 / 3


def _dphi(tv):
    """Its derivative: t^2 on [0, 1], 2t - 1 beyond; C^1, with 0 <= phi'' <= 2."""
    tv = max(tv, 0.0)
    return tv * tv if tv <= 1 else 2 * tv - 1


# T = {l_k <= 0}: l_1 = -beta, l_2 = alpha, l_3 = beta - alpha - 1, in coordinates (beta, alpha)
_FACES = (np.array([-1.0, 0.0]), np.array([0.0, 1.0]), np.array([1.0, -1.0]))
_OFFSETS = (0.0, 0.0, -1.0)
HESSIAN_BOUND = 2 * sum(float(n @ n) for n in _FACES)          # = 8


def triangle_potential(P):
    """Psi(P) = sum_k phi(l_k(P)): convex, C^2, zero exactly on the triangle T."""
    return sum(_phi(float(n @ P) + c) for n, c in zip(_FACES, _OFFSETS))


def triangle_gradient(P):
    """d(P) = grad Psi(P): C^1, zero exactly on T, Lipschitz with constant 8.

    By convexity d.(P - v) >= Psi(P) for every v in T, and because the
    Hessian is at most 8, |d|^2 <= 16 Psi.
    """
    return sum(_dphi(float(n @ P) + c) * n for n, c in zip(_FACES, _OFFSETS))


@dataclass
class Contour:
    kappa: float = 2e-8          # longitudinal displacement of the frame fields
    nu: float = 1.5e-7           # transverse displacement of the frame fields
    kappa_far: float = 2e-3      # far field
    delta0: float = 0.025        # width of the neighbourhood of a cone, inside the ball
    radius: float = 8.0          # frame fields inside, far fields outside 2*radius
    tau0: float = 8.0            # time offset of the far-field centres
    angle: float = 1.1
    region: str = "--"           # signs of the diagonals (s, t)
    p: np.ndarray = field(init=False)
    x: np.ndarray = field(init=False)

    def __post_init__(self):
        self.p, self.s, self.t = light_like_square(self.region, self.angle)
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
        d = triangle_gradient(np.array([beta, alpha]))
        # the longitudinal part needs (kappa D) < 0: reverse it where the diagonal is timelike
        return -np.sign(D) * self.kappa * (d[0] * pm + d[1] * pp) + self.nu * yperp

    def frame_bound(self, j, x0):
        """-(kappa |D|/2) Psi(P) - nu rho: upper bound on V_j.(x0 - x_i), i = j-1, j, j+1."""
        pm, pp, D, alpha, beta, yperp, rho = self.frame(j, x0)
        return -(self.kappa * abs(D) / 2) * triangle_potential(np.array([beta, alpha])) - self.nu * rho

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


# --------------------------------------------------------------------------
# Certificate: every constant used in the proof of Proposition 6
# --------------------------------------------------------------------------
def four_cone_gap(contour: Contour):
    """A value eps such that no real point has all four |x0i^2| <= eps.

    In the frame of a cusp, max_i |x0i^2| <= eps forces |beta D|, |alpha D|
    <= 2 eps, |w| <= eps and rho <= eps + 4 eps^2/|D|, hence
        |x0,far^2 - D'| <= 4 eps |D'|/|D| + eps + 2 sqrt((D' u/D)(eps + 4 eps^2/|D|)).
    If the right side is below |D'| - eps this contradicts |x0,far^2| <= eps.
    Returns the largest eps (to 1e-6) for which the contradiction holds in
    some frame.
    """
    C = contour
    u = -C.s - C.t
    best = 0.0
    for j in range(4):
        D = 2 * mdot(C.p[(j - 1) % 4], C.p[j])
        Dp = C.t if abs(D - C.s) < 1e-9 else C.s
        lo, hi = 0.0, abs(Dp)

        def contradiction(eps):
            rho = eps + 4 * eps**2 / abs(D)
            shift = 4 * eps * abs(Dp) / abs(D) + eps + 2 * np.sqrt(Dp * u / D * rho)
            return shift < abs(Dp) - eps

        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if contradiction(mid) else (lo, mid)
        best = max(best, lo)
    return best


def certificate(contour: Contour) -> dict:
    """Evaluate the explicit conditions of GATE_6_CONTOUR section 9.

    Norms are Euclidean; |a.b| <= |a||b| for the Minkowski product. The
    returned dict has the geometric constants, each bound, and a boolean per
    condition. `all` is their conjunction.
    """
    C = contour
    eu = np.linalg.norm
    Delta = max(eu(x - C.centre) for x in C.x)
    Pn = max(eu(p) for p in C.p)
    Ds = [abs(2 * mdot(C.p[(j - 1) % 4], C.p[j])) for j in range(4)]
    A1 = 2 * Pn / min(Ds)                 # |alpha|, |beta| <= A1 |x0 - x_j|
    A2 = 1 + 2 * A1 * Pn                  # |yperp| <= A2 |x0 - x_j|
    # |d| <= 8 |P| and Lip(d) <= 8 (Hessian bound), |P| <= sqrt(2) A1 |x0 - x_j|
    lam0 = 2 * HESSIAN_BOUND * C.kappa * A1 * Pn + C.nu * A2   # |V_j| <= lam0 |x0 - x_j|, Lip(V_j) <= lam0
    R, d0, kf, tau0 = C.radius, C.delta0, C.kappa_far, C.tau0
    out = {"Delta": Delta, "P": Pn, "A1": A1, "A2": A2, "lambda0": lam0}

    # H1: some frame field is fully switched on near every cone
    out["four_cone_gap"] = four_cone_gap(C)
    out["H1_frame_on"] = bool(2 * d0 < out["four_cone_gap"])

    # C0: far centres are in the common past and future with margin
    out["C0_centres"] = bool(tau0 >= 3 * Delta)

    # C1: far field strictly admissible on |x0i^2| <= delta0, r >= R
    lhs = (tau0 - 2 * Delta) * (R / np.sqrt(17) - Delta)
    out["C1_far_admissible"] = bool(R / np.sqrt(17) > Delta and lhs > d0 + Delta * np.sqrt(d0))
    out["C1_margin"] = lhs - (d0 + Delta * np.sqrt(d0))

    # C2: for r >= R, points with |x0i^2| <= delta0 have |time| >= |space|/2, so Theta = 1
    xi0 = 2 * (2 * Delta + np.sqrt(d0))
    out["C2_theta_is_one"] = bool(R >= np.sqrt(5) / 2 * xi0)

    # C3: in the spacelike wedge |time| <= |space|/4, r >= R, no cone is within delta0
    xi_min = R * 4 / np.sqrt(17)
    worst = (xi_min / 4 + Delta) ** 2 - (xi_min - Delta) ** 2
    out["C3_wedge_free"] = bool(xi_min > Delta and worst < -d0)

    # D: excursion of the phases
    reach = 2 * R + Delta
    frame_term = 2 * lam0 * reach**2                       # |2 V_{i+2}.(x0 - x_i)|
    v_inner = 4 * lam0 * reach + kf * (2 * R + tau0 + Delta)
    out["excursion_frames"] = float(np.arctan(frame_term / (d0 - v_inner**2)))
    # far field: three cases (see the proof): Im > 0 with n.u >= 0; with n.u < 0; partial Theta
    xi = R * 2 / np.sqrt(5)                                  # |space| >= this where 0 < Theta < 1
    case_a = 1 / (1 - kf**2)
    case_b = 1 + 2 * Delta / (3 * (R / np.sqrt(17) - Delta))
    case_c = 1 + (np.sqrt(5) / 2 + Delta / xi) * (tau0 + Delta) / (0.75 * (xi - 4 * Delta))
    out["C4_wedge_margin"] = bool(xi > 4 * Delta)
    # In the blend zone the real part is n^2 - V^2 >= |n^2|/2 instead of (1 - kf^2)|n^2|,
    # so the far term is at most twice the pure outer one: 4 kf covers both zones.
    out["excursion_far"] = float(np.arctan(4 * kf * max(1.0, case_b, case_c)))
    out["excursion_bound"] = out["excursion_frames"] + out["excursion_far"]
    out["D_excursion_small"] = bool(v_inner**2 < d0 / 2 and out["excursion_bound"] < 0.2)
    # pure outer zone: |Re| >= |n^2|/2 in each of the three cases, which the factor 4 kf assumes
    out["D_far_real_part"] = bool(
        kf**2 <= 0.5 and kf**2 * (2 * case_b - 1) <= 0.5
        and kf**2 * (np.sqrt(5) / 2 * xi + 2 * Delta + tau0) ** 2 <= 0.375 * xi * (xi - 4 * Delta))

    # E: the map is a local diffeomorphism, ||dV|| < 1
    j_inner = 4 * (3 * reach**2 * lam0 / d0 + lam0)        # sum_j |dm_j||V_j| + m_j |dV_j|
    j_outer = kf * (1 + 7.5 * (1 + (tau0 + Delta) / R))
    j_blend = 1.5 / R * (4 * lam0 * reach + kf * (2 * R + tau0 + Delta))
    out["jacobian_bound"] = j_inner + j_outer + j_blend
    out["E_local_diffeomorphism"] = bool(out["jacobian_bound"] < 1)

    # F: Lemma 6 needs theta <= 1
    out["theta"] = theta_constant(C)
    out["F_theta_small"] = bool(out["theta"] <= 1)

    # S: the strips of Lemma 7 (same widths as `strip_check`)
    Dmin, Dmax = min(Ds), max(Ds)
    u_abs = abs(C.s + C.t)
    b0 = 4 * d0 / Dmin
    a1 = d0 / (8 * Dmax)
    rho1 = (d0 / 8) ** 2 / max(1.0, Dmax * u_abs / Dmin) / 4
    B = max(transverse_norm(C, j) for j in range(4))
    out.update({"b0": b0, "a1": a1, "rho1": rho1})
    out["S0_strips_overlap"] = bool(b0 < 0.5)
    # |x0,j^2| = |alpha beta D - rho| and |x0,j+1^2| = |alpha D (beta - 1) - rho| on the strip
    out["S2_near_invariants"] = bool(a1 * Dmax * (1 + a1) + rho1 <= d0)
    # |x0,j+2^2| >= |D'|(b0 - a1) - a1 |D| - rho1 - 2 sqrt(rho1 |D' u / D|)
    far_lower = Dmin * (b0 - a1) - a1 * Dmax - rho1 - 2 * np.sqrt(rho1 * Dmax * u_abs / Dmin)
    out["S1_far_margin"] = far_lower - 2 * d0
    out["S1_frame_fully_on"] = bool(far_lower >= 2 * d0)
    # the strips lie in the inner zone
    out["S3_strips_inside"] = bool(Delta + a1 * Pn + B * np.sqrt(rho1) <= R)

    # I: Lemma 9, the bound at infinity
    omega = 0.0
    for _ in range(50):                                   # fixed point of omega = c_u / 2
        c_u = (tau0 - 2 * Delta) * (1 / np.sqrt(17) - Delta / (2 * R)) - Delta * np.sqrt(omega / (2 * R))
        omega = c_u / 2
    out["c_u"], out["omega"] = float(c_u), float(omega)
    out["I0_c_u_positive"] = bool(c_u > 0)
    r0 = 2 * R
    # kappa'^2 (2 |n||u| + |u|^2) <= omega r / 4 for r >= 2R (the left side grows linearly, so check at 2R
    # and compare slopes)
    lhs = kf**2 * (2 * (r0 + Delta) * (tau0 + Delta) + (tau0 + Delta) ** 2)
    out["I1_real_part"] = bool(lhs <= omega * r0 / 4 and 2 * kf**2 * (tau0 + Delta) <= omega / 4)
    # |n^2| < omega r forces |time| >= |space|/2: 2(2 Delta + sqrt(omega r)) <= (2/sqrt 5) r for r >= 2R
    out["I2_theta_is_one"] = bool(2 * (2 * Delta + np.sqrt(omega * r0)) <= 2 / np.sqrt(5) * r0
                                  and r0 > 5 * omega / 4)

    # Q: the corner arg z -> 3 pi (Lemma Q). Near a cone, |V^2| <= mu_Q |Im x0i^2|.
    mu = max(64 * C.kappa * Pn**2 / Dmin, C.nu * B**2)        # |V_j|^2 <= mu (kappa |D| Psi + 2 nu rho)
    M = out["C1_margin"]
    K = kf * (2 * R + tau0 + Delta)                           # |V_out| in the ball of radius 2R
    if M > 0:
        mu_far = kf * (1 + (1.5 * d0 + (tau0 + Delta) ** 2) / (2 * M))
        mu_q = max(3 * mu + np.sqrt(3 * mu), mu_far + np.sqrt(3 * mu) * K**2 / (2 * kf * M))
    else:
        mu_q = np.inf
    gamma = PHASE_CORNER_MARGIN + 3 * out["excursion_bound"]
    out["mu_Q"], out["corner_angle"] = float(mu_q), float(gamma)
    out["Q1_corner_excluded"] = bool(gamma < np.pi / 2 and np.tan(gamma) * mu_q < 1)
    out["Q3_corner_excluded_far"] = bool(
        (2 * R / np.sqrt(17) - Delta) * (tau0 - 2 * Delta) > 2 * kf**2 * Delta**2 / (1 - kf**2))

    # Q4: used in Lemma 10(c): |V_0| < |n_0| and |n| large enough in the outer zone
    n0_min = R / np.sqrt(17) - Delta
    out["Q4_outer_zone"] = bool(kf * (tau0 + Delta) < (1 - kf) * n0_min
                                and 2 * n0_min > 8 * kf**2 * Delta / (3 * (1 - kf**2)))

    out["all"] = all(v for k, v in out.items() if isinstance(v, bool))
    return out


PHASE_CORNER_MARGIN = 0.3      # with a timelike diagonal, |arg z| <= 3 pi - this on the contour


def corner_ratio(contour: Contour, n: int, seed: int = 0) -> float:
    """Largest |V^2| / |Im x0i^2| over sampled points with |x0i^2| <= delta0/2
    inside the ball of radius 2R. Lemma Q bounds it by mu_Q."""
    C = contour
    worst = 0.0
    for x0 in sample_points(C, n, seed):
        if np.linalg.norm(x0 - C.centre) > 2 * C.radius:
            continue
        real = C.invariants(x0)
        v = C.field(x0)
        v2 = abs(mdot(v, v))
        for i in range(4):
            if abs(real[i]) <= C.delta0 / 2:
                im = abs(2 * mdot(v, x0 - C.x[i]))
                if im > 0:
                    worst = max(worst, v2 / im)
                elif v2 > 0:
                    return np.inf
    return float(worst)


# --------------------------------------------------------------------------
# Constants chosen from the geometry
# --------------------------------------------------------------------------
def tuned_contour(region: str = "--", angle: float = 1.1, scale: float = 1.0,
                  jacobian_target: float = 0.3) -> Contour:
    """A contour whose constants are derived from the square, in this order.

    1. tau0 = 3 Delta + 5 Delta: the far centres, condition (C0).
    2. R: the smallest radius meeting (C1), (C2), (C3) and the wedge margin,
       with delta0 provisionally bounded by the four-cone gap.
    3. delta0 = four_cone_gap / 2.5: condition (H1).
    4. kappa' : small enough for v^2 < delta0/2 in (D).
    5. kappa = nu: small enough that the bound (E) on ||dV|| is
       `jacobian_target`, which also makes the frame excursion in (D) small.

    `scale` multiplies every length; the constants scale with it so that the
    conditions, which are homogeneous, are unaffected.
    """
    base = Contour(region=region, angle=angle)
    eu = np.linalg.norm
    Delta = max(eu(x - base.centre) for x in base.x)
    Pn = max(eu(p) for p in base.p)
    Ds = [abs(2 * mdot(base.p[(j - 1) % 4], base.p[j])) for j in range(4)]
    gap = four_cone_gap(base)
    d0 = gap / 2.5
    tau0 = 8 * Delta
    R = max(8.0 * Delta, np.sqrt(5) * (2 * Delta + np.sqrt(d0)) * 1.05, 2.5 * np.sqrt(5) * Delta)
    while True:
        ok1 = R / np.sqrt(17) > Delta and \
            (tau0 - 2 * Delta) * (R / np.sqrt(17) - Delta) > 1.5 * (d0 + Delta * np.sqrt(d0))
        xi_min = R * 4 / np.sqrt(17)
        ok3 = xi_min > Delta and (xi_min / 4 + Delta) ** 2 - (xi_min - Delta) ** 2 < -1.5 * d0
        if ok1 and ok3:
            break
        R *= 1.25
    reach = 2 * R + Delta
    kf = min(2e-3, 0.4 * np.sqrt(d0 / 2) / (2 * R + tau0 + Delta))
    A1 = 2 * Pn / min(Ds)
    A2 = 1 + 2 * A1 * Pn
    j_rest = kf * (1 + 7.5 * (1 + (tau0 + Delta) / R)) + 1.5 / R * kf * (2 * R + tau0 + Delta)
    budget = max(jacobian_target - j_rest, 0.05)
    lam0 = budget / (4 * (3 * reach**2 / d0 + 1) + 6 * reach / R)
    amp = lam0 / (2 * HESSIAN_BOUND * A1 * Pn + A2)
    out = Contour(kappa=amp, nu=amp, kappa_far=kf, delta0=d0, radius=R, tau0=tau0,
                  angle=angle, region=region)
    if scale != 1.0:
        out.p = out.p * scale
        out.x = out.x * scale
        out.centre = out.centre * scale
        out.c_past, out.c_future = out.c_past * scale, out.c_future * scale
        out.s, out.t = out.s * scale**2, out.t * scale**2
        out.delta0 *= scale**2
        out.radius *= scale
        out.tau0 *= scale
    return out


# --------------------------------------------------------------------------
# Absolute integrability: constants and sampled checks (GATE_6_CONTOUR sec. 11)
# --------------------------------------------------------------------------
def transverse_norm(contour: Contour, j: int) -> float:
    """B_j: the largest Euclidean norm of a transverse vector with yperp^2 = -1."""
    pm, pp = contour.p[(j - 1) % 4], contour.p[j]
    D = 2 * mdot(pm, pp)
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
    gram = np.array([[a @ b for b in basis] for a in basis])
    return float(np.sqrt(np.linalg.eigvalsh(gram).max()))


def theta_constant(contour: Contour) -> float:
    """theta with |V^2| <= (theta/2) |Im x0i^2| wherever every active frame controls i.

    |V|^2 <= 4 sum_j m_j (4 kappa^2 P^2 |d_j|^2 + 2 nu^2 B^2 rho_j)
          <= theta sum_j m_j b_j,   b_j = (kappa |D_j|/2)|d_j|^2 + nu rho_j,
    With the smooth potential, |d_j|^2 <= 16 Psi_j and b_j = (kappa |D_j|/2) Psi_j + nu rho_j,
    so theta = max(512 kappa P^2 / min|D|, 8 nu B^2), and |Im x0i^2| >= 2 sum_j m_j b_j.
    """
    C = contour
    Pn = max(np.linalg.norm(p) for p in C.p)
    Dmin = min(abs(2 * mdot(C.p[(j - 1) % 4], C.p[j])) for j in range(4))
    B = max(transverse_norm(C, j) for j in range(4))
    return float(max(512 * C.kappa * Pn**2 / Dmin, 8 * C.nu * B**2))


def lemma6_check(contour: Contour, n: int, seed: int = 0) -> dict:
    """Sampled check of Lemma 6 in the inner zone.

    For r <= R and |x0i^2| <= delta0:
        |V^2| <= (theta/2) |Im x0i^2|   and   |x0i^2 displaced| >= (1/2) max(|x0i^2|, |Im|).
    Returns the largest |V^2| / |Im| and the smallest |displaced| / max(|real|, |Im|).
    """
    C = contour
    worst_v, worst_ratio, count = 0.0, np.inf, 0
    for x0 in sample_points(C, n, seed):
        if np.linalg.norm(x0 - C.centre) > C.radius or C.distance_to_edges(x0) < 1e-9:
            continue
        real = C.invariants(x0)
        displaced, v = C.displaced_invariants(x0)
        v2 = abs(mdot(v, v))
        for i in range(4):
            if abs(real[i]) <= C.delta0:
                im = abs(displaced[i].imag)
                if im == 0.0:
                    continue
                count += 1
                worst_v = max(worst_v, v2 / im)
                worst_ratio = min(worst_ratio, abs(displaced[i]) / max(abs(real[i]), im))
    return {"count": count, "max_v2_over_im": worst_v, "min_ratio": worst_ratio,
            "theta": theta_constant(C)}


def strip_check(contour: Contour, n: int, seed: int = 0) -> dict:
    """Sampled check of the bounds used on the strips around the edges.

    In the frame of x_j, on the strip around the edge to x_{j+1}
    (-a1 <= beta <= 1 - b0, |alpha| <= a1, rho <= rho1):
        |x0,j^2|   >= nu rho
        |x0,j+1^2| >= (nu/4) max(rho, |alpha D| b0 / 2)
        |x0,j-1^2| >= (nu/4) max(rho, |beta D| / 2)
        |x0,j+2^2| >= delta0
    Returns the smallest ratio of each left side to its right side.
    """
    rng = np.random.default_rng(seed)
    C = contour
    Dall = [abs(2 * mdot(C.p[(j - 1) % 4], C.p[j])) for j in range(4)]
    b0 = 4 * C.delta0 / min(Dall)
    a1 = C.delta0 / (8 * max(Dall))
    rho1 = (C.delta0 / 8) ** 2 / max(1.0, max(Dall) * abs(C.s + C.t) / min(Dall)) / 4
    out = {"b0": b0, "a1": a1, "rho1": rho1, "j": np.inf, "next": np.inf, "prev": np.inf,
           "far": np.inf, "frame_on": 1.0}
    for _ in range(n):
        j = int(rng.integers(0, 4))
        pm, pp = C.p[(j - 1) % 4], C.p[j]
        D = 2 * mdot(pm, pp)
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
        beta = rng.uniform(-a1, 1 - b0) if rng.random() < 0.6 else rng.choice([-1, 1]) * a1 * 10 ** rng.uniform(-6, 0)
        beta = min(max(beta, -a1), 1 - b0)
        alpha = rng.choice([-1, 1]) * a1 * 10 ** rng.uniform(-8, 0)
        rho = rho1 * 10 ** rng.uniform(-10, 0)
        th = rng.uniform(0, 2 * np.pi)
        x0 = (C.x[j] + alpha * pm + beta * pp
              + np.sqrt(rho) * (np.cos(th) * basis[0] + np.sin(th) * basis[1]))
        X, _ = C.displaced_invariants(x0)
        out["frame_on"] = min(out["frame_on"], C.weights(x0)[j])
        out["j"] = min(out["j"], abs(X[j]) / (C.nu * rho))
        out["next"] = min(out["next"], abs(X[(j + 1) % 4])
                          / (C.nu / 4 * max(rho, abs(alpha * D) * b0 / 2)))
        out["prev"] = min(out["prev"], abs(X[(j - 1) % 4])
                          / (C.nu / 4 * max(rho, abs(beta * D) / 2)))
        out["far"] = min(out["far"], abs(X[(j + 2) % 4]) / C.delta0)
    return out


def infinity_check(contour: Contour, n: int, seed: int = 0) -> float:
    """Smallest |x0i^2 displaced| / (kappa' max(r, |(x0 - c)^2|)) for r >= 2R."""
    rng = np.random.default_rng(seed)
    C = contour
    worst = np.inf
    for _ in range(n):
        r = 2 * C.radius * 10 ** rng.uniform(0, 3)
        kind = rng.integers(0, 3)
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        if kind == 0:                                   # generic direction
            v = rng.normal(size=4)
            y = r * v / np.linalg.norm(v)
        else:                                           # close to the asymptotic cone
            sign = rng.choice([-1.0, 1.0])
            off = rng.normal() * 10 ** rng.uniform(-3, 1)
            y = np.array([sign * (r / np.sqrt(2) + off), *(direction * r / np.sqrt(2))])
        x0 = C.centre + y
        rr = np.linalg.norm(y)
        if rr < 2 * C.radius:
            continue
        X, _ = C.displaced_invariants(x0)
        worst = min(worst, np.min(np.abs(X)) / (C.kappa_far * max(rr, abs(mdot(y, y)))))
    return float(worst)


# --------------------------------------------------------------------------
# The contour in D = 4 + k dimensions (extra coordinates e)
# --------------------------------------------------------------------------
# Dimensional regularisation adds extra spatial coordinates e, and every
# invariant becomes x0i^2 - e.e. The contour must displace e as well: with a
# real shift mu = e.e and the four-dimensional field, the surface
# x0i^2 = mu is not passed on the Feynman side once mu exceeds delta0. The
# construction is the same in R^(1, 3+k): the transverse space of each frame
# gains the coordinates e, so rho -> rho + e.e and the frame field gets the
# component nu e; the far fields get -kappa' e.

def invariants_extended(contour: Contour, z, e):
    """x0i^2 - e.e for complex z in C^4 and e in C^k (bilinear)."""
    e = np.atleast_1d(e)
    return contour.invariants(z) - e @ e


def field_extended(contour: Contour, x0, e):
    """The displacement (V4, Ve) at the real point (x0, e) of R^(1, 3+k)."""
    C = contour
    e = np.atleast_1d(np.asarray(e, dtype=float))
    mu = e @ e
    X = C.invariants(x0) - mu
    y = x0 - C.centre
    r = np.sqrt(y @ y + mu)
    g = smoothstep(r / C.radius - 1.0)
    v4, ve = np.zeros(4), np.zeros_like(e)
    if g < 1:
        for j in range(4):
            m = smoothstep(abs(X[(j + 2) % 4]) / C.delta0 - 1.0)
            if m > 0:
                v4 = v4 + (1 - g) * m * C.frame_field(j, x0)
                ve = ve + (1 - g) * m * C.nu * e
    if g > 0:
        space = np.sqrt(y[1:] @ y[1:] + mu)
        theta = y[0] / max(space, 1e-300)
        up = smoothstep((theta - 0.25) / 0.25)
        down = smoothstep((-theta - 0.25) / 0.25)
        v4 = v4 + g * (-C.kappa_far) * (up * (x0 - C.c_past) + down * (x0 - C.c_future))
        ve = ve + g * (-C.kappa_far) * (up + down) * e
    return v4, ve


def displaced_extended(contour: Contour, x0, e, lam: float = 1.0):
    """Invariants at (x0 + i lam V4, e + i lam Ve), and the field."""
    v4, ve = field_extended(contour, x0, e)
    e = np.atleast_1d(np.asarray(e, dtype=float))
    return invariants_extended(contour, x0 + 1j * lam * v4, e + 1j * lam * ve), (v4, ve)


def sample_extended(contour: Contour, n: int, seed: int = 0):
    """Points (x0, e) with one extra coordinate, weighted towards the D-dimensional cones."""
    rng = np.random.default_rng(seed)
    for x0 in sample_points(contour, n, seed):
        X = contour.invariants(x0)
        kind = rng.integers(0, 4)
        if kind == 0:
            e = 10 ** rng.uniform(-6, 1.5)
        elif kind == 1 and np.any(X > 0):                   # exactly on a D-dimensional cone
            e = np.sqrt(rng.choice(X[X > 0]))
        elif kind == 2 and np.any(X > 0):                   # close to one
            e = np.sqrt(max(rng.choice(X[X > 0]) + rng.normal() * 10 ** rng.uniform(-6, -1), 0.0))
        else:
            e = 0.0
        yield x0, np.array([e])


def survey_extended(contour: Contour, n: int, seed: int = 0) -> dict:
    """The properties of `survey`, for the contour in 4 + 1 dimensions."""
    C = contour
    out = {"max_excursion": 0.0, "max_im_near": -np.inf, "near_count": 0, "min_ratio": np.inf,
           "largest_mu_near": 0.0}
    for x0, e in sample_extended(C, n, seed):
        mu = float(e @ e)
        dist = np.sqrt(C.distance_to_edges(x0) ** 2 + mu)
        if dist < 1e-5:
            continue
        real = C.invariants(x0) - mu
        displaced, _ = displaced_extended(C, x0, e)
        r = np.sqrt((x0 - C.centre) @ (x0 - C.centre) + mu)
        tol = C.delta0 / 2 if r < C.radius else 0.1 * r
        for i in range(4):
            out["max_excursion"] = max(out["max_excursion"], C.excursion(displaced[i]))
            if abs(real[i]) <= tol:
                out["near_count"] += 1
                out["max_im_near"] = max(out["max_im_near"], displaced[i].imag)
                out["largest_mu_near"] = max(out["largest_mu_near"], mu)
        out["min_ratio"] = min(out["min_ratio"], np.min(np.abs(displaced)) / min(1.0, dist**2))
    return out


def shifted_surface_sides(contour: Contour, mu: float, n: int, seed: int = 0) -> dict:
    """Points with x0i^2 = mu: how many are displaced to the wrong side.

    'real_shift': the four-dimensional field with the invariants shifted by
    the real number mu. 'extended': the field in 4 + 1 dimensions at
    e = sqrt(mu). Only the second is the Feynman prescription for all mu.
    """
    rng = np.random.default_rng(seed)
    C = contour
    out = {"points": 0, "real_shift": 0, "extended": 0}
    for _ in range(n):
        i = int(rng.integers(0, 4))
        v = rng.normal(size=3) * rng.uniform(0, 2)
        x0 = C.x[i] + np.array([np.sqrt(mu + v @ v) * rng.choice([-1, 1]), *v])
        if np.linalg.norm(x0 - C.centre) > C.radius:
            continue
        out["points"] += 1
        out["real_shift"] += int(C.displaced_invariants(x0)[0][i].imag > 0)
        out["extended"] += int(displaced_extended(C, x0, np.array([np.sqrt(mu)]))[0][i].imag > 0)
    return out


# --------------------------------------------------------------------------
# The phase of z itself on the contour
# --------------------------------------------------------------------------
def log_z_imaginary_part(contour: Contour, displaced):
    """Im log z = arg x02^2 + arg x04^2 - arg x01^2 - arg x03^2 + pi (tau_s - tau_t),
    each argument on the branch continued from the Feynman prescription."""
    a = [contour.branch_argument(v) for v in displaced]
    return a[1] + a[3] - a[0] - a[2] + np.pi * ((contour.s > 0) - (contour.t > 0))


def phase_survey(contour: Contour, n: int, seed: int = 0) -> dict:
    """Range of Im log z over a stress sample: at real points and on the contour.

    At real points it is a multiple of pi in [-2 pi, 2 pi] (Lemma P). On the
    contour it must stay within 2 pi plus a small margin.
    """
    C = contour
    out = {"real_min": np.inf, "real_max": -np.inf, "contour_min": np.inf, "contour_max": -np.inf,
           "points": 0}
    for x0 in sample_points(C, n, seed):
        if C.distance_to_edges(x0) < 1e-5:
            continue
        real = C.invariants(x0)
        if np.min(np.abs(real)) > 1e-12:
            v = log_z_imaginary_part(C, real.astype(complex) - 1e-300j) / np.pi
            out["real_min"], out["real_max"] = min(out["real_min"], v), max(out["real_max"], v)
        displaced, _ = C.displaced_invariants(x0)
        v = log_z_imaginary_part(C, displaced) / np.pi
        out["contour_min"], out["contour_max"] = min(out["contour_min"], v), max(out["contour_max"], v)
        out["points"] += 1
    return {k: float(v) for k, v in out.items()}


def strip_power_integral(eta: float, sigma: float, cutoff: float) -> float:
    """int_cutoff^1 d rho rho^(eta - a_j) g(rho, a_+) g(rho, a_-), with
    g(rho, a) = int_0^1 d alpha max(rho, alpha)^(-a), a_j = 1 - sigma and
    a_+ = a_- = 1 + sigma: the strip integral of Lemma 8 for propagator
    powers 1 -+ sigma, using the lower bounds of Lemma 7. It has a finite
    limit as cutoff -> 0 if and only if |sigma| < eta.
    """
    from scipy.integrate import quad

    def g(rho, a):
        if abs(a - 1) < 1e-12:
            return 1 + np.log(1 / rho)
        return rho ** (1 - a) + (1 - rho ** (1 - a)) / (1 - a)

    a_j, a_n = 1 - sigma, 1 + sigma
    f = lambda x: np.exp(-x * (eta - a_j + 1)) * g(np.exp(-x), a_n) ** 2      # rho = exp(-x)
    return quad(f, 0, np.log(1 / cutoff), epsabs=0, epsrel=1e-10, limit=400)[0]
