"""The cusp-local model of the frozen-loop integral, and the vertex lemma.

Near the cusp x_2 the fourth invariant is frozen at its cusp value t. What
is left is a triangle with two light-like sides, which has a closed form.
Its double pole is the product of the collinear poles of the two edges that
meet at the cusp, with coefficient sin(pi p)/(pi p) for the weight z^p.

docs/CUSP_LOCAL_MODEL.md has the statements and proofs.
"""
from __future__ import annotations

import mpmath as mp
import sympy as sp


def measure_identity():
    """Exact change of variables at a cusp (|s| = 1).

    With x02^2 = w, x01^2 = u1 = w + q1, x03^2 = u3 = w + q3 and
    w = q1 q3 - rho, use (u1, u3, z) as coordinates, z = -w/(u1 u3). Then
        dq1 dq3 drho / (w u1 u3) = (du1/u1) (du3/u3) (dz/z)
    up to sign, and rho = u1 u3 (1 + z + z (u1 + u3) + z^2 u1 u3).
    Returns (jacobian ratio, rho residual); both are exact zeros/ones.
    """
    u1, u3, z = sp.symbols("u1 u3 z")
    w = -z * u1 * u3
    q1, q3 = u1 - w, u3 - w
    rho = q1 * q3 - w
    jac = sp.Matrix([q1, q3, rho]).jacobian([u1, u3, z]).det()
    ratio = sp.simplify(jac / (w * u1 * u3) * (u1 * u3 * z))        # should be +-1
    residual = sp.expand(rho - u1 * u3 * (1 + z + z * (u1 + u3) + z**2 * u1 * u3))
    return ratio, residual


def triangle_closed_form(p, eta):
    """eta^2 times the light-like triangle with powers (1+p, 1-p, 1+p),
    normalised so that the box integral is N_eta: the contribution of one
    cusp to eps^2 s t x (regulated integral of z^p), at S = T = 1.

        eta^2 Gamma(eta)^2 Gamma(1 + p - eta) / (Gamma(1+p)^2 Gamma(1 - p + 2 eta))

    Each Gamma(eta) is the collinear pole of one of the two edges at the
    cusp. The limit eta -> 0 is sin(pi p)/(pi p).
    """
    p, eta = mp.mpmathify(p), mp.mpmathify(eta)
    return (eta**2 * mp.gamma(eta) ** 2 * mp.gamma(1 + p - eta)
            * mp.rgamma(1 + p) ** 2 * mp.rgamma(1 - p + 2 * eta))


def dirichlet_check(p, eta):
    """The Feynman-parameter integral of the triangle, done numerically:
    int over the simplex of a1^(eta-1) a3^(eta-1) a2^(-p), against
    Gamma(eta)^2 Gamma(1-p)/Gamma(1 - p + 2 eta)."""
    p, eta = mp.mpf(p), mp.mpf(eta)
    # a1 = r x, a3 = r (1 - x), a2 = 1 - r
    radial = mp.quad(lambda r: r ** (2 * eta - 1) * (1 - r) ** (-p), [0, 0.5, 1])
    angular = mp.beta(eta, eta)
    closed = mp.gamma(eta) ** 2 * mp.gamma(1 - p) / mp.gamma(1 - p + 2 * eta)
    return radial * angular, closed


def vertex_main_term(p, eta, delta, S=1, T=1):
    """The model integral at the vertex of the Feynman simplex belonging to
    the cusp x_2, over the cube a1, a3 <= delta with a4 running to infinity:

        int (a1 a3)^p a4^(-p) (a1 a3 S + a4 T)^(eta - 2)
          = T^(p-1) S^(eta-1-p) Gamma(1-p) Gamma(1+p-eta)/Gamma(2-eta) * delta^(2 eta)/eta^2.

    Returns (numerical value of the a4 integral at a1 a3 S = c = 0.37, its
    closed form, the full main term)."""
    p, eta, delta = mp.mpmathify(p), mp.mpf(eta), mp.mpf(delta)
    S, T = mp.mpf(S), mp.mpf(T)
    c = mp.mpf("0.37")
    numeric = mp.quad(lambda a: a ** (-p) * (c + a * T) ** (eta - 2), [0, c / T, 1, mp.inf])
    g = mp.gamma(1 - p) * mp.gamma(1 + p - eta) / mp.gamma(2 - eta)
    closed = c ** (eta - 1 - p) * T ** (p - 1) * g
    main = T ** (p - 1) * S ** (eta - 1 - p) * g * delta ** (2 * eta) / eta**2
    return numeric, closed, main


def remainder_ratio(p, eta, ratio=1):
    """(N_eta(p) - 4 sin(pi p)/(pi p)) / eta, which the vertex lemma bounds."""
    from .lorentzian_edges import regulated_double_pole

    p = mp.mpmathify(p)
    limit = 4 * mp.sin(mp.pi * p) / (mp.pi * p) if p != 0 else mp.mpf(4)
    return (regulated_double_pole(p, eta, ratio) - limit) / mp.mpf(eta)


# --------------------------------------------------------------------------
# The exact formula for a general weight
# --------------------------------------------------------------------------
def circle_formula_prefactor(eta):
    """c_eta = Gamma(1+eta)^2 Gamma(1-eta)/Gamma(1+2 eta), which tends to 1."""
    eta = mp.mpf(eta)
    return mp.gamma(1 + eta) ** 2 * mp.gamma(1 - eta) / mp.gamma(1 + 2 * eta)


def smearing_density(x, eta):
    """Probability density on x >= 0:
    exp(-(1-eta) x) (1 - exp(-x))^(eta-1) / (Gamma(eta) Gamma(1-eta)).
    It concentrates at x = 0 as eta -> 0; its mean is of order eta."""
    x, eta = mp.mpf(x), mp.mpf(eta)
    return mp.exp(-(1 - eta) * x) * (-mp.expm1(-x)) ** (eta - 1) / (mp.gamma(eta) * mp.gamma(1 - eta))


def regulated_circle_formula(f, eta, nodes: int = 400):
    """eta^2 x (cusp-local model) for a general weight F(z) = f(log z):

        c_eta / (2 pi) * int_{-pi}^{pi} d phi (1 + e^{-i phi})^(2 eta)
                       * int_0^oo dx kappa_eta(x) f(i phi - x).

    Only the principal sheet |Im L| <= pi and |z| <= 1 enter. For f = exp(p L)
    this is `triangle_closed_form(p, eta)` exactly. `f` must accept numpy
    arrays of complex numbers. Gauss-Legendre quadrature; x = y^(1/eta) near
    0 removes the x^(eta-1) endpoint singularity.
    """
    import numpy as np
    from scipy.special import gamma

    eta = float(eta)
    t, wt = np.polynomial.legendre.leggauss(nodes)
    phi, wphi = np.pi * t, np.pi * wt
    y, wy = 0.5 * (t + 1), 0.5 * wt                      # x in (0, 1): x = y^(1/eta)
    x_near = y ** (1 / eta)
    ratio = np.where(x_near > 1e-12, -np.expm1(-x_near) / np.maximum(x_near, 1e-300), 1.0)
    w_near = wy / eta * np.exp(-(1 - eta) * x_near) * ratio ** (eta - 1)
    x_far = 1 + 0.5 * (t + 1) * 59                       # x in (1, 60)
    w_far = 0.5 * 59 * wt * np.exp(-(1 - eta) * x_far) * (-np.expm1(-x_far)) ** (eta - 1)
    x = np.concatenate([x_near, x_far])
    wx = np.concatenate([w_near, w_far]) / (gamma(eta) * gamma(1 - eta))
    values = f(1j * phi[:, None] - x[None, :])
    smeared = values @ wx
    total = np.sum(wphi * (1 + np.exp(-1j * phi)) ** (2 * eta) * smeared)
    return complex(circle_formula_prefactor(eta)) * total / (2 * np.pi)


def timelike_remainder_ratio(p, eta, modulus=1):
    """(N_eta(p) - 4 sin(pi p)/(pi p)) / eta at a timelike diagonal, arg(S/T) = pi."""
    from .lorentzian_edges import regulated_double_pole_complex

    p = mp.mpmathify(p)
    limit = 4 * mp.sin(mp.pi * p) / (mp.pi * p)
    return (regulated_double_pole_complex(p, eta, modulus, mp.pi) - limit) / mp.mpf(eta)


# --------------------------------------------------------------------------
# The Lorentzian derivation (docs/CUSP_LOCAL_MODEL.md section 8)
# --------------------------------------------------------------------------
def _feynman_log(x):
    """log on the branch continued from the Feynman prescription: arg in (-3pi/2, pi/2]."""
    import numpy as np
    a = np.angle(x)
    return np.log(abs(x)) + 1j * (a - 2 * np.pi if a > np.pi / 2 else a)


def light_cone_integral(f, beta, rho):
    """int over real alpha of f(log z)/(w u1 u3) at fixed beta and (complex) rho.

    |s| = 1: w = -alpha beta - rho, u1 = w - beta, u3 = w + alpha,
    z = -w/(u1 u3). Giving rho a phase e^{i psi}, 0 < psi < pi, is the
    Feynman prescription for all three invariants at once. The result is
    zero unless 0 < beta < 1.
    """
    import numpy as np
    from scipy.integrate import quad

    def g(alpha):
        w = -alpha * beta - rho
        u1, u3 = w - beta, w + alpha
        L = _feynman_log(w) - _feynman_log(u1) - _feynman_log(u3) - 1j * np.pi
        return f(L) / (w * u1 * u3)

    marks = sorted({float(np.real(-rho / beta)), float(np.real(-1 - rho / beta)),
                    float(np.real(rho / (1 - beta))), 0.0})
    edges = [-np.inf] + marks + [np.inf]
    total = 0j
    for lo, hi in zip(edges[:-1], edges[1:]):
        total += quad(lambda a: g(a).real, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
        total += 1j * quad(lambda a: g(a).imag, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
    return total


def discontinuity_integral(f, beta, rho, reach: float = 300.0):
    """The same integral after closing the alpha contour, for 0 < beta < 1:
    (1/(1-beta)) int_0^oo du3 [G(u3) - G(u3 e^{-2 pi i})].

    `f` must decay as Re L -> +oo. `reach` is the range of log u3; slowly
    decaying weights need a large one."""
    import numpy as np
    from scipy.integrate import quad

    def G(u3, turn):
        sigma = beta * u3 + rho
        w = -sigma / (1 - beta)
        u1 = w - beta
        L = _feynman_log(w) - _feynman_log(u1) - (np.log(u3) + turn) - 1j * np.pi
        return f(L) / (w * u1 * u3)

    h = lambda t: (G(np.exp(t), 0.0) - G(np.exp(t), -2j * np.pi)) * np.exp(t)
    re = quad(lambda t: h(t).real, -reach, reach, limit=800, epsabs=0, epsrel=1e-10)[0]
    im = quad(lambda t: h(t).imag, -reach, reach, limit=800, epsabs=0, epsrel=1e-10)[0]
    return (re + 1j * im) / (1 - beta)


def lorentzian_product_formula(f, eta, nodes: int = 160):
    """eta^2 x (cusp-local model), normalised as `triangle_closed_form`, from
    the Lorentzian derivation:

        eta^2/(2 pi Gamma(1+eta)) int_0^1 d beta beta^(eta-1) (1-beta)^eta
          int_0^oo dy y^(eta-1)/(1+y)
          int_{-pi}^{pi} d phi (1 + e^{-i phi})^eta f(i phi - log(1-beta) - log(1+y)).

    beta is the position along the edge, y the collinear scale, phi the path
    between the lips. `f` must accept numpy arrays.
    """
    import numpy as np
    from scipy.special import gamma

    eta = float(eta)
    t, wt = np.polynomial.legendre.leggauss(nodes)
    u, wu = 0.5 * (t + 1), 0.5 * wt                           # nodes on (0, 1)
    beta = u ** (1 / eta)                                     # d beta beta^(eta-1) = du/eta
    w_beta = wu / eta * (1 - beta) ** eta
    a = -np.log1p(-beta)
    y_small = u ** (1 / eta)                                  # y in (0, 1)
    w_small = wu / eta / (1 + y_small)
    v = u ** (1 / (1 - eta))                                  # y = 1/v in (1, oo)
    w_large = wu / (1 - eta) / (1 + v)
    b = np.concatenate([np.log1p(y_small), np.log1p(1 / v)])
    w_y = np.concatenate([w_small, w_large])
    phi, w_phi = np.pi * t, np.pi * wt * (1 + np.exp(-1j * np.pi * t)) ** eta
    shift = (a[:, None] - b[None, :]).ravel()
    weights = (w_beta[:, None] * w_y[None, :]).ravel()
    total = 0j
    for k in range(0, shift.size, 4096):                      # in blocks, to bound memory
        values = f(1j * phi[None, :] + shift[k:k + 4096, None])
        total += np.sum(weights[k:k + 4096, None] * w_phi[None, :] * values)
    return eta**2 / (2 * np.pi * gamma(1 + eta)) * total


# --------------------------------------------------------------------------
# The full four-point integrand (docs/CUSP_LOCAL_MODEL.md section 9)
# --------------------------------------------------------------------------
def full_invariants(alpha, beta, rho, kappa, T, eps):
    """The four invariants at s = -1, t = -T, each with -i eps.

    kappa is the component of y_perp along the transverse part of the far
    edge, |kappa| <= sqrt(rho). All four are linear in alpha.
    """
    import numpy as np
    c = 2 * np.sqrt(T * (1 + T))
    w = -alpha * beta - rho - 1j * eps
    u1 = w - beta
    u3 = w + alpha
    u4 = -T * (1 - beta) - alpha * (T + beta) - rho + c * kappa - 1j * eps
    return w, u1, u3, u4


def _full_integrand(f, alpha, beta, rho, kappa, T, eps, shift=0.0):
    import numpy as np
    w, u1, u3, u4 = full_invariants(alpha, beta, rho, kappa, T, eps)
    L = (_feynman_log(w) + _feynman_log(u4) - _feynman_log(u1) - _feynman_log(u3)
         - np.log(T) + shift)                                  # z = w u4/(u1 u3 T)
    return f(L) / (w * u1 * u3 * u4)


def full_light_cone_integral(f, beta, rho, kappa, T, eps=1e-2):
    """int over real alpha of the full four-point integrand F(z)/(w u1 u3 u4).

    It vanishes unless -T < beta < 1: the integral along the null direction
    p_1 is supported in the slab between the null hyperplanes through x_4
    and x_3. This holds for every eps > 0.
    """
    import numpy as np
    from scipy.integrate import quad

    zeros = full_alpha_zeros(beta, rho, kappa, T, eps)
    marks = sorted({float(np.real(z)) for z in zeros.values()} | {0.0})
    edges = [-np.inf] + marks + [np.inf]
    g = lambda a: _full_integrand(f, a, beta, rho, kappa, T, eps)
    total = 0j
    for lo, hi in zip(edges[:-1], edges[1:]):
        total += quad(lambda a: g(a).real, lo, hi, limit=600, epsabs=0, epsrel=1e-10)[0]
        total += 1j * quad(lambda a: g(a).imag, lo, hi, limit=600, epsabs=0, epsrel=1e-10)[0]
    return total


def full_alpha_zeros(beta, rho, kappa, T, eps):
    """Zeros in alpha of the four invariants (each is linear in alpha)."""
    import numpy as np
    c = 2 * np.sqrt(T * (1 + T))
    r = rho + 1j * eps
    return {"w": -r / beta, "u1": -1 - r / beta, "u3": r / (1 - beta),
            "u4": (-T * (1 - beta) - r + c * kappa) / (T + beta)}


def full_keyhole_integral(f, beta, rho, kappa, T, eps=1e-2):
    """The same integral after closing the alpha contour.

    For 0 < beta < 1 only the zero of u3 is in the upper half plane; for
    -T < beta < 0 only the zero of u4 is in the lower. The real line is
    deformed onto a keyhole around a vertical cut from that zero: the two
    sides of the cut, which differ by 2 pi i in log z, and a circle around
    the zero. Outside -T < beta < 1 there is nothing to hang on and the
    integral is zero.
    """
    import numpy as np
    from scipy.integrate import quad

    zeros = full_alpha_zeros(beta, rho, kappa, T, eps)
    if 0 < beta < 1:
        name, up = "u3", True
    elif -T < beta < 0:
        name, up = "u4", False
    else:
        return 0j
    a0 = zeros[name]
    r0 = 0.5 * min(abs(a0 - v) for k, v in zeros.items() if k != name)
    g = lambda a: _full_integrand(f, a, beta, rho, kappa, T, eps)
    step = 1j if up else -1j
    tiny = 1e-12

    def sides(t):                       # s = r0 e^t along the cut
        sdist = r0 * np.exp(t)
        a = a0 + step * sdist
        return (g(a + tiny) - g(a - tiny)) * step * sdist

    re = quad(lambda t: sides(t).real, 0, 60, limit=800, epsabs=0, epsrel=1e-10)[0]
    im = quad(lambda t: sides(t).imag, 0, 60, limit=800, epsabs=0, epsrel=1e-10)[0]
    lo, hi = (np.pi / 2 + 1e-9, 5 * np.pi / 2 - 1e-9) if up else (3 * np.pi / 2 - 1e-9, -np.pi / 2 + 1e-9)
    circle = lambda chi: g(a0 + r0 * np.exp(1j * chi)) * 1j * r0 * np.exp(1j * chi)
    cre = quad(lambda c: circle(c).real, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
    cim = quad(lambda c: circle(c).imag, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
    return re + 1j * im + cre + 1j * cim


# --------------------------------------------------------------------------
# The far end of an edge piece (docs/CUSP_LOCAL_MODEL.md section 10)
# --------------------------------------------------------------------------
def far_end_closed_form(p, eta):
    """J(p, eta) = int_0^oo dv v^(-1-p) int d^(2+2 eta) n (v + |n|^2)^(p-1) (v + |n - e|^2)^(p-1),
    |e| = 1, in closed form:

        pi^(1+eta) Gamma(-p) Gamma(1 - p - eta) Gamma(eta)^2 / (Gamma(1-p)^2 Gamma(2 eta)),

    valid for eta - 1 < ... Re p < 0. The far end of an edge piece is
    (1/(2 eta)) * 2 i sin(pi p) * J, whose double pole has coefficient
    sin(pi p)/(pi p), like the near end.
    """
    p, eta = mp.mpmathify(p), mp.mpf(eta)
    return (mp.pi ** (1 + eta) * mp.gamma(-p) * mp.gamma(1 - p - eta) * mp.gamma(eta) ** 2
            / (mp.gamma(1 - p) ** 2 * mp.gamma(2 * eta)))


def far_end_transverse(v, p, eta):
    """int d^(2+2 eta) n (v + |n|^2)^(p-1) (v + |n - e|^2)^(p-1) by direct
    two-dimensional quadrature: n = (kappa, y'), y' in 1 + 2 eta dimensions."""
    import numpy as np
    from scipy.integrate import dblquad
    from scipy.special import gamma

    v, p, eta = float(v), float(p), float(eta)
    sphere = 2 * np.pi ** (eta + 0.5) / gamma(eta + 0.5)       # unit sphere in 1 + 2 eta dimensions

    def g(x, k):                                               # s = |y'| = e^x
        s = np.exp(x)
        return s ** (2 * eta + 1) * (v + k * k + s * s) ** (p - 1) * (v + (k - 1) ** 2 + s * s) ** (p - 1)

    return sphere * sum(dblquad(g, lo, hi, -30, 10, epsabs=0, epsrel=1e-8)[0]
                        for lo, hi in ((-np.inf, 0.0), (0.0, 1.0), (1.0, np.inf)))


def far_end_transverse_feynman(v, p, eta):
    """The same after combining the two factors with one Feynman parameter:
    pi^(1+eta) Gamma(1-2p-eta)/Gamma(1-p)^2 int_0^1 d xi (xi(1-xi))^(-p) (v + xi(1-xi))^(2p-1+eta)."""
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import gamma

    v, p, eta = float(v), float(p), float(eta)
    def g(x):                                                  # xi = e^x, to resolve xi ~ v
        xi = np.exp(x)
        m = xi * (1 - xi)
        return xi * m ** (-p) * (v + m) ** (2 * p - 1 + eta)

    knot = min(np.log(v), np.log(0.5) - 1e-9)
    integral = 2 * sum(quad(g, lo, hi, epsabs=0, epsrel=1e-11, limit=200)[0]
                       for lo, hi in ((knot - 60, knot), (knot, np.log(0.5))) if hi > lo)
    return np.pi ** (1 + eta) * gamma(1 - 2 * p - eta) / gamma(1 - p) ** 2 * integral


def far_end_numerical(p, eta):
    """J(p, eta) by quadrature in v of the transverse integral (real p < 0)."""
    import numpy as np
    from scipy.integrate import quad

    p, eta = float(p), float(eta)
    return sum(quad(lambda t: np.exp(-p * t) * far_end_transverse_feynman(np.exp(t), p, eta), lo, hi,
                    epsabs=0, epsrel=1e-9, limit=200)[0]
               for lo, hi in ((-40 / eta, -5), (-5, 0), (0, 5), (5, 80)))


def far_end_coefficient(p, eta, T=1.0):
    """eta^2 times the far-end leading term of an edge piece, in units in
    which the near end is `triangle_closed_form`:

        (T/(1+T))^eta Gamma(1+eta) Gamma(eta)^2 eta^2 Gamma(1-p-eta)
          / (Gamma(1+p) Gamma(1-p)^2 Gamma(1 + 2 eta)).

    It tends to sin(pi p)/(pi p).
    """
    p, eta, T = mp.mpmathify(p), mp.mpf(eta), mp.mpf(T)
    return ((T / (1 + T)) ** eta * mp.gamma(1 + eta) * eta**2 * mp.gamma(eta) ** 2 * mp.gamma(1 - p - eta)
            / (mp.gamma(1 + p) * mp.gamma(1 - p) ** 2 * mp.gamma(1 + 2 * eta)))


# --------------------------------------------------------------------------
# The edge piece as an explicit integral (docs/CUSP_LOCAL_MODEL.md section 11)
# --------------------------------------------------------------------------
def _edge_inner_terms(lbeta, lg, p, eta):
    """int_0^1 (1-xi)^(-p) xi^p (xi beta + g)^(eta-2) d xi
       = g^(eta-2) Gamma(1+p) Gamma(1-p) 2F1(2-eta, 1+p; 2; -beta/g),
    returned as a list of (log|term|, sign), using the connection formula
    when beta/g > 1 so that nothing overflows."""
    import numpy as np
    from scipy.special import gamma, hyp2f1, rgamma

    a, b, c = 2 - eta, 1 + p, 2.0
    lX = lbeta - lg
    pref = gamma(1 + p) * gamma(1 - p)
    if lX <= 0:
        v = pref * hyp2f1(a, b, c, -np.exp(lX))
        return [((eta - 2) * lg + np.log(abs(v)), np.sign(v))]
    inv = -np.exp(-lX)
    c1 = pref * gamma(b - a) * rgamma(b) * rgamma(c - a) * hyp2f1(a, a - c + 1, a - b + 1, inv)
    c2 = pref * gamma(a - b) * rgamma(a) * rgamma(c - b) * hyp2f1(b, b - c + 1, b - a + 1, inv)
    return [((eta - 2) * lbeta + np.log(abs(c1)), np.sign(c1)),
            ((eta - 1 + p) * lg - (1 + p) * lbeta + np.log(abs(c2)), np.sign(c2))]


def edge_cube_integral(p, eta, T, rel: float = 1e-9, sign: int = 1):
    """J = (T+1)^(-1) int_0^1 d beta b^(eta-p) int_0^1 d xi int_0^1 d mu
             (1-xi)^(-p) xi^p mu^(-p) (1-mu)^(eta-1)
             (xi beta + tau b mu)^(eta-2) ((T+1) beta + T b mu)^p,
    b = 1 - beta, tau = T/(T+1); real p with |p| < 1 - eta and p + eta not an
    integer. The xi integral is done in closed form; beta and mu are
    integrated in logit variables.

    `sign = -1` replaces T + 1 by T - 1 (T > 1): the edge (2,3) piece when
    the diagonal t = +T is timelike."""
    import numpy as np
    from scipy.integrate import quad

    p, eta, T = float(p), float(eta), float(T)
    if sign == -1 and T <= 1:
        raise ValueError("sign = -1 needs T > 1")
    T1 = T + sign
    tau = T / T1
    R = 44 / eta

    def integrand(xb, xm):
        lbeta, lb = -np.logaddexp(0, -xb), -np.logaddexp(0, xb)
        lmu, l1mu = -np.logaddexp(0, -xm), -np.logaddexp(0, xm)
        lg = np.log(tau) + lb + lmu
        lc = np.logaddexp(np.log(T1) + lbeta, np.log(T) + lb + lmu)
        base = (eta - p) * lb - p * lmu + (eta - 1) * l1mu + p * lc + lbeta + lb + lmu + l1mu
        return sum(s * np.exp(base + l) for l, s in _edge_inner_terms(lbeta, lg, p, eta))

    def outer(xb):
        pts = sorted({-R, 0.0, R, float(np.clip(-xb, -R, R)), float(np.clip(xb, -R, R))})
        return sum(quad(lambda xm: integrand(xb, xm), lo, hi, epsabs=0, epsrel=rel, limit=400)[0]
                   for lo, hi in zip(pts[:-1], pts[1:]))

    pts = [-R, -30, -5, 0, 5, 30, R] if R > 30 else [-R, -5, 0, 5, R]
    total = sum(quad(outer, lo, hi, epsabs=0, epsrel=20 * rel, limit=400)[0] for lo, hi in zip(pts[:-1], pts[1:]))
    return total / T1


def edge_piece(p, eta, T, sign: int = 1):
    """The edge piece A(p; T) of the full integral, with s = -1, t = -T,
    prefactor s t and measure d alpha d beta d^(2+2 eta) y:

        A = -2 i pi^(2+eta) Gamma(2-eta) T^(1-p) / (Gamma(1-p)^2 Gamma(1+p)^2) * J.

    With `sign = -1` it is the edge (2,3) piece for t = +T > 1.
    """
    import numpy as np
    from scipy.special import gamma

    return (-2j * np.pi ** (2 + eta) * gamma(2 - eta) * T ** (1 - p)
            / (gamma(1 - p) ** 2 * gamma(1 + p) ** 2) * edge_cube_integral(p, eta, T, sign=sign))


def edge_near_term(p, eta):
    import numpy as np
    from scipy.special import gamma
    return -2j * np.pi ** (2 + eta) * gamma(eta) ** 2 * gamma(1 + p - eta) / (gamma(1 + p) ** 2 * gamma(1 + 2 * eta - p))


def edge_far_term(p, eta, T, sign: int = 1):
    import numpy as np
    from scipy.special import gamma
    return (-2j * np.pi ** (2 + eta) * (T / (T + sign)) ** eta * gamma(eta) ** 2 * gamma(1 - p - eta)
            / (gamma(1 + p) * gamma(1 - p) ** 2 * gamma(1 + 2 * eta)))


def box_from_two_edges(p, eta, T):
    """A(p; T) + T^eta A(-p; 1/T): the second piece is the first with the two
    diagonals exchanged and z -> 1/z."""
    return edge_piece(p, eta, T) + T ** eta * edge_piece(-p, eta, 1 / T)


# --------------------------------------------------------------------------
# Any kinematic region: zeros, half planes and keyholes (section 14)
# --------------------------------------------------------------------------
def general_alpha_zeros(beta, rho, kappa, s, t, eps):
    """Zeros in alpha of the four invariants for diagonals s = +-1 and t, each
    with -i eps, and the coefficient of alpha in each:

        w  = alpha beta s - rho,        u1 = w + beta s,
        u3 = w - alpha s,               u4 = t (1 + alpha - beta) + w + c kappa,
        c = 2 sqrt(t u / s),  u = -s - t.

    A zero lies in the upper half plane exactly when its coefficient is positive.
    """
    import numpy as np
    u = -s - t
    if t * u / s <= 0:
        raise ValueError("exactly one of s, t, u must be positive")
    c = 2 * np.sqrt(t * u / s)
    r = rho + 1j * eps
    coeff = {"w": beta * s, "u1": beta * s, "u3": s * (beta - 1), "u4": t + beta * s}
    const = {"w": -r, "u1": -r + beta * s, "u3": -r, "u4": t * (1 - beta) - r + c * kappa}
    return {k: -const[k] / coeff[k] for k in coeff}, coeff, const


def _general_integrand(f, alpha, coeff, const, s, t):
    import numpy as np
    val = {k: coeff[k] * alpha + const[k] for k in coeff}
    L = (_feynman_log(val["w"]) + _feynman_log(val["u4"]) - _feynman_log(val["u1"]) - _feynman_log(val["u3"])
         + _feynman_log(complex(s)) - _feynman_log(complex(t)))
    return f(L) / (val["w"] * val["u1"] * val["u3"] * val["u4"])


def general_light_cone_integral(f, beta, rho, kappa, s, t, eps=2e-2):
    """int over real alpha of F(z)/(x01^2 x02^2 x03^2 x04^2) in any region."""
    import numpy as np
    from scipy.integrate import quad

    zeros, coeff, const = general_alpha_zeros(beta, rho, kappa, s, t, eps)
    g = lambda a: _general_integrand(f, a, coeff, const, s, t)
    marks = sorted({float(np.real(z)) for z in zeros.values()} | {0.0})
    edges = [-np.inf] + marks + [np.inf]
    total = 0j
    for lo, hi in zip(edges[:-1], edges[1:]):
        total += quad(lambda a: g(a).real, lo, hi, limit=600, epsabs=0, epsrel=1e-10)[0]
        total += 1j * quad(lambda a: g(a).imag, lo, hi, limit=600, epsabs=0, epsrel=1e-10)[0]
    return total


def general_keyhole_sum(f, beta, rho, kappa, s, t, eps=2e-2):
    """The same integral as a sum of keyholes around the zeros in the half
    plane that holds fewer of them. Returns (value, names of those zeros,
    'U'/'L' pattern for w, u1, u3, u4)."""
    import numpy as np
    from scipy.integrate import quad

    zeros, coeff, const = general_alpha_zeros(beta, rho, kappa, s, t, eps)
    pattern = "".join("U" if zeros[k].imag > 0 else "L" for k in ("w", "u1", "u3", "u4"))
    upper = [k for k in zeros if zeros[k].imag > 0]
    lower = [k for k in zeros if zeros[k].imag < 0]
    if not upper or not lower:
        return 0j, [], pattern
    group, up = (upper, True) if len(upper) <= len(lower) else (lower, False)
    g = lambda a: _general_integrand(f, a, coeff, const, s, t)
    step, tiny = (1j if up else -1j), 1e-12
    total = 0j
    for name in group:
        a0 = zeros[name]
        r0 = 0.4 * min(abs(a0 - v) for k, v in zeros.items() if k != name)

        def sides(tt, a0=a0, r0=r0):
            dist = r0 * np.exp(tt)
            a = a0 + step * dist
            return (g(a + tiny) - g(a - tiny)) * step * dist

        total += quad(lambda tt: sides(tt).real, 0, 60, limit=800, epsabs=0, epsrel=1e-10)[0]
        total += 1j * quad(lambda tt: sides(tt).imag, 0, 60, limit=800, epsabs=0, epsrel=1e-10)[0]
        lo, hi = (np.pi / 2 + 1e-9, 5 * np.pi / 2 - 1e-9) if up else (3 * np.pi / 2 - 1e-9, -np.pi / 2 + 1e-9)
        circle = lambda chi, a0=a0, r0=r0: g(a0 + r0 * np.exp(1j * chi)) * 1j * r0 * np.exp(1j * chi)
        total += quad(lambda ch: circle(ch).real, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
        total += 1j * quad(lambda ch: circle(ch).imag, lo, hi, limit=400, epsabs=0, epsrel=1e-10)[0]
    return total, group, pattern


# --------------------------------------------------------------------------
# The edge jump as the residue of the collinear pole (section 12)
# --------------------------------------------------------------------------
def collinear_one_sided(f, beta, side, psi: float = 0.7):
    """k_side = int_0^oo dr f(L)/(w u3 u1) in the scaling limit alpha -> 0 with
    rho = |alpha| r, in four dimensions, for the cusp-local model at fixed beta:

        w = -(side*beta + r),  u3 = side*(1 - beta) - r,  u1 = -beta,  z = w/(beta u3),

    side = +1 or -1 being the sign of alpha. The -i0 is a rotation of r by e^{i psi}.
    Then int d alpha K(alpha) has the collinear pole (k_+ + k_-)/eta.
    """
    import numpy as np
    from scipy.integrate import quad

    def g(r):
        rr = r * np.exp(1j * psi)
        w = -(side * beta + rr)
        u3 = side * (1 - beta) - rr
        L = _feynman_log(w) - _feynman_log(u3) - np.log(beta)
        return f(L) / (w * u3 * (-beta)) * np.exp(1j * psi)

    marks = sorted({0.0, beta, 1 - beta, 1.0})
    edges = marks + [np.inf]
    total = 0j
    for lo, hi in zip(edges[:-1], edges[1:]):
        total += quad(lambda r: g(r).real, lo, hi, limit=400, epsabs=0, epsrel=1e-11)[0]
        total += 1j * quad(lambda r: g(r).imag, lo, hi, limit=400, epsabs=0, epsrel=1e-11)[0]
    return total


def collinear_residue(f, beta):
    """k_+ + k_-: the coefficient of 1/eta in the collinear pole at fixed beta,
    computed with the transverse integral first. For the weight F it should be
    (2 pi i / beta) M[F(z/(1 - beta))], the density of section 8."""
    return collinear_one_sided(f, beta, 1) + collinear_one_sided(f, beta, -1)


# --------------------------------------------------------------------------
# The edge piece as a one-dimensional integral (section 11.5, section 14)
# --------------------------------------------------------------------------
def edge_piece_one_dim(p, eta, rho, side: int = 0, dps: int = 25):
    """The edge piece for F = z^p as a single integral, with rho = T/(T+1):

        A = -2 i pi^(2+eta) Gamma(2-eta) B(eta,eta) / (Gamma(1-p) Gamma(1+p))
            * rho^eta * int_0^1 dy (1-y)^(eta-1) (1-(1-rho) y)^(-eta) 2F1(eta, 1+p; 2; y).

    For rho > 0 this equals the cube integral `edge_piece`. For rho < 0 it
    is the boundary value at rho + i0 (side = +1) or rho - i0 (side = -1):
    1 - (1-rho) y changes sign at y = 1/(1-rho), an integrable singularity
    of exponent -eta. Returns an mpmath complex number.
    """
    import mpmath as mp

    with mp.workdps(dps):
        p, eta, rho = mp.mpf(p), mp.mpf(eta), mp.mpf(rho)
        if rho < 0 and side not in (1, -1):
            raise ValueError("for rho < 0 give side = +1 or -1")
        if rho == 0:
            raise ValueError("rho must not vanish")

        def factor(y):
            v = 1 - (1 - rho) * y
            # 1 + rho X with X = y/(1-y) is v/(1-y); below zero its phase is e^{+-i pi}
            return (-v) ** (-eta) * mp.expjpi(-side * eta) if v < 0 else v ** (-eta)

        def integrand(u):
            y = 1 - u ** (1 / eta)                      # (1-y)^(eta-1) dy = du / eta
            return factor(y) * mp.hyp2f1(eta, 1 + p, 2, y) / eta

        points = [0, mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.75"), 1]
        if rho < 0:
            points = sorted(set(points + [(1 - 1 / (1 - rho)) ** eta]))
        integral = mp.quad(integrand, points)
        rho_eta = rho ** eta if rho > 0 else (-rho) ** eta * mp.expjpi(side * eta)
        prefactor = -2j * mp.pi ** (2 + eta) * mp.gamma(2 - eta) / (mp.gamma(1 - p) * mp.gamma(1 + p))
        return prefactor * mp.beta(eta, eta) * rho_eta * integral


def timelike_two_edges(p, eta, modulus, sigma: int = -1, dps: int = 25):
    """A(p;T) + T^eta A(-p;1/T) continued to T = modulus * e^{i sigma pi}.

    One of rho = T/(T+1) and 1 - rho = 1/(T+1) is larger than 1 and the
    other is negative; the negative one is approached from the side fixed
    by sigma. Returns (first piece, second piece including T^eta).
    """
    import mpmath as mp

    with mp.workdps(dps):
        if sigma not in (1, -1):
            raise ValueError("sigma must be +1 or -1")
        modulus = mp.mpf(modulus)
        if modulus <= 0 or modulus == 1:
            raise ValueError("the modulus must be positive and different from 1")
        T = -modulus
        rho, other = T / (1 + T), 1 / (1 + T)
        first = edge_piece_one_dim(p, eta, rho, side=sigma if rho < 0 else 0, dps=dps)
        second = edge_piece_one_dim(-mp.mpf(p), eta, other, side=-sigma if other < 0 else 0, dps=dps)
        return first, modulus ** mp.mpf(eta) * mp.expjpi(sigma * mp.mpf(eta)) * second


# --------------------------------------------------------------------------
# A timelike diagonal, computed directly (section 14)
# --------------------------------------------------------------------------
def _hyp_above_cut(a, b, c, X, delta):
    """2F1(a, b; c; X + i0) for X > 1, with delta = X - 1 given separately."""
    import numpy as np
    from scipy.special import gamma, hyp2f1

    e = c - a - b
    if delta < 0.5 and abs(e - round(e)) > 1e-6:
        return (gamma(c) * gamma(e) / (gamma(c - a) * gamma(c - b)) * hyp2f1(a, b, a + b - c + 1, -delta)
                + delta ** e * np.exp(-1j * np.pi * e) * gamma(c) * gamma(-e) / (gamma(a) * gamma(b))
                * hyp2f1(c - a, c - b, e + 1, -delta))
    z = 1.0 / X
    return (gamma(c) * gamma(b - a) / (gamma(b) * gamma(c - a)) * X ** (-a) * np.exp(1j * np.pi * a)
            * hyp2f1(a, a - c + 1, a - b + 1, z)
            + gamma(c) * gamma(a - b) / (gamma(a) * gamma(c - b)) * X ** (-b) * np.exp(1j * np.pi * b)
            * hyp2f1(b, b - c + 1, b - a + 1, z))


def timelike_direct(p, eta, T, rel=1e-6):
    """s t times the iterated integral of z^p / (x01^2 x02^2 x03^2 x04^2) for
    s = -1, t = T > 1, computed in those kinematics.

    Region I, 0 < beta < 1: the alpha contour is closed around the compact
    cut between the zeros of x02^2 and x01^2. Region II, 1 < beta < T: a
    keyhole around the zero of x04^2. The transverse integral is done with
    Feynman parameters and the remaining ones in closed form, down to a
    single integral for region I and a double integral for region II. No
    contour deformation in x0, no Mellin-Barnes, no continuation in T.
    Returns (region I, region II), each multiplied by s t = -T.
    Parameters with p = eta or p + eta = 1 are degenerate for the
    hypergeometric connection formulas used; shift them slightly.
    """
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import gamma

    if not T > 1:
        raise ValueError("this region needs T > 1")
    if abs(p) < 1e-9:
        raise ValueError("p = 0 is excluded by the formulas used; shift it slightly")

    def cquad(f, lo, hi):
        return (quad(lambda x: f(x).real, lo, hi, limit=60, epsabs=0, epsrel=rel)[0]
                + 1j * quad(lambda x: f(x).imag, lo, hi, limit=60, epsabs=0, epsrel=rel)[0])

    front = np.pi ** (1 + eta) * gamma(2 - eta) * T ** (-p) * np.sin(np.pi * p) / p

    def h(m):
        return quad(lambda v: (m - v ** (1 / eta)) ** (eta - 1) / eta, 0, 1, limit=200)[0]

    def one(u):
        xi = min(u ** (1 / eta), 1 - 1e-15)
        m = xi + (1 - xi) * T
        ratio = (1 - xi) * T / m
        F = _hyp_above_cut(2 - eta, 1 + p, 2, 1 / ratio, xi / ((1 - xi) * T))
        return u ** (1 / eta - 1) / eta * xi ** p * (1 - xi) ** (-p) * m ** (-1 - eta) * ratio ** (eta - 2) * F * h(m)

    region_one = -2j * front * np.exp(1j * np.pi * (eta - p)) * cquad(one, 0, 1)

    def two(beta):
        beta = min(max(beta, 1 + 1e-14), T - 1e-14)

        def inner(u):
            chi = min(max(u ** (1 / eta), 1e-30), 1 - 1e-15)
            lam = (T - 1 + chi) / (T - beta)
            g = min((1 - chi) / lam - beta, -1e-30)
            F = _hyp_above_cut(2 - eta, 1 + p, 2, beta / (-g), ((1 - chi) / lam) / (-g))
            return ((1 - chi) ** p * lam ** (-1 - eta) * (chi * beta + (1 - chi) * (beta - 1)) ** (-p)
                    * (-g) ** (eta - 2) * F / eta)

        return (T - beta) ** (p - 1) * cquad(inner, 0, 1)

    region_two = -2j * front * np.exp(1j * np.pi * (eta - p)) * cquad(two, 1, T)
    return -T * region_one, -T * region_two


def timelike_keyhole_piece(p, eta, T, rel=1e-4):
    """s t times the (4,1) piece for s = -1, t = T > 1: the keyhole around the
    zero of x04^2 laid along x04^2 < 0, over the whole slab 0 < beta < T, with
    every other invariant at its Feynman value.

    After two Feynman parameters, the transverse integral, a Beta function
    and an Euler integral it is a double integral over beta and chi with
    omega = 1 - beta - chi:

        (T-beta)^(p+eta) chi^(eta-1) (1-chi)^p (T-1+chi)^(-1-eta)
            (omega + i0)^(-p) (g + i0)^(eta-2) 2F1(2-eta, 1+p; 2; -beta/(g + i0)),
        g = T omega / (T - 1 + chi).

    The substitution beta' = beta/T, mu = (T-1)(1-chi)/(T-1+chi) turns it
    into the second edge piece at rho = -1/(T-1) - i0 (section 14). The
    quadrature here is crude; expect three to four digits.
    """
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import gamma, hyp2f1

    if not T > 1:
        raise ValueError("this region needs T > 1")

    def cquad(f, lo, hi, points=None):
        kw = dict(limit=80, epsabs=0, epsrel=rel)
        if points:
            kw["points"] = points
        return quad(lambda x: f(x).real, lo, hi, **kw)[0] + 1j * quad(lambda x: f(x).imag, lo, hi, **kw)[0]

    def inner(beta):
        def f(chi, one_minus_chi, omega):
            """The integrand, with 1 - chi and omega = 1 - beta - chi supplied
            separately so that neither is formed by cancellation."""
            common = chi ** (eta - 1) * one_minus_chi ** p * ((T - 1 + chi) / (T - beta)) ** (-1 - eta)
            if omega > 0:
                g = T * omega / (T - 1 + chi)
                X = beta / g
                return common * omega ** (-p) * g ** (eta - 2) * (1 + X) ** (-1 - p) * hyp2f1(eta, 1 + p, 2, X / (1 + X))
            g = T * omega / (T - 1 + chi)
            F = _hyp_above_cut(2 - eta, 1 + p, 2, beta / (-g), (T - beta) * one_minus_chi / (T * (-omega)))
            return (common * (-omega) ** (-p) * np.exp(-1j * np.pi * p) * (-g) ** (eta - 2)
                    * np.exp(1j * np.pi * eta) * F)

        cuts = [0.0] + ([1 - beta] if 0 < 1 - beta < 1 else []) + [1.0]
        total, k = 0j, 1 / eta
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            half = (hi - lo) / 2

            def left(u, lo=lo, half=half):
                d = half * max(u, 1e-12) ** k                      # chi = lo + d
                omega = -d if lo > 0 else 1 - beta - d
                return f(lo + d, 1 - lo - d, omega) * half * k * u ** (k - 1)

            def right(u, hi=hi, half=half):
                d = half * max(u, 1e-12) ** k                      # chi = hi - d
                omega = d if hi < 1 else -beta + d
                one_minus = d if hi == 1 else 1 - hi + d
                return f(hi - d, one_minus, omega) * half * k * u ** (k - 1)

            total += cquad(left, 0, 1) + cquad(right, 0, 1)
        return total

    integral = cquad(lambda beta: (T - min(beta, T - 1e-13)) ** (p - 1) * inner(min(max(beta, 1e-13), T - 1e-13)),
                     0, T, points=[1.0])
    return T * (2j * np.sin(np.pi * p) / p) * np.pi ** (1 + eta) * gamma(2 - eta) * T ** (-p) * integral


# --------------------------------------------------------------------------
# A timelike diagonal with the displacement kept (section 14, "removed last")
# --------------------------------------------------------------------------
def _timelike_tau_integral(p, T, beta, Y0, Y1):
    """Continuation in p of int_0^oo dtau tau^(p-1) D1^(-1-p) (D1-beta)^(p-1) D3^(-1-p),
    D1 = A Y0 - beta tau, D3 = A' Y1 + (1-beta) tau, for complex Y0, Y1 in the
    upper half plane and -1 < Re p < 1, p != 0 (the finite part at tau = 0)."""
    import numpy as np
    from scipy.integrate import quad

    A, Ap = T / (T - beta), (T - 1) / (T - beta)

    def cq(f, a, b):
        kw = dict(limit=200, epsabs=1e-12, epsrel=1e-9)
        return quad(lambda x: f(x).real, a, b, **kw)[0] + 1j * quad(lambda x: f(x).imag, a, b, **kw)[0]

    def h(tau):
        D1 = A * Y0 - beta * tau
        return D1 ** (-1 - p) * (D1 - beta) ** (p - 1) * (Ap * Y1 + (1 - beta) * tau) ** (-1 - p)

    h0, k = h(0.0), 1 + p
    marks = [abs(A * Y0 / beta), abs((A * Y0 - beta) / beta)] + ([abs(Ap * Y1 / (beta - 1))] if beta > 1 else [])

    def inner(x):                                   # x = tau^(1+p) on (0, 1)
        tau = x ** (1 / k)
        return (h(tau) - h0) / tau / k if tau > 1e-14 else 0j

    pts = sorted({0.0, 1.0} | {m ** k for m in marks if m < 1})
    total = sum(cq(inner, lo, hi) for lo, hi in zip(pts[:-1], pts[1:])) + h0 / p
    top = max(max(marks), 1.0)
    pts = sorted({1.0, 2 * top, 8 * top} | {m for m in marks if m > 1})
    outer = lambda t: t ** (p - 1) * h(t)
    return total + sum(cq(outer, lo, hi) for lo, hi in zip(pts[:-1], pts[1:])) + cq(outer, pts[-1], np.inf)


def timelike_alpha_split(p, T, beta, y, eps):
    """At displacement eps > 0 and fixed (beta, y), s = -1, t = T > 1: the
    integral over real alpha of z^p / (x01^2 x02^2 x03^2 x04^2), and the two
    keyholes it splits into. The cut from the zero of x03^2 runs to the right
    (x03^2 - i eps > 0) and the cut from the zero of x04^2 to the left
    (x04^2 - i eps < 0); they lie at different heights and never meet.

    y = (component along the direction of kappa, modulus of the rest).
    Returns (real line, K3, K4, K4 from the tau formula with Y -> Y + i eps);
    K3 is zero for beta > 1.
    """
    import numpy as np
    from scipy.integrate import quad

    if not (T > 1 and 0 < beta < T and beta != 1 and eps > 0):
        raise ValueError("needs T > 1, 0 < beta < T, beta != 1, eps > 0")
    c = 2 * np.sqrt(T * (T - 1))
    rho = y[0] ** 2 + y[1] ** 2 + 1j * eps
    coeff = [-beta, -beta, 1 - beta, T - beta]                            # w, u1, u3, u4
    const = [-rho, -rho - beta, -rho, T * (1 - beta) - rho + c * y[0]]
    zeros = [-b / a for a, b in zip(coeff, const)]
    powers = [p - 1, -1 - p, -1 - p, p - 1]

    def f(alpha):
        total = 0j
        for j in range(4):
            x = alpha - zeros[j]
            if coeff[j] > 0 and j == 3:
                lg = np.log(x)                                             # cut to the left
            elif coeff[j] > 0:
                lg = np.log(-x) - 1j * np.pi                               # cut to the right
            else:
                lg = np.log(-1j * x) - 1j * np.pi / 2                      # zero below; cut downwards
            total += powers[j] * (lg + np.log(abs(coeff[j])))
        return np.exp(-1j * np.pi * p) * T ** (-p) * np.exp(total)

    def cq(g, a, b):
        kw = dict(limit=400, epsabs=1e-12, epsrel=1e-10)
        return quad(lambda x: g(x).real, a, b, **kw)[0] + 1j * quad(lambda x: g(x).imag, a, b, **kw)[0]

    marks = sorted(z.real for z in zeros)
    edges = [-np.inf] + marks + [np.inf]
    line = sum(cq(f, lo, hi) for lo, hi in zip(edges[:-1], edges[1:]))

    def keyhole(j, direction):
        a0 = zeros[j]
        r0 = 0.3 * min(abs(a0 - zeros[k]) for k in range(4) if k != j)
        tiny = 1e-11

        def lips(t):
            r = r0 * np.exp(t)
            a = a0 + direction * r
            return (f(a - 1j * tiny) - f(a + 1j * tiny)) * r

        start = 0.0 if direction > 0 else np.pi
        circle = lambda th: f(a0 + r0 * np.exp(1j * th)) * 1j * r0 * np.exp(1j * th)
        return cq(lips, 0, 50) + cq(circle, start + 1e-9, start + 2 * np.pi - 1e-9)

    K3 = keyhole(2, +1) if beta < 1 else 0j
    K4 = keyhole(3, -1)
    y0 = beta * np.sqrt((T - 1) / T)
    y1 = -(1 - beta) * np.sqrt(T / (T - 1))
    Y0 = (y[0] - y0) ** 2 + y[1] ** 2 + 1j * eps
    Y1 = (y[0] - y1) ** 2 + y[1] ** 2 + 1j * eps
    formula = (-2j * np.sin(np.pi * p) * T ** (-p) * (T - beta) ** (p - 1)
               * _timelike_tau_integral(p, T, beta, Y0, Y1))
    return line, K3, K4, formula


def timelike_keyhole_transverse(p, eta, T, beta, eps, n: int = 100):
    """int d^(2+2 eta) y of the (4,1) keyhole at displacement eps > 0 and fixed
    beta, by a tensor Gauss-Legendre rule in y with the tau integral adaptive.
    Works for both signs of p (finite part at tau = 0). Three to four digits
    at n = 100 for eps of order 0.2; the cost grows like n^2."""
    import numpy as np
    from scipy.special import gamma

    if not (T > 1 and 0 < beta < T and beta != 1 and eps > 0):
        raise ValueError("needs T > 1, 0 < beta < T, beta != 1, eps > 0")
    y0 = beta * np.sqrt((T - 1) / T)
    y1 = -(1 - beta) * np.sqrt(T / (T - 1))
    x, w = np.polynomial.legendre.leggauss(n)
    reach = 4.5
    dim = 1 + 2 * eta
    sphere = 2 * np.pi ** (dim / 2) / gamma(dim / 2)
    total = 0j
    for a, wa in zip(reach * x, reach * w):
        ypar = np.sinh(a)
        for c, wc in zip((x + 1) / 2, w / 2):
            r = c * c / (1 - c)
            if r > 150:
                continue
            Y0 = (ypar - y0) ** 2 + r * r + 1j * eps
            Y1 = (ypar - y1) ** 2 + r * r + 1j * eps
            jac = np.cosh(a) * (2 * c - c * c) / (1 - c) ** 2
            total += wa * wc * jac * sphere * r ** (2 * eta) * _timelike_tau_integral(p, T, beta, Y0, Y1)
    return -2j * np.sin(np.pi * p) * T ** (-p) * (T - beta) ** (p - 1) * total


def timelike_second_piece_at_beta(p, eta, T, beta, eps=0.0, dps: int = 18):
    """The same quantity in closed form: after two Feynman parameters, the
    transverse integral, the tau integral and Euler's integral,

        -(2 i sin(pi p)/p) pi^(1+eta) Gamma(2-eta) T^(-p) (T-beta)^(p+eta)
          int_0^1 dchi chi^(eta-1) (1-chi)^p (T-1+chi)^(-1-eta) (omega+i0)^(-p)
              g^(eta-2) 2F1(2-eta, 1+p; 2; -beta/g),

    omega = 1 - beta - chi, g = T omega/(T-1+chi) + i eps (T-1+chi)/((T-beta) chi).
    The displacement enters only through the imaginary part of g. eps = 0
    gives the boundary value g + i0. Returns an mpmath complex number.
    """
    import mpmath as mp

    with mp.workdps(dps):
        p, eta, T, beta = mp.mpf(p), mp.mpf(eta), mp.mpf(T), mp.mpf(beta)
        eps = mp.mpf(eps) if eps > 0 else mp.mpf(10) ** (-dps + 4)
        if not (T > 1 and 0 < beta < T and beta != 1):
            raise ValueError("needs T > 1, 0 < beta < T, beta != 1")

        def f(chi):
            om = 1 - beta - chi
            g = T * om / (T - 1 + chi) + 1j * eps * (T - 1 + chi) / ((T - beta) * chi)
            om_p = om ** (-p) if om > 0 else (-om) ** (-p) * mp.expjpi(-p)
            return (chi ** (eta - 1) * (1 - chi) ** p * (T - 1 + chi) ** (-1 - eta) * om_p
                    * g ** (eta - 2) * mp.hyp2f1(2 - eta, 1 + p, 2, -beta / g))

        pts = [mp.mpf(0), mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, mp.mpf(1)]
        if 0 < 1 - beta < 1:
            pts = sorted(set(pts + [1 - beta]))
        integral = mp.quad(f, pts)
        return (-(2j * mp.sin(mp.pi * p) / p) * mp.pi ** (1 + eta) * mp.gamma(2 - eta)
                * T ** (-p) * (T - beta) ** (p + eta) * integral)
