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
