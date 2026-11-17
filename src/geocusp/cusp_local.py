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
