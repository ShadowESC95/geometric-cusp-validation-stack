from __future__ import annotations

import numpy as np


def Fbar0(z):
    """Normalized four-point finite function quoted in the source paper's Eq. (2.30)."""
    return np.ones_like(np.asarray(z, dtype=complex)) * (-1.0)


def Fbar1(z):
    z=np.asarray(z,dtype=complex)
    return np.log(z)**2 + np.pi**2


def Fbar2(z, dps: int = 20):
    """Two-loop finite function, arXiv:2112.06956 Eq. (2.23) ("Fweak2").

    Transcribed from the paper's TeX source. It is an input from the
    literature, not something this package derives. Returns an mpmath number.
    """
    import mpmath as mp

    with mp.workdps(dps):
        z = mp.mpmathify(z)
        pi, li = mp.pi, mp.polylog
        log = mp.log(z)
        a, b = 1/(z + 1), z/(z + 1)
        s2 = li(2, a) + li(2, b)
        return (-log**4/2
                + log**2*(mp.mpf(2)/3*s2 - 19*pi**2/9)
                + log*(4*li(3, a) - 4*li(3, b))
                + mp.mpf(2)/3*(s2 - pi**2/6)**2
                + mp.mpf(8)/3*pi**2*(s2 - pi**2/6)
                + 8*li(4, a) + 8*li(4, b)
                - 23*pi**4/18)


def circle_mean_mp(func, dps: int = 20):
    """Unit-circle mean on the principal sheet, by adaptive mpmath quadrature.

    Handles the integrable logarithmic endpoint singularities at z=-1 that
    the equally spaced `contour_mean` does not resolve.
    """
    import mpmath as mp

    with mp.workdps(dps):
        nodes = mp.linspace(-mp.pi, mp.pi, 9)
        total = mp.quad(lambda phi: func(mp.exp(1j*phi)), nodes)
        return complex(total/(2*mp.pi))


def insertion_F(g, z, order=0):
    """Wilson-loop insertion normalization inferred from the source's Eq. (2.34): F = -g^2 Fbar.

    Keep this separate from Fbar because the source uses both closely related conventions.
    """
    base=Fbar0(z) if order==0 else Fbar1(z)
    return -g**2*base


def contour_mean(func, radius=1.0, n=8192, phase_offset=1e-7):
    phi=np.linspace(-np.pi+phase_offset,np.pi+phase_offset,n,endpoint=False)
    z=radius*np.exp(1j*phi)
    vals=np.asarray(func(z),dtype=complex)
    return np.mean(vals), phi, vals


def contour_I_Fbar(g, func, radius=1.0, n=8192):
    """Literal implementation of the source Eq. (6.2) form if `func` is the normalized Fbar."""
    avg,phi,vals=contour_mean(func,radius,n)
    return -g**2*avg, phi, vals


def contour_I_insertion(func, radius=1.0, n=8192):
    """Equivalent contour action on the insertion-normalized object with the explicit g^2 removed."""
    avg,phi,vals=contour_mean(func,radius,n)
    return avg, phi, vals


def cusp_coefficient_from_mean(loop_order: int, mean: complex) -> complex:
    """Convert the unit-circle mean of Fbar^(L-1) to the g^(2L) coefficient.

    Dixon et al., arXiv:2605.28926, Eqs. (2.30)–(2.31), (6.2)–(6.5).
    The contour contributes -g^2 times the mean. The coupling derivative
    multiplies the g^(2L) coefficient by 2L.

    The overall constant is a convention: arXiv:2112.06956 prints
    g d/dg Gamma = -2 I[F], which with its own F^(0) = -1 would give
    Gamma^(1) = 1, not the 4 quoted in the same paper. The factor here (-8 I)
    is the one that returns Gamma^(1) = 4. One loop therefore calibrates the
    rule; two and three loops are checks that could fail.
    """
    if loop_order < 1:
        raise ValueError("loop_order must be positive")
    return -4 * mean / loop_order


def weak_cusp(g):
    return 4*g**2 - (4*np.pi**2/3)*g**4 + (44*np.pi**4/45)*g**6


def weak_coefficients():
    return {1:4.0,2:-4*np.pi**2/3,3:44*np.pi**4/45}
