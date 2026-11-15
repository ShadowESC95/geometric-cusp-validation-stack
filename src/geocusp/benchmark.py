from __future__ import annotations

import numpy as np


def Fbar0(z):
    """Normalized four-point finite function quoted in the source paper's Eq. (2.30)."""
    return np.ones_like(np.asarray(z, dtype=complex)) * (-1.0)


def Fbar1(z):
    z=np.asarray(z,dtype=complex)
    return np.log(z)**2 + np.pi**2


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
    return -avg, phi, vals


def weak_cusp(g):
    return 4*g**2 - (4*np.pi**2/3)*g**4 + (44*np.pi**4/45)*g**6


def weak_coefficients():
    return {1:4.0,2:-4*np.pi**2/3,3:44*np.pi**4/45}
