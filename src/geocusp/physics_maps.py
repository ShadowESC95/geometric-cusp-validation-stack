from __future__ import annotations

import numpy as np


def gamma_cusp_weak(g: float) -> float:
    """Planar N=4 SYM weak-coupling benchmark through g^6 in the package convention."""
    return 4*g**2 - (4*np.pi**2/3)*g**4 + (44*np.pi**4/45)*g**6


def local_mass_from_vev(vev, coupling=1.0, m0=0.0):
    """User-framework phenomenological map m_eff(x)=m0+coupling*<phi(x)>.

    This is a declared *ansatz* for downstream continuum experiments, not a theorem.
    """
    return m0 + coupling*np.asarray(vev)


def local_coupling_from_vev(vev, alpha=1.0, beta=0.0):
    """User-framework ansatz g_eff(x)=alpha*f(<phi>)+beta with linear f for baseline tests."""
    return beta + alpha*np.asarray(vev)


def sigma_ir(x, g, log_ratio, width=1.0, amplitude=1.0):
    """A reproducible 1D source profile for continuum/entropy probes."""
    gamma=gamma_cusp_weak(g)
    x=np.asarray(x)
    return amplitude*gamma*log_ratio*np.exp(-(x/width)**2)
