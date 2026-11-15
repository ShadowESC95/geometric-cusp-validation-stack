from __future__ import annotations
import numpy as np


def exceptional_ratio_curve(radius=1.0, x=1.3, z0=0.9, points=4096, phase0=-np.pi):
    phi=np.linspace(phase0,phase0+2*np.pi,points,endpoint=False)
    eta=radius*np.exp(1j*phi)
    z=eta*z0/x
    return eta,z,phi


def winding(points, point=0j):
    from .contours import winding_number

    return winding_number(points, point)


def pullback_dlog_check(x=1.3,z0=.9, radius=1.0, points=1024):
    eta,z,_=exceptional_ratio_curve(radius,x,z0,points)
    dlog_eta=np.mean(1.0/eta[1:]*np.diff(eta))
    dlog_z=np.mean(1.0/z[1:]*np.diff(z))
    return float(abs(dlog_eta-dlog_z)), float(winding(z))
