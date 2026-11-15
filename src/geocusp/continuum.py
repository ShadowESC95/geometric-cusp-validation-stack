from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from scipy.sparse import diags, eye
from scipy.sparse.linalg import spsolve

@dataclass
class ContinuumRun:
    t: np.ndarray
    S: np.ndarray
    sigma: np.ndarray
    x: np.ndarray


def gamma_cusp_weak(g: float) -> float:
    return 4*g*g - (4*np.pi**2/3)*g**4 + (44*np.pi**4/45)*g**6


def run_entropy_transport(n=96, domain=4.0, D=0.25, tau_max=1.5, nt=120, m0=0.7,
                          vev_amp=0.8, gamma_source=1.0, log_ratio=3.0, seed=20260918):
    """Reduced boundary/entropy PDE probe.

    This is not asserted to be the user's final entropy law. It is a numerical
    falsification sandbox for a proposed source sigma_IR = Gamma_cusp * log(w/lambda).
    """
    x=np.linspace(-domain,domain,n)
    dx=x[1]-x[0]
    t=np.linspace(0,tau_max,nt)
    vev=1.0+vev_amp*np.exp(-(x/(0.7*domain))**2)
    m_eff=m0*np.abs(vev)
    sigma=gamma_source*log_ratio*(1+0.15*np.cos(2*np.pi*x/(2*domain)))
    L = diags([np.ones(n-1), -2*np.ones(n), np.ones(n-1)],[-1,0,1],shape=(n,n))/dx**2
    A = eye(n)-0.5*(t[1]-t[0])*D*L
    B = eye(n)+0.5*(t[1]-t[0])*D*L
    # Neumann-like zero flux by reflecting end points.
    A=A.tolil(); B=B.tolil()
    A[0,0]=1+(t[1]-t[0])*D/dx**2; A[0,1]=-(t[1]-t[0])*D/dx**2
    A[-1,-1]=1+(t[1]-t[0])*D/dx**2; A[-1,-2]=-(t[1]-t[0])*D/dx**2
    B[0,0]=1-(t[1]-t[0])*D/dx**2; B[0,1]=(t[1]-t[0])*D/dx**2
    B[-1,-1]=1-(t[1]-t[0])*D/dx**2; B[-1,-2]=(t[1]-t[0])*D/dx**2
    A=A.tocsr(); B=B.tocsr()
    S=np.zeros((nt,n)); S[0]=0.0
    for j in range(1,nt):
        rhs=B@S[j-1]+(t[1]-t[0])*sigma
        S[j]=spsolve(A,rhs)
    return ContinuumRun(t,S,sigma,x)
