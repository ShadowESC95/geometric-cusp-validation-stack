"""QuTiP open-system entropy probe.

Tests whether a simple driven-dissipative boundary degree of freedom produces a
monotone von Neumann entropy under chosen Lindblad dynamics. This is a probe, not
an identification of S_OB with von Neumann entropy.
"""
from __future__ import annotations
import numpy as np


def run(T=10.0, nt=1000, kappa=0.4, drive=0.8):
    try:
        import qutip as qt
    except Exception as exc:
        raise RuntimeError("QuTiP unavailable") from exc
    sx,sz,sm=qt.sigmax(),qt.sigmaz(),qt.sigmam()
    H=0.5*drive*sx
    rho0=qt.basis(2,0)*qt.basis(2,0).dag()
    t=np.linspace(0,T,nt)
    result=qt.mesolve(H,rho0,t,[np.sqrt(kappa)*sm],[])
    S=np.asarray([qt.entropy_vn(rho,base=np.e) for rho in result.states])
    return t,S,result.states

if __name__=='__main__':
    t,S,_=run(); print('QuTiP entropy range:',float(S.min()),float(S.max()))
