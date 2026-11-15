from __future__ import annotations
import numpy as np


def qubit_entropy(rho: np.ndarray) -> float:
    ev=np.linalg.eigvalsh((rho+rho.conj().T)/2).real
    ev=np.clip(ev,1e-15,1)
    return float(-np.sum(ev*np.log(ev)))


def fallback_lindblad(H, rho0, collapse_ops, tlist, dt=None):
    """Simple RK4 fallback for d rho/dt = -i[H,rho] + sum(L rho L† - 1/2{L†L,rho}).
    For serious runs use QuTiP's mesolve/mcsolve/HEOM backends.
    """
    H=np.asarray(H,dtype=complex); rho=np.asarray(rho0,dtype=complex).copy()
    out=[rho.copy()]
    def f(R):
        dr=-1j*(H@R-R@H)
        for L in collapse_ops:
            LL=L.conj().T@L
            dr += L@R@L.conj().T - 0.5*(LL@R+R@LL)
        return dr
    for i in range(len(tlist)-1):
        h=tlist[i+1]-tlist[i]
        k1=f(rho); k2=f(rho+h*k1/2); k3=f(rho+h*k2/2); k4=f(rho+h*k3)
        rho += h*(k1+2*k2+2*k3+k4)/6
        rho=(rho+rho.conj().T)/2
        rho/=np.trace(rho).real
        out.append(rho.copy())
    return np.asarray(out)
