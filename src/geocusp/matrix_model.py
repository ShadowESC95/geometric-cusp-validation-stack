from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class MatrixRun:
    steps: np.ndarray
    action: np.ndarray
    commutator_density: np.ndarray
    eigenvalues: np.ndarray
    X: np.ndarray


def hermitian_noise(N: int, D: int, rng: np.random.Generator) -> np.ndarray:
    out=[]
    for _ in range(D):
        A=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N))
        out.append((A+A.conj().T)/2.0)
    return np.asarray(out)


def action(X: np.ndarray, g: float) -> float:
    D,N,_=X.shape
    s=0.0
    for mu in range(D):
        for nu in range(mu+1,D):
            C=X[mu]@X[nu]-X[nu]@X[mu]
            s += -np.trace(C@C).real
    return float(N/(2*g*g)*s)


def gradient(X: np.ndarray, g: float) -> np.ndarray:
    D,N,_=X.shape
    grad=np.zeros_like(X,dtype=complex)
    # Gradient under Frobenius/trace pairing for the action above.
    for mu in range(D):
        G=np.zeros((N,N),dtype=complex)
        for nu in range(D):
            if mu==nu: continue
            C=X[mu]@X[nu]-X[nu]@X[mu]
            G += X[nu]@C-C@X[nu]
        G *= -N/(g*g)
        grad[mu]=(G+G.conj().T)/2.0
    return grad


def commutator_density(X: np.ndarray) -> float:
    D,N,_=X.shape; s=0.0
    for mu in range(D):
        for nu in range(mu+1,D):
            C=X[mu]@X[nu]-X[nu]@X[mu]
            s += -np.trace(C@C).real/N
    return float(s)


def remove_trace(X: np.ndarray) -> np.ndarray:
    """Project out the U(1) trace modes.

    The commutator action does not depend on them, so under Langevin noise
    they diffuse freely and smear the eigenvalue distribution.
    """
    N=X.shape[-1]
    tr=np.trace(X,axis1=1,axis2=2)/N
    return X-tr[:,None,None]*np.eye(N)


def langevin(N=8,D=4,g=0.7,step_size=3e-4,thermalization=1000,steps=5000,measure_every=10,seed=20260918):
    rng=np.random.default_rng(seed)
    X=0.2*hermitian_noise(N,D,rng)
    for _ in range(thermalization):
        X -= step_size*gradient(X,g) + np.sqrt(2*step_size)*hermitian_noise(N,D,rng)
        X=remove_trace(0.5*(X+X.conj().transpose(0,2,1)))
    rec=[]; acts=[]; comm=[]
    for s in range(steps):
        X -= step_size*gradient(X,g) + np.sqrt(2*step_size)*hermitian_noise(N,D,rng)
        X=remove_trace(0.5*(X+X.conj().transpose(0,2,1)))
        if s%measure_every==0:
            rec.append(s); acts.append(action(X,g)); comm.append(commutator_density(X))
    eig=np.linalg.eigvalsh(X[0])
    return MatrixRun(np.asarray(rec),np.asarray(acts),np.asarray(comm),eig,X)


def integrated_autocorrelation_time(series: np.ndarray, max_lag: int | None = None) -> float:
    x=np.asarray(series,dtype=float); x=x-x.mean()
    if x.size<3: return 0.5
    max_lag=max_lag or min(x.size//5,500)
    var=np.dot(x,x)/x.size
    if var<=0: return 0.5
    tau=0.5
    for lag in range(1,max_lag):
        ac=np.dot(x[:-lag],x[lag:])/((x.size-lag)*var)
        if ac<=0: break
        tau += ac
    return float(tau)
