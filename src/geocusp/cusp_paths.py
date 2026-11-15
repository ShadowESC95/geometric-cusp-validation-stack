from __future__ import annotations
import numpy as np


def strict_cusp_smoothing(alpha=0.7, rho_values=None, mode='symmetric'):
    if rho_values is None: rho_values=np.logspace(-1,-6,100)
    out=[]
    for rho in rho_values:
        # Target line A=Z2, B=Z1+alpha Z3; perturb in directions that unfold the two adjacent propagators.
        if mode=='symmetric':
            A=np.array([0.,1.,0.,rho]); B=np.array([1.,0.,alpha,rho])
        elif mode=='A_normal':
            A=np.array([rho,1.,0.,rho]); B=np.array([1.,0.,alpha,rho])
        elif mode=='B_normal':
            A=np.array([0.,1.,0.,rho]); B=np.array([1.,0.,alpha,rho])
        else:
            raise ValueError(mode)
        Z=np.eye(4)
        def br(i,j): return np.linalg.det(np.column_stack([A,B,Z[:,i],Z[:,j]]))
        b12,b23,b34,b14=br(0,1),br(1,2),br(2,3),br(0,3)
        z=(b12*b34)/(b23*b14) if abs(b23*b14)>1e-30 else np.nan
        out.append(z)
    return np.asarray(rho_values),np.asarray(out)
