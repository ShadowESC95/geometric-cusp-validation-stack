"""FiPy finite-volume entropy/source probe."""
from __future__ import annotations
import argparse
import numpy as np


def run(nx=400, length=12.0, dt=5e-4, steps=2500, D=0.25, gamma=1.0, log_ratio=3.0, width=1.0):
    try:
        from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm
    except Exception as exc:
        raise RuntimeError("FiPy unavailable") from exc
    mesh=Grid1D(nx=nx, dx=length/nx)
    x=np.asarray(mesh.cellCenters.value[0])
    sigma=gamma*log_ratio*np.exp(-(x/width)**2)
    S=CellVariable(mesh=mesh, name='S_OB_proxy', value=0.0)
    eq=TransientTerm()==DiffusionTerm(coeff=D)+sigma
    for _ in range(steps): eq.solve(var=S, dt=dt)
    return mesh, S, sigma

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--nx',type=int,default=400); a=p.parse_args(); _,S,_=run(nx=a.nx); print('FiPy max S:', float(np.max(S.value)))
