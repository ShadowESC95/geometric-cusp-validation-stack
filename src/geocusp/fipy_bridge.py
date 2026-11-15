from __future__ import annotations
import numpy as np

def run(nx=200, length=8.0, dt=1e-3, steps=1000, D=0.25, source=1.0):
    try:
        from fipy import CellVariable, Grid1D, TransientTerm, DiffusionTerm, ImplicitSourceTerm
    except Exception as e:
        raise RuntimeError('FiPy backend unavailable') from e
    mesh=Grid1D(nx=nx, dx=length/nx); S=CellVariable(mesh=mesh,name='S',value=0.0)
    eq=TransientTerm() == DiffusionTerm(coeff=D) + source
    for _ in range(steps): eq.solve(var=S, dt=dt)
    return mesh,S
