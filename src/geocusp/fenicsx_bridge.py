"""FEniCSx backend for a reduced boundary scalar equation.
Run where dolfinx/ufl/petsc4py are installed.
"""
from __future__ import annotations
import numpy as np

def solve(nx=96, length=8.0, m0=0.7, source_amp=1.0):
    try:
        from mpi4py import MPI
        from dolfinx import mesh, fem
        import ufl
    except Exception as e:
        raise RuntimeError('FEniCSx backend unavailable; install dolfinx/ufl/petsc4py/mpi4py') from e
    domain=mesh.create_interval(MPI.COMM_WORLD,nx,[-length/2,length/2])
    V=fem.functionspace(domain,("Lagrange",1))
    u=ufl.TrialFunction(V); v=ufl.TestFunction(V)
    x=ufl.SpatialCoordinate(domain)[0]
    f=source_amp*(1+0*x)
    a=ufl.inner(ufl.grad(u),ufl.grad(v))*ufl.dx + m0**2*u*v*ufl.dx
    L=f*v*ufl.dx
    problem=fem.petsc.LinearProblem(a,L,bcs=[]); uh=problem.solve();
    return uh
