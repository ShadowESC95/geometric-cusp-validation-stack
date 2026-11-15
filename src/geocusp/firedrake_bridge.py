"""Firedrake/PETSc backend. Intended for HPC parameter sweeps."""
from __future__ import annotations

def solve(nx=128, length=8.0, m0=0.7, source_amp=1.0):
    try:
        from firedrake import IntervalMesh, FunctionSpace, TrialFunction, TestFunction, Function, Constant, SpatialCoordinate, inner, grad, dx, LinearVariationalProblem, LinearVariationalSolver
    except Exception as e:
        raise RuntimeError('Firedrake backend unavailable') from e
    mesh=IntervalMesh(nx,length,0)
    V=FunctionSpace(mesh,'CG',1)
    u=TrialFunction(V); v=TestFunction(V); uh=Function(V)
    a=inner(grad(u),grad(v))*dx + Constant(m0*m0)*u*v*dx
    L=Constant(source_amp)*v*dx
    p=LinearVariationalProblem(a,L,uh); LinearVariationalSolver(p).solve()
    return uh
