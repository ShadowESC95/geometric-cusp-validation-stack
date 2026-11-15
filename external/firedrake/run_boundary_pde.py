"""Firedrake/PETSc boundary-PDE adapter with output to VTK."""
from __future__ import annotations
import argparse


def solve(nx=256, length=12.0, D=0.25, mass=0.7, source_amp=1.0):
    try:
        from firedrake import (
            IntervalMesh, FunctionSpace, TrialFunction, TestFunction, Function,
            Constant, SpatialCoordinate, inner, grad, dx, DirichletBC,
            LinearVariationalProblem, LinearVariationalSolver, exp, File,
        )
    except Exception as exc:
        raise RuntimeError("Firedrake unavailable") from exc
    mesh=IntervalMesh(nx, length)
    V=FunctionSpace(mesh, "CG", 2)
    u=TrialFunction(V); v=TestFunction(V); uh=Function(V)
    x=SpatialCoordinate(mesh)[0]
    sigma=Constant(source_amp)*exp(-((x-length/2)/1.0)**2)
    a=D*inner(grad(u),grad(v))*dx + mass**2*u*v*dx
    L=sigma*v*dx
    bc=DirichletBC(V,0.0,"on_boundary")
    problem=LinearVariationalProblem(a,L,uh,bcs=[bc]); LinearVariationalSolver(problem).solve()
    File("firedrake_boundary_solution.pvd").write(uh)
    return mesh,uh

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--nx',type=int,default=256); a=p.parse_args(); solve(nx=a.nx)
