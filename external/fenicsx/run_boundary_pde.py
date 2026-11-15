"""FEniCSx/DOLFINx executable adapter.

Solves the reduced 1-D boundary equation

    -D d2S/dx2 + m_eff(x)^2 S = sigma_IR(x)

with homogeneous Dirichlet boundary conditions. This is deliberately a
continuum *probe* of the proposed boundary source, not a derivation of S_OB.
"""
from __future__ import annotations

import argparse
import numpy as np


def solve(nx=256, length=12.0, D=0.25, mass=0.7, source_amp=1.0, width=1.0):
    try:
        from mpi4py import MPI
        from dolfinx import mesh, fem
        import ufl
    except Exception as exc:
        raise RuntimeError("Install FEniCSx/DOLFINx + mpi4py + PETSc to run this adapter") from exc

    domain = mesh.create_interval(MPI.COMM_WORLD, nx, [-length / 2, length / 2])
    V = fem.functionspace(domain, ("Lagrange", 2))
    u = ufl.TrialFunction(V); v = ufl.TestFunction(V)
    x = ufl.SpatialCoordinate(domain)[0]
    sigma = source_amp * ufl.exp(-(x / width) ** 2)
    a = D * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx + mass**2 * u * v * ufl.dx
    L = sigma * v * ufl.dx
    # Homogeneous Dirichlet BC.
    bc = fem.dirichletbc(fem.Constant(domain, 0.0), fem.locate_dofs_geometrical(V, lambda X: np.isclose(X[0], -length/2) | np.isclose(X[0], length/2)), V)
    from dolfinx.fem.petsc import LinearProblem
    problem = LinearProblem(a, L, bcs=[bc])
    uh = problem.solve()
    return domain, V, uh


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--nx", type=int, default=256); a=p.parse_args()
    domain, V, uh = solve(nx=a.nx)
    if domain.comm.rank == 0:
        print("FEniCSx boundary PDE completed; solution vector size:", uh.x.array.size)
