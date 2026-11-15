from __future__ import annotations

import importlib.util
import os
import shutil
from dataclasses import dataclass


@dataclass(frozen=True)
class Backend:
    name: str
    python_module: str | None = None
    executable: str | None = None
    role: str = ""


BACKENDS = [
    Backend("fenicsx", "dolfinx", role="finite-element continuum PDEs"),
    Backend("firedrake", "firedrake", role="finite-element HPC/PETSc"),
    Backend("deal_ii", executable="dealii-parameter_handler", role="independent C++ FEM cross-check"),
    Backend("meep", "meep", role="FDTD electromagnetic boundary propagation"),
    Backend("qutip", "qutip", role="open quantum dynamics / entropy"),
    Backend("kwant", "kwant", role="quantum transport / boundary modes"),
    Backend("fipy", "fipy", role="finite-volume entropy/transport PDEs"),
    Backend("petsc4py", "petsc4py", role="HPC linear/nonlinear solvers"),
    Backend("flint", "flint", role="rigorous ball arithmetic"),
    Backend("pysecdec", "pySecDec", role="sector decomposition"),
]

EXECUTABLES = {
    "form": ["form"],
    "openfoam": ["foamExec", "blockMesh", "simpleFoam"],
    "comsol": ["comsol"],
    "vasp": ["vasp_std", "vasp_gam", "vasp_ncl"],
    "sage": ["sage"],
    "macaulay2": ["M2", "macaulay2"],
    "singular": ["Singular"],
    "cadabra": ["cadabra2"],
}


def detect() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for spec in BACKENDS:
        py = bool(spec.python_module and importlib.util.find_spec(spec.python_module))
        exe = bool(spec.executable and shutil.which(spec.executable))
        out[spec.name] = {
            "available": bool(py or exe),
            "python_module": spec.python_module,
            "python_importable": py,
            "executable": spec.executable,
            "executable_found": exe,
            "role": spec.role,
        }
    for name, candidates in EXECUTABLES.items():
        found = next((shutil.which(c) for c in candidates if shutil.which(c)), None)
        out[name] = {"available": found is not None, "command": found, "candidates": candidates}
    out["mpi"] = {"available": bool(os.environ.get("OMPI_COMM_WORLD_SIZE") or shutil.which("mpiexec"))}
    return out
