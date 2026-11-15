# Optional backend installation notes

The package detects native installations but does not attempt privileged system installs.

## Python backends

Install into the active environment as appropriate:

```bash
pip install jax jaxlib plotly pyvista meshio
pip install qutip kwant fipy
pip install python-flint
```

For FEniCSx/Firedrake/PETSc/mpi4py, use the installation method recommended by your HPC/site environment.

## Native tools

Make sure the following executables are on `PATH` when used:

```text
form
sage
M2
Singular
cadabra2
comsol
vasp_std
foamExec / blockMesh / simpleFoam
```

For deal.II and Geant4, compile the C++ examples against your local libraries.

For pySecDec, point the handoff directory at your installed pySecDec source/environment.

The `scripts/backend_report.py` command is always the first diagnostic to run.
