# Scientific software matrix

## Core exact/algebraic layer

**SymPy** — exact rational/algebraic manipulation, determinants, symbolic substitutions, residue limits.

**SageMath / Macaulay2 / Singular** — ideal membership, radicals, singular loci, dimension, blow-up/variety experiments.

**FORM** — high-performance symbolic HEP expressions when the SymPy expression graph becomes too large.

**Cadabra / xAct** — tensor/effective-field-theory algebra once the emergent metric is defined.

## Period/integration layer

**mpmath** — arbitrary-precision complex contour work.

**python-flint/Arb** — rigorous ball arithmetic and certified contour bounds.

**pySecDec** — regulated parameter integrals after a genuine integral representation has been derived.

## Continuum layer

**FEniCSx**, **Firedrake/PETSc**, **deal.II**, **FiPy**, **COMSOL** — PDE/variational/transport cross-checks.

## Field/transport layer

**Meep** — Maxwell/FDTD propagation from derived boundary constitutive data.

**Kwant** — lattice boundary-mode / transport calculation.

**QuTiP** — density-matrix dynamics, dissipation, entropy.

## Microscopic/downstream layer

**VASP** — local potential and charge-density field import for vacuum/VEV parameter studies.

**Geant4** — particle propagation/detector-level downstream consequences.

**OpenFOAM** — coarse-grained transport and irreversibility probes.

## What is deliberately NOT done

No package is used to retroactively force the conjecture to work. Each backend has a failure criterion and its outputs are stored separately under `results/` and `figures/`.
