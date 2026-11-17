# Extended research suite

This document defines the experiments that the Python stack is intended to orchestrate. Each layer has a different falsification target.

## A. Positive Grassmannian / amplituhedron support

For the user's support equation

```math
\Omega^{(k,n)}\sim \int \frac{d^{k\times n}C}{\mathrm{Vol}[GL(k)]}
\,\delta^{(4k)}(Y-CZ)\,I(C),
```

the simulation samples a gauge-fixed positive chart, forms $`Y=CZ`$, computes Plücker coordinates and visualizes the pushforward support.

Run:

```bash
PYTHONPATH=src python scripts/run_amplituhedron_sweep.py --samples 10000
```

Primary failure conditions: negative ordered minors, pushforward inconsistencies, chart-dependent support pathologies, or boundary behavior that cannot be reproduced across independent charts.

## B. Exact one-loop geometry

The code evaluates the momentum-twistor brackets, canonical $`d\log`$ form, adjacent-propagator blow-up, strict-cusp blow-up, cross-ratio map and complex linking path.

```bash
PYTHONPATH=src python -m geocusp.one_loop_cli --contour-points 4096
```

Important outputs include the exceptional divisor residue, strict-cusp radial residue and winding of the induced $`z`$-path.

## C. Two-loop collective singularity

The code uses the explicit two-loop negative-geometry integrand, derives the collective radial scaling, scans all residue orderings, and compares exact symbolic scaling with numerical/JAX scaling.

```bash
PYTHONPATH=src python -m geocusp.two_loop_cli --mode all
PYTHONPATH=src python scripts/run_jax_scaling.py
```

The two-loop test is specifically designed to detect an invalid assumption of independent one-loop residue factorization.

## D. Established cusp benchmark

The benchmark separates the normalized finite function and the insertion-normalized function before applying the contour functional. This prevents an accidental sign or $`g^2`$ convention change from being mistaken for a physical result.

```bash
PYTHONPATH=src python scripts/run_benchmark.py
```

## E. Matrix substrate

An IKKT-like Langevin probe is included for commutator density, eigenvalue spectrum and autocorrelation diagnostics. It is a numerical exploration of the matrix substrate and is not a proof of the emergent QFT dictionary.

```bash
PYTHONPATH=src python -m geocusp.matrix_cli --N 12 --D 4 --g 0.7 --steps 10000 --thermalization 2000
```

## F. Continuum boundary / entropy probe

A reduced source-driven equation is exposed to FEniCSx/Firedrake/deal.II/FiPy/COMSOL. The local Python proxy provides a reference solution and convergence target; external backends provide independent numerical implementations.

## G. Quantum boundary probes

Kwant and QuTiP adapters are intentionally downstream. Their input Hamiltonians/dissipators must be replaced by the actual boundary dictionary before any physical claim is made.

## H. Microscopic and phenomenological probes

VASP, Geant4, Meep and OpenFOAM are downstream validators for specific derived constitutive/particle/transport maps. They are not used to retrofit the geometric conjecture.

## I. Certification

SageMath/Macaulay2/Singular should be used for algebraic-variety and singular-locus calculations. python-flint/Arb should be used where a numerical statement needs a certified enclosure. pySecDec is reserved for a genuinely derived regulated parameter integral.
