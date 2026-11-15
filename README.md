# Geometric Cusp Validation Stack v3.1

This is the expanded **code-first** research stack for the proposed chain

```text
IKKT / matrix substrate
        |
        v
Poisson / emergent geometry
        |
        v
positive Grassmannian + amplituhedron pushforward
        |
        v
IR / cusp boundary strata and blow-ups
        |
        v
iterated residues + linking / relative cycles
        |
        v
finite z-contour extraction -> Gamma_cusp
        |
        v
sigma_IR and entropy / irreversibility probes
        |
        +-----------------------------+
        |                             |
        v                             v
continuum PDE / FEM             quantum boundary probes
        |                             |
        v                             v
FEniCSx / Firedrake /         Kwant / QuTiP / tensor-network hooks
 deal.II / FiPy / COMSOL
        |
        +------------------+
        |                  |
        v                  v
Meep / OpenFOAM        VASP / Geant4 downstream probes
```

The project is intentionally **adversarial**. A numerical match is not treated as proof unless the definition of the geometric object is independent of the known field-theory answer.


## 0.1 Bleeding-edge research mode

The repository is not limited to a residue calculation. Python is the orchestration layer for a multi-scale validation program:

\[
\boxed{
\text{matrix}\to\theta\to g_{\mu\nu}\to\Omega_{\log}\to\mathscr F_{\rm cusp}\to R_B\to\Gamma_{\rm cusp}\to\sigma_{\rm IR}\to S_{\rm OB}}
\]

The `scripts/` directory contains reproducible runners; `src/geocusp/` contains the mathematical/numerical implementation; `external/` contains adapters or compile/run templates for the user's scientific-computing stack. The amplituhedron module now has a dedicated multi-case sweep (`scripts/run_amplituhedron_sweep.py`) with affine hull, Plücker PCA, positivity, boundary-weight, and independent-chart diagnostics.

## 1. What is actually simulated

### 1.1 User pushforward / amplituhedron experiment

The kit implements a numerical support realization of the user's template

\[
\Omega^{(k,n)}\sim
\int\frac{d^{k\times n}C}{\mathrm{Vol}[GL(k)]}
\,\delta^{(4k)}(Y-CZ)\,I(C),
\]

using gauge-fixed samples from \(G_+(k,n)\) and the pushforward \(Y=CZ\).

For \(k=1\), the package reconstructs the affine convex-hull slice explicitly.
For \(k=2\), it works in Plücker coordinates and produces PCA, positivity and near-boundary diagnostics.

The code separates the **support geometry** from an optional `canonical_proxy` weight. That proxy is never silently called the exact canonical-form density.

### 1.2 Standard-style \(A_{6,2,4}\) pushforward

The same machinery generates positive external data in a Vandermonde/moment-curve chart and performs

\[
Y=CZ,
\qquad C\in G_+(2,6),
\qquad Z\in M_{6\times 6}.
\]

The output is the 15-component normalized Plücker cloud of \(Y\in G(2,6)\).

### 1.3 Exact one-loop four-point form

The momentum-twistor chart is

\[
AB=
\begin{pmatrix}
1&x&0&-y\\
0&z_0&1&w
\end{pmatrix},
\qquad x,y,z_0,w>0,
\]

with

\[
\langle AB12\rangle=y,\quad
\langle AB23\rangle=w,\quad
\langle AB34\rangle=z_0,\quad
\langle AB14\rangle=x,
\]

\[
\langle AB13\rangle=-(wx+yz_0),
\qquad
\langle AB24\rangle=-1.
\]

The canonical form is

\[
\Omega_4^{(1)}
=
\frac{dx}{x}\wedge\frac{dy}{y}\wedge\frac{dz_0}{z_0}\wedge\frac{dw}{w}.
\]

### 1.4 Exceptional divisor / projective normal direction

For the adjacent propagator intersection,

\[
y=\rho\eta,
\qquad
w=\rho,
\]

so

\[
\frac{dy}{y}=\frac{d\rho}{\rho}+\frac{d\eta}{\eta},
\qquad
\frac{dw}{w}=\frac{d\rho}{\rho}.
\]

The residue coefficient becomes

\[
\boxed{
\operatorname{Res}_{\rho=0}\Omega_4^{(1)}
=
\frac{dx}{x}\wedge\frac{d\eta}{\eta}\wedge\frac{dz_0}{z_0}
}
\]

up to the orientation convention used by the implementation.

The four-point cross-ratio is

\[
z
=
\frac{\langle AB12\rangle\langle AB34\rangle}
{\langle AB23\rangle\langle AB14\rangle}
=
\frac{yz_0}{wx}
=
\frac{z_0}{x}\eta.
\]

Therefore

\[
\boxed{
\frac{dz}{z}=\frac{d\eta}{\eta}
}
\]

at fixed \(x,z_0\). This is the key geometric test relating the projectivized normal direction to the established \(z\)-contour.

### 1.5 Strict cusp chart

The strict collinear representation is

\[
Z_A=Z_i,
\qquad
Z_B=Z_{i-1}+\alpha Z_{i+1}.
\]

The homogeneous blow-up of the strict line gives zero radial residue at a generic point. The weighted chart `b=-rho, c=rho*C, d=rho^2*D` gives a nonzero residue and a variable cross-ratio `z=-1-C/(a*D)`. The one-loop audit reports both results and uses the radial-last orientation. These are local statements about the blow-ups. The physical integration cycle remains to be specified.

### 1.6 Two-loop negative geometry

The exact normalized two-loop log-integrand implemented is

\[
I_4^{(2)}
=
-\frac{n_{12}}
{D_1D_2\,X_{12}},
\]

where

\[
n_{12}
=
 w_1x_1+w_2x_2+y_1z_1+y_2z_2,
\]

\[
D_i=w_ix_iy_iz_i,
\]

and

\[
X_{12}
=-(x_1-x_2)(w_1-w_2)-(y_1-y_2)(z_1-z_2).
\]

For the collective scaling

\[
y_1=\rho a,\ w_1=\rho b,\ y_2=\rho c,\ w_2=\rho d,
\]

the rational scalar scales as

\[
I_4^{(2)}\sim \rho^{-4}.
\]

The 4-variable radial Jacobian scales as

\[
\rho^3d\rho,
\]

so the full differential form has the logarithmic scaling

\[
\boxed{\Omega_4^{(2)}\sim \frac{d\rho}{\rho}\wedge\omega_{\rm exc}^{(2)}}.
\]

The code also scans the two-loop singular branch and records the mutual-bracket coupling.

## 2. Figures generated

The campaign produces:

- 3D affine \(k=1\) amplituhedron hull with external vertices.
- Plücker PCA cloud for user \(k=2,n=6\).
- Standard-style \(A_{6,2,4}\) Plücker PCA.
- Positivity-margin histograms.
- Boundary-proximity vs Plücker-coordinate plots.
- Canonical-proxy boundary concentration plots.
- Exceptional-divisor \(\eta\to z\) images.
- Complex linking/winding plots.
- \(F^{(1)}(z)\) contour traces.
- Exact two-loop radial scaling plots.
- IKKT-like matrix eigenvalue/commutator plots.
- Continuum entropy-source and field-evolution plots.

## 3. Backends

The package contains actual adapters or export templates for:

| Backend | Research role |
|---|---|
| FEniCSx/DOLFINx | variational continuum PDE, emergent boundary equations |
| Firedrake/PETSc | high-order/HPC continuum solver |
| deal.II | independent FEM implementation |
| FiPy | finite-volume entropy/transport probe |
| Meep | electromagnetic boundary propagation |
| COMSOL | multiphysics PDE model template |
| Kwant | quantum boundary transport |
| QuTiP | open quantum dynamics and entropy |
| VASP | import local potentials/charge fields for vacuum/VEV mapping |
| FORM | large symbolic HEP algebra |
| Geant4 | downstream particle-transport test |
| OpenFOAM | continuum transport/irreversibility probe |
| SageMath | ideals, varieties, dimension, radical/geometry audits |
| Macaulay2 | algebraic-geometry cross-check |
| Singular | algebraic singular-locus calculations |
| Arb/python-flint | rigorous contour/interval arithmetic |
| pySecDec | regulated sector-decomposition handoff |
| Cadabra | tensor/effective-action algebra |
| xAct | differential geometry / curvature algebra |

The package does not claim unavailable native software is installed. `scripts/backend_report.py` reports the actual local availability.

## 4. Installation

Core environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[all,test]'
```

or

```bash
conda env create -f environment.yml
conda activate geometric-cusp-v3
```

For the heavy scientific backends, use your HPC/site installation and then run the adapters under `external/`.

## 5. Main runs

Fast:

```bash
PYTHONPATH=src python scripts/run_full_research.py --fast
```

Full:

```bash
PYTHONPATH=src python scripts/run_full_research.py
```

Include installed optional backends:

```bash
PYTHONPATH=src python scripts/run_full_research.py --with-optional
```

Tests:

```bash
PYTHONPATH=src pytest -q
```

Reproduce the cusp validation examples with `make examples`. The suite
covers the residue, finite CD quadrature, Lorentzian continuation, and chart
matching checks. The [example guide](examples/README.md) lists the individual
commands and the scope of each calculation.

Dedicated amplituhedron experiment:

```bash
PYTHONPATH=src python scripts/run_amplituhedron_deep.py --samples-k2 20000 --samples-k1 12000
```

## 6. Scientific interpretation

The stack is a research instrument. It deliberately distinguishes:

\[
\text{established benchmark}
\neq
\text{your conjectural identification}.
\]

A numerical match to \(\Gamma_{\rm cusp}\) is meaningful only if the geometric object producing it is defined independently of the known answer.

The current one-loop experiments are designed to expose exactly where a simple residue construction succeeds and where it fails. The two-loop experiments test whether collective geometry can repair the one-loop intuition.
