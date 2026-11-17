# Geometric Cusp Validation Stack v3.1

A symbolic and numerical test suite for the four-point cusp geometry of the
loop amplituhedron in planar N=4 super Yang–Mills: blow-ups and residues of
the one-loop canonical form, the two-loop logarithm integrand, and the contour
that extracts the cusp anomalous dimension.

## Status: what is established and what is not

**Reproduced from the literature (not new results).** The chart, the forms,
the finite functions \(F^{(0)}=-1\), \(F^{(1)}=\log^2z+\pi^2\), the weight-four
\(F^{(2)}\) and the contour rule are taken from the papers in
[docs/references.md](docs/references.md). The overall constant of the rule is
calibrated at one loop (\(\Gamma^{(1)}_{\rm cusp}=4\)). With that constant
fixed, the code returns \(\Gamma^{(2)}_{\rm cusp}=-4\pi^2/3\) and
\(\Gamma^{(3)}_{\rm cusp}=44\pi^4/45\); the three-loop value agrees to twelve
digits and was not used to set anything. These are checks of the rule and of
the transcribed functions. They are not derivations: \(F^{(1)}\) is recomputed
here from the finite parameter integral, but \(F^{(2)}\) is not recomputed
from the three-loop integrand.

**Computed here, exact in the stated charts.**

- The blow-up of the adjacent-propagator intersection \(y=w=0\) carries a
  simple radial pole with residue \(1/(x\,\eta\,z_0)\).
- The homogeneous blow-up of the strict collinear line has zero radial
  residue; the weighted blow-up with weights (1, 1, 2) has residue
  \(-1/(a\,z\,D)\) and a variable cross-ratio.
- The two-loop logarithm integrand scales as \(\rho^{-4}\) under the
  collective scaling, and its iterated residues depend on the order taken.
- Local contour analysis of a regulated lift of the finite integral
  ([docs/GATE_6_PROOF.md](docs/GATE_6_PROOF.md)).
- The contour of the Lorentzian chain at a light-like edge
  ([docs/GATE_6_CLOSURE.md](docs/GATE_6_CLOSURE.md)).

**Chain-level comparison (new, unreviewed).**
[docs/GATE_6_CLOSURE.md](docs/GATE_6_CLOSURE.md) and
[docs/GATE_6_BOUNDS.md](docs/GATE_6_BOUNDS.md) derive the unit-circle contour
from the Lorentzian integral itself.

- At every interior point of every light-like edge, the coefficient of the
  collinear singularity is \(-2\pi i\) times the unit-circle mean of \(F\); on
  the rest of the null line it is zero. Computed at all four cusps, on both
  edges of each, in all three real kinematic regions, to \(10^{-8}\).
- The real integration region, with the Feynman prescription, maps to a path
  in the \(\log z\) plane that deforms to the unit circle. It leaves the
  principal sheet, with \(|\arg z|\le2\pi\).
- The two-loop weight is regular at \(z=-1\) on the physical sheet. This is
  proved exactly from the polylogarithm inversion relations and Landen's
  identity. The edge chain with it gives \(-11\pi^4/15\), the value behind
  \(\Gamma^{(3)}_{\rm cusp}\).
- The collinear limit and the effect of the exact fourth propagator are
  proved with rates. The regulated integral of \(z^p\) is bounded uniformly
  in the regulator, which carries the result from powers to general weights.
- Only the four edge segments can pinch the integration contour.
- A hard cutoff on the collinear invariants does **not** reproduce the
  double pole for weights that grow at \(z\to0,\infty\). The double pole is
  taken in dimensional regularisation.

- An explicit integration contour is constructed
  ([docs/GATE_6_CONTOUR.md](docs/GATE_6_CONTOUR.md)): the real domain is
  displaced by an imaginary vector field built from four frame fields and
  two far fields. It stays on the Feynman side of every light cone and meets
  the cones only on the edges.

**Open.** The global properties of that contour follow from one proved
inequality plus a compactness argument, and are verified on a 60 000-point
sample; a proof tracking every constant is not written, and the construction
is done for spacelike diagonals only. None of this has been checked by a
specialist.
The weight \(F\) is still an input from
the Lorentzian integral, and the cycle is selected by the Lorentzian real
structure and the \(-i0\) prescription, not by positivity of the
amplituhedron. The relation \(dz/z=d\eta/\eta\) on the exceptional divisor
follows immediately from \(z=\eta z_0/x\) and is not evidence for any
identification. The regulated-cap programme of
[docs/GATE_6_PROOF.md](docs/GATE_6_PROOF.md) is bypassed, not resolved.

**Speculative and not implemented as a connected calculation.** The project
was motivated by a conjectured chain

```text
matrix model -> emergent geometry -> amplituhedron pushforward
  -> cusp boundary / residues -> Gamma_cusp -> IR entropy source
  -> continuum and quantum probes
```

Only the middle (amplituhedron, residues, contour) is computed. The other
stages are separate sandboxes with no derived link between them:

- `matrix_model.py` is a Langevin simulation of a small bosonic
  Yang–Mills matrix integral (default \(N=8\), \(D=4\)). It has no fermions
  and is not the IKKT model; nothing is passed from it to the geometry.
- `continuum.py` solves a 1D diffusion equation with a source
  \(\sigma=\Gamma_{\rm cusp}\log(w/\lambda)\) inserted by hand as an ansatz.
- `physics_maps.py` holds declared ansätze (\(m_{\rm eff}\), \(g_{\rm eff}\)),
  not results.
- The `canonical_proxy` weight in the pushforward sampler is a visual aid,
  not the canonical-form density.

`scripts/` contains the runners, `src/geocusp/` the implementation, and
`examples/` the 63 worked example groups with derivations.

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
- Bosonic matrix-model eigenvalue/commutator plots.
- Continuum entropy-source and field-evolution plots.

## 3. External adapters (templates, mostly untested)

`external/` holds short adapter or template files for FEniCSx, Firedrake,
deal.II, FiPy, Meep, COMSOL, Kwant, QuTiP, VASP, FORM, Geant4, OpenFOAM,
SageMath, Macaulay2, Singular, Arb/python-flint, pySecDec, Cadabra and xAct.
They are starting points for a user who has those packages. Most are a few
lines, none is exercised by the test suite, and none has been run against
the cusp calculation. The core results above need only NumPy, SciPy, SymPy
and mpmath. `scripts/backend_report.py` reports which packages are present
locally.

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

Run the benchmark with `make test && make examples`. The suite has 155 pytest
tests and 63 example groups (about five minutes for the tests). See the
[benchmark validation record](docs/BENCHMARK_VALIDATION.md) for the environment,
expected output, and established results. The [example guide](examples/README.md)
gives the individual commands and derivations. The remaining comparison to
the Lorentzian contour is tracked in the [Gate 6 proof](docs/GATE_6_PROOF.md).

Dedicated amplituhedron experiment:

```bash
PYTHONPATH=src python scripts/run_amplituhedron_deep.py --samples-k2 20000 --samples-k1 12000
```

## 6. Scientific interpretation

A numerical match to \(\Gamma_{\rm cusp}\) is meaningful only if the geometric
object producing it is defined independently of the known answer. At present
the match comes from the published finite functions and contour rule, so it
validates conventions, not the conjectured chain. The one-loop experiments
show where a simple residue construction succeeds and where it fails; the
two-loop experiments show that independent one-loop residues do not
factorize.
