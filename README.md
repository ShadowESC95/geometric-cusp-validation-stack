# Geometric Cusp Validation Stack

A symbolic and numerical study of the four-point cusp geometry of the loop
amplituhedron in planar N=4 super Yang–Mills: blow-ups and residues of the
canonical form, the two-loop logarithm integrand, and the contour that
extracts the cusp anomalous dimension.

It has two main results.

- **The Lorentzian origin of the contour.** At a light-like edge the real
  integration region, with the Feynman prescription, maps to a path between
  the two lips of $`z=-1`$; the unit circle is one representative.
- **The four-cusp decomposition.** The double pole of the integral is a sum
  of four terms, one per cusp (proved for its Feynman-parameter form). Each
  term is a triangle with two light-like sides, and is derived from the
  real Lorentzian domain as an iterated integral: it is supported on the
  edge segment, and it factorises into the position along the edge, the
  collinear scale and the path between the lips.

Together they explain the known contour rule and give two proofs of it for
the three benchmark weights. One is in position space, for spacelike
diagonals, for the Feynman integral defined as the limit of the real
integral with $`-i\varepsilon`$. The other covers all three kinematic
regions and is conditional on a computer-verified certificate. The conditions are listed under
[Scope and limits](#scope-and-limits).

- **355 tests** and **63 worked example groups**, all passing
- Core dependencies: NumPy, SciPy, SymPy, mpmath
- MIT licensed

## Contents

- [Status at a glance](#status-at-a-glance)
- [Quick start](#quick-start)
- [What is computed](#what-is-computed)
- [The contour result](#the-contour-result)
- [The four-cusp decomposition](#the-four-cusp-decomposition)
- [Solver](#solver)
- [Documentation](#documentation)
- [Repository layout](#repository-layout)
- [Scope and limits](#scope-and-limits)
- [Other components](#other-components)
- [Citation and licence](#citation-and-licence)

## Status at a glance

| Result | Status |
| --- | --- |
| One-loop chart, brackets and canonical form | Exact (computer algebra) |
| Blow-up residues at the cusp, ordinary and weighted | Exact |
| Two-loop integrand: scaling and order dependence of residues | Exact |
| Two-loop integrand derived from the inequalities of its region | Exact (re-derivation of a published form) |
| Sampler tested by the sign-flip characterisation, $`k=1,2`$ | Passes; agrees with linear programming at $`k=1`$ |
| Cusp coefficients $`\Gamma^{(1)}, \Gamma^{(2)}, \Gamma^{(3)}`$ from the contour rule | Reproduced from published weights |
| Lorentzian edge chain equals the unit-circle mean | Computed to $`10^{-8}`$ in every case tested |
| Two-loop weight is regular at $`z=-1`$ | Proved |
| Explicit integration contour, all kinematic regions | Constructed; conditions verified in interval arithmetic over a range of squares |
| Regulated integral of $`z^p`$ equals the Mellin–Barnes form | Written proof, for $`\vert\operatorname{Re}p\vert`$ below the regulator |
| Double pole equals $`4M[F]`$ for the weights $`\bar F^{(0)},\bar F^{(1)},\bar F^{(2)}`$ | Written proof, at every square where the certificate passes; the certificate is verified in interval arithmetic for $`0.023\le s/t\le33.6`$ and the corresponding ranges with a timelike diagonal |
| Cusp-local model: a light-like triangle, double pole $`\sin\pi p/(\pi p)`$ | Exact closed form |
| Double pole is the sum of the four cusp-local models, without Mellin–Barnes | Written proof (parametric integral). All three regions; independently checked |
| Exact formula for the cusp-local model of a general weight: $`F`$ on the path between the lips of $`z=-1`$, smeared over width $`\eta`$ | Written proof for the model defined through its parametric form; independently checked |
| Cusp-local model derived from the Lorentzian domain: support on the edge segment, discontinuity across $`z<0`$, exact factorisation | Written proof for the iterated integral in a fixed order; independently checked |
| Full integrand: support in a null slab and exact split into two edge pieces | Written proof; independently checked |
| Leading term at both ends of an edge piece: one unit of $`M[z^p]`$ each | Closed form; independently checked |
| Each edge piece equals its two leading terms up to a single pole | Written proof for an explicit cube representation, uniform in $`\eta`$; independently checked |
| Double pole of the full integral from the Lorentzian domain, weight $`z^p`$, spacelike diagonals | Written proof; independently checked. Iterated integral in a fixed order, $`-i0`$ removed first |
| The same for the three benchmark weights: double pole $`4M[F]`$, spacelike diagonals | Written proof; independently checked. Same iterated integral; no contour, certificate or Mellin–Barnes |
| The same with the $`-i0`$ removed last: the limit of the absolutely convergent real integral | Written proof; independently checked (one range corrected on review) |
| The iterated integral equals the integral over the global contour | Written proof, conditional on the $`D`$-dimensional properties of the contour, which are asserted and sampled, not proved |
| A timelike diagonal in position space | Exact slab and keyhole structure; one of the two edge pieces carries its two units; the other is **not derived** (covered via the parametric form) |
| Review by a specialist | **Not done** |

The contour argument was written and machine-checked in one working
session. An internal review found one substantial gap and several smaller
ones, which are repaired in the documents. A second internal review of the
repairs found them correct in substance and added qualifications, which are
stated. It has not been read by a specialist. Treat it as a proposed
proof.

## Quick start

Python 3.10 or later.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'

make test        # 355 tests, about twenty minutes
make examples    # 63 example groups, about three minutes
```

Individual pieces:

```bash
PYTHONPATH=src python -m geocusp.solver                       # desktop window
PYTHONPATH=src python -m geocusp.solver functional            # the four functionals on a weight
PYTHONPATH=src python -m geocusp.solver edge --position 0.4   # Lorentzian edge chain
PYTHONPATH=src python -m geocusp.solver two-loop              # three-loop coefficient
PYTHONPATH=src python -m geocusp.solver contour --region=+-   # survey of the contour, with the range of arg z
PYTHONPATH=src python -m geocusp.solver certificate           # the inequalities the contour proof needs
PYTHONPATH=src python -m geocusp.solver threshold --sigma 0.3 # convergence of |z^p|: needs |Re p| < eta
PYTHONPATH=src python -m examples.lorentzian_all_cusps        # one example group
```

## What is computed

### One-loop geometry

In the momentum-twistor chart

```math
AB=\begin{pmatrix}1&x&0&-y\\0&z_0&1&w\end{pmatrix},\qquad x,y,z_0,w>0,
```

the canonical form is

```math
\Omega_4^{(1)}=\frac{dx}{x}\wedge\frac{dy}{y}\wedge\frac{dz_0}{z_0}\wedge\frac{dw}{w},
```

and the four-point cross-ratio is

```math
z=\frac{\langle AB12\rangle\langle AB34\rangle}{\langle AB23\rangle\langle AB14\rangle}=\frac{y\,z_0}{w\,x}.
```

- Blowing up the adjacent-propagator intersection $`y=w=0`$ with
  $`y=\rho\eta`$, $`w=\rho`$ gives a simple radial pole with residue
  $`1/(x\,\eta\,z_0)`$.
- The ordinary blow-up of the strict collinear line has zero radial residue.
  The weighted blow-up with weights (1, 1, 2) has residue $`-1/(a\,z\,D)`$
  and a variable cross-ratio.

### Two-loop logarithm integrand

```math
I_4^{(2)}=-\frac{n_{12}}{D_1D_2\,X_{12}},\qquad
n_{12}=w_1x_1+w_2x_2+y_1z_1+y_2z_2,\qquad D_i=w_ix_iy_iz_i .
```

- This integrand is derived here from the inequalities of the region
  $`\langle AB_1AB_2\rangle<0`$, by cutting it into five generalised
  simplices and adding their forms
  ([FROM_GEOMETRY](docs/FROM_GEOMETRY.md)).
- Under the collective scaling $`y_i,w_i\sim\rho`$ the integrand scales as
  $`\rho^{-4}`$, so the form has a logarithmic radial pole.
- The 24 orders of iterated residues give four distinct values. Independent
  one-loop residues do not factorise at two loops.

### Cusp anomalous dimension

With the finite weights $`\bar F^{(L-1)}(z)`$ and the unit-circle mean

```math
M[F]=\frac{1}{2\pi}\int_{-\pi}^{\pi}F(e^{i\phi})\,d\phi,\qquad
\Gamma^{(L)}_{\rm cusp}=-\frac{4}{L}\,M\big[\bar F^{(L-1)}\big],
```

the code returns

| Loop order | Weight | Mean | Coefficient |
| --- | --- | --- | --- |
| 1 | $`-1`$ | $`-1`$ | $`4`$ |
| 2 | $`\log^2z+\pi^2`$ | $`2\pi^2/3`$ | $`-4\pi^2/3`$ |
| 3 | published two-loop weight | $`-11\pi^4/15`$ | $`44\pi^4/45`$ |

The overall constant of the rule is fixed at one loop. Orders two and three
check the published rule and the transcribed weights against known values. The weights themselves are taken from the
literature; the one-loop-lower weight is also recomputed here from its
parameter integral.

## The contour result

Take the frozen-loop integral over the insertion point $`x_0`$, with every
invariant carrying the Feynman $`-i0`$. Near the light-like edge through
$`x_j`$ and $`x_{j+1}`$, the integral over the directions transverse to the
edge produces a collinear singularity whose coefficient is

```math
\frac{d\beta}{\beta(1-\beta)}\times
\begin{cases}-2\pi i\,M[F], & 0<\beta<1 \quad\text{(on the edge)}\\[2pt] 0, & \text{elsewhere on the null line.}\end{cases}
```

The factor $`M[F]`$ arises because the image of the real integration region
in the $`\log z`$ plane is a path from one lip of $`z=-1`$ to the other,
which can be deformed to the unit circle. The path
leaves the principal sheet, with $`|\arg z|\le2\pi`$, and the weighted
blow-up above is exactly the collinear scaling of the edge.

What supports it:

- exact identities for all four invariants in global coordinates;
- the statement computed at all four cusps, on both edges of each, in all
  three real kinematic regions;
- the two-loop weight continued onto its second sheet, giving
  $`-11\pi^4/15`$;
- an explicit, continuously differentiable integration contour whose
  properties reduce to inequalities evaluated by a certificate, for 90
  squares across the three regions;
- written proofs of absolute convergence (for polynomially bounded weights,
  and for $`z^p`$ with $`|\operatorname{Re}p|`$ below the regulator) and of
  the equality with the Mellin–Barnes form of the regulated integral.

What it does not do:

- The double pole of the full integral is not derived from the edge
  statement in position space. It is proved from the parametric integral:
  each cusp contributes one unit of $`M[z^p]`$, and the two poles belong to
  the two edges at the cusp
  ([CUSP_LOCAL_MODEL](docs/CUSP_LOCAL_MODEL.md)), in all three regions. For
  a general weight the cusp-local model, defined through its parametric
  form, is given exactly by an integral along the path between the lips,
  with an $`O(\eta)`$ correction, and is derived from the real Lorentzian
  domain. What is missing is the position-space proof that the full
  integral reduces to the four cusp terms.
- The three-loop value from the edge chain equals the circle mean by
  Cauchy's theorem, so it is not an independent check.
- With a timelike diagonal, the sharp bound $`|\arg z|\le2\pi`$ on the
  contour is measured, not proved. The weaker bound the argument needs,
  below $`3\pi`$, is proved.

Earlier statements found to be wrong and corrected in the documents: the
exchange of integrals was claimed on lines where the integral does not
converge absolutely; a hard cutoff on the collinear invariants does
not reproduce the double pole for weights that grow at $`z\to0,\infty`$, and
dimensional regularisation cannot be treated as a real shift on a
four-dimensional contour.

## The four-cusp decomposition

Freeze the far propagator at a cusp. What is left is a triangle with two
light-like sides. In light-cone coordinates
$`x_0=x_2+\alpha p_1+\beta p_2+y_\perp`$ all three invariants are linear in
$`\alpha`$, and closing the $`\alpha`$ contour gives, exactly:

- the integral vanishes unless $`0<\beta<1`$: it is supported on the edge
  segment;
- on the segment it is a discontinuity across the negative $`z`$ axis: the
  two sides of the edge are the two lips;
- the rest factorises into three independent integrals,

```math
i\int_0^1d\beta\,\beta^{\eta-1}(1-\beta)^\eta\int_0^\infty\frac{dy\,y^{\eta-1}}{1+y}
\int_{-\pi}^{\pi}d\phi\,(1+e^{-i\phi})^\eta\,F\big(e^{i\phi}/((1-\beta)(1+y))\big),
```

  over the position along the edge, the collinear scale and the path
  between the lips.

Each of the first two gives a pole $`1/\eta`$, and the third gives
$`2\pi M[F]`$. Because the measure is a product, the region where both
collinear scales vanish needs no separate estimate. For $`F=z^p`$,
$`\eta-1<\operatorname{Re}p<1+\eta`$, the result is

```math
2\pi i\,\Gamma(1+\eta)\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1-p+2\eta)}\;\approx\;\frac{2\pi i}{\eta^2}\,\frac{\sin\pi p}{\pi p},
```

and the four cusps add to the known coefficient. The derivation is for the
iterated integral in a fixed order (the light-cone coordinate $`\alpha`$
first); the integral is not absolutely convergent, and other orders are
not treated.

For the full integrand, with the far propagator, the same contour closing
is still exact. All four invariants are linear in $`\alpha`$, and the
integral along that null direction is supported in the slab between the
null hyperplanes through $`x_4`$ and $`x_3`$. Inside the slab it is a
discontinuity across the light cone of $`x_3`$ on one side of the edge
$`(1,2)`$ and across that of $`x_4`$ on the other: one piece for the edge
$`(2,3)`$ and one for the edge $`(4,1)`$.

That the full integral equals the sum of the four cusp terms up to a single
pole is proved for the Feynman-parameter form, by elementary estimates on
the simplex, in all three regions. In position space the leading term at each
end of each edge piece is computed, and each gives one unit: at the near
end one pole from the endpoint and one from a collinear limit, at the far
end half a pole from the endpoint and two collinear limits. The estimate
that the rest is a single pole is proved, in variables in which each edge
piece is a positive integral over a cube. This gives the double pole of
the full integral from the Lorentzian domain for the weight $`z^p`$ and
spacelike diagonals, as an iterated integral in a fixed order with the
$`-i0`$ removed first, and by superposition for the three benchmark
weights: the double pole is $`4M[F]`$. This route uses no global contour,
no certificate and no Mellin–Barnes. A timelike diagonal is covered only
through the parametric form. See
[CUSP_LOCAL_MODEL](docs/CUSP_LOCAL_MODEL.md).

## Solver

`geocusp.solver` evaluates four functionals of a weight
$`F=\sum_kc_k\log^kz+Az^p`$: the circle mean, the period over the capped
cycle on the weighted divisor, the Mellin–Barnes double pole, and the
Lorentzian edge chain. It also handles the two-loop weight, surveys the
integration contour in any of the three kinematic regions, evaluates the
certificate of its conditions, and shows the convergence threshold for
powers. It runs locally, with no server and no network access.
The desktop window needs Tk; every function is also available from the
command line as shown in [Quick start](#quick-start).

## Documentation

| Document | Contents |
| --- | --- |
| [GATE_6_CLOSURE](docs/GATE_6_CLOSURE.md) | The edge statement and the comparison of the two integration chains |
| [GATE_6_BOUNDS](docs/GATE_6_BOUNDS.md) | All cusps and regions, error bounds, the regulator, and a status table |
| [GATE_6_CONTOUR](docs/GATE_6_CONTOUR.md) | The explicit contour, its certificate, convergence, and the regulated integral |
| [CUSP_LOCAL_MODEL](docs/CUSP_LOCAL_MODEL.md) | The cusp-local triangle, and the double pole as a sum over cusps without Mellin–Barnes |
| [GATE_6_FUNCTIONAL](docs/GATE_6_FUNCTIONAL.md) | The functionals on the Mellin family, and the solver |
| [FROM_GEOMETRY](docs/FROM_GEOMETRY.md) | The two-loop form from its inequalities; the membership test for the sampler |
| [GATE_6_PROOF](docs/GATE_6_PROOF.md) | The earlier regulated-cap approach and its obstruction |
| [MATHEMATICAL_MODEL](docs/MATHEMATICAL_MODEL.md) | Charts, forms and integrands used by the code |
| [BENCHMARK_VALIDATION](docs/BENCHMARK_VALIDATION.md) | How to reproduce the suite, and what each check establishes |
| [Example guide](examples/README.md) | The 63 example groups with their derivations |
| [References](docs/references.md) | Papers and software |

Start with GATE_6_CLOSURE, then section 7 of GATE_6_BOUNDS for the status of
every statement.

## Repository layout

```text
src/geocusp/      implementation
  residue.py, negative_geometry.py      one- and two-loop algebra
  benchmark.py, contours.py             weights and the contour rule
  lorentzian_cusp.py, lorentzian_edges.py   the edge chain, all cusps and regions
  feynman_contour.py, schwinger.py      the contour, its certificate, the regulated integral
  solver.py, solver_gui.py              command line and desktop window
examples/         63 worked example groups (python -m examples.run_suite)
tests/            355 tests
docs/             derivations and status
scripts/          runners for the wider campaign
external/         adapter templates for other packages (untested)
```

## Scope and limits

- **The top weight is an input.** $`\bar F^{(2)}`$ comes from the
  literature. At one order lower the chain is computed here: the two-loop
  form from its inequalities, then $`\bar F^{(1)}`$ by quadrature. The same
  has not been done at three loops.
- **The contour is selected by physics, not positivity.** It follows from
  the Lorentzian real structure and the Feynman prescription. It is not
  derived from positivity of the amplituhedron.
- **The circle is already in the literature.** arXiv:2112.06956 eq. (2.19)
  writes the functional as a contour between the two lips of $`z=-1`$. What
  is added here is its Lorentzian origin at a light-like edge, an explicit
  contour, and hypotheses with a proof for a class of weights. See
  [GATE_6_CLOSURE](docs/GATE_6_CLOSURE.md) section 9.
- **Reproduction is not derivation.** The three cusp coefficients test the
  rule and the transcribed weights against known values.
- **The relation $`dz/z=d\eta/\eta`$ is immediate** from $`z=\eta z_0/x`$ and
  is not evidence for any identification.
- **The earlier regulated-cap programme** of GATE_6_PROOF is bypassed, not
  resolved.
- **No specialist review.**

## Other components

The repository began as a wider programme linking a matrix model, emergent
geometry, the amplituhedron and an entropy source. Only the amplituhedron and
contour parts are developed. The rest are separate sandboxes with no derived
link between them:

- `matrix_model.py`: Langevin simulation of a small bosonic Yang–Mills matrix
  integral. It has no fermions and is not the IKKT model.
- `continuum.py`: a one-dimensional diffusion equation with a source inserted
  by hand.
- `amplituhedron.py`: samples the positive Grassmannian and pushes forward.
  The samples are now tested by the sign-flip characterisation, which does
  not use the construction. The `canonical_proxy` weight is a visual aid,
  not the canonical form.
- `physics_maps.py`: declared ansätze, not results.
- `external/`: short adapter templates for FEniCSx, Firedrake, deal.II, Meep,
  QuTiP and others. None is exercised by the tests.

`make full` runs that wider campaign.

## Citation and licence

Released under the [MIT licence](LICENSE). Citation metadata is in
[CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button uses it.

The physics this builds on is in the papers listed in
[docs/references.md](docs/references.md), in particular Arkani-Hamed, Henn
and Trnka, *Nonperturbative Negative Geometries* (arXiv:2112.06956), and
Dixon et al., *Multi-Loop Negative Geometries* (arXiv:2605.28926).
