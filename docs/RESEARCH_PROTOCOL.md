# Research Protocol: adversarial validation

The core hypothesis to test is not merely

`R_B = Gamma_cusp`.

The operational question is whether a geometrically defined period constructed from the user's signed/canonical geometry is mathematically identical to the established contour functional after all required residue and pushforward operations.

## Gate 0 — data provenance

Record:

- exact source formula;
- coupling convention;
- orientation convention;
- external data `Z`;
- gauge chart;
- random seed;
- precision;
- backend versions;
- git commit.

## Gate 1 — amplituhedron support

For each sample `C`:

1. verify every required ordered minor is positive;
2. compute `Y=CZ`;
3. compute normalized Plücker coordinates;
4. compare two independent positive charts;
5. visualize support and boundary concentration.

Failure condition: a purported positive pushforward leaves the positive domain beyond floating-point tolerance.

## Gate 2 — one-loop cusp geometry

Compute exact chart brackets and confirm

`<AB12>=y`, `<AB23>=w`, `<AB34>=z0`, `<AB14>=x`.

Then blow up

`y=rho*eta`, `w=rho`.

Verify symbolically that the exceptional form contains a simple `d rho/rho` pole and that

`z=(z0/x)*eta`.

Failure condition: the residue depends on arbitrary defining-function choices or the `z` map fails to preserve the expected `dlog` class.

## Gate 3 — strict cusp falsification

Blow up the strict line

`ZA=Zi`, `ZB=Z_{i-1}+alpha Z_{i+1}`.

Compute the radial residue exactly.

Current result in the supplied chart: zero. This means the strict line is not itself the required simple logarithmic divisor.

## Gate 4 — contour/linking test

Construct a nontrivial cycle around the deleted singular locus and compute its `z` image.

Use both ordinary winding and, when Arb is installed, interval-enclosed winding/certification.

Failure condition: the cycle is null-homologous in the punctured boundary or the contour functional cannot be recovered from the `z` image.

## Gate 5 — two-loop collective geometry

Use the exact two-loop negative-geometry form, then:

- derive the mutual bracket;
- perform collective blow-up;
- measure the radial power;
- compare symbolic and numerical scaling;
- identify whether an iterated residue is order-independent.

Failure condition: different residue orders yield incompatible forms without a mathematically specified normal-crossing/blow-up resolution.

## Gate 6 — independent cusp extraction

Do not fit a normalization.

Use the established finite function `F(g,z)` and compare:

`I[F]` versus your geometric period `R_B`.

The target relation is the source convention's RG extraction equation, with the sign/orientation dictionary recorded explicitly.

## Gate 7 — continuum layer

Only after the geometric source term is fixed, export it to FEniCSx/Firedrake/deal.II/FiPy/COMSOL.

Test existence, uniqueness, stability and convergence of the boundary PDE.

## Gate 8 — quantum layer

Only after the Hamiltonian/open-system map is explicitly specified, test it in Kwant/QuTiP.

## Gate 9 — phenomenological layer

VASP and Geant4 are downstream consistency checks, never substitutes for the geometric proof.
