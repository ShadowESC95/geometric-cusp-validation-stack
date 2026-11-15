# Cusp validation examples

Run from the repository root with the core Python dependencies installed.

```bash
make examples
```

The examples print their numerical and symbolic checks and raise an error on
a failed assertion. The full research campaign remains available with
`make full`. It covers the repository's sampling, matrix, and continuum runs
and uses the optional JAX dependency.

| Command | Checks |
| --- | --- |
| `PYTHONPATH=src python -m examples.geometry` | One-loop signed periods and winding, homogeneous and weighted residues, two-loop radial scaling and residue orders, complex-ray pole tracking and bypass jump, Lorentzian pole collision and partial mutual-pole integrals, finite-frame propagator signs and an interior Feynman-denominator zero |
| `PYTHONPATH=src python -m examples.cd_integral` | Finite CD integral on the positive axis, near-region recovery at fixed cross-ratio, complex contour extraction of the two-loop reference coefficient, uniform parameter-tail bound, and continuation with two finite cutoff paths |
| `PYTHONPATH=src python -m examples.chart_transport` | Exact transport of the weighted rational form and all five Feynman prescriptions, a family leaving the original affine chart, the external residue in the adjacent chart, and matching that residue across charts at a point and over a transverse patch |

## Finite-cutoff continuation

The continuation example keeps a cutoff on the projective parameter
`alpha4`. With `q=exp(-i*phi)`, the combined denominator is linear in
`q*alpha4`. The omitted piece depends on the endpoint through
`w=q*alpha4_max`. For the rotated cutoff, `alpha4_max=L/q` gives `w=L`.
For the positive-real cutoff, `alpha4_max=L` gives `w=q*L`. The example
compares both tails for three cutoff values and checks a uniform bound
of order `log(L)/L` for the tested positive cross-ratios. This retains the
cutoff dependence that drops out when the parameter is integrated to infinity.
Two direct finite-interval quadratures also check that the retained interval
plus its tail reproduces the full value.

References for the finite CD kernel and the cusp contour normalization are
[Arkani-Hamed, Henn and Trnka, Eq. (2.30)](https://arxiv.org/abs/2112.06956)
and [Dixon et al., Eqs. (6.2)–(6.6)](https://arxiv.org/abs/2605.28926).

## Scope of the checks

The CD parameter calculation reproduces the established two-loop coefficient.
The chart checks establish local matching with a transported prescription at
finite regulator. The parameter tail is an estimate in the reduced integral.
Constructing the global cycle on the weighted twistor divisor and controlling
its remaining chart-edge contributions are open steps in the geometric-period
proposal.
