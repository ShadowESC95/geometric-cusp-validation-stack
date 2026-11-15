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
| `PYTHONPATH=src python -m examples.energy_contour` | Real-energy quadrature versus the sum of five upper-half-plane residues, and the integrated weight of an isolated Feynman pinch in one spatial direction |
| `PYTHONPATH=src python -m examples.spatial_pinch` | Regulator scaling at a merger of spatial pinch roots and convergence of a finite three-dimensional spatial patch to the on-shell surface integral |

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

## Energy contour and spatial pinch

The energy example uses a finite Lorentzian frame at `rho=1/2` and the
combined numerator from Eq. (2.30). At positive regulator it compares direct
quadrature along the real energy axis with the sum of all five upper-half-plane
residues. With `p_x=p_y=3/10`, the `D2` upper pole and `D4` lower pole
meet at `p_z=-3/20`, `E=1/20` as the regulator vanishes. Integrating the
real part of the raw energy integral across this isolated pinch gives the
limiting weight `-250*pi**2/147`. The raw integral omits the overall
`1/(i*pi**2)` loop normalization. The transverse momenta are held fixed
throughout this check.

## Merging spatial roots

In the same Lorentzian frame, set `lambda=p_x**2+p_y**2-1/2` and
`q=p_z+3/4`. The opposite-energy poles of `D2` and `D4` coincide on
`lambda=-8*q**2/9` at `E=1/4-q/3`. At fixed `lambda=0`, the roots
merge, and the integral over `q` at `p_y=p_x>0`, multiplied by
`sqrt(eps)`, approaches `-8*pi**2/9`.

The surface is smooth in the coordinates `(q,lambda)`. Over a local patch
in those coordinates and `p_y/p_x`, the real part of the raw, regulated
energy integral has a finite spatial integral as the regulator tends to
zero. The example compares that convergence with the surface density
obtained from the two propagators on shell, including the Jacobian from
`(p_x,p_y)` to `(lambda,p_y/p_x)`.

## Scope of the checks

The CD parameter calculation reproduces the established two-loop coefficient.
The chart checks establish local matching with a transported prescription at
finite regulator. The parameter tail is an estimate in the reduced integral.
Constructing the global cycle on the weighted twistor divisor and controlling
its remaining chart-edge contributions are open steps in the geometric-period
proposal.
