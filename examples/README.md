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
| `PYTHONPATH=src python -m examples.pinch_intersection` | A local intersection of two timelike pinch surfaces, including the principal values and the contribution supported where both surfaces meet |
| `PYTHONPATH=src python -m examples.closed_pinch_curve` | Integration along the full closed intersection and regulator convergence of the lower `D4` residue over a tube around it |
| `PYTHONPATH=src python -m examples.tube_energy_sum` | All five energy residues over that tube, including a group contour through same-side pole coincidences |
| `PYTHONPATH=src python -m examples.spatial_tail` | Wick rotation outside a large sphere, a uniform bound for the exterior integral, and its first two asymptotic coefficients |
| `PYTHONPATH=src python -m examples.null_pinch` | Local integrability and the leading transverse area term over a segment shared by two projected external null edges |
| `PYTHONPATH=src python -m examples.null_insertion_crossing` | Finite-part and delta-derivative terms where an external null edge crosses the insertion pole |
| `PYTHONPATH=src python -m examples.null_endpoint` | A bound on the endpoint energy density and the normalized shrinking-ball coefficient with and without a scaled regulator |
| `PYTHONPATH=src python -m examples.lower_null_endpoint` | The fourth external endpoint on a lower energy branch, including two terms of the normalized shrinking-ball expansion |

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

## Intersecting timelike pinch surfaces

At `p_z=-1`, the two upper poles of `D2` and the insertion meet the
lower pole of `D4` at `(E,p_x,p_y)=(1/3,-1/2,sqrt(7)/6)`. Closing
below, the two pinch surfaces are `u=0` and `v=0` in the `D4=0`
chart `u=D2`, `v=D_insertion`. The spectator denominators there are
`-2/3` and `-4/3`.

The positive imaginary shifts at the `D4` pole give the product of
boundary values `1/(u+i0)/(v+i0)`. Applying the
[Sokhotski–Plemelj formula](https://scipp.ucsc.edu/~haber/archives/physics215_17/Plemelj17.pdf)
in each coordinate gives a double principal value, two single-surface
terms, and a term proportional to `-pi**2*delta(u)*delta(v)`. The last
term contributes `-9*pi**2/(4*sqrt(7))` to the integral of the `D4`
residue with its spatial Jacobian, before multiplying by `-2*pi*i` for
clockwise closure. The example checks the full complex integral over a
small `(u,v)` rectangle at this fixed `p_z`, including the other four
lower poles.

## Closed timelike pinch curve

With `p_z` free, the intersection `D2=D4=D_insertion=0` is the ellipse
`p_y**2 + (41/36)*(p_z+42/41)**2 = 8/41` on the lower `D4` energy
branch. The angular coordinate covers both signs of `p_y`, including
the endpoints where the fixed-`p_z` coordinates degenerate. On the
neighboring ellipses with `(D2,D_insertion)=(u,v)`, integrating the
smooth factor of the residue over the angle gives `2*pi/(2-u)`.
The term supported on the entire intersection is therefore `-pi**3`
in the residue convention, before the `-2*pi*i` closure factor.

The numerical check uses the lower `D4` residue with a common positive
regulator over the tube `abs(u), abs(v) <= 1/320`. Its limiting real
part is the full-curve intersection term. The single principal value
in `u` supplies the imaginary part. This calculation covers one energy
residue and the bounded tube around the entire ellipse.

## Complete energy contour over the tube

Four lower-half-plane poles lie inside a common energy contour throughout
the tube. Summing their residues by contour integration stays regular at
same-side coincidences, including a `D1`–`D3` coincidence on the ellipse.
The fifth lower pole, from `D4`, supplies the boundary value checked in
the preceding example. Adding the two contributions and applying clockwise
closure gives the regulator limit of the real-axis energy integral over
the bounded spatial tube. The example also checks the closure sign and
normalization against direct energy quadrature at two spatial points.
The values use the raw integrand before the overall `1/(i*pi**2)` loop
normalization.

## Large spatial momentum tail

For `r=|p| >= 6` in the finite frame, the spatial centers have norm
at most `3/2` and the energy centers have magnitude at most `1`.
The energy poles therefore lie outside the first and third quadrants,
allowing the energy contour to rotate without crossing a pole.

After rotation to `E=i*t`, write the propagators as `-A_j`. The four
external factors are even in `p_x`. The insertion has its spatial
center at `p_x=-1/2`, so averaging over `p_x` and `-p_x` gives
`-4*p_x**2/(pi**2*prod(A_j)*A_insertion_reflected)` with the loop
normalization included. Each of the six factors has real part at
least `(r**2+t**2)/2`; this gives the integrated bound `14/R**6`
for the region `r>R>=6`, uniformly in the nonnegative regulator.

At zero regulator the large-`R` expansion starts with
`-7/(32*R**6) - 50809/(40960*R**8)`. The expansion uses scaled energy
`t/r` and `1/r`. Reversing the spatial direction and scaled energy
cancels odd powers of `1/r` after integration. The example compares the
coefficients to quadrature at four radii. The ball `r<=R` contains the
finite-momentum pinch regions.

## External null-edge cylinder

Take `p_z=lambda-1` with `1/4<=lambda<=3/8` and a transverse disk
centered at the origin. In this cylinder the opposite-side `D2`–`D3`
and `D1`–`D4` pairs pinch at `E=1-lambda` and `E=-lambda`,
respectively. The pairs have distinct pinch energies over the interval.
The other upper poles are integrated on contours that stay separated
from the pinching poles.

The energy gap of each pinching pair is proportional to
`r_perp**2-i*eps` up to a smooth factor with nonzero real part. Before
reflection the numerator is proportional to `p_x`, so the local
singularity is bounded in magnitude by a constant times `1/r_perp`.
On this reflection-symmetric disk, averaging `p_x` and `-p_x` gives
a factor `p_x**2` and the reflected insertion denominator. This bounds
the averaged energy integral uniformly and gives a cylinder integral
of order the transverse area. Its coefficient after division by the
disk radius squared is approximately `1.990891125602` in the
normalized convention. The example checks the shrinking-disk limit
and the regulator limit at fixed radius. The axial interval avoids
the null-edge endpoints and intersections with the insertion.

## External null edge crossing the insertion pole

Take `p_z=lambda-3/2` and `u=(5*lambda-1)/4` with `abs(u)<=1/32`.
The `D1`–`D4` null pinch crosses the insertion pole at `u=0`,
`E=3/10`. The `D1`–`D2` null pinch also projects into this cylinder
at a different energy. Closing the energy contour below isolates the
lower `D4` and `D1` poles, and a fixed contour groups the remaining lower
poles. The axial integration can then deform into the upper half-plane
for transverse radii up to `1/64`.

After reflection averaging, the two insertion factors approach
`1/(u+i0)**2` in the shrinking-cylinder coefficient. Differentiating the
Sokhotski–Plemelj identity gives the Hadamard finite part and an
`i*pi*delta'(u)` term. With the loop normalization included, the
delta derivative contributes `15*pi*i/2` and the full coefficient is
approximately `-154.611276292+23.561944902i`. The example takes the
regulator to zero at fixed radius and then tests the coefficient as the
radius shrinks. Its axial interval is centered at one crossing.

## External null-edge endpoint

Consider a spatial ball of radius `R` centered at the second external
propagator, `(p_x,p_y,p_z)=(0,0,-1)`. The two incident null edges meet
at `E=1`. At the upper `D2` pole their combined factor is
`D1*D3=2*(p_perp**2-i*eps)`, with the same positive regulator used in all
propagators. Reflection averaging supplies `p_x**2` in the numerator.
The resulting energy density has an integrable `O(1/r)` bound, where
`r` is distance from the ball center. The projected `D1`–`D4` pinch
at a separate energy has bounded density after reflection averaging.
The remaining poles are integrated in groups on fixed energy contours.

At zero regulator the leading density is
`(16/pi)*cos(phi)**2/r`, giving `16*R**2` for the leading normalized
ball integral. The example also checks the finite-regulator scaling
`eps=eta*R**2` at fixed `eta=0.1` against the rescaled endpoint
kernel, and decreases `eta` at fixed ball radius as a separate
check. This calculation covers the endpoint at the second propagator.

## Lower-branch null-edge endpoint

At the fourth external propagator, the two incident null edges meet at
`E=p_z=-1/2` and `p_x=p_y=0`. At the lower `D4` pole the product of
the two regulated null factors is again `2*(p_perp**2-i*eps)`. The
remaining factors at the center are `D2=2` and `D_insertion=1`.
With clockwise energy closure, the reflection-averaged density has
leading term `cos(phi)**2/(pi*r)`. Its ball integral starts at `R**2`.

The `D2`–`D3` pinch projects to the same ball at a different energy.
Its bounded reflection-averaged density contributes to the next term,
alongside the first-order expansion of the lower `D4` residue. The
zero-regulator expansion starts with `R**2 + (1699/150)*R**3`.
The example checks both coefficients, the scaled regulator
`eps=0.1*R**2`, and the clockwise sign against the upper residue sum
at a finite-regulator point.

## Scope of the checks

The CD parameter calculation reproduces the established two-loop coefficient.
The chart checks establish local matching with a transported prescription at
finite regulator. The parameter tail is an estimate in the reduced integral.
Constructing the global cycle on the weighted twistor divisor and controlling
its remaining chart-edge contributions are open steps in the geometric-period
proposal.
