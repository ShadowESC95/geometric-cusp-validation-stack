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
| `PYTHONPATH=src python -m examples.mixed_null_endpoints` | The first and third external endpoints, where a timelike pinch surface adds an imaginary part to each shrinking-ball coefficient |
| `PYTHONPATH=src python -m examples.pinch_inventory` | A pairwise inventory of the finite-frame pinch sets and the small-ball expansion at the isolated insertion center |
| `PYTHONPATH=src python -m examples.second_null_insertion_crossing` | The second insertion crossing, with a real leading cylinder coefficient from an even axial weight |
| `PYTHONPATH=src python -m examples.full_finite_frame` | The complete normalized finite-frame integral via energy residues plus the Wick-rotated exterior, compared with an independent Feynman-parameter calculation and its regulator limit |
| `PYTHONPATH=src python -m examples.finite_cross_ratio` | Symbolic reduction for positive cross-ratios at fixed rho, plus full-domain numerical comparisons at `z=1/2, 3/2, 2` |
| `PYTHONPATH=src python -m examples.complex_cross_ratio` | Continuation of the transported projective integral off the positive axis, and a cut jump from both the logarithm and the moving pole |
| `PYTHONPATH=src python -m examples.cusp_regulator_scaling` | Three regulator paths toward the cusp in a fixed affine frame and the limiting kernel at fixed `eps/rho**2` |
| `PYTHONPATH=src python -m examples.critical_regulator_kernel` | One-dimensional quadrature of the fixed-ratio kernel and its leading logarithmic regulator corrections |
| `PYTHONPATH=src python -m examples.large_regulator_ratio` | Uniform large-ratio expansion of the critical kernel and the independent fixed-regulator cusp coefficient |
| `PYTHONPATH=src python -m examples.matched_regulator_limit` | A common rescaling for the two regulator limits and a proof of their matching coefficient |
| `PYTHONPATH=src python -m examples.regulator_spectrum` | A positive spectral density for the critical kernel, its moment sum rules, and bounds for all positive regulator ratios |
| `PYTHONPATH=src python -m examples.regulator_spectral_edge` | The two leading low-energy density coefficients and the resulting subleading small-regulator terms |
| `PYTHONPATH=src python -m examples.regulator_upper_edge` | Upper spectral exponent, interior cusp, high-order moment asymptotics, and the crossover as the two edge locations merge |
| `PYTHONPATH=src python -m examples.regulator_moment_crossover` | The high-order moment crossover and its first finite-order correction as the edge locations merge |
| `PYTHONPATH=src python -m examples.complex_regulator` | The regulator-plane cut, its two boundary values, and the complex singularity at the upper endpoint |
| `PYTHONPATH=src python -m examples.regulator_interior_cusp` | The square-root correction at the interior spectral cusp, approached along complex rays and both lips of the cut |
| `PYTHONPATH=src python -m examples.weighted_boundary_cycle` | An explicit capped cycle on the complex weighted divisor, its chart orientations, and its scalar pushforward to the benchmark contour |
| `PYTHONPATH=src python -m examples.cd_cycle_seam` | The CD contour difference at the capped-cycle seam and the meridian retained by resolving the moving pole |
| `PYTHONPATH=src python -m examples.cd_cycle_monodromy` | A winding obstruction for periodic CD paths in the fixed positive-x chart and the period of the continued branch |
| `PYTHONPATH=src python -m examples.cd_moving_contour` | A fixed-corner obstruction for moving CD parameters and the period of the nested meridian at chart infinity |
| `PYTHONPATH=src python -m examples.cd_projective_infinity` | A periodic projective family, its signed infinity crossing, and a residue calculation of the lateral jump |
| `PYTHONPATH=src python -m examples.cd_projective_intersections` | A chart-independent signed crossing count, tested by routing the winding through either infinity divisor |
| `PYTHONPATH=src python -m examples.cd_projective_period` | The CD period evaluated with the circle integral first, an absolute seam bound, and its exact finite-width loss |
| `PYTHONPATH=src python -m examples.cd_resolved_seam` | The normal circle at the polar seam, its local cap obstruction, and the excised disk contribution |

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

## Null-edge endpoints on timelike pinch surfaces

The first and third external endpoints also lie on the timelike
`D2`–`D4` pinch surface. Translate to the vertex and reverse both
energy and the axial coordinate for the third endpoint. In these local
coordinates `a=r*mu`, the soft upper pole yields a `1/r` term after
reflection averaging. At the upper pole of the future neighbor, the
past neighbor supplies the factor `4*a-3*delta`, with `delta` the
shift of the future pole away from its axial energy. The common Feynman
regulator places this axial pole below the contour, allowing deformation
of the `mu` interval into the upper half-plane. At leading order its
angular integral is the boundary value of `1/(mu+i0)`, giving
`-i*pi` on the symmetric interval.

The normalized ball integrals divided by `R**2` approach
`16-8*pi*i` at the first endpoint and `1-pi*i/2` at the third. The
example compares both coefficients to the energy-residue quadrature
and checks regulator convergence at a fixed ball radius.

## Finite momentum pinch inventory and insertion center

The opposing branches of distinct propagators can meet where the sum of
the distances to their spatial centers equals their energy separation.
The triangle inequality gives four null segments and two timelike
ellipsoids, all inside `|p|<=9/4`. The ellipsoids are from `D2`–`D4`
and `D4`–`D_insertion`. They meet the null axis at `(-3/2,0)` and
`(-13/10,-1/4)` in `p_z`, respectively. At `p_z=-13/10` the
insertion crosses the interior of the `D1`–`D4` null edge. At
`p_z=-1/4` it crosses the interior of the `D3`–`D4` edge.

The opposite branches of a single propagator meet at its own spatial
center. At the insertion center, the four external factors are
`-1/4,-1/4,-1,1`, and the normalized upper insertion residue starts at
`32/(pi*r)`. For a ball centered there, the leading term is `64*R**2`.
The grouped external residues give a cubic coefficient of about
`-401.838813876`. Their linear variation averages to zero, so the
quartic coefficient comes entirely from the insertion residue and
equals `726`. The example also checks the rescaled soft kernel at
`eps=0.1*R**2`.

## Insertion crossing on the D3–D4 null edge

Set `u=-4*p_z-1` on `abs(u)<=1/32`. The lower `D4` pole meets the
upper `D3` and insertion poles at `E=p_z=-1/4`. At the same spatial
point the lower `D3` pole also meets upper `D2` at `E=1/4`.
The calculation isolates both lower residues and encloses the other
three in one contour. The regulated axial integral deforms onto an
upper semicircle in `u`, including the absolute axial Jacobian `1/4`.

The leading double-pole weight on the first pinch is proportional to
`1/(u**2-1)`. Its derivative vanishes at the crossing, so the
delta-derivative term in the leading coefficient is zero. The finite
part and the other null pair give a normalized cylinder coefficient
of approximately `-63.874350325` after division by the transverse
radius squared. The example takes the regulator to zero at fixed
radius and checks the shrinking-cylinder coefficient separately.

## Full regulated finite-frame integral

Use the common `+i*eps` regulator in all five momentum-space propagators
at `rho=1/2`, `z=1`. Inside `|p|<=3`, integrate the upper energy
residue sum using scrambled Sobol points. Outside this sphere, all
spatial distances to the propagator centers exceed their respective
absolute energy offsets. This allows the exterior Wick rotation. The
two pieces cover the full spatial domain at positive regulator.

Independently, combine the five propagators, shift the energy and spatial
loop variables by their parameter-weighted centroid, and set the insertion
parameter to one. The normalized parameter integrand on this projective
patch is `-4/(A-B*a4-i*eps*P**2)**3`, where
`A=(a1+a2)/4+(1+2*a1)*a3`, `B=1+2*a2`, and
`P=1+a1+a2+a3+a4`. This expression includes the affine numerator
evaluated at the centroid. Rotate `a4` through `+pi/4`. In the swept
sector the denominator has strictly negative imaginary part at
positive regulator. On the final ray, writing `a4=c+i*c`, the real
part is at least `A-B*c`, and the magnitude of the imaginary part is
at least `B*c`. Splitting at `B*c=A/2` bounds the modulus below by
`(A+B*c)/3` for all nonnegative regulators. The resulting integrable
majorant permits the zero-regulator limit on this ray.

At zero regulator the `a4` and `a3` integrals reduce analytically, and
the last two coordinates give
`4*dx*dy/((1+x)*(1+y)*(x+y))` on the positive quadrant. Integrating
over `x` and folding the `y>1` part by inversion leaves
`8*integral_0^1[-log(y)/(1-y**2)]dy=pi**2`. The full-domain spatial
and parametric quadratures are compared at `eps=0.3` and `0.16`.
The smaller-regulator convergence is checked in parameter space.
Scramble spreads describe variation across the chosen deterministic
seeds.

## Full finite-frame integral at non-unit cross-ratios

At fixed `rho=1/2`, use the insertion line with
`a=-1`, `b=-1/2` and `d=1/(2*(1+z))`. Write
`p_j=-(P_insertion-P_j)**2` in the affine frame. This family has
`p1,p2,p3>0>p4` and `-p1*p3/(p2*p4)=z`.
The symbolic check derives the projective denominator and its numerator
from the center separations. It applies to real `z>0`.

The rotated `a4` ray and the regulator bound work with
`A=p1*a1+p2*a2+(p3+2*a1)*a3` and `B=-p4+2*a2`, which are positive.
At zero regulator, integrate out `a4` and `a3`, then set
`x=2*a1/p3`, `y=2*a2/(-p4)`. The remaining integral is
`2*(1+z)*dx*dy/((1+x)*(1+y)*(z*x+y))` on the positive quadrant.
The `x` integral gives `log(y/z)/(y-z)`. On setting `y=z*t` in
the remaining expression and differentiating with respect to `log(z)`,
the derivative is `-2*integral_0^infinity log(s/z)/(1+s)**2 ds`.
The log of `s` integrates to zero by inversion, so the derivative
equals `2*log(z)`. The unit-ratio value fixes the integration constant
to `pi**2`.

The example compares the full spatial and parametric integrals at
`eps=0.3` for `z=1/2, 3/2, 2`, and checks the boundary limit using
smaller regulators in parameter space.

## Complex cross-ratio and the principal cut

Start from the zero-regulator projective chart selected on the positive
axis. Put `a3=w*A/(p3+2*a1)` and
`a4=exp(i*pi/4)*A*k/B`, where `A=(p1*a1+p2*a2)/(1-w)`,
`0<w<1`, and `k>0`. Use the same `x,y` as in the preceding
example. The pullback factors into the two-parameter kernel
`2*(1+z)/((1+x)*(1+y)*(z*x+y))` and a ray factor
`-2*Q/(1-Q*k)**3` with `Q=exp(i*pi/4)`. The ray factor integrates
to one, and the `w` integral has unit length.

The chart coefficients and pullback extend from the positive axis to
the principal slit plane. On compact subsets of that plane,
`abs(z*x+y)` is bounded below by a positive constant times `x+y`
for `x,y>=0`, since the segment between `1` and `z` avoids zero.
The integrable bound makes the pullback holomorphic there. Agreement
on the positive axis identifies the continuation as
`pi**2+Log(z)**2` on the principal sheet.

At `z=-a` with `a>0`, the two lips are
`log(a)**2 +/- 2*pi*i*log(a)`. The jump is `4*pi*i*log(a)`.
In the two-parameter expression the same jump is obtained by
integrating the delta term supported at `y=a*x`. The program compares
both jump calculations and runs the four-parameter quadrature at
complex points off the cut. The construction continues the
zero-regulator projective integral at fixed `rho`. A global loop
contour in weighted twistor space remains to be constructed.

## Regulator scaling in the cusp limit

For real `z>0`, vary the insertion line through `a=-1`, `b=-rho` and
`d=2*rho**2/(1+z)` while keeping the external affine frame fixed.
Its signed separations `p1,p2` vanish as `rho**2`, with `p3,-p4`
remaining positive. Set `a3=rho**2*b3` and `a4=rho**2*b4`.
The numerator and Jacobian supply `rho**6`, which cancels the same
power from the cubed denominator. The common fixed-frame shift in
that denominator is now `-i*eta*P**2`, where `eta=eps/rho**2` and
`P=1+a1+a2+rho**2*(b3+b4)`.

Use the same rotated `b4` ray and parameter map as at finite `rho`.
On the ray `b4=exp(i*pi/4)*A*k/B`, the denominator divided by `A`
has real part at least `1-k/sqrt(2)` and imaginary part at most
`-k/sqrt(2)`. This bounds its modulus below by
`(1+k/sqrt(2))/3` uniformly in the nonnegative regulator. The remaining
weight is proportional to the integrable positive-`z` kernel in
`x,y`. Thus the pointwise regulator limits also hold after integration.

At fixed positive `eps`, `eta` grows without bound and the integral
tends to zero. If `eta` tends to zero, the limit is
`pi**2+log(z)**2`. If `eta` tends to a finite positive value, the
limit retains the shift `-i*eta*(1+a1+a2)**2`. The example tests
these paths at `z=3/2`, using `eps=0.001`, `eps=0.02*rho**2`, and
`eps=rho**4`. This is an order-of-limits result for the common
regulator in the fixed affine frame. The twistor transport example
uses a moving infinity frame and compares the prescribed bracket
shifts under that transport.

## Reduced critical regulator kernel

At fixed `eta=eps/rho**2`, after taking `rho` to zero, the regulator
depends only on `a1,a2`. The rotated `b4` and positive `b3` integrals
can therefore both be done analytically. With `x=2*a1`, `y=2*a2`, the
remaining kernel is

```text
2*(1+z) / ((1+x)*(1+y)*(z*x+y-i*eta*(1+z)*(1+(x+y)/2)**2))
```

It is invariant after integration under `z -> 1/z` by exchange of
`x,y`. Set `t=1/(x+y+2)` and `q=t+(1-2*t)*x/(x+y)`.
Then `q` runs from `t` to `1-t`, and its integrand has denominator
`q*(1-q)*(1-(1+z)*t+(z-1)*q-i*eta*(1+z)/(4*t))`.
Partial fractions give the two logarithms used in the one-dimensional
quadrature.

For fixed positive `z`, write `L=log(1/eta)` and
`C=(1+z)**2*(1+1/z**2)/4`. The leading imaginary correction is
`C*eta*L**2` and the leading real deficit is `pi*C*eta*L`.
To see the coefficients, set `a=eta*(1+z)/4`. Near `t=0` the two
angular endpoints give `log(1/t)` times the difference
`1/(g-i*a/t)-1/g`, with `g=1,z`, multiplied by `2*(1+z)`.
On `a<t<<1` the imaginary part integrates
`a*log(1/t)/(g**2*t)`, giving the squared logarithm. For the real
deficit, setting `t=a*h/g` gives the integral of `1/(1+h**2)` on
the positive axis, which supplies `pi/2`. The difference between
the angular weight and its endpoint logarithms contributes at most
`O(eta*L)` to the imaginary part and `O(eta)` to the real deficit.
The example checks the slow approach of both scaled ratios, as well
as agreement with the four-parameter quadrature at fixed `eta=0.02`.

## Large regulator ratio and the fixed-regulator path

In the folded critical kernel write the non-angular denominator as
`g-i*A`, with `g=1-(1+z)*t+(z-1)*q` and
`A=eta*(1+z)/(4*t)`. For `t<=q<=1-t`, the two endpoint values of
`g` are `1-2*t` and `z*(1-2*t)`. Hence
`0<=g/A<=1/(2*eta)` across the full region, and for `eta>=1` the
inverse has a uniformly convergent geometric series. Its first terms give

```text
i*(8*log(2)-4)/eta + (4-(16/3)*log(2))/eta**2 - i*C3(z)/eta**3
```

The first two coefficients are independent of `z`. The program derives
the third from the angular moment of `g**2` and checks the real and
imaginary remainders using the bound on `g/A`.

For the separate limit with `eps>0` held fixed, do the projective
integration before taking `rho` to zero. Compactify `a1,a2,a3`
with `a_j=b_j/b0`, `b0=1-b1-b2-b3`, and set
`a4=s/b0` on the rotated ray. The Jacobian is `b0**(-5)`.
After dividing by `rho**2`, the resulting kernel is
`-16*p3**2*b0/(A-B*s-i*eps*(1+s)**2)**3`, with
`A=b0*(p1*b1+p2*b2)+(p3*b0+2*b1)*b3` and
`B=-p4*b0+2*b2`. On `s=c+i*c`, the magnitude of the imaginary
part is at least `eps*(1+2*c)`, an integrable majorant for the
inverse cube at any fixed positive regulator. At `rho=0`,
`p1=p2=0` and `p3=-p4=1`, so the coefficient is independent of `z`.

At `eps=0.3` this coefficient is approximately `2.4111i`, with the
quadrature spread reported by the example. Substitution of
`eta=eps/rho**2` into the reduced critical kernel would give
`5.1506i` at the same regulator. Both describe a contribution of
order `rho**2`, with different coefficients for the two limiting
procedures.

## Matched small-regulator and large-ratio limits

Set `lambda=rho**2/eps=1/eta`, `P0=1+a1+a2`,
`U=p3+2*a1` and `V=-p4+2*a2`. Use
`a3=eps*P0**2*alpha/U` and `a4=eps*P0**2*beta/V`.
Rotate `alpha` through `-pi/4` and `beta` through `+pi/4`.
The denominator divided by `eps*P0**2` is

```text
lambda*(z*x+y)/(D*P0**2) + alpha - beta
    - i*(1+eps*P0*(alpha/U+beta/V))**2
```

Here `x=2*a1`, `y=2*a2` and `D=1+z+rho**2*(z-1)`.
Throughout the ray rotations the real part of the squared factor
is at least one. This excludes poles in both deformations. On the
final rays the magnitude of the imaginary part of the denominator
is at least `1+(abs(alpha)+abs(beta))/sqrt(2)`. Its inverse cube is
integrable over the two ray radii. The remaining weight for `eta*F`
is bounded by a constant times
`1/((1+x)*(1+y)*P0**2)` for fixed positive `z` near the joint limit.

The same bound permits both iterated limits and the joint limit
`(eps,lambda)->(0,0)` under the integral. At `eps=0` the two ray
integrals of the inverse cube equal `-1/(2*(source-i))`, where
`source=lambda*(z*x+y)/((1+z)*P0**2)` and the ray differentials are
included. At `lambda=0` the remaining integral gives
`i*(8*log(2)-4)`. Thus `eps*F/rho**2` approaches this coefficient
when both `eps` and `rho**2/eps` vanish. This also matches the
small-`eps` limit of `eps` times the fixed-regulator coefficient.
The example compares the new coordinates with the compactified
integral at positive regulator and checks both boundary paths and
a path with `rho=eps`.

## Positive spectrum of the critical regulator kernel

For `z>0`, write the folded kernel with
`h=4*t*(1-(1+z)*t+(z-1)*q)/(1+z)` and positive weight
`8*t*dt*dq/(q*(1-q))`, with `0<t<1/2` and `t<q<1-t`.
Pushing this weight forward to `h` gives a positive density
`sigma_z(h)` and the Stieltjes representation

```text
F(z,eta) = integral_0^H sigma_z(h)/(h-i*eta) dh
H = max(1,z)/(2*(1+z)).
```

The support follows by maximizing the endpoints in `q`, where `h`
equals `4*t*(1-2*t)/(1+z)` and `z` times that value.
For `z>1`, solving the affine equation for `q` leaves a
one-dimensional integral in `t`. The two quadratic inequalities from
the endpoints determine its allowed intervals. At `z=1`, `h` is
independent of `q` and its two inverse branches in `t` give the
density directly. The exchange `q -> 1-q` gives the reciprocal case.

The inverse moment is `pi**2+log(z)**2`. The zeroth, first and second
moments give the three large-ratio coefficients in the preceding
examples. Let `M0` and `M1` be the zeroth and first moments.
For every `eta>0`, `Re F` is positive and decreases. Both
`eta*Im F` and `eta**2*Re F` increase, with upper bounds `M0`
and `M1` and lower bounds `M0/(1+(H/eta)**2)` and
`M1/(1+(H/eta)**2)`, respectively. These follow by integrating
`eta**2/(h**2+eta**2)` with the positive density and its
`h`-weighted counterpart. The example checks the density against
the existing folded quadrature and verifies the four moment sum rules.

## Low spectral edge and subleading regulator terms

The spectrum near `h=0` comes from the branches near `t=0` and
`t=1/2`. On the small-`t` branch, invert the defining quadratic in
`t` at fixed `q`. With `g=1+(z-1)*q`, it gives
`t=(1+z)*h/(4*g)+O(h**2)`. The angular endpoints occur at
`q=(1+z)*h/4+O(h**2)` and `1-q=(1+z)*h/(4*z)+O(h**2)`.
The pushforward weight divided by `h` therefore has a factor
`(1+z)**2/2` times `1/(q*(1-q)*g**2)` at leading order.
Subtract the linear interpolant of `1/g**2` to isolate the two
endpoint logarithms and integrate the resulting regular remainder.

On the other branch, set `t=1/2-h*v`. The angular variable tends to
`1/2`, and for `z>1` the width in `v` tends to
`(1+z)*(1-1/z)/4`. Substituting into the pushforward integral gives
an additional `2*(1+z)**2/z` in the coefficient of `h`.
The two inverse branches at `z=1` give the continuous limit of the
same formula. Together,

```text
sigma_z(h) = h*(A*log(1/h)+B) + O(h**2*log(1/h))
A = (1+z)**2*(1+z**(-2))/2
B = (1+z)**2/2 * (
    (1+z**(-2))*log(4/(1+z)) + (2*z**(-2)-1)*log(z)
    + 6/z - 1 - z**(-2)
)
```

In the transform at positive `eta`, put `L=log(1/eta)`.
The imaginary part has expansion `eta*(A*L**2/2+B*L+O(1))).
The real deficit from the boundary value is
`(pi/2)*eta*(A*L+B)+o(eta)`. These follow by inserting the density
expansion into the real and imaginary kernels. The constant
`pi/2` is the integral of `1/(1+s**2)` on the positive half-line,
and the integral of `log(s)/(1+s**2)` vanishes by inversion.
The example compares the predicted `B` against both the density
and the real deficit, and checks that subtracting the imaginary
logarithms leaves a finite remainder.

## Upper spectral edge and its merger at z=1

Write `H_minus=min(1,z)/(2*(1+z))`, `H_plus=max(1,z)/(2*(1+z))`,
and `Delta=H_plus-H_minus`. Use reciprocity to take `z>=1`.
The exact edge coordinates obey
`H_plus-h=H_plus*(4*t-1)**2+8*t*Delta*(1-t-q)`.
For `z>1` the maximum occurs at `t=1/4,q=3/4`. The allowed
interval in `t` has width `sqrt(H_plus-h)/(2*sqrt(H_plus))`,
and the angular weight after integrating the delta function
tends to `16/(3*Delta)`. This gives
`sigma(h)~8*sqrt(H_plus-h)/(3*Delta*sqrt(H_plus))`.
Below `H_minus`, an interval opens in the excluded region near
`t=q=1/4`. Its width gives a one-sided square-root deficit with
coefficient `8/(3*Delta*sqrt(H_minus))`.

At `z=1` the whole angular interval at `t=1/4` attains
`H_plus=1/4`. The two inverse branches give
`sigma(h)~2*log(3)/sqrt(H_plus-h)`. If the density has edge power
`p` and coefficient `C`, its normalized moment
`integral (h/H_plus)**n * sigma(h) dh` is asymptotic to
`C*H_plus**(p+1)*Beta(n+1,p+1)`. Thus the normalized moments
decay as `n**(-3/2)` away from `z=1` and as `n**(-1/2)` at
`z=1`. The program checks this prediction up to moment order 1024.

To connect the two regimes, keep `xi=(H_plus-h)/Delta` fixed while
`Delta` shrinks, and use `v=(4*t-1)*sqrt(H_plus/Delta)`.
The allowed positive interval becomes
`sqrt(max(xi-1,0))<=v<=sqrt(xi)`, and
`q=3/4-(xi-v**2)/2` in the limit. The rescaled density
`sqrt(Delta)*sigma(h)` tends to the integral of
`1/(q*(1-q))` over that interval after including the equal
contribution from negative `v`. The profile behaves as
`(16/3)*sqrt(xi)` at small `xi` and
`2*log(3)/sqrt(xi)` at large `xi`. This recovers the two endpoint
powers and resolves the interior cusp at `xi=1`.

## High-order moments across the merging spectral edges

Write `M_n=integral (h/H_plus)**n * sigma(h) dh` and use
`delta=(H_plus-H_minus)/H_plus`. With `q=1-t-u*(1-2*t)`,
the exact factorization is

```text
h/H_plus = (1-(4*t-1)**2)*(1-delta*u)
8*t*abs(dq/du) = 1-(4*t-1)**2
```

Keep `lambda=n*delta` fixed and set `s=sqrt(n)*(4*t-1)`.
The two factors become `exp(-s**2-lambda*u)` at leading order.
Let `q0=3/4-u/2`, `c=(2*u-1)/4`, and `f(q)=1/(q*(1-q))`.
The resulting expansion is `sqrt(n)*M_n=G(lambda)+J(lambda)/n+o(1/n)`,
where

```text
G(lambda) = sqrt(pi)/4 * integral_0^1 exp(-lambda*u)*f(q0) du
J(lambda) = sqrt(pi)/4 * integral_0^1 exp(-lambda*u) * (
    c**2*f''(q0)/4 - (7/8+(lambda*u)**2/2)*f(q0)
) du
```

The order `1/sqrt(n)` term integrates to zero because it is odd in
`s`. The Gaussian moments of `s**2+s**4/2` give the factor `7/8`
in the correction. The profile satisfies `G(0)=sqrt(pi)*log(3)` and
`lambda*G(lambda)->4*sqrt(pi)/3` at large `lambda`. The first limit
recovers the `n**(-1/2)` law. The large-`lambda) coefficient matches
the `n**(-3/2)` edge law with `delta>0).

The example compares both terms against the original folded-coordinate
integral and checks order 1024 independently against the spectral
density. The expansion with the stated correction holds for fixed
`lambda` as `n` grows.

## Complex regulator cut and endpoint singularity

Set `w=i*eta` and write the critical kernel as the transform
`G(w)=integral_0^H sigma(h)/(h-w) dh`. It is analytic off `[0,H]`,
so the regulator cut in `eta` is the segment from `0` to `-i*H`.
At an interior point `v`, the boundary values in the `w` plane are
the principal value plus or minus `i*pi*sigma(v)`. This uses the
[Sokhotski–Plemelj formula](https://scipp.ucsc.edu/~haber/ph215/Plemelj18.pdf)
with the denominator `h-w`. The example subtracts the value of the
density at the pole and adds the principal value of the constant
term analytically.

For `w=H+delta`, take the principal square root of `delta`.
At `z=1`, the edge coefficient `C=2*log(3)` gives
`G(H+delta)~-pi*C/sqrt(delta)`. For `z!=1`, `G(H)` is finite
and `G(H+delta)-G(H)~pi*C*sqrt(delta)` with `C` from the upper
edge example. The signs follow from the negative denominator outside
the support. In the second case the subtracted kernel is
`delta/((H-h)*(H-h+delta))`. In both cases inserting the leading
density and scaling `H-h` gives the convergent integral of
`s**(-1/2)/(1+s)`, equal to `pi`. The quadrature checks the
coefficients along three rays in the right half-plane of `delta`.

## Interior spectral cusp on the complex-regulator cut

For fixed positive `z!=1`, let `c=H_minus`, and use the transform
`G(w)=integral sigma(h)/(h-w) dh`. The density is continuous at `c`.
The interval excluded just below `c` contributes a deficit
`-C*sqrt(c-h)` with `C=8/(3*(H_plus-H_minus)*sqrt(c))`.
In a neighborhood of `c`, the rest of the density has a smooth
continuation across `c`, and the next term in the deficit is
`O((c-h)**(3/2))`.

With `delta=w-c`, write `s=c-h` in the deficit term.
Subtracting its value at `delta=0` leaves
`-C*delta*integral_0^b ds/(sqrt(s)*(s+delta))` for a fixed
local cutoff `b>0`. Extending this integral to infinity changes
the result by `O(delta)` and gives `-pi*C*sqrt(delta)`.
The other terms change by `O(delta)` on each side separately.
Thus the upper and lower boundary values at `c` are finite,
with the same square-root coefficient relative to their respective
limits. The square root is principal, so the correction is real
on the right lip and purely imaginary with opposite signs on the
two left lips. The example verifies these phases and coefficients
with four complex rays and direct boundary-value quadrature.

## Capped weighted-divisor cycle and scalar pushforward

In the weighted blow-up of `(b,c,d)` with weights `(1,1,2)`, the
`d` chart has the cover `b=r*u, c=r*v, d=r**2` with residual
`(u,v)->(-u,-v)` on the exceptional divisor. The cross-ratio is
`z=u*v/a-1`. In radial-last convention the residue on the cover is
`2*da^du^dv/(a*(a-u*v))`. The code checks it against the `b=-rho`
chart using `D=1/u**2`, including the Jacobian sign. It also checks
the canonical form under the change to columns `(3,2)`. The new
weight-two coordinate is `-z*r**2`, so the weighted chart stays
valid along the chosen contour.

Choose `a=R*exp(i*alpha)` with `R>0` and
`z=exp(i*phi)` with `-pi<=phi<=pi`. For `z!=-1`, take
`u=sqrt(abs(a*(z+1)))*exp(i*theta)` and `v=a*(z+1)/u`.
One turn downstairs uses `0<=theta<=pi`. Both ends in `phi`
send the whole theta circle to `u=v=0` at each `a`. This defines
a closed cycle as the image of the `a` circle times the suspension
of the theta circle. The caps lie on the quotient locus. The
form on the cover is regular there for `a!=0`, and the construction
avoids the polar set `a*(a-u*v)=0`.

With orientation `(alpha,phi,theta)), the pullback on the open
cylinder is `2*da/a^dz/z^du/u`. The factor two cancels the half
turn in the cover. Integrating the two oriented fiber circles with
normalization `(2*pi*i)**(-2)` pushes the cycle to `dz/z`.
The pullback is bounded in angle coordinates, so the caps add no
separate term when integrating by a limit of truncated cylinders.
The principal-sheet CD weight has matching limits at `z=-1`.
The example integrates it using the reduced projective quadrature
and obtains the normalized period `2*pi**2/3`. The conversion to
the two-loop coefficient uses the published contour normalization
in Eqs. (6.2)–(6.6) of Multi-Loop Negative Geometries.

This is an explicit complex cycle for the local weighted residue
and an externally supplied CD scalar weight. On the Lorentzian
real slice `c=-2*conj(b)` with `a,d` real, the cross-ratio is real.
The calculation supplies no deformation of that prescription to
the capped cycle, no transport of the CD integration class along
the cycle, and no accounting of additional chains in such a
deformation.

## CD parameter seam and resolved meridian

On the transported projective domain, the auxiliary ray and interval
integrate to one. In the two remaining CD parameters, set `y=x*t`.
The relative integration form at fixed `z` has coefficient
`2*(1+z)/((1+x)*(1+t*x)*(t+z))` in `dx^dt`.
At `z=-a` with `a>0`, the upper and lower limits bypass the pole
`t=a` on opposite sides. The upper left-to-right bypass minus the
lower bypass is one clockwise meridian. The other pole
`t=-1/x` stays outside that meridian for `x>0).

The jump after dividing out `1+z` is
`-4*pi*i*integral_0^infinity dx/((1+x)*(1+a*x))`.
The integral is `log(a)/(a-1)` with value one at `a=1`.
Thus the scalar jump vanishes at `z=-1`, while the winding is
still minus one and the limiting jump after division is
`-4*pi*i`. The numerical check integrates the two detours with
their orientations and compares their difference with the residue.

To track the pole as the outer cycle enters the cap, write
`s=z+1` and `t=1+s*v`. This chart keeps the pole at `v=-1`.
Meridians of size `O(abs(s))` in the original `t` coordinate
retain a full winding in the resolved coordinate. The pulled-back
relative form is
`2*s*dx^dv/((1+x)*(1+x*(1+s*v))*(1+v))`.
It has zero scalar period in the cap limit, with limiting coefficient
`2/((1+x)**2*(1+v))` after division by `s`. This identifies the
meridian that must be accounted for when gluing the projective
CD parameter chains. A global chain in the original loop variables
remains open.

## CD contour winding obstruction and monodromy

Fix `x>0` and transport the `t` path from zero to infinity as
`z` goes counterclockwise around zero. A periodic family of paths
avoiding `t=-z` and `t=-1/x` would make
`(t+z)/(t+1/x)` a homotopy of loops in the punctured plane.
At `t=0` the loop is `x*z` and winds once. At infinity it is
constant and has winding zero. The ratio extends at infinity, so
this discrepancy rules out such a periodic family. At a collision
of the two poles the ratio is taken on their common complement.
The obstruction applies to the fixed positive-`x` path ansatz.

Across the seam the continued branch carries the clockwise
meridian of the moving pole. On the lower semicircle its scalar
period is

```text
-4*pi*i*(1+z)*integral_0^infinity dx/((1+x)*(1-z*x))
  = 4*pi*i*Log(-z).
```

The principal value of `Log(-z)` in this formula has argument
`arg(z)+pi` on the lower semicircle. Adding this period to the
principal-sheet integral gives the branch where `Log(z)` has
increased by `2*pi*i`. Approaching the starting point `z=1` from
below gives a one-turn shift of `-4*pi**2`. The pole collision at
`x=1,z=1` is treated by this lateral limit.

The principal-sheet contour mean is `2*pi**2/3`. If one starts at
`z=1`, continues counterclockwise through the seam, and keeps the
returned sheet on the lower semicircle, the meridian contributes
`-pi**2` to the mean, giving `-pi**2/3` for that continued path.
This distinguishes the two prescriptions even though the seam
values agree. Any identification with a Lorentzian cycle must
explain the branch reset or supply compensating chains. CD chains in
which the `x` contour also moves remain open.

## Moving CD square and nested meridian

Put `X=x/(1+x)`, `Y=y/(1+y)`. The two-parameter form is
`2*(1+z)*dX^dY/D` with `D=z*X*(1-Y)+Y*(1-X)`. This sends the
positive quadrant to the unit square. Its opposite corners
`(1,0)` and `(0,1)` have `D=z` and `D=1`, respectively.
Suppose the square moves continuously and periodically with `z`
in this affine chart, keeps these corners fixed, and avoids
`D=0` on its interior, which includes a connecting arc between
the two corners. Remove the two other corners if the original
corner singularities are treated as relative boundaries.
The image under `D` of a connecting arc would give a homotopy
of the corner loops in the punctured complex plane. Their
windings are one and zero. This rules out the periodic affine
square with these corners even when both parameters move.

The moving divisor is `Y=-z*X/(1-(1+z)*X)`. For `z` off the
seam, it reaches `Y=infinity` at `X=1/(1+z)`. In `U=1/Y`
the form has a logarithmic factor `dU/U` at chart infinity.
The `Y` residue of the moving divisor is
`-2*dX/(X-1/(1+z))`. A positive `X` meridian about that pole
with a clockwise `Y` meridian about the moving divisor therefore
has period `-8*pi**2`, with the orientation `(X,Y)`.
The example checks this by nested contour quadrature at three
points of the lower semicircle and two radii. It is the constant
shift obtained by adding `2*pi*i` to `Log(-z)` in the
moving-pole period. A fixed integer number of these nested
meridians leaves the angle dependence of that period intact.

A general relative chain can change the boundary construction
or involve contributions at infinity. Accounting for those
contributions and relating the result to the Lorentzian contour
remain open.

## Projective CD family and the infinity crossing

Take `x,u>0` and `y=z*u`, and use projective coordinates at
the ends. This family is periodic as `z` winds and has
`zx+y=z*(x+u)` on its interior. It reaches the polar divisor
`y=-1`, or `Y=infinity`, at `z=-1,u=1` on the unit circle.
On the diagonal connecting the fixed opposite corners use
`x=(1-v)/v` and `u=v/(1-v)`. The affine divisor function
is then `D=z*((1-v)**2+v**2)/(1-v+z*v)`.
At `v=1/2,arg(z)=pi` the chart factor has local derivatives
`-2,-i/2` in the orientation `(v,arg(z))`. Its index is
`+1`, so `D` has a pole of index `-1` there. That pole
accounts for the difference between the corner windings.

The pulled-back CD form is
`2*(1+z)*dx^du/((1+x)*(1+z*u)*(x+u))`.
After integrating over `x`, its coefficient in `du` is
`2*(1+z)*log(u)/((u-1)*(1+z*u))`, with a continuous
logarithmic quotient at `u=1`. This integral agrees off the
cut with the principal-sheet CD parameter value. For `z=-a`
the pole at `u=1/a` approaches from above when `Im(z)>0`.
With `x` fixed, its residue is
`-2*(1-a)/((1+x)*(1+a*x))`. Integrating the residue with
the positive orientation gives the lateral jump
`4*pi*i*log(a)`. The crossing on the unit circle has zero
scalar jump at `a=1`; the chart index is still nonzero.

The periodic geometric family crosses the polar locus.
Turning it into a pole-avoiding relative cycle at the seam
and identifying its Lorentzian class remain open.

## Projective CD intersection balance

In homogeneous coordinates `X=X1/X0`, `Y=Y1/Y0`, write the
moving divisor as the zero set of the section
`B=z*X1*(Y0-Y1)+Y1*(X0-X1)`. The ratio
`D=B/(X0*Y0)` is a globally defined meromorphic function.
The three divisors in this balance are `B=0`, `X0=0`,
and `Y0=0`.

Take a periodic family of connecting arcs oriented from `(1,0)`
to `(0,1)`, with `z=exp(i*phi)` and cylinder orientation
`(arc parameter,phi)`. For transverse isolated crossings away
from divisor intersections, let `I_B,I_X,I_Y` count the signed
crossings of the three divisors. Excise a small disk at each
crossing. On the resulting punctured cylinder, `D` is nonzero
and finite. The two boundary windings are `+1` and `0`,
and a small positively oriented circle about a crossing measures
its zero index for `B` or minus its index for either infinity
divisor. Summing the boundary windings gives
`I_B-I_X-I_Y=-1`. With all three divisors avoided, the
boundary loops would be homotopic through nonzero finite values
of `D`, which is impossible with the fixed endpoints.
This also covers a periodically deformed square if it contains
such a connecting arc.

For examples use `x=2*z**m*(1-v)/v`,
`y=3*z**n*v/(1-v)`, with `0<v<1` and integer `m,n`.
The signed counts are `(n-m-1,-m,n)`. The moving divisor
crosses where the magnitudes of `z*x` and `y` coincide;
factoring out `z**n` leaves relative phase `z**(m+1-n)`.
As `v` increases the magnitude of the first term decreases
relative to the second, giving the minus sign in its index count.
The `x` magnitude decreases and the `y` magnitude increases,
giving the signs of the two infinity counts. The program measures
each index on a small loop and checks their balance for six
families, including routing the crossing exclusively through
either infinity divisor. These indices constrain relative
deformations with the fixed endpoints. A changed boundary or
the original Lorentzian loop cycle still requires a separate
construction.

## Circle-first CD period and seam excision

On the family `y=z*u` with `x,u>0` and `abs(z)=1`,
the normalized integrand is the pulled-back CD form
times `dz/(2*pi*i*z)`. Integrate over `z) first. The scalar factor
`(1+z)/(z*(1+u*z))` has residue one at zero. The second pole at
`-1/u` contributes `1/u-1` when `u>1` and lies outside the circle
when `u<1`. At `u=1` the apparent second pole cancels in this scalar
factor. The resulting circle weight is one on `0<u<1` and
`1/u` on `u>1`.

Integrating `1/((1+x)*(x+u))` in `x` gives
`J(u)=log(u)/(u-1)`. Since `J(1/u)=u*J(u)`, both halves
give the same contribution. The full normalized period is
`4*integral_0^1 J(u) du=2*pi**2/3`, using the series for
`1/(1-u)`. This derivation integrates the base circle before
using the CD scalar benchmark.

For absolute convergence use
`abs(1+u*z)**2=(1-u)**2+u*abs(1+z)**2`.
The ratio of the two moduli is bounded by `1/sqrt(u)` for
the numerator `1+z`. The integral of `J(u)/sqrt(u)` over
positive `u` is `pi**2`, by folding at one and setting
`u=v**2` on the lower interval. Thus the full normalized
absolute integral is at most `2*pi**2`, and omitting the arcs
of angular half-width `delta` at the seam costs at most
`2*pi*delta` in absolute value. This bound justifies the
exchange of integrations and shows the seam excision vanishes.

For a retained arc `-a<=arg(z)<=a` with `a<pi`, the base
weight is `(a+(1-u)*arg(1+u*exp(i*a))/u)/pi`, using the
continuous argument along the arc. Integrating this weight gives
an omitted scalar contribution of `delta**2-delta**3/(3*pi)`
at `a=pi-delta`, in agreement with direct scalar quadrature.
This evaluates the scalar period of the excised family. The
chain at the seam and its identification with the Lorentzian
integration class remain open.

## Resolved CD seam and normal disk

On the real sweep `u>0, |z|=1`, set `f=1+u*z=r*exp(i*theta)`.
In a neighborhood of the crossing, the inverse is
`u=abs(f-1), z=(f-1)/u`. Resolving the normal coordinate retains
the `theta` circle at `r=0` over each positive `x`.
The normal circle has winding `+1` in the positive disk
orientation. The boundary of the complement of the excised
disk has the opposite orientation. The corresponding loop in
`z` has zero winding around zero.

At fixed `z`, the coefficient of `dx^df/f` at `f=0` is
`2*(1+z)/((1+x)*(z*x-1))`. The seam value is zero, and the
coefficient divided by `1+z` tends to `-2/(1+x)**2`.
For a positive normal meridian followed in the coordinate
order `(x,theta)`, this yields the stripped period `-4*pi*i`.

For the excision in the real sweep, `du*dphi=(r/u)*dr*dtheta`
and `1+z=i*r*sin(theta)+O(r**2)`. The angular integral of
`i*sin(theta)*exp(-i*theta)` is `pi`. Including the normalized
base measure and the `x` integral gives a contribution
`R**2/2+o(R**2)` from the disk `abs(f)<R`. This is a normal
disk excision, with its own cutoff distinct from the angular
excision in the preceding example.

The ratio `(1+y)/(zx+y)` has winding one on the normal circle
and extends to one at `y=infinity` for finite `x,z`.
A disk with finite `x,z` avoiding both divisors would extend
that loop in the nonzero complex numbers, which is impossible.
In the limiting slice `z=-1`, the exterior disk
`k=R/(1+y)` contains `y=infinity` at `k=0` and meets
the moving divisor `y=x` at `k=R/(1+x)`. The winding and
the vanishing scalar excision term therefore give separate
information. A relative cap with additional boundary pieces
and its connection to the Lorentzian contour remain open.

## Compactified CD seam and the two base poles

Use homogeneous coordinates `[x0:x1]`, `[y0:y1]`, `[z0:z1]`
with affine ratios `x=x1/x0`, `y=y1/y0`, `z=z1/z0`.
Write `B=z1*x1*y0+z0*x0*y1`, `A=x0+x1`, `C=y0+y1`.
The quotient `D=B/(z0*A*C)` extends as a nonzero function
where its zero and pole divisors are excluded, including across
the parameter infinities `x0=0` and `y0=0`. On the resolved
seam circle it winds `-1`. This rules out an absolute cap in
that complement, even after compactifying `x` and `y`.
The quotient `D/z=B/(z1*A*C)` gives the same winding.

An explicit cap exists when both base poles are allowed.
At any positive `x`, first contract the small loop in `z` to
`z=-1` with `y` fixed. In that slice the coordinate
`w=(y+1)/(y-x)` identifies the twice-punctured `y` sphere with
the nonzero complex numbers. The positive `y=-1` meridian
and the negative `y=x` meridian have the same winding in `w`.
The example interpolates their degree-zero factors in logarithmic
coordinates. It then transfers the negative meridian to the
`z` plane while keeping `z*x+y=q` fixed on the loop. At the
end `y=x` and `z=-1+delta*exp(-i*theta)`, with `0<delta<1`.
The cap `z=-1+delta/k`, `abs(k)<=1`, has homogeneous coordinates
`[k:delta-k]`. Its moving section is the nonzero constant
`x*delta`. It crosses `z=0` at `k=delta` and `z=infinity` at
`k=0`, each with index `+1`.

These signed intersections are forced for any cap that avoids
`B=0`, `A=0`, and `C=0`. The boundary windings of `D` and
`D/z` give `-I_infinity=-1` and `-I_zero=-1`, respectively.
The original meridian also bounds a local disk meeting only
`C=0` once. Gluing a cap to the reversed local disk gives a
closed two-cycle of degrees `(0,-1,1)` in the three projective
factors, since `B` has multidegree `(1,1,1)`, and `A` and `C`
measure the first two degrees. Both base intersections measure
the third degree.

The normalized CD form has genuine base residues, with
the overall `1/(2*pi*i)` omitted here. They are
`2/((1+x)*(1+y)*y)` at zero and
`-2/((1+x)*(1+y)*x)` at infinity in the variable `1/z`,
in the coordinate order `(x,y,base)`. This construction gives
a compactified cap at each finite positive `x`. A relative
three-chain with prescribed base-pole and `x` endpoint data
and its Lorentzian identification require further work.

## Swept CD cap and relative boundary terms

For a closed interval `I=[a,b]` strictly inside the positive
`x` axis, sweep the cap of the preceding example over `I`.
Use `f=R*exp(i*theta)` and `base=(f-1)/abs(f-1)` for its
seam boundary. The disk in `k` is punctured at `k=0` and
`k=delta(x)`, where `delta(x)=R/(1+x)`. Remove disks of radii
less than `delta(x)/3`. Their positive normal circles are
`L_infinity` and `L_zero`. The three collar joins cancel with
the disk junction. For the resulting oriented cap `C_x`,
`d C_x = seam - L_zero - L_infinity`. With `x` first, the
swept boundary is `C_b-C_a-I x seam+I x L_zero+I x L_infinity`.
Thus the cap supplies the seam face with the opposite sign
to the excised sweep, and leaves two pole tubes and two end caps.

Only the first collar contributes to the normalized period.
On the second collar `z=-1`, so `dz` pulls back to zero.
On the third, `y=x+(1-s)*q` and `z=-1+s*q/x`. In the
coordinate order `(x,q,s)`, the Jacobian is `q/x`, which
cancels the moving pole. The remaining coefficient is analytic
in the entire small `q` disk for every `x>0`; its circle
integral vanishes. The punctured disk has `y=x`, hence
`dx^dy=0`, including in the reciprocal base chart.
The form has identically zero pullback to the punctured disk,
and the corresponding residue pullbacks along the base-pole
intersection curves vanish as well. The tubes remain in the
oriented boundary.

The first collar uses `y=f-1`, `z=(1-s)*base-s` in the
order `(x,theta,s)`. Deform the filled normal disk from
`s=0` to `s=1`. Across `f=0` its pullback is locally bounded,
since `1+z=O(abs(f))` and `dy` supplies the polar Jacobian.
For Stokes' theorem, first remove an inner circle. Its
collar integral tends to zero quadratically. The deformed
disk at `s=1` and the two constant-`x` faces have zero
pullback. This gives equality of the first collar period
and the original normal disk period on every finite `I`.
The example also compares independent quadratures of the two.

In the coordinate `t=x/(1+x)`, the first collar density is
`-(1-s)*(1+base)**2/(pi*z*(z*t+(f-1)*(1-t)))`. For `R<=.2`,
`abs(1+base)<=2*R/(1-R)` and the real parts of `z` and
`f-1` are at most `-.5` and `-.8`, respectively. This
bounds the density by `25*(1-s)*R**2/pi`, uniformly to
both ends of the `t` interval. The full scalar cap limit
therefore exists, with period `R**2/2+o(R**2)` as `R`
shrinks. The endpoint chain requires its own treatment.
On the last disk the homogeneous moving section is
`delta*t*(1-t)`, so its closure meets the moving divisor
at the endpoints. The scalar bound does not resolve
those faces or identify this relative chain with the
Lorentzian contour.

## Resolved endpoint faces of the swept CD cap

Use a fixed transferred base radius `delta=R` in the last
two pieces of the cap. This is a deformation of `delta=R/(1+x)`
for every positive `x`. Both lie in `0<delta<=R`, so the
transfer has the same pole avoidance and the disk excisions
can follow the base poles. The first collar is unchanged. In
the last disk this choice keeps the two punctures `k=0` and
`k=delta` distinct as `x` tends to either endpoint.

Blow up `(x,y)=(0,0)` over the base and use the charts
`y=x*a` and `x=y*b`. The strict moving section in the
two charts has affine factor `z+a` and `1+z*b` respectively.
In the order `(x,a,z)`, the Jacobian `x` cancels the radial
factor of the moving divisor. The pulled-back form has no
radial residue at `x=0`. At the other endpoint, set
`X=1/x`, `Y=1/y` and use `Y=X*a`. The analogous strict moving
factor is `1+z*a`. Here the Jacobian cancels `X`, and the
radial residue is again zero. The base divisors are still
punctured throughout this construction.

The first collar approaches the proper endpoint fiber with
`y=f-1` bounded away from zero and infinity. On the logarithmic
collar at `z=-1`, split the interpolation in `w=(1+y)/(y-x)`
into two radial segments and the phase rotation on
`abs(w)=x**(-1/2)`. On that circle both coordinates in the
joining chart tend to zero at rate `O(sqrt(x))` near the
lower endpoint and `O(1/sqrt(x))` near the upper one. For
example, at `x=e**2` the bounds for `y` and `b=x/y` are both
`e*(1+e)/(1-e)`, using `y=(x*w+1)/(w-1)`. At `x=e**(-2)`
the same bounds hold for `Y=1/y` and `X/Y`. The phase
rotation therefore collapses to the common chart point.
The two radial legs limit to a leg on the proper endpoint
fiber and a leg on the exceptional curve. Their homogeneous
coordinates keep `y=-1` and the strict moving divisor excluded.

For the last transfer use `a=1+(1-s)*q`, `z=-1+s*q`,
`q=delta*exp(-i*theta)` at the lower end. Then `z+a=q`.
At the upper end the ratio is reciprocal, giving the
affine factor `1+z/a=q/a`. The last disk has ratio one
and base coordinates `[k:delta-k]`, so its strict moving
section is the constant `delta`. These pieces extend
over both exceptional faces with the base punctures
removed. The lifted faces are relative boundary data
on the resolved space; their matching to the endpoint
and base boundary of the Lorentzian chain remains open.

## Reduced principal contour and the shared CD meridian

Take the positive-parameter representative of the reduced
principal contour, with `x,y>0` and `z=exp(i*phi)`. On any
retained arc `abs(phi)<=pi-delta`, deform `y` to `exp(i*s*phi)*y`
for `0<=s<=1`. The two relevant factors are `1+y*exp(i*s*phi)`
and `z*(x+y*exp(-i*(1-s)*phi))`. The bound
`abs(a+b*exp(i*alpha))>=sin(delta/2)*(a+b)` for positive `a,b`
and `abs(alpha)<=pi-delta` keeps both factors nonzero. This
gives a deformation over the retained arcs, with side faces
over their endpoints.

In the positive-parameter representative, write `u=y/x` and
`h=u+z`. Its seam is `u=1, z=-1`; the Jacobian from
`(Re h, Im h)` to `(u,phi)` is negative. A positively oriented
normal circle is therefore `h=R*exp(-i*theta)`. On the
first cap collar take `z_s=(1-s)*z-s`, `y_s=x*(h-z_s)`.
This keeps the moving factor equal to `x*h`, and `Re(y_s)>0`
for the small radii used here. At `s=1` one obtains `z=-1`,
`y=x*(1+h)`. This is the same clockwise moving meridian
used by the rotated cap with fixed transferred radius `R`.
The shared transfer collar and punctured base disk give the
same two base-pole tubes in these reduced representatives.

The two normal-disk periods agree exactly. In the rotated
representative use coordinates `(x_rot,u_rot,z_rot)` with
`y_rot=z_rot*u_rot`. Under `x_rot=1/x`, `u_rot=u`, `z_rot=1/z`,
the normalized three-form pulls back to the positive-parameter
form in `(x,u,z)`. The two inversions preserve the product
orientation, and the rotated normal coordinate is `f=h/z`.
Since `abs(z)=1`, the map takes the `h` disk to the `f` disk
with the same radius. It interchanges the two `x` endpoints.
The example checks this identity symbolically and compares
independent quadratures of the two disks and the principal
cap collar. The collar equals the removed principal disk by
the same punctured-disk Stokes argument as for the rotated cap.

This comparison is in the reduced CD parameter space. The
projective reduction of the scalar integrals gives these
forms. Lifting the matching to the unreduced Feynman chain,
including the auxiliary parameter boundaries and both resolved
endpoint faces, remains a separate step.

## Unreduced CD auxiliary faces and diagonal-cap corners

At zero regulator the unreduced projective map factors the
five-form into the CD three-form and an auxiliary two-form.
For `C0=p1*a1+p2*a2`, set `t=k/(1+k)`, `r=w*(1-t)`,
`s=t`. Write `A=1-r-s`, `B=r`. The auxiliary domain is the
oriented triangle `r,s>=0, r+s<=1`. The projective coordinates
are `a3=C0*B/(A*(p3+S*a1))` and
`a4=C0*Q*s/(A*(-p4+T*a2))`, with `Q=exp(i*pi/4)`.
The normalized fiber form has density
`-2*Q/(1-(1+Q)*s)**3` in the order `(r,s)`. Its auxiliary
denominator stays at least `sin(pi/8)` away from zero,
and the two-form integrates to one by taking the endpoint
difference of `-(1-s)**2/(1-(1+Q)*s)**2`. The example checks
the uniform bound on errors from trimming all three edges.

The three auxiliary faces correspond generically to `a4=0`,
`a5=0` (the projective infinity edge), and `a3=0`, in the
counterclockwise boundary order. The pulled-back fiber form
vanishes on each one-dimensional face. For a three-chain
`K`, its lift has boundary `dK x triangle - K x d(triangle)`.
The auxiliary faces therefore remain part of the lifted
boundary even though integrating the fiber reproduces the
reduced scalar exactly.

On the diagonal cap `y=x`, put `h=1+z`, `d=5*z+3`.
One homogeneous representative of the five Feynman coordinates is
`[2*h*x*(1+x)*A, 2*h*x*(1+x)*A, h*x*B, h*x*Q*s, d*(1+x)*A]`.
On the punctured base disk `h` is nonzero, and the
representative is valid away from the degenerate frame
`d=0`. At `x=0, A=0`, set `A=x*b`; the lifted corner
maps to `[0:0:h*(1-s):h*Q*s:d*b]`. At the other end set
`X=1/x`, `A=X*b`. Its image is
`[2*h*b:2*h*b:h*(1-s):h*Q*s:0]`. The reciprocal ratio
charts cover `b=infinity`. The dependence on `b` exhibits
the auxiliary boundary data lost by taking separate
endpoint limits. This calculation resolves the two corners
of the diagonal cap in the zero-regulator parameter map.
The other cap pieces and comparison with the original
finite-regulator Lorentzian chain remain open.

## Regulated auxiliary infinity edge at z=1

On the real diagonal `y=x>0` at `z=1`, set
`v=x/(1+x)**2`, `G=1+(Q-1)*s`, and `H=1-(1+Q)*s`.
In the triangle variables of the preceding example,
the finite-regulator density is `-2*Q/D**3`, where
`D=H-4*i*eps*(A+v*(G-A)/4)**2/(v*A)`.
This follows from inserting the original Feynman regulator
`-i*eps*(1+a1+a2+a3+a4)**2` before projective normalization.

The argument of `A+v*(G-A)/4` lies in `[0,pi/4]`, since
it is `A*(1-v/4)+(v/4)*((1-s)+Q*s)`. The regulator shift
therefore lies in the sector `[-pi/2,0]`. Multiplication by
`exp(-i*pi/8)` gives `H` the constant imaginary part
`-sin(pi/8)` and gives the shift a nonpositive imaginary
part. This proves the same denominator gap for all
positive `eps`, all `x>0`, and the triangle interior.

In the joint corner `A=v*b`, for fixed positive `b`, the
regulated denominator tends to `H-4*i*eps*(b+G/4)**2/b`
at either `x` endpoint. At fixed interior `x`, the shift instead
grows like `eps/A` near the triangle edge. The zero-regulator
expansion of the integrated fiber begins with
`1+(i*v/2)*eps*log(1/eps)+O_x(eps)`.
For the coefficient, the edge multiplier is `v*G**2/(4*A)`.
Expanding the inverse cube gives `-6*i*Q/H**4` times
that multiplier, and the integral of `G**2/H**4` over `s`
is `-1/(3*Q)`. To bound the remainder, the strip `A<eps`
has area at most `eps` and the density is bounded by the gap.
Above that strip, the quadratic remainder is
`O_x(eps**2/A**2)` and integrates to `O_x(eps)`. The
`log(1-s)` from the upper edge is integrable. The example
checks the logarithmic coefficient with divided differences
in `log(1/eps)` to remove the order-`eps` constant term.

This estimate uses the sector at real `z=1`. The base-cap
continuation to complex `z` and its finite-regulator matching
to the original Lorentzian chain require additional analysis.

## Regulated poles on the diagonal complex base cap

On the diagonal `y=x`, continue the same finite-frame
projective coefficients and keep the auxiliary triangle fixed.
With `A,s,Q,G,H` as above, the normalized regulated denominator is
`D=H-i*eps*U**2/(4*x*A)`, where
`U=A*((5*z+3)/(1+z)+4*x)+x*(G-A)/(1+x)`.
In the cap chart `z=-1+R/k`, this is quadratic in `k`
because `(5*z+3)/(1+z)=5-2*k/R`. The two base punctures
are at `0,R`. The affine frame is singular at `5*R/2`.

Set `x=eps*chi` with `chi>0` fixed. Then `D` tends to
`H-i*A*(5-2*k/R)**2/(4*chi)`. The exact zeros are
`(R/2)*(5+4*x+x*(G-A)/(A*(1+x))) +/- (R/2)*sqrt(-4*i*chi*H/A)`.
For `R=.06`, `A=.5`, `s=.25`, and `chi=.25` there are two
simple zeros strictly inside the unit disk. They stay away
from the three punctures in the tested range and as
`eps` tends to zero. The example checks the zeros in the
original regulated polynomial before normalization, up to
the finite nonzero factor `A/C0`. It also counts index `+1`
at each zero, zero around each original puncture, and
total winding `+2` on the disk rim. The computation stays
in the punctured frame and uses the polynomial continuation
of the normalized denominator to count the winding.

These are intersections with the finite-regulator polar
divisor of the fixed-triangle product lift. The diagonal
pullback of the top form may vanish by rank, as in the
unregulated cap. The intersections still obstruct using
this lift as a chain in the divisor complement. The
calculation does not settle a lift with modified auxiliary
contours or extra meridian boundaries, or their matching
to the Lorentzian prescription.

## Regulated cap rim crossings and meridian balance

Keep `R=.06`, `A=.5`, and `s=.25` from the previous example.
Writing `chi=x/eps`, the limiting roots are `c +/- m*sqrt(chi)`
with `c=5*R/2` and `m=R*sqrt(-i*H/A)`. For either signed
direction `d=+m` or `d=-m`, the equation for `t=sqrt(chi)`
at the rim is `abs(d)**2*t**2+2*Re(conj(c)*d)*t+abs(c)**2-1=0`.
It has exactly one positive solution, and the radial derivative
is positive there. The two outward crossings are at
`chi=189.052122729...` and `chi=270.946552526...` in the joint
limit. The exact finite-regulator roots give nearby transverse
crossings that converge to these values.

Excise a disk around each regulated zero still in the cap.
In positive circle orientation each new meridian has index
one under the normalized denominator. Its contribution as
an inner boundary has index minus one. The original three
punctures have index zero on the tested slices. Consequently
the outer index drops from two to one to zero, and subtracting
the new meridian indices balances it on each slice. The code
checks this directly at `chi=.25`, `220`, and `320`.

A pole-free filling in a chart with a continuous nonvanishing
normalized denominator must have signed boundary index zero.
Holding the rim and the three original puncture loops fixed
therefore rules out removing a positive index by an interior
auxiliary deformation alone in that chart. The extra meridians
provide the slice-wise boundary correction. Their full sweep
over the auxiliary triangle and matching to the regulated
Lorentzian chain require further work.

## Regulated meridian sweep over the auxiliary triangle

Use `lambda=x/(eps*A)` to resolve the lower endpoint jointly
with the auxiliary infinity edge. For each `lambda>0`, the
roots are `c_eps +/- m(s)*sqrt(lambda)`, where
`c_eps=5*R/2+(R/2)*(4*x+eps*lambda*(G-A)/(1+x))` and
`m(s)=-i*R*sqrt(i*H(s))`. The right-half-plane square root
determines a continuous branch with `Im(m)<0`, including by
continuity at `s=0`. Indeed `i*H` lies in the closed right
half-plane and never vanishes for `0<=s<=1`. The limiting
roots at positive `lambda` avoid all three real punctures.
They can approach the frame puncture as `lambda` tends to zero.

In the limit use `c=5*R/2`. The rim equation has one
positive solution `lambda_+(s)` for `+m` and one solution
`lambda_-(s)` for `-m`, by the same quadratic as in the
fixed-slice check. Thus the number of regulated meridians
on a non-crossing slice is the sum of the two indicators
`lambda<lambda_+(s)` and `lambda<lambda_-(s)`. In the
unscaled coordinate, the crossing surfaces are
`x/eps=A*lambda_+(s)` and `x/eps=A*lambda_-(s)`.
The formula extends to `A=0` in the resolved scale, where
the unscaled thresholds tend to zero. This is a joint
limit of the root formula and does not identify the
original Feynman map at `A=x=0`.

The two thresholds coincide at `s=1/(1+Re(Q))`, since
`Re(m)=0` precisely when `Re(H)=0`. The two equal-scale
crossings are transverse in the `(s,lambda)` threshold
diagram. The roots are `c +/- i*sqrt(1-c**2)` there, two
distinct points of the rim. The example checks the formulas
across the limiting triangle, verifies finite-regulator
crossing convergence on selected slices, and measures the
signed boundary balance for each of the three possible counts.
The construction describes the lower-corner crossings and
excision circles. It leaves the upper endpoint, the limiting
frame-puncture face, and the Lorentzian boundary matching open.

## Regulated diagonal cap at the upper endpoint

For `R=.06` and `x>=14`, the normalized denominator of
the diagonal cap has a uniform gap of `sin(pi/8)` for
every positive regulator and every auxiliary triangle
point with `A>0`. To see this, write its regulator
numerator as `U=A*(L(x)-2*k/R)+x*G/(1+x)`, with
`L(x)=5+4*x-x/(1+x)` increasing. For `abs(k)<=1`, the
disk of possible values of `L-2*k/R` lies in the sector
`[-3*pi/16,5*pi/16]`: its distance from `L` to the nearer
sector boundary is `L*sin(3*pi/16)>2/R` at the cutoff.
The same sector contains `G=(1-s)+Q*s`. It is a convex
cone, so it contains `U`. Rotation by `exp(-i*pi/8)` then
puts the regulator shift `-i*eps*U**2/(4*x*A)` in the
closed lower half-plane. The rotated `H` has constant
imaginary part `-sin(pi/8)`. In particular the normalized
polynomial has zero winding on the rim, and both roots
are outside the disk in this range. The physical cap
still uses the original frame and base punctures.

In the upper chart `X=1/x`, `b=x*A`, the multiplier of
`-i*eps` tends to `(4*b+G)**2/(4*b)` as `X` tends to zero
at fixed `b>0`. If `eps` also tends to zero, the balances
`b=eps*beta` and `b=ell/eps` preserve the shifts
`-i*G**2/(4*beta)` and `-4*i*ell`, respectively. These
limits follow by substitution with `X` tending to zero.
Both are compatible with the uniform gap. The example
checks the projected gap, zero winding, both root
moduli, and the fixed and joint regulator limits.
The bound covers the upper region of this diagonal
product cap. Tracking regulated meridians in the finite
`x` region and near the frame puncture, and identifying
the remaining Lorentzian boundary data, remain open.

## Scope of the checks

The CD parameter calculation reproduces the established two-loop coefficient.
The chart checks establish local matching with a transported prescription at
finite regulator. The parameter tail is an estimate in the reduced integral.
The capped cycle realizes the benchmark scalar contour on the complex
weighted divisor. Identifying it with the Lorentzian contour class, including
the CD cycle and any additional chart-edge contributions, remains open.
