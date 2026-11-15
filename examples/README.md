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

## Scope of the checks

The CD parameter calculation reproduces the established two-loop coefficient.
The chart checks establish local matching with a transported prescription at
finite regulator. The parameter tail is an estimate in the reduced integral.
Constructing the global cycle on the weighted twistor divisor and controlling
its remaining chart-edge contributions are open steps in the geometric-period
proposal.
