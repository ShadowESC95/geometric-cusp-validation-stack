# Gate 6: the Lorentzian chain on the weighted divisor

**Status.** This document gives a chain-level comparison between the
Lorentzian Feynman contour and the capped weighted-divisor cycle, by a route
that does not use the regulated cap of [GATE_6_PROOF.md](GATE_6_PROOF.md).
The local statement (sections 2 to 6) rests on exact algebraic identities and
is confirmed numerically to $`10^{-7}`$ or better. The passage from the local
statement to the double pole of the full integral (section 7) is treated in
[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md), which also corrects one statement of
section 7: the double pole must be taken with an analytic regulator, not a
hard cutoff.
The derivation was produced and machine-checked in one session and has not
been reviewed by a specialist. Section 8 lists exactly what is and is not
established.

Code: `src/geocusp/lorentzian_cusp.py`. Tests: `tests/test_gate6_closure.py`.
Example 61: `examples/lorentzian_cusp_contour.py`.

## 1. Result

Take the frozen-loop integral

```math
x_{13}^2x_{24}^2\int d^4x_0\,
\frac{F(z)}{x_{01}^2x_{02}^2x_{03}^2x_{04}^2},
\qquad
z=\frac{x_{02}^2x_{04}^2x_{13}^2}{x_{01}^2x_{03}^2x_{24}^2},
```

with every $`x_{0i}^2`$ carrying $`-i0`$, on a light-like square with
spacelike diagonals $`s=x_{13}^2<0`$, $`t=x_{24}^2<0`$ (metric $`+,-,-,-`$).
Write $`F(z)=f(\log z)`$.

**Edge statement.** Let $`x_0`$ approach the point $`x_2+\beta p_2`$ of the
null line through $`x_2`$ and $`x_3`$ ($`p_2=x_3-x_2`$). The integral over the
directions transverse to the line produces a collinear logarithm whose density
in $`\beta`$ is

```math
\frac{d\beta}{\beta(1-\beta)}\times
\begin{cases}
-2\pi i\,M[F], & 0<\beta<1,\\[2pt]
0, & \beta<0\ \text{or}\ \beta>1,
\end{cases}
\qquad
M[F]=\frac{1}{2\pi}\int_{-\pi}^{\pi}F(e^{i\phi})\,d\phi ,
```

up to a constant independent of $`F`$. The factor $`M[F]`$ arises as the
integral of $`f(L)\,dL`$ along the image of the real Lorentzian contour in
the $`L=\log z`$ plane. That image is an explicit path from $`L=i\pi`$ to
$`L=-i\pi`$; it is homotopic, with endpoints fixed, to the segment
$`L=i\phi`$, which is the unit circle.

**Consequences.**

1. The remaining $`\beta`$ integral diverges logarithmically at both ends of
   the segment, the two cusps. This is the origin of the double pole, whose
   coefficient is proportional to $`M[F]`$. That last statement is not
   proved from the edge statement: the edge density is not uniform as
   $`\beta\to0`$, and the proof of the double pole
   ([GATE_6_BOUNDS.md](GATE_6_BOUNDS.md) section 5) goes through the
   Mellin–Barnes form instead.
2. In the coordinates of the weighted blow-up (section 3) the Lorentzian
   chain and the capped cycle of
   [example 36](../examples/weighted_boundary_cycle.py) have the same class
   in the $`z`$ fibre, relative to the cap $`z=-1`$. The comparison chain is
   the region of the $`L`$ plane between the two paths (section 6).

## 2. Kinematics at the edge (exact)

Let $`p_1=x_2-x_1`$ and $`p_2=x_3-x_2`$ be the null edges at the cusp
$`x_2`$, with $`2p_1\!\cdot p_2=s`$. Decompose

```math
x_0=x_2+\alpha p_1+\beta p_2+y_\perp,\qquad \rho=|y_\perp|^2 .
```

Then, exactly,

```math
x_{02}^2=w,\qquad x_{01}^2=w+q_1,\qquad x_{03}^2=w+q_3,
\qquad
w=\frac{q_1q_3}{|s|}-\rho,\quad q_1=\beta s,\quad q_3=-\alpha s ,
```

and on the plane $`y_\perp=0`$, exactly,
$`x_{04}^2=t(1+\alpha-\beta)+\alpha\beta s`$. On the edge itself ($`\alpha\to0`$) the
fourth propagator is $`t(1-\beta)`$: spacelike inside the segment and
timelike beyond the far cusp.

The transverse integral at fixed $`(q_1,q_3)`$ is, up to the constant
$`1/x_{04}^2`$,

```math
K_f(q_1,q_3)=\int_0^\infty d\rho\,
\frac{f(L(w))}{w\,(w+q_1)\,(w+q_3)},\qquad w\to w-i0 .
```

The same $`-i0`$ enters all three factors and $`L`$, because all three
depend on $`y^2`$ alone. The integrand is therefore the boundary value of one
function holomorphic in the lower half $`w`$ plane.

## 3. Dictionary with the weighted divisor (exact)

In the chart $`A=e_1+ae_3+be_4`$, $`B=e_2+ce_3+de_4`$ of
`strict_cusp_weighted_blowup`,

```math
\langle AB12\rangle=ad-bc,\quad
\langle AB23\rangle=d,\quad
\langle AB14\rangle=-a,\quad
\langle AB34\rangle=1,\qquad
1+z=\frac{bc}{ad}.
```

| Weighted chart | Position space | Weight |
| --- | --- | --- |
| $`a`$ | $`x_{01}^2\propto\beta`$, position along the edge | 0 |
| $`d`$ | $`x_{03}^2\propto\alpha`$ | 2 |
| $`bc`$ | $`\rho=\vert y_\perp\vert ^2`$ | 1 + 1 |
| $`ad-bc`$ | $`x_{02}^2=w`$ | 2 |
| phase of $`b`$ | angle $`\theta`$ of $`y_\perp`$ | 0 |
| $`1+z=bc/(ad)`$ | $`1+z=\rho\vert s\vert /(q_1q_3)`$ at the cusp | 0 |

The strict cusp line $`b=c=d=0`$ is the null line through $`x_2`$ and
$`x_3`$; $`a=0`$ is the cusp $`x_2`$. The weights $`(1,1,2)`$ of the blow-up
are the collinear scaling $`y_\perp\sim\lambda`$, $`\alpha\sim\lambda^2`$.
The exceptional divisor $`E`$ has coordinates $`(a,z,\theta)`$, and the
radial variable $`\lambda`$ carries the collinear logarithm. The tests check
the bracket identities, the equality of the twistor and position-space cross
ratios, and that the strict line is the edge.

On $`E`$ the cap $`bc=0`$ is exactly $`z=-1`$.

## 4. The measure (exact)

With $`z_{\rm e}=-(|s|/q_1)\,w/(w+q_3)`$,

```math
\frac{dw}{w(w+q_1)(w+q_3)}
=\frac{1}{q_1q_3}\,d\log z_{\rm e}
+\frac{1}{q_1(q_1-q_3)}\Big[\frac{dw}{w+q_1}-\frac{dw}{w+q_3}\Big],
\qquad
z=\frac{z_{\rm e}}{1+w/q_1}.
```

The first term carries the collinear pole $`1/q_3`$. The second has no
$`1/q_3`$ and gives no collinear logarithm. On the range that matters
($`|w|\sim|q_3|\ll|q_1|`$) the factor $`1+w/q_1`$ tends to 1, so
$`d\log z_{\rm e}\to dL`$. Hence

```math
q_1q_3\,K_f(q_1,q_3)\;\longrightarrow\;-\int_{\gamma}f(L)\,dL
\qquad(q_3/q_1\to0),
```

where $`\gamma`$ is the image of the $`w`$ contour in the $`L`$ plane, with
corrections of order $`r\log(1/r)`$, $`r=|q_3/q_1|`$. The $`q_3`$ integral
$`\int dq_3/q_3`$ is a principal value across $`q_3=0`$, so the coefficient
of the collinear logarithm is the jump

```math
\big[q_1q_3K_f\big]_{q_3>0}-\big[q_1q_3K_f\big]_{q_3<0}
=-\int_{\gamma_{+}-\gamma_{-}}f(L)\,dL .
```

## 5. Which sheet (exact bookkeeping)

Each $`x_{0i}^2-i0`$ is continued from the spacelike region through the lower
half plane. Its argument is $`-\pi`$ when spacelike and $`0`$ when timelike.
Since $`s/t>0`$,

```math
\arg z=\pi\,(\tau_2+\tau_4-\tau_1-\tau_3),\qquad
\tau_i=\begin{cases}1,&x_{0i}^2>0\\0,&x_{0i}^2<0.\end{cases}
```

The argument can reach $`\pm2\pi`$: the Lorentzian contour leaves the
principal sheet. Treating $`F`$ as a function of $`z`$ on the principal
sheet with a choice of lip gives wrong, cutoff-dependent answers; this
bookkeeping is what makes the cutoff dependence cancel.

## 6. The four paths and the comparison chain

Near the cusp ($`\tau_4=0`$), as $`\rho`$ runs from 0 to its collinear
cutoff, $`|z|`$ starts at 1 (the cap) and ends at $`\Lambda=|s|/|q_1|`$. The
sequence of values of $`\arg z/\pi`$ is:

| Row | $`x_{01}^2`$ | $`x_{03}^2`$ | $`\arg z/\pi`$ along the path | Path in the $`L`$ plane |
| --- | --- | --- | --- | --- |
| 1 | $`-`$ | $`-`$ | $`1,\;0`$ | from $`i\pi`$ to small $`\vert z\vert `$ on $`\arg=\pi`$; then out to $`\log\Lambda`$ on $`\arg=0`$ |
| 2 | $`-`$ | $`+`$ | $`-1,\;0`$ | from $`-i\pi`$ to large $`\vert z\vert `$ on $`\arg=-\pi`$; then in to $`\log\Lambda`$ on $`\arg=0`$ |
| 3 | $`+`$ | $`+`$ | $`-1,\;-2,\;-1`$ | from $`-i\pi`$ to small $`\vert z\vert `$ on $`\arg=-\pi`$; small to large on $`\arg=-2\pi`$; in to $`\log\Lambda-i\pi`$ |
| 4 | $`+`$ | $`-`$ | $`-1`$ | from $`-i\pi`$ straight to $`\log\Lambda-i\pi`$ |

**Spacelike $`x_{01}^2`$ (inside the edge segment).** Let $`\gamma_\pm`$ be
the paths for $`x_{03}^2\gtrless0`$ (rows 2 and 1). Their difference
$`\gamma_+-\gamma_-`$ is one path from $`L=-i\pi`$ to $`L=+i\pi`$: out along
$`\arg=-\pi`$ past $`z=\infty`$, back along $`\arg=0`$ from large to small
$`|z|`$, past $`z=0`$, and in along $`\arg=\pi`$. The cutoff $`\Lambda`$ has
dropped out. Here "past $`z=0`$" and "past $`z=\infty`$" stand for the
arcs traced when the $`w`$ contour, displaced a finite distance into the
lower half plane, goes below the poles $`w=0`$ and $`w=-q_3`$; the path is
compact and lies in the strip $`|\operatorname{Im}L|\le\pi`$, where $`f`$ is
holomorphic. Define the two-chain

```math
\Sigma=\{\operatorname{Re}L>0,\ -\pi<\operatorname{Im}L<0\}
-\{\operatorname{Re}L<0,\ 0<\operatorname{Im}L<\pi\},
```

truncated by those arcs, each piece with its complex orientation. Then

```math
\partial\Sigma=(\gamma_+-\gamma_-)-\{L=i\phi:\ \phi\ \text{from}\ -\pi\ \text{to}\ \pi\},
```

and by Cauchy's theorem $`\int_{\partial\Sigma}f\,dL=0`$. Therefore

```math
\int_{\gamma_+-\gamma_-}f\,dL=i\int_{-\pi}^{\pi}f(i\phi)\,d\phi=2\pi i\,M[F],
\qquad
\big[q_1q_3K_f\big]_{q_3>0}-\big[q_1q_3K_f\big]_{q_3<0}=-2\pi i\,M[F].
```

The two endpoints of the path are the two lips of the cap $`z=-1`$.
$`\Sigma\times S^1_\theta`$ is the comparison chain on the fibre of $`E`$:
its boundary is the Lorentzian chain minus the unit-circle chain, relative
to the cap, and $`f\,dL`$ integrates to zero over that boundary.

**Timelike $`x_{01}^2`$ (the same null line outside the segment).** The
difference of rows 3 and 4 is a closed loop running around the strip
$`-2\pi\le\operatorname{Im}L\le-\pi`$. It bounds, $`f`$ is holomorphic
inside, and the integral is zero for every $`f`$.

**At a finite point of the edge.** The same analysis holds for any
$`0<\beta<1`$, not only near the cusp, once $`x_{04}^2`$ is given its edge
value $`t(1-\beta)`$: the path then starts exactly at $`z=-1`$ and the
contour is exactly the unit circle. Beyond the far cusp ($`\beta>1`$)
$`\tau_4`$ flips and the loop again closes up. Holding $`x_{04}^2`$ at $`t`$
instead moves the endpoints to $`|z|=1/(1-\beta)`$; this is a different
representative of the same class.

**What the comparison shows.** In the limit the difference
$`\gamma_+-\gamma_-`$ is a path from $`L=-i\pi`$ to $`L=i\pi`$ inside the
domain where $`f`$ is holomorphic. Any such path gives the same integral,
$`g(i\pi)-g(-i\pi)`$ with $`g'=f`$. The unit circle is one representative;
the actual image runs out towards $`z=0`$ and $`z=\infty`$. So the statement
is about the endpoints, the two lips of $`z=-1`$, and the sheet. It concerns
the $`z`$ fibre only.

**Relation to the capped cycle.** The cycle of example 36 is
$`\{|a|=R\}\times\{z=e^{i\phi}\}\times\{\theta\}`$. The Lorentzian chain on
$`E`$ is $`\{a\ \text{real}\}\times\gamma\times S^1_\theta`$. The factor in
$`z`$ is identified above. The circle in $`a`$ is the residue at the cusp
$`a=0`$, which is the coefficient of the second logarithm from
$`\int da/a`$ (section 7); this identification of the $`a`$ factor is an
interpretation, not a proved equality of chains. The $`\theta`$ factor is the compact circle in
both. The overall constant, including the order-two quotient in $`\theta`$,
is fixed by $`F=1`$.

## 7. From the edge statement to the double pole

The edge statement gives the coefficient of the collinear singularity as
$`-2\pi i\,M[F]\,d\beta/(\beta(1-\beta))`$ on the segment, for any regulator
acting on the collinear invariant. The $`\beta`$ integral has a logarithmic
endpoint singularity at each cusp, which supplies the second pole. Near the
cusp $`x_2`$ the sector $`|q_1|\ll|q_3|`$ is the edge statement for the
adjacent edge, and the sum over the four sign quadrants tends to
$`2\pi i\,M[F]`$ in both sectors.

**Correction.** An earlier version of this section inferred a
$`\log^2`$ coefficient with a hard cutoff on the collinear invariants. That
is right only for weights bounded at $`z\to0,\infty`$. Where the two
invariants are comparable, growing weights produce cutoff-dependent terms
larger than $`\log^2`$; they are artefacts of cutting through the light
cones and are absent in an analytic regulator. See
[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md), section 4. The double pole is
therefore defined in dimensional regularisation, where
$`L[F]/L[1]=M[F]`$ follows from the Mellin–Barnes evaluation of
[GATE_6_FUNCTIONAL.md](GATE_6_FUNCTIONAL.md) and Lemma 3 of GATE_6_BOUNDS.
What the present document adds is the reason: the circle is the image of the
real integration region at a light-like edge.

All four cusps, both edges at each, and all three real kinematic regions
are computed in GATE_6_BOUNDS section 2.

## 8. What is established and what is not

**Established by exact algebra (sympy) and explicit vectors.**

- The decomposition of section 2 and the edge value of $`x_{04}^2`$.
- The measure identity and the factorisation $`z=z_{\rm e}/(1+w/q_1)`$.
- The bracket dictionary of section 3 and the equality of the twistor and
  position-space cross ratios.

**Established numerically, to the stated accuracy.**

- The four sequences of $`\arg z`$ in section 6.
- The jump equals $`-2\pi i\,M[F]`$ at $`\beta=0.1,0.5,0.9`$ and vanishes at
  $`\beta=-0.3,1.4`$, to $`10^{-7}`$, for $`F=z^{0.3},z^{-0.6}, \log^2z+\pi^2,\log^3z`$.
- The quadrant sum equals $`2\pi i\,M[F]`$ in both collinear sectors, to
  $`10^{-7}`$; the error is of order $`r\log(1/r)`$ in the ratio $`r`$ of
  scales.
- With $`x_{04}^2`$ held at $`t`$ the contour is the circle
  $`|z|=1/(1-\beta)`$.

**Treated in [GATE_6_BOUNDS.md](GATE_6_BOUNDS.md).**

- The rate of the collinear limit (Lemma 1) and the effect of the exact
  fourth propagator (Lemma 2), both proved.
- All four cusps and all three kinematic regions (computed).
- $`\bar F^{(2)}`$ on its continued sheets: regular at $`z=-1`$ (proved
  exactly), and the edge chain gives $`-11\pi^4/15`$.
- The relation to dimensional regularisation, with a uniform bound on the
  regulated integral, and the failure of hard cutoffs for growing weights.
- The pinch surfaces: only the four edge segments.

**Constructed in [GATE_6_CONTOUR.md](GATE_6_CONTOUR.md).**

- An explicit global integration contour for every light-like square, with
  its properties proved by explicit constants, absolute integrability proved
  line by line, and the parametric representation of the regulated integral
  derived, including the contour-independence and limiting step.

**Not addressed.**

- $`F`$ is still an input computed from the Lorentzian integral. The
  geometry supplies the cycle, not the weight.
- The cycle is selected by the Lorentzian real structure and the Feynman
  $`-i0`$. These are physical inputs. They are not derived from positivity
  of the amplituhedron.
- The regulated complex cap of GATE_6_PROOF.md sections 2 to 5, with its rim
  collisions and detours, is bypassed, not resolved. This route never
  introduces that lift.
- The edge statement is not a step in the proof of the double pole. It
  explains the circle; the proof is the Mellin–Barnes argument of
  GATE_6_BOUNDS section 5, which re-derives the rule of arXiv:2112.06956
  and extends it from powers to a class of weights.
- The three-loop value from the edge chain is, by Cauchy's theorem, the
  same number as the circle mean. It tests the continuation of
  $`\bar F^{(2)}`$ to its second sheet, not an independent quantity.
- Lemmas 1 and 2 are proved for weights of polynomial growth in $`L`$. The
  numerical checks with $`z^{0.3}`$ show the same rate but are not covered
  by those proofs.

## 9. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_closure.py
PYTHONPATH=src python -m examples.lorentzian_cusp_contour
```
