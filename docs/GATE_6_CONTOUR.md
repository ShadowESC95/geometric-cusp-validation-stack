# Gate 6: an explicit integration contour

[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md) left one step open: a global integration
contour on which the regulated integrand is absolutely integrable. It was
needed twice, for the exchange of the Mellin superposition with the loop
integral and for the statement that only the edges produce poles. This
document constructs such a contour, and replaces the appeal to Stirling's
formula in Lemma 4 by elementary bounds.

**Status.** The construction is explicit and implemented, in all three real
kinematic regions and for every light-like square: the constants are derived
from the geometry. The displacement is continuously differentiable. Its
central inequality is proved (Lemma 5). Its global properties, the strips
and the bound at infinity are reduced in section 9 to inequalities between
named constants, each proved, and a certificate evaluates all of them; they
hold across the family of squares. Absolute integrability is proved in
section 10. The regulated integral on the contour is identified with the
Mellin–Barnes form in sections 11 and 12, with proofs, including a direct
proof of the one Fourier inversion that is needed. Nothing here has been
reviewed by a specialist.

Code: `src/geocusp/feynman_contour.py`, `src/geocusp/schwinger.py`. Tests:
`tests/test_gate6_contour.py`, `tests/test_schwinger.py`.
Example 63: `examples/feynman_contour_survey.py`.

## 1. The form of the displacement

Displace the real domain by a purely imaginary vector field:

```math
x_0\ \longmapsto\ x_0+iV(x_0),\qquad V\ \text{real}.
```

Then, exactly,

```math
x_{0i}^2\ \longmapsto\ x_{0i}^2-V^2+2i\,V\!\cdot(x_0-x_i).
```

The imaginary part is linear in $`V`$. Call $`V`$ *admissible for $`i`$ at
$`x_0`$* if $`V\!\cdot(x_0-x_i)\le0`$. Admissible fields form a convex cone,
so fields that are admissible on overlapping regions can be joined by a
partition of unity. This is what makes a global construction possible.

The Feynman prescription is the limit of small admissible $`V`$ near each
cone: the invariant moves into the lower half plane.

## 2. Frame fields

At the cusp $`x_j`$ write $`x_0=x_j+\alpha p_{j-1}+\beta p_j+y_\perp`$,
$`D=2p_{j-1}\!\cdot p_j`$, $`\rho=|y_\perp|^2`$, $`P=(\beta,\alpha)`$. Take
$`D<0`$ first (a spacelike diagonal at this cusp). Let

```math
T=\{\beta\ge0,\ \alpha\le0,\ \beta-\alpha\le1\}
```

be the triangle with vertices $`(0,0)`$, $`(0,-1)`$, $`(1,0)`$, which are
the images of $`x_j`$, $`x_{j-1}`$, $`x_{j+1}`$. Write
$`T=\{\ell_1\le0,\ell_2\le0,\ell_3\le0\}`$ with $`\ell_1=-\beta`$,
$`\ell_2=\alpha`$, $`\ell_3=\beta-\alpha-1`$, and let

```math
\Psi(P)=\sum_{k=1}^3\phi(\ell_k(P)),\qquad
\phi(t)=\begin{cases}0,&t\le0\\ t^3/3,&0\le t\le1\\ t^2-t+\tfrac13,&t\ge1.\end{cases}
```

$`\phi`$ is convex and twice continuously differentiable with
$`0\le\phi''\le2`$. So $`\Psi`$ is convex and $`C^2`$, it vanishes exactly on
$`T`$, and its Hessian is at most $`2\sum_k|\nabla\ell_k|^2=8`$. Let
$`d(P)=\nabla\Psi(P)`$ and define

```math
V_j=\kappa\,(d_\beta\,p_{j-1}+d_\alpha\,p_j)+\nu\,y_\perp,\qquad\kappa,\nu>0 .
```

**Lemma 5.** For $`i\in\{j-1,j,j+1\}`$ and every real $`x_0`$,

```math
V_j\!\cdot(x_0-x_i)\ \le\ -\frac{\kappa|D|}{2}\,\Psi(P)-\nu\rho\ \le\ 0,
\qquad\text{and}\qquad|d(P)|^2\le16\,\Psi(P).
```

*Proof.* A direct computation with $`p_{j-1}^2=p_j^2=0`$ and
$`y_\perp\!\cdot p_{j-1}=y_\perp\!\cdot p_j=0`$ gives

```math
V_j\!\cdot(x_0-x_i)=\frac{\kappa D}{2}\;d\cdot(P-v_i)-\nu\rho ,
```

where $`v_i`$ is the vertex of $`T`$ belonging to $`x_i`$. By convexity
$`\Psi(v)\ge\Psi(P)+d\cdot(v-P)`$, and $`\Psi(v)=0`$ for $`v\in T`$, so
$`d\cdot(P-v)\ge\Psi(P)`$; and $`D<0`$. For the second claim, a non-negative
function whose Hessian is at most $`L`$ satisfies
$`|\nabla\Psi|^2\le2L\Psi`$: otherwise the Taylor bound
$`\Psi(P-\nabla\Psi/L)\le\Psi(P)-|\nabla\Psi|^2/(2L)`$ would make $`\Psi`$
negative. $`\blacksquare`$

An earlier version used $`d=P-\pi_T(P)`$, the vector from the nearest point
of $`T`$. That field is Lipschitz but not differentiable. The present one is
continuously differentiable, and so is the whole displacement $`V`$: the
cutoffs are $`C^1`$ and are locally constant where their arguments are not
smooth. Only the ordinary Stokes theorem for a $`C^1`$ map is needed below.

So $`V_j`$ is admissible for its three invariants everywhere, and strictly
so unless $`\rho=0`$ and $`P\in T`$, that is unless $`x_0`$ lies in the flat
triangle $`\operatorname{conv}(x_{j-1},x_j,x_{j+1})`$. On that triangle the
three invariants are $`\alpha\beta D`$, $`\beta D(1+\alpha)`$ and
$`\alpha D(\beta-1)`$, which vanish only on the two edges at $`x_j`$. The
triangle is exactly the set on which the Landau conditions fail to give a
displacement direction for all three invariants at once.

Near the edge through $`x_j`$ and $`x_{j+1}`$ the longitudinal part
vanishes to second order and the field reduces to the rotation of the
transverse plane used in GATE_6_BOUNDS, with $`\rho\to(1+i\nu)^2\rho`$.

The frame field does not control the fourth invariant $`x_{0,j+2}^2`$.

**Timelike diagonal.** If $`D>0`$ at the cusp, replace $`\kappa`$ by
$`-\kappa`$ in $`V_j`$. The computation in the proof gives
$`(-\kappa D/2)\,d\cdot(P-v_i)-\nu\rho`$, and the bound of Lemma 5 holds with
$`|D|`$. The transverse plane is spacelike for either sign of $`D`$, since
two independent null vectors span a Lorentzian plane. Nothing else changes:
the far fields need only points in the common past and future, which always
exist. In the regions $`(s,t)=(+,-)`$ and $`(-,+)`$ two of the four cusps
have $`D>0`$.

## 3. Far fields

Let $`c_-`$ lie in the common chronological past of the four cusps and
$`c_+`$ in their common future. For $`V=-\kappa'(x_0-c_-)`$,

```math
V\!\cdot(x_0-x_i)=-\kappa'\big[x_{0i}^2+(x_i-c_-)\!\cdot(x_0-x_i)\big].
```

On the future nappe of the cone of $`x_i`$ the first term vanishes and the
second is the product of a future timelike and a future null vector, which
is positive. So this field is admissible near every future nappe at once,
and $`-\kappa'(x_0-c_+)`$ near every past nappe.

A point can be near a future nappe of one cone and a past nappe of another
only in a bounded region: $`x_i-x_j`$ would have to be nearly a difference
of a past and a future null vector, hence nearly causal, and for adjacent
cusps this forces $`x_0`$ to the edge between them.

## 4. The global field

Fix a ball of radius $`R`$ about the centre of the square and a width
$`\delta_0`$.

- Inside the ball: $`V_{\rm in}=\sum_jm_jV_j`$ with
  $`m_j=\chi(|x_{0,j+2}^2|/\delta_0)`$, where $`\chi`$ vanishes on $`[0,1]`$
  and equals 1 on $`[2,\infty)`$. The frame field at $`x_j`$ is switched off
  near the cone of the opposite cusp, the one invariant it does not control.
- Outside $`2R`$: $`V_{\rm out}=-\kappa'\Theta_+(x_0-c_-)-\kappa'\Theta_-(x_0-c_+)`$,
  with $`\Theta_\pm`$ supported where the time coordinate is at least a
  quarter of the spatial distance, in the future or the past.
- In between the two are blended.

**Why something is always switched on.** Near the cone of $`x_i`$ the usable
frames are $`j=i-1,i,i+1`$, with weights
$`\chi(|x_{0,i+1}^2|/\delta_0)`$, $`\chi(|x_{0,i+2}^2|/\delta_0)`$,
$`\chi(|x_{0,i-1}^2|/\delta_0)`$. All three vanish only if $`x_0`$ is within
$`2\delta_0`$ of all four cones. That does not happen:

```math
\min_{x_0\in\mathbb R^4}\ \max_i|x_{0i}^2|\approx0.199
```

for the square used here (numerical minimisation), and $`2\delta_0=0.05`$.
Section 9 proves a weaker bound, 0.0897, which is still enough.
This is the statement that the four cones have no common real point.

**Constants.** $`R`$ is chosen so that the far fields are admissible outside
the ball. $`\delta_0`$ is then chosen so that $`2\delta_0`$ is below the
minimum above. Finally $`\kappa,\nu\ll\delta_0/R^2`$, so that where a frame
field is switched on but does not control an invariant, that invariant is at
least $`\delta_0`$ in modulus and moves by a small fraction of itself. The
values used are $`R=8`$, $`\tau_0=8`$ (the time offset of $`c_\pm`$),
$`\delta_0=0.025`$, $`\kappa=2\times10^{-8}`$, $`\nu=1.5\times10^{-7}`$,
$`\kappa'=2\times10^{-3}`$. Section 9 states the precise conditions.

## 5. Properties

**Proposition 6.** For the field of section 4:

1. near the cone of $`x_i`$, $`\operatorname{Im}x_{0i}^2\le0`$ on the
   contour, and the invariant itself vanishes only on the edges (the
   imaginary part alone also vanishes on the flat triangles);
2. the argument of each invariant stays within a small angle $`e`$ of
   $`[-\pi,0]`$, so on the contour $`|\arg z|\le2\pi+4e`$ for spacelike
   diagonals and $`|\arg z|\le3\pi-0.3`$ when one diagonal is timelike
   (Lemma Q);
3. the invariants vanish on the contour only on the four edge segments;
4. near a cusp, in its frame,
   $`|x_{0,j}^2|\ge c\,\nu\max(\rho,|\alpha\beta D|)`$,
   $`|x_{0,j-1}^2|\ge c\,\nu\max(\rho,|\beta D||1+\alpha|)`$,
   $`|x_{0,j+1}^2|\ge c\,\nu\max(\rho,|\alpha D||1-\beta|)`$;
5. the map $`x_0\mapsto x_0+iV`$ is a local diffeomorphism onto the contour.

*Argument.* (1) and (3): by Lemma 5 and section 3 every field that is
switched on near a cone is admissible for it; by section 4 at least one
frame field has weight 1 there, and its bound vanishes only on its triangle,
which meets the cone only on the edges. (2): by (D) of section 9 each
argument lies in $`[-\pi-e,e]`$; $`\arg z`$ is the sum of two of them minus
the other two, plus $`\pi(\tau_s-\tau_t)`$ from the diagonals.
(4): if $`\rho\ge|\alpha\beta D|/2`$ the imaginary part is at least
$`2\nu\rho`$; otherwise the real part is at least $`|\alpha\beta D|/2`$; the
other two are the same. (5): the Jacobian is $`1+i\,\partial V`$ with
$`\partial V`$ real, whose determinant vanishes only if $`\partial V`$ has
an eigenvalue of modulus 1; $`\|\partial V\|<1`$.

**Lemma P (phases at real points).** At every real point where no
invariant vanishes, with $`\tau_i=1`$ if $`x_{0i}^2>0`$ and 0 otherwise,

```math
\frac{\arg z}{\pi}=\tau_2+\tau_4-\tau_1-\tau_3+\tau_s-\tau_t\in
\begin{cases}[-2,2],&s<0,\ t<0\\ [-1,2],&s>0\\ [-2,1],&t>0.\end{cases}
```

*Proof.* The formula is the Feynman prescription. For $`s,t<0`$ the range
is immediate. Let $`s>0`$, so $`x_3-x_1`$ is timelike; say it points to the
future. The null vectors $`x_2-x_1`$ and $`x_3-x_2`$ have product
$`s/2>0`$, so they have the same time orientation, and their sum points to
the future, so both do. Suppose $`\tau_2=1`$. If $`x_0-x_2`$ points to the
future then $`x_0-x_1=(x_0-x_2)+(x_2-x_1)`$ is a future timelike vector plus
a future null vector, hence timelike: $`\tau_1=1`$. If it points to the
past then $`x_3-x_0=(x_3-x_2)+(x_2-x_0)`$ is timelike: $`\tau_3=1`$. The
same holds with $`x_4`$ in place of $`x_2`$. So
$`\tau_2+\tau_4-\tau_1-\tau_3\le1`$ and $`\arg z\le2\pi`$; the lower bound
$`-\pi`$ is immediate. The case $`t>0`$ is the mirror image. $`\blacksquare`$

On the contour the arguments are not close to their real values near two
cones at once, and the bound $`2\pi`$ is not proved there. It is observed:
`phase_survey` finds $`\arg z/\pi`$ in $`[-2,2]`$, $`[-1,2]`$ and
$`[-2,1]`$ on the contour in the three regions, on 30 000 stress points
each. A search aimed at the cusps gives the same: in the frame of each
cusp, with the two null coordinates and the transverse distance spread
over ten decades and both signs, 6000 points and 25 local maximisations
on each of six squares with a timelike diagonal find
$`\max|\arg z|/\pi`$ between $`1.9999994`$ and $`2.0000000`$. Along the
deformation family with $`\varepsilon`$ up to $`10^{-3}`$ a review found
$`2.000003\pi`$, which may be due to $`\varepsilon`$. Since each argument
can leave its real value by $`e`$, the statement to aim for is
$`2\pi+4e`$. A proof would need, in the quadrant beyond a cusp where the
invariants of that cusp and of both neighbours are negative, an upper
bound $`|\operatorname{Im}x_{0j}^2|\le c\,|\operatorname{Re}x_{0j}^2|`$ in
every zone, and the conditions of section 9 give only lower bounds on the
imaginary parts. I have not found it. What is needed in GATE_6_BOUNDS is
less: a bound strictly below $`3\pi`$. That is proved where the
certificate verifies (Q1) and (Q3).

**Lemma Q (the corner is excluded).** Let one diagonal be timelike and let
conditions (Q1) and (Q3) below hold. Then $`|\arg z|\le3\pi-c_0`$ on the
contour, with $`c_0=0.3`$.

*Proof.* Take $`s>0`$; the case $`t>0`$ is the mirror image. Write
$`a_i\in[-\pi-e,e]`$ for the arguments on the contour, so
$`\arg z=\pi+a_2+a_4-a_1-a_3\ge-\pi-4e`$. Suppose $`\arg z>3\pi-c_0`$. Then
$`a_2>-\gamma`$ and $`a_1,a_3<-\pi+\gamma`$ with $`\gamma=c_0+3e<\pi/2`$, so
on the contour

```math
\operatorname{Re}x_{02}^2>0,\qquad\operatorname{Re}x_{01}^2<0,\qquad\operatorname{Re}x_{03}^2<0,\qquad
|\operatorname{Im}x_{02}^2|\le\tan\gamma\;\operatorname{Re}x_{02}^2 .
```

The real parts are the real invariants minus the same number $`V^2`$. In
the frame of $`x_2`$, where $`D=s>0`$, the real invariants are $`w`$,
$`w+A`$ and $`w+B`$ with $`A=\beta D`$, $`B=-\alpha D`$ and
$`w=-AB/D-\rho`$. The three inequalities give $`A<V^2-w<0`$ and $`B<0`$, so
$`AB>0`$ and $`w<0`$. Hence

```math
V^2<w<0,\qquad|x_{02}^2|<|V^2|,\qquad|\operatorname{Im}x_{02}^2|\le\tan\gamma\,|V^2| .
```

*Claim.* At every real point with $`V^2<0`$ and $`|x_{0i}^2|\le|V^2|`$,

```math
|V^2|\le\mu_Q\,|\operatorname{Im}x_{0i}^2|,\qquad\mu_Q\tan\gamma<1 .
```

Given the claim, $`|V^2|\le\mu_Q\tan\gamma\,|V^2|<|V^2|`$, a contradiction.

*Proof of the claim for $`r\le2R`$.* By (D), $`|V^2|\le v^2<\delta_0/2`$, so
$`|x_{0i}^2|\le\delta_0/2`$ and the frame $`i+2`$ is switched off. Each
remaining frame has $`I_j:=-2V_j\!\cdot n\ge\kappa|D_j|\Psi_j+2\nu\rho_j`$ by
Lemma 5. In Euclidean norm
$`|V_j|^2\le4\kappa^2P^2|d|^2+2\nu^2B^2\rho_j\le64\kappa^2P^2\Psi_j+2\nu^2B^2\rho_j\le\mu I_j`$
with

```math
\mu=\max\Big(\frac{64\kappa P^2}{\min_j|D_j|},\ \nu B^2\Big),
```

where $`|y_\perp|\le B\sqrt\rho`$. At most three frames are on and
$`m_j\le1`$, so $`|V_{\rm in}|^2\le3\mu\,I_{\rm in}`$ with
$`I_{\rm in}=-2V_{\rm in}\!\cdot n=\sum m_jI_j`$. For $`r\le R`$ this is the
claim with $`\mu_Q\ge3\mu`$.

For $`R<r<2R`$, $`V=(1-g)V_{\rm in}+gW`$ with $`W=-\kappa'(x_0-c_\mp)`$,
because $`\Theta=1`$ near a cone by (C2). By (C1),
$`I_W:=-2W\!\cdot n=2\kappa'(n^2+u\!\cdot n)\ge2\kappa'M`$, with $`M`$ the
margin of (C1). Since $`(x_0-c_\mp)^2=n^2+2n\!\cdot u+u^2`$ and
$`2|n\!\cdot u|\le I_W/\kappa'+\delta_0`$,

```math
|W^2|\le\kappa'^2\big(\tfrac32\delta_0+(\tau_0+\Delta)^2\big)+\kappa'I_W\le\mu_WI_W,\qquad
\mu_W=\kappa'\Big(1+\frac{\tfrac32\delta_0+(\tau_0+\Delta)^2}{2M}\Big).
```

For the cross term, $`|W|\le K=\kappa'(2R+\tau_0+\Delta)`$ and
$`2|V_{\rm in}||W|\le\sqrt{3\mu}\,I_{\rm in}+\sqrt{3\mu}\,K^2`$. Adding the
three terms of $`V^2`$,

```math
|V^2|\le\mu_Q\big[(1-g)I_{\rm in}+gI_W\big]=\mu_Q|\operatorname{Im}x_{0i}^2|,\qquad
\mu_Q=\max\Big(3\mu+\sqrt{3\mu},\ \mu_W+\frac{\sqrt{3\mu}\,K^2}{2\kappa'M}\Big).
```

*Proof of the claim for $`r\ge2R`$.* Here $`V=bW`$ with $`0\le b\le1`$, and
$`b>0`$ requires $`\tau>\tfrac14|\xi|`$, so $`n_0\ge2R/\sqrt{17}-\Delta>0`$.
We show that $`V^2<0`$ and $`|n^2|\le|V^2|`$ cannot both hold. Put
$`q=\sqrt{|n^2|}`$, so $`|\vec n|\le n_0+q`$ and
$`n\!\cdot u\ge n_0(u_0-|\vec u|)-q|\vec u|\ge n_0(\tau_0-2\Delta)-q\Delta`$.
If $`n\!\cdot u\ge0`$ then $`|(x_0-c)^2|=-n^2-2n\!\cdot u-u^2\le|n^2|`$, so
$`|n^2|\le\kappa'^2|n^2|`$, $`n^2=0`$ and $`(x_0-c)^2>0`$. If
$`n\!\cdot u<0`$ then $`|n\!\cdot u|\le q\Delta`$ and
$`q^2\le\kappa'^2(q^2+2q\Delta)`$, so $`q\le2\kappa'^2\Delta/(1-\kappa'^2)`$,
while $`n_0(\tau_0-2\Delta)\le q\Delta`$. This contradicts (Q3).
$`\blacksquare`$

*Conditions.* (Q1): $`c_0+3e<\pi/2`$ and $`\mu_Q\tan(c_0+3e)<1`$. (Q3):
$`(2R/\sqrt{17}-\Delta)(\tau_0-2\Delta)>2\kappa'^2\Delta^2/(1-\kappa'^2)`$.
The certificate evaluates both. For the contours of section 9,
$`\mu_Q`$ is between 0.017 and 0.019 and $`\tan(c_0+3e)\le0.62`$. On samples
the ratio $`|V^2|/|\operatorname{Im}x_{0i}^2|`$ near a cone is at most 0.003
(`corner_ratio`).

Measured on stress samples of 30 000 points in each region (generic points
at three scales, points within $`10^{-4}`$ to $`10^{-1}`$ of one cone, points
near two cones, and points near the four null lines on and off the
segments), next to the bounds proved in section 9:

| Property | $`(-,-)`$ | $`(+,-)`$ | $`(-,+)`$ | Proved bound |
| --- | --- | --- | --- | --- |
| Near-cone pairs tested | 19 520 | 19 369 | 19 441 | |
| Largest $`\operatorname{Im}x_{0i}^2`$ near a cone | $`-2\times10^{-17}`$ | $`-2\times10^{-17}`$ | $`-7\times10^{-17}`$ | $`\le0`$ |
| Largest excursion of $`\arg x_{0i}^2`$ beyond $`[-\pi,0]`$ | 0.0039 | 0.0039 | 0.0039 | 0.108, 0.074, 0.074 |
| Smallest $`\min_i\vert x_{0i}^2\vert /\min(1,\mathrm{dist}_{\rm edges}^2)`$ | $`3\times10^{-7}`$ | $`1\times10^{-7}`$ | $`2\times10^{-7}`$ | $`>0`$ |
| Smallest ratio in item 4, units of $`\nu`$ | 4 258 | 6 920 | 2 364 | $`\ge\tfrac12`$ |
| Largest $`\Vert \partial V\Vert `$ | 0.018 | 0.018 | 0.018 | 0.38, 0.18, 0.18 |

The imaginary parts near a cone are never positive; their size is set by
$`\nu`$ times the squared distance from the flat triangles, which is tiny for
points placed within $`10^{-4}`$ of a cone.

## 6. Consequences

**Absolute integrability.** The regulated integrand is absolutely integrable
on the contour for small $`\eta=-\epsilon>0`$, for weights of polynomial
growth in $`\log z`$ and for $`z^p`$ with $`|\operatorname{Re}p|<\eta`$. For
$`|\operatorname{Re}p|\ge\eta`$ it is not. The proof is section 10.

**The Fubini step of GATE_6_BOUNDS section 5.** On the contour
$`|z^p|\le|z|^{\operatorname{Re}p}e^{\Phi|\operatorname{Im}p|}`$ with
$`\Phi`$ the bound of item 2. For $`|\operatorname{Re}p|<\eta`$ the
integrand is absolutely integrable (Lemma 8′), and
$`|c(p)|\le Ce^{-(W-\delta')|\operatorname{Im}p|}`$ by Lemma 3. If
$`W-\delta'>\Phi`$ the double integral over the contour and over $`p`$
converges absolutely and the order can be exchanged. The lines of the
$`p`$ integral are then moved out to fixed positions; see steps A to D of
GATE_6_BOUNDS section 5.

**The displacement in Proposition 5.** Item 3 is that displacement: off the
edges the contour does not meet any light cone, so the integrand is bounded
there uniformly in $`\epsilon`$ near 0, and that part of the integral is
analytic in $`\epsilon`$. Poles come from neighbourhoods of the edges only.

**Equality with the Feynman integral.** The family $`x_0+i\lambda V`$,
$`0<\lambda\le1`$, satisfies items 1 and 3 for every $`\lambda`$, so no
singularity crosses the contour as $`\lambda`$ varies and the integral is
independent of $`\lambda`$; its limit $`\lambda\to0`$ is the $`-i0`$
prescription. This is Cauchy's theorem in several variables; the vanishing
of the contribution at infinity uses the decay $`|x_0|^{-8}`$.

## 7. Lemma 4 without Stirling

Lemma 4 of GATE_6_BOUNDS used Stirling's bound. The following suffices and
is elementary. From the Weierstrass product,

```math
\frac{|\Gamma(a+iv)|^2}{\Gamma(a)^2}=\prod_{n\ge0}\Big(1+\frac{v^2}{(a+n)^2}\Big)^{-1},
\qquad
\prod_{n\ge1}\Big(1+\frac{v^2}{n^2}\Big)=\frac{\sinh\pi v}{\pi v}.
```

Comparing factor by factor: for $`0<a\le1`$, $`a+n\le n+1`$, so

```math
|\Gamma(a+iv)|^2\le\Gamma(a)^2\,\frac{\pi v}{\sinh\pi v};
\qquad
\text{for }1\le a\le2:\quad
|\Gamma(a+iv)|^2\le\Gamma(a)^2\,(1+v^2)\,\frac{\pi v}{\sinh\pi v}.
```

Also $`\Gamma(1+p)\Gamma(1-p)=\pi p/\sin\pi p`$ exactly, and
$`|\Gamma(\tfrac\eta2+iy)|^2=|\Gamma(1+\tfrac\eta2+iy)|^2/(\eta^2/4+y^2)`$.

In the integrand of $`N_\eta(p)`$ the two $`p`$-dependent Gamma functions
have real parts in a compact subset of $`(0,2]`$, so their product is at
most $`C(1+v^2)^2\,\pi|v|/\sinh\pi|v|\le C'(1+|v|)^5e^{-\pi|v|}`$ with
$`v=\tau-y`$. The weight is at most
$`C(1+y^2)^2(\pi y/\sinh\pi y)^2/(\eta^2/4+y^2)^2`$, whose integral against
$`e^{\pi|y|}(1+|y|)^5`$ is at most $`C/\eta^3`$. Finally
$`\eta^2/\Gamma(2\eta)=2\eta^3/\Gamma(1+2\eta)\le2.3\,\eta^3`$. Therefore

```math
|N_\eta(\sigma+i\tau)|\le C\,(1+|\tau|)^3\,e^{\pi|\tau|},
```

uniformly in $`\eta`$. The power of $`|\tau|`$ is worse than in Lemma 4 and
does not matter: $`c(p)`$ decays like $`e^{-(W-\delta')|\tau|}`$ with
$`W-\delta'>\pi`$ (Lemma 3). The
test `test_elementary_gamma_bounds` checks the three inequalities on a grid.

## 8. Kinematic regions

Real null edges summing to zero have exactly one positive invariant among
$`s,t,u`$, so there are three regions. The construction is the same in all
of them, with the sign change of section 2 at the cusps whose diagonal is
timelike. In those regions two opposite cusps are timelike separated, and
the future cone of one meets the past cone of the other in a sphere; that
sphere lies inside the ball, where the two frame fields of the other two
cusps control both invariants. The table of section 5 and the certificate of
section 9 cover all three regions.

## 9. Proof of Proposition 6 with explicit constants

Norms $`|\cdot|`$ are Euclidean in $`\mathbb R^4`$, and
$`|a\!\cdot b|\le|a||b|`$ for the Minkowski product. Let $`c`$ be the centre
of the square, $`r=|x_0-c|`$, and

```math
\Delta=\max_i|x_i-c|,\qquad P=\max_j|p_j|,\qquad
A_1=\frac{2P}{\min_j|D_j|},\qquad A_2=1+2A_1P,\qquad
\lambda_0=16\kappa A_1P+\nu A_2 .
```

In the frame of $`x_j`$: $`|\alpha|,|\beta|\le A_1|x_0-x_j|`$,
$`|y_\perp|\le A_2|x_0-x_j|`$, and $`|d(P)|\le8|P|`$ because $`d`$ vanishes
at the vertex $`(0,0)`$ and has Lipschitz constant 8. Hence

```math
|V_j|\le\lambda_0|x_0-x_j|,\qquad\operatorname{Lip}(V_j)\le\lambda_0,
```

the second for the same reason. Write
$`n=x_0-x_i`$, $`u^\mp=x_i-c_\mp`$, and $`(\tau,\xi)`$ for the time and space
parts of $`x_0-c`$. The cutoff $`\chi`$ is the cubic smoothstep, with
$`|\chi'|\le\tfrac32`$. The inner field is used for $`r\le R`$, the outer for
$`r\ge2R`$, and $`V=(1-g)V_{\rm in}+gV_{\rm out}`$ in between with
$`|\nabla g|\le\tfrac3{2R}`$. $`\Theta_\pm`$ equals 1 for
$`\pm\tau\ge\tfrac12|\xi|`$ and 0 for $`\pm\tau\le\tfrac14|\xi|`$.

**(H1) A frame field is fully on near every cone.** If
$`\max_i|x_{0i}^2|\le\epsilon`$ then, in the frame of any cusp,
$`|\beta D|,|\alpha D|\le2\epsilon`$, $`|w|\le\epsilon`$,
$`\rho\le\epsilon+4\epsilon^2/|D|`$, and by the exact formula for the far
invariant

```math
\big|x_{0,j+2}^2-D'\big|\le\frac{4\epsilon|D'|}{|D|}+\epsilon
+2\sqrt{\frac{D'u}{D}\Big(\epsilon+\frac{4\epsilon^2}{|D|}\Big)} .
```

If the right side is less than $`|D'|-\epsilon`$ this contradicts
$`|x_{0,j+2}^2|\le\epsilon`$. Let $`\epsilon_*`$ be the largest $`\epsilon`$
for which it is. *Condition:* $`2\delta_0<\epsilon_*`$. Then for
$`|x_{0i}^2|\le\delta_0`$ one of the three weights
$`\chi(|x_{0,i\pm1}^2|/\delta_0)`$, $`\chi(|x_{0,i+2}^2|/\delta_0)`$ equals 1.
($`\epsilon_*`$ is weaker than the true minimum 0.199 but is proved.)

**(A) Inner field near a cone.** For $`|x_{0i}^2|\le\delta_0`$ the weight of
the frame $`i+2`$ vanishes, and by Lemma 5

```math
V_{\rm in}\!\cdot n=\sum_{j\ne i+2}m_jV_j\!\cdot n
\le-\sum_{j\ne i+2}m_j\Big(\frac{\kappa|D_j|}{2}\Psi_j+\nu\rho_j\Big)\le0 .
```

This is exact and needs no condition. By (H1) some $`m_j=1`$, so equality
forces $`x_0`$ into the flat triangle of that frame, where its three
invariants vanish only on the edges.

**(C0) Centres.** *Condition:* $`\tau_0\ge3\Delta`$. Then $`u^-`$ is future
timelike with $`u^-_0\ge\tau_0-\Delta`$ and $`|\vec u^-|\le\Delta`$.

**(C2) Near a cone outside the ball, $`\Theta=1`$.** If
$`|n^2|\le\delta_0`$ then $`|n_0|\ge|\vec n|-\sqrt{\delta_0}`$, so
$`|\tau|\ge|\xi|-(2\Delta+\sqrt{\delta_0})`$. With
$`\xi_0=2(2\Delta+\sqrt{\delta_0})`$: if $`|\xi|\ge\xi_0`$ then
$`|\tau|\ge\tfrac12|\xi|`$; if $`|\xi|<\xi_0`$ then
$`|\tau|\ge\sqrt{R^2-\xi_0^2}`$. *Condition:*
$`R\ge\tfrac{\sqrt5}{2}\xi_0`$. Then $`|\tau|\ge\tfrac12|\xi|`$ in both
cases.

**(C1) Outer field near a cone.** Take $`\tau\ge\tfrac12|\xi|`$, $`r\ge R`$,
$`|n^2|\le\delta_0`$. Then $`n_0\ge r/\sqrt{17}-\Delta`$ and
$`|\vec n|\le n_0+\sqrt{\delta_0}`$, so

```math
n^2+u^-\!\cdot n\ \ge\ -\delta_0+(\tau_0-2\Delta)\Big(\frac{R}{\sqrt{17}}-\Delta\Big)-\Delta\sqrt{\delta_0}.
```

*Condition:* the right side is positive. Then
$`V_{\rm out}\!\cdot n=-\kappa'(n^2+u^-\!\cdot n)<0`$; the past case is the
mirror image. With (A), $`\operatorname{Im}x_{0i}^2\le0`$ near every cone in
all three zones, strictly so where $`g>0`$.

**(C3) The spacelike wedge.** For $`|\tau|\le\tfrac14|\xi|`$ and $`r\ge R`$,
$`|\xi|\ge4R/\sqrt{17}`$ and
$`n^2\le(\tfrac14|\xi|+\Delta)^2-(|\xi|-\Delta)^2`$. *Condition:* this is
below $`-\delta_0`$ at $`|\xi|=4R/\sqrt{17}`$. So no cone comes within
$`\delta_0`$ where the outer field is switched off.

**(D) Phases away from the cones.** Frames $`j\ne i+2`$ contribute a
non-positive imaginary part everywhere, which keeps the argument in
$`[-\pi,0]`$. The frame $`i+2`$ acts only where $`|x_{0i}^2|\ge\delta_0`$ and
$`r\le2R`$, and there $`|2V_{i+2}\!\cdot n|\le2\lambda_0(2R+\Delta)^2`$, while
the real part moves by at most $`|V|^2\le v^2`$ with
$`v=4\lambda_0(2R+\Delta)+\kappa'(2R+\tau_0+\Delta)`$. *Condition:*
$`v^2<\delta_0/2`$. The excursion from this term is at most
$`\arctan[2\lambda_0(2R+\Delta)^2/(\delta_0-v^2)]`$.

For the outer field with $`\Theta=1`$,
$`x_{0i}^2\mapsto(1-i\kappa')^2n^2-2i\kappa'(1-i\kappa')\,n\!\cdot u-\kappa'^2u^2`$.
The imaginary part is positive only if $`n^2+n\!\cdot u<0`$. If
$`n\!\cdot u\ge0`$ then $`n^2<0`$, the imaginary part is at most
$`2\kappa'|n^2|`$ and the real part at most $`-(1-\kappa'^2)|n^2|`$: excursion
at most $`\arctan[2\kappa'/(1-\kappa'^2)]`$. If $`n\!\cdot u<0`$ then
$`|\vec n|\ge2n_0`$, $`|n^2|\ge\tfrac34|\vec n|^2`$ and
$`|n\!\cdot u|\le\Delta|\vec n|`$: excursion at most
$`\arctan[2\kappa'(1+2\Delta/(3(R/\sqrt{17}-\Delta)))]`$. Where
$`0<\Theta<1`$, $`|\tau|\le\tfrac12|\xi|`$, so
$`-n^2\ge\tfrac34|\xi|(|\xi|-4\Delta)`$ and $`|n|\le\tfrac{\sqrt5}{2}|\xi|+\Delta`$:
excursion at most
$`\arctan\{2\kappa'[1+(\tfrac{\sqrt5}{2}+\Delta/|\xi|)(\tau_0+\Delta)/(\tfrac34(|\xi|-4\Delta))]\}`$
with $`|\xi|\ge2R/\sqrt5`$. *Condition:* $`2R/\sqrt5>4\Delta`$. In the blend zone the
real part is $`n^2-V^2`$, of modulus at least $`|n^2|/2`$ because
$`v^2<\delta_0/2\le|n^2|/2`$ there, so the far term is at most twice the
pure outer one. The certificate therefore uses
$`\arctan[4\kappa'\max(1,\cdot,\cdot)]`$ for the far term, which covers
both zones; with it $`e\le0.11`$ for the contours evaluated. The total
excursion bound is the sum of the frame term and the largest of the three
far terms.

**(E) Local diffeomorphism.**
$`|\nabla m_j|\le\tfrac{3}{2\delta_0}\cdot2|x_0-x_{j+2}|\le3(2R+\Delta)/\delta_0`$,
so

```math
\|\partial V_{\rm in}\|\le4\Big[\frac{3(2R+\Delta)^2\lambda_0}{\delta_0}+\lambda_0\Big].
```

For the outer field $`|\nabla\Theta_\pm|\le6|\nabla(\tau/|\xi|)|\le7.5/r`$ on
its support, so
$`\|\partial V_{\rm out}\|\le\kappa'[1+7.5(1+(\tau_0+\Delta)/R)]`$. The blend
adds $`\tfrac3{2R}[4\lambda_0(2R+\Delta)+\kappa'(2R+\tau_0+\Delta)]`$.
*Condition:* the sum is below 1.

**(F) Lower bounds near a cusp.** The frame at the cusp has weight 1 there.
If $`\rho\ge\tfrac12|\alpha\beta D|`$ the imaginary part of $`x_{0,j}^2`$ is
at least $`2\nu\rho`$ in modulus; otherwise the real part is at least
$`\tfrac12|\alpha\beta D|-|V|^2`$. The other two invariants are treated in
the same way. This gives item 4 with $`c=\tfrac12`$ up to the negligible
$`|V|^2`$.

**(G) Lemma 6.** *Condition:* $`\theta\le1`$, with $`\theta`$ as in section
10.

**(S) The strips of Lemma 7.** With $`b_0=4\delta_0/\min|D|`$,
$`a_1=\delta_0/(8\max|D|)`$ and
$`\rho_1=\delta_0^2/(256\max(1,\max|D|\,|u|/\min|D|))`$: on a strip
$`|x_{0,j}^2|`$ and $`|x_{0,j+1}^2|`$ are at most
$`a_1\max|D|(1+a_1)+\rho_1`$, and

```math
|x_{0,j+2}^2|\ \ge\ \min|D|\,(b_0-a_1)-a_1\max|D|-\rho_1-2\sqrt{\rho_1\max|D|\,|u|/\min|D|}.
```

*Conditions:* the first is at most $`\delta_0`$; the second is at least
$`2\delta_0`$; $`b_0<\tfrac12`$; and
$`\Delta+a_1P+B\sqrt{\rho_1}\le R`$, so that the strips lie in the inner
zone.

**(I) Lemma 9.** $`\omega=c_u/2`$ is the fixed point of
$`c_u=(\tau_0-2\Delta)(1/\sqrt{17}-\Delta/(2R))-\Delta\sqrt{\omega/(2R)}`$.
*Conditions:* $`c_u>0`$;
$`\kappa'^2[2(2R+\Delta)(\tau_0+\Delta)+(\tau_0+\Delta)^2]\le\omega R/2`$
and $`2\kappa'^2(\tau_0+\Delta)\le\omega/4`$, which give
$`\kappa'^2(2|n||u|+|u|^2)\le\omega r/4`$ for all $`r\ge2R`$;
and $`2(2\Delta+\sqrt{2\omega R})\le4R/\sqrt5`$ with $`2R>5\omega/4`$, which
give $`|\tau|\ge\tfrac12|\xi|`$ whenever $`|n^2|<\omega r`$, $`r\ge2R`$.

**The conditions, evaluated.**

| Condition | $`(-,-)`$ | $`(+,-)`$, $`(-,+)`$ |
| --- | --- | --- |
| (H1) $`\epsilon_*`$ against $`2\delta_0=0.05`$ | 0.0897 | 0.345 |
| (C0) $`\tau_0\ge3\Delta`$, with $`\Delta`$ | 0.988 | 1.034 |
| (C1) margin (must be positive) | 5.55 | 5.19 |
| (C2), (C3), wedge margin | hold | hold |
| (D) excursion bound, rad | 0.108 | 0.074 |
| (E) bound on $`\Vert \partial V\Vert `$ | 0.378 | 0.183 |
| (G) $`\theta`$ | $`1.9\times10^{-5}`$ | $`7\times10^{-6}`$ |
| (S) $`b_0`$; margin of the far invariant over $`2\delta_0`$ | 0.092; 0.043 | 0.034; 0.041 |
| (I) $`c_u`$ | 0.921 | 0.884 |

`certificate` computes these from the geometry and the constants, and
`test_every_explicit_condition_of_the_proof_holds` asserts all of them in
each region. A further test checks that the measured excursion and Jacobian
norm are below the proved bounds.

With an earlier choice of constants ($`\kappa=\nu=3\times10^{-6}`$,
$`\delta_0=0.05`$) the sampled properties also held, but conditions (H1),
(D) and (E) did not: the bounds here are not sharp, and the constants were
reduced until every condition holds.

**Every square.** Up to Poincaré transformations and an overall scale, a
light-like square is determined by its kinematic region and one angle, that
is by the cross-ratio $`s/t`$. `tuned_contour` chooses the constants from
the geometry in a fixed order: $`\tau_0=8\Delta`$; $`\delta_0=\epsilon_*/2.5`$;
the smallest $`R`$ meeting (C1) to (C3); $`\kappa'`$ small enough for (D);
and $`\kappa=\nu`$ such that the bound in (E) equals 0.3. Each step can
always be carried out: $`\epsilon_*>0`$ for a non-degenerate square because
the right side of (H1) tends to zero with $`\epsilon`$; (C1) to (C3) hold for
large $`R`$; (D) and (E) hold for small amplitudes. So the conditions can be
met for every non-degenerate square. Evaluated on 30 angles in each region:

| Region | $`s/t`$ covered | $`\epsilon_*`$ | $`R`$ | $`\kappa=\nu`$ | (E) | (D), rad |
| --- | --- | --- | --- | --- | --- | --- |
| $`(-,-)`$ | 0.0057 to 173 | 0.0030 to 0.10 | 7.0 to 8.9 | $`3\times10^{-11}`$ to $`10^{-7}`$ | 0.300 | 0.057 to 0.069 |
| $`(+,-)`$ | $`-174`$ to $`-1.01`$ | 0.0031 to 0.62 | 8.0 to 8.9 | $`3\times10^{-11}`$ to $`10^{-6}`$ | 0.300 | 0.057 to 0.069 |
| $`(-,+)`$ | $`-0.994`$ to $`-0.0057`$ | 0.0031 to 0.62 | 8.0 to 8.9 | $`3\times10^{-11}`$ to $`10^{-6}`$ | 0.300 | 0.057 to 0.069 |

All conditions, including (G), (S) and (I), hold at all 90 squares, and for the same squares rescaled by
0.01 and by 37. The amplitudes shrink as the square degenerates
($`s/t\to0,\infty`$), as they must: the four-cone gap closes. The sampled
properties were also checked at four extreme squares.

**Interval arithmetic.** The table above is floating point at sample
squares. `certificate_iv.interval_certificate` evaluates the same
conditions with intervals, for all squares whose angle lies in a given
interval, with the contour constants fixed at those of its midpoint. Every
comparison is between the upper end of one interval and the lower end of
the other, so a pass is a statement about every square in the range, with
rounding included. Two bounds are replaced by weaker ones that imply them:
$`\arctan x\le x`$, and $`B^2`$ by the trace of the Gram matrix of the
transverse basis. In (I) any $`\omega`$ with $`2\omega\le c_u(\omega)`$
serves in Lemma 9, and $`\omega`$ is taken 2% below the fixed point.

| Region | Angle range | Boxes | $`s/t`$ covered | Result |
| --- | --- | --- | --- | --- |
| $`(-,-)`$ | 0.055 to 3.087 | 12 128 | 0.00076 to 1341 | all conditions hold on every box |
| $`(+,-)`$ | 0.02 to 3.087 | 12 268 | $`-1342`$ to $`-1.0001`$ | all hold |
| $`(-,+)`$ | 0.02 to 3.087 | 12 268 | $`-0.9999`$ to $`-0.00074`$ | all hold |

On these ranges $`e\le0.093`$ and $`\mu_Q\le0.0201`$, with
$`\tan(c_0+3e)\le0.66`$, so $`\mu_Q\tan(c_0+3e)\le0.014`$. For (H1) the value
of $`\epsilon`$ is taken 1% below the floating-point gap and the
contradiction is then verified in intervals at that value, which is all
(H1) needs. The exact transverse norm is
$`B^2=\max(1,|v|^2/(a\!\cdot b)^2)`$ in the notation below; the bound used
is larger by at most 1. An independent review re-ran the three full
covers, with every box passing, checked that boxes outside the ranges
fail, and compared the interval version with the floating-point
certificate line by line: no condition is missing and every comparison is
in the safe direction. An earlier version
covered only $`0.023\le s/t\le33.6`$ and the corresponding timelike ranges;
it failed further out because of the way a transverse basis was chosen,
not because of the conditions. The bound for $`B`$ is now in closed form:
the edges of this family have no component along one spatial axis, so that
axis is transverse with Euclidean norm 1, and the other transverse
direction is
$`v=(a_1b_3-a_3b_1,\ a_0b_3-a_3b_0,\ 0,\ a_1b_0-a_0b_1)`$ for edges $`a,b`$,
with $`-v^2=(a\!\cdot b)^2`$; hence $`B^2\le1+|v|^2/(a\!\cdot b)^2`$. Closer
to the degenerate squares the four-cone gap of (H1) closes and finer boxes
are needed; the limit itself, $`s/t\to0`$ or $`\infty`$, is not reached by
any finite cover. What interval arithmetic certifies is the list of
inequalities; that those inequalities imply the properties of the contour
is the content of the proofs above.

**Dimension.** By Proposition 7 of section 11 the same constants serve in
every dimension. As a check of that conclusion, 4000 points near the
$`D`$-dimensional cones were sampled in $`4+1`$, $`4+2`$ and $`4+3`$
dimensions in each region: the largest imaginary part near a cone is
negative in every case (between $`-2\times10^{-11}`$ and
$`-2\times10^{-10}`$) and the largest excursion of an argument is 0.0003.

## 10. Absolute integrability, line by line

Throughout, $`\eta=-\epsilon>0`$ is small, the weight satisfies
$`|f(L)|\le A(1+|L|)^k`$, and in $`D=4+2\eta`$ dimensions the transverse space
has dimension $`2+2\eta`$, so the transverse measure is
$`\rho^{\eta}d\rho`$ times a finite angular measure. The three invariants of
a frame depend on the transverse coordinates only through $`\rho`$; the
frame field is $`\nu y_\perp`$ in every transverse direction and Lemma 5 is
unchanged.

**Lemma 6.** Let $`r\le R`$ and $`|x_{0i}^2|\le\delta_0`$. Then

```math
|V^2|\le\frac\theta2\,\big|\operatorname{Im}x_{0i}^{2\,\rm def}\big|,
\qquad
\big|x_{0i}^{2\,\rm def}\big|\ \ge\ \tfrac12\max\Big(|x_{0i}^2|,\ \big|\operatorname{Im}x_{0i}^{2\,\rm def}\big|\Big),
\qquad
\theta=\max\Big(\frac{512\kappa P^2}{\min_j|D_j|},\ 8\nu B^2\Big),
```

where $`B`$ is the largest Euclidean norm of a transverse vector with
$`y_\perp^2=-1`$, and "def" denotes the value on the contour. Moreover
$`|\operatorname{Im}x_{0i}^{2\,\rm def}|\ge2\nu\rho_{j}`$ for any frame $`j`$
with $`m_j=1`$ that controls $`i`$.

*Proof.* Near the cone of $`x_i`$ the frame $`i+2`$ is off, so every active
frame controls $`i`$ and, with
$`b_j=\tfrac{\kappa|D_j|}{2}\Psi_j+\nu\rho_j`$, Lemma 5 gives
$`|\operatorname{Im}x_{0i}^{2\,\rm def}|=2|V\!\cdot n|\ge2\sum_jm_jb_j`$. On
the other hand $`|d_\beta p_{j-1}+d_\alpha p_j|\le\sqrt2P|d_j|`$ and
$`|y_\perp|\le B\sqrt{\rho_j}`$, so
$`|V_j|^2\le4\kappa^2P^2|d_j|^2+2\nu^2B^2\rho_j\le64\kappa^2P^2\Psi_j+2\nu^2B^2\rho_j`$
by the second part of Lemma 5. At most four frames are active and
$`m_j\le1`$, so

```math
|V^2|\le|V|^2\le4\sum_jm_j|V_j|^2
\le\sum_jm_j\Big[\frac{512\kappa P^2}{|D_j|}\cdot\frac{\kappa|D_j|}{2}\Psi_j+8\nu B^2\cdot\nu\rho_j\Big]
\le\theta\sum_jm_jb_j .
```

This is the first claim. The real part on the contour is
$`x_{0i}^2-V^2`$. If $`|\operatorname{Im}|\ge|x_{0i}^2|`$ the modulus is at
least $`|\operatorname{Im}|`$; otherwise it is at least
$`|x_{0i}^2|-\tfrac\theta2|\operatorname{Im}|\ge(1-\tfrac\theta2)|x_{0i}^2|`$.
With $`\theta\le1`$ both are at least half the maximum.
$`\blacksquare`$

For the contours used $`\theta`$ is between $`10^{-6}`$ and $`3\times10^{-4}`$.
On samples the ratio $`|V^2|/|\operatorname{Im}|`$ never exceeds
$`\theta/2`$ and the second ratio is 1.000.

**Strips.** In the frame of $`x_j`$ let

```math
S_j^+=\{-a_1\le\beta\le1-b_0,\ |\alpha|\le a_1,\ \rho\le\rho_1\},
```

and let $`S_j^-`$ be the same with the roles of $`\alpha`$ and $`-\beta`$
exchanged (the strip around the other edge at $`x_j`$). Choose
$`b_0=4\delta_0/\min_j|D_j|`$, and $`a_1,\rho_1`$ small enough that on the
strip $`|x_{0,j}^2|\le\delta_0`$, $`|x_{0,j+1}^2|\le\delta_0`$ and
$`|x_{0,j+2}^2|\ge2\delta_0`$ (conditions (S) of section 9); explicitly
$`a_1=\delta_0/(8\max|D|)`$ and
$`\rho_1=\delta_0^2/(256\max(1,\max|D|\,|u|/\min|D|))`$. Then
$`b_0<\tfrac12`$, so $`S_j^+`$ and $`S_{j+1}^-`$ together cover a
neighbourhood of the edge from $`x_j`$ to $`x_{j+1}`$, and the eight strips
cover a neighbourhood of all four edges.

**Lemma 7.** On $`S_j^+`$: the frame $`j`$ is fully on, and

```math
|x_{0,j}^{2\,\rm def}|\ge\nu\rho,\qquad
|x_{0,j+1}^{2\,\rm def}|\ge\frac\nu4\max\Big(\rho,\frac{|\alpha D|b_0}{2}\Big),\qquad
|x_{0,j-1}^{2\,\rm def}|\ge\frac\nu4\max\Big(\rho,\frac{|\beta D|}{2}\Big),\qquad
|x_{0,j+2}^{2\,\rm def}|\ge\delta_0 .
```

*Proof.* $`m_j=\chi(|x_{0,j+2}^2|/\delta_0)=1`$. The frames $`j-1`$ and
$`j+2`$ are off because $`x_{0,j+1}^2`$ and $`x_{0,j}^2`$ are within
$`\delta_0`$ of zero; the frame $`j+1`$ may be partly on and controls
$`x_{0,j}^2`$ and $`x_{0,j+1}^2`$. So Lemma 6 applies to $`i=j`$ and
$`i=j+1`$ with $`|\operatorname{Im}|\ge2\nu\rho`$. This gives the first
bound. For the second, $`x_{0,j+1}^2=\alpha D(\beta-1)-\rho`$: if
$`\rho\le\tfrac12|\alpha D|(1-\beta)`$ its modulus is at least
$`\tfrac12|\alpha D|b_0`$, otherwise $`2\nu\rho\ge\nu|\alpha D|b_0`$; with
Lemma 6 the modulus on the contour is at least half of
$`\max(|x_{0,j+1}^2|,2\nu\rho)`$. For the third, if
$`|x_{0,j-1}^2|\le\delta_0`$ the same argument applies with
$`x_{0,j-1}^2=\beta D(1+\alpha)-\rho`$; if not, the modulus on the contour is
at least $`\delta_0-v^2\ge\delta_0/2`$, which exceeds the right side because
$`\nu`$ is small. The fourth follows from $`|x_{0,j+2}^2|\ge2\delta_0`$ and
$`v^2<\delta_0/2`$. $`\blacksquare`$

On samples of the strips the four left sides exceed the right sides by
factors between 3.9 and $`5\times10^7`$.

**Lemma 8 (the strip integral).** For $`\eta'>0`$,

```math
\int_{S_j^+}\frac{\rho^{\eta'}\,d\alpha\,d\beta\,d\rho}
{\rho\,\max(\rho,|\alpha|)\,\max(\rho,|\beta|)}
\le C\int_0^{\rho_1}\rho^{\eta'-1}\big(1+\log\tfrac1\rho\big)^2d\rho<\infty .
```

*Proof.* $`\int_{-1}^{1}d\alpha/\max(\rho,|\alpha|)=2(1+\log(1/\rho))`$ for
$`\rho\le1`$, and the same for $`\beta`$. $`\blacksquare`$

**Consequence on the strips.** By Lemma 7 the integrand is at most
$`C\,\nu^{-3}\delta_0^{-1}\rho^{\eta}(1+|L|)^k/[\rho\max(\rho,|\alpha|)\max(\rho,|\beta|)]`$.
Each invariant on the contour lies between $`\nu\rho`$ and a constant in
modulus, so $`|L|\le C(1+\log(1/\rho))`$ and
$`(1+|L|)^k\le C\rho^{-\eta/2}`$. Lemma 8 with $`\eta'=\eta/2`$ gives
convergence.

**Lemma 8′ (powers $`1\pm p`$).** Let $`\sigma=\operatorname{Re}p`$. With
the powers $`1-\sigma`$ on $`x_{0,j}^2`$ and $`1+\sigma`$ on its two
neighbours, or the reverse,

```math
\int_{S_j^+}\frac{\rho^{\eta}\,d\alpha\,d\beta\,d\rho}
{\rho^{1\mp\sigma}\,\max(\rho,|\alpha|)^{1\pm\sigma}\,\max(\rho,|\beta|)^{1\pm\sigma}}
\le C_\sigma\int_0^{\rho_1}\rho^{\eta-1-|\sigma|}d\rho ,
```

which is finite if and only if $`|\sigma|<\eta`$.

*Proof.* For $`a>1`$,
$`\int_{-1}^1d\alpha/\max(\rho,|\alpha|)^a\le2a\rho^{1-a}/(a-1)`$; for
$`a<1`$ it is at most $`2/(1-a)`$. With the upper signs and $`\sigma>0`$
the two neighbours give $`\rho^{-2\sigma}`$ and the first factor
$`\rho^{-1+\sigma}`$. With $`\sigma<0`$ the neighbours give a constant and
the first factor $`\rho^{-1-|\sigma|}`$. The lower signs exchange the two
cases. $`\blacksquare`$

So for $`z^p`$ the strips converge when $`|\operatorname{Re}p|<\eta`$: by
Lemma 7 each factor $`|x^2|^{-(1\pm\sigma)}`$ is bounded by the
corresponding factor above, and $`|z^{i\tau}|\le e^{\Phi|\tau|}`$. The
bound cannot be improved by the sharper lower bounds of Proposition 6: in
the region $`\rho\sim|\alpha|\ll|\beta|\ll1`$ the three moduli are of the
order of $`\rho`$, $`|\alpha|`$ and $`|\beta|`$, and the integral of
$`|z|^{\sigma}`$ behaves like $`\int d\beta\,\beta^{\eta-1-\sigma}`$. The
divergence is real. Elsewhere in the ball and at infinity the argument
below is unchanged, since the powers still sum to 4 and are positive. The
test `test_strip_integral_with_powers_converges_only_below_the_regulator`
checks the threshold.

**The rest of the ball.** Let $`K`$ be the set $`r\le2R`$ with the open
strips removed. On $`K`$ no invariant vanishes: near a cone and inside $`R`$
by (A), where equality in Lemma 5 forces the point onto an edge; near a cone
in the blend zone by (C1); elsewhere because the modulus is at least
$`\delta_0-v^2`$. $`K`$ is compact and the invariants are continuous, so
they are bounded below on $`K`$ and the integrand is bounded.

**Lemma 9 (infinity).** For $`r\ge2R`$ and $`y=x_0-c`$,

```math
|x_{0i}^{2\,\rm def}|\ \ge\ c_\infty\,\kappa'\max\big(r,|y^2|\big).
```

*Proof.* On this zone $`V=-\kappa'\Theta(n+u)`$, so
$`V^2=\kappa'^2\Theta^2(n^2+2n\!\cdot u+u^2)`$ exactly. Let

```math
c_u=(\tau_0-2\Delta)\Big(\frac{1}{\sqrt{17}}-\frac{\Delta}{2R}\Big)-\Delta\sqrt{\frac{\omega}{2R}},
\qquad\omega=\tfrac12c_u ,
```

which is positive by condition (I) of section 9 (lengths are in units in
which $`\Delta`$ is of order one). If $`|n^2|\ge\omega r`$, the real part on the
contour is at least
$`|n^2|(1-\kappa'^2)-\kappa'^2(2|n||u|+|u|^2)\ge\tfrac12|n^2|`$, because
$`\kappa'^2(2|n||u|+|u|^2)\le\tfrac14\omega r`$ for small $`\kappa'`$. If
$`|n^2|<\omega r`$, the argument of (C2) with $`\delta_0`$ replaced by
$`\omega r`$ gives $`\Theta=1`$, and that of (C1) gives
$`n\!\cdot u\ge(\tau_0-2\Delta)(r/\sqrt{17}-\Delta)-\Delta\sqrt{\omega r}\ge c_ur`$,
so $`|\operatorname{Im}|=2\kappa'|n^2+n\!\cdot u|\ge2\kappa'(c_u-\omega)r=\kappa'c_ur`$.
In both cases the modulus is at least a constant times
$`\kappa'\max(r,|n^2|)`$. Finally $`|n^2-y^2|\le3\Delta r`$, so
$`\max(r,|n^2|)\ge\max(r,|y^2|)/(1+3\Delta)`$. $`\blacksquare`$

On samples the ratio is at least 5 in units of $`\kappa'`$.

With $`y=(\tau,\xi)`$, $`|y^2|=\big||\tau|-|\xi|\big|(|\tau|+|\xi|)\ge|\varsigma|\,r`$
where $`\varsigma=|\tau|-|\xi|`$. Hence the integrand is at most
$`C\,r^{-4}\max(1,|\varsigma|)^{-4}`$ times powers of $`\log r`$ and the
factor $`r^{2\eta}`$ from the transverse measure, and

```math
\int_{r\ge2R}\frac{d\tau\,d^3\xi\;r^{2\eta}(1+\log r)^k}{r^4\max(1,|\varsigma|)^4}
\le C\int^\infty\frac{\xi^2\,d\xi}{\xi^{4-2\eta}}(1+\log\xi)^k<\infty
```

for $`\eta<\tfrac12`$, the $`\tau`$ integral at fixed $`\xi`$ being bounded
by $`2\int d\varsigma/\max(1,|\varsigma|)^4`$.

Strips, the rest of the ball and infinity exhaust the contour, and the
Jacobian is bounded by (E). So the regulated integrand is absolutely
integrable. $`\blacksquare`$

## 11. The regulated integral on the contour

**Correction: the extra dimensions are displaced too.** In $`D=4+k`$
dimensions the point is $`(x_0,e)`$ with $`e\in\mathbb R^k`$, and every
invariant is $`X_i=x_{0i}^2-e\!\cdot e`$. An earlier version of this section
treated $`\mu=e\!\cdot e`$ as a real shift and used the four-dimensional
contour for every $`\mu`$. That is wrong once $`\mu>\delta_0`$: the frame
field of the opposite cusp is then switched on near the surface
$`x_{0i}^2=\mu`$ and can move it to the wrong side. On 19 865 points of that
surface at $`\mu=5\delta_0`$, 994 are displaced to the wrong side.

The remedy is to carry out the same construction in
$`\mathbb R^{1,3+k}`$. The transverse space of each frame gains the
coordinates $`e`$, so $`\rho\to\rho+e\!\cdot e`$; the frame field gains the
component $`\nu e`$, the far fields the component $`-\kappa'e`$; the cutoffs
use the $`D`$-dimensional invariants. Lemma 5, Lemmas P and Q of section 5 (with $`\rho`$ replaced by
$`\rho+\mu`$ and $`\vec n`$ including $`e`$) and sections 4, 9 and 10 hold
verbatim, with Euclidean norms taken in $`\mathbb R^{4+k}`$; section 10 was
already written with the transverse measure of this picture. With this
contour no point of any surface $`x_{0i}^2=\mu`$ is displaced to the wrong
side (0 of 19 865 at $`\mu=5\delta_0`$, and 0 at $`\mu=1,4,50`$), and the
survey of section 5 in $`4+1`$ dimensions gives the same results, with
near-cone points at $`\mu`$ up to $`4.5\times10^4`$. Non-integer $`k=2\eta`$
is reached in the usual way: the integrand and the field depend on $`e`$
through $`|e|`$ and a radial factor, so the integral reduces to one over
$`x_0`$ and $`m=|e|`$ with weight $`m^{2\eta-1}`$, which is analytic in
$`\eta`$.

**Proposition 7 (dimension).** Let $`k\ge1`$ and extend the field to
$`\mathbb R^{1,3+k}`$ as above: the frame fields gain $`\nu e`$, the far
fields $`-\kappa'\Theta_\pm e`$, and the cutoffs $`m_j`$, the blend $`g`$ and
$`\Theta_\pm`$ are evaluated with the $`D`$-dimensional invariants
$`x_{0i}^2-\mu`$, $`\mu=e\!\cdot e`$, the $`D`$-dimensional radius and the
space part $`(\xi,e)`$. Then every condition of section 9, (H1), (A),
(C0) to (C3), (D), (E), (F), (G), (S), (I), and (Q1), (Q3), (Q4) of
sections 5 and 11, holds with the same constants as for $`k=0`$.

*Proof.* The derivations of those conditions use only the following, each
of which is unchanged or improved.

- The points $`x_i`$, $`c`$, $`c_\pm`$ have $`e=0`$. So $`\Delta`$, $`P`$,
  $`D_j`$, $`u^\mp`$ and all products of these vectors are those of
  $`k=0`$, and $`n\!\cdot u`$ does not involve $`e`$.
- In the frame of $`x_j`$ the coordinates $`\alpha,\beta`$ do not involve
  $`e`$, and $`\rho`$ is replaced by $`\rho+\mu`$. The three invariants of
  the frame depend on the transverse position only through $`\rho+\mu`$,
  and Lemma 5 holds with $`\rho+\mu`$ by the same computation, since $`e`$
  is orthogonal to $`p_{j-1},p_j`$.
- The far invariant is
  $`D'(1+\alpha-\beta)+w+2\sqrt{\rho D'u/D}\cos\theta`$ with
  $`w=\alpha\beta D-\rho-\mu`$ and the square root involving the
  four-dimensional $`\rho\le\rho+\mu`$. The bound of (H1) and the strip
  bounds (S) use only an upper bound on that square root.
- Euclidean norms, with $`|x_0-x_j|`$ now the norm of $`(x_0-x_j,e)`$ in
  $`\mathbb R^{4+k}`$: $`|(y_\perp,e)|\le A_2|x_0-x_j|`$ because $`A_2\ge1`$, and
  the largest Euclidean norm $`B`$ of a transverse vector of unit Minkowski
  length is unchanged because $`B\ge1`$. Hence $`|V_j|\le\lambda_0|x_0-x_j|`$
  and $`\operatorname{Lip}(V_j)\le\lambda_0`$ still hold: the new block of
  $`\partial V_j`$ is $`\nu`$ times the identity and $`\nu\le\lambda_0`$.
- Gradients of the cutoffs: $`|\nabla x_{0i}^2|=2|x_0-x_i|`$ in every
  dimension, and $`|\nabla(\tau/|\xi|)|`$ is bounded as before with
  $`|\xi|`$ the full space part. So (E) is unchanged.
- The time and space decomposition used in (C0) to (C3), (I) and (Q3),
  (Q4) is that of $`\mathbb R^{1,3+k}`$: every inequality there is between
  $`\tau`$, $`|\xi|`$ (now including $`e`$), $`\Delta`$, $`\tau_0`$ and $`R`$.
- At most three frames are switched on near a cone, as before.

For non-integer $`k=2\eta`$ the integrand and the field depend on $`e`$
only through $`m=|e|`$, and the integral is the one over $`(x_0,m)`$,
$`m\ge0`$, with the weight $`m^{2\eta-1}dm`$. Its field and invariants are
those of $`k=1`$ restricted to $`e\ge0`$, so the case $`k=1`$ gives every
pointwise property, the bound on the Jacobian and the dominating
functions. $`\blacksquare`$

One thing the case $`k=1`$ does not give: an application of Stokes'
theorem in the variables $`(x_0,m)`$ has a boundary face at $`m=0`$, which
has to be shown to contribute nothing. Lemma 15 does this for the
integrand of Proposition 8; for the family used in CUSP_LOCAL_MODEL
section 13 it is done there.

This is a proof by inspection of the derivations in section 9. A second
automated review went through it condition by condition and agreed, and
confirmed that `field_extended` implements what is described. The
certificate evaluates the four-dimensional constants; by the proposition
they are the constants in every dimension. The survey in $`4+1`$
dimensions is a numerical check of the conclusion.

**Proposition 8 (parametric representation).** Let $`a_i`$ be complex with
$`\sum a_i=4`$, in one of two cases: all $`a_i=1`$, or $`a_1=a_3=1+p`$,
$`a_2=a_4=1-p`$ with $`|\operatorname{Re}p|<\eta`$. (Steps 1 to 5 hold for
any $`0<\operatorname{Re}a_i<2`$; step 6 has a dominating function only in
these two cases.) Let $`0<\eta<\tfrac12`$. The regulated integral of
$`\prod(-X_i)^{-a_i}`$ over the contour equals

```math
-\,i\pi^{2+\eta}\,\frac{\Gamma(2-\eta)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)\,
\prod\alpha_i^{a_i-1}\Big(-\sum_{i<k}\alpha_i\alpha_kx_{ik}^2\Big)^{\eta-2}.
```

*Derivation.* For $`0<\varphi\le\pi/2`$ replace every squared time
difference by $`e^{-2i\varphi}`$ times itself, in the integration variable
and in the external points alike. With $`\mu=e\!\cdot e`$,

```math
Y_i(\varphi)=\mu+|\xi-\xi_i|^2-e^{-2i\varphi}(\tau-t_i)^2 .
```

*Step 1 (a proper cone).* Each $`Y_i(\varphi)`$ is a non-negative
combination of $`1`$ and $`e^{i(\pi-2\varphi)}`$, so all four lie in a closed
cone of opening $`\pi-2\varphi<\pi`$, and
$`|Y_i(\varphi)|\ge\sin\varphi\,(\mu+|x_0-x_i|^2)`$. The integral over real
$`(x_0,e)`$ is absolutely convergent: the singularity at $`x_0=x_i`$,
$`\mu=0`$ is integrable because $`2\operatorname{Re}a_i<4`$, and the decay
at infinity is $`|x|^{-8}`$. For complex powers
$`|Y^{-a}|\le|Y|^{-\operatorname{Re}a}e^{\pi|\operatorname{Im}a|}`$ in the
cone, so every bound below holds with $`\operatorname{Re}a_i`$ and a
constant factor. The Feynman-parameter formula of step 2 holds for
$`\operatorname{Re}a_i>0`$.

*Step 2 (Feynman parameters).* For numbers in a common proper cone,
$`\prod Y_i^{-a_i}=\frac{\Gamma(4)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta(1-\sum\alpha)\prod\alpha_i^{a_i-1}(\sum\alpha_iY_i)^{-4}`$,
and $`|\sum\alpha_iY_i|\ge\sin\varphi\sum\alpha_i(\mu+|x_0-x_i|^2)`$, so
the double integral converges absolutely and the order can be exchanged.

*Step 3 (completing the square).*
$`\sum\alpha_iY_i=\mu+A_\alpha(\varphi)+|\xi-\bar\xi|^2-e^{-2i\varphi}(\tau-\bar t)^2`$
with $`A_\alpha(\varphi)=-\sum_{i<k}\alpha_i\alpha_kx_{ik}^2(\varphi)`$ and
$`x_{ik}^2(\varphi)=e^{-2i\varphi}(t_i-t_k)^2-|\xi_i-\xi_k|^2`$.

*Step 4 (rotating the time).* Put $`\tau-\bar t=e^{i\psi}\sigma`$ and turn
$`\psi`$ from $`0`$ to $`\varphi-\pi/2`$. The term
$`-e^{-2i\varphi}(\tau-\bar t)^2`$ has argument $`\pi+2\psi-2\varphi`$, which
moves from $`\pi-2\varphi`$ to $`0`$: it stays in the cone, the denominator
never vanishes, and the rotation is allowed. The integral is then Euclidean,
$`\int d^4y\,(c+y^2)^{-4}=\pi^2/(6c^2)`$, with Jacobian $`-ie^{i\varphi}`$.
Integrating the extra dimensions,
$`\frac{\pi^\eta}{\Gamma(\eta)}\int_0^\infty d\mu\,\mu^{\eta-1}(A+\mu)^{-2}=\pi^\eta\Gamma(2-\eta)A^{\eta-2}`$.
Hence, for $`\varphi>0`$,

```math
G(\varphi):=\int d^Dx\prod Y_i(\varphi)^{-a_i}
=-\,ie^{i\varphi}\,\pi^{2+\eta}\frac{\Gamma(2-\eta)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)
\prod\alpha_i^{a_i-1}A_\alpha(\varphi)^{\eta-2}.
\tag{*}
```

Every integral in steps 1 to 4 is absolutely convergent.

*Step 5 (the right side as $`\varphi\to0`$).* With spacelike diagonals
$`\operatorname{Re}A_\alpha(\varphi)\ge A_\alpha(0)=\alpha_1\alpha_3|s|+\alpha_2\alpha_4|t|\ge0`$,
because $`\cos2\varphi\le1`$. Since $`\eta-2<0`$,
$`|A_\alpha(\varphi)|^{\eta-2}\le A_\alpha(0)^{\eta-2}`$, which is integrable
against $`\prod\alpha_i^{a_i-1}`$ for $`\eta>0`$: this is the convergent
regulated Feynman-parameter integral. By dominated convergence the right
side of $`(*)`$ tends to the stated formula.

*Step 6 (the left side as $`\varphi\to0`$).* Three lemmas. Write
$`n=(x_0-x_i,e)`$, $`V=(V_4,V_e)`$ for the displacement in $`D`$ dimensions,
$`|\cdot|`$ for Euclidean norms and $`V\!\cdot n`$ for the Minkowski product;
$`Y_i(\varphi,\lambda)`$ is $`Y_i(\varphi)`$ at the point $`(x_0,e)+i\lambda V`$.
$`J<1`$ is the bound (E) on $`\|\partial V\|`$.

**Lemma 10 (no zeros, and a bound uniform in $`\varphi`$).**

(a) $`\operatorname{Im}\big(e^{i\varphi}Y_i(\varphi,\lambda)\big) =\sin\varphi\,\big(|n|^2-\lambda^2|V|^2\big)-2\lambda\cos\varphi\;V\!\cdot n`$.

(b) $`|V|\le J|n|`$ for every $`i`$.

(c) For $`0<\varphi\le\varphi_0=\delta_0/(8(2R+\Delta+1)^2)`$ and
$`0\le\lambda\le1`$, no $`Y_i(\varphi,\lambda)`$ vanishes.

(d) For $`0\le\varphi\le\varphi_0`$, at $`\lambda=1`$, in the inner zone
with $`|X_i|\le\delta_0`$:
$`|Y_i(\varphi,1)|\ge c_J\max\big(|X_i|,|\operatorname{Im}X_i^{\rm def}|\big)`$
with $`c_J=\tfrac12\min\big(\tfrac14,\tfrac{1-J^2}{16J}\big)`$.

*Proof.* (a) Expand
$`Y=|\vec n+i\lambda\vec V|^2-e^{-2i\varphi}(n_0+i\lambda V_0)^2`$, with
bilinear squares, multiply by $`e^{i\varphi}`$ and take the imaginary part:
the terms in $`\sin\varphi`$ collect to
$`|\vec n|^2-\lambda^2|\vec V|^2+n_0^2-\lambda^2V_0^2`$ and those in
$`\cos\varphi`$ to $`2\lambda(\vec n\!\cdot\!\vec V-n_0V_0)=-2\lambda V\!\cdot n`$.

(b) $`V`$ vanishes at each cusp $`x_i`$: the three frames that contain
$`x_i`$ in their flat triangle have $`V_j(x_i)=0`$, and the fourth is
switched off there. $`\|\partial V\|\le J`$, so
$`|V(x)|=|V(x)-V(x_i)|\le J|x-x_i|`$.

(c) By (a) and (b), if $`V\!\cdot n\le0`$ then
$`\operatorname{Im}(e^{i\varphi}Y)\ge\sin\varphi\,(1-J^2)|n|^2>0`$. Suppose
$`V\!\cdot n>0`$. In the inner and blend zones this happens only where
$`|X_i|\ge\delta_0`$, since near the cone every active field is admissible.
Write $`Y=-X_i^{\lambda}+(1-e^{-2i\varphi})(n_0+i\lambda V_0)^2`$, where
$`X_i^\lambda`$ is the invariant at the displaced point. The first term has
modulus at least $`\delta_0-v^2\ge\delta_0/2`$; the second at most
$`2\varphi(2R+\Delta+1)^2\le\delta_0/4`$. So $`|Y|\ge\delta_0/4`$. In the
outer zone $`V\!\cdot n=-\kappa'(n^2+n\!\cdot u)>0`$ means
$`n^2<-n\!\cdot u`$. Then
$`\operatorname{Re}Y=-(n^2-\lambda^2V^2)+(1-\cos2\varphi)(n_0^2-\lambda^2V_0^2)-2\lambda\sin2\varphi\,n_0V_0`$.
Here $`\Theta>0`$, so $`n_0\ge R/\sqrt{17}-\Delta>0`$ (the past case is the
mirror image). The last two terms are non-negative because
$`V_0=-\kappa'\Theta(n_0+u_0)`$ has the sign opposite to $`n_0`$ and
$`|V_0|<|n_0|`$ by (Q4). For the first, put
$`\varepsilon=\lambda\kappa'\Theta\le\kappa'`$, so that

```math
n^2-\lambda^2V^2=(1-\varepsilon^2)\,n^2-\varepsilon^2\,(2n\!\cdot u+u^2).
```

If $`n\!\cdot u\ge0`$ then $`n^2<-n\!\cdot u\le0`$ and $`u^2>0`$, so this is
negative. If $`n\!\cdot u<0`$ then $`|\vec n||\vec u|>n_0u_0`$, so
$`|\vec n|>n_0(\tau_0-\Delta)/\Delta\ge2n_0`$ by (C0). Hence
$`n^2\le-\tfrac34|\vec n|^2`$ and $`|n\!\cdot u|\le\Delta|\vec n|`$, and the
expression is at most
$`-\tfrac34(1-\varepsilon^2)|\vec n|^2+2\varepsilon^2\Delta|\vec n|`$, which
is negative because $`|\vec n|\ge2n_0`$ and (Q4). So
$`\operatorname{Re}Y>0`$.

(Q4): $`\kappa'(\tau_0+\Delta)<(1-\kappa')(R/\sqrt{17}-\Delta)`$ and
$`2(R/\sqrt{17}-\Delta)>8\kappa'^2\Delta/(3(1-\kappa'^2))`$; the certificate
evaluates it.

(d) Put $`X=X_i`$, $`I=|\operatorname{Im}X_i^{\rm def}|=2|V\!\cdot n|`$; here
$`V\!\cdot n\le0`$. From (a),
$`\operatorname{Im}(e^{i\varphi}Y)\ge\sin\varphi(1-J^2)|n|^2+\cos\varphi\,I`$.
So $`|Y|\ge\cos\varphi\,I\ge\tfrac12I`$. If $`I\ge\tfrac12|X|`$ this is at
least $`\tfrac14|X|`$. Otherwise use
$`\operatorname{Re}(e^{i\varphi}Y)=-\cos\varphi\operatorname{Re}X_i^{\rm def}+\sin\varphi\operatorname{Im}X_i^{\rm def}-4\sin\varphi\,n_0V_0`$.
By Lemma 6, $`|\operatorname{Re}X_i^{\rm def}|\ge|X|-\tfrac\theta2I\ge\tfrac34|X|`$,
and $`|4n_0V_0|\le4J|n|^2`$. If $`4J\sin\varphi|n|^2\le\tfrac18|X|`$ the real
part is at least
$`\tfrac38|X|-\tfrac14|X|\cdot\tfrac12-\tfrac18|X|=\tfrac18|X|`$ in modulus.
If not, $`\sin\varphi|n|^2>|X|/(32J)`$ and the imaginary part is at least
$`(1-J^2)|X|/(32J)`$. $`\blacksquare`$

**Lemma 11 (independence of $`\lambda`$).** For $`0<\varphi\le\varphi_0`$
the integral of $`\prod Y_i^{-a_i}`$ over $`(x_0,e)+i\lambda V`$ does not
depend on $`\lambda\in[0,1]`$.

*Proof.* The integrand is a holomorphic function $`h`$ on the open set
where no $`Y_i`$ vanishes, with its powers continued from $`\lambda=0`$,
where the $`Y_i`$ lie in the cone of step 1. The form $`h\,dz_1\cdots dz_D`$
is closed. The map $`(\lambda,x)\mapsto x+i\lambda V(x)`$ is continuously
differentiable (section 2) and by Lemma 10(c) takes $`[0,1]\times\mathbb R^D`$ into that set. Apply Stokes'
theorem on $`[0,1]\times\{|x|\le\rho\}`$ with a ball of radius $`r`$
removed around each point $`(x_i,0)`$. There $`Y_i`$ vanishes for every
$`\lambda`$, because $`V(x_i)=0`$; the integrand is of order
$`|n|^{-2\operatorname{Re}a_i}`$, the side of the removed tube has area of
order $`r^{D-1}`$ and the form contributes $`|V|\le Jr`$, so that side
contributes $`O(r^{D-2\operatorname{Re}a_i})\to0`$. The difference of the integrals
at $`\lambda=1`$ and $`\lambda=0`$ over the ball is the integral over the
side $`[0,1]\times\{|x|=\rho\}`$. On the side, for fixed $`\varphi>0`$, every
$`|Y_i|\ge c(\varphi)\rho^2`$: where $`V\!\cdot n\le0`$ by Lemma 10(a) and
(b); where $`V\!\cdot n>0`$, either
$`\sin\varphi(1-J^2)|n|^2\ge4\kappa'|n^2|`$ and (a) still gives
$`\tfrac12\sin\varphi(1-J^2)\rho^2`$, or $`|n^2|`$ is at least a multiple of
$`\sin\varphi\,\rho^2`$ and $`|Y|\ge\operatorname{Re}Y\ge\tfrac12|n^2|`$ by
the proof of (c). The side has area of order $`\rho^{D-1}`$, the form
contributes a factor $`|V|\le\kappa'(\rho+\tau_0+\Delta)`$, and
$`|h|\le C\rho^{-8}`$. The side integral is of order $`\rho^{D-8}`$ and tends
to zero. For non-integer $`D`$ the same argument applies to the reduced
variables $`(x_0,m)`$ with the factor $`m^{2\eta-1}`$, which vanishes with
the boundary at $`m=0`$. $`\blacksquare`$

**Lemma 12 (the limit).** As $`\varphi\to0`$, $`G(\varphi)`$ tends to the
integral of $`\prod(-X_i)^{-a_i}`$ over the contour.

*Proof.* By Lemma 11, $`G(\varphi)`$ equals the integral at $`\lambda=1`$.
There the integrand converges pointwise as $`\varphi\to0`$. It is dominated
uniformly in $`\varphi\le\varphi_0`$: near the cones in the inner zone by
Lemma 10(d), which gives the bound of Lemma 6 with $`c_J`$ in place of
$`\tfrac12`$, so that Lemma 7 holds with a different constant and Lemma 8
(unit powers) or Lemma 8′ (powers $`1\pm p`$, $`|\operatorname{Re}p|<\eta`$)
applies; near the cones in the blend zone, which is compact and contains no
edge, $`V\!\cdot n\le-c_1<0`$: it is continuous and strictly negative
there, by (C1) where $`g>0`$ and at $`r=R`$ by the strict case of (A),
since by (H1) a frame is fully on and its flat triangle lies inside
$`r\le\Delta<R`$. Lemma 10(a) then gives
$`|Y|\ge2\cos\varphi\,c_1`$; away
from the cones in the ball by $`|Y|\ge\delta_0/4`$; in the outer zone by
$`|Y|\ge\cos\varphi\cdot2\kappa'|n^2+n\!\cdot u|`$ where $`V\!\cdot n\le0`$
and by $`|Y|\ge\tfrac12|n^2|`$ otherwise, which reproduce Lemma 9. The
dominating function is integrable by section 10. $`\blacksquare`$

Steps 5 and 6 give the proposition. $`\blacksquare`$

*Timelike diagonal.* With $`s>0`$, $`A_\alpha(0)`$ vanishes inside the
simplex and the dominated-convergence argument of step 5 does not apply.
Lemma 14 of section 12 replaces it: the limit exists and is the
Mellin–Barnes form at $`s\to s-i0`$. Step 6 does not use the sign of the
diagonals.

*Numerical checks.*

| Check | Result |
| --- | --- |
| $`\vert Y_i\vert \ge\sin\varphi(\mu+\vert x_0-x_i\vert ^2)`$, $`0\le\arg Y_i\le\pi-2\varphi`$ | holds at 24 000 points, three angles |
| Both sides of $`(*)`$ at $`\varphi=\pi/2`$, $`\mu=0.5`$, unit powers | $`1.4317\pm0.0008`$ against $`1.43230`$ |
| Both sides at $`\varphi=1`$, unit powers | $`1.0225-1.3751i`$ against $`1.0235-1.3763i`$ |
| Both sides at $`\varphi=1`$, powers $`(1.2,0.8,1.2,0.8)`$ | $`1.1356-1.4812i`$ against $`1.1364-1.4820i`$ |
| Right side at $`\varphi=0.1,0.01,0.001`$ | approaches the Minkowski formula linearly |
| Lemma 10(a), relative error of the identity | $`6\times10^{-12}`$ |
| Lemma 10(b), largest $`\vert V\vert /\vert n\vert `$ | 0.0032, against $`J=0.30`$ |
| Lemma 10(c), smallest $`\vert Y_i\vert /(\sin\varphi\vert n\vert ^2)`$, also for $`\varphi`$ up to 0.3 | 1.00 |
| Lemma 10(d), smallest $`\vert Y_i\vert /\max(\vert X_i\vert ,\vert \operatorname{Im}X_i^{\rm def}\vert )`$ | 1.00, against $`c_J=0.095`$ |

The left side of $`(*)`$ is a four-dimensional quasi-Monte Carlo integral at
fixed $`\mu`$; the agreement is within its error. The checks of Lemma 10 are
in $`4+1`$ dimensions, in all three kinematic regions.

**Feynman parameters to Mellin–Barnes.** The Feynman-parameter integral for
$`z^p`$ (powers $`1\pm p`$) is evaluated directly, as three nested
one-dimensional integrals, and compared with the one-fold Mellin–Barnes
integral used in GATE_6_FUNCTIONAL and in Lemma 4:

| $`p`$ | $`\eta`$ | $`s/t`$ | Feynman parameters | Mellin–Barnes |
| --- | --- | --- | --- | --- |
| 0.2 | 0.4 | 1 | 3.520361875 | 3.520361876 |
| $`-0.3`$ | 0.3 | 2.5 | 3.893274… | 3.893274816 |
| 0 | 0.25 | 1 | 4.004157042 | 4.004157043 |

They agree to $`10^{-7}`$ or better, normalisation included.

So the regulated integral on the contour equals the Mellin–Barnes form.

## 12. Three remaining steps

### A. Feynman parameters to Mellin–Barnes

**Lemma 13.** Let $`a_1=a_3=1+p`$, $`a_2=a_4=1-p`$ with
$`|\operatorname{Re}p|<1`$, $`0<\eta<1-|\operatorname{Re}p|`$, and let
$`S,T`$ be non-zero complex numbers with $`|\arg S-\arg T|<\pi`$. Then

```math
\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)\prod\alpha_i^{a_i-1}
\big(\alpha_1\alpha_3S+\alpha_2\alpha_4T\big)^{\eta-2}
=\frac{1}{\Gamma(2-\eta)\Gamma(2\eta)}\int\frac{dw}{2\pi i}\,
\Gamma(-w)\Gamma(2-\eta+w)\,\Gamma(1+p+w)^2\,\Gamma(-1-p+\eta-w)^2\,S^{w}T^{\eta-2-w},
```

the contour being $`\operatorname{Re}w=-1-\operatorname{Re}p+\eta/2`$. With
the prefactor $`\Gamma(2-\eta)/\prod\Gamma(a_i)`$ this is the Mellin–Barnes
form of GATE_6_FUNCTIONAL with $`\epsilon=-\eta`$.

*Proof.* Fix $`0<c<\lambda`$ and put $`f(v)=e^{cv}(1+e^v)^{-\lambda}`$ for
real $`v`$. It is continuous, bounded and integrable, and the Beta integral
gives its Fourier transform:

```math
\hat f(y)=\int_{-\infty}^\infty f(v)e^{iyv}dv
=\int_0^\infty x^{c+iy-1}(1+x)^{-\lambda}dx=\frac{\Gamma(c+iy)\Gamma(\lambda-c-iy)}{\Gamma(\lambda)} ,
```

which is integrable in $`y`$ by the product bounds of section 7.

*Inversion, proved directly.* For $`\delta>0`$, by Fubini (the integrand is
absolutely integrable) and the Gaussian integral
$`\int e^{-\delta y^2}e^{iyx}dy=\sqrt{\pi/\delta}\,e^{-x^2/4\delta}`$,

```math
\frac1{2\pi}\int\hat f(y)e^{-iyv}e^{-\delta y^2}dy
=\int f(u)\,\frac{e^{-(u-v)^2/4\delta}}{\sqrt{4\pi\delta}}\,du .
```

As $`\delta\to0`$ the left side tends to
$`\frac1{2\pi}\int\hat f(y)e^{-iyv}dy`$ by dominated convergence, since
$`\hat f`$ is integrable. The right side is an average of $`f`$ against a
probability density that concentrates at $`v`$, and tends to $`f(v)`$
because $`f`$ is continuous and bounded. Hence
$`f(v)=\frac1{2\pi}\int\hat f(y)e^{-iyv}dy`$.

With $`e^v=A/B`$ and $`w=-c-iy`$ this reads

```math
(A+B)^{-\lambda}=\frac{1}{\Gamma(\lambda)}\int\frac{dw}{2\pi i}\,\Gamma(-w)\Gamma(\lambda+w)\,A^{w}B^{-\lambda-w},
\qquad-\lambda<\operatorname{Re}w<0,
```

for $`A,B>0`$. Both sides are analytic in $`A/B`$ on
$`|\arg(A/B)|<\pi`$: the right side converges there because
$`|\Gamma(-w)\Gamma(\lambda+w)|`$ decays like $`e^{-\pi|\operatorname{Im}w|}`$
while $`|A^wB^{-\lambda-w}|`$ grows at most like
$`e^{|\arg(A/B)||\operatorname{Im}w|}`$. So the formula holds for
$`|\arg(A/B)|<\pi`$.

Apply it with $`A=\alpha_1\alpha_3S`$, $`B=\alpha_2\alpha_4T`$,
$`\lambda=2-\eta`$. The parameter integral becomes a Dirichlet integral with
exponents $`b_1=b_3=1+p+w`$ and $`b_2=b_4=-1-p+\eta-w`$, whose sum is
$`2\eta`$:

```math
\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)\prod\alpha_i^{b_i-1}
=\frac{\Gamma(1+p+w)^2\,\Gamma(-1-p+\eta-w)^2}{\Gamma(2\eta)} .
```

On the contour all four exponents have real part $`\eta/2>0`$, and
$`-(2-\eta)<\operatorname{Re}w<0`$ because $`|\operatorname{Re}p|<1`$. The
exchange of the $`w`$ and $`\alpha`$ integrals is justified by absolute
convergence: the modulus of the $`\alpha`$ integrand is a Dirichlet
integrand with real exponents $`\eta/2`$, independent of
$`\operatorname{Im}w`$ apart from the factor
$`e^{|\arg(S/T)||\operatorname{Im}w|}`$, which the decay
$`e^{-\pi|\operatorname{Im}w|}`$ dominates. $`\blacksquare`$

The numerical agreement of the two forms to $`10^{-7}`$ (section 11) is now a
check of this lemma, not a substitute for it.

**Continuity up to $`|\arg(S/T)|=\pi`$.** On the contour the two Dirichlet
Gamma functions add a further decay: by the product bounds of section 7,
$`|\Gamma(\eta/2+iy)|^4\le Ce^{-2\pi|y|}(1+|y|)^{2+2\eta}(\eta^2/4+y^2)^{-2}`$.
The whole integrand therefore decays like $`e^{-3\pi|y|}`$ against a growth
$`e^{|\arg(S/T)||y|}`$. The Mellin–Barnes integral converges absolutely for
$`|\arg(S/T)|<3\pi`$ and is analytic there; in particular it is continuous
at $`|\arg(S/T)|=\pi`$, where the parameter integral itself is singular.

### B. The timelike diagonal

Let $`s>0`$, $`t<0`$. Work in the gauge $`\alpha_2+\alpha_4=1`$ (allowed
because the integrand is homogeneous of degree $`-4`$), with
$`X=\alpha_2\alpha_4`$. For $`\varphi>0`$,

```math
A_\alpha(\varphi)=A_\alpha(0)+\zeta\,Q_\alpha,\qquad
A_\alpha(0)=-s\,\alpha_1\alpha_3+|t|\,X,\qquad
Q_\alpha=\sum_{i<k}\alpha_i\alpha_k(t_i-t_k)^2,\qquad
\zeta=1-e^{-2i\varphi},
```

and $`\arg\zeta=\pi/2-\varphi`$. $`A_\alpha(0)`$ changes sign inside the
domain.

**Lemma 14.** As $`\varphi\to0`$ the right side of $`(*)`$ tends to the
Mellin–Barnes form of Lemma 13 at $`S=s\,e^{i\pi}`$, $`T=|t|`$, that is at
$`s\to s-i0`$.

*Proof.* $`Q_\alpha`$ is linear in $`\alpha_1`$: $`Q_\alpha=\alpha_1Q_1+Q_0`$
with $`Q_1>0`$, $`Q_0\ge0`$ independent of $`\alpha_1`$. Hence
$`A_\alpha(\varphi)=b_0-\alpha_1b_1`$ with $`b_0=|t|X+\zeta Q_0`$,
$`b_1=s\alpha_3-\zeta Q_1`$. Rotate the $`\alpha_1`$ contour:
$`\alpha_1=re^{-i\gamma}`$, $`0<\gamma<\pi/4`$. Then

```math
A_\alpha(\varphi)=|t|X+\zeta Q_0+r\,s\,\alpha_3\,e^{i(\pi-\gamma)}+r\,\zeta\,Q_1e^{-i\gamma}.
```

For $`\varphi\le\pi/4`$ the four terms have arguments $`0`$,
$`\pi/2-\varphi`$, $`\pi-\gamma`$ and $`\pi/2-\varphi-\gamma`$, all in
$`[0,\pi-\gamma]`$: a cone of opening less than $`\pi`$. Therefore

```math
|A_\alpha(\varphi)|\ \ge\ \sin(\gamma/2)\,\big(|t|X+r\,s\,\alpha_3\big),
```

uniformly in $`\varphi\in[0,\pi/4]`$. The rotation is legitimate for
$`\varphi>0`$: at every intermediate angle $`\gamma'\in[0,\gamma]`$ the four
terms have arguments in $`[0,\pi]`$ and the term $`r\zeta Q_1e^{-i\gamma'}`$
has a strictly positive imaginary part ($`Q_1>0`$ because the two ends of
the timelike diagonal have different times), so $`A_\alpha`$ does not
vanish; the weight $`\alpha_1^{p}(1+\alpha_1+\alpha_3)^{-2\eta}`$ is analytic in the
sector, and the integrand decays like $`r^{\operatorname{Re}p-2-\eta}`$ at
infinity and vanishes like $`r^{\operatorname{Re}p+1}`$ at the origin.

On the rotated contour the integrand converges pointwise as $`\varphi\to0`$
and is dominated by $`\sin(\gamma/2)^{\eta-2}`$ times the integrand of the
spacelike case with $`|s|`$ in place of $`-s`$, which is integrable for
$`\eta>0`$. So the limit is the integral of
$`\prod\alpha_i^{a_i-1}(|t|X-s\alpha_1\alpha_3)^{\eta-2}`$ over the rotated
contour.

The same rotation applied to the two-term integral of Lemma 13 with
$`S=se^{i\theta}`$, $`\theta<\pi`$, is legitimate for the same reason, and
by the same domination its limit as $`\theta\to\pi`$ is that same rotated
integral. By Lemma 13 and the continuity statement after it, that limit is
the Mellin–Barnes form at $`\theta=\pi`$. $`\blacksquare`$

With Lemmas 10 to 12, which do not use the sign of the diagonals,
Proposition 8 holds in all three regions, the right side being read at
$`s\to s-i0`$ or $`t\to t-i0`$ when that diagonal is timelike.

Numerically, at $`S/T=e^{i\pi}`$, $`p=0.2`$, $`\eta=0.4`$: the rotated
parameter integral is $`2.99934861+2.42990491i`$ with $`\gamma=0.5`$ and
$`2.99934849+2.42990482i`$ with $`\gamma=0.9`$; the Mellin–Barnes form is
$`2.99934849+2.42990482i`$. At $`S/T=2.5\,e^{i\pi}`$, $`p=-0.3`$,
$`\eta=0.3`$ they agree to $`2\times10^{-9}`$. The Mellin–Barnes form at
$`\theta=\pi-0.1,\pi-0.01,\pi-0.001`$ approaches its value at $`\pi`$
linearly. The imaginary part is the absorptive part that appears when a
diagonal is timelike.

### C. Non-integer dimension

For integer $`k`$ the integrand and the displacement are invariant under
rotations of $`e\in\mathbb R^k`$: $`V_4`$ depends on $`e`$ only through
$`m=|e|`$, and $`V_e=\sigma(x_0,m)\,e`$ with

```math
\sigma=(1-g)\,\nu\sum_jm_j-g\,\kappa'(\Theta_++\Theta_-),\qquad|\sigma|\le\max(4\nu,\kappa').
```

The map $`e\mapsto e(1+i\lambda\sigma)`$ preserves rays, so the integral over
$`e`$ reduces to one over $`m`$:

```math
\int d^ke\,(\cdots)=\frac{2\pi^{k/2}}{\Gamma(k/2)}\int_0^\infty(\cdots)\,w^{k-1}\,dw,
\qquad w=m\,(1+i\lambda\sigma).
```

**Definition.** For $`k=2\eta`$ the regulated integral over the contour is
the integral of the holomorphic five-form

```math
\Omega=\frac{2\pi^{\eta}}{\Gamma(\eta)}\;h(z,w)\;w^{2\eta-1}\,dw\wedge dz_1\wedge\cdots\wedge dz_4
```

over the image of $`(x_0,m)\mapsto(x_0+i\lambda V_4,\ m(1+i\lambda\sigma))`$ at
$`\lambda=1`$, where $`h`$ is the integrand with $`e\!\cdot e`$ replaced by
$`w^2`$ and $`w^{2\eta-1}`$ is the principal power on $`\operatorname{Re}w>0`$.
It agrees with the $`D`$-dimensional integral for integer $`k`$.

**Lemma 15.** For $`0<\varphi\le\varphi_0`$ and $`0<\eta<\tfrac12`$ the
integral of $`\Omega`$ does not depend on $`\lambda\in[0,1]`$.

*Proof.* $`\Omega`$ is a holomorphic form of top degree on the open set
$`\{\operatorname{Re}w>0,\ Y_i\ne0\}`$ of $`\mathbb C^5`$, hence closed, and
by Lemma 10(c), which was proved for points $`(x_0,e)`$ and uses $`e`$ only
through $`m`$, the map takes $`[0,1]\times\mathbb R^4\times(0,\infty)`$ into
that set; $`\operatorname{Re}w=m>0`$ because $`\sigma`$ is real. Apply
Stokes' theorem on
$`[0,1]\times\{|x_0|\le\rho\}\times[\varepsilon,M]`$. The boundary consists
of the two ends $`\lambda=0,1`$ and three sides.

- $`|x_0|=\rho`$: as in Lemma 11, $`|h|\le C(\rho^2+m^2)^{-4}`$ for fixed
  $`\varphi>0`$, so the side is at most
  $`C\rho^{3}\cdot\kappa'\rho\int_0^\infty m^{2\eta-1}(\rho^2+m^2)^{-4}dm=O(\rho^{2\eta-4})`$.
- $`m=M`$: at most $`C\,M^{2\eta-1}\cdot\kappa'M\int d^4x\,(|x|^2+M^2)^{-4}=O(M^{2\eta-4})`$.
- $`m=\varepsilon`$: on this side $`dw=i\varepsilon(\sigma\,d\lambda+\lambda\,d\sigma)`$,
  so the form is at most
  $`C\varepsilon^{2\eta-1}\cdot\varepsilon\,(|\sigma|+|\nabla\sigma||V_4|)\,|h|`$.
  For $`\varphi>0`$, $`|h|\le C\prod_i|x_0-x_i|^{-2a_i}`$ near the cusps by
  Lemma 10, which is integrable in four dimensions because $`2a_i<4`$, and
  $`|h|\le C|x|^{-8}`$ at infinity. The side is $`O(\varepsilon^{2\eta})`$.

All three tend to zero as $`\rho,M\to\infty`$ and $`\varepsilon\to0`$.
$`\blacksquare`$

Lemma 12 then applies with the same dominating function, the weight
$`m^{2\eta-1}`$ being the transverse measure already used in section 10. The
test `test_extended_field_is_radial_in_the_extra_dimensions` checks the
rotation invariance and that $`|\sigma|<0.01`$.

## 13. What remains

An internal review of this document found gaps, which are repaired above:
absolute integrability for the powers $`1\pm p`$ (Lemma 8′, and the
restriction $`|\operatorname{Re}p|<\eta`$), complex powers in
Proposition 8, the outer zone in Lemma 10(c), the excision of the cusp
points in Lemma 11, the blend zone in Lemma 12, and the phase bound with a
timelike diagonal (Lemma Q). Still open:

- **The phase bound with a timelike diagonal.** $`|\arg z|\le2\pi`$
  is proved at real points (Lemma P) and observed on the contour. On the
  contour $`3\pi-0.3`$ is proved (Lemma Q), which is all the argument
  uses. A bound close to $`2\pi`$ is now proved in
  [GATE_6_PHASE.md](GATE_6_PHASE.md), zone by zone: exactly $`2\pi`$
  for $`r\ge2R`$; $`2\pi+0.05`$ to $`2\pi+0.1`$ for $`r\le R`$; and
  $`2\pi+0.12`$ (default constants) to $`2\pi+0.19`$ in the blend zone,
  which sets the overall bound. It holds under further explicit
  inequalities between the constants, which are verified in interval
  arithmetic on the same cover as the certificate, with the bound at
  most $`2\pi+0.1908`$ throughout (`geocusp.phase_bound`,
  `geocusp.phase_bound_iv`), and it holds for the contour in $`D`$
  dimensions. The exact value $`2\pi`$ for $`r<2R`$ remains unproved.
  It appears to be true with a margin: where the phase is close to
  $`2\pi`$ the largest value found is $`2\pi-2\nu`$. But it is not true
  field by field (the frame at a cusp with a timelike diagonal pushes
  the wrong way at about one point in a hundred, and is outweighed by
  another frame), so the argument used here, which is linear in the
  field, cannot give it, and the proved constants are limited by the
  crude bound $`|V_j|\le\lambda_0|n_j|`$.
  Two things learned there correct
  earlier statements here. The route proposed before Lemma Q, a bound
  $`|\operatorname{Im}x_{0c}^2|\le c\,|\operatorname{Re}x_{0c}^2|`$
  beyond a cusp, is false: the ratio reaches $`10^{22}`$, harmlessly,
  because the neighbour's invariant turns by the same angle; what holds
  is a comparison of arguments, which is linear in the field. And in
  float64 the survey cannot resolve points within $`10^{-15}`$ of a null
  line through two adjacent cusps. The excursion bounds in the tables of
  this document were out of date and now show the values the certificate
  returns, 0.108 and 0.074. Nothing else changes: the
conditions that involve $`e`$, namely (D) itself ($`e<0.2`$) and (Q1),
are evaluated by the certificate, in floating point and in interval
arithmetic, from the value it computes and not from the tables, and
$`\tan(c_0+3e)=0.575`$ is within the bounds quoted in sections 5 and 9;
the strip condition $`2\pi+4e<W`$ holds with room to spare.
- **The certificate** is verified in interval arithmetic on the ranges of
  section 9 and in floating point outside them. The field itself is
  evaluated in floating point in the surveys: within about $`10^{-9}`$ of an edge the float64 field loses the sign of
  $`\operatorname{Im}x_{0i}^2`$ at the $`10^{-24}`$ level. This is rounding,
  not a property of the contour: `certificate_iv.edge_sign_check`
  re-evaluates the inner field in 100-digit arithmetic with exactly null
  edges at points $`10^{-6}`$ to $`10^{-30}`$ from the edges, and every
  near-cone imaginary part is negative, in all three regions.
- **The $`D`$-dimensional contour.** Proposition 7 of section 11 argues,
  condition by condition, that the certificate's constants are those of
  every dimension, and reduces non-integer $`k`$ to $`k=1`$. It is a proof
  by inspection; one automated review agreed with it. Theorem 13.3 of
  CUSP_LOCAL_MODEL, which identifies the contour integral with the Feynman
  integral, depends on it.
- **Review** by someone who works with these integrals.

What the argument takes from general analysis: Fubini's theorem and
dominated convergence, Cauchy's theorem in one variable, the ordinary
Stokes theorem for a continuously differentiable map, Fourier inversion,
the Beta and Dirichlet integrals, the Weierstrass product for the Gamma
function, and the identity theorem for analytic functions.

## 14. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_contour.py tests/test_schwinger.py
PYTHONPATH=src python -m examples.feynman_contour_survey
PYTHONPATH=src python -m geocusp.solver contour --points 5000
```
