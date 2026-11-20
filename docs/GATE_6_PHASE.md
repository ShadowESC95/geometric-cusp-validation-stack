# The phase bound with a timelike diagonal

[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md) proves $`|\arg z|\le2\pi+4e`$ on the
contour when both diagonals are spacelike, and $`|\arg z|\le3\pi-0.3`$ when
one is timelike (Lemma Q). It leaves open a bound close to $`2\pi`$ in the
second case. This document proves

```math
|\arg z|\ \le\ 2\pi+\max\big(b_{\rm in},\ b_{\rm bl},\ 0\big)
```

with explicit constants, the three entries being the bounds for $`r\le R`$,
$`R<r<2R`$ and $`r\ge2R`$, under explicit inequalities between the constants
of the contour (section 6). The right side is $`2\pi+0.12`$ for the default
constants and between $`2\pi+0.17`$ and $`2\pi+0.19`$ at ten tuned squares
evaluated in floating point; for $`r\ge2R`$ the bound is $`2\pi`$ exactly. In
interval arithmetic: for every square with a timelike diagonal there is a
contour of the construction satisfying all the conditions, with the
zone-by-zone bound at most $`2\pi+0.1908`$ (section 7: covers with the
constants of each box midpoint up to the angle 3.087, and a family of
constants scaled with the square from there to the degenerate end,
section 7.1).
For $`r<2R`$ the bound cannot be improved to $`2\pi`$ under these
conditions: there are admissible constants for which the phase exceeds
$`2\pi`$, by an amount of order $`\nu`$ (section 6).

**Status.** Proved, under the conditions listed in section 6. The
argument was found by an automated search and written out. It has been
through two rounds of automated review, each of which re-derived the steps
and tested them with its own high-precision implementation of the field.
The first round (an earlier version: 90 digits, eight squares, about 140
inequalities per point, no violation) found three small gaps and one wrong
constant, none affecting the result. After it the text was extended: a
uniform estimate in case 3d, the constants $`\Lambda_c,\Lambda_n`$,
Propositions R7 and R8, the zone-by-zone theorem and the interval
certificate. The second round reviewed the present text: no substantive
error; it reproduced the interval cover (24 536 boxes, no failure), checked
Proposition R8 at 90 digits on twelve squares (no violation), and asked for
the corrections that are now in the text (a trivial condition for
Proposition R7, the wording of what the interval certificate covers, one
number).

After the second round the blend zone was sharpened: the facts (F5) and
(F6), condition (B6), the constants $`s_0,s_1,s_2`$ and $`b_{\rm bl}`$,
Case N0 with its three sub-cases, the terms of the frames $`c+1`$ and
$`c+2`$ in N1 and N2, the end of Step 4, the zone-by-zone interval bound
and its cover. A further round of the same review covered exactly these
passages: proved as written; it re-derived each step, reproduced the full
cover (24 536 boxes, no failure, largest bound $`2\pi+0.1908`$), tested the
new inequalities at 90 digits on seven squares in four dimensions and on
five squares for $`k=1`$ (no violation), and asked for two corrections of
wording, now made. A separate round reviewed section 9 (the contour in
$`D`$ dimensions): correct as written, with its own $`D`$-dimensional field
and no violation on eleven runs.

After that the interval cover was extended to the angles from 0 to
$`\pi-10^{-7}`$ and Proposition R9 was added. A sixth round reviewed these:
the extended cover is correct and was reproduced; $`u=0`$ is a regular
point; Proposition R9 is correct as stated; but the inequality (6.1), which
an earlier version offered as the remaining step towards
$`|\arg z|\le2\pi`$ for $`r<2R`$, is false for admissible constants, and two
numbers were wrong. The counterexample and the corrections are in section
6. An earlier round re-derived directly the eight statements taken from
GATE_6_CONTOUR (list below): all hold; two written arguments there were
completed. Section 9 (the contour in $`D`$ dimensions) is by a separate
pass: 9.1 to 9.7 were rerun against the present constants and 9.8 proves
in $`D`$ dimensions the conditions taken from GATE_6_CONTOUR; a further
review found 9.8 correct and 9.3, 9.6, 9.7 consistent, with corrections
of wording. Section 7.1 (the degenerate squares, with constants scaled
with the square) was reviewed separately: proved as written, conditional
on the theorems; its corrections are made.

The statements taken from GATE_6_CONTOUR (list below) are no longer used
on trust: a review re-derived all eight in four dimensions, and 9.8
proves them in $`D`$ dimensions. The heading below is kept because this
document does not repeat those proofs.

Nothing here has been reviewed by a specialist.

**What is taken from GATE_6_CONTOUR without re-derivation:** (D)
($`\arg x_{0i}^2\in[-\pi-e,e]`$ on the contour); Lemma 5; Lemma 6
($`|V^2|\le\tfrac\theta2|\operatorname{Im}x_{0i}^2|`$ for $`r\le R`$,
$`|x_{0i}^2|\le\delta_0`$); (H1); (C1) with its margin $`M`$; (C2);
Proposition 6(3); and the inequality $`|V_{\rm in}|^2\le3\mu I_{\rm in}`$ from
the proof of Lemma Q.

The contour in $`D`$ dimensions is treated in section 9, by an audit of
the proof for its dependence on the dimension; that section has had no
second review.

**Range.** Every square with a timelike diagonal is covered, $`|s/t|=1`$
included (section 7, and 7.1 for the squares near the degenerate ones),
with the bound at most $`2\pi+0.1908`$ for a suitable contour at each
square.

**Not covered:** the degenerate squares themselves ($`|s/t|=0,\infty`$,
which are not squares); a bound on the integral uniform as the square
degenerates (the constants of the contour vanish there); the region
$`(-,-)`$, where $`2\pi+4e`$ is already proved.

Code: `src/geocusp/phase_bound.py`, `src/geocusp/phase_bound_iv.py`,
`src/geocusp/scaled_limit.py`.
Tests: `tests/test_phase_bound.py`.

## 0. Notation and conventions

Constants and conditions are those of GATE_6_CONTOUR section 9:
$`\Delta,P,A_1,A_2,\lambda_0,B,\theta,\mu,e`$, the conditions (H1), (C0)–(C3),
(D), (G). Norms $`|\cdot|`$ of vectors are Euclidean, and
$`|a\!\cdot b|\le|a||b|`$ for the Minkowski product.

Let $`c`$ be a cusp whose diagonal is timelike: $`D:=D_c=2p_{c-1}\!\cdot p_c>0`$.
Then $`D':=D_{c\pm1}<0`$ is the other diagonal and $`u=-D-D'\le0`$. The two
edges $`p_{c-1},p_c`$ at $`x_c`$ have the same time orientation, because
their product is positive. Indices are taken mod 4. The value $`u=0`$
occurs at $`|s/t|=1`$: there the quadrilateral is a parallelogram in a
plane, the four transverse planes coincide and $`B_j=1`$; nothing in the
proofs divides by $`u`$.

In the frame of $`x_c`$,

```math
y:=x_0-x_c=\alpha p_{c-1}+\beta p_c+y_\perp,\qquad \rho=-y_\perp^2\ge0,\qquad n_i:=x_0-x_i ,
```
```math
r_c=\alpha\beta D-\rho,\qquad r_{c-1}=r_c+A,\qquad r_{c+1}=r_c+B,\qquad A=\beta D,\quad B=-\alpha D ,
```

where $`r_i=x_{0i}^2`$ are the real invariants. Note $`AB=-D(r_c+\rho)`$.

On the contour

```math
X_i=r_i-V^2-iJ_i,\qquad J_i:=-2V\!\cdot n_i ,
```

and $`u_i:=\pi+\arg X_i`$ with the argument on the Feynman branch, so
$`u_i\in[-e,\pi+e]`$ by (D). A real spacelike invariant has $`u_i=0`$, a real
timelike one $`u_i=\pi`$. Put

```math
K_c:=u_c-u_{c-1}-u_{c+1}.
```

A field *controls* $`i`$ if it is a frame field $`V_j`$ with
$`j\in\{i-1,i,i+1\}`$; then its contribution to $`J_i`$ is $`\ge0`$ (Lemma 5).
The weight of $`V_j`$ is $`m_j=\chi(|r_{j+2}|/\delta_0)`$: the frame $`j`$ is
*off* if $`|r_{j+2}|\le\delta_0`$ and has weight 1 if $`|r_{j+2}|\ge2\delta_0`$.

The letter $`u`$ without index is the invariant $`-D-D'`$; $`u_i`$ with an
index are the angles just defined.

Throughout, $`\phi'(t)=0`$ for $`t\le0`$, $`t^2`$ on $`[0,1]`$, $`2t-1`$ for
$`t\ge1`$; so $`0\le\phi'(t)\le\min(t^2,2t)`$ for $`t\ge0`$. In the frame of a
cusp $`j`$ with coordinates $`(\alpha,\beta)`$, with $`\ell_3=\beta-\alpha-1`$,

```math
d_\beta=\phi'(\ell_3)-\phi'(-\beta),\qquad d_\alpha=\phi'(\alpha)-\phi'(\ell_3).
```

## 1. Reduction

**Lemma R1.** If $`s>0`$ then

```math
\arg z=\pi+K_2+u_4=\pi+K_4+u_2\ \le\ 2\pi+e+\min(K_2,K_4),\qquad \arg z\ge-\pi-4e .
```

If $`t>0`$ then $`-\arg z=\pi+K_1+u_3=\pi+K_3+u_1\le2\pi+e+\min(K_1,K_3)`$ and
$`\arg z\le\pi+4e`$. At real points where no invariant vanishes,
$`K_c\le0`$.

*Proof.* $`\arg z=\pi(\tau_s-\tau_t)+a_2+a_4-a_1-a_3`$ and $`a_i=u_i-\pi`$, so
for $`s>0`$: $`\arg z=\pi+u_2+u_4-u_1-u_3`$. Use $`u_i\le\pi+e`$ and
$`u_i\ge-e`$. For $`t>0`$ the sign of the constant and the roles of the two
pairs are exchanged. At real points $`K_c=\pi(\tau_c-\tau_{c-1}-\tau_{c+1})`$;
if $`\tau_c=1`$ then $`r_c>0`$, so $`AB<0`$, so one of $`r_{c\pm1}`$ exceeds
$`r_c>0`$. $`\blacksquare`$

So an upper bound on $`K_c`$ at **one** cusp with a timelike diagonal, at
every point of the contour, gives the phase bound.

## 2. The pair criterion

**Lemma R2.** Put

```math
N_-:=2(y\!\cdot p_{c-1})\,y-y^2\,p_{c-1},\qquad N_+:=y^2\,p_c-2(y\!\cdot p_c)\,y .
```

1. $`N_-=r_{c-1}n_c-r_cn_{c-1}`$ and $`N_+=r_{c+1}n_c-r_cn_{c+1}`$.
2. $`N_\pm`$ are null. $`N_-`$ is zero or has the time orientation of
   $`p_{c-1}`$; $`N_+`$ is zero or has the orientation opposite to $`p_c`$.
3. For every real vector $`V`$ (with $`J_i=-2V\!\cdot n_i`$ and
   $`X_i=r_i-V^2-iJ_i`$):

```math
r_{c\mp1}J_c-r_cJ_{c\mp1}=-2V\!\cdot N_\mp,\qquad
\operatorname{Im}\big(X_{c\mp1}\overline{X_c}\big)=(r_{c\mp1}-V^2)J_c-(r_c-V^2)J_{c\mp1}
=-2V\!\cdot N_\mp+V^2\,(J_{c\mp1}-J_c).
```

4. Let $`Z_1=\varrho_1-i\jmath_1`$, $`Z_2=\varrho_2-i\jmath_2`$ be non-zero
   with $`\jmath_1,\jmath_2\ge0`$ and
   $`\operatorname{Im}(Z_1\overline{Z_2})=\varrho_1\jmath_2-\varrho_2\jmath_1\ge0`$.
   Then $`\arg Z_1\ge\arg Z_2`$ (arguments in $`[-\pi,0]`$), unless
   $`\jmath_1=\jmath_2=0`$.

*Proof.* (1) $`n_{c-1}=y+p_{c-1}`$, $`r_{c-1}=y^2+2y\!\cdot p_{c-1}`$, so
$`r_{c-1}y-y^2(y+p_{c-1})=2(y\!\cdot p_{c-1})y-y^2p_{c-1}`$. Likewise
$`n_{c+1}=y-p_c`$, $`r_{c+1}=y^2-2y\!\cdot p_c`$.
(2) $`N_-^2=4(y\!\cdot p_{c-1})^2y^2-4(y\!\cdot p_{c-1})^2y^2=0`$, and
$`N_-\!\cdot p_{c-1}=2(y\!\cdot p_{c-1})^2\ge0`$. A non-zero null vector whose
product with a null vector $`p`$ is positive has the orientation of $`p`$.
If $`y\!\cdot p_{c-1}=0`$ then $`y`$ is orthogonal to a null vector, so
$`y^2\le0`$ and $`N_-=|y^2|p_{c-1}`$. The same for $`N_+`$ with
$`N_+\!\cdot p_c=-2(y\!\cdot p_c)^2\le0`$.
(3) from (1) and $`J_i=-2V\!\cdot n_i`$.
(4) both arguments lie in $`[-\pi,0]`$ and
$`\sin(\arg Z_1-\arg Z_2)\ge0`$ with $`\arg Z_1-\arg Z_2\in[-\pi,\pi]`$; the
value $`-\pi`$ requires $`Z_1<0<Z_2`$ real. $`\blacksquare`$

The quantity $`-2V\!\cdot N_\mp`$ is linear in $`V`$. So the statement
"$`u_c\le u_{c\mp1}`$" can be checked field by field, like admissibility.

**Lemma R3 (one frame field, closed form).** Let $`V=V_j`$ be the frame
field of the cusp $`j`$ alone, $`\sigma=\operatorname{sign}D_j`$,
$`(\alpha,\beta,\rho)`$ the coordinates in that frame. Then
$`J_i=\kappa|D_j|\,d\!\cdot(P-v_i)+2\nu\rho`$ for $`i=j-1,j,j+1`$, with
$`v_j=(0,0)`$, $`v_{j-1}=(0,-1)`$, $`v_{j+1}=(1,0)`$ in coordinates
$`(\beta,\alpha)`$. Hence

```math
J_{j-1}=J_j+\kappa|D_j|\,d_\alpha,\qquad J_{j+1}=J_j-\kappa|D_j|\,d_\beta ,
```
```math
r_{j-1}J_j-r_jJ_{j-1}=|D_j|\Big[\sigma\kappa|D_j|\,\beta^2d_\beta+\rho\,\big(2\nu\sigma\beta+\kappa d_\alpha\big)\Big],
```
```math
r_{j+1}J_j-r_jJ_{j+1}=-|D_j|\Big[\sigma\kappa|D_j|\,\alpha^2d_\alpha+\rho\,\big(2\nu\sigma\alpha+\kappa d_\beta\big)\Big].
```

Moreover $`V_j^2=\kappa^2D_j\,d_\beta d_\alpha-\nu^2\rho`$.

*Proof.* $`V_j=-\sigma\kappa(d_\beta p_{j-1}+d_\alpha p_j)+\nu y_\perp`$, and by
the computation in Lemma 5,
$`V_j\!\cdot(x_0-x_i)=-\sigma\kappa\tfrac{D_j}2d\!\cdot(P-v_i)-\nu\rho`$.
This gives $`J_i`$. With $`r_j=\alpha\beta D_j-\rho`$,
$`r_{j-1}=r_j+\beta D_j`$:

```math
r_{j-1}J_j-r_jJ_{j-1}=\beta D_jJ_j-r_j\kappa|D_j|d_\alpha
=\beta D_j\big(\kappa|D_j|(\beta d_\beta+\alpha d_\alpha)+2\nu\rho\big)-(\alpha\beta D_j-\rho)\kappa|D_j|d_\alpha ,
```

and the terms in $`\alpha\beta d_\alpha`$ cancel. With $`r_{j+1}=r_j-\alpha D_j`$:

```math
r_{j+1}J_j-r_jJ_{j+1}=-\alpha D_jJ_j+r_j\kappa|D_j|d_\beta
=-\alpha D_j\big(\kappa|D_j|(\beta d_\beta+\alpha d_\alpha)+2\nu\rho\big)+(\alpha\beta D_j-\rho)\kappa|D_j|d_\beta ,
```

and the terms in $`\alpha\beta d_\beta`$ cancel. The square uses
$`p_{j-1}^2=p_j^2=0`$, $`2p_{j-1}\!\cdot p_j=D_j`$, $`y_\perp^2=-\rho`$.
$`\blacksquare`$

## 3. The outer zone

**Proposition R4.** Assume (C0) ($`\tau_0\ge3\Delta`$),
$`2R>\sqrt2\,\tau_0`$, $`e<\tfrac\pi2`$ and $`\epsilon_\kappa<1`$ (defined below). Then
for $`r\ge2R`$

```math
K_c\le\max(e,\,3e_{\rm sp}),\qquad
\tan e_{\rm sp}=\frac{2\kappa'\big(1+\sqrt2\,\tfrac{\tau_0+\Delta}{\tau_0-2\Delta}\big)}{1-\epsilon_\kappa},\qquad
\epsilon_\kappa=\kappa'^2\Big(1+2\sqrt2\,\tfrac{\tau_0+\Delta}{\tau_0-2\Delta}+\tfrac{(\tau_0+\Delta)^2}{(\tau_0-2\Delta)^2}\Big).
```

*Proof.* Write $`x_0-c=(\tau,\vec\xi)`$, and $`\xi_-:=x_0-c_-`$,
$`\xi_+:=x_0-c_+`$. For $`r\ge2R`$ the field is
$`V=-\kappa'\Theta_+\xi_--\kappa'\Theta_-\xi_+`$, and $`\Theta_+\Theta_-=0`$.

*(i) $`x_0\in J^+(c_-)`$, i.e. $`\tau+\tau_0\ge|\vec\xi|`$.* If $`\tau<0`$
then $`|\tau|+|\vec\xi|\le\tau_0`$ and $`r\le\tau_0<2R`$. So $`\tau\ge0`$,
$`\Theta_-=0`$ and $`V=-\kappa'\vartheta\,\xi_-`$ with $`\vartheta=\Theta_+\in[0,1]`$ and
$`\xi_-`$ future causal. Then $`V^2=\kappa'^2\vartheta^2\xi_-^2\ge0`$.

Suppose first that the edges at $`c`$ are future directed, and take the
neighbour $`{\rm nb}=c-1`$. By Lemma R2(2) $`N_-`$ is future null or zero,
so $`-2V\!\cdot N_-=2\kappa'\vartheta\,\xi_-\!\cdot N_-\ge0`$. Also
$`J_{c-1}-J_c=-2V\!\cdot(n_{c-1}-n_c)=-2V\!\cdot p_{c-1}=2\kappa'\vartheta\,\xi_-\!\cdot p_{c-1}\ge0`$.
By Lemma R2(3),

```math
\operatorname{Im}\big(X_{c-1}\overline{X_c}\big)\ge0 .
```

If the edges at $`c`$ are past directed take $`{\rm nb}=c+1`$: $`N_+`$ is
future null or zero, and
$`J_{c+1}-J_c=2V\!\cdot p_c=-2\kappa'\vartheta\,\xi_-\!\cdot p_c\ge0`$. Again
$`\operatorname{Im}(X_{c+1}\overline{X_c})\ge0`$.

*Signs of the imaginary parts.* If $`r_i>0`$ then $`n_i`$ is future
timelike. Otherwise $`x_0\in J^-(x_i)\cap J^+(c_-)`$; then
$`|\vec\xi|\le\tau+\tau_0`$ and $`|\vec\xi|-\Delta\le\Delta-\tau`$, so
$`-\tau_0\le\tau\le2\Delta\le\tau_0`$ and $`|\vec\xi|\le\tau_0/2+\Delta\le\tau_0`$,
hence $`r\le\sqrt2\tau_0<2R`$. So $`J_i=2\kappa'\vartheta\,n_i\!\cdot\xi_-\ge0`$
whenever $`r_i>0`$. If $`r_i\le0`$ then $`\operatorname{Re}X_i=r_i-V^2\le0`$.
In both cases $`a_i:=\arg X_i\in[-\pi-e,0]`$.

*Conclusion.* $`\sin(a_{\rm nb}-a_c)\ge0`$ with
$`a_{\rm nb}-a_c\in[-\pi-e,\pi+e]`$, so either $`a_{\rm nb}\ge a_c`$ or
$`a_c\ge a_{\rm nb}+\pi`$. In the second case $`a_c\ge-e`$ and
$`a_{\rm nb}\le-\pi`$. Then $`\operatorname{Re}X_c>0`$, so $`r_c>V^2\ge0`$,
$`n_c=y`$ is future timelike, and
$`r_{\rm nb}=r_c+2|y\!\cdot p|\ge r_c`$ ($`p=p_{c-1}`$ future or $`p=p_c`$
past). So $`\operatorname{Re}X_{\rm nb}\ge\operatorname{Re}X_c>0`$ and
$`\operatorname{Im}X_{\rm nb}\le0`$, i.e. $`a_{\rm nb}\in(-\tfrac\pi2,0]`$, a
contradiction. Hence $`u_c\le u_{\rm nb}`$, and with the third index
$`u\ge-e`$: $`K_c\le e`$. No $`X_i`$ vanishes here (Proposition 6(3)).

*(ii) $`x_0\in J^-(c_+)`$.* The mirror image under time reversal: the field
is $`-\kappa'\vartheta\,\xi_+`$ with $`\xi_+`$ past causal, $`V^2\ge0`$, and the
neighbour is $`c+1`$ for future edges, $`c-1`$ for past edges (the one for
which $`N_\pm`$ is past null). The same computation gives $`K_c\le e`$.

*(iii) Neither.* Then $`\tau+\tau_0<|\vec\xi|`$ and $`\tau_0-\tau<|\vec\xi|`$,
i.e. $`|\tau|<|\vec\xi|-\tau_0`$. For every $`i`$,
$`|n_i^0|\le|\tau|+\Delta<|\vec n_i|-(\tau_0-2\Delta)`$. Hence

```math
-r_i=(|\vec n_i|-|n_i^0|)(|\vec n_i|+|n_i^0|)\ \ge\ (\tau_0-2\Delta)\,|\vec n_i|,\qquad|\vec n_i|\ge\tau_0-2\Delta,\qquad|n_i|\le\sqrt2|\vec n_i| .
```

The field is $`\vartheta W`$ with $`W=-\kappa'(n_i+w_i)`$, $`w_i:=x_i-c_\mp`$,
$`|w_i|\le\tau_0+\Delta`$, $`0\le\vartheta\le1`$. So

```math
|\operatorname{Im}X_i|=2|V\!\cdot n_i|\le2\kappa'\big(|r_i|+\sqrt2|\vec n_i|(\tau_0+\Delta)\big)\le2\kappa'|r_i|\Big(1+\sqrt2\,\frac{\tau_0+\Delta}{\tau_0-2\Delta}\Big),
```
```math
|V^2|\le\kappa'^2\big(|r_i|+2|n_i||w_i|+|w_i|^2\big)\le\epsilon_\kappa|r_i| .
```

So $`\operatorname{Re}X_i\le-(1-\epsilon_\kappa)|r_i|<0`$ and
$`|u_i|\le e_{\rm sp}`$ for all $`i`$; $`K_c\le3e_{\rm sp}`$.
$`\blacksquare`$

**Proposition R7.** Assume (R4), $`4e_{\rm sp}<\pi`$ and

```math
\text{(R7)}\qquad\Big(R-\frac{\tau_0}2\Big)(\tau_0-2\Delta)\ \ge\ 2\Delta\,\tau_0 .
```

Then for $`r\ge2R`$: $`\arg z\le2\pi`$ if $`s>0`$, and $`\arg z\ge-2\pi`$ if
$`t>0`$.

*Proof.* Take $`s>0`$; let $`c`$ and $`c+2`$ be the cusps with the timelike
diagonal, $`m`$ and $`m'`$ the other two, so that
$`\arg z=\pi+u_c+u_{c+2}-u_m-u_{m'}`$.

*In $`J^+(c_-)`$.* As in R4(i), $`V=-\kappa'\vartheta\,\xi_-`$,
$`V^2\ge0`$, and every $`a_i=u_i-\pi\le0`$. One of the two cusps, say $`c`$,
has future directed edges; then $`p_{c+1}`$ is past directed (its product
with $`p_c`$ is $`D'/2<0`$), so the edges at $`c+2`$ are past directed. By
R4(i), $`u_c\le u_{c-1}`$ and $`u_{c+2}\le u_{c+3}`$. These are the same
neighbour: put $`m=c-1`$. Hence

```math
\arg z\le\pi+u_m-u_{m'} .
```

Write $`w_i=x_i-c_-=(\tau_0+t_i,\vec x_i)`$ and $`n_i=(\tau-t_i,\vec\xi-\vec x_i)`$.
From $`r\le\tau+|\vec\xi|\le2\tau+\tau_0`$: $`\tau\ge R-\tau_0/2`$. Then, with
$`|t_i|,|\vec x_i|\le\Delta`$ and $`|\vec\xi|\le\tau+\tau_0`$,

```math
n_i\!\cdot w_i\ \ge\ (\tau-\Delta)(\tau_0-\Delta)-(\tau+\tau_0+\Delta)\Delta=\tau(\tau_0-2\Delta)-2\Delta\tau_0\ \ge0
```

by (R7). Since $`J_i=2\kappa'\vartheta\,n_i\!\cdot\xi_-=2\kappa'\vartheta\,(r_i+n_i\!\cdot w_i)`$:

- (a) if $`J_i<0`$ then $`r_i<0`$ and $`|J_i|\le2\kappa'\vartheta|r_i|\le2\kappa'\vartheta\,|\operatorname{Re}X_i|`$,
  so $`u_i\ge-\arctan(2\kappa'\vartheta)`$;
- (b) if $`\operatorname{Re}X_i>0`$ then $`r_i>V^2\ge0`$ and
  $`J_i\ge2\kappa'\vartheta\,r_i\ge2\kappa'\vartheta\operatorname{Re}X_i`$, so
  $`u_i\le\pi-\arctan(2\kappa'\vartheta)`$.

If $`u_{m'}\ge0`$: $`\arg z\le\pi+u_m\le2\pi`$. If $`u_{m'}<0`$: by (a)
$`u_{m'}\ge-\arctan(2\kappa'\vartheta)`$. If $`\operatorname{Re}X_m>0`$, (b)
gives $`u_m\le\pi-\arctan(2\kappa'\vartheta)`$ and $`\arg z\le2\pi`$. If
$`\operatorname{Re}X_m\le0`$ then $`u_m\le\pi/2`$ and
$`\arg z\le\tfrac32\pi+\arctan(2\kappa')<2\pi`$.

*In $`J^-(c_+)`$.* The mirror image under time reversal, with $`m`$ and
$`m'`$ exchanged.

*Elsewhere.* By R4(iii) all $`|u_i|\le e_{\rm sp}`$, so
$`\arg z\le\pi+4e_{\rm sp}<2\pi`$, because $`4e_{\rm sp}<\pi`$.

For $`t>0`$ exchange the two pairs of cusps. $`\blacksquare`$

(R7) holds whenever $`\tau_0=8\Delta`$ and $`R\ge\tau_0`$, as for
`tuned_contour`, and for the default constants.

*What this does not give.* For $`r<2R`$ the frame fields are present. A
frame that does not control an invariant moves it by an amount whose sign
is not fixed, of relative size up to $`\iota/\delta_0`$ (0.007 to 0.012;
$`\iota`$ is defined in section 4), and
$`V^2`$ shifts the real parts. So the method gives $`2\pi+O(\iota/\delta_0)`$
in the inner zone (Proposition R8), not $`2\pi`$. Where $`K_c>0`$ is observed
(generic spacelike points, $`K_c\approx2\kappa'`$) all four $`u_i`$ are within
$`e_{\rm sp}`$ of 0 and $`\arg z\approx\pi`$; where $`\arg z`$ is close to
$`2\pi`$, $`u_{c+2}\approx\pi`$ and $`K_c\le0`$ up to the same small terms. In
the blend zone the constant is set by the sub-case N0nn of Proposition R6.

## 4. The inner zone

**Conditions.** In addition to section 9:

```math
\text{(G1)}\quad4\delta_0\,|u|<D\,|D'| ,\qquad\qquad\text{(G2)}\quad4\delta_0<D .
```

Put $`\gamma_1:=1-2\sqrt{\delta_0|u|/(D|D'|)}\in(0,1]`$.

**Constants.**

```math
\iota:=2\lambda_0(R+\Delta)^2,\qquad v_R:=4\lambda_0(R+\Delta),\qquad e_4:=\arctan\frac{4\iota}{\delta_0-v_R^2},
```
```math
\beta_c:=\frac{2P(R+\Delta)}{D},\qquad G:=\frac{2D}{\delta_0}\max(1,\beta_c),\qquad e_{\rm fr}:=\arctan\frac{\iota}{\delta_0-v_R^2},
```
```math
c_{3d}:=G\Big[2\nu\Big(1+\beta_c\frac{|u|}{|D'|}\Big)+2\kappa\Big(\beta_c\frac{D}{|D'|}+\frac{|u|}{D}+\frac{u^2\beta_c}{|D'|\,D}\Big)\Big].
```

For $`r\le R`$: every frame has $`|V_j|\le\lambda_0(R+\Delta)`$, so
$`|J_i^{(j)}|:=|2V_j\!\cdot n_i|\le\iota`$ for all $`i,j`$, and
$`|V^2|\le|V|^2\le v_R^2`$. In the frame of $`c`$ the coordinates satisfy
$`|\alpha|,|\beta|\le\beta_c`$, because
$`|\beta|=|2y\!\cdot p_{c-1}|/D\le2|y|P/D`$ and $`|y|\le R+\Delta`$; likewise
for $`\alpha`$.

**Lemma R5.1 (two geometric facts).** In the frame of $`c-1`$ write
$`(\alpha',\beta',\rho')`$, in the frame of $`c+1`$ write
$`(\alpha'',\beta'',\rho'')`$. Then

```math
\alpha'=\frac{A}{D'},\qquad\beta'=1+\alpha+\beta\frac{u}{D'}+\frac{2y_\perp\!\cdot p_{c-2}}{D'},\qquad
\beta''=-\frac{B}{D'},\qquad\alpha''=-1+\beta+\alpha\frac{u}{D'}+\frac{2y_\perp\!\cdot p_{c+1}}{D'},
```

and $`|2y_\perp\!\cdot p_{c-2}|,\ |2y_\perp\!\cdot p_{c+1}|\le2\sqrt{\rho\,uD'/D}`$.
Consequently:

- (a) if $`\alpha\ge0`$, $`\beta\ge0`$ and $`r_c\ge-\delta_0`$, then
  $`\beta'\ge\gamma_1>0`$;
- (b) if $`\alpha\le0`$, $`\beta\le0`$ and $`r_c\ge-\delta_0`$, then
  $`\alpha''\le-\gamma_1<0`$.

*Proof.* In the frame of $`j`$: $`r_{j-1}=r_j+\beta_jD_j`$ and
$`r_{j+1}=r_j-\alpha_jD_j`$. For $`j=c-1`$: $`r_c=r_{c-1}-\alpha'D'`$, so
$`\alpha'=A/D'`$; and $`\beta'D'=r_{c-2}-r_{c-1}=2n_{c-1}\!\cdot p_{c-2}`$
because $`n_{c-2}=n_{c-1}+p_{c-2}`$. With $`n_{c-1}=y+p_{c-1}`$,
$`2p_{c-1}\!\cdot p_{c-2}=D'`$, $`2p_c\!\cdot p_{c-2}=u`$:
$`\beta'D'=D'+\alpha D'+\beta u+2y_\perp\!\cdot p_{c-2}`$. For $`j=c+1`$:
$`r_c=r_{c+1}+\beta''D'`$, so $`\beta''=-B/D'`$; and
$`\alpha''D'=r_{c+1}-r_{c+2}=2n_{c+1}\!\cdot p_{c+1}`$ because
$`n_{c+2}=n_{c+1}-p_{c+1}`$. With $`n_{c+1}=y-p_c`$,
$`2p_c\!\cdot p_{c+1}=D'`$, $`2p_{c-1}\!\cdot p_{c+1}=u`$:
$`\alpha''D'=-D'+\alpha u+\beta D'+2y_\perp\!\cdot p_{c+1}`$.

Decompose $`p_{c-2}=\tfrac uDp_{c-1}+\tfrac{D'}Dp_c+q`$ with $`q`$ transverse;
$`0=p_{c-2}^2=\tfrac{uD'}{D}+q^2`$, so $`-q^2=uD'/D>0`$ and
$`|y_\perp\!\cdot q|\le\sqrt{\rho\,uD'/D}`$ (Cauchy–Schwarz in the
transverse plane). The same for $`p_{c+1}=\tfrac{D'}Dp_{c-1}+\tfrac uDp_c+q'`$.

(a) $`\rho=\alpha\beta D-r_c\le\alpha\beta D+\delta_0`$, so
$`2\sqrt{\rho uD'/D}\le2\sqrt{\alpha|D'|\,\beta|u|}+2\sqrt{\delta_0uD'/D}\le\alpha|D'|+\beta|u|+2\sqrt{\delta_0uD'/D}`$.
Since $`D'<0`$ and $`u<0`$,
$`\beta'=1+\big(\alpha|D'|+\beta|u|-2y_\perp\!\cdot p_{c-2}\big)/|D'|\ \ge1-2\sqrt{\delta_0uD'/D}/|D'|=\gamma_1`$.

(b) Now $`\alpha u=|\alpha||u|`$, $`\beta D'=|\beta||D'|`$,
$`\rho\le|\alpha||\beta|D+\delta_0`$, and
$`\alpha''=-1-\big(|\alpha||u|+|\beta||D'|+2y_\perp\!\cdot p_{c+1}\big)/|D'|\le-1+2\sqrt{\delta_0uD'/D}/|D'|=-\gamma_1`$.
$`\blacksquare`$

**Proposition R5.** Assume the conditions of section 9, (G1), (G2), and
$`v_R^2<\delta_0`$, $`\theta\le2`$ (so that $`e_4`$ and $`\varepsilon_V`$ are defined). Then
for $`r\le R`$, off the four edges,

```math
K_c\ \le\ K_{\rm in}:=\max\Big(2e+e_4,\ \ e+2\varepsilon_V,\ \ \arctan\!\Big(2\nu+2c_{3d}+\frac{2\iota}{\delta_0}\Big)+2\arctan\frac{\iota}{\delta_0}+3\varepsilon_V\Big),
\qquad\varepsilon_V:=\max\Big(\arcsin\frac\theta2,\ \arcsin\frac{v_R^2}{\delta_0}\Big).
```

*Proof.* The field is $`V=\sum_jm_jV_j`$.

**Case 1: $`r_c\le-\delta_0`$.** $`\operatorname{Re}X_c=r_c-V^2\le-(\delta_0-v_R^2)<0`$
and $`|J_c|\le4\iota`$, so $`u_c\le e_4`$. With $`u_{c\pm1}\ge-e`$:
$`K_c\le e_4+2e`$.

**Case 2: $`r_c\ge\delta_0`$.** Then $`AB<0`$, so one of $`r_{c\pm1}`$, say
$`r_{\rm nb}`$, is $`\ge r_c\ge\delta_0`$. Then
$`\operatorname{Re}X_{\rm nb}\ge\delta_0-v_R^2`$ and $`|J_{\rm nb}|\le4\iota`$,
so $`u_{\rm nb}\ge\pi-e_4`$. With $`u_c\le\pi+e`$ and the third $`u\ge-e`$:
$`K_c\le2e+e_4`$.

**Case 3: $`|r_c|<\delta_0`$.** The frame $`c+2`$ is off. The frame $`c-1`$
is on only if $`|r_{c+1}|>\delta_0`$; the frame $`c+1`$ only if
$`|r_{c-1}|>\delta_0`$. By Lemma 6 (applied to $`i=c`$),
$`|V^2|\le\tfrac\theta2J_c`$ and $`J_c\ge0`$.

*Linearised numbers.* Put $`\dot X_i:=r_i-iJ_i`$ and
$`\dot u_i:=\pi+\arg\dot X_i`$. If $`|r_i|\le\delta_0`$ then
$`|X_i-\dot X_i|=|V^2|\le\tfrac\theta2J_i\le\tfrac\theta2|\dot X_i|`$ (Lemma
6); if $`|r_i|>\delta_0`$ then $`|X_i-\dot X_i|\le v_R^2\le\tfrac{v_R^2}{\delta_0}|\dot X_i|`$.
A perturbation of relative size $`\epsilon`$ turns a complex number by at
most $`\arcsin\epsilon`$. So

```math
|u_i-\dot u_i|\le\varepsilon_V\qquad(i=c-1,c,c+1).
```

**3a: $`\beta\ge0`$ and $`\alpha\ge0`$.** Then $`A\ge0\ge B`$, so
$`r_{c-1}\ge r_c\ge r_{c+1}`$.

If $`r_{c-1}\ge\delta_0`$: as in Case 2, $`K_c\le2e+e_4`$.

Otherwise $`|r_{c-1}|<\delta_0`$ (because $`r_{c-1}\ge r_c>-\delta_0`$). The
frames $`c+1`$ and $`c+2`$ are off; only the frames $`c`$ and $`c-1`$ can be
on, and both control $`c`$ and $`c-1`$, so $`J_c,J_{c-1}\ge0`$.

- Frame $`c`$ ($`\sigma=+1`$): $`\beta D=A=r_{c-1}-r_c<2\delta_0`$, so
  $`\ell_3=\beta-\alpha-1<2\delta_0/D-1<0`$ by (G2). Hence $`d_\beta=0`$ and
  $`d_\alpha=\phi'(\alpha)`$, and by Lemma R3

```math
r_{c-1}J_c^{(c)}-r_cJ_{c-1}^{(c)}=D\,\rho\,\big(2\nu\beta+\kappa\phi'(\alpha)\big)\ \ge0 .
```

- Frame $`c-1`$ ($`\sigma=-1`$, coordinates $`\alpha',\beta'`$). The pair
  $`(c-1,c)`$ is its "next" pair, so by Lemma R3 with $`j=c-1`$

```math
r_{c-1}J_c^{(c-1)}-r_cJ_{c-1}^{(c-1)}=-\big(r_cJ_{c-1}^{(c-1)}-r_{c-1}J_c^{(c-1)}\big)
=-|D'|\Big[\kappa|D'|\alpha'^2d'_\alpha+\rho'\big(2\nu\alpha'-\kappa d'_\beta\big)\Big].
```

  By Lemma R5.1, $`\alpha'=A/D'\le0`$ and $`\beta'\ge\gamma_1>0`$. So
  $`\phi'(\alpha')=0=\phi'(-\beta')`$, $`d'_\alpha=-\phi'(\ell_3')\le0`$,
  $`d'_\beta=\phi'(\ell_3')\ge0`$, and the right side equals
  $`|D'|\big[\kappa|D'|\alpha'^2\phi'(\ell_3')+\rho'(2\nu|\alpha'|+\kappa\phi'(\ell_3'))\big]\ge0`$.

Adding with the weights, $`r_{c-1}J_c-r_cJ_{c-1}\ge0`$, i.e.
$`\operatorname{Im}(\dot X_{c-1}\overline{\dot X_c})\ge0`$. If
$`J_c=J_{c-1}=0`$ then $`V^2=0`$ by Lemma 6, the $`X_i`$ are the real
invariants, and Lemma R1 gives $`K_c\le0`$. Otherwise Lemma R2(4) gives
$`\dot u_c\le\dot u_{c-1}`$, so $`u_c\le u_{c-1}+2\varepsilon_V`$, and with
$`u_{c+1}\ge-e`$: $`K_c\le e+2\varepsilon_V`$.

**3b: $`\beta\le0`$ and $`\alpha\le0`$.** Then $`B\ge0\ge A`$, so
$`r_{c+1}\ge r_c\ge r_{c-1}`$.

If $`r_{c+1}\ge\delta_0`$: $`K_c\le2e+e_4`$ as in Case 2.

Otherwise $`|r_{c+1}|<\delta_0`$. The frames $`c-1`$ and $`c+2`$ are off;
only $`c`$ and $`c+1`$ can be on, and both control $`c`$ and $`c+1`$.

- Frame $`c`$: $`|\alpha|D=B<2\delta_0`$, so
  $`\ell_3=\beta-\alpha-1\le2\delta_0/D-1<0`$ by (G2). With $`\alpha\le0`$:
  $`d_\alpha=0`$, $`d_\beta=-\phi'(-\beta)`$. Lemma R3 gives

```math
r_{c+1}J_c^{(c)}-r_cJ_{c+1}^{(c)}=-D\,\rho\,\big(2\nu\alpha+\kappa d_\beta\big)=D\,\rho\,\big(2\nu|\alpha|+\kappa\phi'(|\beta|)\big)\ \ge0 .
```

- Frame $`c+1`$ ($`\sigma=-1`$, coordinates $`\alpha'',\beta''`$). The pair
  $`(c+1,c)`$ is its "previous" pair, so by Lemma R3 with $`j=c+1`$

```math
r_{c+1}J_c^{(c+1)}-r_cJ_{c+1}^{(c+1)}=-\big(r_cJ_{c+1}^{(c+1)}-r_{c+1}J_c^{(c+1)}\big)
=|D'|\Big[\kappa|D'|\beta''^2d''_\beta+\rho''\big(2\nu\beta''-\kappa d''_\alpha\big)\Big].
```

  By Lemma R5.1, $`\beta''=-B/D'=B/|D'|\ge0`$ and $`\alpha''\le-\gamma_1<0`$.
  So $`\phi'(-\beta'')=0=\phi'(\alpha'')`$, $`d''_\beta=\phi'(\ell_3'')\ge0`$,
  $`d''_\alpha=-\phi'(\ell_3'')\le0`$, and the right side is
  $`|D'|\big[\kappa|D'|\beta''^2\phi'(\ell_3'')+\rho''(2\nu\beta''+\kappa\phi'(\ell_3''))\big]\ge0`$.

As in 3a: $`\dot u_c\le\dot u_{c+1}`$, and $`K_c\le e+2\varepsilon_V`$.

**3c: $`\beta\ge0`$ and $`\alpha\le0`$.** Then $`A,B\ge0`$, so
$`r_{c\pm1}\ge r_c>-\delta_0`$. If one of them is $`\ge\delta_0`$:
$`K_c\le2e+e_4`$ as in Case 2. Otherwise all three of
$`|r_{c-1}|,|r_c|,|r_{c+1}|`$ are below $`\delta_0`$. Then the frames
$`c-1,c+1,c+2`$ are off, and by (H1) $`|r_{c+2}|>2\delta_0`$, so
$`V=V_c`$. Also $`\beta D,|\alpha|D<2\delta_0`$, so
$`\ell_3=\beta+|\alpha|-1<4\delta_0/D-1<0`$ by (G2). Hence $`d=0`$ and
$`V=\nu y_\perp`$, $`J_{c-1}=J_c=J_{c+1}=2\nu\rho`$, $`V^2=-\nu^2\rho`$.
If $`\rho=0`$ the point is real and $`K_c\le0`$ by Lemma R1. If
$`\rho>0`$ the three $`X_i`$ have the same negative imaginary part and
$`\operatorname{Re}X_{c\pm1}\ge\operatorname{Re}X_c`$, so
$`u_c\le u_{c\pm1}`$, and $`u_{c\mp1}\ge0`$: $`K_c\le0`$.

**3d: $`\beta<0`$ and $`\alpha>0`$** (the quadrant beyond the cusp). Put
$`a=\alpha`$, $`b=|\beta|`$. Then $`A=-bD<0`$, $`B=-aD<0`$, and all three
invariants are negative:

```math
|r_c|=abD+\rho,\qquad|r_{c-1}|=|r_c|+bD,\qquad|r_{c+1}|=|r_c|+aD,\qquad\rho\le|r_c|<\delta_0 .
```

The frames $`c-1,c,c+1`$ can be on; $`J_c\ge0`$. We bound
$`J_c^{(j)}/|r_c|`$ frame by frame.

*Frame $`c`$.* $`\ell_3=-b-a-1<0`$, so $`d_\beta=-\phi'(b)`$,
$`d_\alpha=\phi'(a)`$, and

```math
J_c^{(c)}=\kappa D\big(a\phi'(a)+b\phi'(b)\big)+2\nu\rho,\qquad J_{c-1}^{(c)}=J_c^{(c)}+\kappa D\phi'(a),\qquad J_{c+1}^{(c)}=J_c^{(c)}+\kappa D\phi'(b).
```

Since $`a\big((1+a)bD+\rho\big)\le(1+a)\big(abD+\rho\big)`$,

```math
\frac{\kappa D\,a\phi'(a)}{|r_c|}\le\frac{\kappa D(1+a)\phi'(a)}{|r_{c-1}|}\le\frac{J_{c-1}^{(c)}}{|r_{c-1}|},\qquad
\frac{\kappa D\,b\phi'(b)}{|r_c|}\le\frac{\kappa D(1+b)\phi'(b)}{|r_{c+1}|}\le\frac{J_{c+1}^{(c)}}{|r_{c+1}|},\qquad\frac{2\nu\rho}{|r_c|}\le2\nu .
```

Hence

```math
\frac{J_c^{(c)}}{|r_c|}\le\frac{J_{c-1}^{(c)}}{|r_{c-1}|}+\frac{J_{c+1}^{(c)}}{|r_{c+1}|}+2\nu. \tag{3d.1}
```

*Frame $`c-1`$* (on only if $`|r_{c+1}|>\delta_0`$). Here
$`\alpha'=A/D'=bD/|D'|>0`$, $`\alpha'|D'|=|A|=bD\le|r_{c-1}|`$, and
$`\rho'=\alpha'\beta'D'-r_{c-1}=|r_{c-1}|-\alpha'\beta'|D'|`$. All invariants
being negative, Lemma R3 gives

```math
J_c^{(c-1)}|r_{c-1}|-J_{c-1}^{(c-1)}|r_c|=|D'|\Big[\kappa|D'|\alpha'^2d'_\alpha+\rho'\big(2\nu\alpha'-\kappa d'_\beta\big)\Big]=:E' ,
```

with $`d'_\alpha\le\phi'(\alpha')\le2\alpha'`$ and
$`-d'_\beta\le\phi'(\beta'_-)\le2\beta'_-`$, $`\beta'_-:=\max(-\beta',0)`$.

(i) $`bD/|r_c|\le G`$. If $`|r_c|\le\delta_0/2`$ then
$`aD=|r_{c+1}|-|r_c|>\delta_0/2`$ and $`bD/|r_c|\le bD/(abD)=1/a<2D/\delta_0`$;
if $`|r_c|>\delta_0/2`$ then $`bD/|r_c|<2bD/\delta_0\le2\beta_cD/\delta_0`$.

(ii) $`\beta'_-\le b\,|u|/|D'|`$. By Lemma R5.1,
$`\beta'\ge1+a-b\tfrac{|u|}{|D'|}-2\sqrt{\rho|u|/(D|D'|)}`$, and
$`2\sqrt{\rho|u|/(D|D'|)}\le1-\gamma_1\le1`$ because $`\rho<\delta_0`$. So
$`-\beta'\le b\tfrac{|u|}{|D'|}-a-\gamma_1\le b\tfrac{|u|}{|D'|}`$.

Now $`\rho'\le|r_{c-1}|+\alpha'\beta'_-|D'|`$, so

```math
E'\le|D'|\Big[2\kappa|D'|\alpha'^3+\big(|r_{c-1}|+\alpha'\beta'_-|D'|\big)\big(2\nu\alpha'+2\kappa\beta'_-\big)\Big].
```

Divide by $`|r_c||r_{c-1}|`$ and use $`\alpha'|D'|=bD\le|r_{c-1}|`$,
$`b\le\beta_c`$, (i), (ii), term by term:

```math
\frac{2\kappa|D'|^2\alpha'^3}{|r_c||r_{c-1}|}\le2\kappa\alpha'\frac{bD}{|r_c|}\le2\kappa\,\beta_c\frac{D}{|D'|}\,G,\qquad
\frac{|r_{c-1}|\,2\nu\alpha'|D'|}{|r_c||r_{c-1}|}=2\nu\frac{bD}{|r_c|}\le2\nu G,
```
```math
\frac{|r_{c-1}|\,2\kappa\beta'_-|D'|}{|r_c||r_{c-1}|}\le\frac{2\kappa\,b|u|}{|r_c|}\le2\kappa\frac{|u|}{D}\,G,\qquad
\frac{\alpha'\beta'_-|D'|^2\,2\nu\alpha'}{|r_c||r_{c-1}|}=\frac{2\nu(bD)^2\beta'_-}{|r_c||r_{c-1}|}\le2\nu\beta'_-\frac{bD}{|r_c|}\le2\nu\,\beta_c\frac{|u|}{|D'|}\,G,
```
```math
\frac{\alpha'\beta'_-|D'|^2\,2\kappa\beta'_-}{|r_c||r_{c-1}|}=\frac{2\kappa\,bD\,|D'|\,\beta_-'^2}{|r_c||r_{c-1}|}\le\frac{2\kappa|D'|\beta_-'^2}{|r_c|}\le\frac{2\kappa\,b^2u^2}{|D'|\,|r_c|}\le2\kappa\frac{u^2\beta_c}{|D'|D}\,G .
```

The sum is $`c_{3d}`$. Hence

```math
\frac{J_c^{(c-1)}}{|r_c|}\le\frac{J_{c-1}^{(c-1)}}{|r_{c-1}|}+c_{3d}. \tag{3d.2}
```

*Frame $`c+1`$* (on only if $`|r_{c-1}|>\delta_0`$). Here
$`\beta''=-B/D'=-aD/|D'|<0`$, $`|\beta''||D'|=|B|=aD\le|r_{c+1}|`$, and
$`\rho''=\alpha''\beta''D'-r_{c+1}=|r_{c+1}|+\alpha''|\beta''||D'|`$. Lemma R3
gives

```math
J_c^{(c+1)}|r_{c+1}|-J_{c+1}^{(c+1)}|r_c|=|D'|\Big[-\kappa|D'|\beta''^2d''_\beta+\rho''\big(2\nu|\beta''|+\kappa d''_\alpha\big)\Big]=:E'' ,
```

with $`-d''_\beta\le\phi'(|\beta''|)\le2|\beta''|`$ and
$`d''_\alpha\le\phi'(\alpha''_+)\le2\alpha''_+`$, $`\alpha''_+:=\max(\alpha'',0)`$.

(i) $`aD/|r_c|\le G`$: if $`|r_c|\le\delta_0/2`$ then
$`bD=|r_{c-1}|-|r_c|>\delta_0/2`$ and $`aD/|r_c|\le1/b<2D/\delta_0`$; otherwise
$`aD/|r_c|<2aD/\delta_0\le2\beta_cD/\delta_0`$.

(ii) $`\alpha''_+\le a\,|u|/|D'|`$: by Lemma R5.1,
$`\alpha''\le-1-b+a\tfrac{|u|}{|D'|}+2\sqrt{\rho|u|/(D|D'|)}\le a\tfrac{|u|}{|D'|}-b-\gamma_1`$.

Then $`\rho''\le|r_{c+1}|+\alpha''_+|\beta''||D'|`$ and

```math
E''\le|D'|\Big[2\kappa|D'||\beta''|^3+\big(|r_{c+1}|+\alpha''_+|\beta''||D'|\big)\big(2\nu|\beta''|+2\kappa\alpha''_+\big)\Big],
```

and the same five estimates, with $`(\alpha',\beta'_-,b,r_{c-1})`$ replaced by
$`(|\beta''|,\alpha''_+,a,r_{c+1})`$, give

```math
\frac{J_c^{(c+1)}}{|r_c|}\le\frac{J_{c+1}^{(c+1)}}{|r_{c+1}|}+c_{3d}. \tag{3d.3}
```

(G1) enters these two paragraphs only through $`\gamma_1\ge0`$.

*Sum.* $`J_i=\sum_{j\in\{c-1,c,c+1\}}m_jJ_i^{(j)}`$. Multiply (3d.1),
(3d.2), (3d.3) by $`m_c,m_{c-1},m_{c+1}\in[0,1]`$ and add:

```math
\frac{J_c}{|r_c|}\ \le\ \frac{J_{c-1}-m_{c+1}J_{c-1}^{(c+1)}}{|r_{c-1}|}+\frac{J_{c+1}-m_{c-1}J_{c+1}^{(c-1)}}{|r_{c+1}|}+2\nu+2c_{3d} .
```

The two subtracted terms are the contributions of the frames that do not
control that invariant. They are present only where the invariant exceeds
$`\delta_0`$ in modulus, and are at most $`\iota`$ in modulus. Therefore

```math
\frac{J_c}{|r_c|}\ \le\ T_{c-1}+T_{c+1}+2\nu+2c_{3d}+\frac{2\iota}{\delta_0},\qquad T_i:=\frac{J_i}{|r_i|}\ \ge-\frac{\iota}{\delta_0}\quad(i=c\pm1). \tag{3d.4}
```

(The lower bound on $`T_i`$: the controlling frames contribute $`\ge0`$ to
$`J_i`$.) Now $`\dot u_i=\arctan T_i`$ for all three indices, because the
real parts are negative. By (3d.4), monotonicity and subadditivity of
$`\arctan`$ on $`[0,\infty)`$,

```math
\dot u_c\le\arctan(T_{c-1}^+)+\arctan(T_{c+1}^+)+\arctan\!\Big(2\nu+2c_{3d}+\frac{2\iota}{\delta_0}\Big),\qquad
\arctan(T_i^+)\le\dot u_i+\arctan\frac{\iota}{\delta_0} .
```

With $`|u_i-\dot u_i|\le\varepsilon_V`$ this is the third entry of
$`K_{\rm in}`$. All cases are covered, since the four sign combinations of
$`(\alpha,\beta)`$ are 3a–3d. $`\blacksquare`$

**Proposition R8 (the inner zone with the frame part of the excursion).**
Assume the hypotheses of Proposition R5. Then for $`r\le R`$, off the four
edges,

```math
u_i\in[-e_{\rm fr},\pi+e_{\rm fr}]\quad(i=1,\dots,4),\qquad
K_c\le K_{\rm in}^{\rm fr}:=\max\Big(2e_{\rm fr}+e_4,\ e_{\rm fr}+2\varepsilon_V,\ \arctan\!\Big(2\nu+2c_{3d}+\frac{2\iota}{\delta_0}\Big)+2\arctan\frac\iota{\delta_0}+3\varepsilon_V\Big),
```

and therefore $`|\arg z|\le2\pi+b_{\rm in}`$ with
$`b_{\rm in}:=e_{\rm fr}+K_{\rm in}^{\rm fr}`$.

*Proof.* For $`r\le R`$, $`J_i=\sum_jm_jJ_i^{(j)}`$. The frames
$`j\ne i+2`$ control $`i`$, so their terms are $`\ge0`$. If $`J_i\ge0`$ then
$`\operatorname{Im}X_i\le0`$ and $`u_i\in[0,\pi]`$. If $`J_i<0`$ then
$`m_{i+2}>0`$, so $`|r_i|>\delta_0`$, and
$`0<\operatorname{Im}X_i\le m_{i+2}|J_i^{(i+2)}|\le\iota`$, while
$`|\operatorname{Re}X_i|\ge\delta_0-v_R^2`$. On the Feynman branch this
gives $`u_i\in[-e_{\rm fr},0)`$ if $`\operatorname{Re}X_i<0`$ and
$`u_i\in(\pi,\pi+e_{\rm fr}]`$ if $`\operatorname{Re}X_i>0`$.

The proof of Proposition R5 uses (D) only through the inequalities
$`u_i\ge-e`$ and $`u_i\le\pi+e`$ at the point considered: in Case 1
($`u_{c\pm1}\ge-e`$), in Case 2 and its repetitions in 3a, 3b, 3c
($`u_c\le\pi+e`$ and one $`u\ge-e`$), and at the end of 3a and 3b
($`u_{c\pm1}\ge-e`$). Case 3d does not use it. Repeating the proof with
$`e_{\rm fr}`$ gives $`K_c\le K^{\rm fr}_{\rm in}`$. Finally
$`\arg z=\pi+K_c+u_{c+2}\le2\pi+e_{\rm fr}+K_c`$ as in Lemma R1.
$`\blacksquare`$

$`b_{\rm in}`$ is 0.046 for the default constants and between 0.084 and
0.097 at the tuned squares of the table in section 6.

**Proposition R5′ (sharp constant next to the cusp).** Assume (H1), (G2),
$`\nu<1`$ and $`4\kappa^2<1`$.
At points with $`r\le R`$ and $`|r_{c-1}|,|r_c|,|r_{c+1}|\le\delta_0`$,

```math
K_c\le\frac{2\nu}{1-\nu^2}+\frac{4\kappa^2}{1-4\kappa^2}.
```

*Proof.* As in 3c the field is exactly $`V_c`$, and
$`|\alpha|,|\beta|\le2\delta_0/D`$ gives $`\ell_3<0`$. All three $`J_i\ge0`$.
By Lemma R3, $`V^2=\kappa^2Dd_\beta d_\alpha-\nu^2\rho`$.

*$`\alpha,\beta\ge0`$:* $`d_\beta=0`$, $`d_\alpha=\phi'(\alpha)`$,
$`V^2=-\nu^2\rho`$, $`J_{c-1}-J_c=\kappa D\phi'(\alpha)`$, so by Lemma R2(3)
and R3,
$`\operatorname{Im}(X_{c-1}\overline{X_c})=D\rho\big(2\nu\beta+(1-\nu^2)\kappa\phi'(\alpha)\big)\ge0`$.
Then $`u_c\le u_{c-1}`$ (Lemma R2(4); if $`J_c=J_{c-1}=0`$ then
$`\rho=0`$ and $`\alpha=0`$, so $`V=0`$ and the point is real) and
$`u_{c+1}\ge0`$: $`K_c\le0`$.

*$`\alpha,\beta\le0`$:* $`d_\alpha=0`$, $`d_\beta=-\phi'(|\beta|)`$,
$`V^2=-\nu^2\rho`$, $`J_{c+1}-J_c=\kappa D\phi'(|\beta|)`$, and
$`\operatorname{Im}(X_{c+1}\overline{X_c})=D\rho\big(2\nu|\alpha|+(1-\nu^2)\kappa\phi'(|\beta|)\big)\ge0`$:
$`K_c\le0`$.

*$`\alpha\le0\le\beta`$:* case 3c, $`K_c\le0`$.

*$`\alpha=a>0>\beta=-b`$:* $`V^2=-\kappa^2D\phi'(a)\phi'(b)-\nu^2\rho`$, so
$`q_i:=-\operatorname{Re}X_i`$ are

```math
q_c=D\big(ab-\kappa^2\phi'(a)\phi'(b)\big)+(1-\nu^2)\rho\ \ge\ (1-4\kappa^2)abD+(1-\nu^2)\rho>0,\qquad q_{c-1}=q_c+bD,\quad q_{c+1}=q_c+aD .
```

With $`\gamma:=1/(1-4\kappa^2)`$: $`a\,q_{c-1}=aq_c+abD\le aq_c+\gamma q_c\le\gamma(1+a)q_c`$, so
$`\kappa Da\phi'(a)/q_c\le\gamma\kappa D(1+a)\phi'(a)/q_{c-1}\le\gamma J_{c-1}/q_{c-1}`$;
likewise with $`b`$ and $`c+1`$; and $`2\nu\rho/q_c\le2\nu/(1-\nu^2)`$. So
$`\tan u_c\le\gamma\tan u_{c-1}+\gamma\tan u_{c+1}+\tfrac{2\nu}{1-\nu^2}`$ with
all $`u_i\in[0,\tfrac\pi2)`$. Since
$`\arctan(\gamma x)\le\arctan x+\tfrac{\gamma-1}{2}`$ for $`x\ge0`$
($`\partial_\gamma\arctan(\gamma x)=x/(1+\gamma^2x^2)\le1/(2\gamma)`$), the
claim follows. $`\blacksquare`$

**Remark.** $`2\nu`$ is the supremum of $`K_c`$ in the region of Proposition
R5′, and it is attained there to five digits. It is not the maximum over the
inner zone: at a point of Case 1 where the frame $`c+2`$ is just switching
on, an independent check found $`K_c=2.25\times10^{-6}=7.5\times2\nu`$ (default
constants). The proved constant $`K_{\rm in}`$ is five orders of magnitude
above either.

## 5. The blend zone

For $`R<r<2R`$ the field is $`V=(1-g)V_{\rm in}+gV_{\rm out}`$ with
$`g\in(0,1)`$, $`V_{\rm in}=\sum_jm_jV_j`$ and
$`V_{\rm out}=-\kappa'\Theta_+\xi_--\kappa'\Theta_-\xi_+`$. The frame fields
are not switched off here: their weights $`m_j`$ depend only on
$`|r_{j+2}|/\delta_0`$. Near a cone they are not small compared with the far
field either; their share of $`\operatorname{Im}X_i`$ can exceed the far
field's for $`g`$ below about 0.5 with the default constants and below about
0.95 with the tuned ones. So the pair criterion of Lemma R2 is checked field
by field.

**Further constants.** With $`\mu=\max(64\kappa P^2/\min_j|D_j|,\nu B^2)`$ and
$`M`$ the margin of (C1), as in Lemma Q; with $`B_c`$ the transverse norm of
the frames of the two cusps with the timelike diagonal and $`B_n`$ that of
the other two frames ($`B=\max(B_c,B_n)`$):

```math
K:=\kappa'(2R+\tau_0+\Delta),\quad v_{\rm in}:=4\lambda_0(2R+\Delta),\quad v:=v_{\rm in}+K,\quad
\iota_2:=2\lambda_0(2R+\Delta)^2,\quad\pi_0:=2\lambda_0(2R+\Delta)P,
```
```math
\mu_-:=3\mu+K\sqrt{\frac{3\mu}{2\kappa'M}},\qquad\mu_+:=\mu_-+\frac{K^2}{2\kappa'M},\qquad\omega_0:=v_{\rm in}(2K+v_{\rm in}),\qquad\delta_v:=\delta_0-v^2,
```
```math
\varepsilon:=\frac{\sqrt2\,P}{R-\Delta},\qquad b_0:=\frac{\varepsilon\,\delta_0}{\min_j|D_j|},\qquad\delta_1:=\delta_0(1-\varepsilon),
```
```math
\Lambda_c:=\frac{R-\Delta-B_c\sqrt{\delta_0+b_0^2D}}{P},\qquad
\Lambda_n:=\frac{R-\Delta-B_n\sqrt{\delta_0(1+\varepsilon)+2b_0^2|D'|}}{P},\qquad a_{\min}:=\Lambda_c-b_0,
```
```math
\varsigma:=\Big(2\nu+\frac{\kappa b_0^2}{a_{\min}(1+a_{\min})}\Big)+\frac{2\nu+\kappa b_0^2}{a_{\min}},\qquad F:=\frac{2KP\mu_-}{\delta_v}.
```

```math
u_2:=\frac{\iota_2}{\delta_v},\qquad h:=0.52,\qquad
\mathfrak a:=\frac{2P(2R+\Delta)}{D},\qquad\mathfrak a':=\frac{2P(2R+\Delta)}{|D'|},
```
```math
\varsigma_c:=2\nu+\frac{2D}{\delta_0}\Big(\kappa\mathfrak a^2+\nu\mathfrak a^2+\kappa\mathfrak a\max(\mathfrak a+1,2)\Big),\qquad
\varsigma_n:=\frac{2|D'|}{\delta_0}\Big(\kappa\mathfrak a'^2+\nu(\mathfrak a'^2+\mathfrak a')+\kappa\mathfrak a'(\mathfrak a'+1)\Big),
```
```math
\tau_1:=\frac{2\Delta\tau_0}{\tau_0-2\Delta},\qquad W_1:=\Delta(\tau_0+\Delta)+\Delta(\tau_1+\tau_0+\Delta),\qquad
\rho_1:=\Big(\sqrt{R^2-\tau_1^2}-\Delta\Big)^2-(\tau_1+\Delta)^2,
```
```math
\tan e^+_{\rm bl}:=h\,u_2,\qquad\tan e^-_{\rm bl}:=h\,u_2+\frac{2\kappa'(1+W_1/\rho_1)}{1-v^2/\delta_0},\qquad\chi(t)=t^2(3-2t)\ \ (0\le t\le1).
```

For $`r<2R`$: $`|V_j|\le\lambda_0(2R+\Delta)`$, so
$`|J^{(j)}_i|\le\iota_2`$, $`|2V_j\!\cdot p|\le\pi_0`$ for an edge $`p`$,
$`|V_{\rm in}|\le v_{\rm in}`$, $`|V_{\rm out}|\le K`$, $`|V^2|\le v^2`$. In the frame of $`c`$
the coordinates satisfy $`|\alpha|,|\beta|\le\mathfrak a`$, in the frame of
$`c\mp1`$ they satisfy $`|\alpha'|,|\beta'|\le\mathfrak a'`$ (as in section 4,
with $`|n_j|\le2R+\Delta`$).

**Conditions (B).** Besides section 9 of GATE_6_CONTOUR ((C0), (C1), (C2),
(D), in particular $`v^2<\delta_0/2`$, which also gives
$`\omega_0<\delta_0`$):

- (B1) $`R\ge\tau_0`$ and $`R>\sqrt5\,\Delta`$;
- (B2) $`\kappa\le\nu`$;
- (B3) $`b_0\le\tfrac12`$, $`\Lambda_c\ge2`$, $`\Lambda_n\ge2`$, $`4\delta_0\le\min_j|D_j|`$;
- (B4) $`\delta_1>v^2`$ (in particular $`\varepsilon<1`$);
- (B5) $`s<1`$ and $`e+\arcsin s+\arctan\mu_-<\tfrac\pi2`$, with $`s`$ defined
  in the statement;
- (B6) $`\tau_1<R`$, $`\sqrt{R^2-\tau_1^2}-\Delta>\tau_1+\Delta`$ (so that
  $`\rho_1>0`$);
- and $`v^2<(\tau_0-2\Delta)^2`$, $`4e'_{\rm sp}<\pi`$.

**Proposition R6.** Assume the conditions of section 9 and (B). Then for
$`R<r<2R`$

```math
K_c\le K_{\rm bl}:=\max\big(e^-_{\rm bl}+\arcsin s,\ 3e_{\rm sp}'\big),\qquad
|\arg z|\le2\pi+b_{\rm bl},\quad b_{\rm bl}:=e^+_{\rm bl}+e^-_{\rm bl}+\arcsin s,
```
```math
\tan e'_{\rm sp}=\frac{\dfrac{4\iota_2}{(\tau_0-2\Delta)^2}+2\kappa'\Big(1+\sqrt2\dfrac{\tau_0+\Delta}{\tau_0-2\Delta}\Big)}{1-v^2/(\tau_0-2\Delta)^2},
```

where $`s=\max(s_0,s_1,s_2,s_3)`$ and, with
$`T_V:=2v^2\pi_0/\delta_v^2`$ and $`T_F:=2KP\,\omega_0/\delta_v^2`$,

```math
s_0=\max\Big(3hu_2+T_V+T_F,\ \ 2hu_2+T_F,\ \ (1+2h)u_2+(\varsigma_c+\varsigma_n)\frac{\delta_0^2}{\delta_v^2}+T_V+T_F\Big),
```
```math
s_1=\max\Big(hu_2+\frac{2\mu_+\pi_0}{\delta_v}+F,\ \ (1+\mu_+)\frac{\delta_0}{\delta_v}\,\varsigma+\frac{2\mu_-\pi_0}{\delta_v}+\chi(\varepsilon)\,\iota_2\Big(\frac1{\delta_v}+\frac1{\delta_1-v^2}\Big)+F\Big),
```
```math
s_2=hu_2+\frac{2\mu_+\pi_0}{\delta_v}+F,\qquad s_3=\mu_-+(1+\mu_+)^2\varsigma .
```

*Proof.* Write $`x_0-c=(\tau,\vec\xi)`$.

**Step 0: three regions.** (a) If $`x_0\notin J^+(c_-)\cup J^-(c_+)`$ then,
as in R4(iii), $`-r_i\ge(\tau_0-2\Delta)|\vec n_i|\ge(\tau_0-2\Delta)^2`$ and
$`|n_i|\le\sqrt2|\vec n_i|`$ for all $`i`$. The frames contribute at most
$`4\iota_2`$ to $`|\operatorname{Im}X_i|`$, the far field at most
$`2\kappa'|r_i|(1+\sqrt2\tfrac{\tau_0+\Delta}{\tau_0-2\Delta})`$, and
$`|V^2|\le v^2`$. So $`|u_i|\le e'_{\rm sp}`$, $`K_c\le3e'_{\rm sp}`$ and
$`|\arg z|\le\pi+4e'_{\rm sp}<2\pi`$.

(b) If $`x_0\in J^+(c_-)`$: $`\tau+\tau_0\ge|\vec\xi|`$. If $`\tau\le0`$ then
$`|\tau|+|\vec\xi|\le\tau_0`$ and $`r\le\tau_0\le R`$. So $`\tau>0`$,
$`\Theta_-=0`$, and

```math
V=(1-g)V_{\rm in}+g\,\vartheta\,W,\qquad W=-\kappa'\xi_-,\quad\vartheta=\Theta_+\in[0,1],\quad\xi_-\ \text{future causal},\quad W^2\ge0 .
```

(c) If $`x_0\in J^-(c_+)`$: the mirror image under time reversal.

**Symmetries.** The construction commutes with time reversal
$`x^0\mapsto-x^0`$ (which exchanges $`c_\pm`$ and $`\Theta_\pm`$) and with
the relabelling $`x_{c+k}\mapsto x_{c-k}`$, under which
$`p_{c-1}\leftrightarrow-p_c`$, $`(\alpha,\beta)\mapsto(-\beta,-\alpha)`$,
the triangle $`T`$ and $`\Psi`$ are invariant, the neighbours $`c\mp1`$ are
exchanged and the orientation of the edges at $`c`$ is reversed. So it is
enough to treat: $`x_0\in J^+(c_-)`$, the edges $`p_{c-1},p_c`$ **future
directed**, and the neighbour $`c-1`$. In the other three combinations the
neighbour is the one for which $`N_\mp`$ is future null in $`J^+(c_-)`$ and
past null in $`J^-(c_+)`$; the numerical checks cover all four. From here
on the standard combination is assumed. Then $`p_{c-2}`$ is past directed,
because $`2p_{c-2}\!\cdot p_{c-1}=D'<0`$.

**Step 1: facts in $`J^+(c_-)`$, $`R<r<2R`$.**

(F1) If $`r_i>0`$ then $`n_i`$ is future timelike. Otherwise
$`x_0\in J^-(x_i)`$, so with $`x_i-c=(t_i,\vec x_i)`$:
$`t_i-\tau\ge|\vec x_i-\vec\xi|\ge|\vec\xi|-|\vec x_i|`$, hence
$`\tau+|\vec\xi|\le t_i+|\vec x_i|\le\sqrt2\Delta`$ and, since $`\tau>0`$,
$`r\le\sqrt2\Delta<R`$ by (B1).

(F2) Let $`n=n_i`$ with $`n^0>0`$, $`|n|\ge R-\Delta`$ and
$`|n^2|\le\delta`$. Then
$`|\vec n|-n^0\le\delta/(|\vec n|+n^0)\le\delta/(R-\Delta)`$, so for a future
null edge $`p`$ (with $`p^0=|p|/\sqrt2`$):
$`n\!\cdot p\ge-p^0(|\vec n|-n^0)\ge-\tfrac{P\delta}{\sqrt2(R-\Delta)}`$, and
for a past null edge $`n\!\cdot p\le\tfrac{P\delta}{\sqrt2(R-\Delta)}`$. If
$`|r_i|\le\delta_0`$ ("$`i`$ is near"), then $`\Theta_+=1`$ by (C2),
$`\tau\ge|\vec\xi|/2`$, so $`n_i^0\ge R/\sqrt5-\Delta>0`$ by (B1), and the
estimate holds with $`\delta=\delta_0`$: a coordinate
$`2n_i\!\cdot p/D_j`$ along a future edge is $`\ge-b_0`$ if $`D_j>0`$, and a
coordinate along a past edge with $`D_j<0`$ is $`\ge-b_0`$ too. Moreover
$`J^W_i:=-2W\!\cdot n_i\ge2\kappa'M`$ by (C1).

(F3) *Norm.* If $`n=a_1p+a_2p'+y_\perp`$ with $`p,p'`$ the two edges of a
frame and $`y_\perp`$ transverse to it, $`\rho=-y_\perp^2`$, then
$`|n|\le(|a_1|+|a_2|)P+B_\bullet\sqrt\rho`$, where $`B_\bullet`$ is the
transverse norm of that frame: $`B_c`$ for the frame of $`c`$, $`B_n`$ for
the frame of $`c-1`$. For $`n=n_i`$: $`|n_i|\ge R-\Delta`$.

(F4) *Lemma V.* If $`|r_i|\le\delta_0`$ then, with
$`V^2_-:=\max(-V^2,0)`$,

```math
V^2_-\le\mu_-J_i,\qquad|V^2|\le\mu_+J_i .
```

*Proof of (F4).* The frame $`i+2`$ is off, the others control $`i`$, so
$`J^{\rm in}_i:=-2V_{\rm in}\!\cdot n_i\ge0`$ and
$`|V_{\rm in}|^2\le3\mu J^{\rm in}_i`$ (proof of Lemma Q). By (F2),
$`J_i=(1-g)J_i^{\rm in}+gJ_i^W\ge(1-g)J_i^{\rm in}+2\kappa'Mg`$. Now
$`V^2=(1-g)^2V_{\rm in}^2+2g(1-g)V_{\rm in}\!\cdot W+g^2W^2`$ with
$`0\le W^2\le K^2`$. First,
$`(1-g)^2|V_{\rm in}^2|\le3\mu(1-g)J_i^{\rm in}\le3\mu J_i`$. Second,
$`2g(1-g)|V_{\rm in}||W|\le K\sqrt{\tfrac{3\mu}{2\kappa'M}}\cdot2\sqrt{2\kappa'Mg\cdot g(1-g)^2J_i^{\rm in}}\le K\sqrt{\tfrac{3\mu}{2\kappa'M}}\,J_i`$.
Third, $`g^2W^2\le\tfrac{K^2}{2\kappa'M}\cdot2\kappa'Mg\le\tfrac{K^2}{2\kappa'M}J_i`$.
The first two bound $`V^2_-`$, all three bound $`|V^2|`$. $`\square`$

Without the hypothesis $`|r_i|\le\delta_0`$:
$`V^2_-\le(1-g)\big(2g|V_{\rm in}||W|+(1-g)|V_{\rm in}|^2\big)\le\omega_0`$.

(F5) *The frame $`i+2`$ against $`\operatorname{Re}X_i`$.* For every $`i`$,

```math
m_{i+2}\,\big|J_i^{(i+2)}\big|\ \le\ h\,u_2\,|R_i| .
```

*Proof.* If $`m_{i+2}>0`$ then $`|r_i|=\delta_0(1+t)`$ with $`t>0`$ and
$`m_{i+2}=\chi(t)`$ ($`=1`$ for $`t\ge1`$). The polynomial
$`2t^3-3t^2+0.52t+0.52`$ has its minimum on $`[0,1]`$ at
$`t=\tfrac12+\tfrac1{12}\sqrt{36-12.48}=0.9041\ldots`$, where it is larger than
$`0.015`$; so $`\chi(t)\le h(1+t)`$ for all $`t\ge0`$. Since
$`|R_i|\ge\delta_0(1+t)-v^2\ge(1+t)(\delta_0-v^2)`$ and
$`|J_i^{(i+2)}|\le\iota_2`$, the claim follows. $`\square`$

(F6) *Excursions in $`J^+(c_-)`$.* For every $`i`$:
$`-e^-_{\rm bl}\le u_i\le\pi+e^+_{\rm bl}`$.

*Proof.* $`J_i=(1-g)\sum_jm_jJ_i^{(j)}+g\vartheta J_i^W`$, with
$`J_i^W=2\kappa'(r_i+n_i\!\cdot w_i)`$, $`w_i=x_i-c_-`$. If $`J_i\ge0`$ then
$`u_i\in[0,\pi]`$. Let $`J_i<0`$. Then $`|r_i|>\delta_0`$: if $`i`$ is near,
the frame $`i+2`$ is off, the other frames contribute $`\ge0`$ and
$`J_i^W>0`$ by (F2). So $`|R_i|\ge|r_i|(1-v^2/\delta_0)`$, and $`R_i`$ has the
sign of $`r_i`$. The frames $`j\ne i+2`$ contribute $`\ge0`$, and the frame
$`i+2`$ at least $`-hu_2|R_i|`$ by (F5).

If $`r_i>\delta_0`$: $`n_i`$ is future timelike (F1), so $`J_i^W\ge0`$ and
$`0<\operatorname{Im}X_i\le hu_2\,R_i`$: $`u_i\le\pi+e^+_{\rm bl}`$.

If $`r_i<-\delta_0`$: I bound the negative part of $`J_i^W`$. With
$`w_i=(\tau_0+t_i,\vec x_i)`$, $`n_i=(\tau-t_i,\vec\xi-\vec x_i)`$,
$`|t_i|,|\vec x_i|\le\Delta`$, $`|\vec\xi|\le\tau+\tau_0`$ and $`\tau>0`$:
for $`\tau\ge\Delta`$,
$`n_i\!\cdot w_i\ge(\tau-\Delta)(\tau_0-\Delta)-(\tau+\tau_0+\Delta)\Delta=\tau(\tau_0-2\Delta)-2\Delta\tau_0`$,
which is $`\ge0`$ for $`\tau\ge\tau_1`$. So $`n_i\!\cdot w_i<0`$ forces
$`\tau<\tau_1`$ (note $`\tau_1\ge2\Delta`$), and then
$`n_i\!\cdot w_i\ge-\Delta(\tau_0+\Delta)-(\tau_1+\tau_0+\Delta)\Delta=-W_1`$,
while $`|\vec n_i|\ge|\vec\xi|-\Delta\ge\sqrt{R^2-\tau_1^2}-\Delta`$ and
$`|n_i^0|\le\tau_1+\Delta`$ give $`-r_i\ge\rho_1>0`$ by (B6). Hence, whenever
$`J_i^W<0`$: $`|J_i^W|\le2\kappa'(|r_i|+W_1)\le2\kappa'|r_i|(1+W_1/\rho_1)`$ if
$`n_i\!\cdot w_i<0`$, and $`|J_i^W|\le2\kappa'|r_i|`$ otherwise. So
$`0<\operatorname{Im}X_i\le\big(hu_2+\tfrac{2\kappa'(1+W_1/\rho_1)}{1-v^2/\delta_0}\big)|R_i|`$
and $`u_i\ge-e^-_{\rm bl}`$. $`\square`$

**Step 2: decomposition.** Put $`R_i:=r_i-V^2=\operatorname{Re}X_i`$,
$`\Delta^W:=-2\vartheta W\!\cdot p_{c-1}=2\kappa'\vartheta\,\xi_-\!\cdot p_{c-1}\in[0,2KP]`$,
$`\Delta^{(j)}:=J^{(j)}_{c-1}-J^{(j)}_c=-2V_j\!\cdot p_{c-1}`$, and
$`E_j:=r_{c-1}J^{(j)}_c-r_cJ^{(j)}_{c-1}`$. By Lemma R2(3) and linearity,

```math
Q:=\operatorname{Im}\big(X_{c-1}\overline{X_c}\big)=g\Big[-2\vartheta W\!\cdot N_-+V^2\Delta^W\Big]+(1-g)\sum_jm_j\Big[E_j+V^2\Delta^{(j)}\Big],
\qquad E_j+V^2\Delta^{(j)}=R_{c-1}J^{(j)}_c-R_cJ^{(j)}_{c-1}. \tag{5.1}
```

$`-2\vartheta W\!\cdot N_-=2\kappa'\vartheta\,\xi_-\!\cdot N_-\ge0`$ because
$`N_-`$ is future null or zero. So the far bracket is
$`\ge-V^2_-\Delta^W`$. I show

```math
Q\ \ge\ -s\,|X_c|\,|X_{c-1}| \tag{5.2}
```

in every case. Two bounds on a frame bracket are used repeatedly:

- (crude) $`|R_{c-1}J^{(j)}_c-R_cJ^{(j)}_{c-1}|\le\iota_2(|R_c|+|R_{c-1}|)`$;
- (signed) if the frame $`j`$ controls $`i\in\{c,c-1\}`$ then
  $`J^{(j)}_i\ge0`$.

*Closed forms for the frames $`c`$ and $`c-1`$* (Lemma R3;
$`(\alpha',\beta',\rho')`$ are the coordinates in the frame of $`c-1`$,
$`\alpha'=A/D'`$):

```math
E_c=D\big[\kappa D\beta^2d_\beta+\rho\,(2\nu\beta+\kappa d_\alpha)\big],\qquad
E_{c-1}=-|D'|\big[\kappa|D'|\alpha'^2d'_\alpha+\rho'\,(2\nu\alpha'-\kappa d'_\beta)\big],
```
```math
\Delta^{(c)}=\kappa D\,d_\alpha,\qquad\Delta^{(c-1)}=\kappa|D'|\,d'_\beta .
```

(E+) *If $`\beta>0`$, $`\kappa\phi'(\ell_3)\le2\nu\beta`$, and
$`\beta'\ge0`$ or $`\kappa\phi'(-\beta')\le2\nu|\alpha'|`$, then
$`E_c\ge0`$ and $`E_{c-1}\ge0`$.* Indeed $`d_\beta=\phi'(\ell_3)\ge0`$,
$`d_\alpha\ge-\phi'(\ell_3)`$; and $`\alpha'=\beta D/D'<0`$ gives
$`d'_\alpha=-\phi'(\ell_3')\le0`$,
$`2\nu\alpha'-\kappa d'_\beta=-2\nu|\alpha'|-\kappa\phi'(\ell'_3)+\kappa\phi'(-\beta')\le0`$.

**Step 3: the cases.** If $`i`$ is not near, $`|X_i|\ge|R_i|\ge\delta_v`$ and
$`|r_i|\le\tfrac{\delta_0}{\delta_v}|R_i|`$. If $`i`$ is near,
$`|X_i|\ge J_i`$ and $`|r_i|\le|R_i|+|V^2|\le(1+\mu_+)|X_i|`$ by (F4).

In (5.1) the bracket of the frame $`j`$ carries the weight $`m_j`$. For
a frame that does not control $`i`$, i.e. $`j=i+2`$, (F5) gives
$`m_j|J_i^{(j)}|\le hu_2|R_i|`$; for any frame $`m_j|J_i^{(j)}|\le\iota_2`$.

**Case N0: neither $`c`$ nor $`c-1`$ is near.** Then $`|R_c|,|R_{c-1}|\ge\delta_v`$
and $`R_i`$ has the sign of $`r_i`$. Far bracket:
$`\ge-\omega_0\cdot2KP\ge-T_F|X_c||X_{c-1}|`$. The frame $`c+1`$ controls
$`c`$ but not $`c-1`$ (so $`m_{c+1}|J^{(c+1)}_{c-1}|\le hu_2|R_{c-1}|`$), the
frame $`c+2`$ controls $`c-1`$ but not $`c`$ (so
$`m_{c+2}|J^{(c+2)}_c|\le hu_2|R_c|`$ and
$`0\le m_{c+2}J^{(c+2)}_{c-1}\le\chi\big(\tfrac{|r_c|}{\delta_0}-1\big)\iota_2\le hu_2|R_c|`$,
as in the proof of (F5)).

*N0pp: $`r_c>\delta_0`$.* Then $`y`$ is future timelike (F1), $`A>0`$ and
$`R_{c-1}>R_c>0`$. Frames $`c-1,c`$: as in N1b below, $`E_c,E_{c-1}\ge0`$,
and $`|V^2\Delta^{(j)}|\le v^2\pi_0`$, so their brackets are
$`\ge-T_V|X_c||X_{c-1}|`$ together. Frame $`c+1`$:
$`R_{c-1}J_c^{(c+1)}\ge0`$ and
$`-R_c\,m_{c+1}J^{(c+1)}_{c-1}\ge-hu_2R_cR_{c-1}`$. Frame $`c+2`$:
$`-R_c\,m_{c+2}J^{(c+2)}_{c-1}\ge-hu_2R_c^2\ge-hu_2R_cR_{c-1}`$ and
$`R_{c-1}\,m_{c+2}J_c^{(c+2)}\ge-hu_2R_{c-1}R_c`$. Total: (5.2) with
$`3hu_2+T_V+T_F`$.

*N0mix: $`r_{c-1}>\delta_0`$, $`r_c<-\delta_0`$.* Then $`R_{c-1}>0>R_c`$.
Frames $`c-1,c`$: both $`J`$'s are $`\ge0`$, the brackets are $`\ge0`$.
Frame $`c+1`$: $`R_{c-1}J_c^{(c+1)}\ge0`$ and
$`|R_c|\,m_{c+1}J^{(c+1)}_{c-1}\ge-hu_2|R_c|R_{c-1}`$. Frame $`c+2`$:
$`|R_c|J^{(c+2)}_{c-1}\ge0`$ and
$`R_{c-1}\,m_{c+2}J_c^{(c+2)}\ge-hu_2R_{c-1}|R_c|`$. Total: $`2hu_2+T_F`$.

*N0nn: $`r_c<-\delta_0`$, $`r_{c-1}<-\delta_0`$.* (The fourth combination,
$`r_c>\delta_0`$ with $`r_{c-1}<0`$, does not occur.) Here the brackets are
$`-|R_{c-1}|J_c^{(j)}+|R_c|J_{c-1}^{(j)}`$. Frame $`c+1`$:
$`-|R_{c-1}|\,m_{c+1}J_c^{(c+1)}\ge-\iota_2|R_{c-1}|\ge-u_2|R_c||R_{c-1}|`$ and
$`|R_c|\,m_{c+1}J^{(c+1)}_{c-1}\ge-hu_2|R_c||R_{c-1}|`$. Frame $`c+2`$:
$`-|R_{c-1}|\,m_{c+2}J_c^{(c+2)}\ge-hu_2|R_c||R_{c-1}|`$ and
$`|R_c|J^{(c+2)}_{c-1}\ge0`$. Frames $`c`$ and $`c-1`$:
$`|V^2\Delta^{(j)}|\le v^2\pi_0`$, and the closed forms satisfy

```math
E_c\ge-\varsigma_c\,|r_c||r_{c-1}|,\qquad E_{c-1}\ge-\varsigma_n\,|r_c||r_{c-1}| , \tag{5.3}
```

proved below. With $`|r_i|\le\tfrac{\delta_0}{\delta_v}|R_i|`$ the total is
(5.2) with $`(1+2h)u_2+(\varsigma_c+\varsigma_n)\delta_0^2/\delta_v^2+T_V+T_F`$.

*Proof of (5.3) for $`E_c=D[\kappa D\beta^2d_\beta+\rho(2\nu\beta+\kappa d_\alpha)]`$.*
Recall $`\rho=\alpha\beta D+|r_c|`$ and $`|r_c|,|r_{c-1}|>\delta_0`$.

If $`\beta\ge0`$: $`d_\beta=\phi'(\ell_3)\ge0`$, so
$`E_c\ge D\rho\,(2\nu\beta-\kappa\phi'(\ell_3^+))`$. If $`\alpha\ge-1`$ then
$`\ell_3\le\beta`$, $`\phi'(\ell_3^+)\le2\beta`$ and $`E_c\ge0`$ by (B2). If
$`\alpha<-1`$ then $`\alpha\beta\le0`$, $`\rho\le|r_c|`$,
$`\phi'(\ell_3^+)\le2(\beta+|\alpha|)\le4\mathfrak a`$ and
$`E_c\ge-4\kappa D\mathfrak a|r_c|\ge-\tfrac{4\kappa D\mathfrak a}{\delta_0}|r_c||r_{c-1}|`$.

If $`\beta=-b<0`$: $`|r_{c-1}|=|r_c|+bD\ge bD`$ and
$`\rho\le\alpha_-bD+|r_c|`$ with $`\alpha_-:=\max(-\alpha,0)`$. Dropping the
non-negative terms,
$`E_c\ge-D\big[\kappa Db^2\phi'(b)+\rho\,(2\nu b+\kappa\phi'(\ell_3^+))\big]`$,
and $`\ell_3^+\le\alpha_-`$. Term by term:

```math
\kappa D^2b^2\phi'(b)\le2\kappa Db^2\cdot bD\le\frac{2\kappa D\mathfrak a^2}{\delta_0}|r_c||r_{c-1}|,\qquad
2\nu\,bD\,\rho\le2\nu|r_{c-1}|\big(\alpha_-bD+|r_c|\big)\le2\nu\Big(\frac{\mathfrak a^2D}{\delta_0}+1\Big)|r_c||r_{c-1}|,
```
```math
\kappa D\rho\,\phi'(\ell_3^+)\le2\kappa D\alpha_-\big(\alpha_-bD+|r_c|\big)\le2\kappa D\alpha_-\Big(\frac{\alpha_-}{|r_c|}+\frac1{|r_{c-1}|}\Big)|r_c||r_{c-1}|\le\frac{2\kappa D\mathfrak a(\mathfrak a+1)}{\delta_0}|r_c||r_{c-1}| .
```

The sum, and the bound of the case $`\beta\ge0`$, are at most
$`\varsigma_c|r_c||r_{c-1}|`$.

*Proof of (5.3) for $`E_{c-1}=-|D'|[\kappa|D'|\alpha'^2d'_\alpha+\rho'(2\nu\alpha'-\kappa d'_\beta)]`$.*
Dropping the non-negative terms,
$`E_{c-1}\ge-|D'|\big[\kappa|D'|\alpha'^2\phi'(\alpha'_+)+\rho'\,(2\nu\alpha'_++\kappa\phi'(\beta'_-))\big]`$.
Recall $`r_c=r_{c-1}+\alpha'|D'|`$ and $`\rho'=\alpha'\beta'D'+|r_{c-1}|`$.

If $`\alpha'>0`$: $`|r_{c-1}|=|r_c|+\alpha'|D'|\ge\alpha'|D'|`$ and
$`\rho'\le\alpha'\beta'_-|D'|+|r_{c-1}|`$. Then

```math
\kappa D'^2\alpha'^2\phi'(\alpha')\le2\kappa|D'|\alpha'^2\cdot\alpha'|D'|\le\frac{2\kappa|D'|\mathfrak a'^2}{\delta_0}|r_c||r_{c-1}|,\qquad
2\nu\alpha'|D'|\rho'\le2\nu\alpha'|D'|\big(\alpha'\beta'_-|D'|+|r_{c-1}|\big)\le\frac{2\nu|D'|(\mathfrak a'^2+\mathfrak a')}{\delta_0}|r_c||r_{c-1}|,
```
```math
\kappa|D'|\rho'\phi'(\beta'_-)\le2\kappa|D'|\beta'_-\big(\alpha'\beta'_-|D'|+|r_{c-1}|\big)\le\frac{2\kappa|D'|\mathfrak a'(\mathfrak a'+1)}{\delta_0}|r_c||r_{c-1}| ,
```

using $`\alpha'|D'|\le|r_{c-1}|`$ and $`|r_c|>\delta_0`$ in each.

If $`\alpha'\le0`$: only the last term is present; if $`\beta'<0`$ then
$`\alpha'\beta'D'\le0`$, $`\rho'\le|r_{c-1}|`$, and it is at most
$`2\kappa|D'|\mathfrak a'|r_{c-1}|\le\tfrac{2\kappa|D'|\mathfrak a'}{\delta_0}|r_c||r_{c-1}|`$.
In both cases the sum is at most $`\varsigma_n|r_c||r_{c-1}|`$. $`\square`$

So (5.2) holds in Case N0 with $`s_0`$.

**Case N1: $`c`$ near, $`c-1`$ not.** The frame $`c+2`$ is off; $`c-1`$ and
$`c`$ control both invariants, $`c+1`$ controls $`c`$ only. By (F2),
$`\alpha,\beta\ge-b_0`$. Far bracket:
$`\ge-V^2_-\Delta^W\ge-\mu_-J_c\cdot2KP\ge-F|X_c||X_{c-1}|`$.

*N1a: $`r_{c-1}>\delta_0`$ and $`R_c\le0`$.* Then $`R_{c-1}>0\ge R_c`$. For
$`j=c-1,c`$ both $`J^{(j)}`$ are $`\ge0`$, so the bracket is $`\ge0`$. For
$`j=c+1`$: $`R_{c-1}J^{(c+1)}_c\ge0`$ and
$`-R_c\,m_{c+1}J^{(c+1)}_{c-1}\ge-hu_2|R_c||R_{c-1}|\ge-hu_2|X_c||X_{c-1}|`$ by (F5).

*N1b: $`r_{c-1}>\delta_0`$ and $`R_c>0`$.* If $`r_c>0`$: by (F1) $`y`$ is
future timelike, so $`\alpha,\beta>0`$ and $`\rho<\alpha\beta D`$. Then
$`\ell_3<\beta`$, so $`\kappa\phi'(\ell_3)\le2\kappa\beta\le2\nu\beta`$ by
(B2); and by Lemma R5.1 with $`\rho\le\alpha\beta D`$,
$`\beta'\ge1+\big(\alpha|D'|+\beta|u|-2\sqrt{\alpha|D'|\beta|u|}\big)/|D'|\ge1`$.
So (E+) gives $`E_c,E_{c-1}\ge0`$. If $`r_c\le0`$: $`r_{c-1}>0\ge r_c`$ and
$`E_c,E_{c-1}\ge0`$ because the $`J`$'s are $`\ge0`$. In both sub-cases
$`|V^2\Delta^{(j)}|\le\mu_+J_c\pi_0\le\tfrac{\mu_+\pi_0}{\delta_v}|X_c||X_{c-1}|`$
for $`j=c-1,c`$. For $`j=c+1`$: $`R_{c-1}J^{(c+1)}_c\ge0`$ and
$`-R_c\,m_{c+1}J^{(c+1)}_{c-1}\ge-hu_2R_cR_{c-1}\ge-hu_2|X_c||X_{c-1}|`$ by (F5).

*N1c: $`r_{c-1}<-\delta_0`$.* Then $`A=r_{c-1}-r_c<0`$, so $`\beta<0`$. By
(F2) $`|\beta|D\le\varepsilon\delta_0`$, i.e.
$`|A|\le\varepsilon\delta_0`$. Hence $`r_c\le-\delta_1`$ and
$`|r_{c-1}|\le\delta_0(1+\varepsilon)`$.

*$`\alpha>0`$ and $`\alpha\ge a_{\min}`$.* If $`\alpha\le0`$ then
$`\rho=\alpha\beta D-r_c\le b_0^2D+\delta_0`$ and (F3) gives
$`R-\Delta\le2b_0P+B_c\sqrt{\delta_0+b_0^2D}<2P+B_c\sqrt{\delta_0+b_0^2D}`$,
contradicting $`\Lambda_c\ge2`$. So $`\alpha>0`$,
$`\rho=|r_c|-\alpha|\beta|D\le\delta_0`$, and (F3) gives
$`\alpha\ge\Lambda_c-b_0=a_{\min}`$ (the frame of $`c`$, with
$`\rho\le\delta_0\le\delta_0+b_0^2D`$).

*$`\beta'\ge0`$.* Here $`c-1`$ is not near, so (F2) is applied with a larger
$`\delta`$. First $`n^0_{c-1}=n_c^0+p^0_{c-1}>0`$, because $`c`$ is near and
$`p_{c-1}`$ is future directed. With $`\delta=\delta_0(1+\varepsilon)`$ and
$`p_{c-2}`$ past directed, $`D'<0`$:
$`\beta'=2n_{c-1}\!\cdot p_{c-2}/D'\ge-b_0(1+\varepsilon)=:-b_0'`$, and
$`b_0'<2b_0\le1`$ by (B4), (B3). Also
$`\alpha'=A/D'=|A|/|D'|\in(0,b_0]`$. Suppose $`\beta'<0`$. Then

```math
\rho'=\alpha'\beta'D'-r_{c-1}=\alpha'|\beta'||D'|+|r_{c-1}|\le b_0b_0'|D'|+\delta_0(1+\varepsilon)\le\delta_0(1+\varepsilon)+2b_0^2|D'| ,
```

and (F3) in the frame of $`c-1`$ gives
$`R-\Delta\le(b_0+b_0')P+B_n\sqrt{\delta_0(1+\varepsilon)+2b_0^2|D'|}<2P+B_n\sqrt{\delta_0(1+\varepsilon)+2b_0^2|D'|}`$,
contradicting $`\Lambda_n\ge2`$. So $`\beta'\ge0`$.

*The brackets.* With $`a=\alpha`$, $`h=|\beta|`$:
$`|r_c|=ahD+\rho`$, $`|r_{c-1}|=(1+a)hD+\rho`$, $`\ell_3<0`$,
$`d_\beta=-\phi'(h)`$, $`d_\alpha=\phi'(a)`$,
$`d'_\beta=\phi'(\ell_3')\ge0`$, $`d'_\alpha\le\phi'(\alpha')\le\alpha'^2`$,
$`\rho'=|r_{c-1}|-\alpha'\beta'|D'|\le|r_{c-1}|`$,
$`\alpha'|D'|=|A|=hD\le|r_c|/a`$, $`h\le b_0`$, $`\alpha'\le b_0`$:

```math
E_c\ge-D\big[\kappa Dh^4+2\nu\rho\,h\big]\ge-\Big(\frac{\kappa h^2}{a(1+a)}+2\nu\Big)|r_c||r_{c-1}|,\qquad
E_{c-1}\ge-|A|\,|r_{c-1}|\big(\kappa\alpha'^2+2\nu\big)\ge-\frac{2\nu+\kappa b_0^2}{a}|r_c||r_{c-1}| .
```

So $`E_c+E_{c-1}\ge-\varsigma|r_c||r_{c-1}|\ge-(1+\mu_+)\tfrac{\delta_0}{\delta_v}\varsigma\,|X_c||X_{c-1}|`$.
Also $`\Delta^{(c)}=\kappa D\phi'(a)\ge0`$ and
$`\Delta^{(c-1)}=\kappa|D'|\phi'(\ell_3')\ge0`$, both at most $`\pi_0`$, so
$`V^2\Delta^{(j)}\ge-V^2_-\pi_0\ge-\tfrac{\mu_-\pi_0}{\delta_v}|X_c||X_{c-1}|`$
for $`j=c-1,c`$. For $`j=c+1`$: its weight is
$`m_{c+1}=\chi(|r_{c-1}|/\delta_0-1)\le\chi(\varepsilon)`$, because
$`|r_{c-1}|\le\delta_0(1+\varepsilon)`$; the crude bound with
$`|R_c|\ge\delta_1-v^2`$, $`|R_{c-1}|\ge\delta_v`$ gives
$`\ge-\chi(\varepsilon)\iota_2\big(\tfrac1{\delta_v}+\tfrac1{\delta_1-v^2}\big)|X_c||X_{c-1}|`$.

The case $`r_{c-1}<-\delta_0`$ with $`r_c>0`$ does not occur ($`r_c>0`$
gives $`r_{c-1}>r_c`$). So in Case N1 (5.2) holds with $`s_1`$.

**Case N2: $`c-1`$ near, $`c`$ not.** The frame $`c+1`$ is off; $`c-1,c`$
control both, $`c+2`$ controls $`c-1`$ only. By (F2) applied to
$`n_{c-1}`$: $`A=2n_{c-1}\!\cdot p_{c-1}\ge-b_0D`$,
$`1+\alpha=2n_{c-1}\!\cdot p_c/D\ge-b_0`$, and
$`\beta'=2n_{c-1}\!\cdot p_{c-2}/D'\ge-b_0`$. If $`r_c>\delta_0`$ then
$`y`$ is future timelike, $`A>0`$ and $`r_{c-1}>r_c>\delta_0`$, impossible.
So $`r_c<-\delta_0`$, $`A>0`$, $`\beta>0`$, $`\alpha'<0`$. Far bracket:
$`V^2_-\le\mu_-J_{c-1}`$, so $`\ge-F|X_c||X_{c-1}|`$. For $`j=c+2`$:
$`-R_cJ^{(c+2)}_{c-1}\ge0`$ and
$`R_{c-1}\,m_{c+2}J_c^{(c+2)}\ge-hu_2|R_{c-1}||R_c|\ge-hu_2|X_c||X_{c-1}|`$ by (F5).

For $`j=c-1,c`$: if $`R_{c-1}\ge0`$ the brackets are $`\ge0`$ (signs). Let
$`R_{c-1}<0`$. If $`r_{c-1}\ge0`$ then $`E_j\ge0`$ by the signs. If
$`r_{c-1}<0`$ I verify (E+):

- $`\kappa\phi'(\ell_3)\le2\nu\beta`$. Here $`\ell_3=\beta-(1+\alpha)\le\beta+b_0`$.
  If $`1+\alpha\ge0`$: $`\phi'(\ell_3^+)\le2\ell_3^+\le2\beta`$. If
  $`-b_0\le1+\alpha<0`$: write
  $`n_{c-1}=(1+\alpha)p_{c-1}+\beta p_c+y_\perp`$, with
  $`\rho=(1+\alpha)\beta D-r_{c-1}\le\delta_0`$; (F3) applied to this
  decomposition (the frame of $`c`$) gives
  $`R-\Delta\le(b_0+\beta)P+B_c\sqrt{\delta_0}`$, so
  $`\beta\ge\Lambda_c-b_0\ge1`$. Then $`\ell_3\ge1`$ and
  $`\phi'(\ell_3)=2\ell_3-1\le2\beta+2b_0-1\le2\beta`$ by (B3). Then use (B2).
- If $`\beta'<0`$: $`\rho'=\alpha'\beta'D'-r_{c-1}\le\delta_0`$ (because
  $`\alpha'\beta'\ge0`$, $`D'<0`$), and (F3) in the frame of $`c-1`$ gives
  $`|\alpha'|\ge\Lambda_n-b_0\ge1`$ (with $`|\beta'|\le b_0`$); so
  $`\kappa\phi'(-\beta')\le\kappa b_0^2\le2\nu|\alpha'|`$ by (B2).

So $`E_c,E_{c-1}\ge0`$. In all sub-cases
$`|V^2\Delta^{(j)}|\le\mu_+J_{c-1}\pi_0\le\tfrac{\mu_+\pi_0}{\delta_v}|X_c||X_{c-1}|`$.
So (5.2) holds with $`s_2`$.

**Case N3: both near.** Only the frames $`c-1,c`$ can be on; both control
both invariants, so $`J_c,J_{c-1}\ge0`$. $`|A|\le2\delta_0`$, so
$`|\beta|\le2\delta_0/D\le\tfrac12`$ and
$`|\alpha'|\le2\delta_0/|D'|\le\tfrac12`$ by (B3); and
$`\alpha,\beta,\beta'\ge-b_0`$ by (F2).

*$`\beta\ge0`$.* $`\ell_3=\beta-\alpha-1\le\tfrac12+b_0-1\le0`$ by (B3), so
$`d_\beta=0`$, $`d_\alpha=\phi'(\alpha)\ge0`$ and
$`E_c=D\rho(2\nu\beta+\kappa\phi'(\alpha))\ge0`$. If $`\beta'<0`$ then, as in
N2, $`|\alpha'|\ge\Lambda_n-b_0>1`$, a contradiction; so $`\beta'\ge0`$,
$`E_{c-1}\ge0`$ by (E+), and $`\Delta^{(c)},\Delta^{(c-1)}\ge0`$. Then by
(5.1)

```math
Q\ \ge\ V^2\Big(g\Delta^W+(1-g)\sum_jm_j\Delta^{(j)}\Big)=V^2\,(J_{c-1}-J_c)\ \ge\ -V^2_-\,J_{c-1}\ \ge\ -\mu_-J_cJ_{c-1} .
```

*$`\beta<0`$.* As in N1c, now with $`c-1`$ near so that (F2) applies
directly: $`\alpha\ge a_{\min}`$, $`\beta'\ge0`$,
$`E_c+E_{c-1}\ge-\varsigma|r_c||r_{c-1}|\ge-(1+\mu_+)^2\varsigma|X_c||X_{c-1}|`$,
$`\Delta^{(c)},\Delta^{(c-1)}\ge0`$, and the $`V^2`$ terms are
$`\ge-\mu_-J_cJ_{c-1}`$ as above. So (5.2) holds with $`s_3`$.

**Step 4: from (5.2) to the angles.** With $`a_i=\arg X_i\in[-\pi-e,e]`$,
(5.2) says $`\sin(a_{c-1}-a_c)\ge-s`$, and
$`a_{c-1}-a_c\in[-\pi-2e,\pi+2e]`$. So either
$`a_{c-1}-a_c\ge-\arcsin s`$, or $`a_{c-1}-a_c\le-\pi+\arcsin s`$. In the
second case $`a_c\ge-e-\arcsin s>-\tfrac\pi2`$ and
$`a_{c-1}\le-\pi+e+\arcsin s`$, so $`R_c>0>R_{c-1}`$. If $`r_c>0`$ then
$`y`$ is future timelike (F1) and $`R_{c-1}=R_c+2y\!\cdot p_{c-1}>R_c`$, a
contradiction. If $`r_c\le0`$ then $`V^2<0`$ and
$`|r_c|<V^2_-\le\omega_0<\delta_0`$, so $`c`$ is near, and by (F4)
$`0<R_c\le\mu_-J_c`$; then $`a_c\le-\tfrac\pi2+\arctan\mu_-`$, which
contradicts $`a_c\ge-e-\arcsin s`$ by (B5). Hence
$`u_c\le u_{c-1}+\arcsin s`$, and with $`u_{c+1}\ge-e^-_{\rm bl}`$ (F6):
$`K_c\le e^-_{\rm bl}+\arcsin s`$. Finally
$`\arg z=\pi+K_c+u_{c+2}\le2\pi+e^+_{\rm bl}+K_c\le2\pi+b_{\rm bl}`$ by (F6) (for
$`s>0`$; for $`t>0`$ the same for $`-\arg z`$), and the other side is Lemma
R1. $`\blacksquare`$

**Remarks.**

1. The constant is set by the sub-case N0nn: $`s\approx(1+2h)u_2+(\varsigma_c+\varsigma_n)`$,
   that is 0.08 for the default constants and 0.12 to 0.13 for the tuned
   ones; the other entries are at most 0.09. An earlier version used
   $`|J_i^{(j)}|\le\iota_2`$ for all four frames in Case N0 and $`e`$ for the
   excursions, which gave $`s\approx8u_2`$ and a blend-zone bound of 0.36
   to 0.57. What is left in N0nn is one full unit $`u_2`$, from the frame
   $`c+1`$ acting on $`X_c`$ where $`\beta<0`$; the two half units are
   frames that do not control an invariant. The measured minimum of
   $`Q/(|X_c||X_{c-1}|)`$ in the blend zone is $`-1.6\times10^{-6}`$, the
   measured maximum of $`K_c`$ is 0.0027, and $`m_{i+2}|J_i^{(i+2)}|/|R_i|`$
   never exceeds $`0.13\,hu_2`$: the unit $`u_2`$, which comes from
   $`|V_j|\le\lambda_0|n_j|`$, is itself far from sharp.
2. (B2) is used only in N0pp, N0nn, N1b and N2. It holds for the default constants and
   with equality for `tuned_contour`.

## 6. The whole contour

**Theorem.** Let one diagonal be timelike. Assume the conditions of
GATE_6_CONTOUR section 9 and

- (R4) $`\tau_0\ge3\Delta`$, $`2R>\sqrt2\,\tau_0`$;
- (R7) $`(R-\tfrac{\tau_0}2)(\tau_0-2\Delta)\ge2\Delta\tau_0`$;
- (G1) $`4\delta_0|u|<D|D'|`$; (G2) $`4\delta_0<D`$;
- (B1) $`R\ge\tau_0`$, $`R>\sqrt5\Delta`$; (B2) $`\kappa\le\nu`$;
  (B3) $`b_0\le\tfrac12`$, $`\Lambda_c\ge2`$, $`\Lambda_n\ge2`$,
  $`4\delta_0\le\min_j|D_j|`$; (B4) $`\delta_1>v^2`$; (B5) $`s<1`$,
  $`e+\arcsin s+\arctan\mu_-<\tfrac\pi2`$; (B6) $`\tau_1<R`$ and
  $`\sqrt{R^2-\tau_1^2}-\Delta>\tau_1+\Delta`$;
- the trivial ones: $`e\le\tfrac\pi4`$ (it follows from (D), which has
  $`e<0.2`$), $`4e_{\rm sp}<\pi`$, $`\epsilon_\kappa<1`$,
  $`v_R^2<\delta_0`$, $`\theta\le2`$, $`\nu<1`$, $`4\kappa^2<1`$,
  $`v^2<\delta_0/2`$, $`v^2<(\tau_0-2\Delta)^2`$, $`4e'_{\rm sp}<\pi`$.

Then on the contour, off the four edges,

```math
|\arg z|\le2\pi+b_{\rm in}\quad(r\le R),\qquad|\arg z|\le2\pi+b_{\rm bl}\quad(R<r<2R),\qquad|\arg z|\le2\pi\quad(r\ge2R),
```

with $`b_{\rm in}=e_{\rm fr}+K^{\rm fr}_{\rm in}`$ (Proposition R8) and
$`b_{\rm bl}=e^+_{\rm bl}+e^-_{\rm bl}+\arcsin s`$ (Proposition R6). In particular
$`|\arg z|\le2\pi+\max(b_{\rm in},b_{\rm bl})`$ everywhere.

*Proof.* Take $`s>0`$. For $`r\le R`$: Proposition R8. For $`R<r<2R`$:
Proposition R6 at one cusp with a timelike diagonal. For
$`r\ge2R`$: Proposition R7. The other one-sided bound,
$`\arg z\ge-\pi-4e\ge-2\pi`$, is Lemma R1 with $`e\le\tfrac\pi4`$. For $`t>0`$ the signs are exchanged.
$`\blacksquare`$

Without (R7) the outer zone still has $`|\arg z|\le2\pi+e+\max(e,3e_{\rm sp})`$
by Proposition R4. All the conditions are inequalities between the
constants that the certificate already computes. (B3) implies
$`4\delta_0\le D`$, the non-strict form of (G2), which suffices in the proofs.

**Values.** `phase_bound.phase_bound(contour)` evaluates every constant and
condition in floating point. All conditions hold at the squares below.

| Square | $`\lvert s/t\rvert`$ | $`u_2`$ | $`s`$ | $`r\le R`$: $`b_{\rm in}`$ | $`R<r<2R`$: $`b_{\rm bl}`$ | $`r\ge2R`$ | $`c_{3d}`$ | $`\min(\Lambda_c,\Lambda_n)`$ | bound $`-2\pi`$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $`(+,-)`$, angle 1.1, default constants | 1.3759 | 0.026 | 0.082 | 0.046 | 0.117 | 0.000 | 0.0029 | 4.76 | 0.117 |
| $`(-,+)`$, angle 1.1, default constants | 0.7268 | 0.026 | 0.082 | 0.046 | 0.117 | 0.000 | 0.0029 | 4.76 | 0.117 |
| $`(+,-)`$, angle 0.02, tuned | 1.0001 | 0.043 | 0.117 | 0.084 | 0.169 | 0.000 | 0.0029 | 4.54 | 0.169 |
| $`(-,+)`$, angle 0.5, tuned | 0.93879 | 0.043 | 0.117 | 0.084 | 0.169 | 0.000 | 0.0031 | 4.60 | 0.169 |
| $`(+,-)`$, angle 0.7, tuned | 1.1332 | 0.043 | 0.116 | 0.084 | 0.169 | 0.000 | 0.0033 | 4.64 | 0.169 |
| $`(-,+)`$, angle 1.1, tuned | 0.7268 | 0.044 | 0.115 | 0.084 | 0.168 | 0.000 | 0.0038 | 4.72 | 0.168 |
| $`(+,-)`$, angle 2.0, tuned | 3.4255 | 0.045 | 0.115 | 0.084 | 0.169 | 0.000 | 0.0060 | 4.91 | 0.169 |
| $`(+,-)`$, angle 2.5, tuned | 10.058 | 0.049 | 0.122 | 0.088 | 0.178 | 0.000 | 0.0076 | 4.98 | 0.178 |
| $`(-,+)`$, angle 2.9, tuned | 0.014521 | 0.052 | 0.129 | 0.094 | 0.186 | 0.000 | 0.0088 | 5.01 | 0.186 |
| $`(+,-)`$, angle 3.087, tuned | 1342.5 | 0.054 | 0.133 | 0.097 | 0.190 | 0.000 | 0.0092 | 5.01 | 0.190 |

**Scope.**

- The constants are not sharp. The blend zone decides, through the
  sub-case N0nn of Proposition R6. The observed maxima of $`K_c`$ are about
  $`10^{-6}`$, 0.003 and 0.004 in the three zones. For the constants of
  the table the phase observed on the contour stays below $`2\pi`$; for
  other admissible constants it does not (see below).
- At the ten squares of the table (floating point, constants of
  `tuned_contour` or the default ones) the bound is between $`2\pi+0.12`$
  and $`2\pi+0.19`$; section 7 gives at most $`2\pi+0.1908`$ for every
  square with a timelike diagonal, with a suitable contour at each
  square. Earlier versions had $`2\pi+0.36`$ to $`2\pi+0.57`$ (Case N0
  with the crude bound on all four frames), and before that a constant
  that was not uniform in the square.
- *The exact value.* For $`r\ge2R`$, $`|\arg z|\le2\pi`$ exactly
  (Proposition R7). For $`r<2R`$ the statement $`|\arg z|\le2\pi`$ is
  **false in general** under the conditions of the theorem.
  - *Counterexample.* Take the default contour of the region $`(+,-)`$
    with $`\delta_0`$, $`\kappa`$, $`\nu`$ multiplied by 0.4
    ($`\delta_0=0.01`$, $`\kappa=8\times10^{-9}`$, $`\nu=6\times10^{-8}`$). The
    certificate of GATE_6_CONTOUR and all conditions of the theorem hold,
    $`b_{\rm in}=0.0458`$, and at
    $`x_0=(0.04233112,\,-0.09263057,\,-6.5\times10^{-7},\,0.05679225)`$
    the phase on the contour is $`\arg z=2\pi+3.54\times10^{-7}=2\pi+5.9\,\nu`$.
    There $`r_i/\delta_0=(-1.0014,\,1.8908,\,382.07,\,1.8909)`$ and the
    weights are $`(1,\,0.967,\,0,\,0.967)`$. The value was confirmed with
    the float64 field of the repository and in multiprecision by an
    independent check.
  - *Other constants.* Largest $`\arg z-2\pi`$ found on the set of
    Proposition R9, in units of $`\nu`$, and the largest ratio of the two
    sides of (6.1) where it was recorded (constants multiplied means
    $`\delta_0,\kappa,\nu`$ multiplied by the factor):

    | Contour | $`(\arg z-2\pi)/\nu`$ | ratio in (6.1) |
    | --- | --- | --- |
    | $`(+,-)`$ 1.1, tuned, as is | $`-2.00`$ | 0.278 |
    | the same, constants $`\times0.3`$ | $`-2.00`$ | 0.691 |
    | the same, $`\times0.2`$ | $`-0.91`$ | |
    | the same, $`\times0.15`$ | $`+0.62`$ | |
    | the same, $`\times0.1`$ | $`+3.29`$ | 1.355 |
    | $`(-,+)`$ 1.1, tuned, $`\times0.1`$ | $`+3.29`$ | |
    | $`(+,-)`$ 2.0, tuned, $`\times0.3`$ | $`+1.84`$ | |
    | $`(+,-)`$ 1.1, default, as is | $`-0.33`$ | 0.963 |
    | the same, $`\times0.4`$ | $`+5.90`$ | |
    | the same, $`\times0.25`$ | $`+10.5`$ | |

  - *Mechanism.* At the worst points $`|r_{m'}|\approx\delta_0`$, so the
    frame of $`m`$ is off, and $`r_c=r_{c+2}\approx1.88\,\delta_0`$, so the
    frames of $`c`$ and of $`c+2`$ have weight about 0.96 and each acts on
    the other's invariant, which it does not control. By Proposition
    R9(4) that term contains $`2\nu\,y_\perp\!\cdot p_{j+1}`$, of size up to
    $`2\nu\sqrt{\rho\,|u||D'|/D}`$ with $`\rho\sim\delta_0`$, over an
    invariant of size $`2\delta_0`$. So the uncontrolled terms are of order
    $`\nu\sqrt{|u||D'|/(D\delta_0)}`$ while the slacks $`|C^{(j)}|`$ are of
    order $`\nu`$: the ratio grows like $`\delta_0^{-1/2}`$ as $`\delta_0`$
    is reduced, and vanishes at $`u=0`$.
  - *Consequences.* The inequality (6.1) fails for admissible constants.
    An exact bound $`2\pi`$ for $`r<2R`$ could hold only under an
    additional condition of the form $`\delta_0\ge{\rm const}\cdot|u||D'|/D`$,
    or for specific constants. Numerically it holds for the default
    constants, with margin $`0.33\,\nu`$ (not $`2\nu`$), and for
    `tuned_contour` at the squares tried (ratios in (6.1) between 0.02
    and 0.59). Where it fails, the excess is of order $`\nu`$, far below
    $`b_{\rm in}`$; the theorem is not affected.
  - Proposition R9 below gives the first-order structure behind this: the
    controlled parts are non-positive frame by frame, and the excess comes
    from the four uncontrolled terms.
- *The limits $`|s/t|\to1,0,\infty`$.* Section 7 now certifies the angles
  from 0 to $`\pi-10^{-7}`$.
  - $`|s/t|=1`$ (angle 0) is a regular point: $`u=0`$ there, no constant
    has $`u`$ in a denominator, and interval boxes containing the angle 0
    pass. So this end is covered, endpoint included.
  - Towards the angle $`\pi`$ ($`|s/t|\to\infty`$ in $`(+,-)`$, $`\to0`$ in
    $`(-,+)`$) the square degenerates: with $`d=\pi-\text{angle}`$,
    $`|D'|=d^2(1+O(d^2))`$, and for `tuned_contour`
    $`\delta_0\to0.05359\,|D'|`$, $`\kappa=\nu\to9.61\times10^{-7}\,\delta_0|D'|`$,
    $`\kappa'\to2.343\times10^{-3}\,d`$, $`B_n^2|D'|\to8`$, while
    $`R=\tau_0\to8.94427`$ and $`\Delta\to1.11803`$. Every constant of this
    document converges ($`u_2\to0.0543`$, $`c_{3d}\to0.0093`$,
    $`b_{\rm in}\to0.0980`$, $`b_{\rm bl}\to0.1910`$), so the conditions
    hold with a margin that does not shrink. The cover stops where
    float64 stops: `tuned_contour` computes the constants in floating
    point, and for $`d<10^{-7}`$ ($`|D'|<10^{-14}`$) their rounding errors
    show in the interval bound; at $`d=1.1\times10^{-8}`$ the certificate of
    GATE_6_CONTOUR itself (H1) fails for that reason.
  - The degenerate end is settled in section 7.1: with constants
    prescribed as functions of $`d`$ ($`\delta_0=\tfrac3{64}|D'|`$,
    $`\kappa=\nu=2^{-20}\delta_0|D'|`$, $`\kappa'=2^{-9}|D'|^{1/2}`$,
    $`R=\tau_0=9`$), every condition, written in scaled variables, is
    verified in interval arithmetic for all $`0<d\le0.1`$, with the bound
    at most $`2\pi+0.18826`$, and $`2\pi+0.18745`$ for $`d\le10^{-6}`$. The
    float64 constants of `tuned_contour` were not the reason the graded
    cover stopped; the way (H1) is checked unscaled was (section 7.1).
- The contour in $`D`$ dimensions (GATE_6_CONTOUR section 11) is treated
  in section 9.
- The region $`(-,-)`$ is not treated because nothing is needed there:
  $`|\arg z|\le2\pi+4e`$ follows from (D) alone.

**Proposition R9 (the exact value, to first order: the controlled parts).**
Let $`s>0`$ (for $`t>0`$ exchange the two pairs of cusps), $`r\le R`$, and
let $`x_0`$ be a real point with

```math
r_c>0,\qquad r_{c+2}>0,\qquad r_m>0>r_{m'} ,
```

where $`c,c+2`$ are the cusps with the timelike diagonal and $`m,m'`$ the
other two. (This is the set where the real phase is $`2\pi`$.) Assume (B2).
For the frame $`j`$, with $`J_i=J_i^{(j)}`$, put

```math
L^{(j)}:=\frac{J_m}{r_m}-\frac{J_c}{r_c}-\frac{J_{c+2}}{r_{c+2}}-\frac{J_{m'}}{|r_{m'}|}=C^{(j)}+U^{(j)},
```

where $`U^{(j)}`$ is the term of the invariant $`j+2`$, which the frame does
not control, and $`C^{(j)}`$ the sum of the other three. Then:

1. $`\arg z-2\pi=\varphi_m-\varphi_c-\varphi_{c+2}-\varphi_{m'}`$ with
   $`\varphi_i=\arctan(J_i/|\operatorname{Re}X_i|)`$, provided
   $`\operatorname{Re}X_i`$ has the sign of $`r_i`$. To first order in $`V`$
   this is $`\sum_jm_jL^{(j)}`$.
2. $`r_m>r_c`$ and $`r_m>r_{c+2}`$.
3. $`C^{(j)}\le0`$ for each of the four frames. More precisely
   $`C^{(m')}\le-(\kappa|D'|\Psi_{m'}+2\nu\rho_{m'})\big(\tfrac1{r_c}+\tfrac1{r_{c+2}}+\tfrac1{|r_{m'}|}\big)`$;
   $`C^{(m)}\le-\big(2\nu\rho_m+\kappa|D'|\phi'(\ell_3)(\ell_3-1)\big)/r_m`$,
   where $`\ell_3>1`$ in the frame of $`m`$; and
   $`C^{(c)}\le-J^{(c)}_{m'}/|r_{m'}|`$, likewise for $`c+2`$.
4. The uncontrolled term is given by
   $`J^{(j)}_{j+2}=J^{(j)}_j+\kappa|D_{j+1}|\,(d_\alpha-d_\beta)+2\nu\,y_\perp\!\cdot p_{j+1}`$
   (coordinates and $`d`$ of the frame $`j`$), and
   $`|y_\perp\!\cdot p_{j+1}|\le\sqrt{\rho_j|u||D_{j+1}|/|D_j|}`$.

*Proof.* (1) $`X_i=\operatorname{Re}X_i-iJ_i`$; for the three positive
invariants $`u_i=\pi-\varphi_i`$, for $`m'`$: $`u_{m'}=\varphi_{m'}`$; and
$`\arg z=\pi+u_c+u_{c+2}-u_m-u_{m'}`$.

(2) $`r_c>0`$ gives $`AB<0`$ (section 0), so exactly one neighbour of $`c`$
has a larger invariant; the neighbours are $`m`$ and $`m'`$ and
$`r_{m'}<0<r_c`$. The same at $`c+2`$.

(3) *Frame $`m'`$:* it controls $`c,m',c+2`$; all three terms are
$`-J_i/|r_i|`$ with $`J_i\ge\kappa|D'|\Psi_{m'}+2\nu\rho_{m'}`$ (Lemma 5).

*Frame $`m`$* ($`D_m=D'<0`$; it controls $`m-1,m,m+1`$, which are $`c`$ and
$`c+2`$ and $`m`$). By (2), $`r_{m+1}=r_m+\alpha|D'|<r_m`$ and
$`r_{m-1}=r_m-\beta|D'|<r_m`$, so $`\alpha=-a<0`$, $`\beta=b>0`$,
$`r_m=ab|D'|-\rho`$, and $`r_{m+1}=a(b-1)|D'|-\rho>0`$,
$`r_{m-1}=b(a-1)|D'|-\rho>0`$ give $`a,b>1`$. Then $`\ell_1,\ell_2<0`$,
$`\ell_3=a+b-1>1`$, $`d_\beta=\phi'(\ell_3)=-d_\alpha`$, and by Lemma R3
$`J_{m\pm1}=\kappa|D'|\phi'(\ell_3)\ell_3+2\nu\rho=:J_n`$,
$`J_m=J_n+\kappa|D'|\phi'(\ell_3)`$. Since $`r_{m\pm1}<r_m`$,

```math
C^{(m)}=\frac{J_m}{r_m}-\frac{J_n}{r_{m-1}}-\frac{J_n}{r_{m+1}}\le\frac{J_m-2J_n}{r_m}=-\frac{2\nu\rho+\kappa|D'|\phi'(\ell_3)(\ell_3-1)}{r_m} .
```

*Frame $`c`$* (it controls $`m,c,m'`$). If $`m=c-1`$: $`r_m=r_c+\beta D>r_c`$
gives $`\beta>0`$, and $`r_c>0`$ then gives $`\alpha>0`$, so
$`\ell_3<\beta`$, $`d_\beta=\phi'(\ell_3)\ge0`$, $`d_\alpha\ge-2\beta`$, and
by Lemma R3 and (B2)
$`r_mJ_c-r_cJ_m=D[\kappa D\beta^2d_\beta+\rho(2\nu\beta+\kappa d_\alpha)]\ge0`$,
i.e. $`J_m/r_m\le J_c/r_c`$. If $`m=c+1`$: $`\alpha<0`$, $`\beta<0`$,
$`\ell_3<|\alpha|`$, and
$`r_mJ_c-r_cJ_m=D[\kappa D\alpha^2\phi'(\ell_3)+\rho(2\nu|\alpha|-\kappa\phi'(\ell_3)+\kappa\phi'(|\beta|))]\ge0`$.
In both cases $`C^{(c)}\le-J_{m'}/|r_{m'}|\le0`$. The frame $`c+2`$ is the
same.

(4) $`n_{j+2}=n_j-p_j-p_{j+1}`$, so
$`J_{j+2}=J_j+2V_j\!\cdot(p_j+p_{j+1})`$. With
$`V_j=-\sigma\kappa(d_\beta p_{j-1}+d_\alpha p_j)+\nu y_\perp`$,
$`2p_{j-1}\!\cdot p_j=D_j`$, $`2p_{j-1}\!\cdot p_{j+1}=u`$,
$`2p_j\!\cdot p_{j+1}=D_{j+1}`$ and $`D_j+u=-D_{j+1}`$:
$`2V_j\!\cdot(p_j+p_{j+1})=\sigma\kappa D_{j+1}(d_\beta-d_\alpha)+2\nu y_\perp\!\cdot p_{j+1}`$,
and $`\sigma D_{j+1}=-|D_{j+1}|`$. The bound on $`y_\perp\!\cdot p_{j+1}`$ is
as in Lemma R5.1. $`\blacksquare`$

*The uncontrolled terms.* By Proposition R9, $`\arg z\le2\pi`$ to first
order on this set is equivalent to

```math
\sum_jm_j\,U^{(j)}\ \le\ \sum_jm_j\,\big|C^{(j)}\big| . \tag{6.1}
```

Each $`U^{(j)}`$ is explicit by (4). It can be positive and exceed the
slack $`|C^{(j)}|`$ of its own frame (on samples by a factor up to 13 for a
frame at a cusp with the timelike diagonal). The sum (6.1) is **not** a
consequence of the conditions of the theorem: it fails for the constants
of the counterexample in the scope item above (ratio of the two sides up
to 1.355 in the table there). For the default constants the largest ratio
is 0.963, at $`|r_{m'}|\approx0.96\,\delta_0`$,
$`r_c=r_{c+2}\approx1.876\,\delta_0`$; for the tuned squares tried it is
between 0.02 and 0.59. (An earlier version of this text gave 0.67 for the
default constants, from a sample that missed that region.) So Proposition
R9 stands as a first-order description of the phase near $`2\pi`$, and
(6.1) is not a statement to aim for without a further condition on
$`\delta_0`$.

## 7. The conditions in interval arithmetic

`phase_bound_iv.interval_phase_bound(region, angle_lo, angle_hi)` evaluates
the conditions (R4), (R7), (G1), (G2), (B1)–(B6) and the trivial ones, and an upper bound for the constant, with mpmath intervals, for all
squares whose angle lies in the interval, with the constants of the
midpoint square; `cover` runs it on adjacent intervals. It uses the helpers
of `certificate_iv` and calls `certificate_iv.interval_certificate` on the
same box, so a pass includes the conditions of GATE_6_CONTOUR section 9.

What a pass on a box says: for every square whose angle lies in the box,
the conditions hold for the contour built with the constants of the
midpoint square of that box. This is the convention of `certificate_iv`.
It is legitimate because the theorem holds for any constants that satisfy
the conditions; it is not a statement about `tuned_contour` evaluated at
every angle.

How each quantity is enclosed:

- $`\Delta,P,D,D',u,A_1,\lambda_0`$ and the margin $`M`$ of (C1): interval
  evaluation of their formulas on the interval square, as in
  `certificate_iv`.
- $`e`$: the upper bound of `certificate_iv` ($`\arctan x\le x`$).
- $`B_c,B_n,B`$: `certificate_iv._transverse_bound_planar`, the root of the
  trace of the Gram matrix, per frame. It is an upper bound of the
  transverse norm and is valid for this family, whose edges have no
  component along one spatial axis. $`\theta`$ and $`\mu`$ are formed from
  $`B`$.
- $`b_0,\varepsilon,\delta_1,\Lambda_c,\Lambda_n,a_{\min},\mu_\pm,\omega_0,c_{3d},u_2,\mathfrak a,\mathfrak a',\varsigma_c,\varsigma_n,\tau_1,W_1,\rho_1,s_0,\dots,s_3`$,
  $`e_{\rm fr}`$, $`e_4`$, $`\varepsilon_V`$: interval evaluation of their formulas.
- $`\arctan x`$ is replaced by $`x`$, and $`\arcsin x`$ by
  $`t-t^3/3+t^5/5`$ with $`t=x/\sqrt{1-x^2}\le1`$ (by $`t`$ if $`t>1`$).

The proof uses $`e`$, $`B`$, $`\theta`$, $`\mu`$ only as upper bounds, and
every constant is monotone in them in the direction that keeps the
conclusion: larger $`B`$ lowers $`\Lambda_c,\Lambda_n,a_{\min}`$ and raises
$`\mu_\pm`$, $`\varsigma`$ and the $`s_k`$. The interval version bounds the
three entries of the theorem: $`b_{\rm in}`$ (Proposition R8), $`b_{\rm bl}`$
(Proposition R6), and for $`r\ge2R`$ the value 0 if (R7) holds on the box
and $`e+\max(e,3e_{\rm sp})`$ otherwise. `bound_upper` is the largest of
the three.

| Region | Angle range | Boxes | $`s/t`$ covered | Result | Largest bound $`-2\pi`$ |
| --- | --- | --- | --- | --- | --- |
| $`(+,-)`$ | 0.02 to 3.087 | 12 268 | $`-1342`$ to $`-1.0001`$ | all conditions and (R7) hold on every box | 0.1908 |
| $`(-,+)`$ | 0.02 to 3.087 | 12 268 | $`-0.9999`$ to $`-0.00074`$ | all hold | 0.1908 |

On these ranges $`b_{\rm in}\le0.0977`$, $`b_{\rm bl}\le0.1908`$, $`s\le0.134`$,
$`e\le0.093`$, $`c_{3d}\le0.0093`$, $`\min(\Lambda_c,\Lambda_n)\ge4.36`$. The
largest bound is in the last box, $`[3.08675,3.087]`$; on the first three
quarters of the range it is at most 0.176. The cover was run as four
consecutive pieces of 3 067 boxes per region (24 536 boxes in all).

**Beyond the range of GATE_6_CONTOUR.** `phase_bound_iv.cover` on
$`[0,0.02]`$ and `phase_bound_iv.graded_cover` from 3.087 towards $`\pi`$
(boxes of width 0.004 times the distance from $`\pi`$) extend the cover:

| Region | Angles | Boxes | $`s/t`$ covered | Result | Largest bound $`-2\pi`$ |
| --- | --- | --- | --- | --- | --- |
| $`(+,-)`$ | 0 to 0.02 | 80 | $`-1.0001`$ to $`-1`$ | all conditions and (R7) hold on every box | 0.1692 |
| $`(-,+)`$ | 0 to 0.02 | 80 | $`-1`$ to $`-0.9999`$ | all hold | 0.1692 |
| $`(+,-)`$ | 3.087 to $`\pi-10^{-7}`$ | 3 296 | $`-4.0\times10^{14}`$ to $`-1342`$ | all hold | 0.1938 |
| $`(-,+)`$ | 3.087 to $`\pi-10^{-7}`$ | 3 296 | $`-0.00074`$ to $`-2.5\times10^{-15}`$ | all hold | 0.1938 |

Each pass includes `certificate_iv.interval_certificate` on the same box.
With the cover above: for every square with angle in
$`[0,\pi-10^{-7}]`$, in both regions with a timelike diagonal (15 644 boxes
per region), the conditions of the theorem hold for the contour with the
constants of the box midpoint, and the bound is at most $`2\pi+0.1938`$
(0.1920 for distances from $`\pi`$ of at least $`10^{-6}`$). The graded cover
continues to pass down to distance $`1.14\times10^{-8}`$
($`|s/t|=3\times10^{16}`$), with a bound that degrades to 0.33, and then
fails in (H1) of `certificate_iv`, for the reason given at the end of 7.1.
Section 7.1 replaces it near the end with a better family of constants.

*Caveat near the angle $`\pi`$.* The interval statement concerns the
mathematical contour with the constants of the box midpoint, and is
rigorous there. The float64 implementation of the field is not an accurate
realisation of that contour at the far end of the range: its edges are
null only to $`10^{-16}`$, while $`|D'|\approx10^{-14}`$ at distance $`10^{-7}`$
from $`\pi`$. Surveys of the field in float64 mean little there.

The certificate of GATE_6_CONTOUR section 9 alone, in the region $`(-,-)`$
(`graded_cover(..., certificate_only=True)`, 3 296 and 3 298 boxes),
passes from 3.087 to $`\pi-10^{-7}`$ and from 0.055 down to $`10^{-7}`$:
$`s/t`$ from $`2.5\times10^{-15}`$ to $`4.0\times10^{14}`$. It also continues
to distance $`1.1\times10^{-8}`$ at both ends and then fails in (H1). I did
not change `certificate_iv`.

Section 7.1 covers the rest of the range, up to the degenerate squares.

### 7.1 The degenerate squares in scaled variables

This closes the gap left by the covers near the degenerate end. Code:
`src/geocusp/scaled_limit.py`; tests: `tests/test_scaled_limit.py`.

**Statement.** Let $`d`$ be the distance of the angle from a degenerate end:
$`d=\pi-\text{angle}`$ in $`(+,-)`$, $`(-,+)`$ and $`(-,-)`$, and
$`d=\text{angle}`$ at the other end of $`(-,-)`$. In each of these four
cases, and for every square with $`0<d\le0.1`$, the contour with the
constants below satisfies every condition checked by
`certificate_iv.interval_certificate` (GATE_6_CONTOUR section 9). With a
timelike diagonal it also satisfies every condition checked by
`phase_bound_iv.interval_phase_bound`, together with (R7), which is what
`graded_cover` requires. On that range the bound of the theorem of section 6 is

```math
|\arg z|\le2\pi+0.18826\quad(0<d\le0.1),\qquad |\arg z|\le2\pi+0.18745\quad(0<d\le10^{-6}).
```

The range $`d\le0.1`$ overlaps the uniform covers of section 7 and
GATE_6_CONTOUR section 9, which run up to the angle 3.087
($`d\ge0.0546`$) and, in $`(-,-)`$, start at the angle 0.055. Together with
them, and with the 80-box cover of the angles $`[0,0.02]`$ in section 7
(which the uniform covers do not reach), every square with a timelike
diagonal and angle in $`[0,\pi)`$, and every square of $`(-,-)`$ with angle
in $`(0,\pi)`$, has a contour that satisfies all the conditions. The graded
cover towards $`\pi`$ is then no longer needed. Taking, for each square, the
contour of the uniform cover or of this family, the bound of the theorem
is at most $`2\pi+0.1908`$ for every square with a timelike diagonal: 0.1692
on $`[0,0.02]`$, 0.1908 on $`[0.02,3.087]`$, 0.18826 for $`d\le0.1`$.

**The family.** Use the normalisation of `light_like_square`, and put
$`\sigma=\sin(d/2)`$. The small diagonal is then exactly
$`|D'|=4\sigma^2`$: it is $`|t|`$ in $`(+,-)`$ and at the end $`\pi`$ of
$`(-,-)`$, and $`|s|`$ in $`(-,+)`$ and at the end 0 of $`(-,-)`$. The
constants are

```math
\delta_0=\tfrac3{64}\,|D'|,\qquad \kappa=\nu=2^{-20}\,\delta_0|D'|,\qquad
\kappa'=2^{-9}\,|D'|^{1/2},\qquad R=\tau_0=9 .
```

(H1) is verified at $`\epsilon=\tfrac18|D'|`$ and (I) at $`\omega=\tfrac14`$.
Both are free choices, as in section 9. The numbers are dyadic, so
$`\kappa=\nu`$ is an exact interval and (B2) is decided exactly. This is
not the family of `tuned_contour`, which is not needed: the theorems hold
for any constants that satisfy the conditions. Evaluated at
$`d=10^{-4}`$, `tuned_contour` has $`\delta_0/|D'|=0.05359`$,
$`\kappa/(\delta_0|D'|)=9.61\times10^{-7}`$, $`\kappa'/d=2.343\times10^{-3}`$
and $`R=\tau_0=8.944`$, and the phase-bound constants are $`e=0.05430`$,
$`c_{3d}=0.0093`$, $`b_{\rm in}=0.0981`$ and $`b_{\rm bl}=0.1910`$. This is
consistent with the limits quoted in section 6.

**Scaled variables.** Each quantity $`Q`$ of the two interval modules is
written

```math
Q=m(\sigma)\,\sigma^{k},
```

where $`k`$ is an explicit exponent (a multiple of $`\tfrac12`$) and $`m`$
is an explicit expression with no division by $`\sigma`$, continuous on
$`[0,\sigma_0]`$. The code carries an interval enclosure of $`m`$ over a
box of $`\sigma`$, together with $`k`$:

- A product adds the exponents.
- A sum or a maximum is formed at the smaller exponent:
  $`m_1\sigma^{k_1}+m_2\sigma^{k_2}=(m_1+m_2\sigma^{k_2-k_1})\sigma^{k_1}`$
  with $`\sigma^{k_2-k_1}\in[\sigma_{\rm lo}^{k_2-k_1},\sigma_{\rm hi}^{k_2-k_1}]`$.
- A minimum of positive terms with different exponents is formed at the
  largest exponent. The other terms then carry a factor
  $`\sigma^{-n}\ge\sigma_{\rm hi}^{-n}`$, whose upper end is infinite at
  $`\sigma=0`$.
- Every comparison is made at the common exponent $`\min(k_1,k_2)`$, that
  is after dividing both sides by $`\sigma^k>0`$, between the upper end of
  one interval and the lower end of the other.

A box may therefore contain $`\sigma=0`$, and a pass is a statement about
every $`\sigma>0`$ in it. The functions `scaled_certificate` and
`scaled_phase_bound` are line-by-line transcriptions of
`interval_certificate` and `interval_phase_bound`, with the same variable
names and the same branches.

*Geometry.* $`\Delta`$, $`P`$, the timelike $`D=2x\!\cdot y=4`$ (pairs
$`\{a,b\}`$, $`\{m,q\}`$) and their transverse bounds are computed from the
edge vectors in intervals, with $`\sin(\text{angle})=2\sigma\cos(d/2)`$ and
$`\cos(\text{angle})=\mp(1-2\sigma^2)`$. Two closed forms replace interval
evaluation, because evaluating them from the edge vectors loses the order in
$`\sigma`$ through the cancellation in $`1\pm\cos(\text{angle})`$. Take a
pair of edges other than $`\{a,b\}`$, $`\{m,q\}`$. Then $`2x\!\cdot y=-2f`$
with $`f=1\pm\cos(\text{angle})`$, which is $`2\sigma^2`$ or
$`2\cos^2(d/2)`$. The vector $`v`$ of `_transverse_bound_planar` has
$`|v|^2=2\sin^2(\text{angle})+f^2=2f(2-f)+f^2=4f-f^2`$, so its bound is

```math
B^2\le1+\frac{|v|^2}{(x\!\cdot y)^2}=\frac4f .
```

For the frames of $`D'`$ this gives $`B_n=\sqrt2\,\sigma^{-1}`$, so
$`B_n^2|D'|=8`$ exactly. Both identities are checked symbolically and against
the interval code.

*Exponents.* The exponents $`k`$ of the derived quantities, in
$`(+,-)`$:

| $`k`$ | quantities |
| --- | --- |
| $`-2`$ | $`A_1`$, $`A_2`$, $`G`$, $`w=\lvert u\rvert/\lvert D'\rvert`$, $`\mathfrak a'`$ |
| $`-1`$ | $`B`$, $`B_n`$ |
| 0 | $`\Delta`$, $`P`$, $`D`$, $`\lvert u\rvert`$, $`M`$, $`e`$, the bound of (E), $`b_0`$, $`\Lambda_c`$, $`\Lambda_n`$, $`c_{3d}`$, $`\mathfrak a`$, $`\varsigma_n`$, $`s_0,s_1,s_2`$, $`b_{\rm in}`$, $`b_{\rm bl}`$ |
| $`\tfrac12`$ | $`F`$ |
| 1 | $`\kappa'`$, $`v`$, $`K`$, $`\mu_Q`$, $`\mu_+`$, $`e_{\rm sp}`$, the outer and blend parts of (E) |
| $`\tfrac32`$ | $`\mu_-`$, $`s_3`$ |
| 2 | $`\delta_0`$, $`\delta_1`$, $`\lambda_0`$, $`\theta`$, $`\mu`$, $`v^2`$, $`a_1`$, $`\iota`$, $`\iota_2`$, $`\pi_0`$, $`\varsigma_c`$, the far lower bound of (S1), $`\delta_0-v^2`$ |
| 3 | $`\omega_0`$ |
| 4 | $`\kappa=\nu`$, $`v_R^2`$, $`\varsigma`$ |
| 6 | $`\rho_1`$ of (S) |

*(H1).* Write $`\epsilon=e|D'|`$. In a frame whose own diagonal is the large
one ($`D=4`$), the contradiction of (H1), divided by $`|D'|`$, reads

```math
\frac{4e|D'|}{D}+e+2\sqrt{\frac{|u|}{D}\Big(e+\frac{4e^2|D'|}{D}\Big)}\;<\;1-e .
```

This is continuous at $`d=0`$. There $`|u|/D=1`$ and it becomes
$`2e+2\sqrt e<1`$, that is $`e<e_*=1-\tfrac{\sqrt3}2=0.133975`$. This is why
`tuned_contour` has $`\delta_0/|D'|\to e_*/2.5=0.05359`$. `scaled_gap` is
the scaled form of `four_cone_gap`: a bisection in $`e`$ at fixed exponent.
It gives $`e_*=0.133975`$ at $`d=0`$ in all four cases, 0.1338 at
$`d=0.1`$ in $`(+,-)`$, and 0.1334 at $`d=0.1`$ in $`(-,-)`$. At
$`e=\tfrac18`$ the two sides of the inequality above are 0.832 and 0.875.
The code divides by $`\sigma^2`$ rather than by $`|D'|=4\sigma^2`$, which
gives the 3.328 and 3.5 of the table. In the frames whose own diagonal is
$`D'`$, both sides are of order 1 and the contradiction fails in the limit
(5.46 against 4), but (H1) needs only one frame. The condition $`2\delta_0<\epsilon`$ is
$`\tfrac3{32}<\tfrac18`$.

**Conditions.** The table below is for $`(+,-)`$ on the single box
$`0\le d\le10^{-6}`$, that is $`\sigma\in[0,\sin(5\times10^{-7})]`$. It
gives the exponents of the two sides, the power of $`\sigma`$ divided out,
and the scaled sides (upper end on the left, lower end on the right).
Where the left side has the higher exponent its scaled value vanishes in
the limit; in the other rows both sides are of the same order.

<details><summary>All 63 labelled comparisons, grouped</summary>

| Condition | $`k`$ left, right | divided by | left | right |
| --- | --- | --- | --- | --- |
| (H1) contradiction, frames with $`D=4`$ | 2, 2 | $`\sigma^2`$ | 3.328 | 3.5 |
| (H1) contradiction, frames with $`D=\lvert D'\rvert`$ (not needed) | 0, 0 | 1 | 5.464 | 4 |
| (H1) $`2\delta_0<\epsilon`$ | 2, 2 | $`\sigma^2`$ | 0.375 | 0.5 |
| (C0) $`3\Delta\le\tau_0`$ | 0, 0 | 1 | 3.354 | 9 |
| (C1) $`\Delta<R/\sqrt{17}`$ | 0, 0 | 1 | 1.118 | 2.183 |
| (C1) margin $`>0`$ | 0, 0 | 1 | 0 | 7.202 |
| (C2) | 0, 0 | 1 | 5.000 | 9 |
| (C3) $`\Delta<\xi_{\min}`$; worst $`<-\delta_0`$ | 0, 0; 0, 2 | 1 | 1.118; $`-47.07`$ | 8.731; $`-4.7\times10^{-14}`$ |
| (C4) | 0, 0 | 1 | 4.472 | 8.050 |
| (D) $`v^2<\delta_0/2`$ | 2, 2 | $`\sigma^2`$ | 0.01206 | 0.09375 |
| (D) $`e<0.2`$ | 0, 0 | 1 | 0.05365 | 0.2 |
| (D) far real part, three inequalities | 2, 0 | 1 | $`\le1.6\times10^{-15}`$ | $`\ge0.5`$ |
| (E) $`\Vert\partial V\Vert<1`$ | 0, 0 | 1 | 0.3012 | 1 |
| (F) $`\theta\le1`$ | 2, 0 | 1 | $`4.6\times10^{-17}`$ | 1 |
| (S0) $`b_0<\tfrac12`$ | 0, 0 | 1 | 0.1875 | 0.5 |
| (S2) | 2, 2 | $`\sigma^2`$ | 0.02344 | 0.1875 |
| (S1) far invariant $`\ge2\delta_0`$ | 2, 2 | $`\sigma^2`$ | 0.375 | 0.7031 |
| (S3) | 0, 0 | 1 | 1.118 | 9 |
| (I0) $`c_u>0`$, $`2\omega\le c_u`$ | 0, 0 | 1 | 0; 0.5 | 1.089 |
| (I1) two inequalities | 2, 0 | 1 | $`\le1.9\times10^{-15}`$ | $`\ge0.0625`$ |
| (I2) two inequalities | 0, 0 | 1 | 8.715; 0.3125 | 16.10; 18 |
| (Q1) $`\gamma<\pi/2`$; $`\mu_Q\sin\gamma<\cos\gamma`$ | 0, 0; 1, 0 | 1 | 0.4609; $`7.0\times10^{-9}`$ | 1.571; 0.8956 |
| (Q3), (Q4) | 2 or 1, 0 | 1 | $`\le2.0\times10^{-8}`$ | $`\ge1.065`$ |
| (R4) $`3\Delta\le\tau_0`$; $`\sqrt2\tau_0<2R`$ | 0, 0 | 1 | 3.354; 12.73 | 9; 18 |
| (R4) $`\epsilon_\kappa<1`$; $`e\le\pi/4`$; $`4e_{\rm sp}<\pi`$ | 2, 0; 0, 0; 1, 0 | 1 | $`2.8\times10^{-17}`$; 0.05365; $`4.9\times10^{-8}`$ | 1; 0.7854; 3.142 |
| (R7) | 0, 0 | 1 | 20.12 | 30.44 |
| (G1) $`4\delta_0\lvert u\rvert<D\lvert D'\rvert`$ | 2, 2 | $`\sigma^2`$ | 3 | 16 |
| (G2) $`4\delta_0<D`$ | 2, 0 | 1 | $`1.9\times10^{-13}`$ | 4 |
| (R5) $`v_R^2<\delta_0`$ | 4, 2 | $`\sigma^2`$ | $`6.8\times10^{-20}`$ | 0.1875 |
| (R5) $`\theta\le2`$, $`\nu<1`$, $`4\kappa^2<1`$ | 2, 4, 8; 0 | 1 | $`\le4.6\times10^{-17}`$ | $`\ge1`$ |
| (B1) $`\tau_0\le R`$; $`\sqrt5\Delta<R`$ | 0, 0 | 1 | 9; 2.5 | 9; 9 |
| (B2) $`\kappa\le\nu`$ | 4, 4 | $`\sigma^4`$ | $`7.153\times10^{-7}`$ | $`7.153\times10^{-7}`$ |
| (B3) $`b_0\le\tfrac12`$; $`\Lambda_c\ge2`$; $`\Lambda_n\ge2`$ | 0, 0 | 1 | 0.01189; 2; 2 | 0.5; 5.550; 5.087 |
| (B3) $`4\delta_0\le\min\lvert D\rvert`$ | 2, 2 | $`\sigma^2`$ | 0.75 | 4 |
| (B4) $`v^2<\delta_1`$ | 2, 2 | $`\sigma^2`$ | 0.01206 | 0.1399 |
| (R6) $`v^2<\delta_0/2`$; $`v^2<(\tau_0-2\Delta)^2`$ | 2, 2; 2, 0 | $`\sigma^2`$; 1 | 0.01206; $`3.0\times10^{-15}`$ | 0.09375; 45.75 |
| (B6) three inequalities | 0, 0 | 1 | 2.975; 4.093; 0 | 9; 7.376; 37.65 |
| (B5) $`s<1`$; $`e+\arcsin s+\mu_-<\pi/2`$; $`4e'_{\rm sp}<\pi`$ | 0, 0; 0, 0; 1, 0 | 1 | 0.1313; 0.1853; $`4.9\times10^{-8}`$ | 1; 1.571; 3.142 |

</details>

**Results.**

| Case | $`d`$ | Boxes in $`\sigma`$ | Result | $`e`$ | (E) | $`\mu_Q`$ | $`b_{\rm in}`$ | bound $`-2\pi`$ ($`=b_{\rm bl}`$) | $`\min(\Lambda_c,\Lambda_n)`$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $`(+,-)`$, $`(-,+)`$ | $`(0,10^{-6}]`$ | 1 | all hold | 0.05365 | 0.3012 | $`1.6\times10^{-8}`$ | 0.09842 | 0.18745 | 5.087 |
| $`(+,-)`$, $`(-,+)`$ | $`(0,0.1]`$ | 100 | all hold | 0.05815 | 0.3055 | 0.0016 | 0.09844 | 0.18826 | 5.087 |
| $`(-,-)`$, both ends | $`(0,10^{-6}]`$ | 1 | all hold | 0.05365 | 0.3012 | $`1.6\times10^{-8}`$ | | | |
| $`(-,-)`$, both ends | $`(0,0.1]`$ | 100 | all hold | 0.05814 | 0.3055 | 0.0016 | | | |

Upper ends are rounded up and $`\min(\Lambda_c,\Lambda_n)`$ down.

At $`d=0`$ itself the scaled bound is $`b_{\rm bl}=0.187448`$, with $`b_{\rm in}=0.098414`$,
$`s=0.1313`$ and $`c_{3d}=0.00937`$. The regions $`(+,-)`$ and $`(-,+)`$ are
mirror images. In `certificate_iv._square` the edges of $`(-,+)`$ are
$`[q,a,b,m]`$, a cyclic shift of $`[a,b,m,q]`$, so the square is the same up
to translation and relabelling, and every quantity in both modules is
symmetric under relabelling. Both regions were evaluated anyway; they agree
to $`10^{-13}`$ on the box $`d\le10^{-6}`$, and on the 100-box cover the two
bounds are 0.1882549 and 0.1882560. The two ends of $`(-,-)`$ are congruent, by $`x_3\mapsto-x_3`$ composed
with $`x\mapsto-x`$ and a relabelling, and both were evaluated (same
numbers). With the family fixed, the conditions
continue to hold up to $`d=0.5`$ in $`(+,-)`$ and up to $`d=0.35`$ in
$`(-,-)`$. Beyond that (H1) fails at $`e=\tfrac18`$ once $`e_*(d)`$ drops
below $`\tfrac18`$ (near $`d=0.77`$ and $`d=0.40`$). The family is meant for
the ends only.

**Checks** (`test_scaled_limit.py`: 29 tests, about 10 s).

- *The closed forms.* Checked symbolically with sympy for all six pairs of
  edges, and against `2*_mdot` and `_transverse_bound_planar` on the interval
  square at $`d=0.3,10^{-2},10^{-5}`$ (overlap, and agreement to
  $`10^{-25}`$ at 50 digits).
- *The rewriting against the unscaled code.* The repository's
  `interval_certificate` and `interval_phase_bound` are run, unmodified, on
  the one-point box of the angle, at $`d=0.3,10^{-3},10^{-6}`$, in each of
  the four cases. For the run, `tuned_contour` and `certificate` are replaced
  in memory so that the constants are the family's intervals and the
  repository's $`\epsilon`$ and $`\omega`$ enclose the family's. Their local
  variables are captured with `sys.setprofile`, and every one of the 65
  (certificate) or 166 (both) named intermediates is compared with the
  scaled value $`m\sigma^k`$. At 50 digits they agree to below $`10^{-30}`$
  relative. Every interval-valued local of the repository functions has a
  counterpart except `angle`. Three differences are expected. The arcsin
  bounds are points at an upper end in both codes, so they need not
  overlap, although they agree to $`10^{-30}`$. The repository's (B2) fails
  in this mode only because one non-degenerate interval is compared with
  itself. And `all_conditions` fails with it.
- *A real contour.* The same comparison was made with a `Contour` whose
  constants are the family's, rounded to float64. The repository's own
  checks pass at $`d=0.3,10^{-2},10^{-4},10^{-6}`$, and all intermediates
  agree to $`3\times10^{-16}`$.
- *The bookkeeping.* On the box $`\sigma\in[0,0.05]`$, the scaled interval of
  every named quantity contains its scaled value at points of the box.
- *Not vacuous.* Changing one constant makes the expected condition fail:
  $`\kappa=\nu`$ times 3.5 fails (E) alone, $`e=0.14`$ fails (H1), $`\delta_0=0.07|D'|`$ fails
  $`2\delta_0<\epsilon`$, $`\omega=0.6`$ fails (I0), and
  $`\tau_0=3`$ fails (C0).

**What is rigorous and what is not.**

- *Rigorous.* The interval evaluation uses outward rounding (mpmath `iv`,
  30 digits), and every comparison is between an upper and a lower end. The
  only binary numbers involved are the family's dyadic constants and the
  constants that the repository code itself uses (0.3, 0.52, 0.2, 0.375,
  ...), all converted exactly. The $`\sigma`$-box is an interval enclosure
  of $`\sin(d/2)`$. No floating-point constants of `tuned_contour` enter.
- *The two closed forms* are elementary identities, proved above and
  checked with sympy. They are not interval-evaluated.
- *The transcription* is a copy of the repository functions, and its
  agreement with them is tested, not proved. It is line by line with the
  same names, and agrees at twelve squares to $`10^{-30}`$. A reader should
  diff the two pairs of functions.
- *Not re-examined.* That the conditions imply the properties of the
  contour is the content of GATE_6_CONTOUR and of sections 1 to 6. The
  uniform covers that the range $`d\le0.1`$ is joined to are taken as
  reported.

**Not covered.**

- $`d=0`$ itself, which is not a square. The statement is for every
  $`d>0`$ in the range, uniformly.
- What is uniform is the list of inequalities and the phase constant. The
  constants themselves vanish ($`\delta_0\propto d^2`$,
  $`\kappa=\nu\propto d^4`$, $`\kappa'\propto d`$), and with them the
  quantitative lower bounds that scale with them: Lemma 6, the strips of
  Lemma 7, and Lemma 9. Nothing here bounds the integral uniformly as
  $`d\to0`$.
- The float64 implementation of the field is still not an accurate
  realisation of the contour near the end (section 7). The statement
  concerns the mathematical contour.
- $`D`$ dimensions. By section 9.4 the conditions are the same numbers. The
  scaled evaluation was not repeated separately for them.
- *Why the graded cover stopped at $`d=1.1\times10^{-8}`$.* With the
  family's constants, `graded_cover` (ratio 0.004) stops at the same box,
  $`d\in[1.1328,1.1373]\times10^{-8}`$, in (H1), so the cause is not the
  constants of `tuned_contour`. The repository checks (H1) unscaled, at
  $`\epsilon=0.99`$ times the float64 four-cone gap. On a box of that width
  the interval check passes for $`\epsilon`$ up to 0.99139 times the exact
  midpoint gap, so the float64 gap may be at most 0.14% too large. On the
  failing box it is 0.51% too large ($`1.73472\times10^{-17}`$ against
  $`1.72602\times10^{-17}`$ at 50 digits). The main reason is that the small
  diagonal in two of the frames is quantised in float64 ($`-2.22\times10^{-16}`$
  against the exact $`-1.288\times10^{-16}`$), so that frame's value exceeds
  the correct one of the large frame and is taken as the maximum; the step
  of the bisection is secondary. With a 50-digit gap the graded cover
  passes on to $`d=1.054\times10^{-8}`$, where $`1+\cos(\text{angle})`$ is
  exactly 0 in float64 and the float certificate fails. In scaled form
  (H1) has the margin 3.328 against 3.5.

**Review.** An independent review checked every operation of the scaled
arithmetic as an enclosure, diffed the transcriptions against
`interval_certificate` and `interval_phase_bound` (only the intended
differences: $`\epsilon`$, $`\omega`$, the two closed forms, the scaling),
re-derived the closed forms from the edge vectors, looked for conditions
that the theorems use but the certificate does not check (two smallness
steps in Lemmas 7 and 9 of GATE_6_CONTOUR, both implied by listed
conditions), and ran the unmodified repository functions at 60 digits with
the family's constants at ten values of $`d`$ from 0.1 to $`10^{-12}`$ in all
four cases: all pass, with bounds equal to the scaled ones to the last
digit. Verdict: proved as written, conditional on the theorems, which it
did not re-derive. The corrections it asked for (the numbers in the last
paragraph, the cover of $`[0,0.02]`$, rounding) are made.

```bash
PYTHONPATH=src python -m geocusp.scaled_limit --extend 0.1 --boxes 100
PYTHONPATH=src pytest -q tests/test_scaled_limit.py
```

## 8. Numerical checks

The field was re-implemented in 50 to 70 digit arithmetic with exactly null
edges and compared with the float64 field (relative difference at most
$`3\times10^{-13}`$). Points were aimed at the null lines through adjacent
cusps on and beyond the segments, at offsets from $`10^{-30}`$ to 1, at the
cones, and at generic points in all three zones. Unless stated otherwise
the squares are $`(+,-)`$ at angles 1.1 (default constants) and 0.7 and
$`(-,+)`$ at 1.1.

| Statement | Points per square | Violations |
| --- | --- | --- |
| Lemma R2(3) and Lemma R3 (identities) | 16 000 | 0; agreement to $`10^{-37}`$ |
| Signs of the closed forms in cases 3a, 3b, 3c; (3d.1) | 16 000 | 0 |
| Lemma R5.1(a), (b) | 15 900 | 0 |
| (3d.2), (3d.3) with $`c_{3d}`$; steps (i), (ii); squares at angles 1.1, 2.5, 2.9, 3.05, 0.05 | 8 400 to 12 700 | 0; measured at most 0.0037 against $`c_{3d}\le0.0092`$ |
| Proposition R5′: field is $`V_c`$, $`K_c\le`$ bound | 15 700 | 0; maximum equals the bound to five digits |
| Proposition R4: all statements of the proof | 32 000 | 0 |
| Proposition R7: $`n_i\!\cdot w_i\ge0`$, (a), (b), the two orderings, $`\lvert\arg z\rvert\le2\pi`$; squares at angles 1.1, 1.1, 2.5 | 29 000 | 0; largest $`\lvert\arg z\rvert/\pi`$: 1.991 |
| Proposition R6: (F1), (F2), (F4), the signs of $`E_c,E_{c-1}`$, (5.2) by case (earlier constants); also at angle 2.9 | 7 100 to 9 200 | 0 |
| Proposition R6 with the present constants, squares at angles 1.1 (default), 1.1, 2.5, 3.05, 0.05 (tuned): (F5); (F6) with its two sub-claims; (5.2) by case N0pp, N0mix, N0nn, N1, N2, N3; (5.3) for $`E_c`$ and $`E_{c-1}`$; coordinate bounds; $`K_c\le K_{\rm bl}`$; $`\lvert\arg z\rvert\le2\pi+b_{\rm bl}`$ | 9 300 to 9 700 point–cusp pairs, of which N0: 6 400 to 7 100 | 0; smallest $`Q/(\lvert X_c\rvert\lvert X_{\rm nb}\rvert)`$: $`-1.6\times10^{-6}`$; $`E_c`$, $`E_{c-1}`$ use at most $`2.5\times10^{-5}`$ of their bounds; (F5) at most 0.13 of its bound; (F6) at most 0.13 |
| Proposition R9 (float64; squares at 1.1 default, 1.1, 2.5 and 0.3 tuned; $`r\le R`$): $`C^{(j)}\le0`$ for each frame; the closed form (4) | 6 800 to 16 800 | 0 for $`C^{(j)}\le0`$ (about 130 000 frame evaluations); (4) to $`10^{-12}`$ |
| (6.1) and the phase on the set of Proposition R9 | targeted searches by an independent check, and the samples above | (6.1) fails for admissible constants (section 6); default constants: ratio at most 0.963, $`\arg z-2\pi\le-0.33\,\nu`$; tuned squares tried: ratio 0.02 to 0.59. My own samples (8 200 to 14 700 points on three squares) reached only ratio 0.67 and $`-2\nu`$: they missed the worst region |
| $`K_c`$ and $`\lvert\arg z\rvert`$ over all zones | 53 000 | largest $`K_c`$: $`5\times10^{-6}`$, 0.0027, 0.0039; largest $`\lvert\arg z\rvert/\pi`$: 1.9999999976 |

Proposition R8 was checked at 90 digits by the second review on twelve
squares (no violation). The survey in the tests also compares the sampled
phase for $`r\le R`$ with $`2\pi+b_{\rm in}`$ and $`K_c`$ with
$`K^{\rm fr}_{\rm in}`$, in float64.

In float64 the real parts of the invariants are lost within about
$`10^{-15}`$ of a null line through two adjacent cusps, and a float64 survey
shows spurious values of $`K_c`$ up to 0.06 there. The surveys in
`tests/test_phase_bound.py` therefore skip points where an invariant is
within rounding error of zero.

```bash
PYTHONPATH=src pytest -q tests/test_phase_bound.py
PYTHONPATH=src python -c "from geocusp import phase_bound_iv as p; print(p.cover('+-', 0.02, 3.087, 12268))"
```

## 9. The contour in $`D`$ dimensions

Sections 1 to 7 treat the four-dimensional contour. This section goes
through them for the contour of GATE_6_CONTOUR section 11 in
$`\mathbb R^{1,3+k}`$, $`k\ge1`$ (`feynman_contour.field_extended`), and for
its radial reduction used for $`D=4+2\eta`$.

**Result.** The theorem of section 6 holds for the $`D`$-dimensional
contour, zone by zone, with the same conditions and the same constants,
each of them the same number as for $`k=0`$. The proof is the proof of
sections 1 to 5 read in $`\mathbb R^{1,3+k}`$ with the substitutions (9.1)
below. No inequality changes and no new condition appears. Three steps
need a remark, because their four-dimensional wording uses the dimension:
the Cauchy–Schwarz inequality in Lemma R5.1; the transverse norms
$`B_c,B_n,B`$ in (F3) and in the bounds
$`|V^2|\le\tfrac\theta2|{\rm Im}|`$, $`|V_{\rm in}|^2\le3\mu I_{\rm in}`$; and
the results taken from GATE_6_CONTOUR. They are 9.2 (a), (b), (c). The
interval certificate of section 7 transfers as it stands (9.4).

**Status.** Derived line by line by an automated pass and tested as
described in 9.6. The statements (D), (C1), (C2) and Proposition 6(3) of
GATE_6_CONTOUR, which an earlier version of this section used in $`D`$
dimensions on the strength of Proposition 7 of that document (a proof by
inspection), are proved in 9.8. The facts (F5), (F6), condition (B6) and
the sharpened Case N0 with (5.3) were added to section 5 after the first
version of this section; their rows in 9.3 have been audited by the pass
that wrote this section, and 9.6 was run again against the present
constants. A referee read the first version of 9.1 to 9.7. A second
review checked 9.8 line by line, with its own 60-digit $`D`$-dimensional
field on four contours and two stress contours for $`k=1,2`$ (571 000
points, no violation), and checked the present 9.3, 9.6 and 9.7 for
consistency with sections 1 to 7; it found no mathematical error and
asked for corrections of wording, now made. The 130-digit campaign of
9.6 has not been rerun by a reviewer. Nothing here has been reviewed by a
specialist.

### 9.1 Setting and statement

GATE_6_CONTOUR section 11 writes $`e`$ for the extra coordinates and
$`\mu=e\!\cdot e`$. Both letters are in use here ($`e`$ is the excursion
bound of (D), $`\mu`$ the constant of Lemma Q), so in this section the extra
coordinates are $`\zeta\in\mathbb R^k`$ and $`m^2:=\zeta\!\cdot\zeta`$.

A point is $`\hat x=(x_0,\zeta)`$, the product is
$`(a,\zeta)\!\cdot(b,\zeta')=a\!\cdot b-\zeta\!\cdot\zeta'`$, and Euclidean
norms are those of $`\mathbb R^{4+k}`$. The cusps $`x_i`$, the edges $`p_j`$,
the centre $`c`$ and $`c_\pm`$ have $`\zeta=0`$. With
$`n_i:=\hat x-x_i=(x_0-x_i,\zeta)`$ the real invariants are

```math
r_i=n_i^2=x_{0i}^2-m^2 .
```

In the frame of a cusp $`j`$, $`\hat x-x_j=\alpha p_{j-1}+\beta p_j+\hat y_\perp`$
with the same $`\alpha,\beta`$ as in four dimensions and

```math
\hat y_\perp=(y_\perp,\zeta),\qquad\hat\rho:=-\hat y_\perp^2=\rho+m^2 . \tag{9.1}
```

The field is that of section 4 of GATE_6_CONTOUR built in
$`\mathbb R^{1,3+k}`$:

```math
V_j=-\sigma\kappa\,(d_\beta\,p_{j-1}+d_\alpha\,p_j)+\nu\,\hat y_\perp,\qquad
V_{\rm out}=-\kappa'\Theta_+(\hat x-c_-)-\kappa'\Theta_-(\hat x-c_+),
```

$`V_{\rm in}=\sum_jm_jV_j`$ with $`m_j=\chi(|r_{j+2}|/\delta_0)`$,
$`V=(1-g)V_{\rm in}+gV_{\rm out}`$; the weights use the $`D`$-dimensional
invariants $`r_i`$, the blend $`g`$ the radius $`r=|\hat x-c|`$ in
$`\mathbb R^{4+k}`$, and $`\Theta_\pm`$ the ratio of the time $`\tau`$ to the
length of the full space part $`(\vec\xi,\zeta)`$. So the extra components
of the field are $`\nu\zeta`$ for each frame and
$`-\kappa'(\Theta_++\Theta_-)\zeta`$ for the far field. This is what
`field_extended` computes (compared in 9.6). The zones are $`r\le R`$,
$`R<r<2R`$, $`r\ge2R`$ with this $`r`$. On the contour
$`X_i=r_i-V^2-iJ_i`$, $`J_i=-2V\!\cdot n_i`$, as before, and $`u_i`$, $`K_c`$,
$`z`$ are defined as in section 0.

**Theorem D.** Let one diagonal be timelike and let $`k\ge1`$. Assume the
conditions of the theorem of section 6, evaluated with the
four-dimensional constants
$`\Delta,P,A_1,\lambda_0,B_c,B_n,B,\theta,\mu,M,e`$ and the constants derived
from them. Then on the $`D`$-dimensional contour, off the four edges (which
lie in $`\zeta=0`$),

```math
|\arg z|\le2\pi+b_{\rm in}\quad(r\le R),\qquad|\arg z|\le2\pi+b_{\rm bl}\quad(R<r<2R),\qquad|\arg z|\le2\pi\quad(r\ge2R),
```

with the same $`b_{\rm in}=e_{\rm fr}+K^{\rm fr}_{\rm in}`$ and $`b_{\rm bl}`$ as
in Propositions R8 and R6. Propositions R4, R5, R5′, R6, R7 and R8 hold
as stated, with $`r`$ the $`D`$-dimensional radius; in particular, for
$`r\le R`$, $`u_i\in[-e_{\rm fr},\pi+e_{\rm fr}]`$ for all four $`i`$.

The same bounds hold on the contour in the variables $`(x_0,m)`$, $`m\ge0`$,
used for non-integer $`D=4+2\eta`$ (see 9.5).

### 9.2 The three steps that use the dimension

**(a) The transverse space.** The orthogonal complement of
$`\operatorname{span}(p_{j-1},p_j)`$ in $`\mathbb R^{1,3+k}`$ is the sum of the
four-dimensional transverse plane and the $`\zeta`$-space. It is negative
definite of dimension $`2+k`$, so $`\hat\rho\ge0`$ with equality only for
$`\hat y_\perp=0`$, and the Cauchy–Schwarz inequality
$`|a\!\cdot b|\le\sqrt{(-a^2)(-b^2)}`$ holds in it. Everything the proof says
about "the transverse plane" uses only this. In Lemma R5.1 the vectors
$`q,q'`$ (the transverse parts of $`p_{c-2}`$, $`p_{c+1}`$) have no
$`\zeta`$-component, so $`\hat y_\perp\!\cdot q=y_\perp\!\cdot q`$ and

```math
|2\hat y_\perp\!\cdot p_{c-2}|,\ |2\hat y_\perp\!\cdot p_{c+1}|\ \le\ 2\sqrt{\rho\,uD'/D}\ \le\ 2\sqrt{\hat\rho\,uD'/D} .
```

The proofs of R5.1 and of steps (ii) of case 3d use the right-hand form,
with $`\hat\rho=\alpha\beta D-r_c`$.

**(b) The transverse norms.** $`B_j`$ is the largest Euclidean norm of a
vector of Minkowski square $`-1`$ transverse to the frame $`j`$, in four
dimensions; $`B_c`$, $`B_n`$ are its values for the two kinds of frame and
$`B`$ their maximum. For such a vector $`v=(v^0,\vec v)`$:
$`|v|^2=1+2(v^0)^2\ge1`$, so $`B_j\ge1`$ for every frame. Hence, in the frame
$`j`$,

```math
|\hat y_\perp|^2=|y_\perp|^2+m^2\le B_j^2\rho+m^2\le B_j^2(\rho+m^2)=B_j^2\hat\rho ,
```

and $`B_j`$ is also the largest norm in $`2+k`$ transverse dimensions (it is
attained at $`\zeta=0`$). So $`|\hat y_\perp|\le B_j\sqrt{\hat\rho}`$ with the
same numbers $`B_c`$, $`B_n`$, $`B`$. This is what (F3) uses frame by frame,
and what Lemma 6 and $`|V_{\rm in}|^2\le3\mu I_{\rm in}`$ use with $`B`$; so
$`\Lambda_c`$, $`\Lambda_n`$, $`a_{\min}`$, $`\theta`$ and $`\mu`$ are unchanged.

**(c) What is taken from GATE_6_CONTOUR.** Its Proposition 7 states that
every condition of its section 9 holds in $`\mathbb R^{1,3+k}`$ with the
constants of $`k=0`$. The items used here:

- *Geometric constants.* $`\Delta`$, $`P`$, $`D_j`$, $`u`$ are built from
  vectors with $`\zeta=0`$. $`|\alpha|,|\beta|\le A_1|\hat x-x_j|`$ because
  $`\alpha=2n_j\!\cdot p_j/D_j`$ and $`|a\!\cdot b|\le|a||b|`$; for the same
  reason $`|\alpha|,|\beta|\le\beta_c`$ in the frame of $`c`$ for $`r\le R`$;
  $`|\hat y_\perp|\le A_2|\hat x-x_j|`$ by the triangle inequality; hence
  $`|V_j|\le\lambda_0|\hat x-x_j|`$, and $`\iota`$, $`v_R`$, $`\iota_2`$,
  $`\pi_0`$, $`v_{\rm in}`$, $`K`$, $`v`$ bound the same quantities. The bound
  on $`|V_j|`$ includes the component $`\nu\zeta`$.
- *Lemma 5* with $`\hat\rho`$: $`x_i-x_j`$ lies in the span of the two edges,
  so $`n_i=(\alpha-v_\alpha)p_{j-1}+(\beta-v_\beta)p_j+\hat y_\perp`$ and
  $`V_j\!\cdot n_i=-\sigma\kappa\tfrac{D_j}{2}d\!\cdot(P-v_i)-\nu\hat\rho`$,
  using only $`p^2=0`$, $`\hat y_\perp\!\cdot p=0`$,
  $`\hat y_\perp^2=-\hat\rho`$.
- *Lemma 6* and *$`|V_{\rm in}|^2\le3\mu I_{\rm in}`$*: their proofs use Lemma
  5, $`|d|^2\le16\Psi`$, $`|d_\beta p_{j-1}+d_\alpha p_j|\le\sqrt2P|d|`$ and
  $`|\hat y_\perp|\le B\sqrt{\hat\rho}`$, i.e. (b).
- *(H1)*: the far invariant is
  $`r_{j+2}=D'(1+\alpha-\beta)+(\alpha\beta D-\hat\rho)-2\hat y_\perp\!\cdot q'`$
  (from Lemma R5.1), and $`\max_i|r_i|\le\epsilon`$ gives
  $`\rho\le\hat\rho\le\epsilon+4\epsilon^2/|D|`$, so the bound of (H1) holds
  by (a).
- *(C1) with its margin $`M`$, (C2), (D), Proposition 6(3)*: proved in 9.8
  for every $`k`$, with the constants of $`k=0`$. Their derivations are
  inequalities between $`\tau`$, the length of the space part (now
  including $`\zeta`$), $`\Delta`$, $`\tau_0`$, $`R`$, and products with the
  vectors $`x_i-c_\mp`$, which have $`\zeta=0`$.

### 9.3 Sections 1 to 6, statement by statement

"Unchanged" means: the text holds word for word with vectors in
$`\mathbb R^{1,3+k}`$, $`y_\perp\to\hat y_\perp`$, $`\rho\to\hat\rho`$ (also
$`\rho',\rho''`$ in the neighbouring frames), norms in $`\mathbb R^{4+k}`$
and $`r`$ the $`D`$-dimensional radius.

| Statement | In $`D`$ dimensions | Reason |
| --- | --- | --- |
| Section 0 | (9.1): $`r_c=\alpha\beta D-\hat\rho`$, $`r_{c\mp1}=r_c+A,\ r_c+B`$, $`AB=-D(r_c+\hat\rho)`$ | $`n_{c-1}=y+p_{c-1}`$, $`n_{c+1}=y-p_c`$ and $`\hat y_\perp\!\cdot p=0`$ |
| Lemma R1 | unchanged | arithmetic of the $`u_i`$ and (D); at real points $`r_c>0`$ gives $`AB<0`$ since $`\hat\rho\ge0`$ |
| Lemma R2 (1), (3) | unchanged | bilinear identities |
| Lemma R2 (2) | unchanged | $`N_\mp^2=0`$ is an identity; a non-zero null vector with positive product with a null vector $`p`$ has the orientation of $`p`$, and a vector orthogonal to a null vector has square $`\le0`$, in $`\mathbb R^{1,n}`$ for every $`n`$ |
| Lemma R2 (4) | unchanged | complex numbers only |
| Lemma R3 | holds with $`\hat\rho`$: $`J_i=\kappa\lvert D_j\rvert d\!\cdot(P-v_i)+2\nu(\rho+m^2)`$, $`V_j^2=\kappa^2D_jd_\beta d_\alpha-\nu^2(\rho+m^2)`$, and $`\rho\to\rho+m^2`$ in both closed forms | Lemma 5 in the form of 9.2 (c); the square uses $`\hat y_\perp^2=-\hat\rho`$ |
| R4, preamble and (i) | unchanged | $`J^\pm(c_\mp)`$ are the cones of $`\mathbb R^{1,3+k}`$; $`r\le\lvert\tau\rvert+\lvert(\vec\xi,\zeta)\rvert`$; $`\hat x-c_-`$ future causal, so $`V^2\ge0`$; Lemma R2 |
| R4, signs of the imaginary parts | unchanged | triangle inequality for the space part, with $`\vec x_i`$ padded by zeros |
| R4, conclusion; (ii) | unchanged | $`r_{\rm nb}=r_c+2\lvert y\!\cdot p\rvert`$; time reversal does not touch $`\zeta`$ |
| R4 (iii) | unchanged | $`-r_i=(\lvert\vec n_i\rvert-\lvert n_i^0\rvert)(\lvert\vec n_i\rvert+\lvert n_i^0\rvert)`$ with $`\vec n_i\in\mathbb R^{3+k}`$; $`w_i=x_i-c_\mp`$ has $`\zeta=0`$ |
| R7, orderings $`u_c\le u_m`$, $`u_{c+2}\le u_m`$ | unchanged | R4 (i); orientations of the edges |
| R7, $`\tau\ge R-\tau_0/2`$ | unchanged | $`r\le\tau+\lvert(\vec\xi,\zeta)\rvert\le2\tau+\tau_0`$ |
| R7, $`n_i\!\cdot w_i\ge0`$ | unchanged | $`w_i=(\tau_0+t_i,\vec x_i,0)`$ has no $`\zeta`$-component, so $`n_i\!\cdot w_i=(\tau-t_i)(\tau_0+t_i)-(\vec\xi-\vec x_i)\!\cdot\vec x_i`$ contains only the four-dimensional part $`\vec\xi`$ of the space part, and $`\lvert\vec\xi\rvert\le\lvert(\vec\xi,\zeta)\rvert\le\tau+\tau_0`$; then (R7) |
| R7, $`J_i=2\kappa'\vartheta(r_i+n_i\!\cdot w_i)`$, (a), (b) | unchanged, with $`r_i=x_{0i}^2-m^2`$ | $`\hat x-c_-=n_i+w_i`$ as vectors of $`\mathbb R^{1,3+k}`$, so $`n_i\!\cdot(\hat x-c_-)=n_i^2+n_i\!\cdot w_i`$. This is where the component $`-\kappa'\Theta\zeta`$ of the far field enters: it supplies the term $`-m^2`$ of $`r_i`$. (a) and (b) use $`V^2\ge0`$ |
| R7, elsewhere | unchanged | R4 (iii) |
| Section 4, constants $`\iota,v_R,e_4,\beta_c,G,e_{\rm fr},c_{3d}`$ | unchanged | 9.2 (c) |
| Lemma R5.1, formulas | unchanged, with $`2y_\perp\!\cdot p`$ read as $`2\hat y_\perp\!\cdot p`$ (the same number) | differences of invariants |
| Lemma R5.1, inequality, (a), (b) | holds with $`\hat\rho`$ | 9.2 (a); (a), (b) use $`\hat\rho=\alpha\beta D-r_c\le\lvert\alpha\beta\rvert D+\delta_0`$ |
| R5, cases 1, 2 | unchanged | $`\lvert V^2\rvert\le v_R^2`$, $`\lvert J_c\rvert\le4\iota`$, $`AB<0`$ |
| R5, case 3, linearised numbers | unchanged | Lemma 6 |
| R5, 3a, 3b | hold with $`\hat\rho,\hat\rho',\hat\rho''`$ | the closed forms are $`\ge0`$ because $`\hat\rho,\hat\rho',\hat\rho''\ge0`$; Lemma R5.1; if $`J_c=J_{\rm nb}=0`$, Lemma 5 forces $`\Psi_j=0`$ and $`\hat\rho_j=0`$ for every frame that is on, hence $`V=0`$ |
| R5, 3c | holds: $`V=\nu\hat y_\perp`$, $`J_{c-1}=J_c=J_{c+1}=2\nu\hat\rho`$, $`V^2=-\nu^2\hat\rho`$ | (H1) by 9.2 (c); $`\hat\rho=0`$ means $`y_\perp=0`$ and $`\zeta=0`$ |
| R5, 3d, frame $`c`$, (3d.1) | holds with $`\lvert r_c\rvert=abD+\hat\rho`$, $`\hat\rho\le\lvert r_c\rvert<\delta_0`$ | Lemma R3 |
| R5, 3d, frame $`c-1`$: step (i) $`bD/\lvert r_c\rvert\le G`$ | unchanged | uses $`aD=\lvert r_{c+1}\rvert-\lvert r_c\rvert`$, $`\lvert r_c\rvert\ge abD`$ and $`b\le\beta_c`$ |
| step (ii) $`\beta'_-\le b\lvert u\rvert/\lvert D'\rvert`$ | holds with $`\hat\rho`$ | Lemma R5.1 and $`2\sqrt{\hat\rho\lvert u\rvert/(D\lvert D'\rvert)}\le1-\gamma_1`$ because $`\hat\rho<\delta_0`$ |
| the five estimates and (3d.2) | hold with $`\hat\rho'=\lvert r_{c-1}\rvert-\alpha'\beta'\lvert D'\rvert\le\lvert r_{c-1}\rvert+\alpha'\beta'_-\lvert D'\rvert`$ | algebra on Lemma R3, (i), (ii) |
| R5, 3d, frame $`c+1`$, (3d.3); sum, (3d.4) | hold with $`\hat\rho''`$ | the same |
| Proposition R8, first part | unchanged: still only the frame $`i+2`$ can make $`J_i<0`$ | the frames $`j\ne i+2`$ contribute $`\kappa\lvert D_j\rvert d\!\cdot(P-v_i)+2\nu\hat\rho_j\ge0`$ by Lemma 5; their $`\zeta`$-component $`\nu\zeta`$ adds $`2\nu m^2\ge0`$. For the frame $`i+2`$, $`\lvert J_i^{(i+2)}\rvert\le2\lvert V_{i+2}\rvert\lvert n_i\rvert\le\iota`$ with norms in $`\mathbb R^{4+k}`$, which includes its $`\zeta`$-component, and it is on only where $`\lvert r_i\rvert>\delta_0`$ |
| Proposition R8, second part | unchanged | the proof of R5 with $`e_{\rm fr}`$ in place of $`e`$ |
| Proposition R5′ | holds with $`\hat\rho`$; same constant | as 3c and 3d |
| Section 5, constants | unchanged, including $`B_c,B_n,\Lambda_c,\Lambda_n,a_{\min}`$ | 9.2 (b), (c) |
| R6, step 0, symmetries | unchanged | as R4; the relabelling $`x_{c+k}\mapsto x_{c-k}`$ does not touch $`\zeta`$ |
| (F1) | unchanged | as in R4 |
| (F2) | unchanged | $`\lvert\vec n\rvert-n^0\le\delta/(\lvert\vec n\rvert+n^0)`$ and $`\lvert\vec n\rvert+n^0\ge\lvert n\rvert`$; a null edge has $`p^0=\lvert p\rvert/\sqrt2`$ and $`\vec n\!\cdot\vec p`$ involves only the four-dimensional part of $`\vec n`$; (C1), (C2) by 9.2 (c) |
| (F3), per frame | holds with $`\hat\rho`$ and the same $`B_c`$, $`B_n`$ | triangle inequality and 9.2 (b) |
| (F4), Lemma V | unchanged | $`\lvert V_{\rm in}\rvert^2\le3\mu J_i^{\rm in}`$ by 9.2 (c), $`J_i^W\ge2\kappa'M`$ by (F2), $`0\le W^2\le K^2`$ |
| (5.1), (5.2), crude and signed bounds, (E+) | unchanged, closed forms with $`\hat\rho,\hat\rho'`$ | Lemma R2(3), Lemma R3 |
| (F5) | unchanged | the weight is $`m_{i+2}=\chi(\lvert r_i\rvert/\delta_0-1)`$ with the $`D`$-dimensional $`r_i`$; $`\chi(t)\le h(1+t)`$ is a statement about the cutoff; $`\lvert J_i^{(i+2)}\rvert\le2\lvert V_{i+2}\rvert\lvert n_i\rvert\le\iota_2`$ and $`\lvert R_i\rvert\ge\lvert r_i\rvert-v^2`$ with norms in $`\mathbb R^{4+k}`$ |
| (F6): $`J_i<0`$ forces $`\lvert r_i\rvert>\delta_0`$; the case $`r_i>\delta_0`$ | unchanged | (F2), Lemma 5 with $`\hat\rho`$, (F5); $`n_i`$ and $`\hat x-c_-`$ are future causal in $`\mathbb R^{1,3+k}`$, so $`J_i^W\ge0`$ |
| (F6): $`J_i^W=2\kappa'(r_i+n_i\!\cdot w_i)`$ | unchanged, with $`r_i=x_{0i}^2-m^2`$ | as in R7: $`\hat x-c_-=n_i+w_i`$; the component $`-\kappa'\Theta\zeta`$ of the far field supplies $`-m^2`$ |
| (F6): $`n_i\!\cdot w_i<0`$ forces $`\tau<\tau_1`$, and then $`n_i\!\cdot w_i\ge-W_1`$ | unchanged; same $`\tau_1,W_1`$ | $`w_i`$ has no $`\zeta`$-component, so $`n_i\!\cdot w_i=(\tau-t_i)(\tau_0+t_i)-(\vec\xi-\vec x_i)\!\cdot\vec x_i`$ contains only $`\vec\xi`$, and $`\lvert\vec\xi\rvert\le\lvert(\vec\xi,\zeta)\rvert\le\tau+\tau_0`$ in $`J^+(c_-)`$ |
| (F6): $`-r_i\ge\rho_1`$ for $`\tau<\tau_1`$; the bound $`e^-_{\rm bl}`$ | unchanged; same $`\rho_1`$, $`e^\pm_{\rm bl}`$; (B6) is a condition on the same numbers | $`\lvert(\vec\xi,\zeta)\rvert^2=r^2-\tau^2\ge R^2-\tau_1^2`$, the space part of $`n_i`$ in $`\mathbb R^{3+k}`$ has length $`\ge\lvert(\vec\xi,\zeta)\rvert-\Delta`$, and $`\lvert n_i^0\rvert\le\tau_1+\Delta`$ |
| N0, preamble | unchanged | (F5); $`0\le m_{c+2}J^{(c+2)}_{c-1}`$ by Lemma 5; $`\lvert R_i\rvert\ge\delta_v`$ |
| N0pp | unchanged | (F1): $`y`$ is future timelike in $`\mathbb R^{1,3+k}`$, so $`\alpha,\beta>0`$ and $`\hat\rho<\alpha\beta D`$; then (E+) as in N1b |
| N0mix | unchanged | signs of Lemma 5; (F5) |
| N0nn, (5.3) | hold with $`\hat\rho=\alpha\beta D+\lvert r_c\rvert`$, $`\hat\rho'=\alpha'\beta'D'+\lvert r_{c-1}\rvert`$ and the same $`\mathfrak a,\mathfrak a',\varsigma_c,\varsigma_n`$ | algebra on the closed forms of Lemma R3 with $`\hat\rho,\hat\rho'\ge0`$; it uses neither $`J^\pm(c_\mp)`$ nor the orientation of the edges; $`\lvert\alpha\rvert,\lvert\beta\rvert\le2\lvert n_c\rvert P/D\le\mathfrak a`$ with $`\lvert n_c\rvert\le2R+\Delta`$ in $`\mathbb R^{4+k}`$, likewise $`\mathfrak a'`$ |
| N1a, N1b | unchanged | $`y`$ future timelike gives $`\alpha,\beta>0`$ and $`\hat\rho<\alpha\beta D`$; Lemma R5.1 with $`\hat\rho`$; (F5) for the frame $`c+1`$ |
| N1c | holds with $`\hat\rho\le\delta_0+b_0^2D`$ and $`\hat\rho'\le\delta_0(1+\varepsilon)+2b_0^2\lvert D'\rvert`$; $`\Lambda_c`$, $`\Lambda_n`$ as in four dimensions; $`m_{c+1}\le\chi(\varepsilon)`$ | (F2), (F3) in the frames of $`c`$ and $`c-1`$; the weight uses the $`D`$-dimensional $`\lvert r_{c-1}\rvert\le\delta_0(1+\varepsilon)`$ |
| N2, N3 | hold with $`\hat\rho,\hat\rho'\le\delta_0`$ | (F2), (F3), (E+); (F5) for the frame $`c+2`$ in N2 |
| Step 4 | unchanged | complex numbers, (F1), (F4), (B5); the last lines use (F6) |
| Theorem, zone by zone; list of conditions | unchanged; every constant is the number computed for $`k=0`$ | R8, R1 with R6, R7; the other side $`\arg z\ge-\pi-4e`$ is R1 |

Answers to the questions a reader will ask:

- The extra component of $`V_j`$ is $`\nu\zeta`$; it changes $`J_i`$ by
  $`2\nu m^2`$ for the three invariants of the frame and $`V_j^2`$ by
  $`-\nu^2m^2`$, i.e. $`\rho\to\rho+m^2`$ in both, and in nothing else.
- The far field acts on the extra coordinates as $`-\kappa'\Theta_\pm\zeta`$:
  it is the $`(4+k)`$-vector $`-\kappa'\Theta_\pm(\hat x-c_\mp)`$. So $`V`$ is
  proportional to $`\hat x-c_-`$, which is causal in $`J^+(c_-)`$ by
  definition; this gives $`V^2\ge0`$ in R4 and R7,
  $`W^2\ge0`$ in R6, and $`J_i=2\kappa'\vartheta(r_i+n_i\!\cdot w_i)`$ with the
  $`D`$-dimensional $`r_i`$ in R7. With the extra components dropped that
  identity fails (9.6).
- $`B_c,B_n\ge1`$ is needed, and is true, for (F3) and for the constants
  $`\theta`$, $`\mu`$. It is not the only point: Lemma R5.1 needs 9.2 (a).

### 9.4 The interval certificate of section 7

A pass of `interval_phase_bound` on a box is a statement about numbers:
for every square with angle in the box, the conditions of the theorem hold
for the constants of the midpoint, and the constant is at most the value
reported. By 9.2 the constants and conditions of Theorem D are these same
numbers, and the conditions of GATE_6_CONTOUR section 9 that the pass
includes are the conditions of the $`D`$-dimensional contour: for (H1),
(C0) to (C3) and (D) by 9.8, for the others by Proposition 7 of that
document (they are not used in sections 1 to 6). The enclosure of $`B_c,B_n,B`$ (root of the trace of the Gram
matrix of the four-dimensional transverse basis) is an upper bound of the
four-dimensional norms, hence of the $`D`$-dimensional ones, which are
equal to them. So the cover of section 7 holds for the $`D`$-dimensional
contour without recomputation: for every square with angle in
$`[0,\pi-10^{-7}]`$, $`|\arg z|\le2\pi+0.1938`$ on the contour in
$`\mathbb R^{1,3+k}`$ built with the constants of the midpoint of its box,
and (R7) holds on every box. The interval version now bounds the three
zone entries $`b_{\rm in}`$, $`b_{\rm bl}`$ and 0. The cover was not rerun
in $`D`$ dimensions; nothing in it depends on $`k`$.

### 9.5 All $`k`$ at once, and non-integer dimension

The field is radial in $`\zeta`$: $`V_4`$ depends on $`\zeta`$ through
$`m^2`$ only, and the extra component is $`\lambda\,\zeta`$ with the scalar
$`\lambda=(1-g)\nu\sum_jm_j-g\kappa'(\Theta_++\Theta_-)`$, a function of
$`(x_0,m^2)`$. Hence $`(\zeta+i\lambda\zeta)^2=m^2(1+i\lambda)^2`$ and the
four numbers $`X_i`$ at $`(x_0,\zeta)`$ are those of $`k=1`$ at
$`(x_0,|\zeta|)`$. So Theorem D for $`k=1`$ implies it for every $`k`$,
independently of 9.3. For $`D=4+2\eta`$ the contour of GATE_6_CONTOUR
section 11 is the $`k=1`$ contour restricted to $`m\ge0`$, and Theorem D is a
pointwise statement, so it holds there. Nothing is said here about the
integral over $`m`$ or about the face $`m=0`$ in Stokes' theorem; those are
the subject of GATE_6_CONTOUR section 11.

### 9.6 Numerical checks

The checks of this subsection were run against the present text of
sections 1 to 6 and the present constants of
`phase_bound.phase_bound` ($`s_0,\dots,s_3`$, $`e^\pm_{\rm bl}`$,
$`K_{\rm bl}`$, $`b_{\rm bl}`$, (B6)). They replace an earlier campaign
(1 406 096 points, no violation) that was run before the blend zone was
sharpened.

The $`(4+k)`$-dimensional field was implemented again in 130-digit decimal
arithmetic with exactly null edges (Python `decimal`). The comparisons of
this implementation with others were made for the earlier campaign and
were not repeated, the field being unchanged: an mpmath evaluation at 160
digits, written separately, agrees with it to $`10^{-86}`$ relative on about
5 960 invariants for each of two squares, down to moduli of $`10^{-45}`$;
it agrees with `field_extended` in float64 to $`1.3\times10^{-12}`$ relative
in the field and $`2.5\times10^{-13}`$ in the invariants (400 generic points
for each of six squares and $`k=1,2,3`$); and for $`k=2,3`$ the float64
invariants at $`(x_0,\zeta)`$ equal those of $`k=1`$ at $`(x_0,|\zeta|)`$ to
$`6\times10^{-14}`$ (9.5).

Points: generic points at scales $`10^{-1}`$ to $`10^{1.6}`$; points on the
null lines through adjacent cusps, on the segments, at the cusps and up to
16 edge lengths beyond, with four-dimensional offsets from $`10^{-30}`$ to
$`0.3`$ or zero; points on the $`D`$-dimensional cones
$`x_{0i}^2=m^2`$ at distances $`10^{-3}`$ to $`300`$ from the cusp, with the
null direction generic, almost along $`\zeta`$, or almost four-dimensional,
and offsets $`10^{-30}`$ to $`0.3`$; points in the frame of a cusp with
$`\alpha,\beta`$ and the transverse offset spread over twelve decades and
both signs, and the same with lengths scaled by $`\delta_0`$ so that
$`|r_c|\sim\delta_0`$; points at $`r=R(1\pm\epsilon)`$ and $`2R(1\pm\epsilon)`$
and at $`|r_i|=\delta_0(1\pm\epsilon)`$, $`2\delta_0(1\pm\epsilon)`$ with
$`\epsilon`$ from $`10^{-25}`$ to $`10^{-1}`$, reached by rescaling or by
choosing $`m`$; points aimed at the thin shell of case N1c; points in the
wedge $`0.15\le|\tau|/|(\vec\xi,\zeta)|\le0.6`$; and, new in this campaign,
two families aimed at Case N0 and at the weight transition: points with
$`r_i=\pm\delta_0(1+t)`$, $`t`$ uniform in $`(0,1.3)`$ or between $`10^{-8}`$
and 1, at a prescribed distance from $`x_i`$ along a null direction of
$`\mathbb R^{1,3+k}`$ (the frame $`i+2`$ partly on; these give N0pp and
N0mix); and points far out on the null line through two adjacent cusps
with $`\hat\rho=\delta_0(1+t)`$, split in any proportion between the
four-dimensional transverse plane and $`\zeta`$, so that both invariants
are $`-\delta_0(1+t)`$ up to small longitudinal offsets (N0nn with two
weights in transition). The extra coordinate is zero, between $`10^{-30}`$
and $`10^{-10}`$, comparable with the four-dimensional offset within a
factor $`10^{\pm3}`$, or large: $`m^2`$ ranges up to $`10^9`$.

The squares include three extreme tuned ones: $`(-,+)`$ at angle 2.9
($`|s/t|=0.0145`$, $`\delta_0=0.0031`$, $`\kappa=\nu=1.7\times10^{-10}`$,
$`B_n=11.7`$), $`(+,-)`$ at angle 3.05 ($`\delta_0=4.5\times10^{-4}`$,
$`\kappa=\nu=3.6\times10^{-12}`$, $`B_n=30.9`$), and $`(+,-)`$ at angle 0.05
($`|s/t|`$ close to 1, $`\delta_0=0.26`$). All conditions of the theorem,
including (R7) and (B6), hold at the seven squares used.

| Square | $`k`$ | Points | inner / blend / outer | Largest $`K_c`$: inner, blend, outer | Largest $`\lvert\arg z\rvert/\pi`$: inner, blend, outer | $`b_{\rm in}`$; $`b_{\rm bl}`$; $`s`$ |
| --- | --- | --- | --- | --- | --- | --- |
| $`(+,-)`$, 1.1, default | 1 | 81 246 | 35 815 / 24 308 / 21 123 | $`1.7\times10^{-6}`$, 0.0027, 0.0040 | 1.9999999045, 1.9999942, 1.9795 | 0.046; 0.117; 0.082 |
| $`(-,+)`$, 1.1, default | 1 | 76 647 | 23 232 / 37 321 / 16 094 | $`1.1\times10^{-6}`$, 0.0027, 0.0039 | 1.9999999045, 1.9999945, 1.9805 | 0.046; 0.117; 0.082 |
| $`(+,-)`$, 0.7, tuned | 1 | 80 208 | 24 763 / 39 528 / 15 917 | $`1.4\times10^{-6}`$, 0.0027, 0.0040 | 1.9999996151, 1.9999361, 1.9755 | 0.084; 0.169; 0.116 |
| $`(-,+)`$, 0.5, tuned | 1 | 51 463 | 19 069 / 22 588 / 9 806 | $`1.5\times10^{-6}`$, 0.0027, 0.0040 | 1.9999995339, 1.9999260, 1.9753 | 0.084; 0.169; 0.117 |
| $`(-,+)`$, 2.9, tuned | 1 | 82 273 | 38 046 / 25 569 / 18 658 | $`1.0\times10^{-5}`$, 0.00076, 0.0011 | 1.9999999999, 1.9999995, 1.9965 | 0.094; 0.186; 0.129 |
| $`(+,-)`$, 0.05, tuned | 1 | 60 435 | 17 777 / 29 951 / 12 707 | $`2.0\times10^{-6}`$, 0.0027, 0.0040 | 1.9999993570, 1.9998895, 1.9702 | 0.084; 0.169; 0.117 |
| $`(+,-)`$, 3.05, tuned | 1 | 125 900 | 55 184 / 43 082 / 27 634 | $`8.7\times10^{-6}`$, 0.00029, 0.00042 | $`2-2.3\times10^{-12}`$, 1.9999999, 1.9985 | 0.097; 0.189; 0.132 |
| $`(+,-)`$, 1.1, default | 2 | 54 700 | 19 174 / 24 755 / 10 771 | $`9.0\times10^{-7}`$, 0.0027, 0.0039 | 1.9999999045, 1.9999942, 1.9785 | 0.046; 0.117; 0.082 |
| $`(-,+)`$, 2.9, tuned | 2 | 54 364 | 21 505 / 24 843 / 8 016 | $`8.9\times10^{-6}`$, 0.00076, 0.0011 | 1.9999999999, 1.9999995, 1.9962 | 0.094; 0.186; 0.129 |
| $`(+,-)`$, 3.05, tuned (blend families only) | 2 | 26 850 | 5 842 / 17 058 / 3 950 | $`4\times10^{-10}`$, 0.00029, 0.00029 | 1.9999999998, 1.9999999, 1.9983 | 0.097; 0.189; 0.132 |
| $`(+,-)`$, 0.05, tuned (transition families only) | 2 | 27 278 | 4 237 / 16 299 / 6 742 | $`<0`$, 0.00096, 0.0037 | 1.9999923, 1.9998572, 1.9688 | 0.084; 0.169; 0.117 |
| $`(-,+)`$, 1.1, default | 3 | 51 068 | 17 285 / 23 474 / 10 309 | $`7.9\times10^{-7}`$, 0.0026, 0.0040 | 1.9999999045, 1.9999937, 1.9797 | 0.046; 0.117; 0.082 |
| $`(+,-)`$, 0.05, tuned | 3 | 52 118 | 17 298 / 23 944 / 10 876 | $`2.0\times10^{-6}`$, 0.0027, 0.0040 | 1.9999993570, 1.9998897, 1.9704 | 0.084; 0.169; 0.117 |
| $`(-,+)`$, 2.9, tuned (transition families only) | 3 | 26 277 | 5 513 / 16 340 / 4 424 | $`<0`$, $`<0`$, $`<0`$ | 1.9999986, 1.9996180, 1.5135 | 0.094; 0.186; 0.129 |
| $`(+,-)`$, 0.7, tuned (blend families only) | 3 | 24 648 | 3 741 / 15 815 / 5 092 | $`9.3\times10^{-7}`$, 0.0027, 0.0026 | 1.9999996038, 1.9999316, 1.9721 | 0.084; 0.169; 0.116 |

In all, 875 475 points (308 481 inner, 384 875 blend, 182 119 outer;
558 172 with $`k=1`$, 163 192 with $`k=2`$, 154 111 with $`k=3`$); 2 118 of
them lie within $`10^{-40}`$ of an edge and are excluded from the
statements about phases. Every statement below was tested at each point
where its hypotheses hold, at both cusps with a timelike diagonal. There
is no violation of any of them. The largest $`|\arg z|`$ found is
$`2\pi-7\times10^{-12}`$, in the inner zone of the square at angle 3.05; in
the blend zone it is $`2\pi-2.5\times10^{-7}`$ and in the outer zone
$`2\pi-0.0047`$.

| Statement, in the form of 9.3 | Tests | Remarks |
| --- | --- | --- |
| $`\hat\rho=\rho+m^2`$; $`r_j=\alpha\beta D-\hat\rho`$, $`r_{j\mp1}`$; $`AB=-D(r_c+\hat\rho)`$ | 3 501 900 | identities, to $`10^{-95}(1+r^2)^2`$ absolute |
| Lemma R2 (1), (2), (3), with the actual field | 1 750 950 | identities; orientation of $`N_\mp`$ |
| Lemma R3: $`J_i`$ with $`2\nu\hat\rho`$ (three per frame), both closed forms, $`V_j^2`$ | 10 505 700; 3 501 900 each | identities |
| Lemma 5 with $`\hat\rho`$ | 10 505 700 | |
| $`\lvert\hat y_\perp\rvert\le B_j\sqrt{\hat\rho}`$ and (F3) with the norm of each frame; $`\lvert\alpha\rvert,\lvert\beta\rvert\le A_1\lvert\hat x-x_j\rvert`$, $`\lvert V_j\rvert\le\lambda_0\lvert\hat x-x_j\rvert`$ | 3 501 900 each | four-dimensional $`B_c,B_n,A_1,\lambda_0`$ |
| (D): $`u_i\in[-e,\pi+e]`$ | 3 493 428 | largest excursion 0.0040; see 9.8 |
| Proposition 6(3): no $`X_i`$ vanishes off the edges | 875 475 | see 9.8 |
| Lemma R5.1: formulas; $`\lvert2\hat y_\perp\!\cdot p\rvert\le2\sqrt{\rho uD'/D}\le2\sqrt{\hat\rho uD'/D}`$ | 1 750 950 | |
| Lemma R5.1 (a), (b) | 315 861; 316 370 | smallest margin 0.0011 |
| Lemma 6 and $`\lvert V\rvert^2\le3\mu J_i`$, $`r\le R`$, $`\lvert r_i\rvert\le\delta_0`$ | 599 620 | largest $`\lvert V^2\rvert/J_i`$: $`0.075\cdot\theta/2`$ |
| $`\lvert\alpha\rvert,\lvert\beta\rvert\le\beta_c`$ in the frame of $`c`$, $`r\le R`$ | 612 726 | |
| R5 cases 1, 2 | 363 880; 106 467 | |
| 3a, 3b: frames off, $`\ell_3<0`$, closed form of $`E_c`$ and its sign, sign of $`E_{c\mp1}`$, $`u_c\le u_{c\mp1}+2\varepsilon_V`$ | 19 120; 17 347 | |
| 3c: $`V=\nu\hat y_\perp`$, $`J=2\nu\hat\rho`$, $`V^2=-\nu^2\hat\rho`$, $`K_c\le0`$ | 9 364 | |
| 3d: $`\lvert r_c\rvert=abD+\hat\rho`$; $`2\sqrt{\hat\rho\lvert u\rvert/(D\lvert D'\rvert)}\le1-\gamma_1`$; (3d.1); (3d.4) with $`c_{3d}`$ | 17 193 | |
| 3d, frame $`c-1`$: (i); (ii); $`\hat\rho'\le\lvert r_{c-1}\rvert+\alpha'\beta'_-\lvert D'\rvert`$; (3d.2) with $`c_{3d}`$ | 8 472 | left side of (i) at most $`0.46\,G`$; of (3d.2) at most $`0.38\,c_{3d}`$ |
| 3d, frame $`c+1`$: (i); (ii); $`\hat\rho''`$; (3d.3) with $`c_{3d}`$ | 8 139 | at most $`0.48\,G`$; $`0.40\,c_{3d}`$ |
| R8: frames $`j\ne i+2`$ contribute $`\ge0`$ to $`J_i`$; $`u_i\in[-e_{\rm fr},\pi+e_{\rm fr}]`$ | 1 225 452 | |
| R8: $`J_i<0`$ implies $`m_{i+2}>0`$, $`\lvert r_i\rvert>\delta_0`$, $`\lvert J_i\rvert\le m_{i+2}\lvert J_i^{(i+2)}\rvert\le\iota`$, $`\lvert\operatorname{Re}X_i\rvert\ge\delta_0-v_R^2`$ | 82 324 | |
| R8: $`K_c\le K^{\rm fr}_{\rm in}`$; R5: $`K_c\le K_{\rm in}`$ | 612 726 | |
| R5′: field is $`V_c`$; $`K_c\le`$ bound | 33 289 | largest $`K_c`$ equals the bound to five digits or better |
| R4 (i), (ii): $`V^2\ge0`$, $`\operatorname{Im}(X_{\rm nb}\overline{X_c})\ge0`$, signs of $`J_i`$, $`u_c\le u_{\rm nb}`$ | 305 406 | all four combinations of zone and orientation |
| R4 (iii); $`K_c\le\max(e,3e_{\rm sp})`$ | 58 832; 364 238 | |
| R7 in $`J^\pm(c_\mp)`$: $`\lvert\tau\rvert\ge R-\tau_0/2`$; $`n_i\!\cdot w_i\ge0`$ (four $`i`$); $`J_i=2\kappa'\vartheta(r_i+n_i\!\cdot w_i)`$; all $`u_i\le\pi`$; (a); (b); the common neighbour $`m`$ with $`u_c,u_{c+2}\le u_m`$ | 152 703 each | |
| R6 step 0 (a): $`\lvert u_i\rvert\le e'_{\rm sp}`$, $`\lvert\arg z\rvert\le\pi+4e'_{\rm sp}`$; (b), (c) | 50 786; 718 964 | |
| (F1) | 995 000 | |
| (F2): $`\Theta=1`$, $`\lvert n_i^0\rvert\ge R/\sqrt5-\Delta`$, $`J_i^W\ge2\kappa'M`$; products with the four edges | 327 262 | |
| (F4) Lemma V; $`\lvert V_{\rm in}\rvert^2\le3\mu J_i^{\rm in}`$; $`V^2_-\le\omega_0`$ | 327 262; 718 964 | largest $`V^2_-/(\mu_-J_i)`$: 0.0005 |
| (F5): $`m_{i+2}\lvert J_i^{(i+2)}\rvert\le hu_2\lvert R_i\rvert`$, all four $`i`$, the whole blend zone | 1 539 500 | 1 375 869 with $`m_{i+2}>0`$, of which 227 880 with $`0<t<1`$; largest ratio of the two sides 0.41 |
| (F5), steps: $`\chi(t)\le h(1+t)`$, $`\lvert J_i^{(i+2)}\rvert\le\iota_2`$, $`\lvert R_i\rvert\ge(1+t)(\delta_0-v^2)`$ | 1 375 869 | smallest $`h(1+t)-\chi(t)`$: 0.01596 (the minimum of the polynomial is 0.0159) |
| (F6): $`-e^-_{\rm bl}\le u_i\le\pi+e^+_{\rm bl}`$ in $`J^\pm(c_\mp)`$; $`J_i^W=2\kappa'(r_i+n_i\!\cdot w_i)`$ | 1 437 928 | both orientations; smallest distance of $`u_i`$ from the ends 0.013 |
| (F6): $`n_i\!\cdot w_i<0`$ implies $`\lvert\tau\rvert<\tau_1`$, $`n_i\!\cdot w_i\ge-W_1`$, $`-r_i\ge\rho_1`$ | 6 740 | smallest $`-r_i/\rho_1`$: 1.62 |
| (F6): $`J_i<0`$ implies $`\lvert r_i\rvert>\delta_0`$, $`R_i`$ has the sign of $`r_i`$, $`\lvert R_i\rvert\ge\lvert r_i\rvert(1-v^2/\delta_0)`$; for $`r_i<-\delta_0`$ the bound on the negative part of $`J_i^W`$ and $`\operatorname{Im}X_i\le\tan e^-_{\rm bl}\lvert R_i\rvert`$ | 81 405 | all with $`r_i<-\delta_0`$; largest $`\operatorname{Im}X_i/(\tan e^-_{\rm bl}\lvert R_i\rvert)`$: 0.13 |
| (F6), the case $`r_i>\delta_0`$ with $`J_i<0`$ | 0 | never reached, also not by 154 089 further points placed just outside $`r=R`$ |
| N0, preamble: $`\lvert R_c\rvert,\lvert R_{\rm nb}\rvert\ge\delta_v`$ with the signs of $`r`$; the three uses of (F5); $`0\le m_{c+2}J^{(c+2)}_{\rm nb}\le\chi(\lvert r_c\rvert/\delta_0-1)\iota_2`$ | 593 078 | all four combinations |
| N0pp: $`r_{\rm nb}>r_c`$, $`R_{\rm nb}>R_c>0`$, $`E_c,E_{\rm nb}\ge0`$; (5.2) with $`3hu_2+T_V+T_F`$ | 231 135 | 103 561 with a weight in transition; largest $`-Q/(\lvert X_c\rvert\lvert X_{\rm nb}\rvert)`$: $`4\times10^{-6}`$ of the entry |
| N0mix: $`R_{\rm nb}>0>R_c`$, brackets of the frames $`c`$, nb $`\ge0`$; (5.2) with $`2hu_2+T_F`$ | 149 494 | 94 254 with a weight in transition; $`Q>0`$ at every point |
| N0nn: both below $`-\delta_0`$ (no fourth combination); (5.2) with the third entry of $`s_0`$ | 212 449 | 100 003 with a weight in transition; largest $`-Q/(\lvert X_c\rvert\lvert X_{\rm nb}\rvert)`$: $`1.5\times10^{-5}`$ of the entry |
| (5.3) for both neighbours wherever $`r_c,r_{\rm nb}<-\delta_0`$, $`R<r<2R`$; $`\lvert\alpha\rvert,\lvert\beta\rvert\le\mathfrak a`$, $`\lvert\alpha'\rvert,\lvert\beta'\rvert\le\mathfrak a'`$ | 613 199 | largest $`-E_c/(\varsigma_c\lvert r_c\rvert\lvert r_{\rm nb}\rvert)`$: 0.25; largest $`-E_{\rm nb}/(\varsigma_n\lvert r_c\rvert\lvert r_{\rm nb}\rvert)`$: 0.79 |
| (5.2) with $`s_1,s_2,s_3`$ in cases N1, N2, N3 | 59 280; 38 884; 27 722 | all four combinations; largest $`-Q/(\lvert X_c\rvert\lvert X_{\rm nb}\rvert)`$ over the blend zone: $`1.7\times10^{-6}`$ |
| Step 4: $`u_c\le u_{\rm nb}+\arcsin s`$; $`K_c\le K_{\rm bl}`$ | 718 964; 769 750 | largest $`K_c`$ 0.0027 against $`K_{\rm bl}\ge0.10`$ |
| N1: $`\alpha,\beta\ge-b_0`$; N1b: $`E_c,E_{c-1}\ge0`$, $`\beta'\ge1`$ | 14 558; 827 | standard combination only |
| N1c: $`\beta<0`$, $`\lvert A\rvert\le\varepsilon\delta_0`$, $`\alpha\ge a_{\min}`$, $`\beta'\ge0`$, $`E_c+E_{c-1}\ge-\varsigma\lvert r_c\rvert\lvert r_{c-1}\rvert`$; $`\lvert r_{c-1}\rvert\le\delta_0(1+\varepsilon)`$, $`m_{c+1}\le\chi(\varepsilon)`$ | 2 464 | standard combination only |
| N2: the listed inequalities; $`E_c,E_{c-1}\ge0`$; the two hypotheses of (E+) when $`r_{c-1}<0`$ | 9 964; 7 357 | standard combination only |
| N2: $`1+\alpha<0`$ implies $`\beta\ge\Lambda_c-b_0`$; $`\beta'<0`$ implies $`\lvert\alpha'\rvert\ge\Lambda_n-b_0`$ | 118; 1 974 | standard combination only |
| N3, $`\beta\ge0`$ and $`\beta<0`$ | 5 142; 1 738 | standard combination only |
| Theorem D: $`\lvert\arg z\rvert\le2\pi+b_{\rm in}`$; $`2\pi+b_{\rm bl}`$; $`2\pi`$ | 306 363; 384 875; 182 119 | by zone; no phase exceeds $`2\pi`$; smallest margin in the outer zone 0.0047 |
| The other side: $`\pm\arg z\ge-\pi-4e`$ | 873 357 | |

Remarks on the tests.

- Sign statements are tested with an absolute slack of
  $`10^{-95}(1+r^2)^2`$, and comparisons with a constant computed in float64
  with a relative slack of $`10^{-12}`$. Coordinates and invariants below
  the rounding level are not classified as non-zero, and $`\gamma_1`$ is
  formed in decimal arithmetic; these precautions come from three
  artefacts of rounding met in the earlier campaign at 110 digits (points
  exactly on a null line or cone; points $`10^{-30}`$ beyond a cusp at the
  square of angle 3.05, where the true
  $`|X_i|\approx\kappa|D|\Psi\approx10^{-101}`$; $`\hat\rho=\delta_0(1-10^{-25})`$).
  No statement was flagged in the present campaign.
- For Proposition R5′ the ratio of the largest $`K_c`$ to the bound is
  1.00001 at the square of angle 3.05, where the bound is
  $`7\times10^{-12}`$ and the angles are formed in float64; the test allows
  $`10^{-15}`$ absolute.
- The detailed statements of Case N0 (the signs in N0pp and N0mix, (5.3),
  the entry of $`s_0`$ by sub-case) were tested in all four combinations of
  zone and orientation, with the neighbour of section 5; (5.3) also for the
  other neighbour.
- As a control, the same program was run with the extra components of the
  field removed (the "real shift" of GATE_6_CONTOUR section 11), on 60 000
  points of the square at angle 2.9 and 60 000 blend-zone points of the
  default square. Then the identities of Lemma R3 fail at 88 to 90% of
  the points, Lemma 5 at 26 to 35%, the identity
  $`J_i=2\kappa'\vartheta(r_i+n_i\!\cdot w_i)`$ of R7 at 95 to 96%, and the closed
  forms in 3a, 3b, 3c at 50 to 83%. No bound on $`K_c`$, on $`\arg z`$ or in
  (5.2), (F5), (F6) is violated in those samples, so a survey of the phase
  alone would not distinguish the two fields.

A survey in float64 with `feynman_contour.displaced_extended` itself
($`k=1`$; the points above together with `sample_extended`; points where a
real or displaced invariant is below $`10^{-9}(1+r^2)`$ skipped, about 17%)
gives, on 94 501, 93 290 and 94 750 points for the squares $`(+,-)`$ 1.1
default, $`(-,+)`$ 2.9 tuned and $`(+,-)`$ 0.05 tuned: no violation of
$`K_c\le K^{\rm fr}_{\rm in},K_{\rm bl},\max(e,3e_{\rm sp})`$ with the present
$`K_{\rm bl}`$, or of the three zone bounds $`b_{\rm in}`$, $`b_{\rm bl}`$, 0
on $`|\arg z|-2\pi`$; largest $`|\arg z|/\pi`$ by zone at most
1.9999999999, 1.9999994, 1.9960.

The programs for these checks are not part of the repository.

### 9.7 What is not covered

- The sub-claim of (F6) for $`r_i>\delta_0`$ with $`J_i<0`$ was never
  reached by a sample point, in 9.6 or in a further run aimed at it. It
  rests on its proof: (F1), $`J_i^W\ge0`$ and (F5).
- The detailed sign statements of cases N1, N2, N3 were tested only in the
  standard combination of section 5 (future edges at $`c`$,
  $`\hat x\in J^+(c_-)`$, neighbour $`c-1`$); in the other three
  combinations only (F1), (F2), (F4), (F5), (F6), (5.2) by case and the
  conclusions were tested. The branch "$`1+\alpha<0`$" of N2 was reached at
  only 118 points. Case N0 was tested in all four combinations.
- The interval cover of section 7 (angles from 0 to $`\pi-10^{-7}`$) was
  not rerun; 9.4 is an argument that it applies, not a computation. As in
  four dimensions the limits $`|s/t|\to0,\infty`$ and the region $`(-,-)`$
  are not treated.
- Proposition R9 and the counterexample to $`|\arg z|\le2\pi`$ for
  $`r<2R`$ (section 6) concern four dimensions; Proposition R9 has not
  been checked for $`k\ge1`$. The observation in 9.6 that no sampled phase
  exceeds $`2\pi`$ is about the seven squares and constants used there.
- The numerical checks use $`k=1,2,3`$ and seven squares. Non-integer
  dimension is covered only through the reduction of 9.5, which is
  pointwise; nothing here concerns the integrals.
- Of the conditions of GATE_6_CONTOUR section 9, 9.8 proves those that
  sections 1 to 6 use: (H1), (A), (C0) to (C2), (D), Proposition 6(3), with
  Lemma 5, and also (C3), which is not used; Lemma 6 and $`|V_{\rm in}|^2\le3\mu I_{\rm in}`$ are in 9.2 (c). The
  remaining ones, (E), (F), (S), (I) and (Q), are not used here and still
  rest on Proposition 7 of that document.
- This section audits the proof of sections 1 to 6 for its dependence on
  the dimension. It is not a further review of that proof.

### 9.8 The conditions of GATE_6_CONTOUR in $`D`$ dimensions, written out

Section 9 uses (D), (C1), (C2) and Proposition 6(3) of GATE_6_CONTOUR for
the contour in $`\mathbb R^{1,3+k}`$. Proposition 7 of that document asserts
them by inspection of the four-dimensional derivations. This subsection
proves them from the definitions, for every $`k\ge0`$, with the constants
that `feynman_contour.certificate` computes. The proofs also need (H1) and
(A), which are included; (C3) is proved as well, for completeness, but no
proof here or in sections 1 to 6 uses it ((C2) already excludes the wedge
near a cone). The case $`k=0`$ is the four-dimensional
statement; nothing below uses $`k`$.

**Result.** All of them hold in $`\mathbb R^{1,3+k}`$ as stated, and every
constant is literally the four-dimensional number (table at the end).
Nothing fails. Two steps of the four-dimensional text of (D) are
treated here in more detail; neither depends on the dimension. First, how
the frame term and the far term combine: the proof below gives the convex
combination $`(1-g)\tan e_{\rm frames}+g\tan e_{\rm far}`$. The sum used in
the four-dimensional text is also valid, since the imaginary part is at
most $`(\tan e_{\rm frames}+\tan e_{\rm far})|\operatorname{Re}X_i|`$ and
$`\arctan`$ is subadditive on $`[0,\infty)`$; the convex combination is
only sharper. Second, the real part in the zone $`r\ge2R`$ when
$`n\!\cdot u<0`$ or $`0<\Theta<1`$, where $`V^2`$ shifts it (Lemma F, last
claim). Here this uses $`\tan e_{\rm far}\le1`$, which follows from the
condition $`e<0.2`$ that the certificate checks under (D) (it is also named
in GATE_6_CONTOUR section 13). GATE_6_CONTOUR section 9 now closes the same
step by another route, the paragraph beginning "In the pure outer zone"
with the certificate condition `D_far_real_part`. Both routes prove the
same inequality, $`|\operatorname{Re}X_i|\ge\tfrac12|r_i|`$ wherever the far
term is positive and $`r\ge2R`$, and give the same $`e`$: on the route
here $`e<0.2`$ implies the first two inequalities of `D_far_real_part`,
and the third is not needed because the exact $`V_{\rm out}^2`$ is used;
on the other route $`e<0.2`$ is not needed for this step. A review
confirmed both routes and that the second also holds verbatim in
$`\mathbb R^{1,3+k}`$.

**Status.** Derived by the automated pass that wrote section 9 and tested
as described at the end. Reviewed once (see the status of 9.1); not by a
specialist.

#### Setting

Notation of 9.1: $`\hat x=(x_0,\zeta)`$, $`\zeta\in\mathbb R^k`$,
$`m^2=\zeta\!\cdot\zeta`$; the cusps, the edges, $`c`$ and $`c_\pm=c\pm\tau_0e_t`$
have $`\zeta=0`$. Write

```math
\hat x-c=(\tau,\hat\xi),\quad\hat\xi=(\vec\xi,\zeta)\in\mathbb R^{3+k},\quad r=|\hat x-c|,\qquad
n=n_i=\hat x-x_i=(n^0,\hat n),\qquad x_i-c=(t_i,\vec x_i,0),
```

so $`t_i^2+|\vec x_i|^2\le\Delta^2`$, $`|n^0-\tau|\le\Delta`$,
$`\big||\hat n|-|\hat\xi|\big|\le\Delta`$, and $`r_i=n^2=(n^0)^2-|\hat n|^2`$. Put
$`u^\mp=u_i^\mp:=x_i-c_\mp=(t_i\pm\tau_0,\vec x_i,0)`$; then
$`\hat x-c_\mp=n+u^\mp`$. For two vectors $`a\!\cdot b`$ is the Minkowski
product of $`\mathbb R^{1,3+k}`$, $`|a\!\cdot b|\le|a||b|`$ and
$`|a^2|\le|a|^2`$ with Euclidean norms.

The field (GATE_6_CONTOUR sections 4, 9, 11; `field_extended`):
$`V=(1-g)V_{\rm in}+gV_{\rm out}`$, $`g=0`$ for $`r\le R`$, $`g=1`$ for
$`r\ge2R`$, $`0<g<1`$ in between;
$`V_{\rm in}=\sum_jm_jV_j`$, $`m_j=\chi(|r_{j+2}|/\delta_0)\in[0,1]`$, which
is 0 for $`|r_{j+2}|\le\delta_0`$ and 1 for $`|r_{j+2}|\ge2\delta_0`$;
$`V_{\rm out}=-\kappa'\Theta_+(n+u^-)-\kappa'\Theta_-(n+u^+)`$, where
$`\Theta_\pm=1`$ for $`\pm\tau\ge\tfrac12|\hat\xi|`$, $`\Theta_\pm=0`$ for
$`\pm\tau\le\tfrac14|\hat\xi|`$, $`0<\Theta_\pm<1`$ in between; so
$`\Theta_+\Theta_-=0`$. On the contour
$`X_i=r_i-V^2+2i\,V\!\cdot n_i`$.

Constants (all computed by the certificate from the four-dimensional
square):

```math
A_1=\frac{2P}{\min_j|D_j|},\quad A_2=1+2A_1P,\quad\lambda_0=16\kappa A_1P+\nu A_2,\quad
\iota_2=2\lambda_0(2R+\Delta)^2,\quad v=4\lambda_0(2R+\Delta)+\kappa'(2R+\tau_0+\Delta),
```
```math
M=(\tau_0-2\Delta)\Big(\frac R{\sqrt{17}}-\Delta\Big)-\delta_0-\Delta\sqrt{\delta_0},\qquad\xi_0=2(2\Delta+\sqrt{\delta_0}),\qquad\xi_c=\frac{2R}{\sqrt5},
```
```math
c_b=1+\frac{2\Delta}{3(R/\sqrt{17}-\Delta)},\qquad c_c=1+\Big(\frac{\sqrt5}2+\frac\Delta{\xi_c}\Big)\frac{\tau_0+\Delta}{\tfrac34(\xi_c-4\Delta)},
```
```math
\tan e_{\rm frames}=\frac{\iota_2}{\delta_0-v^2},\qquad\tan e_{\rm far}=4\kappa'\max(1,c_b,c_c),\qquad e=e_{\rm frames}+e_{\rm far},
```

($`e_{\rm frames}`$ and $`e_{\rm far}`$ are the certificate's
`excursion_frames` and `excursion_far`. The first is not the constant
$`e_{\rm fr}`$ of section 4 and Proposition R8, which uses
$`\iota=2\lambda_0(R+\Delta)^2`$ for $`r\le R`$ in place of $`\iota_2`$ and is
smaller by about a factor of four.)

and $`\epsilon_*`$ of (H1). Conditions used, each one a line of the
certificate: (H1) $`2\delta_0<\epsilon_*`$; (C0) $`\tau_0\ge3\Delta`$; (C1)
$`R/\sqrt{17}>\Delta`$ and $`M>0`$; (C2) $`R\ge\tfrac{\sqrt5}2\xi_0`$; (C3)
$`\xi_m:=4R/\sqrt{17}>\Delta`$ and
$`(\tfrac14\xi_m+\Delta)^2-(\xi_m-\Delta)^2<-\delta_0`$; (D) $`v^2<\delta_0/2`$,
$`\xi_c>4\Delta`$ and $`e<0.2`$.

#### Lemma 5, norms, (A), (H1)

**(L5) Lemma 5 with $`\hat\rho`$.** In the frame of $`j`$,
$`n_j=\alpha p_{j-1}+\beta p_j+\hat y_\perp`$ with
$`\hat y_\perp=(y_\perp,\zeta)`$ orthogonal to both edges and
$`\hat y_\perp^2=-\hat\rho`$, $`\hat\rho=\rho+m^2\ge0`$. Since
$`n_{j-1}=n_j+p_{j-1}`$ and $`n_{j+1}=n_j-p_j`$, the coordinates of $`n_i`$
($`i=j-1,j,j+1`$) are $`P-v_i`$ in the order $`(\beta,\alpha)`$, with the
same transverse part. With
$`V_j=-\sigma\kappa(d_\beta p_{j-1}+d_\alpha p_j)+\nu\hat y_\perp`$,
$`\sigma=\operatorname{sign}D_j`$, $`p^2=0`$ and $`2p_{j-1}\!\cdot p_j=D_j`$:

```math
V_j\!\cdot n_i=-\sigma\kappa\frac{D_j}2\,d\!\cdot(P-v_i)-\nu\hat\rho=-\frac{\kappa|D_j|}2\,d\!\cdot(P-v_i)-\nu\hat\rho\ \le\ -\frac{\kappa|D_j|}2\Psi_j-\nu\hat\rho_j\ \le0 ,
```

by convexity of $`\Psi`$ and $`\Psi(v_i)=0`$. If the right side is 0 then
$`\Psi_j=0`$ and $`\hat\rho_j=0`$, hence $`d=\nabla\Psi=0`$ (a minimum of a
$`C^1`$ function), $`\hat y_\perp=0`$ and $`V_j=0`$; and $`\hat x`$ lies in the
flat triangle $`\operatorname{conv}(x_{j-1},x_j,x_{j+1})`$, in particular
$`\zeta=0`$.

**(N) Norms.** $`|\alpha|=|2n_j\!\cdot p_j|/|D_j|\le A_1|n_j|`$, likewise
$`\beta`$; $`|\hat y_\perp|\le|n_j|+(|\alpha|+|\beta|)P\le A_2|n_j|`$;
$`|d|\le8|P|\le8\sqrt2A_1|n_j|`$ because $`d`$ vanishes at the vertex
$`(0,0)`$ and has Lipschitz constant 8; and
$`|d_\beta p_{j-1}+d_\alpha p_j|\le\sqrt2P|d|`$. So
$`|V_j|\le\lambda_0|n_j|`$. For $`r\le2R`$: $`|n_j|\le2R+\Delta`$,
$`|V_{\rm in}|\le4\lambda_0(2R+\Delta)`$,
$`|V_{\rm out}|\le\kappa'(r+\tau_0)`$, hence

```math
|V|\le v,\qquad|V^2|\le v^2,\qquad m_j\,|2V_j\!\cdot n_i|\le\iota_2\quad\text{for all }i,j .
```

**(A).** If $`|r_i|\le\delta_0`$ then $`m_{i+2}=0`$, and by (L5)
$`2V_{\rm in}\!\cdot n_i\le-\sum_{j\ne i+2}m_j(\kappa|D_j|\Psi_j+2\nu\hat\rho_j)\le0`$.

**(H1).** *If $`\max_i|r_i|\le\epsilon`$ then $`\epsilon>\epsilon_*`$.* Fix a
frame $`j`$, $`D=D_j`$, $`D'=D_{j+1}`$, $`u=-D-D'=2p_{j-1}\!\cdot p_{j+1}`$.
From $`r_{j-1}=r_j+\beta D`$, $`r_{j+1}=r_j-\alpha D`$,
$`r_j=\alpha\beta D-\hat\rho`$: $`|\beta D|,|\alpha D|\le2\epsilon`$ and
$`\hat\rho\le\epsilon+4\epsilon^2/|D|`$. Decompose
$`p_{j+1}=\tfrac{D'}Dp_{j-1}+\tfrac uDp_j+q'`$; $`q'`$ is transverse, has no
$`\zeta`$-component, and $`0=p_{j+1}^2`$ gives $`-q'^2=uD'/D\ge0`$. With
$`n_{j+2}=n_j-p_j-p_{j+1}`$,

```math
r_{j+2}=r_j-\alpha D-(\alpha u+\beta D'+2\hat y_\perp\!\cdot q')+D'=r_j+D'(1+\alpha-\beta)-2\hat y_\perp\!\cdot q' ,
```

and $`|\hat y_\perp\!\cdot q'|=|y_\perp\!\cdot q'|\le\sqrt{\rho\,uD'/D}\le\sqrt{\hat\rho\,uD'/D}`$
(Cauchy–Schwarz in the negative definite transverse space). So

```math
\big|r_{j+2}-D'\big|\le\frac{4\epsilon|D'|}{|D|}+\epsilon+2\sqrt{\frac{D'u}D\Big(\epsilon+\frac{4\epsilon^2}{|D|}\Big)} .
```

The right side increases with $`\epsilon`$ and $`|D'|-\epsilon`$ decreases;
$`\epsilon_*`$ is the largest $`\epsilon`$ for which the right side is below
$`|D'|-\epsilon`$ in some frame, and then $`|r_{j+2}|>\epsilon`$, a
contradiction, for every $`\epsilon\le\epsilon_*`$. $`\square`$

*Consequence.* If $`|r_i|\le\delta_0`$, one of the frames $`i-1,i,i+1`$ has
weight 1: their weights are $`\chi(|r_{i+1}|/\delta_0)`$,
$`\chi(|r_{i+2}|/\delta_0)`$, $`\chi(|r_{i-1}|/\delta_0)`$, and if none is 1
all four $`|r|`$ are below $`2\delta_0<\epsilon_*`$.

#### (C2), (C1), (C3)

**(C2).** *If $`r\ge R`$ and $`|r_i|\le\delta_0`$ then
$`|\tau|\ge\tfrac12|\hat\xi|`$; so $`\Theta=1`$ for the sign of $`\tau`$ and
the other $`\Theta`$ is 0.*
$`(n^0)^2\ge|\hat n|^2-\delta_0`$ gives $`|n^0|\ge|\hat n|-\sqrt{\delta_0}`$
(trivially if $`|\hat n|^2<\delta_0`$). Hence
$`|\tau|\ge|n^0|-\Delta\ge|\hat\xi|-\tfrac12\xi_0`$. If $`|\hat\xi|\ge\xi_0`$
this is $`\ge\tfrac12|\hat\xi|`$. If $`|\hat\xi|<\xi_0`$ then
$`|\tau|=\sqrt{r^2-|\hat\xi|^2}>\sqrt{R^2-\xi_0^2}\ge\tfrac12\xi_0>\tfrac12|\hat\xi|`$
by the condition. $`\square`$

**(C1).** *If $`r\ge R`$, $`\tau\ge\tfrac14|\hat\xi|`$ and $`|r_i|\le\delta_0`$,
then $`n^2+u^-\!\cdot n\ge M>0`$.* $`r^2\le17\tau^2`$, so
$`n^0\ge\tau-\Delta\ge r/\sqrt{17}-\Delta>0`$, and
$`|\hat n|\le\sqrt{(n^0)^2+\delta_0}\le n^0+\sqrt{\delta_0}`$. By (C0),
$`u^-_0=\tau_0+t_i\ge\tau_0-\Delta`$ and $`|\vec u^-|=|\vec x_i|\le\Delta`$.
Because $`u^-`$ has no $`\zeta`$-component,
$`u^-\!\cdot n=u^-_0n^0-\vec x_i\!\cdot\vec n\ge(\tau_0-\Delta)n^0-\Delta|\hat n|`$,
where $`\vec n`$ is the four-dimensional space part of $`n`$ and
$`|\vec n|\le|\hat n|`$. So

```math
n^2+u^-\!\cdot n\ \ge\ -\delta_0+(\tau_0-2\Delta)n^0-\Delta\sqrt{\delta_0}\ \ge\ M .\qquad\square
```

The mirror statement holds for $`\tau\le-\tfrac14|\hat\xi|`$ with $`u^+`$:
$`n^2+u^+\!\cdot n\ge M`$ (time reversal). *Consequence, with (C2) and (A):*
if $`r>R`$ and $`|r_i|\le\delta_0`$ then
$`2V_{\rm out}\!\cdot n=-2\kappa'(n^2+u^\mp\!\cdot n)\le-2\kappa'M`$ and

```math
\operatorname{Im}X_i=(1-g)\,2V_{\rm in}\!\cdot n+g\,2V_{\rm out}\!\cdot n\ \le\ -2g\kappa'M<0 .
```

**(C3).** *If $`r\ge R`$ and $`|\tau|\le\tfrac14|\hat\xi|`$ then
$`r_i<-\delta_0`$ for all $`i`$ (and $`V_{\rm out}=0`$).*
$`r^2\le\tfrac{17}{16}|\hat\xi|^2`$, so $`s:=|\hat\xi|\ge\xi_m`$.
$`|n^0|\le\tfrac14s+\Delta`$ and $`|\hat n|\ge s-\Delta>0`$, so
$`r_i\le f(s):=(\tfrac14s+\Delta)^2-(s-\Delta)^2=\tfrac54s(2\Delta-\tfrac34s)`$.
The condition $`f(\xi_m)<-\delta_0<0`$ gives $`\xi_m>\tfrac83\Delta`$, and
$`f'(s)=\tfrac52\Delta-\tfrac{15}8s<0`$ for $`s>\tfrac43\Delta`$. So
$`f(s)\le f(\xi_m)<-\delta_0`$. $`\square`$

#### The far field away from the cones

**Lemma F.** Let $`r>R`$ and $`\vartheta:=\Theta_+>0`$ (so $`\Theta_-=0`$,
$`\tau>\tfrac14|\hat\xi|`$); put $`a=n^2=r_i`$, $`b=n\!\cdot u^-`$. Then
$`2V_{\rm out}\!\cdot n=-2\kappa'\vartheta(a+b)`$. Suppose $`a+b\le0`$. Then

1. $`|a|>\delta_0`$ and $`a<0`$;
2. $`0\le2V_{\rm out}\!\cdot n\le2\kappa'\vartheta\,c\,|a|`$, and
   $`|b|\le(c-1)|a|`$ if $`b<0`$, where $`c=1`$ if $`b\ge0`$, $`c=c_b`$ if
   $`b<0`$ and $`\vartheta=1`$, $`c=c_c`$ if $`b<0`$ and $`\vartheta<1`$;
3. $`a-V_{\rm out}^2\le-|a|\big(1-\kappa'^2(2c-1)\big)\le-\tfrac12|a|`$.

The same holds for $`\Theta_->0`$ with $`u^+`$.

*Proof.* The formula is $`V_{\rm out}=-\kappa'\vartheta(n+u^-)`$.

(1) If $`|a|\le\delta_0`$, (C1) gives $`a+b\ge M>0`$. So $`|a|>\delta_0`$.
The sign is shown with (2).

(2) *$`b\ge0`$:* $`a\le-b\le0`$, so $`a<0`$ and
$`-(a+b)\le|a|`$.

*$`b<0`$ and $`\vartheta=1`$:* $`\tau\ge\tfrac12|\hat\xi|`$, so
$`\tau\ge r/\sqrt5`$ and $`n^0\ge r/\sqrt5-\Delta\ge R/\sqrt{17}-\Delta>0`$.
$`b=u^-_0n^0-\vec x_i\!\cdot\vec n<0`$ gives
$`(\tau_0-\Delta)n^0<\Delta|\hat n|`$, so $`|\hat n|>2n^0`$ by (C0). Then
$`a\le-\tfrac34|\hat n|^2<0`$, $`|b|\le\Delta|\hat n|`$, and
$`|b|/|a|\le\tfrac{4\Delta}{3|\hat n|}\le\tfrac{2\Delta}{3(R/\sqrt{17}-\Delta)}=c_b-1`$.

*$`b<0`$ and $`\vartheta<1`$:* $`|\tau|<\tfrac12|\hat\xi|`$, so
$`r^2<\tfrac54|\hat\xi|^2`$ and $`s:=|\hat\xi|>2r/\sqrt5\ge\xi_c>4\Delta`$.
$`|n^0|\le\tfrac12s+\Delta`$, $`|\hat n|\ge s-\Delta`$, so
$`a\le(\tfrac12s+\Delta)^2-(s-\Delta)^2=-\tfrac34s(s-4\Delta)<0`$. Also
$`|n|\le r+\Delta\le\tfrac{\sqrt5}2s+\Delta`$ and
$`|u^-|\le\tau_0+\Delta`$, so

```math
\frac{|b|}{|a|}\le\Big(\frac{\sqrt5}2+\frac\Delta s\Big)\frac{\tau_0+\Delta}{\tfrac34(s-4\Delta)}\le c_c-1 ,
```

the middle expression being decreasing in $`s`$. In both cases
$`-(a+b)=|a|+|b|\le c|a|`$.

(3) $`V_{\rm out}^2=\kappa'^2\vartheta^2(a+2b+u^2)`$ with $`u^2=(u^-)^2>0`$
(C0). So $`a-V_{\rm out}^2=a(1-\kappa'^2\vartheta^2)-\kappa'^2\vartheta^2(2b+u^2)`$.
If $`2b+u^2\ge0`$ this is $`\le a(1-\kappa'^2)`$. Otherwise $`b<0`$ and
$`0<-(2b+u^2)<2|b|\le2(c-1)|a|`$, so it is
$`\le-|a|(1-\kappa'^2)+2\kappa'^2(c-1)|a|`$. Finally
$`4\kappa'c\le\tan e_{\rm far}\le\tan0.2<1`$ gives $`\kappa'\le\kappa'c\le\tfrac14`$
and $`\kappa'^2(2c-1)\le2\kappa'\cdot\kappa'c\le\tfrac18`$. $`\square`$

*Remark.* The sub-case "$`b<0`$, $`\vartheta=1`$" is empty when
$`R/\sqrt5\ge\Delta\tau_0/(\tau_0-3\Delta)`$: there
$`(\tau-\Delta)(\tau_0-\Delta)\le n^0u^-_0<\Delta|\hat n|\le\Delta(2\tau+\Delta)`$,
i.e. $`\tau<\Delta\tau_0/(\tau_0-3\Delta)`$, while $`\tau\ge R/\sqrt5`$. This
holds for the default constants and for `tuned_contour` ($`\tau_0=8\Delta`$,
$`R\ge8\Delta`$). A review showed more: the sub-case is empty for every
contour that satisfies (C0) and (C1). There $`b<0`$ forces
$`\tau(\tau_0+t_i-2|\vec x_i|)<t_i(\tau_0+t_i)-|\vec x_i|^2`$, and the largest
such $`\tau`$ over $`t_i^2+|\vec x_i|^2\le\Delta^2`$ with $`\tau_0\ge3\Delta`$
is $`1.102\Delta`$ (computed numerically), whereas $`\vartheta=1`$, $`r>R`$
and (C1) give $`\tau>R/\sqrt5>\sqrt{17/5}\,\Delta=1.844\Delta`$. This is why
the sub-case was never reached. The proof does not use it.

#### (D)

**Proposition (D).** At every real point $`\hat x`$ and for every $`i`$:
if $`\operatorname{Im}X_i>0`$ then $`|r_i|>\delta_0`$, $`\operatorname{Re}X_i\ne0`$ and

```math
\operatorname{Im}X_i\ \le\ \big[(1-g)\tan e_{\rm frames}+g\tan e_{\rm far}\big]\,|\operatorname{Re}X_i|\ \le\ \tan e\;|\operatorname{Re}X_i| .
```

Consequently, wherever $`X_i\ne0`$, its argument on the Feynman branch
lies in $`[-\pi-e,e]`$.

*Proof.* By (L5) the frames $`j\ne i+2`$ contribute $`\le0`$ to
$`\operatorname{Im}X_i=2V\!\cdot n_i`$, so

```math
\operatorname{Im}X_i\ \le\ (1-g)\,m_{i+2}\big|2V_{i+2}\!\cdot n_i\big|+g\,\big(2V_{\rm out}\!\cdot n_i\big)^+ ,\tag{9.2}
```

where the first term is absent for $`r\ge2R`$ and the second for $`r\le R`$.

*$`|r_i|\le\delta_0`$.* Then $`m_{i+2}=0`$; and if $`g>0`$,
$`2V_{\rm out}\!\cdot n_i<0`$ by (C2), (C1). So $`\operatorname{Im}X_i\le0`$.
This is the first claim.

*$`|r_i|>\delta_0`$, first term.* It is present only for $`r<2R`$, where it
is at most $`(1-g)\iota_2`$ by (N), and
$`|\operatorname{Re}X_i|\ge|r_i|-v^2>\delta_0-v^2`$. So it is at most
$`(1-g)\tan e_{\rm frames}|\operatorname{Re}X_i|`$.

*$`|r_i|>\delta_0`$, second term.* Let $`g>0`$ and
$`2V_{\rm out}\!\cdot n_i>0`$. Then one of $`\Theta_\pm`$ is positive and
Lemma F applies with $`a+b<0`$:
$`2V_{\rm out}\!\cdot n_i\le2\kappa'c|r_i|`$ with
$`c\le\max(1,c_b,c_c)`$. For $`r<2R`$:
$`|\operatorname{Re}X_i|\ge|r_i|-v^2\ge\tfrac12|r_i|`$ because
$`v^2<\delta_0/2<|r_i|/2`$. For $`r\ge2R`$: $`V=V_{\rm out}`$ and
$`|\operatorname{Re}X_i|\ge\tfrac12|r_i|`$ by Lemma F(3). In both zones the
second term of (9.2) is at most
$`g\cdot4\kappa'c|\operatorname{Re}X_i|\le g\tan e_{\rm far}|\operatorname{Re}X_i|`$.

Adding the two terms gives the inequality, with
$`\operatorname{Re}X_i\ne0`$ in each case; and
$`(1-g)\tan e_{\rm frames}+g\tan e_{\rm far}\le\max(\tan e_{\rm frames},\tan e_{\rm far})\le\tan e`$.

*The argument.* If $`\operatorname{Im}X_i\le0`$ the argument is in
$`[-\pi,0]`$. If $`\operatorname{Im}X_i>0`$ and $`\operatorname{Re}X_i>0`$ it
is in $`(0,e]`$. If $`\operatorname{Im}X_i>0`$ and
$`\operatorname{Re}X_i<0`$ the principal value is in $`[\pi-e,\pi)`$, which on
the branch continued from the Feynman prescription (values in
$`(-\tfrac32\pi,\tfrac12\pi]`$; $`e<\tfrac\pi2`$) is $`[-\pi-e,-\pi)`$.
$`\blacksquare`$

The proof gives slightly more than (D): the two excursions do not add,
the bound is the convex combination with the blend $`g`$.

#### Proposition 6(3)

**Proposition.** If $`X_i=0`$ at a real point $`\hat x`$ then $`\zeta=0`$
and $`x_0`$ lies on one of the two edge segments at $`x_i`$.

*Proof.* *$`|r_i|>\delta_0`$, $`r<2R`$:*
$`|\operatorname{Re}X_i|\ge\delta_0-v^2>0`$.

*$`|r_i|>\delta_0`$, $`r\ge2R`$:* $`V=V_{\rm out}`$. If
$`\Theta_+=\Theta_-=0`$ then $`V=0`$ and $`X_i=r_i\ne0`$. Otherwise, with
$`a,b`$ as in Lemma F: if $`a+b>0`$ then $`\operatorname{Im}X_i<0`$; if
$`a+b\le0`$ then $`\operatorname{Re}X_i\le-\tfrac12|r_i|<0`$ by Lemma F(3).

*$`|r_i|\le\delta_0`$, $`r>R`$:* $`\operatorname{Im}X_i\le-2g\kappa'M<0`$
(consequence of (C1)).

*$`|r_i|\le\delta_0`$, $`r\le R`$:* $`V=V_{\rm in}`$ and by (A)
$`\operatorname{Im}X_i\le-\sum_{j\ne i+2}m_j(\kappa|D_j|\Psi_j+2\nu\hat\rho_j)`$.
If $`\operatorname{Im}X_i=0`$, every frame with $`m_j>0`$ has $`\Psi_j=0`$
and $`\hat\rho_j=0`$, hence $`V_j=0`$ by (L5); so $`V=0`$ and $`X_i=r_i`$. By
(H1) some $`j\in\{i-1,i,i+1\}`$ has $`m_j=1`$, and for it $`\hat x`$ lies in
the flat triangle of the frame $`j`$: $`\zeta=0`$, $`y_\perp=0`$,
$`\beta\ge0`$, $`\alpha\le0`$, $`\beta-\alpha\le1`$. There
$`r_j=\alpha\beta D_j`$, $`r_{j-1}=\beta D_j(1+\alpha)`$,
$`r_{j+1}=\alpha D_j(\beta-1)`$. $`r_j=0`$ means $`\alpha=0`$ or $`\beta=0`$;
$`r_{j-1}=0`$ means $`\beta=0`$, or $`\alpha=-1`$ and then $`\beta=0`$;
$`r_{j+1}=0`$ means $`\alpha=0`$, or $`\beta=1`$ and then $`\alpha=0`$. In each
case $`x_0=x_j+\beta p_j`$ with $`\beta\in[0,1]`$ or
$`x_0=x_j+\alpha p_{j-1}`$ with $`\alpha\in[-1,0]`$: a point of an edge
segment, and since $`r_i=0`$ there, of a segment ending at $`x_i`$.
$`\blacksquare`$

#### The constants

| Statement | Holds in $`\mathbb R^{1,3+k}`$ | Constants | Where the extra coordinates enter |
| --- | --- | --- | --- |
| Lemma 5 | yes, with $`\hat\rho=\rho+m^2`$ | $`\kappa,\nu,D_j`$: the same | the component $`\nu\zeta`$ of $`V_j`$ gives $`-\nu m^2`$ |
| $`\lvert V_j\rvert\le\lambda_0\lvert n_j\rvert`$ | yes | $`A_1,A_2,\lambda_0`$: the same numbers | norms in $`\mathbb R^{4+k}`$ |
| (H1) | yes | $`\epsilon_*`$: the same number | $`\rho\le\hat\rho`$; $`q'`$ has no $`\zeta`$-component |
| (C0) | yes | $`\Delta,\tau_0`$: the same | not at all |
| (C1) | yes | $`M`$: the same number | $`u^\mp`$ has no $`\zeta`$-component; $`\lvert\vec n\rvert\le\lvert\hat n\rvert`$; $`V_{\rm out}\!\cdot n=-\kappa'(n^2+u\!\cdot n)`$ needs the component $`-\kappa'\Theta\zeta`$ |
| (C2) | yes | $`\xi_0`$: the same number | $`\hat\xi`$ and $`\hat n`$ in $`\mathbb R^{3+k}`$; $`\Theta`$ uses $`\lvert\hat\xi\rvert`$ |
| (C3) | yes | $`\xi_m`$: the same number | the same |
| (D) | yes | $`e_{\rm frames},e_{\rm far},e`$, $`c_b,c_c`$, $`v`$, $`\iota_2`$: the same numbers | through (C1), (C2), Lemma F and the norms |
| Proposition 6(3) | yes | none of its own | $`\hat\rho_j=0`$ forces $`\zeta=0`$ |

#### Numerical tests

The statements were tested with the 130-digit implementation of 9.6, in
two runs, at every point where their hypotheses hold and for all four
$`i`$. Run 1: the 875 475 points of 9.6 (seven squares with a timelike
diagonal, $`k=1,2,3`$). Run 2: 565 122 further points on six further
contours, in seven runs: the region $`(-,-)`$ with default constants
($`k=1`$ and $`k=3`$, one contour) and tuned at angles 0.3 ($`k=2`$) and 2.8
($`k=1`$); and, in each of the three regions, the default square with
$`\tau_0=3.2`$ instead of 8 (close to the limit $`3\Delta`$ of (C0), which
is about 3.10 with a timelike diagonal and 2.96 in $`(-,-)`$), which is the
only way the sub-case $`n\!\cdot u<0`$ of Lemma F was reached. The
certificate passes at all thirteen contours. Points of run 2:
generic; on the null lines; on the $`D`$-dimensional cones at all
distances; at $`r=R,2R`$ and at $`|r_i|=\delta_0,2\delta_0`$ to within
$`10^{-25}`$; in the wedge $`0.15\le|\tau|/|\hat\xi|\le0.6`$, including
$`|\tau|/|\hat\xi|=\tfrac14,\tfrac12`$ to within $`10^{-20}`$; and with the
space direction along $`\pm(x_i-c)`$.

| Statement | Tests (run 1 + run 2) | Violations | Remarks |
| --- | --- | --- | --- |
| (L5): $`-2V_j\!\cdot n_i\ge\kappa\lvert D_j\rvert\Psi_j+2\nu\hat\rho_j`$ | 17 287 164 | 0 | |
| (N): $`\lvert\alpha\rvert,\lvert\beta\rvert\le A_1\lvert n_j\rvert`$, $`\lvert V_j\rvert\le\lambda_0\lvert n_j\rvert`$; $`\lvert V\rvert\le v`$ for $`r<2R`$ | 5 762 388; 1 085 282 | 0 | |
| (A): $`\lvert r_i\rvert\le\delta_0`$: frame $`i+2`$ off, $`-2V_{\rm in}\!\cdot n_i\ge\sum m_j(\kappa\lvert D_j\rvert\Psi_j+2\nu\hat\rho_j)`$; one of the frames $`i-1,i,i+1`$ has weight 1 | 724 418 | 0 | |
| (H1): $`\max_i\lvert r_i\rvert>\epsilon_*`$ | 1 437 380 | 0 | largest $`\epsilon_*/\max_i\lvert r_i\rvert`$: 0.63 |
| (C2): $`\lvert n^0\rvert\ge\lvert\hat n\rvert-\sqrt{\delta_0}`$, $`\lvert\tau\rvert\ge\tfrac12\lvert\hat\xi\rvert`$, $`\Theta=1`$ and the other 0 | 265 255 | 0 | largest $`\lvert\hat\xi\rvert/(2\lvert\tau\rvert)`$: 0.63 |
| (C1): $`\pm n^0\ge r/\sqrt{17}-\Delta`$, $`\lvert\hat n\rvert\le\lvert n^0\rvert+\sqrt{\delta_0}`$, $`n^2+u\!\cdot n\ge M`$, $`V_{\rm out}\!\cdot n=-\kappa'(n^2+u\!\cdot n)`$, $`\operatorname{Im}X_i\le-2g\kappa'M`$ | 265 255 | 0 | smallest $`(n^2+u\!\cdot n)/M`$: 6.6 |
| (C3): $`\lvert\hat\xi\rvert\ge\xi_m`$, $`\Theta_\pm=0`$, all $`r_i<-\delta_0`$ | 68 581 | 0 | largest $`r_i/\delta_0`$: $`-162`$ |
| Lemma F: $`2V_{\rm out}\!\cdot n=-2\kappa'\vartheta(a+b)`$, the other $`\Theta`$ is 0 | 3 698 052 | 0 | identity |
| Lemma F (1), (2), $`\lvert\operatorname{Re}X_i\rvert\ge\tfrac12\lvert a\rvert`$, where $`a+b\le0`$ | 1 228 643 | 0 | largest $`2V_{\rm out}\!\cdot n/(4\kappa'c\lvert\operatorname{Re}X_i\rvert)`$: 0.50 |
| of these: $`b\ge0`$; $`b<0`$ and $`\vartheta<1`$ (with its three inequalities); $`b<0`$ and $`\vartheta=1`$ | 1 227 838; 805; 0 | 0 | the last sub-case was never reached (see the remark after Lemma F) |
| Lemma F (3) for $`r\ge2R`$ | 458 734 | 0 | |
| (D), where $`\operatorname{Im}X_i>0`$: $`\lvert r_i\rvert>\delta_0`$; (9.2); the bound with $`(1-g)\tan e_{\rm frames}+g\tan e_{\rm far}`$ | 1 318 658 (135 387 inner, 724 537 blend, 458 734 outer) | 0 | largest ratio to that bound: 0.16; largest $`\operatorname{Im}X_i/(\tan e\,\lvert\operatorname{Re}X_i\rvert)`$: 0.077 |
| (D): argument in $`[-\pi-e,e]`$ | 5 749 520 | 0 | |
| 6(3): $`X_i\ne0`$ off the edges (distance above $`10^{-40}`$) | 1 440 597 points | 0 | |
| 6(3), steps: $`\lvert\operatorname{Re}X_i\rvert\ge\delta_0-v^2`$ for $`\lvert r_i\rvert>\delta_0`$, $`r<2R`$; $`\operatorname{Im}X_i\le-\sum m_j(\dots)`$ for $`\lvert r_i\rvert\le\delta_0`$, $`r\le R`$ | 3 692 819; 459 163 | 0 | |

(H1) was also tested by minimisation: $`\max_i|x_{0i}^2-m^2|`$ was
minimised over $`(x_0,m)`$ in float64 (SLSQP, up to 400 starts) for six
of the squares of run 1 (all but $`(-,+)`$ at angle 0.5) and for the
default square in $`(-,-)`$. The minimum found with $`m`$ free equals the
four-dimensional one ($`m=0`$) to six digits in every case (0.7268, 0.7268,
0.8824, 0.01452, 0.9994, 0.002096, 0.1986) and is between 1.5 and 2.2
times $`\epsilon_*`$.

Control: with the extra components of the field removed, on 60 000 points
each of two squares, the identity $`V_{\rm out}\!\cdot n=-\kappa'(n^2+u\!\cdot n)`$
of (C1) fails at 79% of its points, the identity of Lemma F at 89 to 90%,
(A) at 11 to 27%, and Lemma F(2) at 32 to 37%.

*Not checked.* The sub-case "$`b<0`$, $`\vartheta=1`$" of Lemma F was not
reached by any point; it rests on its three-line proof. The other
sub-case with $`b<0`$ was reached only on the three contours with
$`\tau_0=3.2`$ (805 points). The minimisation for (H1) is a local method
with random starts, not a certified global minimum. Lipschitz bounds, (E),
(F), (G), (S), (I) and (Q) of GATE_6_CONTOUR section 9 are not treated
here; (G) (Lemma 6) is in 9.2 (c).
