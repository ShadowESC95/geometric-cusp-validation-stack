# The phase bound with a timelike diagonal

[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md) proves $`|\arg z|\le2\pi+4e`$ on the
contour when both diagonals are spacelike, and $`|\arg z|\le3\pi-0.3`$ when
one is timelike (Lemma Q). It leaves open a bound close to $`2\pi`$ in the
second case. This document proves

```math
|\arg z|\ \le\ 2\pi+e+\max\big(K_{\rm in},K_{\rm bl},e,3e_{\rm sp}\big)
```

with explicit constants, under explicit inequalities between the constants
of the contour (section 6). For the default constants the right side is
$`2\pi+0.36`$. The constants are not sharp and not uniform in the square;
see the scope remarks in section 6.

**Status.** Proved, under the conditions listed in section 6. The
argument was found by an automated search and written out; a separate
automated review then re-derived every step, tested about 140
inequalities per point with its own 90-digit implementation of the field
on eight squares (no violation), and found three small gaps and one wrong
constant, none affecting the result. Those are repaired in this text, and
the repaired passages (the definition of $`K_{\rm in}`$, (B3), (F1), and
cases N1c, N2 and N3 of section 5) have not been reviewed again. Nothing
here has been reviewed by a specialist.

**What is taken from GATE_6_CONTOUR without re-derivation:** (D)
($`\arg x_{0i}^2\in[-\pi-e,e]`$ on the contour); Lemma 5; Lemma 6
($`|V^2|\le\tfrac\theta2|\operatorname{Im}x_{0i}^2|`$ for $`r\le R`$,
$`|x_{0i}^2|\le\delta_0`$); (H1); (C1) with its margin $`M`$; (C2);
Proposition 6(3); and the inequality $`|V_{\rm in}|^2\le3\mu I_{\rm in}`$ from
the proof of Lemma Q.

**Not covered:** the contour in $`D`$ dimensions; interval arithmetic (the
conditions are evaluated in floating point at individual squares); the
region $`(-,-)`$, where $`2\pi+4e`$ is already proved.

Code: `src/geocusp/phase_bound.py`. Tests: `tests/test_phase_bound.py`.

## 0. Notation and conventions

Constants and conditions are those of GATE_6_CONTOUR section 9:
$`\Delta,P,A_1,A_2,\lambda_0,B,\theta,\mu,e`$, the conditions (H1), (C0)–(C3),
(D), (G). Norms $`|\cdot|`$ of vectors are Euclidean, and
$`|a\!\cdot b|\le|a||b|`$ for the Minkowski product.

Let $`c`$ be a cusp whose diagonal is timelike: $`D:=D_c=2p_{c-1}\!\cdot p_c>0`$.
Then $`D':=D_{c\pm1}<0`$ is the other diagonal and $`u=-D-D'<0`$. The two
edges $`p_{c-1},p_c`$ at $`x_c`$ have the same time orientation, because
their product is positive. Indices are taken mod 4.

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
$`\xi_-`$ future causal. Then $`V^2=\kappa'^2b^2\xi_-^2\ge0`$.

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
\beta_{\max}:=A_1(R+\Delta),\qquad\alpha'_{\max}:=\beta_{\max}\frac{D}{|D'|},\qquad\beta_*:=\frac{|D'|}{|u|}\,\gamma_1,
```
```math
c_1:=(2\nu+2\kappa\alpha'_{\max})\,\frac{2D}{\delta_0}\max(1,\beta_{\max}),\qquad
c_1':=\frac{2|D'|\big[\kappa\alpha_{\max}'^2+(1+\beta_{\max})(\nu\alpha'_{\max}+\kappa\beta_{\max})\big]}{\tfrac12\delta_0\min(1,\beta_*)},\qquad
c_*:=\max(c_1,c_1').
```

For $`r\le R`$: every frame has $`|V_j|\le\lambda_0(R+\Delta)`$, so
$`|J_i^{(j)}|:=|2V_j\!\cdot n_i|\le\iota`$ for all $`i,j`$, and
$`|V^2|\le|V|^2\le v_R^2`$. In every frame the coordinates satisfy
$`|\alpha|,|\beta|\le\beta_{\max}`$ (section 9).

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
K_c\ \le\ K_{\rm in}:=\max\Big(2e+e_4,\ \ e+2\varepsilon_V,\ \ \arctan\!\Big(2\nu+2c_*+\frac{2\iota}{\delta_0}\Big)+2\arctan\frac{\iota}{\delta_0}+3\varepsilon_V\Big),
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
$`\alpha'=A/D'=bD/|D'|>0`$, $`\alpha'|D'|=|A|\le|r_{c-1}|`$, and
$`\rho'=\alpha'\beta'D'-r_{c-1}=|r_{c-1}|-\alpha'\beta'|D'|`$. All
invariants being negative, Lemma R3 (as in 3a) gives

```math
J_c^{(c-1)}|r_{c-1}|-J_{c-1}^{(c-1)}|r_c|=r_cJ_{c-1}^{(c-1)}-r_{c-1}J_c^{(c-1)}
=|D'|\Big[\kappa|D'|\alpha'^2d'_\alpha+\rho'\big(2\nu\alpha'-\kappa d'_\beta\big)\Big]=:E' ,
```

with $`d'_\alpha\le\phi'(\alpha')\le2\alpha'`$ and
$`-d'_\beta\le\phi'((-\beta')^+)\le2(-\beta')^+`$.

If $`\beta'\ge0`$: $`\rho'\le|r_{c-1}|`$ and $`-d'_\beta\le0`$, so
$`E'\le|D'|\big[\kappa|D'|\alpha'^2\phi'(\alpha')+2\nu\alpha'|r_{c-1}|\big]\le|A|\,|r_{c-1}|\,\big(\kappa\phi'(\alpha')+2\nu\big)`$.
Also $`|A|/|r_c|\le\tfrac{2D}{\delta_0}\max(1,b)`$: if
$`|r_c|\le\delta_0/2`$ then $`aD=|r_{c+1}|-|r_c|>\delta_0/2`$ and
$`|A|/|r_c|\le bD/(abD)=1/a<2D/\delta_0`$; if $`|r_c|>\delta_0/2`$ then
$`|A|/|r_c|<2bD/\delta_0`$. With $`\alpha'\le\alpha'_{\max}`$,
$`b\le\beta_{\max}`$:

```math
\frac{J_c^{(c-1)}}{|r_c|}-\frac{J_{c-1}^{(c-1)}}{|r_{c-1}|}=\frac{E'}{|r_c||r_{c-1}|}\le c_1 .
```

If $`\beta'<0`$: by Lemma R5.1,
$`0>\beta'\ge1+a-b\tfrac{|u|}{|D'|}-2\sqrt{\rho|u|/(D|D'|)}`$ and
$`\rho<\delta_0`$, so $`b>\tfrac{|D'|}{|u|}\gamma_1=\beta_*`$. Then
$`|r_c|\ge\tfrac12\delta_0\min(1,\beta_*)`$: if $`|r_c|\le\delta_0/2`$ then
$`a>\delta_0/(2D)`$ and $`|r_c|\ge abD>\tfrac12\delta_0\beta_*`$. Further
$`\rho'=|r_{c-1}|+\alpha'|\beta'||D'|\le|r_{c-1}|(1+|\beta'|)`$ and
$`\kappa|D'|^2\alpha'^2\cdot2\alpha'=2\kappa\alpha'^2|D'|\,|A|\le2\kappa\alpha'^2|D'||r_{c-1}|`$,
so $`E'\le|r_{c-1}||D'|\big[2\kappa\alpha'^2+(1+|\beta'|)(2\nu\alpha'+2\kappa|\beta'|)\big]`$,
and with $`|\beta'|\le\beta_{\max}`$:

```math
\frac{E'}{|r_c||r_{c-1}|}\le c_1' .
```

In both cases

```math
\frac{J_c^{(c-1)}}{|r_c|}\le\frac{J_{c-1}^{(c-1)}}{|r_{c-1}|}+c_*. \tag{3d.2}
```

*Frame $`c+1`$* (on only if $`|r_{c-1}|>\delta_0`$). Here
$`\beta''=-B/D'=-aD/|D'|<0`$, $`|\beta''||D'|=|B|\le|r_{c+1}|`$, and
$`\rho''=\alpha''\beta''D'-r_{c+1}=|r_{c+1}|+\alpha''|\beta''||D'|`$. Lemma R3
(as in 3b) gives

```math
J_c^{(c+1)}|r_{c+1}|-J_{c+1}^{(c+1)}|r_c|=r_cJ_{c+1}^{(c+1)}-r_{c+1}J_c^{(c+1)}
=|D'|\Big[-\kappa|D'|\beta''^2d''_\beta+\rho''\big(2\nu|\beta''|+\kappa d''_\alpha\big)\Big]=:E'' ,
```

with $`-d''_\beta\le\phi'(|\beta''|)\le2|\beta''|`$ and
$`d''_\alpha\le\phi'((\alpha'')^+)\le2(\alpha'')^+`$.

If $`\alpha''\le0`$: $`\rho''\le|r_{c+1}|`$ and $`d''_\alpha\le0`$, so
$`E''\le|B|\,|r_{c+1}|\,\big(\kappa\phi'(|\beta''|)+2\nu\big)`$. And
$`|B|/|r_c|\le\tfrac{2D}{\delta_0}\max(1,a)`$: if $`|r_c|\le\delta_0/2`$ then
$`bD=|r_{c-1}|-|r_c|>\delta_0/2`$ and $`|B|/|r_c|\le1/b<2D/\delta_0`$;
otherwise $`|B|/|r_c|<2aD/\delta_0`$. With
$`|\beta''|\le\alpha'_{\max}`$, $`a\le\beta_{\max}`$:
$`E''/(|r_c||r_{c+1}|)\le c_1`$.

If $`\alpha''>0`$: by Lemma R5.1,
$`0<\alpha''\le-1-b+a\tfrac{|u|}{|D'|}+2\sqrt{\rho|u|/(D|D'|)}`$, so
$`a>\beta_*`$, and as before $`|r_c|\ge\tfrac12\delta_0\min(1,\beta_*)`$ (if
$`|r_c|\le\delta_0/2`$ then $`b>\delta_0/(2D)`$ and
$`|r_c|\ge abD>\tfrac12\delta_0\beta_*`$). Further
$`\rho''\le|r_{c+1}|(1+\alpha'')`$ and
$`E''\le|r_{c+1}||D'|\big[2\kappa\beta''^2+(1+\alpha'')(2\nu|\beta''|+2\kappa\alpha'')\big]`$,
so $`E''/(|r_c||r_{c+1}|)\le c_1'`$. In both cases

```math
\frac{J_c^{(c+1)}}{|r_c|}\le\frac{J_{c+1}^{(c+1)}}{|r_{c+1}|}+c_*. \tag{3d.3}
```

*Sum.* $`J_i=\sum_{j\in\{c-1,c,c+1\}}m_jJ_i^{(j)}`$. Multiply (3d.1),
(3d.2), (3d.3) by $`m_c,m_{c-1},m_{c+1}\in[0,1]`$ and add:

```math
\frac{J_c}{|r_c|}\ \le\ \frac{J_{c-1}-m_{c+1}J_{c-1}^{(c+1)}}{|r_{c-1}|}+\frac{J_{c+1}-m_{c-1}J_{c+1}^{(c-1)}}{|r_{c+1}|}+2\nu+2c_* .
```

The two subtracted terms are the contributions of the frames that do not
control that invariant. They are present only where the invariant exceeds
$`\delta_0`$ in modulus, and are at most $`\iota`$ in modulus. Therefore

```math
\frac{J_c}{|r_c|}\ \le\ T_{c-1}+T_{c+1}+2\nu+2c_*+\frac{2\iota}{\delta_0},\qquad T_i:=\frac{J_i}{|r_i|}\ \ge-\frac{\iota}{\delta_0}\quad(i=c\pm1). \tag{3d.4}
```

(The lower bound on $`T_i`$: the controlling frames contribute $`\ge0`$ to
$`J_i`$.) Now $`\dot u_i=\arctan T_i`$ for all three indices, because the
real parts are negative. By (3d.4), monotonicity and subadditivity of
$`\arctan`$ on $`[0,\infty)`$,

```math
\dot u_c\le\arctan(T_{c-1}^+)+\arctan(T_{c+1}^+)+\arctan\!\Big(2\nu+2c_*+\frac{2\iota}{\delta_0}\Big),\qquad
\arctan(T_i^+)\le\dot u_i+\arctan\frac{\iota}{\delta_0} .
```

With $`|u_i-\dot u_i|\le\varepsilon_V`$ this is the third entry of
$`K_{\rm in}`$. All cases are covered, since the four sign combinations of
$`(\alpha,\beta)`$ are 3a–3d. $`\blacksquare`$

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
$`M`$ the margin of (C1), as in Lemma Q:

```math
K:=\kappa'(2R+\tau_0+\Delta),\quad v_{\rm in}:=4\lambda_0(2R+\Delta),\quad v:=v_{\rm in}+K,\quad
\iota_2:=2\lambda_0(2R+\Delta)^2,\quad\pi_0:=2\lambda_0(2R+\Delta)P,
```
```math
\mu_-:=3\mu+K\sqrt{\frac{3\mu}{2\kappa'M}},\qquad\mu_+:=\mu_-+\frac{K^2}{2\kappa'M},\qquad\omega_0:=v_{\rm in}(2K+v_{\rm in}),\qquad\delta_v:=\delta_0-v^2,
```
```math
\varepsilon:=\frac{\sqrt2\,P}{R-\Delta},\qquad b_0:=\frac{\varepsilon\,\delta_0}{\min_j|D_j|},\qquad\Lambda_b:=\frac{R-\Delta-B\sqrt{\delta_0+2b_0\max_j|D_j|}}{P},\qquad a_{\min}:=\Lambda_b-b_0,\qquad
\delta_1:=\delta_0(1-\varepsilon),
```
```math
\varsigma:=\Big(2\nu+\frac{\kappa b_0^2}{a_{\min}(1+a_{\min})}\Big)+\frac{2\nu+\kappa b_0^2}{a_{\min}},\qquad F:=\frac{2KP\mu_-}{\delta_v}.
```

For $`r<2R`$: $`|V_j|\le\lambda_0(2R+\Delta)`$, so
$`|J^{(j)}_i|\le\iota_2`$, $`|2V_j\!\cdot p|\le\pi_0`$ for an edge $`p`$,
$`|V_{\rm in}|\le v_{\rm in}`$, $`|V_{\rm out}|\le K`$, $`|V^2|\le v^2`$.

**Conditions (B).** Besides section 9 of GATE_6_CONTOUR ((C0), (C1), (C2),
(D), in particular $`v^2<\delta_0/2`$, which also gives
$`\omega_0<\delta_0`$):

- (B1) $`R\ge\tau_0`$ and $`R>\sqrt5\,\Delta`$;
- (B2) $`\kappa\le\nu`$;
- (B3) $`b_0\le\tfrac12`$, $`\Lambda_b\ge2`$, $`4\delta_0\le\min_j|D_j|`$;
- (B4) $`\delta_1>v^2`$ (in particular $`\varepsilon<1`$);
- (B5) $`s<1`$ and $`e+\arcsin s+\arctan\mu_-<\tfrac\pi2`$, with $`s`$ defined
  in the statement;
- and $`v^2<(\tau_0-2\Delta)^2`$.

**Proposition R6.** Assume the conditions of section 9 and (B). Then for
$`R<r<2R`$

```math
K_c\le K_{\rm bl}:=\max\big(e+\arcsin s,\ 3e_{\rm sp}'\big),\qquad
\tan e'_{\rm sp}=\frac{\dfrac{4\iota_2}{(\tau_0-2\Delta)^2}+2\kappa'\Big(1+\sqrt2\dfrac{\tau_0+\Delta}{\tau_0-2\Delta}\Big)}{1-v^2/(\tau_0-2\Delta)^2},
```

where $`s=\max(s_0,s_1,s_2,s_3)`$,

```math
s_0=\frac{8\iota_2}{\delta_v}+\frac{2KP\,\omega_0}{\delta_v^2},\qquad
s_2=\frac{\iota_2+2\mu_+\pi_0}{\delta_v}+F,\qquad s_3=\mu_-+(1+\mu_+)^2\varsigma,
```
```math
s_1=\max\Big(\frac{\iota_2+2\mu_+\pi_0}{\delta_v}+F,\ \ (1+\mu_+)\frac{\delta_0}{\delta_v}\,\varsigma+\frac{2\mu_-\pi_0}{\delta_v}+\iota_2\Big(\frac1{\delta_v}+\frac1{\delta_1-v^2}\Big)+F\Big).
```

*Proof.* Write $`x_0-c=(\tau,\vec\xi)`$.

**Step 0: three regions.** (a) If $`x_0\notin J^+(c_-)\cup J^-(c_+)`$ then,
as in R4(iii), $`-r_i\ge(\tau_0-2\Delta)|\vec n_i|\ge(\tau_0-2\Delta)^2`$ and
$`|n_i|\le\sqrt2|\vec n_i|`$ for all $`i`$. The frames contribute at most
$`4\iota_2`$ to $`|\operatorname{Im}X_i|`$, the far field at most
$`2\kappa'|r_i|(1+\sqrt2\tfrac{\tau_0+\Delta}{\tau_0-2\Delta})`$, and
$`|V^2|\le v^2`$. So $`|u_i|\le e'_{\rm sp}`$ and $`K_c\le3e'_{\rm sp}`$.

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
$`|n|\le(|a_1|+|a_2|)P+B\sqrt\rho`$. For $`n=n_i`$: $`|n_i|\ge R-\Delta`$.

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

**Case N0: neither $`c`$ nor $`c-1`$ is near.** Far bracket:
$`\ge-\omega_0\cdot2KP\ge-\tfrac{2KP\omega_0}{\delta_v^2}|X_c||X_{c-1}|`$.
Frames: by the crude bound for the four frames,
$`\ge-4\iota_2(|R_c|+|R_{c-1}|)\ge-\tfrac{8\iota_2}{\delta_v}|X_c||X_{c-1}|`$.
So (5.2) holds with $`s_0`$.

**Case N1: $`c`$ near, $`c-1`$ not.** The frame $`c+2`$ is off; $`c-1`$ and
$`c`$ control both invariants, $`c+1`$ controls $`c`$ only. By (F2),
$`\alpha,\beta\ge-b_0`$. Far bracket:
$`\ge-V^2_-\Delta^W\ge-\mu_-J_c\cdot2KP\ge-F|X_c||X_{c-1}|`$.

*N1a: $`r_{c-1}>\delta_0`$ and $`R_c\le0`$.* Then $`R_{c-1}>0\ge R_c`$. For
$`j=c-1,c`$ both $`J^{(j)}`$ are $`\ge0`$, so the bracket is $`\ge0`$. For
$`j=c+1`$: $`R_{c-1}J^{(c+1)}_c\ge0`$ and
$`-R_cJ^{(c+1)}_{c-1}\ge-|R_c|\iota_2\ge-\tfrac{\iota_2}{\delta_v}|X_c||X_{c-1}|`$.

*N1b: $`r_{c-1}>\delta_0`$ and $`R_c>0`$.* If $`r_c>0`$: by (F1) $`y`$ is
future timelike, so $`\alpha,\beta>0`$ and $`\rho<\alpha\beta D`$. Then
$`\ell_3<\beta`$, so $`\kappa\phi'(\ell_3)\le2\kappa\beta\le2\nu\beta`$ by
(B2); and by Lemma R5.1 with $`\rho\le\alpha\beta D`$,
$`\beta'\ge1+\big(\alpha|D'|+\beta|u|-2\sqrt{\alpha|D'|\beta|u|}\big)/|D'|\ge1`$.
So (E+) gives $`E_c,E_{c-1}\ge0`$. If $`r_c\le0`$: $`r_{c-1}>0\ge r_c`$ and
$`E_c,E_{c-1}\ge0`$ because the $`J`$'s are $`\ge0`$. In both sub-cases
$`|V^2\Delta^{(j)}|\le\mu_+J_c\pi_0\le\tfrac{\mu_+\pi_0}{\delta_v}|X_c||X_{c-1}|`$
for $`j=c-1,c`$. For $`j=c+1`$: $`R_{c-1}J^{(c+1)}_c\ge0`$ and
$`-R_cJ^{(c+1)}_{c-1}\ge-R_c\iota_2\ge-\tfrac{\iota_2}{\delta_v}|X_c||X_{c-1}|`$.

*N1c: $`r_{c-1}<-\delta_0`$.* Then $`A=r_{c-1}-r_c<0`$, so $`\beta<0`$. By
(F2) $`|\beta|D\le\varepsilon\delta_0`$, i.e.
$`|A|\le\varepsilon\delta_0`$. Hence $`r_c\le-\delta_1`$ and
$`|r_{c-1}|\le\delta_0(1+\varepsilon)`$.

*$`\alpha>0`$ and $`\alpha\ge a_{\min}`$.* If $`\alpha\le0`$ then
$`\rho=\alpha\beta D-r_c\le b_0^2D+\delta_0`$ and (F3) gives
$`R-\Delta\le2b_0P+B\sqrt{\delta_0+b_0^2D}<2P+B\sqrt{\delta_0+2b_0\max|D|}`$,
contradicting $`\Lambda_b\ge2`$. So $`\alpha>0`$,
$`\rho=|r_c|-\alpha|\beta|D\le\delta_0`$, and (F3) gives
$`\alpha\ge\Lambda_b-b_0=a_{\min}`$.

*$`\beta'\ge0`$.* Here $`c-1`$ is not near, so (F2) is applied with a larger
$`\delta`$. First $`n^0_{c-1}=n_c^0+p^0_{c-1}>0`$, because $`c`$ is near and
$`p_{c-1}`$ is future directed. With $`\delta=\delta_0(1+\varepsilon)`$ and
$`p_{c-2}`$ past directed, $`D'<0`$:
$`\beta'=2n_{c-1}\!\cdot p_{c-2}/D'\ge-b_0(1+\varepsilon)=:-b_0'`$, and
$`b_0'<2b_0\le1`$ by (B4), (B3). Also
$`\alpha'=A/D'=|A|/|D'|\in(0,b_0]`$. Suppose $`\beta'<0`$. Then

```math
\rho'=\alpha'\beta'D'-r_{c-1}=\alpha'|\beta'||D'|+|r_{c-1}|\le b_0b_0'|D'|+\delta_0+b_0\min|D|\le\delta_0+2b_0\max|D| ,
```

and (F3) gives
$`R-\Delta\le(b_0+b_0')P+B\sqrt{\delta_0+2b_0\max|D|}<2P+B\sqrt{\delta_0+2b_0\max|D|}`$,
contradicting $`\Lambda_b\ge2`$. So $`\beta'\ge0`$.

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
for $`j=c-1,c`$. For $`j=c+1`$ the crude bound with
$`|R_c|\ge\delta_1-v^2`$, $`|R_{c-1}|\ge\delta_v`$ gives
$`\ge-\iota_2\big(\tfrac1{\delta_v}+\tfrac1{\delta_1-v^2}\big)|X_c||X_{c-1}|`$.

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
$`R_{c-1}J_c^{(c+2)}\ge-|R_{c-1}|\iota_2\ge-\tfrac{\iota_2}{\delta_v}|X_c||X_{c-1}|`$.

For $`j=c-1,c`$: if $`R_{c-1}\ge0`$ the brackets are $`\ge0`$ (signs). Let
$`R_{c-1}<0`$. If $`r_{c-1}\ge0`$ then $`E_j\ge0`$ by the signs. If
$`r_{c-1}<0`$ I verify (E+):

- $`\kappa\phi'(\ell_3)\le2\nu\beta`$. Here $`\ell_3=\beta-(1+\alpha)\le\beta+b_0`$.
  If $`1+\alpha\ge0`$: $`\phi'(\ell_3^+)\le2\ell_3^+\le2\beta`$. If
  $`-b_0\le1+\alpha<0`$: write
  $`n_{c-1}=(1+\alpha)p_{c-1}+\beta p_c+y_\perp`$, with
  $`\rho=(1+\alpha)\beta D-r_{c-1}\le\delta_0`$; (F3) applied to this
  decomposition gives $`R-\Delta\le(b_0+\beta)P+B\sqrt{\delta_0}`$, so
  $`\beta\ge\Lambda_b-b_0\ge1`$. Then $`\ell_3\ge1`$ and
  $`\phi'(\ell_3)=2\ell_3-1\le2\beta+2b_0-1\le2\beta`$ by (B3). Then use (B2).
- If $`\beta'<0`$: $`\rho'=\alpha'\beta'D'-r_{c-1}\le\delta_0`$ (because
  $`\alpha'\beta'\ge0`$, $`D'<0`$), and (F3) in the frame of $`c-1`$ gives
  $`|\alpha'|\ge\Lambda_b-b_0\ge1`$; so
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
N2, $`|\alpha'|\ge\Lambda_b-b_0>1`$, a contradiction; so $`\beta'\ge0`$,
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
$`u_c\le u_{c-1}+\arcsin s`$, and with $`u_{c+1}\ge-e`$:
$`K_c\le e+\arcsin s`$. $`\blacksquare`$

**Remarks.**

1. The constant is set by the crude Case N0: $`s\approx8\iota_2/\delta_v`$.
   For the default constants $`s_0=0.209`$, $`s_1=0.067`$, $`s_2=0.030`$,
   $`s_3=0.0006`$. The measured minimum of $`Q/(|X_c||X_{c-1}|)`$ in the
   blend zone is $`-7\times10^{-7}`$ and the measured maximum of $`K_c`$ is
   $`0.0027`$. A sharper Case N0 would need the closed forms of the frames
   $`c,c-1`$ at points where both invariants are spacelike and beyond
   $`\delta_0`$.
2. (B2) is used only in N1b and N2. It holds for the default constants and
   with equality for `tuned_contour`.

## 6. The whole contour

**Theorem.** Let one diagonal be timelike. Assume the conditions of
GATE_6_CONTOUR section 9 and

- (R4) $`\tau_0\ge3\Delta`$, $`2R>\sqrt2\,\tau_0`$;
- (G1) $`4\delta_0|u|<D|D'|`$; (G2) $`4\delta_0<D`$;
- (B1) $`R\ge\tau_0`$, $`R>\sqrt5\Delta`$; (B2) $`\kappa\le\nu`$;
  (B3) $`b_0\le\tfrac12`$, $`\Lambda_b\ge2`$, $`4\delta_0\le\min_j|D_j|`$;
  (B4) $`\delta_1>v^2`$; (B5) $`s<1`$, $`e+\arcsin s+\arctan\mu_-<\tfrac\pi2`$;
- the trivial ones: $`e<\tfrac\pi2`$, $`\epsilon_\kappa<1`$,
  $`v_R^2<\delta_0`$, $`\theta\le2`$, $`\nu<1`$, $`4\kappa^2<1`$,
  $`v^2<\delta_0/2`$, $`v^2<(\tau_0-2\Delta)^2`$.

Then on the whole contour, off the four edges,

```math
|\arg z|\ \le\ 2\pi+e+\max\big(K_{\rm in},\ K_{\rm bl},\ e,\ 3e_{\rm sp}\big).
```

*Proof.* Lemma R1 with Propositions R5 ($`r\le R`$), R6 ($`R<r<2R`$) and R4
($`r\ge2R`$), applied at one cusp with a timelike diagonal. The other
one-sided bound, $`\arg z\ge-\pi-4e`$ for $`s>0`$ and $`\arg z\le\pi+4e`$ for
$`t>0`$, is Lemma R1. $`\blacksquare`$

All the new conditions are inequalities between the constants that the
certificate already computes. (G2) is implied by (B3).

**Values.** `phase_bound.phase_bound(contour)` evaluates every constant and
condition. All conditions hold at the squares below.

| Square | $`e`$ | $`K_{\rm in}`$ | $`K_{\rm bl}`$ | $`\max(e,3e_{\rm sp})`$ | bound $`-2\pi`$ |
| --- | --- | --- | --- | --- | --- |
| $`(+,-)`$, angle 1.1, default constants | 0.074 | 0.174 | 0.285 | 0.074 | 0.359 |
| $`(-,+)`$, angle 1.1, default constants | 0.074 | 0.174 | 0.285 | 0.074 | 0.359 |
| $`(+,-)`$, angle 0.3, tuned | 0.089 | 0.227 | 0.443 | 0.089 | 0.533 |
| $`(-,+)`$, angle 0.5, tuned | 0.090 | 0.227 | 0.444 | 0.090 | 0.533 |
| $`(+,-)`$, angle 0.7, tuned | 0.090 | 0.227 | 0.445 | 0.090 | 0.534 |
| $`(-,+)`$, angle 1.1, tuned | 0.090 | 0.228 | 0.447 | 0.090 | 0.537 |
| $`(+,-)`$, angle 2.0, tuned | 0.091 | 0.363 | 0.463 | 0.091 | 0.554 |
| $`(+,-)`$, angle 2.5, tuned | 0.082 | 1.495 | 0.485 | 0.082 | 1.577 |

**Scope.**

- The bound is not uniform in the square. $`K_{\rm in}`$ is $`2e+e_4`$ for
  most squares, but its entry from case 3d contains
  $`c_1'\sim1/(\delta_0\min(1,\beta_*))`$, which grows when $`|D'|`$ is small;
  at angle 2.5 it gives $`2\pi+1.58`$. Every value in the table is below
  $`3\pi-0.3`$, the bound of Lemma Q, and the phase observed on the contour
  never exceeds $`2\pi`$ by more than rounding. The constants are not sharp:
  the observed maxima of $`K_c`$ are about $`10^{-6}`$ in the inner zone and
  $`2\kappa'`$ in the other two.
- The conditions are evaluated in floating point at individual squares,
  not in interval arithmetic over a range of squares.
- The contour in $`D`$ dimensions (GATE_6_CONTOUR section 11) is not
  treated.
- The region $`(-,-)`$ is not treated because nothing is needed there:
  $`|\arg z|\le2\pi+4e`$ follows from (D) alone.

## 7. Numerical checks

The field was re-implemented in 50 to 70 digit arithmetic with exactly null
edges and compared with the float64 field (relative difference at most
$`3\times10^{-13}`$). On three squares ($`(+,-)`$ at angles 1.1 and 0.7,
$`(-,+)`$ at 1.1), with points aimed at the null lines through adjacent
cusps on and beyond the segments, at offsets from $`10^{-30}`$ to 1, at the
cones, and at generic points in all three zones:

| Statement | Points per square | Violations |
| --- | --- | --- |
| Lemma R2(3) and Lemma R3 (identities) | 16 000 | 0; agreement to $`10^{-37}`$ |
| Signs of the closed forms in cases 3a, 3b, 3c; (3d.1) | 16 000 | 0 |
| Lemma R5.1(a), (b) | 15 900 | 0 |
| (3d.2), (3d.3) against $`c_1`$, $`c_1'`$ | 2 500 to 3 400 | 0; measured at most 0.0012 |
| Proposition R5′: field is $`V_c`$, $`K_c\le`$ bound | 15 700 | 0; maximum equals the bound to five digits |
| Proposition R4: all statements of the proof | 32 000 | 0 |
| Proposition R6: (F1), (F2), (F4), the signs of $`E_c,E_{c-1}`$, (5.2) by case | 9 200 | 0 |
| $`K_c`$ and $`\lvert\arg z\rvert`$ over all zones | 53 000 | largest $`K_c`$: $`5\times10^{-6}`$, 0.0027, 0.0039; largest $`\lvert\arg z\rvert/\pi`$: 1.9999999976 |

An independent re-derivation of the constants and an independent
re-implementation of the checks agreed with these.

In float64 the real parts of the invariants are lost within about
$`10^{-15}`$ of a null line through two adjacent cusps, and a float64 survey
shows spurious values of $`K_c`$ up to 0.06 there. The surveys in
`tests/test_phase_bound.py` therefore skip points where an invariant is
within rounding error of zero.

```bash
PYTHONPATH=src pytest -q tests/test_phase_bound.py
```
