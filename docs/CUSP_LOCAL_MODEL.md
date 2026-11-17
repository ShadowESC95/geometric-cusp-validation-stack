# The cusp-local model and a proof of the double pole without Mellin–Barnes

The question: does the local structure at a cusp determine the double pole
of the whole integral directly?

**Answer.** Largely, with one step still missing. Proved here, for the
weight $`z^p`$ and then for the class of weights of GATE_6_BOUNDS Lemma 3
(items 1 to 3 for $`|\operatorname{Re}p|<1`$, item 4 for
$`\eta-1<\operatorname{Re}p<1+\eta`$):

1. The cusp-local model, in which the far propagator is frozen at its cusp
   value, is a triangle with two light-like sides. It has a closed form. Its
   double pole is the product of the collinear poles of the two edges at the
   cusp, with coefficient $`\sin\pi p/(\pi p)=M[z^p]`$ (section 2).
2. The full parametric integral equals the sum of the four cusp-local
   models up to a single pole, by elementary estimates on the Feynman
   simplex, without Mellin–Barnes (section 3). This now covers a timelike
   diagonal as well (section 6).
3. For a general weight the regulated cusp-local model, defined through
   its parametric form, is given exactly by an integral of $`F`$ along the
   path between the two lips of $`z=-1`$, smeared radially by an explicit
   probability density of width $`\eta`$ (section 7). So for the model
   $`\eta^2\times`$(model)$`\to M[F]`$ with an $`O(\eta)`$ correction, for
   every weight polynomially bounded in $`\log z`$ on the principal sheet.

4. The cusp-local model is derived directly from the real Lorentzian
   domain (section 8). Closing the contour in one light-cone coordinate
   shows that the integral is supported on the edge segment and equals a
   discontinuity across the negative $`z`$ axis. Three elementary integrals
   remain, one for the position along the edge, one for the collinear
   scale and one along the path between the lips. They factorise exactly,
   so the region where both collinear scales vanish is under complete
   control.

What is **not** proved: the same for the full four-point integral. There
the far propagator is not frozen, and the statement that the full integral
differs from the four cusp-local models by a single pole is proved only for
the parametric integral (sections 3 and 6).

Code: `src/geocusp/cusp_local.py`. Tests: `tests/test_cusp_local.py`.

## 1. Exact variables at a cusp

Set $`|s|=1`$. At the cusp $`x_2`$ (GATE_6_CLOSURE section 2),
$`x_{02}^2=w`$, $`x_{01}^2=u_1=w+q_1`$, $`x_{03}^2=u_3=w+q_3`$,
$`w=q_1q_3-\rho`$. With the far propagator frozen,
$`z=-w/(u_1u_3)`$. Use $`(u_1,u_3,z)`$ as coordinates. Then, exactly,

```math
\frac{dq_1\,dq_3\,d\rho}{w\,u_1u_3}=-\frac{du_1}{u_1}\,\frac{du_3}{u_3}\,\frac{dz}{z},
\qquad
\rho=u_1u_3\,\big[1+z+z(u_1+u_3)+z^2u_1u_3\big].
```

So the measure of the cusp-local model is a product of three logarithmic
forms, not only in the soft limit, and the regulator is
$`\rho^\eta=|u_1u_3|^\eta\,|1+z+\dots|^\eta`$. As $`u_1,u_3\to0`$ the
bracket tends to $`1+z`$: the $`z`$ integral starts at $`z=-1`$. This is
the exact form of the statement that the two logarithms and the contour in
$`z`$ separate. (`measure_identity`, symbolic.)

## 2. The cusp-local model is a light-like triangle

Freezing $`x_{04}^2=t`$ in $`\prod(-x_{0i}^2)^{-a_i}`$, with
$`a_1=a_3=1+p`$, $`a_2=a_4=1-p`$, leaves $`T^{p-1}`$ ($`T=-t`$) times

```math
\mathcal T_2=\int d^Dx_0\,(-x_{01}^2)^{-1-p}(-x_{02}^2)^{-1+p}(-x_{03}^2)^{-1-p},
\qquad x_{12}^2=x_{23}^2=0,\quad x_{13}^2=s .
```

By Feynman parameters, with $`S=-s`$ and $`D=4+2\eta`$, the parametric form
is

```math
-i\pi^{2+\eta}\,\frac{\Gamma(1+p-\eta)}{\Gamma(1+p)^2\Gamma(1-p)}\,S^{\eta-1-p}
\int d^3\alpha\,\delta\big(1-\textstyle\sum\alpha\big)\,\alpha_1^{\eta-1}\alpha_3^{\eta-1}\alpha_2^{-p}
=-i\pi^{2+\eta}\,S^{\eta-1-p}\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1-p+2\eta)} ,
```

the integral being Dirichlet's. Only $`\alpha_1\alpha_3S`$ appears because
the other two separations are null.

- The two factors $`\Gamma(\eta)`$ come from $`\alpha_1\to0`$ and
  $`\alpha_3\to0`$. They are $`\Gamma(D/2-a_2-a_3)`$ and
  $`\Gamma(D/2-a_1-a_2)`$: the collinear poles of the edges $`(2,3)`$ and
  $`(1,2)`$, each present whenever the powers of the two propagators at
  that edge sum to 2.
- The coefficient of the double pole is
  $`\Gamma(1+p)/(\Gamma(1+p)^2\Gamma(1-p))=\sin\pi p/(\pi p)`$.

In the normalisation in which the full integral is $`N_\eta(p)`$,

```math
\eta^2\times(\text{cusp-local model})=\eta^2\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1-p+2\eta)}
\ \longrightarrow\ \frac{\sin\pi p}{\pi p}\qquad(S=T=1).
```

At the cusps $`x_1,x_3`$ the roles of $`1\pm p`$ are exchanged, which
changes $`p\to-p`$ and leaves the limit the same. The four cusps add to
$`4\sin\pi p/(\pi p)`$, the known value. (`triangle_closed_form`.)

The step from the position-space triangle to its parametric form is the
same as Proposition 8 of GATE_6_CONTOUR. It has not been re-proved for the
triangle, whose contour at infinity differs. Section 3 does not need it.

## 3. The vertex lemma

Let $`S,T>0`$, $`p=\sigma+i\tau`$ with $`|\sigma|<1`$, and

```math
B_\eta(p)=\int_\Delta\alpha_1^{p}\alpha_3^{p}\alpha_2^{-p}\alpha_4^{-p}\,A^{\eta-2},
\qquad A=\alpha_1\alpha_3S+\alpha_2\alpha_4T ,
```

over the simplex $`\Delta=\{\alpha_i\ge0,\ \sum\alpha_i=1\}`$. By
Proposition 8 the regulated integral of $`z^p`$ is a constant times
$`\Gamma(2-\eta)B_\eta(p)/(\Gamma(1+p)^2\Gamma(1-p)^2)`$.

**Theorem.** For $`0<\eta\le\eta_0`$ with $`\eta_0<1-|\sigma|`$,

```math
B_\eta(p)=\frac{4\,\Gamma(1+p)\Gamma(1-p)}{\eta^2}\,S^{-1-p}\,T^{p-1}+R_\eta(p),
\qquad|R_\eta(p)|\le\frac{C\,(1+|\tau|)}{\eta},
```

with $`C`$ depending on $`\sigma`$, $`S`$, $`T`$ and $`\eta_0`$ only. Each
vertex of the simplex contributes one quarter of the leading term.

*Proof.* Fix $`\delta\le1/6`$. For each $`j`$ let $`V_j`$ be the set where
the three parameters other than $`\alpha_j`$ are at most $`\delta`$; then
$`\alpha_j\ge1/2`$ on $`V_j`$. All bounds are on moduli, so $`p`$ may be
replaced by $`\sigma`$ in them.

*The vertex $`V_2`$.* Use $`(\alpha_1,\alpha_3,\alpha_4)`$ and put
$`c=\alpha_1\alpha_3S`$,
$`m=(\alpha_1\alpha_3)^p\alpha_4^{-p}(c+\alpha_4T)^{\eta-2}`$. The integrand
is $`m\,\alpha_2^{-p}\big[(c+\alpha_2\alpha_4T)/(c+\alpha_4T)\big]^{\eta-2}`$,
and the quotient in brackets lies in $`[\alpha_2,1]`$. Hence

```math
|\text{integrand}-m|\le C(1+|\tau|)\,(\alpha_1+\alpha_3+\alpha_4)\,|m| .
```

(i) For every $`c>0`$,

```math
\int_0^\infty\alpha_4^{-p}(c+\alpha_4T)^{\eta-2}d\alpha_4
=c^{\,\eta-1-p}\,T^{p-1}\,\frac{\Gamma(1-p)\Gamma(1+p-\eta)}{\Gamma(2-\eta)} .
```

The part with $`\alpha_4>\delta`$ is at most
$`\int_\delta^\infty\alpha_4^{-\sigma}(\alpha_4T)^{\eta-2}`$, a constant, and
contributes $`O(1)`$ after integration against $`(\alpha_1\alpha_3)^\sigma`$.
So

```math
\int_{V_2}m=S^{\eta-1-p}T^{p-1}\,\frac{\Gamma(1-p)\Gamma(1+p-\eta)}{\Gamma(2-\eta)}
\Big(\int_0^\delta\alpha^{\eta-1}d\alpha\Big)^2+O(1)
=\frac{\Gamma(1-p)\Gamma(1+p)}{\eta^2}S^{-1-p}T^{p-1}+O\Big(\frac1\eta\Big),
```

because $`\delta^{2\eta}`$, $`S^\eta`$ and the Gamma quotient are
$`1+O(\eta)`$.

(ii) The error. $`\int|m|\alpha_1`$: by (i) with $`\sigma`$ it is at most
$`C\int\!\!\int\alpha_1(\alpha_1\alpha_3)^{\eta-1}=O(1/\eta)`$; the same for
$`\alpha_3`$. $`\int|m|\alpha_4`$: use
$`(c+x)^{\eta-2}\le c^{-\theta}x^{\eta-2+\theta}`$ with
$`\theta=(1+\sigma)/2`$. Then the $`\alpha_4`$ integrand is
$`\alpha_4^{\eta-1+(1-\sigma)/2}`$ and the rest is
$`(\alpha_1\alpha_3)^{(\sigma-1)/2}`$, both integrable: $`O(1)`$.

The vertex $`V_4`$ is the same with $`\alpha_2\leftrightarrow\alpha_4`$. At
$`V_1`$ and $`V_3`$ the roles of the two pairs are exchanged, giving
$`\Gamma(1+p)\Gamma(1-p-\eta)/\Gamma(2-\eta)\cdot T^{\eta-1+p}S^{-p-1}`$: the
same leading term.

*The rest, $`K=\Delta\setminus\bigcup V_j`$.* Let
$`\varepsilon=\delta/2`$. If $`\alpha_1,\alpha_2\le\varepsilon`$ and the
point is in neither $`V_3`$ nor $`V_4`$, then
$`\alpha_3,\alpha_4>\delta`$, so $`A\ge\delta(\alpha_1S+\alpha_2T)`$. With
$`\alpha_1=\lambda x`$, $`\alpha_2=\lambda(1-x)`$ that part is at most

```math
C\int_0^{2\varepsilon}\lambda^{\eta-1}d\lambda\int_0^1x^\sigma(1-x)^{-\sigma}\big(xS+(1-x)T\big)^{\eta-2}dx=O\Big(\frac1\eta\Big).
```

The pairs $`(1,4)`$, $`(3,2)`$, $`(3,4)`$ are the same. These are the four
edges of the simplex on which $`A=0`$: the collinear single poles. At every
other point of $`K`$ either $`\alpha_1,\alpha_3>\varepsilon`$ or
$`\alpha_2,\alpha_4>\varepsilon`$, so $`A\ge\varepsilon^2\min(S,T)`$ and the
integral is $`O(1)`$. $`\blacksquare`$

**Corollary.** With $`N_\eta(p)`$ the normalised regulated integral,

```math
\Big|N_\eta(p)-\frac{4\sin\pi p}{\pi p}\Big|\le C\,\eta\,(1+|\tau|)\,e^{2\pi|\tau|},
```

uniformly for $`\sigma`$ in a compact subset of $`(-1,1)`$; the constant
grows as $`|\sigma|\to1`$. The same bounds give absolute convergence of
$`B_\eta(p)`$ locally uniformly in $`p`$, so $`B_\eta`$ and $`N_\eta`$ are
holomorphic in $`p`$ for $`|\sigma|<1-\eta`$.

This replaces Lemma 4, Lemma 13 and the residue computation of
GATE_6_FUNCTIONAL in the chain of GATE_6_BOUNDS section 5 when both
diagonals are spacelike. Proposition 8 is still needed for step B. The
growth $`e^{2\pi|\tau|}`$ is worse than Lemma 4's $`e^{\pi|\tau|}`$ (and is
loose: numerically the growth is about $`e^{\pi|\tau|}`$). It is enough
provided $`\delta'`$ is chosen with $`W-\delta'>2\pi`$ in steps C and D,
which is possible for the three weights: $`W=3\pi-\delta`$ for
$`\bar F^{(2)}`$ and any $`W`$ for the others.

Remarks on the proof. The inequality used in (ii) needs
$`0\le\theta\le2-\eta`$, which holds. The Gamma quotient in (i) is
$`1+O(\eta\log(2+|\tau|))`$, absorbed by the factor $`1+|\tau|`$. The
parametric $`N_\eta`$ and the Mellin–Barnes one of GATE_6_FUNCTIONAL differ
by the factor $`(S/T)^{\eta/2}=1+O(\eta)`$ when $`S\ne T`$.

With a timelike diagonal $`A`$ changes sign on the simplex and this proof
does not apply as it stands. Section 6 extends it.

**Review.** An independent review re-derived sections 1 and 2, checked the
proof of the theorem line by line, and computed $`B_\eta(p)`$ by its own
direct integration over the simplex for six choices of $`p,S,T`$, including
$`S\ne T`$ and complex $`p`$: $`(\eta^2B_\eta-\text{leading term})/\eta`$
stays bounded in every case (about 5 at $`p=0.3`$, 85 at $`p=0.85`$). It
found no error and asked for the remarks above.

## 4. Numerical checks

| Check | Result |
| --- | --- |
| Measure identity and the formula for $`\rho`$ | exact (symbolic) |
| Dirichlet integral at $`p=0.3`$, $`\eta=0.25`$ | 18.5837284983 both ways |
| $`\alpha_4`$ integral of (i), real and complex $`p`$ | agree to 12 digits |
| $`\eta^2\times`$ triangle at $`p=0.3`$, $`\eta=10^{-1},10^{-2},10^{-3}`$ | 0.9655, 0.8707, 0.8596 → 0.8584 |
| $`(N_\eta-4\sin\pi p/\pi p)/\eta`$ at $`p=0.3`$, $`\eta=0.2\ldots0.025`$ | 0.65, 1.41, 1.87, 2.12: bounded |
| The same at $`p=0.2+1.5i`$ | bounded, about $`-18+16i`$ |

## 5. What this does and does not answer

- **Each cusp contributes $`M[z^p]`$, and the contributions add.** That is
  now a theorem at the level of the parametric integral.
- **The double pole is the product of two edge poles.** In the triangle the
  two $`\Gamma(\eta)`$ belong to the two edges at the cusp. This is the
  precise sense in which the edges determine the pole.
- **The vertex of the Feynman simplex with $`\alpha_j\to1`$ is the cusp
  $`x_j`$**, in the sense that the vertex model and the parametric form of
  the triangle at that cusp have the same double pole. They differ at order
  $`1/\eta`$. This identification is by comparison of the two formulas.
- **The exact formula of section 7** gives the complete effect of the
  regulator at a cusp for a general weight: a smearing of width $`\eta`$ and
  a factor $`(1+1/z)^{2\eta}`$ on the path between the lips. For the model
  the correction to $`M[F]`$ is $`O(\eta)`$. If the model equals the frozen
  position-space integral, which is not proved, the mismatch between
  quadrants cannot promote to a higher pole.
- **Section 8 derives the cusp-local model in position space**, by a
  different organisation of the integral from the one sketched next. In the
  exact variables of section 1 the local integral is
  $`\int\frac{du_1}{u_1}\frac{du_3}{u_3}|u_1u_3|^\eta\,G(u_1,u_3)`$, and its
  double pole is the sum over the four sign quadrants of the limits of
  $`G`$, which is the quadrant sum $`2\pi i\,M[F]`$ of GATE_6_CLOSURE. Two
  obstacles remain. On the real section the iterated integrals are not
  absolutely convergent, so the order cannot be exchanged without the
  contour. And for a weight growing like $`\log^kz`$ each quadrant has a
  $`\log^{k+1}`$ that cancels only in the sum, while the regulator differs
  between quadrants at $`|z|\sim1/|u|`$. A proof has to show that this
  mismatch does not promote to a higher pole. Section 7 shows that it does
not for the parametric model, by a different route. A direct derivation would integrate
$`u_1,u_3`$ over the real chain at fixed $`z`$ and must reproduce the kernel
of section 7. Section 8 does not follow this route: it integrates the
light-cone coordinate $`\alpha`$ first, by closing its contour, which
removes the sum over quadrants altogether.

## 6. A timelike diagonal

Let $`s>0`$, $`t<0`$, and $`T=|t|`$. By Lemma 14 of GATE_6_CONTOUR the
parametric integral is then the limit $`\theta\to\pi`$ of the integral with
$`S=|s|e^{i\theta}`$, and equals the integral over the rotated contour
$`\alpha_1=re^{-i\gamma}`$. Lemma 14 proves this for $`\gamma<\pi/4`$. At
$`\varphi=0`$ the rotation extends to any $`\gamma<\pi`$: at fixed
$`(\alpha_3,\alpha_2)`$ the zero of $`A`$ in $`\alpha_1`$ is on the positive
real axis and that of $`1+\alpha_1+\alpha_3`$ on the negative real axis,
both outside the sector swept; the arc at radius $`R`$ contributes
$`R^{\sigma-1-\eta}\to0`$; and both rays are absolutely integrable by
$`|A|\ge\sin(\gamma/2)(\dots)`$ and
$`|1+re^{-i\gamma}+\alpha_3|\ge\cos(\gamma/2)(1+r+\alpha_3)`$. The integrand is homogeneous of degree
$`-4`$, so the constraint may be put on the simplex in
$`(r,\alpha_2,\alpha_3,\alpha_4)`$:

```math
B_\eta(p)=e^{-i\gamma(1+p)}\int_\Delta r^{p}\alpha_3^{p}\alpha_2^{-p}\alpha_4^{-p}\,
h(r)^{-2\eta}\,\big(r\alpha_3S'+\alpha_2\alpha_4T\big)^{\eta-2},
\qquad S'=|s|e^{i(\pi-\gamma)},\quad h(r)=1-r(1-e^{-i\gamma}).
```

**Theorem (timelike).** With $`S=|s|e^{i\pi}`$,

```math
B_\eta(p)=\frac{4\,\Gamma(1+p)\Gamma(1-p)}{\eta^2}\,S^{-1-p}\,T^{p-1}+R_\eta(p),
\qquad|R_\eta(p)|\le\frac{C_\gamma\,(1+|\tau|)\,e^{\gamma\tau}}{\eta}.
```

*Proof.* The proof of section 3 goes through with $`S'`$ in place of $`S`$,
with these changes. (a) For $`c=r\alpha_3S'`$ and $`x\ge0`$,
$`|c+xT|\ge\sin(\gamma/2)\,(|c|+xT)`$, because the two terms make an angle
$`\pi-\gamma`$. Every bound on a modulus therefore holds with an extra
factor $`\sin(\gamma/2)^{\eta-2}`$. (b) The quotient
$`(c+\alpha_2\alpha_4T)/(c+\alpha_4T)`$ equals
$`1-(1-\alpha_2)\alpha_4T/(c+\alpha_4T)`$, which differs from 1 by at most
$`(1-\alpha_2)/\sin(\gamma/2)`$; take $`\delta\le\sin(\gamma/2)/6`$. (c) The
integral (i) holds for complex $`c`$ with $`|\arg c|<\pi`$ by analytic
continuation, and gives $`c^{\eta-1-p}`$. (d)
$`\cos(\gamma/2)\le|h|\le1`$, so $`|h^{-2\eta}-1|\le C_\gamma\eta`$, and the
factor changes the integral by $`O(\eta)\cdot O(1/\eta^2)`$. The leading
terms of the four vertices add to
$`e^{-i\gamma(1+p)}\,4\Gamma(1+p)\Gamma(1-p)\,S'^{-1-p}T^{p-1}/\eta^2`$, which
is the stated one and does not depend on $`\gamma`$. The factor
$`e^{\gamma\tau}`$ is the modulus of $`e^{-i\gamma(1+p)}`$. $`\blacksquare`$

**Corollary.** In the normalisation $`N_\eta`$, where
$`|(s/t)^p|=|s/t|^\sigma e^{-\pi\tau}`$,

```math
\Big|N_\eta(p)-\frac{4\sin\pi p}{\pi p}\Big|\le C_\gamma\,\eta\,(1+|\tau|)\,e^{(\gamma-\pi)\tau+2\pi|\tau|}.
```

Taking $`\gamma`$ small for $`\tau>0`$ and $`\gamma=\pi-\gamma'`$ for
$`\tau<0`$ gives the bound $`C\eta(1+|\tau|)e^{(2\pi+\gamma')|\tau|}`$ for any
$`\gamma'>0`$. Steps C and D of GATE_6_BOUNDS section 5 then need
$`W-\delta'>2\pi+\gamma'`$, which holds for the three weights. Step C also
needs $`N_\eta`$ holomorphic in $`p`$ between the lines: the rotated
integrand is entire in $`p`$ and is dominated, locally uniformly in $`p`$,
by $`\sin(\gamma/2)^{\eta-2}\cos(\gamma/2)^{-2\eta}`$ times the spacelike
integrand, so $`B_\eta`$ is holomorphic for $`|\sigma|<1-\eta`$. So Lemma 4
and Lemma 13 are no longer needed in any region; Lemma 14 is used only for
its first half, the passage to the rotated contour.

Numerically, at $`\arg(S/T)=\pi`$, $`(N_\eta-4\sin\pi p/\pi p)/\eta`$ stays
bounded as $`\eta\to0`$: about $`1.8+5.7i`$ at $`p=0.3`$, $`9+82i`$ at
$`p=0.2+1.5i`$ and $`-51+48i`$ at $`p=0.2-1.5i`$.

## 7. The exact formula for a general weight

Write the normalised cusp-local model of section 2 as

```math
\mathcal T_\eta(p)=\eta^2\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1-p+2\eta)} .
```

Two classical integrals,

```math
\frac{\Gamma(1+2\eta)}{\Gamma(1+p)\Gamma(1-p+2\eta)}=\frac1{2\pi}\int_{-\pi}^{\pi}e^{ip\phi}\,\big(1+e^{-i\phi}\big)^{2\eta}\,d\phi,
\qquad
\frac{\Gamma(1+p-\eta)}{\Gamma(1+p)}=\frac{1}{\Gamma(\eta)}\int_0^\infty e^{-px}\,e^{-(1-\eta)x}\big(1-e^{-x}\big)^{\eta-1}dx ,
```

the first a Beta integral on the unit circle and the second Euler's, give
for $`\operatorname{Re}p>\eta-1`$

```math
\mathcal T_\eta(p)=\frac{c_\eta}{2\pi}\int_{-\pi}^{\pi}d\phi\,\big(1+e^{-i\phi}\big)^{2\eta}\int_0^\infty dx\,\kappa_\eta(x)\,e^{p(i\phi-x)},
```

```math
c_\eta=\frac{\Gamma(1+\eta)^2\Gamma(1-\eta)}{\Gamma(1+2\eta)}\to1,\qquad
\kappa_\eta(x)=\frac{e^{-(1-\eta)x}\,(1-e^{-x})^{\eta-1}}{\Gamma(\eta)\Gamma(1-\eta)},\qquad\int_0^\infty\kappa_\eta=1 .
```

**Theorem (regulated circle).** Let $`0<\eta<1`$, let $`f(L)=F(e^L)`$
satisfy the hypothesis of GATE_6_BOUNDS Lemma 3 with $`W>\pi`$, and
define the cusp-local model of
$`F`$ by linearity, $`\mathcal T_\eta[F]=\int c(p)\,\mathcal T_\eta(p)\,dp`$.
Then, exactly,

```math
\mathcal T_\eta[F]=\frac{c_\eta}{2\pi}\int_{-\pi}^{\pi}d\phi\,\big(1+e^{-i\phi}\big)^{2\eta}\int_0^\infty dx\,\kappa_\eta(x)\,F\big(e^{i\phi-x}\big).
```

At the cusps $`x_1,x_3`$ the same holds with $`F(e^{-i\phi+x})`$.

*Proof.* Insert the formula for $`\mathcal T_\eta(p)`$ and exchange the
integrals. On the lines $`\operatorname{Re}p=\mp\sigma`$,
$`|e^{p(i\phi-x)}|\le e^{\pi|\tau|+\sigma x}`$, $`\kappa_\eta`$ decays like
$`e^{-(1-\eta)x}`$ and $`|c(p)|\le Ce^{-(W-\delta')|\tau|}`$ with
$`W-\delta'>\pi`$, so the triple integral converges absolutely for
$`\sigma<\min(1/m,1-\eta)`$. The points $`L=i\phi-x`$ have
$`|\operatorname{Im}L|\le\pi<W`$, inside the strip of Lemma 3. By Lemma 3, $`\int c(p)e^{pL}dp=f(L)`$. $`\blacksquare`$

**What it says.**

- The cusp-local model is an integral of $`F`$ over the path between the
  two lips of $`z=-1`$, here the unit circle, moved to $`|z|=e^{-x}`$ at the
  cusps $`x_2,x_4`$ and to $`|z|=e^{x}`$ at $`x_1,x_3`$, with $`x`$
  distributed by $`\kappa_\eta`$. Only the principal sheet enters.
- $`\kappa_\eta`$ has mean $`\psi(1)-\psi(1-\eta)\approx\eta\pi^2/6`$, and
  $`(1+e^{-i\phi})^{2\eta}\to1`$ with an integrable logarithm at
  $`\phi=\pm\pi`$. If $`|f|+|f'|\le A(1+|L|)^k`$ on
  $`|\operatorname{Im}L|\le\pi`$, $`\operatorname{Re}L\le0`$ (or $`\ge0`$ at
  $`x_1,x_3`$), and $`\eta\le\eta_0<1`$, then

```math
\big|\mathcal T_\eta[F]-M[F]\big|\le C_{k,\eta_0}\,A\,\eta .
```

  The constant is not computed; it grows quickly with $`k`$.
- So $`\eta^2\times`$(model) is bounded and tends to $`M[F]`$ for every such
  weight, whatever its growth in $`\log z`$. "No pole beyond the second" is
  meant in this sense; for meromorphy in $`\eta`$ see section 8.
- For the full integral, summing the vertex theorems against $`c(p)`$ gives
  $`\epsilon^2st\,I_\eta[F]=\sum_j\mathcal T^{(j)}_\eta[F]+O(\eta)=4M[F]+O(\eta)`$.
  This needs more than the theorem above: $`W>\Phi`$ for step A of
  GATE_6_BOUNDS section 5, $`W-\delta'>2\pi`$ (spacelike) or
  $`2\pi+\gamma'`$ (timelike), a square where the certificate passes,
  Proposition 8 and the first half of Lemma 14. The three weights
  $`\bar F^{(0)},\bar F^{(1)},\bar F^{(2)}`$ satisfy these. The difference
  between the full integral and the sum of the four models is a genuine
  single pole: $`(N_\eta-\sum_j\mathcal T^{(j)}_\eta)/\eta\to-0.86`$ at
  $`p=0.3`$.

**Review.** An independent review re-derived sections 6 and 7, found no
error, and confirmed with its own code the rotated integral at
$`\gamma=0.7,1.6,2.4`$ against Mellin–Barnes to $`10^{-6}`$, the double
integral of this section against the closed form to $`10^{-15}`$, and the
$`O(\eta)`$ rates for $`L,L^2+\pi^2,L^3,L^4,L^8`$. The qualifications it
asked for are those stated here.

**Status.** The formula is exact for the cusp-local model defined through
its parametric form. Section 8 derives the model from the Lorentzian
domain and shows that the two agree.

Numerical checks (`regulated_circle_formula`):

| Weight | $`\eta`$ | Formula | Reference |
| --- | --- | --- | --- |
| $`z^{0.3}`$ | 0.2 | 1.04665719 | 1.04665717 (closed form) |
| $`z^{0.2+i}`$ | 0.15 | $`2.9136875-0.7260480i`$ | $`2.9136869-0.7260478i`$ |
| $`z^{-0.4}`$ | 0.05 | 0.77929882 | 0.77929878 |
| $`\log^2z+\pi^2`$ | 0.2, 0.05, 0.0125 | 8.7060, 7.1154, 6.7163 | same from $`\partial_p^2`$ of the closed form; $`M[F]=6.5797`$ |

For the last weight $`(\mathcal T_\eta[F]-M[F])/\eta`$ is 10.6, 10.7, 10.9.

## 8. Derivation from the Lorentzian domain

Set $`s=-1`$ (so $`2p_1\!\cdot p_2=-1`$); the dimension is $`4+2\eta`$. With
$`x_0=x_2+\alpha p_1+\beta p_2+y_\perp`$ and
$`\rho=|y_\perp|^2`$,

```math
w=x_{02}^2=-\alpha\beta-\rho,\qquad u_1=x_{01}^2=w-\beta,\qquad u_3=x_{03}^2=w+\alpha,\qquad z=-\frac{w}{u_1u_3},
```

each with $`-i0`$, and $`L=\log z`$ on the branch fixed by the
prescription. The cusp-local model in position space is the iterated
integral

```math
\mathcal T^{\rm L}_\eta[F]=\int_{\mathbb R}d\beta\int_0^\infty d\rho\,\rho^\eta\int_{\mathbb R}d\alpha\;\frac{f(L)}{w\,u_1\,u_3}.
```

The measure is $`d^Dx_0=\tfrac12\,\pi^{1+\eta}/\Gamma(1+\eta)\;d\alpha\,d\beta\,\rho^\eta d\rho`$;
the constant is omitted. The $`-i0`$ on all three invariants is the limit
$`\psi\to0`$ of rotating the $`\rho`$ ray to $`\rho=e^{i\psi}r`$,
$`0<\psi<\pi`$, measure included (a factor $`e^{i\psi(1+\eta)}`$). The
integral does not depend on $`\psi`$, and for $`\psi>0`$ the $`\alpha`$
integrand is non-singular on the real line. The order of integration is
part of the definition: $`\alpha`$ innermost, then $`\rho`$, then
$`\beta`$. The integral is not absolutely convergent on $`\mathbb R^3`$,
since its vanishing off the edge segment is a cancellation, and other
orders are not treated here.

**Theorem (Lorentzian derivation).** Let $`0<\eta<1`$ and let $`f`$ be
holomorphic on $`-2\pi-\delta<\operatorname{Im}L<\pi+\delta`$ with
$`|f(L)|\le Ae^{\sigma_-|\operatorname{Re}L|}`$ for
$`\operatorname{Re}L\le0`$ and
$`|f(L)|\le Ae^{\sigma_+\operatorname{Re}L}`$ for $`\operatorname{Re}L\ge0`$,
where $`\sigma_-<1-\eta`$ and $`\sigma_+<1+\eta`$. This includes
polynomially bounded $`f`$ and $`e^{pL}`$ with
$`\eta-1<\operatorname{Re}p<1+\eta`$. Then

```math
\mathcal T^{\rm L}_\eta[F]=i\int_0^1d\beta\,\beta^{\eta-1}(1-\beta)^{\eta}\int_0^\infty\frac{dy\,y^{\eta-1}}{1+y}
\int_{-\pi}^{\pi}d\phi\,\big(1+e^{-i\phi}\big)^{\eta}\,f\big(i\phi-\log(1-\beta)-\log(1+y)\big).
```

*Proof.* **Step 1: support on the edge segment.** All three invariants are
linear in $`\alpha`$. Their zeros, displaced by the prescription, are

```math
\alpha_w=-\frac{\rho+i0}{\beta},\qquad\alpha_1=-1-\frac{\rho+i0}{\beta},\qquad\alpha_3=\frac{\rho+i0}{1-\beta}.
```

For $`\beta<0`$ all three lie in the upper half plane, for $`\beta>1`$ all
three in the lower. The integrand is holomorphic in the other half plane
and decays at least like $`|\alpha|^{-3+\max(\sigma_-,\sigma_+-2)}`$ up to
logarithms, so the $`\alpha`$ integral vanishes. In that half plane all
three invariants have negative imaginary part, so
$`\operatorname{Im}L\in(-2\pi,\pi)`$ and $`f`$ is evaluated inside the
strip of the hypothesis. (On the segment, in step 2, only
$`|\operatorname{Im}L|<\pi`$ occurs.) Hence only
$`0<\beta<1`$ contributes: the integral is supported on the edge segment.
This is the exact form of the vanishing of the edge jump off the segment.

**Step 2: the discontinuity.** For $`0<\beta<1`$, $`\alpha_w`$ and
$`\alpha_1`$ are in the lower half plane and $`\alpha_3`$ in the upper.
Close the contour in the upper half plane around the cut from
$`\alpha_3`$. With $`u_3=(1-\beta)(\alpha-\alpha_3)`$ as the variable,

```math
\int_{\mathbb R}d\alpha\,\frac{f(L)}{w\,u_1u_3}=-\frac{1}{1-\beta}\oint_H\frac{du_3}{u_3}\,\frac{f(L)}{w\,u_1},
\qquad w=-\frac{\beta u_3+\rho}{1-\beta},\quad u_1=w-\beta,
```

where $`H`$ comes in along $`\arg u_3=0`$, circles the origin and goes out
along $`\arg u_3=-2\pi`$. On the two sides $`w`$ and $`u_1`$ are negative
and $`\arg z=\mp\pi`$: the two sides of $`H`$ are the two lips of the
negative $`z`$ axis. When $`f`$ decays at $`\operatorname{Re}L\to+\infty`$
the circle can be shrunk, and the $`\alpha`$ integral equals
$`\frac{1}{1-\beta}\int_0^\infty\frac{du_3}{u_3}\frac{f(\ell-i\pi)-f(\ell+i\pi)}{|w||u_1|}`$,
the discontinuity of $`F`$ across $`z<0`$ against a positive measure.

**Step 3: three independent variables.** Put $`u_3=\rho t`$, so
$`w=-\rho(1+\beta t)/(1-\beta)`$, with $`t`$ on a Hankel contour whose
circle has a radius below $`1/\beta`$ that does not depend on $`\rho`$.
Then $`\operatorname{Re}(1+\beta t)>0`$ on the contour. Exchange the
$`\rho`$ and $`t`$ integrals: on the product the integrand is at most
$`\rho^{\eta-1}|f|/(|t||1+\beta t||1+y|)`$ with $`t`$ bounded away from 0,
which is integrable under the hypothesis (this is where
$`\sigma_-<1-\eta`$ is used, at $`y\to\infty`$; no decay at
$`\operatorname{Re}L\to+\infty`$ is needed because the circle is kept).
At fixed $`t`$ rotate the $`\rho`$ ray so that
$`y=\rho(1+\beta t)/(\beta(1-\beta))`$ is real and positive. The rotation is
by less than $`\pi`$, crosses neither $`1+\beta t=0`$ nor $`y=-1`$, and the
arc decays. Then exchange back, so that $`t`$ is integrated at fixed
$`y`$. The outer $`\beta`$ integral converges at $`\beta\to1`$ because
$`\sigma_+<1+\eta`$.
Then $`u_1=-\beta(1+y)`$ and

```math
z=-\frac{1+\beta t}{\beta t}\cdot\frac{1}{(1-\beta)(1+y)} .
```

With $`x=\beta t/(1+\beta t)`$ the measure becomes
$`\beta^{\eta-1}(1-\beta)^\eta d\beta\cdot y^{\eta-1}(1+y)^{-1}dy\cdot(1-x)^\eta dx/x`$,
and $`z=-1/(x(1-\beta)(1+y))`$. The contour in $`x`$ runs around the
segment $`(0,1)`$ from $`x=1`$ back to $`x=1`$; the integrand is
holomorphic in the disc $`|x|<1`$ cut along that segment, so the contour
may be taken to be the unit circle. Writing $`x=-e^{-i\phi}`$ gives the
formula. $`\blacksquare`$

**The three variables.**

| Variable | Meaning | Measure | Effect |
| --- | --- | --- | --- |
| $`\beta\in(0,1)`$ | position along the edge | $`\beta^{\eta-1}(1-\beta)^\eta d\beta`$ | pole $`1/\eta`$ from the end $`\beta\to0`$, the cusp |
| $`y\in(0,\infty)`$ | collinear scale, $`(\beta u_3+\rho)/(\beta(1-\beta))`$ | $`y^{\eta-1}(1+y)^{-1}dy`$ | pole $`1/\eta`$ from $`y\to0`$, the collinear limit |
| $`\phi\in(-\pi,\pi)`$ | the path between the lips | $`(1+e^{-i\phi})^\eta d\phi`$ | the functional $`M[F]`$ |

The argument of $`f`$ is the sum of the three contributions. The measure is
a product, so the region where both collinear scales vanish,
$`\beta\to0`$ and $`y\to0`$, contributes exactly the product of the two
poles: $`B(\eta,1+\eta)\,B(\eta,1-\eta)`$. Nothing is left to estimate
there.

**Consequences.**

1. For $`F=z^p`$, $`\eta-1<\operatorname{Re}p<1+\eta`$, the three integrals
   are Beta integrals and

```math
\mathcal T^{\rm L}_\eta[z^p]=2\pi i\,\Gamma(1+\eta)\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1+2\eta-p)}
=\frac{2\pi i\,\Gamma(1+\eta)}{\eta^2}\,\mathcal T_\eta(p).
```

   This is the closed form of section 2, obtained here without Feynman
   parameters; with the omitted constant of the measure it reproduces the
   prefactor $`-i\pi^{2+\eta}`$ of section 2 exactly. The range of $`p`$ is
   that of the iterated integral; the closed form continues beyond it. The
   theorem gives a general weight directly. The parametric model of
   section 7 is defined by linearity from powers, so the two agree for the
   weights of Lemma 3 (with $`W>2\pi`$ so that $`f`$ exists on the strip,
   and lines $`\sigma<\min(1/m,1-\eta)`$), by Fubini on the product
   formula. The formula of section 7 is then a second expression for the
   same quantity.
2. Normalising by $`2\pi i\,\Gamma(1+\eta)/\eta^2`$, the two outer measures
   become probability measures concentrated at $`\beta=0`$ and $`y=0`$ with
   widths of order $`\eta`$, up to a factor tending to 1, and

```math
\frac{\eta^2\,\mathcal T^{\rm L}_\eta[F]}{2\pi i\,\Gamma(1+\eta)}\;\longrightarrow\;\frac1{2\pi}\int_{-\pi}^{\pi}f(i\phi)\,d\phi=M[F],
```

   with an $`O(\eta)`$ error for $`|f|+|f'|\le A(1+|L|)^k`$ on
   $`|\operatorname{Im}L|\le\pi`$.
3. The right side has the form
   $`\int\!\!\int\beta^{\eta-1}y^{\eta-1}h(\beta,y,\eta)\,d\beta\,dy`$ with
   $`h`$ jointly analytic at $`\beta=y=0`$, because $`f`$ is holomorphic
   near $`|\operatorname{Im}L|\le\pi`$, and with uniform convergence at
   $`\beta\to1`$ and $`y\to\infty`$ for small $`|\eta|`$. Subtracting
   $`h(0,y)`$, $`h(\beta,0)`$ and $`h(0,0)`$ shows that the right side
   continues meromorphically in $`\eta`$ to a neighbourhood of 0, with a
   pole of order at most two. This is a statement about the cusp-local
   model, not about the full integral.
4. Step 2 at fixed $`\beta`$ is the mechanism of the edge statement of
   GATE_6_CLOSURE: the two sides of the edge are the two lips. The leading
   density is $`2\pi i\,M[F]\,d\beta/\beta`$ per $`1/\eta`$, which matches
   the edge statement once the far propagator is restored. The two are
   identified as mechanisms; the edge statement itself is a jump in
   $`q_3`$ of a transverse integral in four dimensions, and their identity
   is not proved as a formula. With that understood, for the cusp-local
   model the edge discontinuity, integrated along the edge with its
   endpoint at the cusp, gives the double pole.

**What this does not cover.** The full four-point integral. With the far
propagator present, $`x_{04}^2`$ is not linear in $`\alpha`$ alone, since it
depends on the angle of $`y_\perp`$, and step 1 fails as it stands. The
passage from the full integral to the four cusp-local models is proved for
the parametric integral (sections 3 and 6), not in position space. The
iterated integral used here is the natural real-domain definition; its
equality with the integral over the global contour of GATE_6_CONTOUR is
not proved for the triangle.

**Review.** An independent review re-derived every step, found no wrong
sign or constant, and asked for the hypotheses and the exchanges of order
now stated. With its own code it confirmed: the $`\alpha`$ integral below
$`3\times10^{-12}`$ off the segment; the discontinuity form to 8 to 15
digits; the product formula against the closed form to $`10^{-11}`$ at six
values of $`(p,\eta)`$; independence of $`\psi`$; and the original triple
integral in $`(\beta,\rho,\alpha)`$ against the closed form,
$`27.8301i`$ against $`27.8280i`$ at $`p=0.3`$, $`\eta=0.5`$ and
$`-0.5606+43.0209i`$ against $`-0.5578+43.0210i`$ at $`p=-0.2+0.3i`$,
$`\eta=0.4`$, limited by its quadrature.

**Numerical checks.**

| Check | Result |
| --- | --- |
| $`\alpha`$ integral off the segment, $`\beta=-0.4,1.3,1.6`$ | $`10^{-12}`$ |
| $`\alpha`$ integral against the discontinuity form, $`z^{-0.35+0.4i}`$, $`\beta=0.3,0.7`$ | agree to 6 digits |
| The same for $`(\log^2z+\pi^2)z^{-0.2}`$, $`\beta=0.4`$ | $`186.272665+36.524833i`$ both ways |
| Product formula against the closed form, $`z^{0.3}`$, $`z^{0.2+i}`$, $`z^{-0.4}`$ | agree to $`10^{-4}`$ (quadrature) |
| Product formula against section 7 for $`\log^2z+\pi^2`$ | 8.7058 against 8.7060; 7.1151 against 7.1154 |
