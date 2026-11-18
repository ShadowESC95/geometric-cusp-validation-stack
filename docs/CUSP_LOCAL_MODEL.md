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

5. For the full four-point integrand the same contour closing is exact
   (section 9): integrated first along the null direction of an edge, the
   integral is supported in the slab between two null hyperplanes and
   splits into two pieces, one for each of two opposite edges.

6. Each of the two pieces is an explicit integral over a cube and equals
   its two end terms up to a single pole (sections 10 and 11). With
   section 9 this gives the double pole of the full integral,
   $`4\sin\pi p/(\pi p)`$, from the Lorentzian domain, for the weight
   $`z^p`$.

7. The same holds for every weight holomorphic and polynomially bounded on
   a strip $`|\operatorname{Im}\log z|<W`$ with $`W>2\pi`$, in particular
   for the three benchmark weights: the double pole is $`4M[F]`$
   (section 12).

Conditions and what is **not** proved: the result is for the iterated
integral in a fixed order with the displacement removed first, and for
spacelike diagonals. A timelike diagonal is not treated in position space,
and the equality of the iterated integral with the integral over the
global contour is not proved.

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
propagator present step 1 still holds, in a stronger form (section 9), but
the product structure of step 3 is lost. The passage from the full
integral to the four cusp-local models is proved for the parametric
integral (sections 3 and 6), not in position space. The
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

## 9. The full integrand: exact reduction to two edges

Keep the far propagator. Set $`s=-1`$, $`t=-T`$, and let $`\kappa`$ be the
component of $`y_\perp`$ along the transverse part of the far edge, so
$`\kappa^2\le\rho`$. By the exact formula for the far invariant,

```math
x_{04}^2=-T(1-\beta)-\alpha\,(T+\beta)-\rho+c\,\kappa,\qquad c=2\sqrt{T(1+T)},
\qquad z=\frac{w\,x_{04}^2}{u_1u_3\,T}.
```

So all four invariants are linear in $`\alpha`$ at fixed $`\beta`$ and fixed
transverse position. Their zeros, displaced by the prescription, lie in
these half planes:

| | $`\beta<-T`$ | $`-T<\beta<0`$ | $`0<\beta<1`$ | $`\beta>1`$ |
| --- | --- | --- | --- | --- |
| $`x_{02}^2`$, $`x_{01}^2`$ | upper | upper | lower | lower |
| $`x_{03}^2`$ | upper | upper | upper | lower |
| $`x_{04}^2`$ | upper | lower | lower | lower |

**Theorem (two edges).** Let $`f`$ be holomorphic and polynomially bounded
on $`|\operatorname{Im}L|<2\pi+\delta`$. At fixed $`\beta`$, $`\rho`$ and
$`\kappa`$,

```math
\int_{\mathbb R}d\alpha\,\frac{f(L)}{x_{01}^2x_{02}^2x_{03}^2x_{04}^2}=
\begin{cases}
0,&\beta<-T\ \text{or}\ \beta>1,\\[2pt]
\text{keyhole around the zero of }x_{03}^2,&0<\beta<1,\\[2pt]
\text{keyhole around the zero of }x_{04}^2,&-T<\beta<0 .
\end{cases}
```

*Proof.* As in step 1 of section 8: where all four zeros lie in one half
plane the contour closes in the other; where exactly one zero lies apart,
the contour is deformed onto a keyhole around a cut from that zero. Where
the integral vanishes, the four invariants have negative imaginary part in
the closing half plane, so $`|\operatorname{Im}L|<2\pi`$. In the keyhole
ranges the invariant whose zero is encircled changes sign, and
$`\operatorname{Im}L`$, which is harmonic and bounded, stays within
$`[-2\pi,2\pi]`$ by its values on the real line and on the two lips. The
integrand decays like $`|\alpha|^{-4}`$; no growth condition on $`f`$ is
needed for this, because $`z`$ tends to a finite non-zero constant as
$`|\alpha|\to\infty`$, so the theorem also covers $`z^p`$. $`\blacksquare`$

**Geometry.** The surfaces $`\beta=`$ const are the null hyperplanes
orthogonal to $`p_1`$. The one with $`\beta=1`$ passes through $`x_3`$, the
one with $`\beta=0`$ contains the edge from $`x_1`$ to $`x_2`$, and the one
with $`\beta=-T`$ passes through $`x_4`$. So the integral along the null
direction $`p_1`$ is supported in the slab between the hyperplanes through
$`x_4`$ and $`x_3`$. The hyperplane containing the edge $`(1,2)`$ divides
the slab in two. On one side the integral is a discontinuity across the
light cone of $`x_3`$ and belongs to the edge $`(2,3)`$; on the other it is
a discontinuity across the light cone of $`x_4`$ and belongs to the edge
$`(4,1)`$. For $`F=1`$, the box, each keyhole is a residue.

This is exact for the full integrand, for every weight in the class and
for every value of the displacement $`\varepsilon>0`$ in $`-i\varepsilon`$.

**Numerical check** (`full_light_cone_integral`, `full_keyhole_integral`),
at $`T=1.7`$, $`\rho=0.35`$, $`\kappa=0.3`$, $`\varepsilon=0.02`$, for the
weights $`z^{-0.35+0.4i}`$, $`\log^2z+\pi^2`$ and 1:

| $`\beta`$ | $`\alpha`$ integral | Keyhole |
| --- | --- | --- |
| $`-2.5,\ -1.9,\ 1.2,\ 2.0`$ | below $`10^{-12}`$ in every case | 0 |
| $`-1.0`$, $`\log^2z+\pi^2`$ | $`-0.0681280-1.729232i`$ | the same to 7 digits |
| $`-0.3`$, $`z^{-0.35+0.4i}`$ | $`-6.201512-21.85390i`$ | the same |
| $`0.4`$, $`F=1`$ | $`-2.205243-13.70115i`$ | the same |
| $`0.8`$, $`z^{-0.35+0.4i}`$ | $`0.7404036-0.2133773i`$ | the same |

**What remains.** Each edge piece has two ends, and the double pole comes
from the ends.

- *Near end* ($`\beta\to0^+`$ for the edge $`(2,3)`$: the cusp $`x_2`$).
  There $`x_{04}^2\to t(1-\beta)`$ and the integrand tends to that of the
  cusp-local model of section 8. The corrections are of relative order
  $`\beta`$, $`y`$ and $`\sqrt{\beta y}`$ in the variables of section 8,
  each of which removes at least one of the two poles.
- *Far end* ($`\beta\to1`$: the cusp $`x_3`$). With $`b=1-\beta`$ and the
  scaling $`\rho,u_3\sim b^2`$, $`\kappa\sim b`$, the measure gives
  $`\int db\,b^{2\eta-1}=1/(2\eta)`$, and the scaled integral has two
  collinear pinches, one for each edge at $`x_3`$, each giving $`1/\eta`$:
  in total $`\tfrac{1}{2\eta}(\tfrac1\eta+\tfrac1\eta)=1/\eta^2`$, one unit.

The far-end coefficient is computed in section 10: it is one unit, as at
the near end. The estimate that turns the two leading terms into "edge
piece $`=`$ two end terms $`+\,O(1/\eta)`$" is section 11.

## 10. The far end of an edge piece

Take the piece of section 9 for the edge $`(2,3)`$, $`0<\beta<1`$, and its
far end $`b=1-\beta\to0`$, the cusp $`x_3`$. After closing the $`\alpha`$
contour, with $`u_3`$ on the keyhole,

```math
x_{02}^2=-\frac{\beta u_3+\rho}{b},\qquad
x_{04}^2=-\frac{(T+1)\,|y_\perp-y_*|^2+(T+\beta)\,u_3}{b},\qquad |y_*|^2=b^2\,\frac{T}{1+T},
```

where the second is exact: $`x_{04}^2`$ is a negative definite quadratic in
the transverse position, vanishing only at $`y_\perp=y_*`$, $`u_3=0`$, which
is the edge $`(3,4)`$. Scale $`y_\perp=b\sqrt{T/(1+T)}\,n`$ and
$`u_3=b^2\tfrac{T}{1+T}v`$. As $`b\to0`$, $`x_{01}^2\to s`$ and

```math
z\to-\frac{(v+|n|^2)(v+|n-e|^2)}{v},\qquad|e|=1,
```

with no dependence on $`T`$. The $`b`$ integral gives
$`\int_0db\,b^{2\eta-1}=1/(2\eta)`$, and what multiplies it is the universal
integral

```math
\Sigma_\eta[F]=\int d^{2+2\eta}n\oint_H\frac{dv}{v}\,\frac{F(z)}{(v+|n|^2)(v+|n-e|^2)} .
```

The two factors in the denominator are the two edges at $`x_3`$: the
integrand is symmetric under $`n\to e-n`$.

**Proposition (far-end coefficient).** For $`F=z^p`$,
$`\eta-1<\operatorname{Re}p<0`$,

```math
\Sigma_\eta[z^p]=2i\sin(\pi p)\,J,\qquad
J=\int_0^\infty dv\,v^{-1-p}\int d^{2+2\eta}n\,(v+|n|^2)^{p-1}(v+|n-e|^2)^{p-1}
=\pi^{1+\eta}\,\frac{\Gamma(-p)\,\Gamma(1-p-\eta)\,\Gamma(\eta)^2}{\Gamma(1-p)^2\,\Gamma(2\eta)} .
```

*Proof.* The keyhole gives $`2i\sin\pi p`$ times the integral over
$`v>0`$, as in section 8. Combine the two factors with one Feynman
parameter $`\xi`$:
$`\xi|n|^2+(1-\xi)|n-e|^2=|n-(1-\xi)e|^2+\xi(1-\xi)`$. The integral over
$`n`$ gives
$`\pi^{1+\eta}\Gamma(1-2p-\eta)/\Gamma(2-2p)\,(v+\xi(1-\xi))^{2p-1+\eta}`$,
the integral over $`v`$ a Beta function,
$`(\xi(1-\xi))^{p-1+\eta}\Gamma(-p)\Gamma(1-p-\eta)/\Gamma(1-2p-\eta)`$, and
the integral over $`\xi`$ is
$`\int_0^1(\xi(1-\xi))^{\eta-1}d\xi=\Gamma(\eta)^2/\Gamma(2\eta)`$.
$`\blacksquare`$

The two poles in $`\Gamma(\eta)^2/\Gamma(2\eta)\approx2/\eta`$ come from
$`\xi\to0`$ and $`\xi\to1`$: the two collinear limits, one for each edge at
$`x_3`$. With the $`1/(2\eta)`$ from $`b`$,

```math
\frac{1}{2\eta}\,\Sigma_\eta[z^p]\;\longrightarrow\;-\frac{2\pi^2i}{\eta^2}\,\frac{\sin\pi p}{\pi p} .
```

**Comparison with the near end.** In the same units, with the transverse
measure $`d^{2+2\eta}y`$ and the prefactor $`st`$, the two leading terms of
the edge piece are

```math
\text{near end }(x_2):\quad-2\pi^{2+\eta}i\,\frac{\Gamma(\eta)^2\,\Gamma(1+p-\eta)}{\Gamma(1+p)^2\,\Gamma(1+2\eta-p)},
\qquad
\text{far end }(x_3):\quad-2\pi^{2+\eta}i\,\Big(\frac{T}{1+T}\Big)^{\eta}\frac{\Gamma(\eta)^2\,\Gamma(1-p-\eta)}{\Gamma(1+p)\,\Gamma(1-p)^2\,\Gamma(1+2\eta)} .
```

Both have the double pole $`-2\pi^2i\,\sin\pi p/(\pi p)\,/\eta^2`$, with the
same sign. So each edge piece carries two units, one per end, and the two
edge pieces carry the four units of the known result. The mechanism at the
two ends is different:

| | Endpoint | Collinear | Product |
| --- | --- | --- | --- |
| Near end | $`\int\beta^{\eta-1}d\beta=1/\eta`$ | one edge, $`1/\eta`$ | $`1/\eta^2`$ |
| Far end | $`\int b^{2\eta-1}db=1/(2\eta)`$ | two edges, $`2/\eta`$ | $`1/\eta^2`$ |

For a general weight the far-end integral is treated in the same way as
section 8 and gives $`M[F]`$; this is not written out.

**Numerical checks** (`far_end_transverse`, `far_end_numerical`):

| Check | Result |
| --- | --- |
| Transverse integral against its Feynman-parameter form, $`v=0.3`$, $`p=-0.3`$, $`\eta=0.5`$ | 8.5202446630 against 8.5202446627 |
| The same at $`v=0.02`$, $`p=-0.6`$, $`\eta=0.3`$ | agree to $`10^{-13}`$ |
| $`J`$ against the closed form at $`(p,\eta)=(-0.3,0.5),(-0.6,0.3),(-0.2,0.1)`$ | agree to $`10^{-10}`$ or better |
| $`\eta^2\times`$ far-end term at $`p=0.3`$, $`T=1.7`$, $`\eta=10^{-1},10^{-2},10^{-3}`$ | 0.8818, 0.8600, 0.8585 → 0.8584 |

**Status.** The far-end leading term is exact for the scaled integrand.
That the edge piece minus its two leading terms is $`O(1/\eta)`$ is proved
in section 11, where the far-end term is also re-derived by a second
route.

## 11. The edge piece: exact reduction and the remainder

Units as in section 9: $`s=-1`$, $`t=-T`$, $`b=1-\beta`$,
$`\varrho=T/(T+1)`$; $`p=\sigma+i\tau`$. The edge piece for the edge
$`(2,3)`$ is

```math
A(p;T)=st\int_0^1d\beta\int d^{2+2\eta}y_\perp\int_{\mathbb R}d\alpha\,\frac{z^p}{x_{01}^2x_{02}^2x_{03}^2x_{04}^2},
```

and the full integral with the same prefactor and measure is $`A`$ plus the
piece for the edge $`(4,1)`$. The order is part of the definition: the
$`\alpha`$ integral first, with the displacement $`\varepsilon\to0`$ taken
at fixed $`(\beta,y_\perp)`$, then the transverse integral, then
$`\beta`$.

### 11.1 Reduction to an integral over a cube

**Proposition.** For $`\eta-1<\operatorname{Re}p<0`$,

```math
A(p;T)=-2i\pi^{2+\eta}\,\frac{\Gamma(2-\eta)\,T^{1-p}}{\Gamma(1-p)^2\,\Gamma(1+p)^2}\;\mathcal J(p;T),
```

```math
\mathcal J=\frac{1}{T+1}\int_0^1d\beta\,b^{\eta-p}\int_0^1d\xi\int_0^1d\mu\;
(1-\xi)^{-p}\,\xi^{p}\,\mu^{-p}\,(1-\mu)^{\eta-1}\,\big(\xi\beta+\varrho b\mu\big)^{\eta-2}\,\big((T+1)\beta+Tb\mu\big)^{p}.
```

*Proof.* By section 9 the $`\alpha`$ integral is a keyhole around the zero
of $`x_{03}^2`$. For $`\operatorname{Re}p<0`$ the circle shrinks and the
keyhole is $`2i\sin\pi p`$ times an integral over $`u_3>0`$ of a positive
function (on the cut $`x_{02}^2,x_{01}^2,x_{04}^2<0`$ and $`\arg z=\mp\pi`$):

```math
A=2iT^{1-p}\sin(\pi p)\int_0^1d\beta\,b^{2-p}\int_0^\infty\frac{du_3}{u_3^{1+p}}\int d^{2+2\eta}y\;
(\rho+\beta u_3)^{p-1}(\rho+\beta u_3+\beta b)^{-1-p}\big((T+1)|y-y_*|^2+(T+\beta)u_3\big)^{p-1}.
```

Combine the three factors with Feynman parameters $`\xi_1,\xi_2,\xi_3`$,
with exponents $`1-p`$, $`1+p`$, $`1-p`$. The quadratic form in $`y`$ is
$`\Lambda|y-\bar y|^2+(\xi_1+\xi_2)\xi_3Tb^2/\Lambda`$ with
$`\Lambda=\xi_1+\xi_2+\xi_3(T+1)`$, and the $`y`$ integral gives
$`\pi^{1+\eta}\Gamma(2-p-\eta)/\Gamma(3-p)\,\Lambda^{-1-\eta}M^{p-2+\eta}`$,
where $`M`$ is linear in $`u_3`$. The $`u_3`$ integral is a Beta function,
$`\Gamma(-p)\Gamma(2-\eta)/\Gamma(2-p-\eta)`$ times powers. The integrand is
homogeneous of degree $`-3`$ in $`\xi`$, so put $`\xi_1+\xi_2=1`$,
$`\xi=\xi_2`$, $`\lambda=\xi_3`$, and then
$`\mu=\lambda(T+1)/(1+\lambda(T+1))`$. Finally
$`\sin(\pi p)\Gamma(-p)=-\pi/\Gamma(1+p)`$. $`\blacksquare`$

The Feynman parameters here are an integration technique for the
transverse and $`u_3`$ integrals, applied after the contour has been closed
in position space. The cube has a direct reading:

| Limit | Meaning | Behaviour |
| --- | --- | --- |
| $`\beta\to0`$ | near end, the cusp $`x_2`$ | $`\beta^{\eta-1}`$ |
| $`b\to0`$ | far end, the cusp $`x_3`$ | $`b^{2\eta-1}`$ |
| $`\xi,\mu\to0`$ | collinear to the edge $`(2,3)`$ | pole $`1/\eta`$ |
| $`\mu\to1`$ | collinear to the edge $`(3,4)`$ | pole $`1/\eta`$ |

$`\mathcal J`$ converges and is analytic in $`p`$ for
$`-1<\operatorname{Re}p<1-\eta`$; from here on $`A`$ denotes the right side
on that strip. The integrand is positive for real $`p`$; for complex $`p`$
its modulus is the integrand at $`\sigma`$.

### 11.2 The remainder

**Theorem (edge piece).** Let $`p=\sigma+i\tau`$ with $`\sigma`$ in a compact
subset of $`(-1+\eta_0,1-\eta_0)`$ and $`0<\eta\le\eta_0<\tfrac12`$. Then

```math
\mathcal J=\frac{T^{p-1}\,\Gamma(1+p)\Gamma(1-p)}{\eta^2}
+\frac{(T+1)^{p-1}\varrho^{\eta-1+p}}{2\eta}\,\frac{\Gamma(1+p)\,\Gamma(1-p-\eta)\,\Gamma(\eta)^2}{\Gamma(2-\eta)\,\Gamma(2\eta)}
+R,\qquad|R|\le\frac{C\,(1+|\tau|)}{\eta},
```

with $`C`$ depending on the compact set, $`T`$ and $`\eta_0`$ only.
Equivalently

```math
A(p;T)=A_{\rm near}+A_{\rm far}+O\Big(\frac{(1+|\tau|)\,e^{2\pi|\tau|}}{\eta}\Big),
```

with the two terms of section 10, each equal to
$`-2\pi^2i\,\frac{\sin\pi p}{\pi p}\,\eta^{-2}`$ up to $`O(1/\eta)`$. The
bound is at fixed $`T`$ and is weak in $`\tau`$: the leading terms grow
like $`e^{\pi|\tau|}`$ only.

*Proof.* Write $`\mathcal J=(T+1)^{-1}\int_0^1b^{\eta-p}G(\beta)\,d\beta`$
with $`G`$ the integral over $`\xi,\mu`$, and split at $`\beta=\tfrac12`$.
All bounds are on moduli; differences of complex powers cost a factor
$`1+|\tau|`$.

*Near half, $`\beta\le\tfrac12`$.* Put $`\mu=\beta m/(\varrho b)`$,
$`0<m<M=\varrho b/\beta`$. Then $`\xi\beta+\varrho b\mu=\beta(\xi+m)`$,
$`(T+1)\beta+Tb\mu=(T+1)\beta(1+m)`$, and

```math
b^{\eta-p}G(\beta)=\beta^{\eta-1}\,\varrho^{p-1}(T+1)^p\,b^{\eta-1}\,\hat g(\beta),\qquad
\hat g=\int_0^1d\xi\int_0^Mdm\;k(\xi,m)\,\Big(1-\frac mM\Big)^{\eta-1},
```

```math
k=(1-\xi)^{-p}\xi^{p}m^{-p}(\xi+m)^{\eta-2}(1+m)^{p},\qquad g_0=\int_0^1d\xi\int_0^\infty dm\,k .
```

$`g_0`$ converges ($`k\sim m^{\eta-2}`$ at infinity) and
$`g_0=\Gamma(1+p)\Gamma(1-p)/\eta+O(1)`$: replacing
$`(1-\xi)^{-p}(1+m)^p`$ by 1 gives exactly
$`\int_0^1\xi^{\eta-1}d\xi\cdot B(1-p,1+p-\eta)`$, and the difference gains
a power of $`\xi`$ or $`m`$ and is bounded. Next,

```math
|\hat g-g_0|\le\int_0^1\!\!\int_M^\infty|k|+\int_0^1\!\!\int_0^M|k|\,\Big[\Big(1-\frac mM\Big)^{\eta-1}-1\Big].
```

The first term is at most $`CM^{\eta-1}`$. In the second, for $`m\le M/2`$
the bracket is at most $`4m/M`$, which gives $`C/M`$ from $`m\le1`$ and
$`CM^{\eta-1}/\eta`$ from $`m\ge1`$; for $`M/2<m<M`$,
$`|k|\le CM^{\eta-2}`$ and $`\int(1-m/M)^{\eta-1}dm\le M/\eta`$. Since
$`M\ge\varrho/(2\beta)`$,

```math
|\hat g-g_0|\le\frac C\eta\,\beta^{1-\eta}.
```

With $`b^{\eta-1}=1+O(\beta)`$ and $`g_0=O(1/\eta)`$,
$`b^{\eta-p}G=\beta^{\eta-1}\varrho^{p-1}(T+1)^pg_0+E`$ with
$`|E|\le C\eta^{-1}(1+\beta^{\eta})`$, and

```math
\frac1{T+1}\int_0^{1/2}b^{\eta-p}G\,d\beta=\frac{T^{p-1}\Gamma(1+p)\Gamma(1-p)}{\eta^2}+O\Big(\frac1\eta\Big).
```

*Far half, $`b\le\tfrac12`$.* Put $`\xi=\varrho bx/\beta`$,
$`0<x<X=\beta/(\varrho b)`$. Then $`\xi\beta+\varrho b\mu=\varrho b(x+\mu)`$ and

```math
b^{\eta-p}G(\beta)=\varrho^{\eta-1+p}\,\beta^{-1-p}\,b^{2\eta-1}\,h(\beta),\qquad
h=\int_0^Xdx\int_0^1d\mu\;x^{p}\mu^{-p}(1-\mu)^{\eta-1}(x+\mu)^{\eta-2}\Big(1-\frac xX\Big)^{-p}\big((T+1)\beta+Tb\mu\big)^{p}.
```

Let

```math
h_0=\int_0^\infty dx\int_0^1d\mu\;x^{p}\mu^{-p}(1-\mu)^{\eta-1}(x+\mu)^{\eta-2}
=B(1+p,1-p-\eta)\,B(\eta,\eta)
=\frac{\Gamma(1+p)\Gamma(1-p-\eta)}{\Gamma(2-\eta)}\cdot\frac{\Gamma(\eta)^2}{\Gamma(2\eta)} ,
```

by the $`x`$ integral first, which gives $`\mu^{p+\eta-1}`$. This is the
far-end coefficient of section 10 obtained independently:
$`B(\eta,\eta)`$ is the two collinear limits $`\mu\to0`$ and $`\mu\to1`$.
The differences between $`h`$ and $`(T+1)^ph_0`$ are: the tail $`x>X`$, at
most $`C\eta^{-1}X^{\sigma+\eta-1}`$; the factor $`(1-x/X)^{-p}-1`$, at most
$`C(1+|\tau|)\eta^{-1}X^{-\theta}`$ with
$`\theta=\min(1,1-\sigma-\eta)-\epsilon`$ for a small fixed $`\epsilon>0`$
(for $`x\le X/2`$ it is at most $`C|p|x/X`$, for $`x>X/2`$ the integrand is
at most $`CX^{\sigma+\eta-2}`$ times an integrable function); and
$`((T+1)\beta+Tb\mu)^p-(T+1)^p`$ and $`\beta^{-1-p}-1`$, each $`O(|p|b)`$
times a quantity of order $`1/\eta`$. In each case the $`\mu`$ integral
contributes $`B(1-\sigma,\eta)\le C/\eta`$. As $`X\ge1/(2\varrho b)`$,

```math
b^{\eta-p}G=\varrho^{\eta-1+p}(T+1)^ph_0\,b^{2\eta-1}+E',\qquad|E'|\le\frac{C(1+|\tau|)}{\eta}\,b^{2\eta-1+\theta},
```

and $`\int_0^{1/2}b^{2\eta-1}db=2^{-2\eta}/(2\eta)`$ with
$`h_0=O(1/\eta)`$ gives the second term of the theorem plus
$`O(1/\eta)`$. $`\blacksquare`$

### 11.3 The other edge and the full integral

The relabelling $`x_1\leftrightarrow x_2`$, $`x_3\leftrightarrow x_4`$
exchanges $`s`$ and $`t`$, sends $`z\to1/z`$, and reverses $`p_1`$. It
preserves the foliation by null hyperplanes orthogonal to $`p_1`$ and
exchanges the two parts of the slab. So the piece for the edge $`(4,1)`$ is
the piece for $`(2,3)`$ of the relabelled square: in the present units it
is $`T^{\eta}A(-p;1/T)`$, and

```math
\text{full integral}=A(p;T)+T^{\eta}\,A(-p;1/T).
```

**Corollary (double pole from the Lorentzian domain).** By the lemma
below, for the iterated integral defined at the start of this section,
spacelike diagonals and the weight $`z^p`$,
$`|\operatorname{Re}p|<1-\eta`$,

```math
st\int\frac{z^p}{x_{01}^2x_{02}^2x_{03}^2x_{04}^2}
=-2\pi^2i\,\frac{4}{\eta^2}\,\frac{\sin\pi p}{\pi p}+O\Big(\frac1\eta\Big):
```

two units from each edge piece, one from each of its ends.

**What the corollary uses.** The proposition identifies $`A`$ with the
position-space edge piece for $`\operatorname{Re}p<0`$. The second piece
needs $`A(-p)`$, that is $`\operatorname{Re}p>0`$. So at a given $`p`$ only
one of the two pieces is identified directly, and the corollary needs the
following lemma for the other. Without it the statement is: each edge
piece carries two units on its own half strip.

**Lemma (continuation).** The position-space edge piece is analytic in
$`p`$ for $`-1<\operatorname{Re}p<1-\eta`$ and equals $`A(p;T)`$ there.

*Proof.* At fixed $`(\beta,y)`$ write the keyhole integral as
$`\oint_Hu_3^{-1-p}\varphi(u_3)\,du_3`$ times a phase, with

```math
\varphi(u_3)=(\rho+\beta u_3)^{p-1}(\rho+\beta u_3+\beta b)^{-1-p}\big((T+1)|y-y_*|^2+(T+\beta)u_3\big)^{p-1},
```

which is analytic near $`u_3\ge0`$ for $`\rho>0`$, $`y\ne y_*`$. Integrate
by parts on $`H`$. The ends at infinity do not contribute, since
$`u_3^{-p}\varphi\sim u_3^{-3}`$. So for $`p\ne0`$

```math
\oint_Hu_3^{-1-p}\varphi\,du_3=\frac1p\oint_Hu_3^{-p}\varphi'\,du_3,
```

and for $`\operatorname{Re}p<1`$ the circle on the right shrinks:

```math
\text{keyhole}=\frac{2i\sin\pi p}{p}\int_0^\infty u_3^{-p}\,\varphi'(u_3)\,du_3 .
```

For $`\operatorname{Re}p<0`$ this is the formula used in the proposition,
by an ordinary integration by parts. $`\varphi'`$ is a sum of three terms,
each $`\varphi`$ with one exponent lowered by one and a factor
$`(p-1)\beta`$, $`-(1+p)\beta`$ or $`(p-1)(T+\beta)`$. Take moduli and
apply to each term the reduction of the proposition: Feynman parameters
with exponents $`(2-\sigma,1+\sigma,1-\sigma)`$, $`(1-\sigma,2+\sigma,1-\sigma)`$
or $`(1-\sigma,1+\sigma,2-\sigma)`$, the $`y`$ integral, and
$`\int_0^\infty u_3^{-\sigma}(M_0+c_1u_3)^{\sigma-3+\eta}du_3=M_0^{\eta-2}c_1^{\sigma-1}B(1-\sigma,2-\eta)`$,
finite for $`\sigma<1`$. Here
$`c_1=\beta(\xi_1+\xi_2)+\xi_3(T+\beta)`$. The result is the integrand of
$`\mathcal J`$ at $`\sigma`$, in its homogeneous form, times
$`\beta\xi_1/c_1`$, $`\beta\xi_2/c_1`$ or $`(T+\beta)\xi_3/c_1`$
respectively. These three factors are non-negative and sum to 1. Since
$`\mathcal J`$ converges for $`-1<\sigma<1-\eta`$, the integral of
$`|u_3^{-p}\varphi'|`$ over $`(\beta,y,u_3)`$ is finite. It is locally
bounded in $`p`$: the modulus of the integrand is the exponential of a
function linear in $`\sigma`$, hence convex in $`\sigma`$ and bounded on a
compact interval by the sum of its values at the two ends. By Tonelli and
Fubini the iterated integral in the defined order exists and equals this
absolutely convergent integral (the null set $`\rho=0`$ or $`y=y_*`$ is
excluded), and an integral of a function holomorphic in $`p`$ with a
locally bounded majorant is holomorphic. At $`p=0`$ the keyhole at fixed
$`(\beta,y)`$ is entire in $`p`$, so the formula extends by continuity; it
gives the residue. The piece equals $`A`$ for
$`\eta-1<\operatorname{Re}p<0`$, and $`A`$ is analytic on the strip, so
they agree on the strip. $`\blacksquare`$

With the lemma the corollary holds for $`|\operatorname{Re}p|<1-\eta`$.
The lemma was added after the first independent review of this section and
was then refereed separately: the integration by parts, the phases, the
three reductions and the logic were re-derived and found correct, and the
keyhole identity was confirmed numerically for $`\operatorname{Re}p>0`$ to
between $`10^{-14}`$ and $`10^{-27}`$. Independently of it, the identity
$`A(p;T)+T^\eta A(-p;1/T)=`$ (the box) is confirmed numerically to 12
digits below, where the box is computed by Mellin–Barnes.

### 11.4 Numerical checks

`edge_cube_integral`, `edge_piece`, `box_from_two_edges`. The box value is
$`-2i\pi^{2+\eta}\,T^{\eta}N_\eta(p)/\eta^2`$ with $`N_\eta`$ the
Mellin–Barnes form at $`s/t=1/T`$.

| $`p`$ | $`\eta`$ | $`T`$ | (two edge pieces)/(box) |
| --- | --- | --- | --- |
| 0.30 | 0.25 | 1 | 1.000000000000 |
| 0.15 | 0.10 | 1 | 1.000000000000 |
| 0.40 | 0.05 | 1 | 1.000000000000 |
| 0.30 | 0.20 | 3 | 1.000000 after the factor $`T^\eta`$ |
| $`-0.25`$ | 0.10 | 0.4 | 1.000000 after the factor $`T^\eta`$ |
| 0.35 | 0.10 | 2.5 | 1.000000 after the factor $`T^\eta`$ |

The remainder, $`\eta\,(A-A_{\rm near}-A_{\rm far})/(-2\pi^2i)`$, with the
closed forms of section 10:

| $`p`$ | $`T`$ | $`\eta=0.2`$ | 0.1 | 0.05 | 0.025 | 0.0125 |
| --- | --- | --- | --- | --- | --- | --- |
| $`-0.3`$ | 1 | $`-0.584`$ | $`-0.055`$ | 0.177 | 0.287 | 0.341 |
| 0.3 | 1 | $`-0.551`$ | $`-0.023`$ | 0.198 | 0.299 | 0.347 |
| $`-0.3`$ | 3 | $`-0.825`$ | $`-0.346`$ | $`-0.142`$ | $`-0.046`$ | 0.000 |
| 0.5 | 0.4 | $`-0.873`$ | $`-0.183`$ | 0.080 | 0.197 | 0.253 |
| $`-0.6`$ | 2 | $`-1.885`$ | $`-0.958`$ | $`-0.671`$ | $`-0.551`$ | $`-0.495`$ |

It converges in every row: the remainder is a single pole.

### 11.5 What is and is not established

Established, for the weight $`z^p`$, spacelike diagonals, and the iterated
integral in the order $`\alpha`$, $`y_\perp`$, $`\beta`$ with
$`\varepsilon\to0`$ innermost:

- the integral is exactly two edge pieces (section 9);
- each piece is given by the cube integral, positive for real $`p`$, on
  the whole strip $`|\operatorname{Re}p|<1-\eta`$ (proposition and
  continuation lemma);
- the cube function equals its two end terms up to a single pole, with a
  bound uniform in $`\eta`$ at fixed $`p`$ and $`T`$ (theorem);
- each end term carries $`\sin\pi p/(\pi p)`$ (sections 8 and 10);
- hence the double pole of the full integral is
  $`4\sin\pi p/(\pi p)`$ (corollary).

No Mellin–Barnes representation and no global contour is used. An
independent review re-derived sections 9 to 11 and, separately, the
continuation lemma, found no wrong sign, constant or exponent, and
confirmed with its own code the two-edge theorem, the reduction at fixed
$`\beta`$ against a three-dimensional position-space integral (ratio
1.0000001), the remainder table with a fourth row at $`p=0.7`$, and the
box identity at five points.

Not established:

- (the benchmark weights are treated in section 12;)
- a timelike diagonal: the table of half planes in section 9 and the
  positivity used here change;
- the limit $`\varepsilon\to0`$ taken after the outer integrals, and other
  orders of integration;
- the equality of this iterated integral with the integral over the global
  contour of GATE_6_CONTOUR.

## 12. General weights from the Lorentzian domain

**Theorem (benchmark weights).** Let both diagonals be spacelike, and let
$`f(L)=F(e^L)`$ be holomorphic on $`|\operatorname{Im}L|<W`$ with
$`W>2\pi`$ and $`|f(L)|\le A(1+|L|)^k`$ there. For the iterated integral of
section 11 (the $`\alpha`$ integral first, the displacement removed at
fixed $`(\beta,y_\perp)`$, then $`y_\perp`$, then $`\beta`$),

```math
st\int\frac{F(z)}{x_{01}^2x_{02}^2x_{03}^2x_{04}^2}=-2\pi^2i\,\frac{4\,M[F]}{\eta^2}+O\Big(\frac1\eta\Big),\qquad\eta\to0^+ .
```

The hypothesis holds for $`\bar F^{(0)}`$ and $`\bar F^{(1)}`$ with any
$`W`$ and for $`\bar F^{(2)}`$ with $`W=3\pi-\delta`$ (GATE_6_BOUNDS,
Lemma 3′).

*Proof.* By GATE_6_BOUNDS Lemma 3, for $`0<\sigma<1/m`$,
$`f(L)=\int_{\operatorname{Re}p=-\sigma}c_+(p)e^{pL}dp+\int_{\operatorname{Re}p=\sigma}c_-(p)e^{pL}dp`$
on the strip, with $`|c_\pm|\le Ce^{-(W-\delta')|\tau|}`$; choose
$`\delta'`$ with $`W-\delta'>2\pi`$. Consider the edge piece for the edge
$`(2,3)`$.

*Step 1, the $`\alpha`$ integral.* Fix $`(\beta,y)`$ with $`\rho>0`$,
$`y\ne y_*`$, and $`\varepsilon>0`$. On the real $`\alpha`$ line
$`|\operatorname{Im}L|\le2\pi`$ (section 9), so
$`|z^p|\le|z|^{\pm\sigma}e^{2\pi|\tau|}`$; the integrand has no zero on the
line and decays like $`|\alpha|^{-4}`$. With the decay of $`c_\pm`$ the
double integral over $`\alpha`$ and $`p`$ converges absolutely, and the
order can be exchanged. For each $`p`$ the $`\alpha`$ integral is the
keyhole of section 9. To take $`\varepsilon\to0`$ under the $`p`$ integral
use a contour that does not move with $`\varepsilon`$: in the variable
$`u_3`$, two lips $`\operatorname{Im}u_3=\pm\delta_0`$,
$`\operatorname{Re}u_3\ge0`$, joined by a semicircle on the left. For
$`\rho>0`$, $`y\ne y_*`$ the other three zeros have
$`\operatorname{Re}u_3<0`$ at a distance independent of
$`\varepsilon`$, so for small $`\delta_0`$ the contour encloses only the
branch point. On it $`\operatorname{Im}L=\mp\pi+O(\delta_0)`$ on the lips
and interpolates on the arc, $`|z|^{\pm\sigma}`$ is bounded and the tails
decay like $`|\alpha|^{-4}`$, all uniformly in small $`\varepsilon`$. So the
integral is bounded by $`C(\beta,y)e^{2\pi|\tau|}`$ uniformly, and by
dominated convergence

```math
\lim_{\varepsilon\to0}\int d\alpha\,\frac{F(z)}{\prod x_{0i}^2}=\sum_\pm\int c_\pm(p)\,K_p(\beta,y)\,dp ,
```

with $`K_p`$ the keyhole integral for $`z^p`$ at $`\varepsilon=0`$.

*Step 2, the outer integrals.* By the proof of the continuation lemma,
$`\int d\beta\int d^{2+2\eta}y\,|K_p|\le C(1+|\tau|)\,|\sin(\pi p)/p|\le C(1+|\tau|)e^{\pi|\tau|}`$
on the two lines. So the $`p`$ integral can be exchanged with the $`y`$ and
$`\beta`$ integrals, and the edge piece for $`F`$ is

```math
A[F]=\sum_\pm\int c_\pm(p)\,A(p;T)\,dp .
```

*Step 3, the limit.* By the theorem of section 11,
$`A(p;T)=A_{\rm near}(p)+A_{\rm far}(p)+R_p`$ with
$`|R_p|\le C(1+|\tau|)e^{2\pi|\tau|}/\eta`$ on the two lines, uniformly in
$`\eta\le\eta_0`$. Since $`W-\delta'>2\pi`$, $`\int|c_\pm||R_p|\le C/\eta`$.
The near term is exactly $`-2\pi^{2+\eta}i\,\eta^{-2}\,\mathcal T_\eta(p)`$.
Section 7 needs $`|f|+|f'|`$ polynomially bounded on
$`|\operatorname{Im}L|\le\pi`$, which follows from the hypothesis by
Cauchy's estimate on discs of radius $`(W-\pi)/2`$. So, with
$`\pi^\eta=1+O(\eta)`$,

```math
\sum_\pm\int c_\pm A_{\rm near}=-\frac{2\pi^{2+\eta}i}{\eta^2}\,\mathcal T_\eta[F]=-\frac{2\pi^2i}{\eta^2}\big(M[F]+O(\eta)\big).
```

The far term is
$`-2\pi^{2+\eta}i\,\varrho^\eta\eta^{-2}\,\Phi_\eta(p)`$ with
$`\Phi_\eta(p)=\frac{\sin\pi p}{\pi p}\cdot\frac{\eta^2\Gamma(\eta)^2}{\Gamma(1+2\eta)}\cdot\frac{\Gamma(1-p-\eta)}{\Gamma(1-p)}`$.
On the lines the last two factors differ from 1 by at most
$`C\eta(1+\log(2+|\tau|))`$, and $`\int c(p)\frac{\sin\pi p}{\pi p}dp=M[F]`$
by Lemma 3 on the unit circle. So the far term also gives
$`-2\pi^2i\,\eta^{-2}(M[F]+O(\eta))`$, and

```math
A[F]=-2\pi^2i\,\frac{2\,M[F]}{\eta^2}+O\Big(\frac1\eta\Big).
```

*Step 4, the other edge.* The relabelling of section 11.3 is a change of
variables that does not depend on the weight, so steps 1 and 2 apply to
the piece for the edge $`(4,1)`$: it is the first piece with the diagonals
exchanged, the weight $`F(1/z)`$ and a factor $`T^\eta=1+O(\eta)`$. That
weight is $`f(-L)`$, with the same strip, bound and mean, and section 11
applies at $`1/T`$. $`\blacksquare`$

The constants depend on $`T`$ and are not uniform as $`T\to0`$ or
$`\infty`$. The statement is a bound as $`\eta\to0^+`$; it is not a
statement about the pole structure in $`\eta`$.

**What this replaces.** For spacelike diagonals, and for the iterated
integral, this derives the coefficient $`4M[\bar F^{(L-1)}]`$ for the
benchmark weights without the global contour and its certificate, without
Mellin–Barnes and without the phase bound on the contour. It uses the
contour closing in $`\alpha`$, elementary integrals, and Lemma 3. That the
iterated integral equals the integral over the global contour, or the
parametric integral, is not proved; for $`z^p`$ the two agree numerically
(section 11.4).

**Numerical check.** For $`F=\log^2z+\pi^2`$ the edge piece is
$`\partial_p^2A(p)|_{p=0}+\pi^2A(0)`$. With `edge_piece`,
$`\eta^2A[F]/(-2\pi^2i)`$ at $`\eta=0.2,0.1,0.05,0.025`$ is
17.88, 15.70, 14.48, 13.83 at $`T=1`$ and 18.10, 15.76, 14.50, 13.84 at
$`T=2.5`$, against $`2M[F]=13.16`$; the difference divided by $`\eta`$ is
between 24 and 27.

**Remarks on what is still fixed by the definition.**

- *The other null direction.* The cyclic relabelling $`x_i\to x_{i+1}`$
  exchanges the roles of $`\alpha`$ and $`\beta`$, sends $`z\to1/z`$ and
  exchanges $`s`$ and $`t`$. So the theorem holds equally with the $`\beta`$
  integral first; the two pieces are then the edges $`(1,2)`$ and
  $`(3,4)`$.
- *$`\alpha`$ and $`\rho`$ in the cusp-local model.* With $`\rho`$ rotated
  by a phase the model integrand is absolutely integrable in
  $`(\alpha,\rho)`$ at fixed $`\beta`$: it decays like the inverse cube at
  infinity and is at most $`|(\alpha,\rho)|^{-2}`$ at the one point where
  two invariants vanish, against the measure $`\rho^\eta d\rho\,d\alpha`$;
  at the isolated zero of $`x_{01}^2`$ it is of order one over the distance,
  which is integrable for the weights considered. So in the model those two
  integrals may be exchanged, with $`\rho`$ on the rotated ray. With the
  far propagator the rotation is not available and this is not shown.
- *The edge statement as a formula.* By the previous remark the collinear
  pole of the model at fixed $`\beta`$ can be computed with $`\rho`$ first.
  The $`\rho`$ integral then behaves like $`|\alpha|^{\eta-1}`$ times the
  two one-sided limits of GATE_6_CLOSURE, and the pole coefficient is their
  difference, the edge jump $`-2\pi i\,M[F]`$; with $`\alpha`$ first it is
  the coefficient $`2\pi i\,M[F]/\beta`$ of section 8. If the scaling limit
  in the $`\rho`$ integral is justified, this identifies the edge jump with
  the residue of the collinear pole. That limit is not written out, so this
  is an outline; for $`z^p`$ the two one-sided integrals have been checked
  in closed form to give $`2\pi i\,\frac{\sin\pi p}{\pi p}(1-\beta)^{-p}/\beta`$.

Not treated: a timelike diagonal; removing the displacement after the
outer integrals; the order with the transverse integral first for the full
integrand; and the equality of the iterated integral with the integral
over the global contour.

An independent review re-derived this section, found no wrong constant,
sign or exponent, and asked for the three additions now made (the contour
uniform in $`\varepsilon`$, the Cauchy estimate, the transfer of the
relabelling) and for the qualifier on "What this replaces". It reproduced
the numerical check with an independent three-dimensional quadrature:
17.8808, 15.7035, 14.4813 at $`T=1`$.
