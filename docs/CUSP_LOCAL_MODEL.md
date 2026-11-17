# The cusp-local model and a proof of the double pole without Mellin–Barnes

The question: does the local structure at a cusp determine the double pole
of the whole integral directly?

**Answer.** Partly. Two things are proved here, for spacelike diagonals and
the weight $`z^p`$, $`|\operatorname{Re}p|<1`$:

1. The cusp-local model, in which the far propagator is frozen at its cusp
   value, is a triangle with two light-like sides. It has a closed form. Its
   double pole is the product of the collinear poles of the two edges at the
   cusp, with coefficient $`\sin\pi p/(\pi p)=M[z^p]`$.
2. The full parametric integral equals the sum of the four cusp-local
   models up to a single pole. So the double pole is
   $`4\sin\pi p/(\pi p)`$, one unit of $`M[z^p]`$ per cusp. The proof is
   elementary real analysis on the Feynman simplex. It does not use
   Mellin–Barnes.

What is **not** proved: the same statement carried out entirely in position
space, with the Lorentzian paths of GATE_6_CLOSURE as the source of
$`M[F]`$ and no Mellin transform in $`z`$. Section 5 says where that stands.

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

uniformly for $`\sigma`$ in a compact subset of $`(-1,1)`$. This replaces
Lemma 4, Lemma 13 and the residue computation of GATE_6_FUNCTIONAL in the
chain of GATE_6_BOUNDS section 5 when both diagonals are spacelike. The
growth $`e^{2\pi|\tau|}`$ is worse than Lemma 4's $`e^{\pi|\tau|}`$ and is
still enough: steps C and D there need $`W-\delta'>2\pi`$, and
$`W>\Phi\ge2\pi`$ holds already.

With a timelike diagonal $`A`$ changes sign on the simplex and this proof
does not apply; Lemma 14 and the Mellin–Barnes form are still needed there.

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
  $`x_j`$**, since the vertex model is the parametric form of the triangle
  at that cusp. This identification is by comparison of the two formulas.
- **Not done: a position-space proof with the Lorentzian paths.** In the
  exact variables of section 1 the local integral is
  $`\int\frac{du_1}{u_1}\frac{du_3}{u_3}|u_1u_3|^\eta\,G(u_1,u_3)`$, and its
  double pole is the sum over the four sign quadrants of the limits of
  $`G`$, which is the quadrant sum $`2\pi i\,M[F]`$ of GATE_6_CLOSURE. Two
  obstacles remain. On the real section the iterated integrals are not
  absolutely convergent, so the order cannot be exchanged without the
  contour. And for a weight growing like $`\log^kz`$ each quadrant has a
  $`\log^{k+1}`$ that cancels only in the sum, while the regulator differs
  between quadrants at $`|z|\sim1/|u|`$. A proof has to show that this
  mismatch does not promote to a higher pole. That is the open theorem.
