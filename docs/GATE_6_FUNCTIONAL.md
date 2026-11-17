# Gate 6 at the level of functionals

**Status.** The scalar form of Gate 6 holds on the Mellin family $`z^p`$,
$`|\Re p|<1`$, and on every polynomial in $`\log z`$. This covers
$`\bar F^{(0)}`$ and $`\bar F^{(1)}`$ exactly. The chain-level identity in
[GATE_6_PROOF.md](GATE_6_PROOF.md) is a stronger statement and stays open.
The argument below is the standard Mellin derivation of the contour rule in
[Arkani-Hamed, Henn and Trnka, arXiv:2112.06956](https://arxiv.org/abs/2112.06956),
written out in this repository's conventions and checked numerically in
[the Mellin example](../examples/mellin_functional.py).

## 1. Three linear functionals

For a finite weight $`F(z)`$ define

1. the circle mean $`M[F]=(2\pi)^{-1}\int_{-\pi}^{\pi}F(e^{i\phi})\,d\phi`$
   on the principal sheet;
2. the geometric period $`R_B[F]=(2\pi i)^{-3}\int_{C}F(z)\,\omega`$, where
   $`\omega=2\,da\wedge du\wedge dv/[a(a-uv)]`$ is the weighted-divisor
   residue and $`C`$ is the capped cycle of example 36, with $`z=uv/a-1`$;
3. the Lorentzian double pole $`L[F]`$, the coefficient of $`4/\epsilon^2`$
   in $`x_{13}^2x_{24}^2\int d^Dx_0\,F(z)/(i\pi^{D/2}\prod_j x_{0j}^2)`$
   with $`D=4-2\epsilon`$ and
   $`z=x_{02}^2x_{04}^2x_{13}^2/(x_{01}^2x_{03}^2x_{24}^2)`$.

All three are linear in $`F`$. The claim is $`R_B=M=L`$ on the stated family.

## 2. Geometric period equals the circle mean

On the cycle $`a=Re^{i\alpha}`$, $`z=e^{i\phi}`$,
$`u=\sqrt{|a(z+1)|}\,e^{i\theta}`$, $`v=a(z+1)/u`$, one has $`a-uv=-az`$.
The radial derivative of $`u`$ cancels in the Jacobian, leaving
$`\partial(a,u,v)/\partial(\alpha,\phi,\theta)=i\,a^2z`$. The pulled-back
coefficient is therefore the constant $`-2i`$. With $`0\le\alpha<2\pi`$,
$`-\pi<\phi<\pi`$ and $`0\le\theta<\pi`$,

```math
R_B[F]=\frac{-2i\cdot2\pi\cdot\pi}{(2\pi i)^3}\int_{-\pi}^{\pi}F(e^{i\phi})\,d\phi=M[F].
```

The example integrates $`\omega`$ on the cycle with a finite-difference
Jacobian of the parameterization, so it does not assume this reduction.
In particular $`R_B[z^p]=\sin(\pi p)/(\pi p)`$.

## 3. Lorentzian double pole of $`z^p`$

Multiplying by $`z^p`$ shifts the propagator powers to $`1+p`$ on legs 1, 3
and $`1-p`$ on legs 2, 4. Their sum is 4, so with $`s=x_{13}^2`$,
$`t=x_{24}^2`$ the Feynman-parameter form is

```math
\frac{\Gamma(2+\epsilon)}{\Gamma(1+p)^2\Gamma(1-p)^2}
\int d^4\alpha\,\delta(1-\textstyle\sum\alpha)\,
\frac{(\alpha_1\alpha_3)^{p}(\alpha_2\alpha_4)^{-p}}
{(\alpha_1\alpha_3s+\alpha_2\alpha_4t)^{2+\epsilon}} .
```

Splitting the denominator with one Mellin–Barnes variable $`w`$ and doing
the parameter integrals gives

```math
\frac{1}{\Gamma(1+p)^2\Gamma(1-p)^2\,\Gamma(-2\epsilon)}
\int\frac{dw}{2\pi i}\,
\Gamma(-w)\Gamma(2+\epsilon+w)\Gamma(1+p+w)^2\Gamma(-1-p-\epsilon-w)^2\,
s^{w}t^{-2-\epsilon-w}.
```

The contour separates the double pole at $`w=-1-p`$ from the double pole at
$`w=-1-p-\epsilon`$. It exists for $`\epsilon<0`$, the infrared side. As
$`\epsilon\to0`$ the two pinch. The residue at the right pole has leading
term $`2g(0)/\epsilon^3`$ with
$`g(0)=\Gamma(1+p)\Gamma(1-p)\,s^{-1-p}t^{-1+p}`$. Closing to the right
contributes a minus sign and $`1/\Gamma(-2\epsilon)=-2\epsilon+O(\epsilon^2)`$.
After the factors $`st`$ and $`(s/t)^p`$ from the definition of $`z`$,

```math
\epsilon^2\,x_{13}^2x_{24}^2\!\int\! z^p\;\longrightarrow\;
\frac{4}{\Gamma(1+p)\Gamma(1-p)}=4\,\frac{\sin\pi p}{\pi p},
\qquad\text{so}\qquad L[z^p]=\frac{\sin\pi p}{\pi p}=M[z^p].
```

At $`p=0`$ this is the massless box and reproduces its $`4/\epsilon^2`$
and $`-2\log(st)/\epsilon`$ terms. The example evaluates the
Mellin–Barnes integral at $`\epsilon=-0.1,\dots,-0.002`$ and recovers the
closed form to $`2\times10^{-4}`$ after one Richardson step.

## 4. Extension and scope

Both sides are analytic in $`p`$ on $`|\Re p|<1`$. Differentiating at
$`p=0`$ gives $`L[\log^kz]=M[\log^kz]`$ for every $`k`$. Hence
$`L[\bar F^{(0)}]=-1`$ and $`L[\bar F^{(1)}]=2\pi^2/3`$, and
$`\Gamma^{(L)}=-4M_L/L`$ follows in the normalization of
[GATE_6_PROOF.md](GATE_6_PROOF.md), section 1. For a general $`F`$ the
identity follows wherever $`F`$ has a Mellin representation
$`F(z)=\int c(p)\,z^p\,dp`$ on a contour inside the strip and the
$`p`$ integral commutes with the regulated loop integral.

What this does not establish:

- It compares numbers produced by two functionals. It gives no chain whose
  boundary is the difference of the capped cycle and the Feynman contour.
  The rim collisions and relative faces of GATE_6_PROOF.md are untouched.
  [GATE_6_CLOSURE.md](GATE_6_CLOSURE.md) supplies such a chain on the fibre
  of the weighted divisor by a different route.
- It uses dimensional regularization. The finite-frame $`i\epsilon P^2`$
  regulator of the examples is tied to it only through regulator
  independence of the leading double logarithm.
- $`F`$ is taken in four dimensions. Its $`O(\epsilon)`$ corrections affect
  subleading poles only.
- The weight $`F`$ itself still comes from the Lorentzian integral. The
  geometry supplies the cycle and the functional, not $`F`$.

## 5. Solver

`geocusp.solver` evaluates four functionals for a weight
$`F=\sum_k c_k\log^kz+Az^p`$: the circle mean, the capped-cycle period, the
Mellin–Barnes double pole, and the Lorentzian edge chain of
[GATE_6_CLOSURE.md](GATE_6_CLOSURE.md). It also handles the two-loop weight
$`\bar F^{(2)}`$, surveys the integration contour of
[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md), integrates the CD weight at a complex
$`z`$, and reports the roots and meridian periods on a rim slice. It runs
locally in one process with no server and no network access.

```bash
PYTHONPATH=src python -m geocusp.solver                 # desktop window
PYTHONPATH=src python -m geocusp.solver functional --logs 9.8696044 0 1 --loop-order 2
PYTHONPATH=src python -m geocusp.solver edge --cusp 3 --edge prev --position 0.7 --region +-
PYTHONPATH=src python -m geocusp.solver edge --two-loop --position 0.4
PYTHONPATH=src python -m geocusp.solver two-loop        # about a minute; --no-edge is quick
PYTHONPATH=src python -m geocusp.solver contour --points 5000
PYTHONPATH=src python -m geocusp.solver cd 0.35+0.8j
PYTHONPATH=src python -m geocusp.solver rim --t 1 --lam 470 --s 0.27
```

`edge --position x` evaluates the collinear singularity at a point of the
null line through a cusp; the edge segment is `0 < x < 1`. The window has a
"Lorentzian edge chain" tab with the same controls, a scan along the line,
and a button for the two-loop weight.
