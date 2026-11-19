# The three-loop weight

The results in this repository were developed on three weights whose
coefficients were known in advance. This note records one test beyond
them, and says plainly what it is and is not.

## The function

Henn, Korchemsky and Mistlberger (arXiv:1911.10174) give the function
$`F`$ to three loops, as harmonic polylogarithms of $`x=z`$ with letters
$`0`$ and $`-1`$, of weight up to six. At three loops there is a planar
part and a non-planar part,

```math
F^{(3)}=F^{(3)}_{\rm P}+\frac{1}{N^2}\Big(\frac{1}{1+x}\,F^{(3)}_{\rm NP1}+F^{(3)}_{\rm NP2}\Big).
```

In the normalisation used here $`\bar F^{(L)}=4^{L+1}F^{(L)}`$. The one-
and two-loop expressions are in `geocusp.hpl`. The three-loop files are
the paper's ancillary files and are not distributed here; with a copy,

```bash
PYTHONPATH=src python -m geocusp.hpl /path/to/anc
HKM_ANC=/path/to/anc python -m pytest tests/test_hpl.py
```

`geocusp.hpl` evaluates the harmonic polylogarithms by a series near the
origin and continuation along a path, on any sheet.

## The circle mean

The rule $`\Gamma^{(L)}=-(4/L)\,M[\bar F^{(L-1)}]`$, with $`M`$ the mean
over the unit circle on the principal sheet, gives, in powers of
$`\alpha_s/\pi`$ times the colour factor:

| Order | From the circle mean | Known value | Difference |
| --- | --- | --- | --- |
| Two loops | $`-0.8224670334241`$ | $`-\zeta_2/2`$ | $`10^{-14}`$ |
| Three loops | $`+1.4881944463528`$ | $`11\zeta_4/8`$ | $`3\times10^{-14}`$ |
| Four loops, planar | $`-3.6618383900323`$ | $`-(\zeta_3^2/8+219\zeta_6/64)`$ | $`6\times10^{-14}`$ |
| Four loops, non-planar | $`-9.1644233443`$ | $`-\tfrac32(\tfrac32\zeta_3^2+\tfrac{31}{8}\zeta_6)`$ | $`5\times10^{-11}`$ |

The planar value is $`-8(73\pi^6/630+4\zeta_3^2)`$ in powers of $`g^2`$,
the known four-loop coefficient. The non-planar value is the coefficient
of $`N^2`$ that comes from the quartic Casimir term of that paper.

This is a check, not a prediction: the four-loop values are published,
and the circle rule itself is from the literature. What it adds is one
more order, a function of weight six with 46 harmonic polylogarithms in
place of a short expression, and a non-planar piece.

## Does the theorem apply to it?

The double-pole theorem (GATE_6_BOUNDS, steps A to D, and
CUSP_LOCAL_MODEL section 12) needs the weight, as a function of
$`L=\log z`$, to be holomorphic on a strip $`|\operatorname{Im}L|<W`$ wider
than the range of $`\arg z`$ on the contour, with polynomial growth. For
the three-loop weight this was not obvious in advance, for two reasons:
the only possible singular points in the strip are
$`L=\pm i\pi`$, that is $`z=-1`$, which is exactly where the two sides of
an edge sit; and the non-planar part carries an explicit factor
$`1/(1+x)`$.

It does apply. An earlier version of this note had numerical evidence
only (the monodromy around $`x=-1`$ to 43 digits). The statement is now
proved, by an exact computation. It is a computer-assisted proof: it
rests on the module described below and on the published expressions as
transcribed from the paper's files.

**Proposition (holomorphy at $`z=-1`$).** Let
$`F(x)=\sum_wc_wH_w(x)`$ be a finite combination of harmonic
polylogarithms with letters $`0`$ and $`-1`$. Put $`y=1+x`$ and write
$`\omega_0=dy/y`$, $`\omega_1=dy/(y-1)`$. Words are written as the
indices of $`H`$ are, the outermost integration first. For a word $`v`$
in $`\omega_0,\omega_1`$ let

```math
T_\pm(v)=\sum_{v=v_2\,\omega_1^{m}}(-1)^{n_0(v_2)}\,\zeta^{\rm reg}(\bar v_2)\,\frac{(\pm i\pi)^m}{m!},
```

where the sum is over the ways of writing $`v`$ as a word $`v_2`$
followed by $`m\ge0`$ letters $`\omega_1`$, $`n_0`$ counts the letters
$`\omega_0`$, $`\bar v_2`$ is $`v_2`$ with $`\omega_1\to du/u`$ and
$`\omega_0\to du/(1-u)`$, and $`\zeta^{\rm reg}`$ is the shuffle-regularised
multiple zeta value. Let $`\hat w`$ be $`w`$ with $`0\to\omega_1`$,
$`-1\to\omega_0`$, and

```math
d^\pm_u=\sum_{w:\ \hat w=u\,v}c_w\,T_\pm(v).
```

If $`d^\pm_u=0`$ for every word $`u`$ whose last letter is $`\omega_0`$,
then the continuation of $`F`$ from $`x>0`$ through the upper (lower)
half plane to a neighbourhood of $`x=-1`$ extends to a holomorphic
function on $`|1+x|<1`$, with value $`d^\pm_\emptyset`$ at $`x=-1`$ and
derivative $`-\sum_{k\ge0}d^\pm_{\omega_0^k\omega_1}`$.

*Proof.* Since $`dx/x=\omega_1`$ and $`dx/(1+x)=\omega_0`$, $`H_w(x)`$ is
the iterated integral of $`\hat w`$ from the tangential base point at
$`y=1`$ (the regularisation that makes $`H_{0\dots0}=\log^nx/n!`$). The
path to a point near $`y=0`$ is a half circle about $`y=1`$, then the
segment from 1 to 0, then a path from 0 to $`y`$. Iterated integrals of a
composite path are sums over the ways of cutting the word (Chen), also
with tangential base points. On the half circle only the words
$`\omega_1^m`$ survive, with value $`(\pm i\pi)^m/m!`$ (this is
$`\log x=\log(1-y)\pm i\pi`$). On the segment, $`y=1-u`$ turns the word
into one in $`du/u`$, $`-du/(1-u)`$ on $`[0,1]`$, a regularised multiple
zeta value. So $`F=\sum_ud_uI_u(y)`$, where $`I_u`$ is the iterated
integral of $`u`$ from the tangential base point at $`y=0`$. By the
shuffle relations, $`I_{u'\omega_0^k}`$ is a polynomial in $`\log y`$
whose coefficients are $`I_{u''}`$ with $`u''`$ ending in $`\omega_1`$ or
empty, and those are power series in $`y`$ with radius of convergence 1
that vanish at $`y=0`$ (the only other singular point of the forms is
$`y=1`$). If the coefficients $`d_u`$ with $`u`$ ending in $`\omega_0`$
vanish, no logarithm is present. The derivative is read off from
$`I_{\omega_0^k\omega_1}=-\mathrm{Li}_{k+1}(y)=-y+O(y^2)`$, all other
non-empty words being $`O(y^2)`$. $`\blacksquare`$

**The computation** (`geocusp.hpl_exact`). The $`d_u`$ lie in
$`\mathbb Q[i\pi,\zeta_3,\zeta_5]`$ once the multiple zeta values of
weight up to six are reduced to $`\zeta_2=\pi^2/6`$, $`\zeta_3`$,
$`\zeta_5`$. The module does not assume a table for that. It writes
down, in each weight, the shuffle products, the stuffle products and
Hoffman's relation (the double shuffle relations of Ihara, Kaneko and
Zagier, all theorems), and solves them in rational arithmetic; the
solver checks that the system has full rank, so that every value is
determined, and that it is consistent (1, 2, 4, 8, 16 unknowns in
weights 2 to 6). It reproduces $`\zeta(2,1)=\zeta_3`$,
$`\zeta(3,1)=\pi^4/360`$, $`\zeta(3,2)=3\zeta_2\zeta_3-\tfrac{11}2\zeta_5`$
and the others tested. With that table the coefficients of the published
expressions, which are rational multiples of even powers of $`\pi`$,
$`\zeta_3`$ and $`\zeta_5`$, give:

| Function | Conditions $`d_u=0`$ | Result | Value at $`x=-1`$ |
| --- | --- | --- | --- |
| $`F^{(2)}`$ | 3 | all hold | $`11\pi^4/2880`$ |
| $`F^{(3)}_{\rm P}`$ | 15 | all hold | $`\zeta_3^2/32-67\pi^6/51840=-1.19738`$ |
| $`F^{(3)}_{\rm NP1}`$ | 4 | all hold | $`0`$ |
| $`F^{(3)}_{\rm NP2}`$ | 7 | all hold | $`21\zeta_3^2/8+49\pi^6/8640=9.24529`$ |

on both lips, as identities in $`\mathbb Q[i\pi,\zeta_3,\zeta_5]`$: every
one of the coefficients is the zero polynomial. The two-loop row is a
second proof of GATE_6_BOUNDS section 2 ($`64\times11\pi^4/2880=11\pi^4/45`$),
by a different method. The zero of $`F^{(3)}_{\rm NP1}`$ cancels the pole
of $`1/(1+x)`$, and the non-planar part has the value

```math
\Big[\frac{F^{(3)}_{\rm NP1}}{1+x}+F^{(3)}_{\rm NP2}\Big]_{x=-1}
=\frac{21}{8}\zeta_3^2+\frac{3\pi^2}{8}\zeta_3-\frac{15}{16}\zeta_5+\frac{49\pi^6}{8640}=12.7221\ldots
```

the same on both lips. The finiteness is not an accident: the paper
imposes it at $`x=-1`$ as the boundary condition that fixes the
integration constants of its master integrals (its discussion of the
differential equations: "impose the finiteness condition at $`x=-1`$").
What is shown here is more than finiteness: there is no logarithm at
all, so no branch point. As a check of the whole
expansion and not only of the vanishing coefficients, the series
$`\sum_ud_uI_u(y)`$ summed at $`|y|=0.3`$ agrees with the independent
numerical continuation of `geocusp.hpl` to about $`2\times10^{-13}`$ for
all four functions on both lips.

**Corollary (the strip).** For the planar part and for the non-planar
part of the three-loop weight, $`f(L)=F(e^L)`$ is holomorphic on
$`|\operatorname{Im}L|<3\pi`$, and for every $`\delta>0`$
$`|f(L)|\le C_\delta(1+|L|)^6`$ on $`|\operatorname{Im}L|\le3\pi-\delta`$.

*Proof.* *Holomorphy.* $`F`$ is holomorphic on the plane cut along
$`x\le0`$, that is on $`|\operatorname{Im}L|<\pi`$, and continues along
every path that avoids $`x\in\{0,-1,\infty\}`$. In the variable $`L`$ the
only such points in the strip are $`L=\pm i\pi`$. By the proposition the
germ that arrives at $`\pm i\pi`$ from the principal sheet is holomorphic
there, so its continuation round either point returns to itself; every
path in the strip from the principal sheet therefore leads to the same
germ, and $`f`$ is single valued and holomorphic on the strip.

*Growth.* For $`\operatorname{Re}L\le-\log2`$: each $`H_w(x)`$ is a
polynomial in $`\log x`$ of degree at most its weight whose coefficients
are power series in $`x`$ of radius 1 (Remiddi and Vermaseren). That
expression is entire in $`L`$ on this half strip, which is simply
connected and meets the principal sheet, so it is the continuation, and
it is bounded by $`C(1+|L|)^6`$ there; $`1/(1+x)`$ is bounded by 2. For
$`\operatorname{Re}L\ge\log2`$: with $`t=1/x`$, $`dx/x=-dt/t`$ and
$`dx/(1+x)=dt/(1+t)-dt/t`$, so by the composition of paths (from
$`t=\infty`$ to $`0`$ along the positive axis, then to $`t`$) $`F`$ is a
combination with constant coefficients of harmonic polylogarithms of
$`t`$ with the same letters and weight at most six, and the same
argument applies with $`\log t=-L`$; $`1/(1+x)`$ is bounded by 1. The
rest of $`|\operatorname{Im}L|\le3\pi-\delta`$ is compact, and $`f`$ is
continuous on it. $`\blacksquare`$

So the hypotheses of Lemma 3 of GATE_6_BOUNDS hold for the three-loop
weight with $`W=3\pi-\delta`$ and $`k=6`$, exactly as for the two-loop
weight. The theorems of this repository that are stated for a weight
satisfying that hypothesis therefore apply to it, for the planar and
the non-planar part separately: the double pole $`4M[\bar F^{(3)}]`$ of the
regulated integral, in position space for spacelike and for timelike
diagonals (CUSP_LOCAL_MODEL sections 12 to 14) and by the contour route
where the certificate passes.

An independent review re-derived the proposition, read the module line
by line, and checked the result with its own code, which does not use
this repository's evaluator: all 126 regularised multiple zeta values of
weight up to six against 40-digit iterated integrals ($`2\times10^{-40}`$);
the four functions continued to a circle round $`x=-1`$ on both lips, with
monodromy below $`10^{-38}`$, values and slopes equal to the exact ones to
$`10^{-22}`$, and the module's series equal to its continuation to
$`10^{-39}`$. It found one wrong sign in the table above, corrected, and
asked for the qualification "computer-assisted". It did not check that
the files are the published functions. On that point: the files used
are the ancillary files of the arXiv submission itself, downloaded
twice with identical contents, and the circle means computed from them
reproduce the published cusp anomalous dimension through four loops
(table above), which a wrong coefficient would not. The theorems
applied to the weight are stated for the class of Lemma 3 of
GATE_6_BOUNDS (steps A to D and Theorem M there with $`W>\Phi`$;
CUSP_LOCAL_MODEL section 12 with $`W>2\pi`$, section 13 with
$`W>2\pi+4e`$, section 14 with $`W>5\pi/2`$), and use nothing else about
the weight; $`W=3\pi-\delta`$ satisfies each.

What remains numerical: that the strip cannot be made wider. After one
turn around the origin the loop around $`x=-1`$ changes both parts by an
amount of order $`10^3`$ to $`10^4`$ (`monodromy_at_minus_one`), so there are
singularities at $`L=\pm3\pi i`$, as for the two-loop weight. Nothing
depends on this. And the published expressions themselves are taken from
the paper's files.

The one place where "which region contributes" was a real question is
the factor $`1/(1+x)`$: a pole at $`z=-1`$ would have put a new
contribution on the edges themselves, where the theorem's superposition
does not apply. The zero of $`F^{(3)}_{\rm NP1}`$, now exact, removes it.

## What this does not show

- It is not a case where the coefficient was unknown. The method reduces
  the coefficient to a circle mean of the weight; the difficulty of a new
  order is in obtaining the weight, which this method does not do. The
  four-loop weight, which would give the five-loop coefficient, is not
  available to me.
- It is not a different integral. The contour, the regions and the
  theorem are the same; only the weight changed. Nothing here tests the
  method on more than four cusps, on a weight that depends on more than
  one cross-ratio, or on a contour that is not light-like.
