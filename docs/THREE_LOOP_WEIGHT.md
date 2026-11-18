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

Numerically:

- **At $`z=-1`$ on the principal sheet.** Continued once around
  $`x=-1`$, both parts return to themselves, to $`10^{-12}`$: there is no
  branch point. The planar part tends to $`-1.19738`$. The factor
  $`F^{(3)}_{\rm NP1}`$ vanishes linearly ($`3\times10^{-6}`$ at distance
  $`10^{-6}`$), so the pole of $`1/(1+x)`$ is cancelled and the non-planar
  part tends to $`12.722`$.
- **On the next sheet.** After one turn around the origin the same loop
  around $`x=-1`$ changes both parts by an amount of order $`10^3`$. So
  the nearest singularities are at $`L=\pm3\pi i`$, as for the two-loop
  weight, and the strip is $`W=3\pi-\delta`$.
- **Growth.** $`|F|/(1+|L|)^6`$ is below $`10^{-3}`$ for
  $`|L|`$ up to 40 on the real axis.

So the hypotheses hold as far as these checks go, and the theorem then
gives the double pole $`4M[\bar F^{(3)}]`$ of the regulated integral, in
position space for spacelike diagonals and by the contour route in all
three regions, for the planar and the non-planar part separately. This
is numerical evidence for the hypotheses, not a proof of them: for the
two-loop weight holomorphy on the strip is proved (GATE_6_BOUNDS section
2), and the same argument has not been written for weight six.

The one place where "which region contributes" was a real question is
the factor $`1/(1+x)`$: a pole at $`z=-1`$ would have put a new
contribution on the edges themselves, where the theorem's superposition
does not apply. The zero of $`F^{(3)}_{\rm NP1}`$ removes it.

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
