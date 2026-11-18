# What is derived from the geometry

Two statements that use the positive geometry itself, not a formula taken
from the literature. Code: `src/geocusp/negative_geometry.py`,
`src/geocusp/amplituhedron.py`. Tests: `tests/test_from_geometry.py`.

## 1. The two-loop form from its inequalities

In the chart $`AB_i=\begin{pmatrix}1&x_i&0&-y_i\\0&z_i&1&w_i\end{pmatrix}`$
the two-loop negative geometry is the region

```math
x_i,\,y_i,\,z_i,\,w_i>0,\qquad
Q=(x_1-x_2)(w_1-w_2)+(y_1-y_2)(z_1-z_2)>0 ,
```

that is $`\langle AB_1AB_2\rangle<0`$. Its canonical form is computed here
from these inequalities alone.

Put $`A=x_1-x_2`$ and $`c=(y_1-y_2)(z_1-z_2)`$, so $`Q=A(w_1-w_2)+c`$ is
linear in $`w_1`$ and in $`w_2`$. The region is the disjoint union of five
pieces. On each, $`w_1`$ and then $`w_2`$ run over an interval whose ends are
rational in the other variables, over a base that is a product of
simplices:

| Piece | $`w_2`$ | $`w_1`$ |
| --- | --- | --- |
| $`A>0,\ c>0`$ | $`0<w_2<c/A`$ | $`w_1>0`$ |
| $`A>0,\ c>0`$ | $`w_2>c/A`$ | $`w_1>w_2-c/A`$ |
| $`A>0,\ c<0`$ | $`w_2>0`$ | $`w_1>w_2-c/A`$ |
| $`A<0,\ c>0`$ | $`w_2>0`$ | $`0<w_1<w_2-c/A`$ |
| $`A<0,\ c<0`$ | $`w_2>c/A`$ | $`0<w_1<w_2-c/A`$ |

The form of an interval $`(a,\infty)`$ is $`dw/(w-a)`$ and of $`(0,b)`$ is
$`b\,dw/(w(b-w))`$. The base forms are those of the orderings, for example
$`dx_1dx_2/(x_2(x_1-x_2))`$ for $`x_1>x_2>0`$; $`c>0`$ is the union of
$`\{y_1>y_2,z_1>z_2\}`$ and $`\{y_1<y_2,z_1<z_2\}`$. Each piece is a
generalised simplex, so its form is the product of these.

**Result.** The sum of the five forms is, exactly,

```math
\frac{w_1x_1+w_2x_2+y_1z_1+y_2z_2}{D_1D_2\,Q}
=-\frac{n_{12}}{D_1D_2X_{12}},\qquad D_i=w_ix_iy_iz_i,\quad X_{12}=-Q .
```

All spurious poles, at $`x_1=x_2`$, $`y_1=y_2`$, $`z_1=z_2`$,
$`Aw_2=c`$, cancel. This is the two-loop logarithm integrand used
throughout the repository; it was previously entered as a formula.

Checks: the identity is verified symbolically; on 5000 random points every
point with $`Q>0`$ lies in exactly one piece and no point with $`Q<0`$ lies
in any; each piece's form is positive on its piece.

**Assumptions.** That the sum is the canonical form assumes that the region
is a positive geometry and that forms add under triangulation, which is
standard. Checked here are the tiling, the cancellation of spurious poles,
positivity, and one level of the residue recursion: the residue of the form
at the facet $`w_1=0`$ equals, exactly, the canonical form of that facet
computed independently by the same interval rule in $`w_2`$
(`two_loop_boundary_check`). The form, $`Q`$ and the region are invariant
under the exchange of the two loops and, applied to both loops at once,
under $`x\leftrightarrow w`$, $`y\leftrightarrow z`$ and
$`(x,w)\leftrightarrow(y,z)`$. These are even permutations, so they
preserve orientation, and together they reach all eight coordinate facets,
where the same therefore holds. At the ninth facet, $`Q=0`$, eliminate
$`w_1=w_2-c/A`$: the residue equals, exactly, the form of the facet
computed by the interval rule in $`w_2`$, with the pieces having $`A<0`$
taken with the opposite sign because $`dw_1=dQ/A`$ reverses the induced
orientation there (`two_loop_mutual_facet_check`). So the first level of
the recursion holds at every facet. The sign rule is the induced boundary
orientation, not an adjustment: an independent review found the residue
negative at every one of 2385 sampled facet points with $`A<0`$, and
reproduced it in a second chart. At the second level, the further residues of
the facet $`w_1=0`$ at $`w_2=0`$, $`x_1=0`$ and $`x_2=0`$ equal, exactly, the
forms of those three faces computed from their inequalities
(`two_loop_second_level_check`). The remaining faces of codimension two
and all deeper levels are not checked. The inequalities themselves are
taken from the literature. An independent review rebuilt the five forms and
confirmed the identity, the tiling on 20 000 points, and that
$`\langle AB_1AB_2\rangle=-Q`$.

**What this is.** A re-derivation of a published integrand
(arXiv:2112.06956) by the standard method of cutting the region into
simplices. It is not a new result. It does make the chain

```text
inequalities  ->  two-loop form  ->  weight log^2 z + pi^2  ->  contour  ->  -4 pi^2/3
```

computed inside the repository at this order: the second arrow is the
quadrature of `examples/cd_integral.py`, the third the contour rule.

**What it is not.** The same has not been done at three loops. The chain
graph would go the same way with more pieces. The triangle graph, with
three mutual conditions, does not reduce to intervals in this simple way,
and the two remaining loop integrals that give $`\bar F^{(2)}`$ are not
attempted. So $`\bar F^{(2)}`$ remains an input.

## 2. An independent membership test for the sampler

The sampler builds $`Y=CZ`$ with $`C`$ in the positive Grassmannian, so
positivity of $`C`$ is a construction, not a test. The sign-flip
characterisation of arXiv:1704.05069, eq. (5.13), tests $`Y`$ without any
$`C`$. The paper states one condition $`[Y\,i\,i{+}1\,j\,j{+}1]>0`$ with the
convention $`Z_{n+1}=(-1)^{k-1}Z_1`$; written out, that is:

```math
\langle Y\,i\,i{+}1\,j\,j{+}1\rangle>0,\qquad
(-1)^{k-1}\langle Y\,i\,i{+}1\,n\,1\rangle>0,\qquad
\{\langle Y\,1\,2\,3\,i\rangle\}_{i=4}^{n}\ \text{has exactly }k\text{ sign flips.}
```

| $`k`$ | $`n`$ | Sampled points passing | Random $`C`$ rejected |
| --- | --- | --- | --- |
| 1 | 6 | 400 of 400 | 291 of 300 |
| 1 | 8 | 400 of 400 | 288 of 300 |
| 2 | 6 | 400 of 400 | 300 of 300 |
| 2 | 8 | 400 of 400 | 300 of 300 |

For $`k=1`$ the amplituhedron is a cyclic polytope and membership can be
decided by linear programming. On 400 mixed points the sign-flip test and
the linear programme agree every time. So the test can fail, and the
sampler passes it. This checks a published characterisation numerically; it
proves nothing new.
