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
reproduced it in a second chart. At the second level, the further residues
of the facet $`w_1=0`$ at each of the other seven coordinates and at
$`Q=0`$ equal, exactly, the forms of those eight faces computed from their
inequalities (`two_loop_second_level_check`); at $`Q=0`$ the piece with
$`A>0`$ enters with a minus sign, by the same orientation rule. The
symmetries above carry $`w_1=0`$ to every coordinate facet, and residues
at two facets can be taken in either order up to sign, so these eight
represent every face of codimension two among the nine finite facets
(the 36 pairs of facets fall into eight orbits, each containing one of
these): the second level holds on all of them. The form also has poles
at infinity in this chart: at $`w_1=\infty`$, at $`w_2=\infty`$ on the
facet $`w_1=0`$, and at $`x_1=x_2`$ on the facet $`Q=0`$, where $`w_1`$
runs to infinity. At each of these three the residue equals the form of
the face cut out by the inequalities there, up to the orientation of the
upper end of an interval (`two_loop_infinity_check`): for example at
$`w_1=\infty`$ the condition $`Q>0`$ forces $`x_1>x_2`$ and leaves the other
variables free. The remaining pairs involving a face at infinity are not
checked. The face $`w_1=w_2=0`$ is the union of two
products of orthants and ordered pairs that meet in codimension two; its
form is the sum of two known forms, given additivity.

**Every depth.** `region_form` computes the form of a region
$`\{v_i>0,\ c_k>0\}`$ from its inequalities, for the regions that occur
here: a constraint that factors is split into the sign patterns of its
factors, and a variable occurring in one constraint only,
$`\alpha v+\beta>0`$, is integrated out by the interval rule ($`1/v`$,
$`\alpha/(\alpha v+\beta)`$ or $`\beta/(v(\alpha v+\beta))`$ according to
the signs of $`\alpha`$ and $`\beta`$). It is not a general algorithm: it
declines regions in which every variable occurs in two constraints.
Applied to $`Q>0`$ it reproduces the two-loop integrand.

*Coordinate faces.* `coordinate_recursion_check` starts from the full
region and visits every face reached by setting coordinates to zero in
any order. On each, the residue of its form at every remaining coordinate
is compared with the form of the next face. All agree exactly. The
traversal has 399 states and 1360 comparisons. Of the states, 199 are
proper faces, 100 are whole orthants and 100 are empty; of the
comparisons, 720 are between two non-zero forms, 160 test that a proper
face has no pole where it should have none, and 480 are trivial.

One thing had to be got right for this. On some steps the inequality
vanishes identically on the next face: for instance on $`x_1=x_2=y_1=0`$
the quadric is $`-y_2(z_1-z_2)`$, which is zero at $`y_2=0`$. The face
there is not empty and not the whole orthant. It is the limit of the
region, $`z_1<z_2`$, and the residue is its form,
$`1/(z_1(z_2-z_1)w_1w_2)`$. Reached in the other order, through
$`x_1=x_2=y_2=0`$ and then $`y_1=0`$, the same coordinate plane carries
the other half, $`z_1>z_2`$; reached through $`y_1=y_2`$ it carries a
third object, the signed union. Three divisors pass through that plane,
so the boundary is not a normal crossing there, and a face is defined by
the chain of residues that reaches it. The rule (`face_constraint`) is:
the inequality on the face is the first non-vanishing coefficient of the
constraint in the vanishing variable. With the naive restriction the
recursion fails at every such step that has a non-zero residue (112
steps in a run over subsets of coordinates); with the rule it holds.

*Zeros and infinities.* `full_recursion_check` treats the chart as eight
intervals $`[0,\infty]`$ and sends each coordinate to either end, in every
order. At infinity the inequality on the face is the leading coefficient
of the constraint. This gives 1323 states and 10 080 comparisons, 6096 of
them with a non-zero residue, and all agree exactly.

*The quadric.* `quadric_full_recursion_check` takes every one of those
faces that still carries an inequality $`c>0`$ and every chart in which
$`c=\alpha u+\beta`$ is linear in an eliminated variable $`u`$. The
residue of the face's form at $`c=0`$ equals the form of
$`\{\alpha>0,\beta<0\}`$ minus that of $`\{\alpha<0,\beta>0\}`$, the
orientation rule above (1112 residues). Each such quadric face is then
followed down: at each end of each remaining coordinate $`v`$ the residue
equals the same expression with $`\alpha`$ and $`\beta`$ replaced by their
leading coefficients in $`v`$ at that end, which say where $`u`$ stays
positive in the limit. This one rule covers the steps at which the chart
degenerates ($`u\to0`$ or $`u\to\infty`$ on the deeper face), which an
earlier version skipped. Followed to the bottom there are 4328 quadric
faces in charts and 29 920 comparisons, 22 688 with a non-zero residue,
and all agree exactly.

**What this shows.** In this chart, compactified as $`[0,\infty]^8`$, the
form derived from the inequalities has, at every boundary reached by
coordinates going to zero or infinity and by the quadric, in every order
and to every depth checked, a residue equal to the form derived from the
inequalities of that boundary. A check on coordinate residues alone could
not see a term such as $`k/Q`$; the quadric residues exclude it. There
are no stray poles (`pole_check`): the denominator of each of the 1323
face forms consists of remaining coordinates and factors of that face's
own inequality, each to the first power, and the same holds for the 4328
quadric-face forms with the factors of $`\alpha`$ and $`\beta`$.

It is not a proof that the form is the canonical form of the region.
That needs one of two things. Either each interval piece is a positive
geometry whose form is the product form and forms add over subdivisions,
in which case `region_form` is canonical by construction; or the
boundary divisors of the region in its natural compactification, the
product of two Grassmannians, are all among those checked here, in which
case uniqueness gives it. I have not shown either. Two boundaries of a
quadric face in a chart are not coordinates of that chart, $`u=0`$ and
$`u=\infty`$; they are coordinates of other charts and are reached there.

An independent review of the routine and of the coordinate and first
quadric traversals wrote its own version, which eliminates a different
variable and builds each interval from its end points. It agreed on all
396 non-trivial faces, matched 98 random polytopes and cones against the
dual-volume formula, reproduced all the counts, and supplied the
breakdown into trivial and non-trivial comparisons quoted above. The
traversals with infinities and the following of quadric faces to every
depth were written after that review. The inequalities themselves are
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
