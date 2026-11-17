# Gate 6: an explicit integration contour

[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md) left one step open: a global integration
contour on which the regulated integrand is absolutely integrable. It was
needed twice, for the exchange of the Mellin superposition with the loop
integral and for the statement that only the edges produce poles. This
document constructs such a contour, and replaces the appeal to Stirling's
formula in Lemma 4 by elementary bounds.

**Status.** The construction is explicit and implemented. Its central
inequality is proved. Its global properties follow from that inequality, a
compactness argument and a choice of constants; they are verified on a
stress sample of 60 000 points, not proved estimate by estimate. It is
carried out for spacelike diagonals (\(s,t<0\)). Nothing here has been
reviewed by a specialist.

Code: `src/geocusp/feynman_contour.py`. Tests: `tests/test_gate6_contour.py`.
Example 63: `examples/feynman_contour_survey.py`.

## 1. The form of the displacement

Displace the real domain by a purely imaginary vector field:

\[
x_0\ \longmapsto\ x_0+iV(x_0),\qquad V\ \text{real}.
\]

Then, exactly,

\[
x_{0i}^2\ \longmapsto\ x_{0i}^2-V^2+2i\,V\!\cdot(x_0-x_i).
\]

The imaginary part is linear in \(V\). Call \(V\) *admissible for \(i\) at
\(x_0\)* if \(V\!\cdot(x_0-x_i)\le0\). Admissible fields form a convex cone,
so fields that are admissible on overlapping regions can be joined by a
partition of unity. This is what makes a global construction possible.

The Feynman prescription is the limit of small admissible \(V\) near each
cone: the invariant moves into the lower half plane.

## 2. Frame fields

At the cusp \(x_j\) write \(x_0=x_j+\alpha p_{j-1}+\beta p_j+y_\perp\),
\(D=2p_{j-1}\!\cdot p_j<0\), \(\rho=|y_\perp|^2\), \(P=(\beta,\alpha)\). Let

\[
T=\{\beta\ge0,\ \alpha\le0,\ \beta-\alpha\le1\}
\]

be the triangle with vertices \((0,0)\), \((0,-1)\), \((1,0)\), which are
the images of \(x_j\), \(x_{j-1}\), \(x_{j+1}\). Let
\(d(P)=P-\pi_T(P)\), where \(\pi_T\) is the nearest-point projection onto
\(T\). Define

\[
V_j=\kappa\,(d_\beta\,p_{j-1}+d_\alpha\,p_j)+\nu\,y_\perp,\qquad\kappa,\nu>0 .
\]

**Lemma 5.** For \(i\in\{j-1,j,j+1\}\) and every real \(x_0\),

\[
V_j\!\cdot(x_0-x_i)\ \le\ -\frac{\kappa|D|}{2}\,|d(P)|^2-\nu\rho\ \le\ 0 .
\]

*Proof.* A direct computation with \(p_{j-1}^2=p_j^2=0\) and
\(y_\perp\!\cdot p_{j-1}=y_\perp\!\cdot p_j=0\) gives

\[
V_j\!\cdot(x_0-x_i)=\frac{\kappa D}{2}\;d\cdot(P-v_i)-\nu\rho ,
\]

where \(v_i\) is the vertex of \(T\) belonging to \(x_i\). For the
projection onto a closed convex set, \((P-\pi_T(P))\cdot(v-\pi_T(P))\le0\)
for every \(v\in T\). Hence
\(d\cdot(P-v)=|d|^2-d\cdot(v-\pi_T(P))\ge|d|^2\), and \(D<0\).
\(\blacksquare\)

So \(V_j\) is admissible for its three invariants everywhere, and strictly
so unless \(\rho=0\) and \(P\in T\), that is unless \(x_0\) lies in the flat
triangle \(\operatorname{conv}(x_{j-1},x_j,x_{j+1})\). On that triangle the
three invariants are \(\alpha\beta D\), \(\beta D(1+\alpha)\) and
\(\alpha D(\beta-1)\), which vanish only on the two edges at \(x_j\). The
triangle is exactly the set on which the Landau conditions fail to give a
displacement direction for all three invariants at once.

Near the edge through \(x_j\) and \(x_{j+1}\) this field reduces to the
rotation of the transverse plane used in GATE_6_BOUNDS, with
\(\rho\to(1+i\nu)^2\rho\).

The frame field does not control the fourth invariant \(x_{0,j+2}^2\).

## 3. Far fields

Let \(c_-\) lie in the common chronological past of the four cusps and
\(c_+\) in their common future. For \(V=-\kappa'(x_0-c_-)\),

\[
V\!\cdot(x_0-x_i)=-\kappa'\big[x_{0i}^2+(x_i-c_-)\!\cdot(x_0-x_i)\big].
\]

On the future nappe of the cone of \(x_i\) the first term vanishes and the
second is the product of a future timelike and a future null vector, which
is positive. So this field is admissible near every future nappe at once,
and \(-\kappa'(x_0-c_+)\) near every past nappe.

A point can be near a future nappe of one cone and a past nappe of another
only in a bounded region: \(x_i-x_j\) would have to be nearly a difference
of a past and a future null vector, hence nearly causal, and for adjacent
cusps this forces \(x_0\) to the edge between them.

## 4. The global field

Fix a ball of radius \(R\) about the centre of the square and a width
\(\delta_0\).

- Inside the ball: \(V_{\rm in}=\sum_jm_jV_j\) with
  \(m_j=\chi(|x_{0,j+2}^2|/\delta_0)\), where \(\chi\) vanishes on \([0,1]\)
  and equals 1 on \([2,\infty)\). The frame field at \(x_j\) is switched off
  near the cone of the opposite cusp, the one invariant it does not control.
- Outside \(2R\): \(V_{\rm out}=-\kappa'\Theta_+(x_0-c_-)-\kappa'\Theta_-(x_0-c_+)\),
  with \(\Theta_\pm\) supported where the time coordinate is at least a
  quarter of the spatial distance, in the future or the past.
- In between the two are blended.

**Why something is always switched on.** Near the cone of \(x_i\) the usable
frames are \(j=i-1,i,i+1\), with weights
\(\chi(|x_{0,i+1}^2|/\delta_0)\), \(\chi(|x_{0,i+2}^2|/\delta_0)\),
\(\chi(|x_{0,i-1}^2|/\delta_0)\). All three vanish only if \(x_0\) is within
\(2\delta_0\) of all four cones. That does not happen:

\[
\min_{x_0\in\mathbb R^4}\ \max_i|x_{0i}^2|\approx0.199
\]

for the square used here (numerical minimisation), and \(2\delta_0=0.1\).
This is the statement that the four cones have no common real point.

**Constants.** \(R\) is chosen so that the far fields are admissible outside
the ball. \(\delta_0\) is then chosen below the minimum above. Finally
\(\kappa,\nu\ll\delta_0/R^2\), so that where a frame field is switched on but
does not control an invariant, that invariant is at least \(\delta_0\) in
modulus and moves by a small fraction of itself. The values used are
\(R=8\), \(\delta_0=0.05\), \(\kappa=\nu=3\times10^{-6}\),
\(\kappa'=2\times10^{-3}\).

## 5. Properties

**Proposition 6.** For the field of section 4:

1. near the cone of \(x_i\), \(\operatorname{Im}x_{0i}^2\le0\) on the
   contour, with equality only on the edges;
2. away from the cones each invariant stays within a small angle of its
   real value, so \(|\arg z|\le2\pi+\delta''\) on the contour with
   \(\delta''\) as small as desired;
3. the invariants vanish on the contour only on the four edge segments;
4. near a cusp, in its frame,
   \(|x_{0,j}^2|\ge c\,\nu\max(\rho,|\alpha\beta D|)\),
   \(|x_{0,j-1}^2|\ge c\,\nu\max(\rho,|\beta D||1+\alpha|)\),
   \(|x_{0,j+1}^2|\ge c\,\nu\max(\rho,|\alpha D||1-\beta|)\);
5. the map \(x_0\mapsto x_0+iV\) is a local diffeomorphism onto the contour.

*Argument.* (1) and (3): by Lemma 5 and section 3 every field that is
switched on near a cone is admissible for it; by section 4 at least one
frame field has weight 1 there, and its bound vanishes only on its triangle,
which meets the cone only on the edges. (2): by the choice of constants.
(4): if \(\rho\ge|\alpha\beta D|/2\) the imaginary part is at least
\(2\nu\rho\); otherwise the real part is at least \(|\alpha\beta D|/2\); the
other two are the same. (5): the Jacobian is \(1+i\,\partial V\) with
\(\partial V\) real, whose determinant vanishes only if \(\partial V\) has
an eigenvalue of modulus 1; \(\|\partial V\|<1\).

Measured on the sample (60 000 points: generic points at three scales,
points within \(10^{-4}\) to \(10^{-1}\) of one cone, points near two cones,
and points near the four null lines on and off the segments):

| Property | Measured |
| --- | --- |
| Largest \(\operatorname{Im}x_{0i}^2\) over 42 716 near-cone pairs | \(-5.5\times10^{-16}\) (never positive) |
| Largest excursion of \(\arg x_{0i}^2\) beyond \([-\pi,0]\) | \(3.9\times10^{-3}\) rad |
| Smallest \(\min_i|x_{0i}^2|/\min(1,\mathrm{dist}_{\rm edges}^2)\) | \(1.8\times10^{-6}\) |
| Smallest ratio in item 4, in units of \(\nu\), at the four cusps | 213 to 1152 |
| Largest \(\|\partial V\|\) over 8 000 points | 0.19 |

## 6. Consequences

**Absolute integrability.** By item 4 the regulated integrand near a cusp
is bounded by a constant times

\[
\frac{\rho^{\eta}\,d\alpha\,d\beta\,d\rho}
{\max(\rho,|\alpha\beta|)\,\max(\rho,|\beta|)\,\max(\rho,|\alpha|)},
\qquad\eta=-\epsilon>0,
\]

times a polynomial in logarithms. For \(\rho<|\alpha\beta|\) the \(\rho\)
integral gives \(|\alpha\beta|^{\eta-1}\); for
\(|\alpha\beta|<\rho<\min(|\alpha|,|\beta|)\) it gives
\(\min(|\alpha|,|\beta|)^\eta/(\eta|\alpha\beta|)\); both are integrable in
\((\alpha,\beta)\) for \(\eta>0\), and \(\rho>\min(|\alpha|,|\beta|)\) is
easier. Along the interior of an edge the same bounds apply with one of
\(\alpha,\beta\) of order one. Away from the edges the integrand is bounded
by item 3, and at infinity it decays as \(|x_0|^{-8}\) with the cones
resolved by the far field. So the regulated integrand is absolutely
integrable on the contour for small \(\eta>0\).

**The Fubini step of GATE_6_BOUNDS section 5.** On the contour
\(|z^p|\le|z|^{\operatorname{Re}p}e^{(2\pi+\delta'')|\operatorname{Im}p|}\)
by item 2, the integrand is absolutely integrable, and
\(|c(p)|\le Ce^{-(2\pi+\delta')|\operatorname{Im}p|}\) by Lemma 3. Taking
\(\delta''<\delta'\), the double integral over the contour and over \(p\)
converges absolutely, and the order can be exchanged.

**The displacement in Proposition 5.** Item 3 is that displacement: off the
edges the contour does not meet any light cone, so the integrand is bounded
there uniformly in \(\epsilon\) near 0, and that part of the integral is
analytic in \(\epsilon\). Poles come from neighbourhoods of the edges only.

**Equality with the Feynman integral.** The family \(x_0+i\lambda V\),
\(0<\lambda\le1\), satisfies items 1 and 3 for every \(\lambda\), so no
singularity crosses the contour as \(\lambda\) varies and the integral is
independent of \(\lambda\); its limit \(\lambda\to0\) is the \(-i0\)
prescription. This is Cauchy's theorem in several variables; the vanishing
of the contribution at infinity uses the decay \(|x_0|^{-8}\).

## 7. Lemma 4 without Stirling

Lemma 4 of GATE_6_BOUNDS used Stirling's bound. The following suffices and
is elementary. From the Weierstrass product,

\[
\frac{|\Gamma(a+iv)|^2}{\Gamma(a)^2}=\prod_{n\ge0}\Big(1+\frac{v^2}{(a+n)^2}\Big)^{-1},
\qquad
\prod_{n\ge1}\Big(1+\frac{v^2}{n^2}\Big)=\frac{\sinh\pi v}{\pi v}.
\]

Comparing factor by factor: for \(0<a\le1\), \(a+n\le n+1\), so

\[
|\Gamma(a+iv)|^2\le\Gamma(a)^2\,\frac{\pi v}{\sinh\pi v};
\qquad
\text{for }1\le a\le2:\quad
|\Gamma(a+iv)|^2\le\Gamma(a)^2\,(1+v^2)\,\frac{\pi v}{\sinh\pi v}.
\]

Also \(\Gamma(1+p)\Gamma(1-p)=\pi p/\sin\pi p\) exactly, and
\(|\Gamma(\tfrac\eta2+iy)|^2=|\Gamma(1+\tfrac\eta2+iy)|^2/(\eta^2/4+y^2)\).

In the integrand of \(N_\eta(p)\) the two \(p\)-dependent Gamma functions
have real parts in a compact subset of \((0,2]\), so their product is at
most \(C(1+v^2)^2\,\pi|v|/\sinh\pi|v|\le C'(1+|v|)^5e^{-\pi|v|}\) with
\(v=\tau-y\). The weight is at most
\(C(1+y^2)^2(\pi y/\sinh\pi y)^2/(\eta^2/4+y^2)^2\), whose integral against
\(e^{\pi|y|}(1+|y|)^5\) is at most \(C/\eta^3\). Finally
\(\eta^2/\Gamma(2\eta)=2\eta^3/\Gamma(1+2\eta)\le2.3\,\eta^3\). Therefore

\[
|N_\eta(\sigma+i\tau)|\le C\,(1+|\tau|)^3\,e^{\pi|\tau|},
\]

uniformly in \(\eta\). The power of \(|\tau|\) is worse than in Lemma 4 and
does not matter: \(c(p)\) decays like \(e^{-(2\pi+\delta')|\tau|}\). The
test `test_elementary_gamma_bounds` checks the three inequalities on a grid.

## 8. What remains

- The global properties in Proposition 6 are derived from Lemma 5 by a
  compactness argument and a choice of constants, and checked by sampling.
  A proof that tracks every constant, including the blend region between the
  ball and the far field, is not written.
- The construction is for \(s,t<0\). In the other two kinematic regions one
  diagonal is timelike, the signs in Lemma 5 change at two of the cusps, and
  the far field must be adapted. The local results of GATE_6_BOUNDS hold in
  all three regions.
- Dimensional regularisation is applied on the contour by continuing the
  transverse measure. That the result agrees with the Mellin–Barnes
  definition of the regulated integral is standard and not re-derived.
- Review by someone who works with these integrals.

## 9. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_contour.py
PYTHONPATH=src python -m examples.feynman_contour_survey
```
