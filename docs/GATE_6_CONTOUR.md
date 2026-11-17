# Gate 6: an explicit integration contour

[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md) left one step open: a global integration
contour on which the regulated integrand is absolutely integrable. It was
needed twice, for the exchange of the Mellin superposition with the loop
integral and for the statement that only the edges produce poles. This
document constructs such a contour, and replaces the appeal to Stirling's
formula in Lemma 4 by elementary bounds.

**Status.** The construction is explicit and implemented, in all three real
kinematic regions and for every light-like square: the constants are derived
from the geometry. Its central inequality is proved (Lemma 5). Its global
properties are reduced in section 9 to inequalities between named constants,
each proved, and a certificate evaluates them; they hold across the family of
squares. Absolute integrability is proved line by line in section 10. The
parametric representation of the regulated integral is derived in section
11, including the contour-independence and limiting step (Lemmas 10 to 12).
Section 11 also corrects an earlier statement: the contour must displace the
extra dimensions of dimensional regularisation, not treat them as a real
shift. Nothing here has been reviewed by a specialist.

Code: `src/geocusp/feynman_contour.py`, `src/geocusp/schwinger.py`. Tests:
`tests/test_gate6_contour.py`, `tests/test_schwinger.py`.
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
\(D=2p_{j-1}\!\cdot p_j\), \(\rho=|y_\perp|^2\), \(P=(\beta,\alpha)\). Take
\(D<0\) first (a spacelike diagonal at this cusp). Let

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

**Timelike diagonal.** If \(D>0\) at the cusp, replace \(\kappa\) by
\(-\kappa\) in \(V_j\). The computation in the proof gives
\((-\kappa D/2)\,d\cdot(P-v_i)-\nu\rho\), and the bound of Lemma 5 holds with
\(|D|\). The transverse plane is spacelike for either sign of \(D\), since
two independent null vectors span a Lorentzian plane. Nothing else changes:
the far fields need only points in the common past and future, which always
exist. In the regions \((s,t)=(+,-)\) and \((-,+)\) two of the four cusps
have \(D>0\).

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

for the square used here (numerical minimisation), and \(2\delta_0=0.05\).
Section 9 proves a weaker bound, 0.0897, which is still enough.
This is the statement that the four cones have no common real point.

**Constants.** \(R\) is chosen so that the far fields are admissible outside
the ball. \(\delta_0\) is then chosen so that \(2\delta_0\) is below the
minimum above. Finally \(\kappa,\nu\ll\delta_0/R^2\), so that where a frame
field is switched on but does not control an invariant, that invariant is at
least \(\delta_0\) in modulus and moves by a small fraction of itself. The
values used are \(R=8\), \(\tau_0=8\) (the time offset of \(c_\pm\)),
\(\delta_0=0.025\), \(\kappa=\nu=1.5\times10^{-7}\),
\(\kappa'=2\times10^{-3}\). Section 9 states the precise conditions.

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

Measured on stress samples of 30 000 points in each region (generic points
at three scales, points within \(10^{-4}\) to \(10^{-1}\) of one cone, points
near two cones, and points near the four null lines on and off the
segments), next to the bounds proved in section 9:

| Property | \((-,-)\) | \((+,-)\) | \((-,+)\) | Proved bound |
| --- | --- | --- | --- | --- |
| Near-cone pairs tested | 19 520 | 19 369 | 19 441 | |
| Largest \(\operatorname{Im}x_{0i}^2\) near a cone | \(-7\times10^{-17}\) | \(-2\times10^{-17}\) | \(-9\times10^{-17}\) | \(\le0\) |
| Largest excursion of \(\arg x_{0i}^2\) beyond \([-\pi,0]\) | 0.0039 | 0.0039 | 0.0039 | 0.083, 0.049, 0.049 |
| Smallest \(\min_i|x_{0i}^2|/\min(1,\mathrm{dist}_{\rm edges}^2)\) | \(4\times10^{-7}\) | \(5\times10^{-7}\) | \(5\times10^{-7}\) | \(>0\) |
| Smallest ratio in item 4, units of \(\nu\) | 4 258 | 6 920 | 2 364 | \(\ge\tfrac12\) |
| Largest \(\|\partial V\|\) | 0.018 | 0.018 | 0.018 | 0.37, 0.18, 0.18 |

The imaginary parts near a cone are never positive; their size is set by
\(\nu\) times the squared distance from the flat triangles, which is tiny for
points placed within \(10^{-4}\) of a cone.

## 6. Consequences

**Absolute integrability.** The regulated integrand is absolutely integrable
on the contour for small \(\eta=-\epsilon>0\). The proof is section 10.

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

## 8. Kinematic regions

Real null edges summing to zero have exactly one positive invariant among
\(s,t,u\), so there are three regions. The construction is the same in all
of them, with the sign change of section 2 at the cusps whose diagonal is
timelike. In those regions two opposite cusps are timelike separated, and
the future cone of one meets the past cone of the other in a sphere; that
sphere lies inside the ball, where the two frame fields of the other two
cusps control both invariants. The table of section 5 and the certificate of
section 9 cover all three regions.

## 9. Proof of Proposition 6 with explicit constants

Norms \(|\cdot|\) are Euclidean in \(\mathbb R^4\), and
\(|a\!\cdot b|\le|a||b|\) for the Minkowski product. Let \(c\) be the centre
of the square, \(r=|x_0-c|\), and

\[
\Delta=\max_i|x_i-c|,\qquad P=\max_j|p_j|,\qquad
A_1=\frac{2P}{\min_j|D_j|},\qquad A_2=1+2A_1P,\qquad
\lambda_0=2\kappa A_1P+\nu A_2 .
\]

In the frame of \(x_j\): \(|\alpha|,|\beta|\le A_1|x_0-x_j|\),
\(|y_\perp|\le A_2|x_0-x_j|\), and \(|d(P)|\le|P|\) because the vertex
\((0,0)\) belongs to \(T\). Hence

\[
|V_j|\le\lambda_0|x_0-x_j|,\qquad\operatorname{Lip}(V_j)\le\lambda_0,
\]

the second because \(P\mapsto d(P)\) is 1-Lipschitz. Write
\(n=x_0-x_i\), \(u^\mp=x_i-c_\mp\), and \((\tau,\xi)\) for the time and space
parts of \(x_0-c\). The cutoff \(\chi\) is the cubic smoothstep, with
\(|\chi'|\le\tfrac32\). The inner field is used for \(r\le R\), the outer for
\(r\ge2R\), and \(V=(1-g)V_{\rm in}+gV_{\rm out}\) in between with
\(|\nabla g|\le\tfrac3{2R}\). \(\Theta_\pm\) equals 1 for
\(\pm\tau\ge\tfrac12|\xi|\) and 0 for \(\pm\tau\le\tfrac14|\xi|\).

**(H1) A frame field is fully on near every cone.** If
\(\max_i|x_{0i}^2|\le\epsilon\) then, in the frame of any cusp,
\(|\beta D|,|\alpha D|\le2\epsilon\), \(|w|\le\epsilon\),
\(\rho\le\epsilon+4\epsilon^2/|D|\), and by the exact formula for the far
invariant

\[
\big|x_{0,j+2}^2-D'\big|\le\frac{4\epsilon|D'|}{|D|}+\epsilon
+2\sqrt{\frac{D'u}{D}\Big(\epsilon+\frac{4\epsilon^2}{|D|}\Big)} .
\]

If the right side is less than \(|D'|-\epsilon\) this contradicts
\(|x_{0,j+2}^2|\le\epsilon\). Let \(\epsilon_*\) be the largest \(\epsilon\)
for which it is. *Condition:* \(2\delta_0<\epsilon_*\). Then for
\(|x_{0i}^2|\le\delta_0\) one of the three weights
\(\chi(|x_{0,i\pm1}^2|/\delta_0)\), \(\chi(|x_{0,i+2}^2|/\delta_0)\) equals 1.
(\(\epsilon_*\) is weaker than the true minimum 0.199 but is proved.)

**(A) Inner field near a cone.** For \(|x_{0i}^2|\le\delta_0\) the weight of
the frame \(i+2\) vanishes, and by Lemma 5

\[
V_{\rm in}\!\cdot n=\sum_{j\ne i+2}m_jV_j\!\cdot n
\le-\sum_{j\ne i+2}m_j\Big(\frac{\kappa|D_j|}{2}|d_j|^2+\nu\rho_j\Big)\le0 .
\]

This is exact and needs no condition. By (H1) some \(m_j=1\), so equality
forces \(x_0\) into the flat triangle of that frame, where its three
invariants vanish only on the edges.

**(C0) Centres.** *Condition:* \(\tau_0\ge3\Delta\). Then \(u^-\) is future
timelike with \(u^-_0\ge\tau_0-\Delta\) and \(|\vec u^-|\le\Delta\).

**(C2) Near a cone outside the ball, \(\Theta=1\).** If
\(|n^2|\le\delta_0\) then \(|n_0|\ge|\vec n|-\sqrt{\delta_0}\), so
\(|\tau|\ge|\xi|-(2\Delta+\sqrt{\delta_0})\). With
\(\xi_0=2(2\Delta+\sqrt{\delta_0})\): if \(|\xi|\ge\xi_0\) then
\(|\tau|\ge\tfrac12|\xi|\); if \(|\xi|<\xi_0\) then
\(|\tau|\ge\sqrt{R^2-\xi_0^2}\). *Condition:*
\(R\ge\tfrac{\sqrt5}{2}\xi_0\). Then \(|\tau|\ge\tfrac12|\xi|\) in both
cases.

**(C1) Outer field near a cone.** Take \(\tau\ge\tfrac12|\xi|\), \(r\ge R\),
\(|n^2|\le\delta_0\). Then \(n_0\ge r/\sqrt{17}-\Delta\) and
\(|\vec n|\le n_0+\sqrt{\delta_0}\), so

\[
n^2+u^-\!\cdot n\ \ge\ -\delta_0+(\tau_0-2\Delta)\Big(\frac{R}{\sqrt{17}}-\Delta\Big)-\Delta\sqrt{\delta_0}.
\]

*Condition:* the right side is positive. Then
\(V_{\rm out}\!\cdot n=-\kappa'(n^2+u^-\!\cdot n)<0\); the past case is the
mirror image. With (A), \(\operatorname{Im}x_{0i}^2\le0\) near every cone in
all three zones, strictly so where \(g>0\).

**(C3) The spacelike wedge.** For \(|\tau|\le\tfrac14|\xi|\) and \(r\ge R\),
\(|\xi|\ge4R/\sqrt{17}\) and
\(n^2\le(\tfrac14|\xi|+\Delta)^2-(|\xi|-\Delta)^2\). *Condition:* this is
below \(-\delta_0\) at \(|\xi|=4R/\sqrt{17}\). So no cone comes within
\(\delta_0\) where the outer field is switched off.

**(D) Phases away from the cones.** Frames \(j\ne i+2\) contribute a
non-positive imaginary part everywhere, which keeps the argument in
\([-\pi,0]\). The frame \(i+2\) acts only where \(|x_{0i}^2|\ge\delta_0\) and
\(r\le2R\), and there \(|2V_{i+2}\!\cdot n|\le2\lambda_0(2R+\Delta)^2\), while
the real part moves by at most \(|V|^2\le v^2\) with
\(v=4\lambda_0(2R+\Delta)+\kappa'(2R+\tau_0+\Delta)\). *Condition:*
\(v^2<\delta_0/2\). The excursion from this term is at most
\(\arctan[2\lambda_0(2R+\Delta)^2/(\delta_0-v^2)]\).

For the outer field with \(\Theta=1\),
\(x_{0i}^2\mapsto(1-i\kappa')^2n^2-2i\kappa'(1-i\kappa')\,n\!\cdot u-\kappa'^2u^2\).
The imaginary part is positive only if \(n^2+n\!\cdot u<0\). If
\(n\!\cdot u\ge0\) then \(n^2<0\), the imaginary part is at most
\(2\kappa'|n^2|\) and the real part at most \(-(1-\kappa'^2)|n^2|\): excursion
at most \(\arctan[2\kappa'/(1-\kappa'^2)]\). If \(n\!\cdot u<0\) then
\(|\vec n|\ge2n_0\), \(|n^2|\ge\tfrac34|\vec n|^2\) and
\(|n\!\cdot u|\le\Delta|\vec n|\): excursion at most
\(\arctan[2\kappa'(1+2\Delta/(3(R/\sqrt{17}-\Delta)))]\). Where
\(0<\Theta<1\), \(|\tau|\le\tfrac12|\xi|\), so
\(-n^2\ge\tfrac34|\xi|(|\xi|-4\Delta)\) and \(|n|\le\tfrac{\sqrt5}{2}|\xi|+\Delta\):
excursion at most
\(\arctan\{2\kappa'[1+(\tfrac{\sqrt5}{2}+\Delta/|\xi|)(\tau_0+\Delta)/(\tfrac34(|\xi|-4\Delta))]\}\)
with \(|\xi|\ge2R/\sqrt5\). *Condition:* \(2R/\sqrt5>4\Delta\). The total
excursion bound is the sum of the frame term and the largest of the three
far terms.

**(E) Local diffeomorphism.**
\(|\nabla m_j|\le\tfrac{3}{2\delta_0}\cdot2|x_0-x_{j+2}|\le3(2R+\Delta)/\delta_0\),
so

\[
\|\partial V_{\rm in}\|\le4\Big[\frac{3(2R+\Delta)^2\lambda_0}{\delta_0}+\lambda_0\Big].
\]

For the outer field \(|\nabla\Theta_\pm|\le6|\nabla(\tau/|\xi|)|\le7.5/r\) on
its support, so
\(\|\partial V_{\rm out}\|\le\kappa'[1+7.5(1+(\tau_0+\Delta)/R)]\). The blend
adds \(\tfrac3{2R}[4\lambda_0(2R+\Delta)+\kappa'(2R+\tau_0+\Delta)]\).
*Condition:* the sum is below 1.

**(F) Lower bounds near a cusp.** The frame at the cusp has weight 1 there.
If \(\rho\ge\tfrac12|\alpha\beta D|\) the imaginary part of \(x_{0,j}^2\) is
at least \(2\nu\rho\) in modulus; otherwise the real part is at least
\(\tfrac12|\alpha\beta D|-|V|^2\). The other two invariants are treated in
the same way. This gives item 4 with \(c=\tfrac12\) up to the negligible
\(|V|^2\).

**The conditions, evaluated.**

| Condition | \((-,-)\) | \((+,-)\), \((-,+)\) |
| --- | --- | --- |
| (H1) \(\epsilon_*\) against \(2\delta_0=0.05\) | 0.0897 | 0.345 |
| (C0) \(\tau_0\ge3\Delta\), with \(\Delta\) | 0.988 | 1.034 |
| (C1) margin (must be positive) | 5.55 | 5.19 |
| (C2), (C3), wedge margin | hold | hold |
| (D) excursion bound, rad | 0.083 | 0.049 |
| (E) bound on \(\|\partial V\|\) | 0.368 | 0.179 |

`certificate` computes these from the geometry and the constants, and
`test_every_explicit_condition_of_the_proof_holds` asserts all of them in
each region. A further test checks that the measured excursion and Jacobian
norm are below the proved bounds.

With an earlier choice of constants (\(\kappa=\nu=3\times10^{-6}\),
\(\delta_0=0.05\)) the sampled properties also held, but conditions (H1),
(D) and (E) did not: the bounds here are not sharp, and the constants were
reduced until every condition holds.

**Every square.** Up to Poincaré transformations and an overall scale, a
light-like square is determined by its kinematic region and one angle, that
is by the cross-ratio \(s/t\). `tuned_contour` chooses the constants from
the geometry in a fixed order: \(\tau_0=8\Delta\); \(\delta_0=\epsilon_*/2.5\);
the smallest \(R\) meeting (C1) to (C3); \(\kappa'\) small enough for (D);
and \(\kappa=\nu\) such that the bound in (E) equals 0.3. Each step can
always be carried out: \(\epsilon_*>0\) for a non-degenerate square because
the right side of (H1) tends to zero with \(\epsilon\); (C1) to (C3) hold for
large \(R\); (D) and (E) hold for small amplitudes. So the conditions can be
met for every non-degenerate square. Evaluated on 30 angles in each region:

| Region | \(s/t\) covered | \(\epsilon_*\) | \(R\) | \(\kappa=\nu\) | (E) | (D), rad |
| --- | --- | --- | --- | --- | --- | --- |
| \((-,-)\) | 0.0057 to 173 | 0.0030 to 0.10 | 7.0 to 8.9 | \(10^{-10}\) to \(4\times10^{-7}\) | 0.300 | 0.057 to 0.069 |
| \((+,-)\) | \(-174\) to \(-1.01\) | 0.0031 to 0.62 | 8.0 to 8.9 | \(10^{-10}\) to \(4\times10^{-6}\) | 0.300 | 0.057 to 0.069 |
| \((-,+)\) | \(-0.994\) to \(-0.0057\) | 0.0031 to 0.62 | 8.0 to 8.9 | \(10^{-10}\) to \(4\times10^{-6}\) | 0.300 | 0.057 to 0.069 |

All conditions hold at all 90 squares, and for the same squares rescaled by
0.01 and by 37. The amplitudes shrink as the square degenerates
(\(s/t\to0,\infty\)), as they must: the four-cone gap closes. The sampled
properties were also checked at four extreme squares.

## 10. Absolute integrability, line by line

Throughout, \(\eta=-\epsilon>0\) is small, the weight satisfies
\(|f(L)|\le A(1+|L|)^k\), and in \(D=4+2\eta\) dimensions the transverse space
has dimension \(2+2\eta\), so the transverse measure is
\(\rho^{\eta}d\rho\) times a finite angular measure. The three invariants of
a frame depend on the transverse coordinates only through \(\rho\); the
frame field is \(\nu y_\perp\) in every transverse direction and Lemma 5 is
unchanged.

**Lemma 6.** Let \(r\le R\) and \(|x_{0i}^2|\le\delta_0\). Then

\[
|V^2|\le\frac\theta2\,\big|\operatorname{Im}x_{0i}^{2\,\rm def}\big|,
\qquad
\big|x_{0i}^{2\,\rm def}\big|\ \ge\ \tfrac12\max\Big(|x_{0i}^2|,\ \big|\operatorname{Im}x_{0i}^{2\,\rm def}\big|\Big),
\qquad
\theta=\max\Big(\frac{32\kappa P^2}{\min_j|D_j|},\ 8\nu B^2\Big),
\]

where \(B\) is the largest Euclidean norm of a transverse vector with
\(y_\perp^2=-1\), and "def" denotes the value on the contour. Moreover
\(|\operatorname{Im}x_{0i}^{2\,\rm def}|\ge2\nu\rho_{j}\) for any frame \(j\)
with \(m_j=1\) that controls \(i\).

*Proof.* Near the cone of \(x_i\) the frame \(i+2\) is off, so every active
frame controls \(i\) and, with
\(b_j=\tfrac{\kappa|D_j|}{2}|d_j|^2+\nu\rho_j\), Lemma 5 gives
\(|\operatorname{Im}x_{0i}^{2\,\rm def}|=2|V\!\cdot n|\ge2\sum_jm_jb_j\). On
the other hand \(|d_\beta p_{j-1}+d_\alpha p_j|\le\sqrt2P|d_j|\) and
\(|y_\perp|\le B\sqrt{\rho_j}\), so
\(|V_j|^2\le4\kappa^2P^2|d_j|^2+2\nu^2B^2\rho_j\). At most four frames are
active and \(m_j\le1\), so

\[
|V^2|\le|V|^2\le4\sum_jm_j|V_j|^2
\le\sum_jm_j\Big[\frac{32\kappa P^2}{|D_j|}\cdot\frac{\kappa|D_j|}{2}|d_j|^2+8\nu B^2\cdot\nu\rho_j\Big]
\le\theta\sum_jm_jb_j .
\]

This is the first claim. The real part on the contour is
\(x_{0i}^2-V^2\). If \(|\operatorname{Im}|\ge|x_{0i}^2|\) the modulus is at
least \(|\operatorname{Im}|\); otherwise it is at least
\(|x_{0i}^2|-\tfrac\theta2|\operatorname{Im}|\ge(1-\tfrac\theta2)|x_{0i}^2|\).
With \(\theta\le1\) both are at least half the maximum.
\(\blacksquare\)

For the contours used \(\theta\) is between \(10^{-6}\) and \(3\times10^{-5}\).
On samples the ratio \(|V^2|/|\operatorname{Im}|\) never exceeds
\(\theta/2\) and the second ratio is 1.000.

**Strips.** In the frame of \(x_j\) let

\[
S_j^+=\{-a_1\le\beta\le1-b_0,\ |\alpha|\le a_1,\ \rho\le\rho_1\},
\]

and let \(S_j^-\) be the same with the roles of \(\alpha\) and \(-\beta\)
exchanged (the strip around the other edge at \(x_j\)). Choose
\(b_0=4\delta_0/\min_j|D_j|\), and \(a_1,\rho_1\) small enough that on the
strip \(|x_{0,j}^2|\le\delta_0\), \(|x_{0,j+1}^2|\le\delta_0\) and
\(|x_{0,j+2}^2|\ge2\delta_0\); explicitly
\(a_1=\delta_0/(8\max|D|)\) and
\(\rho_1=\delta_0^2/(256\max(1,\max|D|\,|u|/\min|D|))\). Then
\(b_0<\tfrac12\), so \(S_j^+\) and \(S_{j+1}^-\) together cover a
neighbourhood of the edge from \(x_j\) to \(x_{j+1}\), and the eight strips
cover a neighbourhood of all four edges.

**Lemma 7.** On \(S_j^+\): the frame \(j\) is fully on, and

\[
|x_{0,j}^{2\,\rm def}|\ge\nu\rho,\qquad
|x_{0,j+1}^{2\,\rm def}|\ge\frac\nu4\max\Big(\rho,\frac{|\alpha D|b_0}{2}\Big),\qquad
|x_{0,j-1}^{2\,\rm def}|\ge\frac\nu4\max\Big(\rho,\frac{|\beta D|}{2}\Big),\qquad
|x_{0,j+2}^{2\,\rm def}|\ge\delta_0 .
\]

*Proof.* \(m_j=\chi(|x_{0,j+2}^2|/\delta_0)=1\). The frames \(j-1\) and
\(j+2\) are off because \(x_{0,j+1}^2\) and \(x_{0,j}^2\) are within
\(\delta_0\) of zero; the frame \(j+1\) may be partly on and controls
\(x_{0,j}^2\) and \(x_{0,j+1}^2\). So Lemma 6 applies to \(i=j\) and
\(i=j+1\) with \(|\operatorname{Im}|\ge2\nu\rho\). This gives the first
bound. For the second, \(x_{0,j+1}^2=\alpha D(\beta-1)-\rho\): if
\(\rho\le\tfrac12|\alpha D|(1-\beta)\) its modulus is at least
\(\tfrac12|\alpha D|b_0\), otherwise \(2\nu\rho\ge\nu|\alpha D|b_0\); with
Lemma 6 the modulus on the contour is at least half of
\(\max(|x_{0,j+1}^2|,2\nu\rho)\). For the third, if
\(|x_{0,j-1}^2|\le\delta_0\) the same argument applies with
\(x_{0,j-1}^2=\beta D(1+\alpha)-\rho\); if not, the modulus on the contour is
at least \(\delta_0-v^2\ge\delta_0/2\), which exceeds the right side because
\(\nu\) is small. The fourth follows from \(|x_{0,j+2}^2|\ge2\delta_0\) and
\(v^2<\delta_0/2\). \(\blacksquare\)

On samples of the strips the four left sides exceed the right sides by
factors between 3.9 and \(5\times10^7\).

**Lemma 8 (the strip integral).** For \(\eta'>0\),

\[
\int_{S_j^+}\frac{\rho^{\eta'}\,d\alpha\,d\beta\,d\rho}
{\rho\,\max(\rho,|\alpha|)\,\max(\rho,|\beta|)}
\le C\int_0^{\rho_1}\rho^{\eta'-1}\big(1+\log\tfrac1\rho\big)^2d\rho<\infty .
\]

*Proof.* \(\int_{-1}^{1}d\alpha/\max(\rho,|\alpha|)=2(1+\log(1/\rho))\) for
\(\rho\le1\), and the same for \(\beta\). \(\blacksquare\)

**Consequence on the strips.** By Lemma 7 the integrand is at most
\(C\,\nu^{-3}\delta_0^{-1}\rho^{\eta}(1+|L|)^k/[\rho\max(\rho,|\alpha|)\max(\rho,|\beta|)]\).
Each invariant on the contour lies between \(\nu\rho\) and a constant in
modulus, so \(|L|\le C(1+\log(1/\rho))\) and
\((1+|L|)^k\le C\rho^{-\eta/2}\). Lemma 8 with \(\eta'=\eta/2\) gives
convergence.

**The rest of the ball.** Let \(K\) be the set \(r\le2R\) with the open
strips removed. On \(K\) no invariant vanishes: near a cone and inside \(R\)
by (A), where equality in Lemma 5 forces the point onto an edge; near a cone
in the blend zone by (C1); elsewhere because the modulus is at least
\(\delta_0-v^2\). \(K\) is compact and the invariants are continuous, so
they are bounded below on \(K\) and the integrand is bounded.

**Lemma 9 (infinity).** For \(r\ge2R\) and \(y=x_0-c\),

\[
|x_{0i}^{2\,\rm def}|\ \ge\ c_\infty\,\kappa'\max\big(r,|y^2|\big).
\]

*Proof.* On this zone \(V=-\kappa'\Theta(n+u)\), so
\(V^2=\kappa'^2\Theta^2(n^2+2n\!\cdot u+u^2)\) exactly. Let

\[
c_u=(\tau_0-2\Delta)\Big(\frac{1}{\sqrt{17}}-\frac{\Delta}{2R}\Big)-\Delta\sqrt{\frac{\omega}{2R}},
\qquad\omega=\tfrac12c_u ,
\]

which is positive for the constants used (lengths are in units in which
\(\Delta\) is of order one). If \(|n^2|\ge\omega r\), the real part on the
contour is at least
\(|n^2|(1-\kappa'^2)-\kappa'^2(2|n||u|+|u|^2)\ge\tfrac12|n^2|\), because
\(\kappa'^2(2|n||u|+|u|^2)\le\tfrac14\omega r\) for small \(\kappa'\). If
\(|n^2|<\omega r\), the argument of (C2) with \(\delta_0\) replaced by
\(\omega r\) gives \(\Theta=1\), and that of (C1) gives
\(n\!\cdot u\ge(\tau_0-2\Delta)(r/\sqrt{17}-\Delta)-\Delta\sqrt{\omega r}\ge c_ur\),
so \(|\operatorname{Im}|=2\kappa'|n^2+n\!\cdot u|\ge2\kappa'(c_u-\omega)r=\kappa'c_ur\).
In both cases the modulus is at least a constant times
\(\kappa'\max(r,|n^2|)\). Finally \(|n^2-y^2|\le3\Delta r\), so
\(\max(r,|n^2|)\ge\max(r,|y^2|)/(1+3\Delta)\). \(\blacksquare\)

On samples the ratio is at least 5 in units of \(\kappa'\).

With \(y=(\tau,\xi)\), \(|y^2|=\big||\tau|-|\xi|\big|(|\tau|+|\xi|)\ge|\varsigma|\,r\)
where \(\varsigma=|\tau|-|\xi|\). Hence the integrand is at most
\(C\,r^{-4}\max(1,|\varsigma|)^{-4}\) times powers of \(\log r\) and the
factor \(r^{2\eta}\) from the transverse measure, and

\[
\int_{r\ge2R}\frac{d\tau\,d^3\xi\;r^{2\eta}(1+\log r)^k}{r^4\max(1,|\varsigma|)^4}
\le C\int^\infty\frac{\xi^2\,d\xi}{\xi^{4-2\eta}}(1+\log\xi)^k<\infty
\]

for \(\eta<\tfrac12\), the \(\tau\) integral at fixed \(\xi\) being bounded
by \(2\int d\varsigma/\max(1,|\varsigma|)^4\).

Strips, the rest of the ball and infinity exhaust the contour, and the
Jacobian is bounded by (E). So the regulated integrand is absolutely
integrable. \(\blacksquare\)

## 11. The regulated integral on the contour

**Correction: the extra dimensions are displaced too.** In \(D=4+k\)
dimensions the point is \((x_0,e)\) with \(e\in\mathbb R^k\), and every
invariant is \(X_i=x_{0i}^2-e\!\cdot e\). An earlier version of this section
treated \(\mu=e\!\cdot e\) as a real shift and used the four-dimensional
contour for every \(\mu\). That is wrong once \(\mu>\delta_0\): the frame
field of the opposite cusp is then switched on near the surface
\(x_{0i}^2=\mu\) and can move it to the wrong side. On 19 865 points of that
surface at \(\mu=5\delta_0\), 994 are displaced to the wrong side.

The remedy is to carry out the same construction in
\(\mathbb R^{1,3+k}\). The transverse space of each frame gains the
coordinates \(e\), so \(\rho\to\rho+e\!\cdot e\); the frame field gains the
component \(\nu e\), the far fields the component \(-\kappa'e\); the cutoffs
use the \(D\)-dimensional invariants. Lemma 5 and sections 4, 9 and 10 hold
verbatim, with Euclidean norms taken in \(\mathbb R^{4+k}\); section 10 was
already written with the transverse measure of this picture. With this
contour no point of any surface \(x_{0i}^2=\mu\) is displaced to the wrong
side (0 of 19 865 at \(\mu=5\delta_0\), and 0 at \(\mu=1,4,50\)), and the
survey of section 5 in \(4+1\) dimensions gives the same results, with
near-cone points at \(\mu\) up to \(4.5\times10^4\). Non-integer \(k=2\eta\)
is reached in the usual way: the integrand and the field depend on \(e\)
through \(|e|\) and a radial factor, so the integral reduces to one over
\(x_0\) and \(m=|e|\) with weight \(m^{2\eta-1}\), which is analytic in
\(\eta\).

**Proposition 8 (parametric representation).** Let \(a_i>0\) with
\(\sum a_i=4\) and \(0<\eta<\tfrac12\). The regulated integral of
\(\prod(-X_i)^{-a_i}\) over the contour equals

\[
-\,i\pi^{2+\eta}\,\frac{\Gamma(2-\eta)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)\,
\prod\alpha_i^{a_i-1}\Big(-\sum_{i<k}\alpha_i\alpha_kx_{ik}^2\Big)^{\eta-2}.
\]

*Derivation.* For \(0<\varphi\le\pi/2\) replace every squared time
difference by \(e^{-2i\varphi}\) times itself, in the integration variable
and in the external points alike. With \(\mu=e\!\cdot e\),

\[
Y_i(\varphi)=\mu+|\xi-\xi_i|^2-e^{-2i\varphi}(\tau-t_i)^2 .
\]

*Step 1 (a proper cone).* Each \(Y_i(\varphi)\) is a non-negative
combination of \(1\) and \(e^{i(\pi-2\varphi)}\), so all four lie in a closed
cone of opening \(\pi-2\varphi<\pi\), and
\(|Y_i(\varphi)|\ge\sin\varphi\,(\mu+|x_0-x_i|^2)\). The integral over real
\((x_0,e)\) is absolutely convergent: the singularity at \(x_0=x_i\),
\(\mu=0\) is integrable because \(2a_i<4\), and the decay at infinity is
\(|x|^{-8}\).

*Step 2 (Feynman parameters).* For numbers in a common proper cone,
\(\prod Y_i^{-a_i}=\frac{\Gamma(4)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta(1-\sum\alpha)\prod\alpha_i^{a_i-1}(\sum\alpha_iY_i)^{-4}\),
and \(|\sum\alpha_iY_i|\ge\sin\varphi\sum\alpha_i(\mu+|x_0-x_i|^2)\), so
the double integral converges absolutely and the order can be exchanged.

*Step 3 (completing the square).*
\(\sum\alpha_iY_i=\mu+A_\alpha(\varphi)+|\xi-\bar\xi|^2-e^{-2i\varphi}(\tau-\bar t)^2\)
with \(A_\alpha(\varphi)=-\sum_{i<k}\alpha_i\alpha_kx_{ik}^2(\varphi)\) and
\(x_{ik}^2(\varphi)=e^{-2i\varphi}(t_i-t_k)^2-|\xi_i-\xi_k|^2\).

*Step 4 (rotating the time).* Put \(\tau-\bar t=e^{i\psi}\sigma\) and turn
\(\psi\) from \(0\) to \(\varphi-\pi/2\). The term
\(-e^{-2i\varphi}(\tau-\bar t)^2\) has argument \(\pi+2\psi-2\varphi\), which
moves from \(\pi-2\varphi\) to \(0\): it stays in the cone, the denominator
never vanishes, and the rotation is allowed. The integral is then Euclidean,
\(\int d^4y\,(c+y^2)^{-4}=\pi^2/(6c^2)\), with Jacobian \(-ie^{i\varphi}\).
Integrating the extra dimensions,
\(\frac{\pi^\eta}{\Gamma(\eta)}\int_0^\infty d\mu\,\mu^{\eta-1}(A+\mu)^{-2}=\pi^\eta\Gamma(2-\eta)A^{\eta-2}\).
Hence, for \(\varphi>0\),

\[
G(\varphi):=\int d^Dx\prod Y_i(\varphi)^{-a_i}
=-\,ie^{i\varphi}\,\pi^{2+\eta}\frac{\Gamma(2-\eta)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)
\prod\alpha_i^{a_i-1}A_\alpha(\varphi)^{\eta-2}.
\tag{$*$}
\]

Every integral in steps 1 to 4 is absolutely convergent.

*Step 5 (the right side as \(\varphi\to0\)).* With spacelike diagonals
\(\operatorname{Re}A_\alpha(\varphi)\ge A_\alpha(0)=\alpha_1\alpha_3|s|+\alpha_2\alpha_4|t|\ge0\),
because \(\cos2\varphi\le1\). Since \(\eta-2<0\),
\(|A_\alpha(\varphi)|^{\eta-2}\le A_\alpha(0)^{\eta-2}\), which is integrable
against \(\prod\alpha_i^{a_i-1}\) for \(\eta>0\): this is the convergent
regulated Feynman-parameter integral. By dominated convergence the right
side of \((*)\) tends to the stated formula.

*Step 6 (the left side as \(\varphi\to0\)).* Three lemmas. Write
\(n=(x_0-x_i,e)\), \(V=(V_4,V_e)\) for the displacement in \(D\) dimensions,
\(|\cdot|\) for Euclidean norms and \(V\!\cdot n\) for the Minkowski product;
\(Y_i(\varphi,\lambda)\) is \(Y_i(\varphi)\) at the point \((x_0,e)+i\lambda V\).
\(J<1\) is the bound (E) on \(\|\partial V\|\).

**Lemma 10 (no zeros, and a bound uniform in \(\varphi\)).**

(a) \(\operatorname{Im}\big(e^{i\varphi}Y_i(\varphi,\lambda)\big)
=\sin\varphi\,\big(|n|^2-\lambda^2|V|^2\big)-2\lambda\cos\varphi\;V\!\cdot n\).

(b) \(|V|\le J|n|\) for every \(i\).

(c) For \(0<\varphi\le\varphi_0=\delta_0/(8(2R+\Delta+1)^2)\) and
\(0\le\lambda\le1\), no \(Y_i(\varphi,\lambda)\) vanishes.

(d) For \(0\le\varphi\le\varphi_0\), at \(\lambda=1\), in the inner zone
with \(|X_i|\le\delta_0\):
\(|Y_i(\varphi,1)|\ge c_J\max\big(|X_i|,|\operatorname{Im}X_i^{\rm def}|\big)\)
with \(c_J=\tfrac12\min\big(\tfrac14,\tfrac{1-J^2}{16J}\big)\).

*Proof.* (a) Expand
\(Y=|\vec n+i\lambda\vec V|^2-e^{-2i\varphi}(n_0+i\lambda V_0)^2\), with
bilinear squares, multiply by \(e^{i\varphi}\) and take the imaginary part:
the terms in \(\sin\varphi\) collect to
\(|\vec n|^2-\lambda^2|\vec V|^2+n_0^2-\lambda^2V_0^2\) and those in
\(\cos\varphi\) to \(2\lambda(\vec n\!\cdot\!\vec V-n_0V_0)=-2\lambda V\!\cdot n\).

(b) \(V\) vanishes at each cusp \(x_i\): the three frames that contain
\(x_i\) in their flat triangle have \(V_j(x_i)=0\), and the fourth is
switched off there. \(V\) is \(J\)-Lipschitz, so
\(|V(x)|=|V(x)-V(x_i)|\le J|x-x_i|\).

(c) By (a) and (b), if \(V\!\cdot n\le0\) then
\(\operatorname{Im}(e^{i\varphi}Y)\ge\sin\varphi\,(1-J^2)|n|^2>0\). Suppose
\(V\!\cdot n>0\). In the inner and blend zones this happens only where
\(|X_i|\ge\delta_0\), since near the cone every active field is admissible.
Write \(Y=-X_i^{\lambda}+(1-e^{-2i\varphi})(n_0+i\lambda V_0)^2\), where
\(X_i^\lambda\) is the invariant at the displaced point. The first term has
modulus at least \(\delta_0-v^2\ge\delta_0/2\); the second at most
\(2\varphi(2R+\Delta+1)^2\le\delta_0/4\). So \(|Y|\ge\delta_0/4\). In the
outer zone \(V\!\cdot n=-\kappa'(n^2+n\!\cdot u)>0\) means
\(n^2<-n\!\cdot u\). Then
\(\operatorname{Re}Y=-(n^2-\lambda^2V^2)+(1-\cos2\varphi)(n_0^2-\lambda^2V_0^2)-2\lambda\sin2\varphi\,n_0V_0\).
The last two terms are non-negative because \(V_0=-\kappa'(n_0+u_0)\) has
the sign opposite to \(n_0\) and \(|V_0|<|n_0|\). For the first,
\(V^2=\kappa'^2(n^2+2n\!\cdot u+u^2)\): if \(n\!\cdot u\ge0\) then
\(n^2<0\) and \(n^2-\lambda^2V^2<0\) once \(|n^2|\ge n\!\cdot u\ge c_ur\)
exceeds \(\kappa'^2(2n\!\cdot u+u^2)\); if \(n\!\cdot u<0\) then
\(n^2\le-\tfrac34|\vec n|^2\) and the same holds. So \(\operatorname{Re}Y>0\).

(d) Put \(X=X_i\), \(I=|\operatorname{Im}X_i^{\rm def}|=2|V\!\cdot n|\); here
\(V\!\cdot n\le0\). From (a),
\(\operatorname{Im}(e^{i\varphi}Y)\ge\sin\varphi(1-J^2)|n|^2+\cos\varphi\,I\).
So \(|Y|\ge\cos\varphi\,I\ge\tfrac12I\). If \(I\ge\tfrac12|X|\) this is at
least \(\tfrac14|X|\). Otherwise use
\(\operatorname{Re}(e^{i\varphi}Y)=-\cos\varphi\operatorname{Re}X_i^{\rm def}+\sin\varphi\operatorname{Im}X_i^{\rm def}-4\sin\varphi\,n_0V_0\).
By Lemma 6, \(|\operatorname{Re}X_i^{\rm def}|\ge|X|-\tfrac\theta2I\ge\tfrac34|X|\),
and \(|4n_0V_0|\le4J|n|^2\). If \(4J\sin\varphi|n|^2\le\tfrac18|X|\) the real
part is at least
\(\tfrac38|X|-\tfrac14|X|\cdot\tfrac12-\tfrac18|X|=\tfrac18|X|\) in modulus.
If not, \(\sin\varphi|n|^2>|X|/(32J)\) and the imaginary part is at least
\((1-J^2)|X|/(32J)\). \(\blacksquare\)

**Lemma 11 (independence of \(\lambda\)).** For \(0<\varphi\le\varphi_0\)
the integral of \(\prod Y_i^{-a_i}\) over \((x_0,e)+i\lambda V\) does not
depend on \(\lambda\in[0,1]\).

*Proof.* The integrand is a holomorphic function \(h\) on the open set
where no \(Y_i\) vanishes, with its powers continued from \(\lambda=0\),
where the \(Y_i\) lie in the cone of step 1. The form \(h\,dz_1\cdots dz_D\)
is closed. The map \((\lambda,x)\mapsto x+i\lambda V(x)\) is Lipschitz and
by Lemma 10(c) takes \([0,1]\times\mathbb R^D\) into that set. Apply Stokes'
theorem on \([0,1]\times\{|x|\le\rho\}\): the difference of the integrals
at \(\lambda=1\) and \(\lambda=0\) over the ball is the integral over the
side \([0,1]\times\{|x|=\rho\}\). On the side, for fixed \(\varphi>0\), every
\(|Y_i|\ge c(\varphi)\rho^2\): where \(V\!\cdot n\le0\) by Lemma 10(a) and
(b); where \(V\!\cdot n>0\), either
\(\sin\varphi(1-J^2)|n|^2\ge4\kappa'|n^2|\) and (a) still gives
\(\tfrac12\sin\varphi(1-J^2)\rho^2\), or \(|n^2|\) is at least a multiple of
\(\sin\varphi\,\rho^2\) and \(|Y|\ge\operatorname{Re}Y\ge\tfrac12|n^2|\) by
the proof of (c). The side has area of order \(\rho^{D-1}\), the form
contributes a factor \(|V|\le\kappa'(\rho+\tau_0+\Delta)\), and
\(|h|\le C\rho^{-8}\). The side integral is of order \(\rho^{D-8}\) and tends
to zero. For non-integer \(D\) the same argument applies to the reduced
variables \((x_0,m)\) with the factor \(m^{2\eta-1}\), which vanishes with
the boundary at \(m=0\). \(\blacksquare\)

**Lemma 12 (the limit).** As \(\varphi\to0\), \(G(\varphi)\) tends to the
integral of \(\prod(-X_i)^{-a_i}\) over the contour.

*Proof.* By Lemma 11, \(G(\varphi)\) equals the integral at \(\lambda=1\).
There the integrand converges pointwise as \(\varphi\to0\). It is dominated
uniformly in \(\varphi\le\varphi_0\): near the cones in the inner zone by
Lemma 10(d), which gives the bound of Lemma 6 with \(c_J\) in place of
\(\tfrac12\), so that Lemmas 7 and 8 apply with a different constant; away
from the cones in the ball by \(|Y|\ge\delta_0/4\); in the outer zone by
\(|Y|\ge\cos\varphi\cdot2\kappa'|n^2+n\!\cdot u|\) where \(V\!\cdot n\le0\)
and by \(|Y|\ge\tfrac12|n^2|\) otherwise, which reproduce Lemma 9. The
dominating function is integrable by section 10. \(\blacksquare\)

Steps 5 and 6 give the proposition. \(\blacksquare\)

*Timelike diagonal.* With \(s>0\), \(A_\alpha(0)\) can vanish inside the
simplex and step 5 is a boundary value: the right side at \(\varphi\to0^+\)
is the formula with \(s\to s-i0\). Step 6 does not use the sign of the
diagonals. The dominated-convergence argument of step 5 is not adapted to
that case here.

*Numerical checks.*

| Check | Result |
| --- | --- |
| \(|Y_i|\ge\sin\varphi(\mu+|x_0-x_i|^2)\), \(0\le\arg Y_i\le\pi-2\varphi\) | holds at 24 000 points, three angles |
| Both sides of \((*)\) at \(\varphi=\pi/2\), \(\mu=0.5\), unit powers | \(1.4317\pm0.0008\) against \(1.43230\) |
| Both sides at \(\varphi=1\), unit powers | \(1.0225-1.3751i\) against \(1.0235-1.3763i\) |
| Both sides at \(\varphi=1\), powers \((1.2,0.8,1.2,0.8)\) | \(1.1356-1.4812i\) against \(1.1364-1.4820i\) |
| Right side at \(\varphi=0.1,0.01,0.001\) | approaches the Minkowski formula linearly |
| Lemma 10(a), relative error of the identity | \(6\times10^{-12}\) |
| Lemma 10(b), largest \(|V|/|n|\) | 0.0032, against \(J=0.30\) |
| Lemma 10(c), smallest \(|Y_i|/(\sin\varphi|n|^2)\), also for \(\varphi\) up to 0.3 | 1.00 |
| Lemma 10(d), smallest \(|Y_i|/\max(|X_i|,|\operatorname{Im}X_i^{\rm def}|)\) | 1.00, against \(c_J=0.095\) |

The left side of \((*)\) is a four-dimensional quasi-Monte Carlo integral at
fixed \(\mu\); the agreement is within its error. The checks of Lemma 10 are
in \(4+1\) dimensions, in all three kinematic regions.

**Feynman parameters to Mellin–Barnes.** The Feynman-parameter integral for
\(z^p\) (powers \(1\pm p\)) is evaluated directly, as three nested
one-dimensional integrals, and compared with the one-fold Mellin–Barnes
integral used in GATE_6_FUNCTIONAL and in Lemma 4:

| \(p\) | \(\eta\) | \(s/t\) | Feynman parameters | Mellin–Barnes |
| --- | --- | --- | --- | --- |
| 0.2 | 0.4 | 1 | 3.520361875 | 3.520361876 |
| \(-0.3\) | 0.3 | 2.5 | 3.893274… | 3.893274816 |
| 0 | 0.25 | 1 | 4.004157042 | 4.004157043 |

They agree to \(10^{-7}\) or better, normalisation included.

So the regulated integral on the contour equals the Mellin–Barnes form.

## 12. What remains

- With a timelike diagonal, step 5 of Proposition 8 is a boundary value in
  the Mandelstam invariant, and its dominated-convergence argument is not
  adapted to that case. Steps 1 to 4 and 6 hold in all regions.
- The equality of the Feynman-parameter and Mellin–Barnes forms is checked
  numerically to \(10^{-7}\); it is a standard identity and is not proved
  here.
- In section 10 the constant \(c_u>0\) of Lemma 9 and the choice of
  \(a_1,\rho_1\) in the strips are verified for the contours used, by the
  sampled checks, and are not added to the certificate.
- The non-integer dimension is handled by the usual reduction to the radial
  variable of the extra dimensions; the corresponding Stokes argument is
  indicated in Lemma 11, not written separately.
- Review by someone who works with these integrals.

## 13. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_contour.py tests/test_schwinger.py
PYTHONPATH=src python -m examples.feynman_contour_survey
PYTHONPATH=src python -m geocusp.solver contour --points 5000
```
