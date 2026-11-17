# Gate 6: an explicit integration contour

[GATE_6_BOUNDS.md](GATE_6_BOUNDS.md) left one step open: a global integration
contour on which the regulated integrand is absolutely integrable. It was
needed twice, for the exchange of the Mellin superposition with the loop
integral and for the statement that only the edges produce poles. This
document constructs such a contour, and replaces the appeal to Stirling's
formula in Lemma 4 by elementary bounds.

**Status.** The construction is explicit and implemented, in all three real
kinematic regions. Its central inequality is proved (Lemma 5). Its global
properties are reduced in section 9 to a list of inequalities between named
constants, each proved, and the list is evaluated for the contour used
(`certificate`): every condition holds. The properties are also measured on
stress samples and lie within the proved bounds. Section 10 identifies the
regulated integral on the contour with the Mellin–Barnes form; one textbook
representation is cited there and not re-derived. Nothing here has been
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

## 10. The regulated integral on the contour

**Dimensional regularisation as a shift.** In \(D=4+2\eta\) dimensions with
\(\eta>0\) the integrand depends on the \(2\eta\) extra components of
\(x_0\) only through their squared length \(\mu\), and every invariant
becomes \(x_{0i}^2-\mu\). Integrating the extra directions,

\[
\int d^Dx_0\,(\cdots)=\frac{\pi^{\eta}}{\Gamma(\eta)}\int_0^\infty d\mu\,\mu^{\eta-1}
\int d^4x_0\,(\cdots)\big|_{x_{0i}^2\to x_{0i}^2-\mu}.
\]

The shift is real, so the contour of this document serves for every
\(\mu\ge0\): imaginary parts are unchanged, and the bounds of item 4 hold
with \(x_{0i}^2-\mu\) in place of \(x_{0i}^2\) by the same two-case
argument. For \(\mu>0\) the shifted invariants also vanish at real points
inside the flat triangles, where \(V=0\); these are isolated curves on which
a single invariant vanishes, with an integrable singularity.

**From the shifted integral to Feynman parameters.** The four-dimensional
integral with shifted invariants and powers \(a_i\) (\(\sum a_i=4\)) has the
Schwinger representation

\[
\int\frac{d^4x_0}{i\pi^2}\prod_i\frac{1}{(-x_{0i}^2+\mu-i0)^{a_i}}
=\frac{\Gamma(2)}{\prod\Gamma(a_i)}\int d^4\alpha\,\delta\big(1-\textstyle\sum\alpha\big)
\frac{\prod\alpha_i^{a_i-1}}{(\alpha_1\alpha_3|s|+\alpha_2\alpha_4|t|+\mu)^2}.
\]

This is the textbook parametric representation of a one-loop integral with
a common mass shift. **It is cited, not re-derived here**, and it is the one
input of this section taken from the literature.

**Extra dimensions.** With \(A=\alpha_1\alpha_3|s|+\alpha_2\alpha_4|t|\),

\[
\frac{1}{\Gamma(\eta)}\int_0^\infty d\mu\,\mu^{\eta-1}\frac{\Gamma(2)}{(A+\mu)^2}
=\Gamma(2-\eta)\,A^{\eta-2},
\]

which turns the right side into the \(D\)-dimensional Feynman-parameter
integral with exponent \(2-\eta=2+\epsilon\). Checked numerically.

**Feynman parameters to Mellin–Barnes.** The Feynman-parameter integral for
\(z^p\) (powers \(1\pm p\)) is evaluated directly, as three nested
one-dimensional integrals, and compared with the one-fold Mellin–Barnes
integral used in GATE_6_FUNCTIONAL and in Lemma 4:

| \(p\) | \(\eta\) | \(s/t\) | Feynman parameters | Mellin–Barnes |
| --- | --- | --- | --- | --- |
| 0.2 | 0.4 | 1 | 3.520361875 | 3.520361876 |
| \(-0.3\) | 0.3 | 2.5 | 3.893274… | 3.893274816 |
| 0 | 0.25 | 1 | 4.004157042 | 4.004157043 |

They agree to \(10^{-7}\) or better, normalisation included
(`test_feynman_parameter_form_equals_the_mellin_barnes_form`).

So the regulated integral on the contour equals the Mellin–Barnes form,
given the Schwinger representation.

## 11. What remains

- The Schwinger representation in section 10 is cited. Deriving it for the
  \(-i0\) integral means exchanging the parameter integrals with the
  \(x_0\) integral; on the contour of this document the integrand is
  absolutely integrable, which is what that exchange needs, but the
  derivation is not written out.
- The bounds of section 9 are proved for the stated construction. The
  absolute integrability in section 6 uses item 4 and a power-counting
  estimate that is given in outline.
- The certificate is evaluated for one square in each region. The
  conditions are inequalities between geometric constants and can be
  re-evaluated for any other square with `certificate`.
- Review by someone who works with these integrals.

## 12. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_contour.py
PYTHONPATH=src python -m examples.feynman_contour_survey
PYTHONPATH=src python -m geocusp.solver contour --points 5000
```
