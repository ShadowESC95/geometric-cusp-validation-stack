# Gate 6: bounds, regulators and the remaining cases

This document supplies the three items left open by
[GATE_6_CLOSURE.md](GATE_6_CLOSURE.md), section 8:

1. error bounds for the local statement;
2. the relation between the local statement and the double pole in
   dimensional regularisation;
3. the cusps, kinematic regions and weight that were covered there by
   symmetry only.

It also **corrects** one statement of that document (section 4 below).

**Status.** Item 3 is complete: it is computed, not inferred, and the
regularity of the two-loop weight at \(z=-1\) is proved exactly. Item 1 has
two lemmas with proofs and identifies a limitation of hard cutoffs. Item 2
rests on the published Mellin–Barnes evaluation, a uniform bound (Lemma 4)
and a pinch analysis (Proposition 5). The integration contour these need is
constructed in [GATE_6_CONTOUR.md](GATE_6_CONTOUR.md), which also gives an
elementary proof of Lemma 4. Nothing here has been reviewed by a specialist.

Code: `src/geocusp/lorentzian_edges.py`. Tests: `tests/test_gate6_bounds.py`
(71 tests). Example 62: `examples/lorentzian_all_cusps.py`.

## 1. Setting

Metric \((+,-,-,-)\); edges \(p_i=x_{i+1}-x_i\); \(s=x_{13}^2\),
\(t=x_{24}^2\), \(u=-s-t\). Near the cusp \(x_i\) write
\(x_0=x_i+\alpha p_{i-1}+\beta p_i+y_\perp\), \(\rho=|y_\perp|^2\), and let
\(D=2p_{i-1}\!\cdot p_i\) (equal to \(s\) at cusps 2, 4 and \(t\) at cusps
1, 3) and \(D'\) the other diagonal. Then **exactly**

\[
x_{0,i}^2=w,\qquad
x_{0,i-1}^2=w+c_-,\qquad
x_{0,i+1}^2=w+c_+,\qquad
x_{0,i+2}^2=D'(1+\alpha-\beta)+w+2\sqrt{\rho D'u/D}\,\cos\theta ,
\]

\[
w=\alpha\beta D-\rho,\qquad c_-=\beta D,\qquad c_+=-\alpha D,
\]

where \(\theta\) is the angle of \(y_\perp\) and \(D'u/D>0\). The coordinates
\((\alpha,\beta,\rho,\theta)\) are global, so this is the exact integrand
everywhere, not an expansion. The test
`test_local_formulas_are_exact_at_every_cusp` checks all four formulas
against explicit vectors at every cusp in every kinematic region.

Every invariant, including \(s\) and \(t\), carries \(-i0\). Hence

\[
\arg z=\pi\,(\tau_2+\tau_4+\tau_s-\tau_1-\tau_3-\tau_t),
\qquad \tau=1\ \text{if timelike},\ 0\ \text{if spacelike}.
\]

The kernel is
\(K_f(c_-,c_+)=\int_0^\infty d\rho\,f(L)/[w(w+c_-)(w+c_+)]\) with
\(w\to w-i0\), \(J=c_-c_+K_f\), and on the edge through \(x_i\) and
\(x_{i+1}\) the collinear invariant is \(c_+\) at fixed \(c_-=\beta D\).

## 2. Item 3: every cusp, every region, and the two-loop weight

**Kinematic regions.** Real null edges summing to zero have exactly one
positive invariant among \(s,t,u\). The regions are therefore
\((s,t)=(-,-),(+,-),(-,+)\).

**Result.** At each of the four cusps, on each of its two edges, in each of
the three regions:

- at a point of the edge segment the jump of \(J\) across the collinear
  invariant is \(-2\pi i\,M[F]\);
- on the same null line outside the segment the jump is zero;
- the constant weight gives exactly \(-2\pi i\) on the segment.

This is 24 edge/region combinations, each checked on and off the segment to
\(10^{-8}\) with the weight \(z^{0.3}\), which is not symmetric under
\(z\to1/z\) and so distinguishes the odd cusps.

The sheets visited differ between cusps and regions:

| Region | Cusps 2, 4 | Cusps 1, 3 |
| --- | --- | --- |
| \((-,-)\) | \(\arg z/\pi\in\{-2,-1,0,1\}\) | \(\{-1,0,1,2\}\) |
| \((+,-)\) | \(\{-1,0,1\}\) | \(\{-1,0,1,2\}\) |
| \((-,+)\) | \(\{-2,-1,0,1\}\) | \(\{-1,0,1\}\) |

In every case \(|\arg z|\le2\pi\). This bound is used in section 5.

**The two-loop weight is regular at \(z=-1\) (Proposition 1).** Put
\(e=1+z\) with \(\operatorname{Im}e>0\) and \(\lambda=\log e\). The inversion
relations

\[
\mathrm{Li}_2(y)+\mathrm{Li}_2(1/y)=-\tfrac{\pi^2}{6}-\tfrac12\log^2(-y),\quad
\mathrm{Li}_3(y)-\mathrm{Li}_3(1/y)=-\tfrac{\pi^2}{6}\log(-y)-\tfrac16\log^3(-y),
\]
\[
\mathrm{Li}_4(y)+\mathrm{Li}_4(1/y)=-\tfrac{7\pi^4}{360}-\tfrac{\pi^2}{12}\log^2(-y)-\tfrac1{24}\log^4(-y),
\]

applied to \(y=1/(1+z)\) and \(y=z/(1+z)\), express \(\bar F^{(2)}\) through
\(a_n=\mathrm{Li}_n(e)\), \(b_n=\mathrm{Li}_n(-e/(1-e))\) and
\(c=\log(1-e)\), all analytic at \(e=0\), and through \(\lambda\), with
\(\log(-1/e)=i\pi-\lambda\), \(\log((1-e)/e)=c-\lambda\) and
\(\log z=i\pi+c\). Expanding (`f2_expansion_at_minus_one`):

- the coefficients of \(\lambda^4\) and \(\lambda^3\) vanish identically;
- the coefficient of \(\lambda^2\) is \(\tfrac23(2a_2+2b_2+c^2)\) and that of
  \(\lambda\) is \(-\tfrac23(c+i\pi)(2a_2+2b_2+c^2)\);
- \(2a_2+2b_2+c^2=2\,[\mathrm{Li}_2(e)+\mathrm{Li}_2(-e/(1-e))+\tfrac12\log^2(1-e)]=0\)
  by Landen's identity (its derivative vanishes and it is zero at \(e=0\)).

So \(\bar F^{(2)}\) equals the \(\lambda^0\) coefficient, which is analytic at
\(e=0\), with value \(11\pi^4/45\). This is exact. It extends
\(\bar F^{(2)}\) through \(z=-1\); for \(\operatorname{Im}e<0\) the extension
is the sheet with \(\operatorname{Im}L>\pi\), which the test compares with
the independent continuation formula to \(10^{-18}\). By reflection the same
holds on the lower lip. The singularities over \(z=-1\) therefore lie at
\(\operatorname{Im}L=\pm3\pi\) and beyond, and
\(f^{(2)}(L)=\bar F^{(2)}(e^L)\) is holomorphic for
\(|\operatorname{Im}L|<3\pi\).

As an independent numerical check, the continuation across \(-1<z<0\) and
the continuation across \(z<-1\) agree to \(10^{-23}\) at every point tested.

With this continuation the edge chain gives

\[
\frac{\text{jump}}{-2\pi i}=-71.43333342=-\frac{11\pi^4}{15}
\quad\text{on the segment},\qquad
\approx10^{-7}\quad\text{off it},
\]

at cusps 1 and 2. This is the value that yields
\(\Gamma^{(3)}_{\rm cusp}=44\pi^4/45\). The three-loop coefficient therefore
comes out of the Lorentzian edge chain as well as out of the circle mean.

## 3. Item 1: bounds

Assume \(f\) is holomorphic on the strip \(|\operatorname{Im}L|<2\pi+\delta\)
with \(|f|+|f'|+|f''|\le A(1+|L|)^k\) there.

**Lemma 1 (collinear limit).** Fix \(\eta>0\). For
\(\eta\le\beta\le1-\eta\) and \(r=|c_+/c_-|\le\tfrac12\),

\[
\big|J(c_-,c_+)-J(c_-,0^\pm)\big|\le C\,A\,r\,(1+\log(1/r))^{k+1},
\]

with \(C\) depending on \(\eta\), \(k\), the tilt of the contour and the
ratios of \(s,t\) only.

*Proof.* Put \(w=c_+u\). Then

\[
J=-\int_\Gamma\frac{du}{u(u+1)}\;\frac{f\big(L_0(u)-\sigma\log(1+ru)\big)}{1+ru},
\]

where \(L_0\) depends on \(c_+\) only through its sign, and \(\Gamma\) is a
ray from \(u_0=-\beta\) at a fixed angle to the real axis. On \(\Gamma\):
the distance to the poles \(u=0,-1\) is at least \(d=\eta\sin\psi\);
\(|1+ru|\ge m=\tfrac12\sin\psi\); and
\(|L_0(u)|\le C_0+2\log(1+|u|)+2\log(1/d)\). The limit \(J(c_-,0^\pm)\) is
the same integral with \(r=0\). For \(|u|\le1/r\),

\[
\Big|\frac{f(L_0-\sigma\log(1+ru))}{1+ru}-f(L_0)\Big|
\le\frac{1}{m}\Big(\sup|f'|\,|\log(1+ru)|+|f|\,r|u|\Big)
\le C_1A\,(1+\log(1+|u|))^k\,r|u| ,
\]

and for \(|u|>1/r\) each term is at most
\(C_2A(1+\log(1+|u|))^k\). Since \(|u(u+1)|^{-1}\le C_3\min(1,|u|^{-2})\),

\[
|J-J_0|\le C\,A\Big[r\int_1^{1/r}\frac{(1+\log u)^k}{u}\,du
+\int_{1/r}^\infty\frac{(1+\log u)^k}{u^2}\,du\Big]
\le C\,A\,r\,(1+\log(1/r))^{k+1}.\qquad\blacksquare
\]

Consequently, at each interior point of an edge, \(J\) has one-sided limits
as the collinear invariant tends to zero, approached at the stated rate.
For **any** regulator acting on the collinear invariant alone, the leading
singular term is the jump \(J(c_-,0^+)-J(c_-,0^-)\) times the regulator's
logarithm or pole. With the result of section 2 this jump is
\(-2\pi i\,M[F]\) on the segment and \(0\) off it. The test
`test_lemma_1_error_is_bounded_by_r_times_a_power_of_log` checks the bound at
\(\beta=0.05,0.4,0.8\) over four decades of \(r\).

The constant degrades as \(\beta\to0\): the bound becomes
\(C\,A\,r\,(1+\log(1/r)+\log(1/\beta))^{k+1}\). This is what one expects, and
section 4 shows it is not an artefact of the proof.

**Lemma 2 (exact fourth propagator).** Let \(J^{\rm ex}\) be defined with
the exact \(x_{0,i+2}^2\) of section 1, averaged over \(\theta\), on the
contour \(\rho=re^{i\psi}\), \(0\le r\le\rho_{\max}\), with \(\rho_{\max}\)
small enough that \(|\delta|\le\tfrac12\) below. Then for \(\beta\) in the
segment,

\[
\big|J^{\rm ex}-J\big|\le
\frac{C\,A}{(1-\beta)^3\sin^3\psi}\,|c_+|\,
\Big(1+\log\frac{\rho_{\max}}{|c_+|}+\log\frac1\beta\Big)^{k+1}.
\]

*Proof.* The ratio of the exact propagator to its edge value is \(1+\delta\),
\(\delta=\delta_0+\delta_1\cos\theta\), with

\[
\delta_0=\frac{D'\alpha+w}{D'(1-\beta)},\qquad
\delta_1=\frac{2\sqrt{\rho D'u/D}}{D'(1-\beta)} .
\]

It enters through \(G(\delta)=f(L+\sigma\log(1+\delta))/(1+\delta)\). By
Taylor's theorem \(G(\delta)=f(L)+G'(0)\delta+R\) with
\(G'(0)=\sigma f'(L)-f(L)\) and \(|R|\le\tfrac12\sup|G''|\,|\delta|^2\). The
angular average of \(\delta\) is \(\delta_0\), and of \(|\delta|^2\) at most
\(2|\delta_0|^2+|\delta_1|^2\). Since \(|\alpha|=|c_+/D|\),
\(|w|\le\beta|c_+|+r\) and \(|\delta_1|^2\) is linear in \(r\),

\[
\big|\langle G(\delta)\rangle-f(L)\big|\le
\frac{C_1A\,(1+|L|)^k}{(1-\beta)^2}\,\big(|c_+|+r\big).
\]

On the contour the distance from \(re^{i\psi}\) to a real point \(a\) is at
least \(\sin\psi\,\max(r,|a|)\). With \(w_0=-\beta c_+\) this gives
\(|w|\ge\sin\psi\max(r,\beta|c_+|)\),
\(|w+c_+|\ge\sin\psi\max(r,(1-\beta)|c_+|)\) and
\(|w+c_-|\ge\tfrac12\sin\psi\,|c_-|\). Hence

\[
|J^{\rm ex}-J|\le|c_-c_+|\int_0^{\rho_{\max}}
\frac{C_1A(1+|L|)^k(|c_+|+r)}{(1-\beta)^2\,|w||w+c_-||w+c_+|}\,dr
\le\frac{C_2A\,|c_+|}{(1-\beta)^3\sin^3\psi}\int_0^{\rho_{\max}}
\frac{(1+|L|)^k}{\max(r,\beta|c_+|)}\,dr ,
\]

using \((|c_+|+r)/\max(r,(1-\beta)|c_+|)\le2/(1-\beta)\). On the contour
\(|L|\le C_0+\log(1/\beta)+2|\log(r/|c_+|)|\), and the integral is bounded by
\((1+\log(\rho_{\max}/(\beta|c_+|)))^{k+1}\) up to a constant.
\(\blacksquare\)

The test `test_lemma_2_...` evaluates the two-dimensional integral with the
exact propagator. The difference falls as \(c\log(1/c)\): \(7.9\times10^{-2}\),
\(1.5\times10^{-3}\), \(2.2\times10^{-5}\) at \(c=10^{-3},10^{-5},10^{-7}\).
The jump with the exact propagator still tends to \(-2\pi i\,M[F]\).

**Together.** Lemmas 1 and 2 make the edge statement of GATE_6_CLOSURE a
statement about the exact integrand: at each interior point of each edge the
coefficient of the collinear singularity is \(-2\pi i\,M[F]\), and it is
zero on the rest of the null line.

## 4. Correction: a hard cutoff is not an admissible regulator

GATE_6_CLOSURE section 7 said that summing the four sign quadrants gives a
density \(2\pi i\,M[F]\) "in both sectors, and each sector contributes
\(\tfrac12\log^2\)". The first half is correct. The inference to a
\(\log^2\) with a **hard cutoff** on the collinear invariants is correct
only for weights that stay bounded as \(z\to0,\infty\).

Where the two collinear invariants are comparable (\(r\sim1\)) the quadrant
sum departs from its collinear limit. For the constant weight the departure
is a fixed function of \(r\) (it equals 2 at \(r=1\)), and integrates to a
single logarithm. For a growing weight it scales with the weight at the
cutoff \(\Lambda=|D|/|c|\):

| Weight | \(\big|\int\frac{dr}{r}\,[Q(r)-Q_{\rm lim}]\big|/2\pi\) at \(c=10^{-12},10^{-16},10^{-20}\) | Scaling |
| --- | --- | --- |
| \(\log z\) | \(44,\ 58,\ 56\) | \(\propto\log\Lambda\) |
| \(\log^2z+\pi^2\) | \(1.2\times10^3,\ 2.2\times10^3,\ 3.1\times10^3\) | \(\propto\log^2\Lambda\) |
| \(z^{0.3}\) | \(6.0\times10^3,\ 9.5\times10^4,\ 1.4\times10^6\) | \(\propto\Lambda^{0.3}\) |

(40-point quadrature in \(\log r\) over 32 units; the third value in each
row is the least accurate.) Integrated over the remaining variable with a hard cutoff these
would give \(\log^{k+1}\) or a power of the cutoff, exceeding \(\log^2\).

They are artefacts of cutting through the light cones. The integrand is a
distribution such as \((x-i0)^{-1-p}\) there; its integral is finite by
contour deformation, and removing a slab \(|x|<\delta\) removes a piece of
the same size as what is kept. In an analytic regulator the same terms are
pure powers \(|c|^{-p}\) with no logarithm, hence no pole unless the exponent
vanishes, and their derivatives at \(p=0\) give at most a single pole.

So:

- the **edge statement** is regulator independent and unaffected;
- the **double pole** must be defined with an analytic regulator, as in
  section 5, not with a hard cutoff on \(x_{0j}^2\).

The test `test_hard_cutoff_is_not_an_admissible_regulator_for_growing_weights`
records this behaviour.

## 5. Item 2: dimensional regularisation

**Local density with the dimensional weight.** In \(D=4-2\epsilon\) the
transverse measure becomes \(\rho^{-\epsilon}d\rho\). For \(F=z^p\) and
\(\epsilon>|\operatorname{Re}p|\) the quadrant sum, rescaled by
\(|c_-c_+/D|^{\epsilon}\), has a limit at the cusp. In the region
\(\rho\sim|c_-c_+/D|\) one has \(z=-1+\rho|D|/(c_-c_+)\), so each quadrant is
a ray integral of \(z^{p-1}|1+z|^{-\epsilon}\) on its sheet, an Euler Beta
function. The sum is

\[
Q_\epsilon(p)=-2\cos(\pi p)\,\Gamma(p)\Big[e^{-i\pi p}\frac{\Gamma(\epsilon-p)}{\Gamma(\epsilon)}
-\frac{\Gamma(1-\epsilon)}{\Gamma(1+p-\epsilon)}\Big]
+2e^{-i\pi p}\frac{\Gamma(1-\epsilon)\Gamma(\epsilon-p)}{\Gamma(1-p)} .
\]

The numerical quadrant sum with the \(\rho^{-\epsilon}\) weight agrees with
this to \(2\times10^{-4}\) at \((p,\epsilon)=(-0.2,0.5)\) and \((-0.4,0.3)\). The
closed form is meromorphic in \(\epsilon\), and

\[
\lim_{\epsilon\to0}Q_\epsilon(p)=2\pi i\,\frac{\sin\pi p}{\pi p}=2\pi i\,M[z^p],
\]

with error linear in \(\epsilon\). This is the value of the unregulated
loop integral. The other regions (\(\rho\sim|c_\pm|\)) scale as pure powers
of \(|c_\mp|\) for \(p\ne0\) and give no double pole.

This is the method of regions. It is consistent with, and explains, the
global Mellin–Barnes evaluation of
[GATE_6_FUNCTIONAL.md](GATE_6_FUNCTIONAL.md), which gives
\(L[z^p]=\sin\pi p/(\pi p)\) for the whole integral. It is not an independent
proof: the domains of convergence of the two steps (\(\epsilon>|{\rm Re}\,p|\)
for the rays, \(\epsilon<0\) for the collinear integrals) do not overlap,
and the continuation between them is the standard dimensional-regularisation
argument.

**Lemma 3 (from the Mellin family to a class of weights).** Suppose
\(f(L)=F(e^L)\) is holomorphic for \(|\operatorname{Im}L|<2\pi+\delta\) and
polynomially bounded there. Then \(F\) is a superposition
\(\int c(p)\,z^p\,dp\) on lines \(\operatorname{Re}p=\pm\sigma\),
\(0<\sigma<1\), with \(|c(p)|\le C\,e^{-(2\pi+\delta')|\operatorname{Im}p|}\),
and the superposition converges absolutely and uniformly on the Lorentzian
integration region.

*Reason.* Split \(f=f\,\chi+f\,(1-\chi)\) with
\(\chi(L)=1/(1+e^{L/4})\), whose poles lie at
\(\operatorname{Im}L=\pm4\pi\). Each piece decays exponentially in one
direction of \(\operatorname{Re}L\), so its two-sided Laplace transform
exists on a strip of \(p\). Holomorphy in a strip of half-width
\(2\pi+\delta\) gives the stated decay of \(c(p)\) by shifting the contour
in \(L\). On the Lorentzian region \(|\arg z|\le2\pi\) (section 2), so
\(|z^p|\le|z|^{\pm\sigma}e^{2\pi|\operatorname{Im}p|}\), and the integrand of
the \(p\) integral is dominated by an integrable function.
\(\blacksquare\)

The hypothesis holds for \(\bar F^{(0)}\) and \(\bar F^{(1)}\) (entire in
\(L\)) and for \(\bar F^{(2)}\) (holomorphic for
\(|\operatorname{Im}L|<3\pi\), section 2). The need for \(2\pi\), not
\(\pi\), is a consequence of the sheet bookkeeping: the Lorentzian contour
leaves the principal sheet.

**Lemma 4 (uniform bound).** Let
\(N_\eta(p)=\epsilon^2st\times\)(regulated integral of \(z^p\)) at
\(\epsilon=-\eta\). For \(|\operatorname{Re}p|\le\sigma_0<1\) and
\(0<\eta\le\eta_0\) with \(\sigma_0+\eta_0<1\),

\[
|N_\eta(\sigma+i\tau)|\le C\,\frac{e^{\pi|\tau|}}{1+|\tau|},
\]

with \(C\) independent of \(\eta\) and \(\tau\).

*Proof.* On the contour \(w=-1-p+\eta/2+iy\) of the Mellin–Barnes integral
two of the Gamma functions lose their dependence on \(p\):

\[
N_\eta(p)=\frac{\eta^2}{\Gamma(2\eta)}\int\frac{dy}{2\pi}\,
\big|\Gamma(\tfrac\eta2+iy)\big|^4\Big(\frac st\Big)^{\eta/2+iy}\,
\frac{\Gamma(1+p-\tfrac\eta2-iy)\,\Gamma(1-p-\tfrac\eta2+iy)}
{\Gamma(1+p)^2\,\Gamma(1-p)^2}.
\]

By Stirling's bound on vertical strips the numerator is at most
\(C_1(1+|\tau-y|)\,e^{-\pi|\tau-y|}\), and
\(|\Gamma(1+p)\Gamma(1-p)|^{-2}=|\sin\pi p/(\pi p)|^2\le C_2e^{2\pi|\tau|}/(1+\tau^2)\).
Since \(e^{-\pi|\tau-y|}\le e^{-\pi|\tau|}e^{\pi|y|}\), the ratio is at most
\(C_3e^{\pi|\tau|}e^{\pi|y|}(1+|y|)/(1+|\tau|)\), using
\(1+|\tau-y|\le(1+|\tau|)(1+|y|)\). The weight satisfies
\(|\Gamma(\eta/2+iy)|^4\le C_4e^{-2\pi|y|}(1+|y|)^{2+2\eta}(\eta^2/4+y^2)^{-2}\)
(it behaves as \((\eta^2/4+y^2)^{-2}\) near \(y=0\) and as
\(4\pi^2|y|^{2\eta-2}e^{-2\pi|y|}\) for large \(|y|\)), so its integral
against \(e^{\pi|y|}(1+|y|)\) is at most \(C_5/\eta^3\). With
\(\eta^2/\Gamma(2\eta)\le3\eta^3\) the factors of \(\eta\) cancel.
\(\blacksquare\)

The test evaluates \(|N_\eta(\sigma+i\tau)|(1+\tau)e^{-\pi\tau}\) for
\(\tau\) up to 8 and \(\eta=0.1,0.02,0.004\): it stays between 0.57 and 3.6.
It also checks \(N_\eta(p)\to4\sin\pi p/(\pi p)\) at a complex \(p\).

**Conclusion of item 2.** Combine Lemma 3
(\(|c(p)|\le Ce^{-(2\pi+\delta')|\tau|}\)) with Lemma 4. The function
\(c(p)N_\eta(p)\) is dominated, uniformly in \(\eta\), by
\(Ce^{-(\pi+\delta')|\tau|}\), which is integrable. By dominated convergence

\[
\lim_{\eta\to0}\int c(p)\,N_\eta(p)\,dp=\int c(p)\,\frac{4\sin\pi p}{\pi p}\,dp=4\,M[F],
\]

the last step being Fubini on the unit circle. So the double pole of the
superposition is \(4M[F]\).

**The Fubini step.** The argument also needs the regulated integral with
weight \(F\) to equal \(\int c(p)\times\)(regulated integral of
\(z^p\))\(\,dp\). This holds on an integration contour on which the
regulated integrand is absolutely integrable and \(|\arg z|\le2\pi+\delta''\)
with \(\delta''<\delta'\). Such a contour is constructed in
[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md) in all three kinematic regions; its
properties are proved there with explicit constants, and section 10 of that
document identifies the regulated integral on it with the Mellin–Barnes
form.

## 6. Pinch surfaces (the regions away from the edges)

**Proposition 5.** The only points at which the integration contour is
pinched are the four closed edge segments.

*Argument.* All four propagators, and the weight through \(L\), carry the
same \(-i0\). A pinch at a real point requires weights \(\lambda_i>0\) on the
invariants that vanish there with
\(\sum_i\lambda_i\,\partial_{x_0}x_{0i}^2=0\), that is
\(\sum_i\lambda_i(x_0-x_i)=0\): \(x_0\) is a convex combination of the
\(x_i\) involved and is null-separated from each. With
\(x_0=x_2-ap_1+(c+d)p_2+dp_3\) the four invariants are explicit polynomials
(`landau_pinch_analysis`):

- two adjacent points: both invariants vanish on the whole segment between
  them, which is the edge;
- two opposite points: on the diagonal the invariants are \(s(a-1)^2\) and
  \(sa^2\), which vanish only at its ends;
- three points: the only solution is the cusp;
- four points: no solution.

Infinity is not a singular point, because \(d^4x_0/\prod x_{0i}^2\) is
invariant under conformal inversion. The weight adds no singular surfaces:
on the sheets reached, \(f\) is singular only over \(z=0,\infty\)
(Proposition 1 removes \(z=-1\)), which are the same light cones.

Where the contour is not pinched it can be displaced into the complex
domain, the integrand is bounded there, and the regulated integral is
analytic in \(\epsilon\) at \(\epsilon=0\). So poles arise only from the
edges: a single pole from the interior of each edge, with the coefficient of
section 3, and a double pole only from the ends. \(\blacksquare\)

This is the standard Landau and Coleman–Norton reasoning. The displacement
of the contour is constructed explicitly in
[GATE_6_CONTOUR.md](GATE_6_CONTOUR.md). The restriction to *convex* combinations is
the same fact as the vanishing of the edge jump off the segment found in
section 2.

## 7. Summary of status

| Statement | Status |
| --- | --- |
| Exact form of the four invariants in \((\alpha,\beta,\rho,\theta)\) | Proved (identity), tested at every cusp and region |
| Edge statement at every cusp, edge and region | Computed to \(10^{-8}\) |
| Collinear limit with rate \(r\log^{k+1}\) (Lemma 1) | Proved for interior points of an edge |
| Exact fourth propagator (Lemma 2) | Proved; rate confirmed numerically |
| \(\bar F^{(2)}\) regular at \(z=-1\), holomorphic for \(|{\rm Im}L|<3\pi\) | Proved (inversion relations and Landen's identity) |
| Edge chain with \(\bar F^{(2)}\) gives \(-11\pi^4/15\) | Computed to \(10^{-8}\) |
| \(|\arg z|\le2\pi\) on the Lorentzian region | From the sheet rule; tested |
| Hard cutoff equals \(\log^2\times M[F]\) | **False** for growing weights (section 4) |
| Local dimensional density \(\to2\pi i\,M[z^p]\) | Closed form, checked numerically |
| \(L[z^p]=M[z^p]\) for the whole integral | Mellin–Barnes (GATE_6_FUNCTIONAL) |
| Uniform bound \(e^{\pi|{\rm Im}\,p|}\) on the regulated integral (Lemma 4) | Proved, with elementary Gamma bounds (GATE_6_CONTOUR section 7) |
| Extension to \(\bar F^{(0,1,2)}\) (Lemmas 3, 4) | Proved, given the contour |
| Only the edges pinch the contour (Proposition 5) | Algebra exact; displacement constructed |
| Global absolutely convergent contour (Proposition 6) | Constructed for every square in all three regions; proved with explicit constants; absolute integrability proved line by line |
| Regulated integral on the contour equals the Mellin–Barnes form | Proved in all three regions (GATE_6_CONTOUR sections 11 and 12) |

**What would complete it.** Review by someone who works with these
integrals.

## 8. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_bounds.py      # about 3.5 minutes
PYTHONPATH=src python -m examples.lorentzian_all_cusps
```
