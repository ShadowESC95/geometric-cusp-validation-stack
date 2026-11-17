# Gate 6: bounds, regulators and the remaining cases

This document supplies the three items left open by
[GATE_6_CLOSURE.md](GATE_6_CLOSURE.md), section 8:

1. error bounds for the local statement;
2. the relation between the local statement and the double pole in
   dimensional regularisation;
3. the cusps, kinematic regions and weight that were covered there by
   symmetry only.

It also **corrects** one statement of that document (section 4 below).

**Status.** Item 3 is complete: it is computed, not inferred. Item 1 gives
one lemma with a full proof, one with a proof sketch and numerical
confirmation of the predicted rate, and identifies a limitation of hard
cutoffs. Item 2 is an argument resting on the published Mellin–Barnes
evaluation plus a lemma whose hypothesis is verified here for the weights of
interest; the interchange it needs is stated, not proved. Nothing here has
been reviewed by a specialist.

Code: `src/geocusp/lorentzian_edges.py`. Tests: `tests/test_gate6_bounds.py`
(64 tests). Example 62: `examples/lorentzian_all_cusps.py`.

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

**The two-loop weight.** \(\bar F^{(2)}\) is given on the principal sheet.
Its continuation across the lower lip can be made across \(-1<z<0\), where
\(\mathrm{Li}_n(1/(1+z))\) crosses its cut, or across \(z<-1\), where
\(\mathrm{Li}_n(z/(1+z))\) does. The two results agree to \(10^{-23}\) at
every point tested. Therefore \(\bar F^{(2)}\) has no monodromy around
\(z=-1\) on the physical sheet: it is regular there, with
\(\bar F^{(2)}(-1)=11\pi^4/45\). Its singularities over \(z=-1\) lie on other
sheets, at \(\operatorname{Im}L=\pm3\pi\) and beyond, so
\(f^{(2)}(L)=\bar F^{(2)}(e^L)\) is holomorphic for
\(|\operatorname{Im}L|<3\pi\).

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
the exact \(x_{0,i+2}^2\) of section 1, averaged over \(\theta\), and
\(\rho\le\rho_{\max}\). Then for fixed \(\beta\) in the segment,

\[
\big|J^{\rm ex}-J\big|\le C\,A\,|c_+|\,\big(1+\log(\rho_{\max}/|c_+|)\big).
\]

*Sketch.* The ratio of the exact propagator to its edge value is
\(1+\delta\) with

\[
\delta=\frac{D'\alpha+w+2\sqrt{\rho D'u/D}\cos\theta}{D'(1-\beta)} .
\]

It enters through \(f(L+\sigma\log(1+\delta))/(1+\delta)\). Expanding to
second order, the term linear in \(\sqrt\rho\cos\theta\) averages to zero.
What remains is \(O(\rho+|c_+|)\) times \(|f|+|f'|+|f''|\). The measure
\(c_-c_+\,d\rho/|w(w+c_-)(w+c_+)|\) behaves as \(|c_+|\,d\rho/\rho^2\) for
\(\rho\gg|c_+|\), so the correction integrates to
\(|c_+|\log(\rho_{\max}/|c_+|)\). \(\blacksquare\)

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

**What is still assumed.** To conclude \(L[F]=\int c(p)L[z^p]\,dp=M[F]\) in
dimensional regularisation one must exchange the \(p\) integral with the
regulated loop integral and with the extraction of the \(1/\epsilon^2\)
coefficient. Lemma 3 gives domination of the integrand on the integration
region. A bound on the regulated integral of \(z^p\) that is uniform in
\(\operatorname{Im}p\) is also required; the closed Mellin–Barnes form
suggests growth no faster than \(e^{\pi|\operatorname{Im}p|}\), which is
enough, but this is not proved here.

## 6. Summary of status

| Statement | Status |
| --- | --- |
| Exact form of the four invariants in \((\alpha,\beta,\rho,\theta)\) | Proved (identity), tested at every cusp and region |
| Edge statement at every cusp, edge and region | Computed to \(10^{-8}\) |
| Collinear limit with rate \(r\log^{k+1}\) (Lemma 1) | Proved for interior points of an edge |
| Exact fourth propagator (Lemma 2) | Proof sketch; rate confirmed numerically |
| \(\bar F^{(2)}\) regular at \(z=-1\), holomorphic for \(|{\rm Im}L|<3\pi\) | Numerical (two continuations agree to \(10^{-23}\)) |
| Edge chain with \(\bar F^{(2)}\) gives \(-11\pi^4/15\) | Computed to \(10^{-8}\) |
| \(|\arg z|\le2\pi\) on the Lorentzian region | From the sheet rule; tested |
| Hard cutoff equals \(\log^2\times M[F]\) | **False** for growing weights (section 4) |
| Local dimensional density \(\to2\pi i\,M[z^p]\) | Closed form, checked numerically |
| \(L[z^p]=M[z^p]\) for the whole integral | Mellin–Barnes (GATE_6_FUNCTIONAL) |
| Extension to \(\bar F^{(0,1,2)}\) (Lemma 3) | Hypothesis verified; one interchange assumed |
| No double pole from regions away from the edges | Not proved; consistent with the Mellin–Barnes total |

**What would complete it.** (a) A written proof of Lemma 2 with constants.
(b) The uniform bound in \(\operatorname{Im}p\) needed in section 5.
(c) An analytic, not numerical, proof that \(\bar F^{(2)}\) is regular at
\(z=-1\), for example from its symbol or from the inversion relations of the
polylogarithms. (d) Review by someone who works with these integrals.

## 7. Reproduce

```bash
PYTHONPATH=src pytest -q tests/test_gate6_bounds.py      # about 3 minutes
PYTHONPATH=src python -m examples.lorentzian_all_cusps
```
