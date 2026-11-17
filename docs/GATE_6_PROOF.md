# Gate 6 contour identity

**Status.** The four-point benchmark is reproduced through two loops. The
global identity between the independently constructed weighted-divisor chain
and the Lorentzian Feynman chain remains open. The scope, baseline, and commands
for the completed suite are recorded in [benchmark validation](BENCHMARK_VALIDATION.md).

**Update.** [GATE_6_CLOSURE.md](GATE_6_CLOSURE.md) compares the two chains by
a different route, in position space at a light-like edge, without the
regulated complex cap used below. It does not resolve the rim collisions of
sections 2 to 5; it avoids introducing that lift. Read its section 8 for the
level of rigor.

## 1. Scalar identity established by the benchmark

Use \(\bar F^{(0)}=-1\) and \(\bar F^{(1)}=\log^2 z+\pi^2\), with the principal
logarithm on the open unit-circle arc \(-\pi<\phi<\pi\), \(z=e^{i\phi}\).
For loop order \(L\), write \(M_L=(2\pi)^{-1}\int_{-\pi}^{\pi}
\bar F^{(L-1)}(e^{i\phi})\,d\phi\). Then \(M_1=-1\) and \(M_2=2\pi^2/3\),
the latter from integrating \(\pi^2-\phi^2\). The implementation uses the
conventions \(\mathcal I=-g^2\) times the circle mean and
\(g\partial_g\Gamma=-8\mathcal I\). Thus \(\Gamma^{(L)}=-4M_L/L\).
[Dixon et al., arXiv:2605.28926, Eqs. (2.30)–(2.31) and (6.2)–(6.5)](https://arxiv.org/abs/2605.28926)

The weighted-divisor construction in [example 36](../examples/weighted_boundary_cycle.py)
uses radial-last orientation and a quotient fiber whose lift traverses half
the phase circle before the order-two identification. Its oriented factor is
\(+1\). The CD period with the circle integrated first is independently
\(4\int_0^1\log u/(u-1)\,du=2\pi^2/3\), by expansion of \((1-u)^{-1}\)
and termwise integration. The bound on the reduced seam excision justifies
that calculation. These are scalar results in the stated contour and
normalization. The original finite CD numerator and loop normalization come
from the combined Lorentzian integral.
[Arkani-Hamed, Henn and Trnka, arXiv:2112.06956, Eqs. (2.30)–(2.34)](https://arxiv.org/abs/2112.06956)

## 2. Transverse obstruction to integrating the meridians separately

The diagonal band of [example 60](../examples/cd_diagonal_meridian_band.py) is
the coefficient of \(dt\) at \(t=y/x=1\). The following local calculation
extends the analysis off the diagonal. It uses the scaled density from
[example 58](../examples/cd_remote_meridian_period.py) and the same cap radius
\(R=0.06=3/50\). Set

\[
Q=e^{i\pi/4},\quad H=1-(1+Q)s,\quad x=\epsilon A\lambda,\quad
K=k/R=5/2+\sqrt\lambda\,u,\quad J=-1+(1-t)K.
\]

Here \(0\le s\le1\), \(0<A\le1-s\), and the base cap has \(|k|\le1\).
At \(\epsilon=0\) **after scaling**, the regulated polar equation in \(K\) is

\[
N(K)=H[-1+(1-t)K]+\frac{i}{\lambda}(K-5/2)^2=0.
\]

Multiplying by \(\lambda/i\) makes the quadratic monic. With \(b=1-t\), its
discriminant is \(\Delta=-\lambda H[\lambda H b^2+i(4-10b)]\) and its root
center is \(K_c=5/2+i\lambda Hb/2\). Set \(s_*=1/(1+\Re Q)=2-\sqrt2\) and
\(a=\Im(Q)s_*=\sqrt2-1\), so \(H(s_*)=-ia\). At a double root, with
\(b\ne0\),

\[
\lambda a b^2=4-10b,\qquad K_c=\frac{2}{b}-\frac52.
\]

Taking \(K_c=\pm1/R\) puts the double root on the rim. The resulting values
are \(t=103/115\approx0.895652174\) and \(t=97/85\approx1.141176471\), both
within the transverse interval used in example 57. In both cases
\(\lambda=(R^{-2}-25/4)/a\approx655.526044\). In contrast, the simultaneous
diagonal rim crossings in example 52 occur at distinct roots. At either
transverse double root, \(J=1-5b/2\ne0\), and \(K_c\) differs from the base
and frame punctures \(0,1,5/2\).

### Persistence and absolute integrability

Fix an interior \(A\), for example \(A=1/4\). Near either rim collision the
three real equations (real and imaginary parts of the discriminant and
\(R^2|K_c|^2-1=0\)) have independent derivatives in \((s,\lambda,t)\).
Indeed, at the displayed solution the real part of
\(\lambda H b^2+i(4-10b)\) has nonzero \(s\) derivative and zero
\(\lambda,t\) derivatives. At \(s=s_*\) its imaginary equation has nonzero
\(\lambda\) derivative \(-ab^2\). Along that discriminant-zero curve the
rim equation has derivative \(2R^2K_c\,dK_c/dt\ne0\). The quadratic in
[example 58](../examples/cd_remote_meridian_period.py) is smooth near these
points at small finite \(\epsilon\) with nonzero leading coefficient. The
implicit function theorem therefore gives nearby rim collisions at small
positive regulator and under small changes of \(A\).

Write the local slice density as \(f(K)/[(K-K_+)^3(K-K_-)^3]\) with analytic
\(f\) and \(f(K_c)\ne0\). This follows from the density in example 58 after
the linear change to \(K\), since \(J\ne0\) and the base factors are nonzero.
With \(d=K_+-K_-\), the residue at \(K_+\) is
\(6f(K_+)/d^5-3f'(K_+)/d^4+f''(K_+)/(2d^3)\). The residue at \(K_-\) has
the opposite leading term. Therefore a single residue has magnitude bounded
below by a constant times \(|\Delta|^{-5/2}\) near the collision. The common
factor \(\epsilon A\) from \(dx\) is nonzero at fixed positive regulator.

Use local real coordinates \((\Re\Delta,\Im\Delta,\rho)\), where
\(\rho=R^2|K_c|^2-1\). For an angular sector of \(\Delta\) where the split
has nonzero radial component, \(|\rho|<c|\Delta|^{1/2}\) with a sufficiently
small \(c>0\) puts one root on each side of the rim. The Jacobian in these
coordinates is nonsingular. Writing \(r=|\Delta|\), a lower bound for the
absolute residue integral over this region contains
\(\int_0^\eta r^{-5/2}r\,r^{1/2}\,dr=\int_0^\eta dr/r\).
The ordinary absolute integral of the isolated one-root contribution thus
fails near the rim collision. The coefficient of the full chain, including
the rim and transported contour, is undecided by this estimate. A grouped
contour around the two roots has a smooth period locally, since the enclosing
loop can be held off all poles even at coincidence. Its boundary relation to
the physical contour must still be specified.

## 3. Local rim detours through a double root

Let \(p=(\lambda,t,A,s)\) and let \(p_*\) denote an interior auxiliary slice of
either collision in section 2. Work in the \(K\) plane and give the disk
\(D_R=\{|K|\le R^{-1}\}\) its complex orientation. Choose a closed disk \(E\)
centered at the collision point \(K_*\), small enough to miss \(0,1,5/2\) and
to meet the rim transversely at two points. On a compact neighborhood of
\(p_*\) and for all sufficiently small \(\epsilon\ge0\) in the scaled variables,
both regulated roots stay strictly inside \(E\). This follows by continuity of
the monic quadratic coefficients, or by applying Rouché's theorem on
\(\partial E\) at the double root. The boundary of \(E\) can stay fixed.

The boundary of \(D_R\cup E\) detours outside the original rim and includes
both roots. The boundary of \(D_R\setminus\operatorname{int}E\) detours inside
the rim and excludes both roots. Keep every original base and frame hole
in the two chains. The outer detour uses the arc of \(\partial E\) outside
\(D_R\) counterclockwise. The inner detour uses the arc inside \(D_R\)
clockwise. Write the resulting oriented rim cycles as \(C_{\mathrm{out}}\)
and \(C_{\mathrm{in}}\). Their difference, as one-chains, is
\(C_{\mathrm{out}}-C_{\mathrm{in}}=\partial E\) counterclockwise.

At a parameter point where the original rim avoids the roots, let \(S\) be
the set of roots in \(D_R\) and \(S^c\) the remaining roots in \(E\). Write
\(M_J\) for the sum of counterclockwise meridians around the roots in \(J\),
with \(M_\varnothing=0\). The punctured regions \(E\cap D_R\) and
\(E\setminus D_R\) supply the respective homologies

\[
C_{\mathrm{rim}}-C_{\mathrm{in}}\sim M_S,\qquad
C_{\mathrm{out}}-C_{\mathrm{rim}}\sim M_{S^c}.
\]

Indeed the positive rim arc occurs in the first region and the negative
rim arc in the second. After subtracting the included meridians, each
difference bounds its region in the complement. These formulas cover each
single-root side where the rim contains one root, and give the periods with
the orientation fixed. For instance, transporting the original rim to the
outer detour adds the meridian of the excluded root. Summing the two
relations yields the pair contour. Each detour extends through the
collision with no pole on its path. The individual meridian terms from the
original rim can be singular even though their grouped sum extends.

### The grouped period and its regulator bound

Write the sliced form in the orientation of example 58 as

\[
\omega_\epsilon=\frac{4QJ^2}{(1+x)(1+tx)K(1-K)N_\epsilon(K)^3}
\frac{dK}{2\pi i},
\qquad N_\epsilon(K)=HJ+i(u+hP_0)^2,
\quad u=\frac{K-5/2}{\sqrt\lambda},\quad h=\epsilon\sqrt\lambda.
\]

Here \(P_0\) is the affine function of \(K\) from example 58, with no
singularity on \(E\). On the compact product of the parameter patch and
\(\partial E\), \(N_\epsilon\) is uniformly bounded away from zero after
shrinking the patch and regulator interval. The base factors are also
bounded away from zero. Consequently the grouped period
\(G_\epsilon(p)=\oint_{\partial E}\omega_\epsilon\) is bounded and smooth in
the real parameters and regulator, even through the double root. For a
fixed compact patch \(V\) of the \(p\) variables, the difference of the two
detour periods is therefore

\[
\epsilon\int_V A\,G_0(p)\,d\lambda\,dt\,dA\,ds+O(\epsilon^2).
\]

This is from \(dx=\epsilon A\,d\lambda\) modulo \(dA\) at fixed regulator.
The error follows from the uniform bound on the first regulator derivative
on the fixed loop. The estimate applies to the local detour difference.

At \(p_*\) and \(\epsilon=0\) in the scaled variables,
\(N_0(K)=i(K-K_*)^2/\lambda_*\). The rational form has no residue at
infinity, and its only poles outside \(E\) are the base poles at \(K=0,1\).
Taking the negative sum of those residues gives

\[
G_0(p_*)=4iQ\lambda_*^3\left(\frac{t_*^2}{(1-K_*)^6}
-\frac{1}{K_*^6}\right).
\]

At \(K_*=+50/3\) this is approximately \(-6.052350324174+6.052350324174i\).
At \(K_*=-50/3\) it is approximately \(3.045968993180-3.045968993180i\).
Direct circle quadrature of the existing slice density agrees with the
remote-residue formula to within \(3\times10^{-14}\) at the two limiting
collisions and at nearby continued collisions for
\(A=1/4\), \(\epsilon=10^{-7},10^{-6}\). Thus the local two-detour ambiguity
has a bounded, generally nonzero coefficient at order \(\epsilon\), and
its period vanishes as \(\epsilon\) tends to zero on this compact patch.

The two detours differ topologically in the regulated divisor complement
at positive regulator. The original rim and the meridians specify the
local relation above. The Lorentzian \(i0\) prescription still has to
select the detour and its attachments on the parameter boundary. A change
of detour over a parameter subregion carries an interface chain with the
pair-meridian boundary, so the local period bound alone does not identify
the full relative cycle.

## 4. Physical prescription with fixed affine Feynman weights

In the finite frame of [the regulated parameter integral](../examples/finite_cross_ratio.py),
set the fifth homogeneous weight to one. Hold the four remaining weights
\(a_1,\dots,a_4\) fixed and define \(\ell=K-5/2\), \(P=1+\sum_j a_j\). This
uses a different slicing from sections 2 and 3, which keep \((\lambda,t,A,s)\)
fixed while varying \(K\). In the fixed-weight slicing, the coefficients of
the terms linear in \(a_j\) are

\[
(p_1,p_2,p_3,p_4)=\left(1+\frac{3}{2\ell},-1-\frac{5}{2\ell},
-\frac{2}{\ell},\frac{2}{\ell}\right).
\]

Thus the affine Feynman denominator with the original regulator is exactly

\[
F_\epsilon(\ell)=C+\frac{L}{\ell}-i\epsilon P^2,
\quad C=2(a_1a_3-a_2a_4)+a_1-a_2,\quad
L=\tfrac32a_1-\tfrac52a_2-2a_3+2a_4.
\]

For positive real weights, \(P>0\). If \(L\ne0\), the unique finite base zero
at \(\epsilon>0\) is simple and lies at

\[
\ell_\epsilon=-\frac{L}{C-i\epsilon P^2},\qquad
\Im\ell_\epsilon=-\frac{\epsilon LP^2}{C^2+\epsilon^2P^4}.
\]

At a real zero \(\ell_0=-L/C\ne0\), the displacement to first order is
\(i\epsilon P^2/F_0'(\ell_0)\), with \(F_0'(\ell_0)=-L/\ell_0^2\).
The oriented real path through that point consequently passes above the
pole for \(L>0\) and below it for \(L<0\). For example, fixed weights
\((1,1/2,1/4,1/4)\) give \(K_0=13/6\) and
\(\lim_{\epsilon\downarrow0}\Im K_\epsilon/\epsilon=-4\). Fixed weights
\((1/4,1/2,1/4,1/2)\) give \(K_0=19/10\) and the slope \(+6\).
Both zeros are on the real frame side \(\ell<0\). When \(L=0\), the
denominator is independent of \(K\) off the frame and has no finite zero
on this positive-weight slice at \(\epsilon>0\). When \(C=0\) and \(L\ne0\),
the zero tends to infinity as the regulator is removed.

### Relation to the rotated ray and the two rim detours

For real \(\ell<0\), keep \(a_1,a_2,a_3\) positive and rotate
\(a_4=r e^{i\theta}\) through \(0\le\theta\le\pi/4\) at fixed \(r>0\).
Write \(P=P_*+c+id\), where \(P_*=1+a_1+a_2+a_3>0\),
\(c=r\cos\theta\), \(d=r\sin\theta\), and \(c\ge d\ge0\).
The imaginary part of the regulated denominator along the rotation is

\[
\Im F_\epsilon=(2/\ell-2a_2)d-\epsilon[(P_*+c)^2-d^2]<0.
\]

This gives a pointwise pole-free connection to the rotated ray on the real
frame side. It transports the local sign of the Lorentzian prescription.
The outer base poles and the frame puncture remain excluded from this
local comparison.

The two collisions in section 2 are in the moving-weight slicing. In the
resolved chart \(w=A\ell/x\), the affine weights obey \(a_1=-x/\ell\)
and \(a_3=BJ/[2(1+x)w]\), \(a_4/Q=sJ/[2(1+tx)w]\), with
\(B=1-A-s\). At the positive limiting rim point, \(\ell>0\) forces
\(a_1<0\) along the real cap with \(x>0\). At the negative limiting rim
point, \(w<0\) and \(J=1-5(1-t_*)/2>0\) force \(a_3<0\) and
\(a_4/Q<0\) in the auxiliary interior. The strict signs persist on
nearby real cap points. The positive rotated ray does not directly meet
those real cap neighborhoods. Accordingly the fixed-weight sign rule
gives the physical local bypass on its own slice, while assigning
\(C_{\mathrm{in}}\) or \(C_{\mathrm{out}}\) requires a comparison that
also transports the weights and the relative faces. A fixed infinitesimal
imaginary shift of the base coordinate contributes \(\Im F_0\) with sign
set by \(F_0'\), which changes between the two positive-weight examples.
It cannot reproduce the sign of \(-i\epsilon P^2\) on both examples.

## 5. Local transport through the positive-weight vertex

There is a local comparison that moves the weights at fixed base
coordinate. Start with a compact slice of the positive rotated ray in
section 4 and scale its four affine weights to zero, keeping \(a_5=1\).
On the real frame side \(\ell<0\), the imaginary-part inequality in section 4
continues to hold, including at the vertex where \(F_\epsilon=-i\epsilon\).
At the vertex the base path can move in the punctured base chart while
the denominator remains \(-i\epsilon\). This moves a local path from the
real frame side to the rim patch, with the base and frame punctures avoided.

In that rim patch, let \(a_j^{\mathrm{cap}}(K)\) denote the affine weights of
the moving-weight chart, with \(x=\epsilon A\lambda\). Take the explicit
interpolation \(a_j=q a_j^{\mathrm{cap}}(K)\), \(0\le q\le1\). The cap
weights are \(O(\epsilon)\) on the patch. More explicitly they are

\[
(a_1,a_2,a_3,a_4)^{\mathrm{cap}}
=\frac{\epsilon\lambda}{\ell}
\left(-A,-At,\frac{BJ}{2(1+x)},\frac{QsJ}{2(1+tx)}\right),
\qquad B=1-A-s.
\]

The exact cap identity \(F_0=C_0 H/A\), with
\(C_0=p_1a_1+p_2a_2=-xJ/\ell^2\), gives the linear contribution in these
weights up to \(O(\epsilon^2)\). Thus along this interpolation

\[
\frac{F_\epsilon(q a^{\mathrm{cap}})}{\epsilon}
=-i-q\frac{\lambda HJ}{\ell^2}+O(\epsilon).
\]

The estimate is uniform on any compact patch separated from the frame
and on \(0\le q\le1\). At either limiting rim collision \(K_*=\pm R^{-1}\),
use \(H=-ia\), \(a=\sqrt2-1\), and
\(\lambda aJ=\ell^2-(K-K_*)^2\). The scaled denominator becomes

\[
-i\left[(1-q)+q\frac{(K-K_*)^2}{\ell^2}\right].
\]

For real \(K\) near \(K_*\), its imaginary part is strictly negative on
each closed interval \(0\le q\le q_0<1\). The sign survives at sufficiently
small positive regulator and sufficiently small perturbations of the
collision parameters. At \(q=1\), use the fixed local detour boundary of
section 3. This is a deformation of a local path. Scaling the weight
coordinates to the vertex also collapses the corresponding integration
directions, so the homotopy's lateral faces must be retained when making
a comparison of full chains.

### Which roots lie inside during this transport

Let \(c=5/2\). After clearing \(\ell^2\), the limiting polynomial along the
interpolation is proportional to

\[
(1-q)(K-c)^2+q(K-K_*)^2.
\]

For \(0<q<1\) its two roots lie on opposite sides of the real axis at

\[
K=(1-q)c+qK_*\ \pm\ i\sqrt{q(1-q)}(K_*-c),
\qquad |K|^2=(1-q)c^2+qR^{-2}<R^{-2}.
\]

Both roots approach the rim from its interior in this scaled transport.
Continuing the rim locally while keeping both enclosed leads to the outer
detour \(C_{\mathrm{out}}\) from section 3. The calculation identifies that
detour for this vertex-based path. It does not assert that every transport
of the physical cycle has the same relative boundary. At \(q=0\), the
cleared polynomial has a double zero at \(K=c\), but the actual denominator
is constant. Tracking the pair all the way from the vertex must therefore
account for the frame puncture, where this normalization changes the
winding. At positive regulator, continuity also preserves strict interior
root placement on compact subintervals of \(0<q<1\) for this family. The
endpoints require the frame and rim charts. Exact polynomial identities
and finite-regulator evaluations on both patches check these statements.

## 6. Conditions for closing Gate 6

1. **Construct and control the regulated comparison through the rim
   collisions.** Fix a common prescription and cutoff chain for the capped
   lift and the Lorentzian integral. Section 3 gives both local rim detours,
   their signed meridian relations, and the bound on their period difference.
   Section 4 fixes the local prescription on the original weight slices and
   connects it to the rotated ray. Section 5 supplies a local vertex-based
   transport and its rim detour. Determine the frame and lateral boundary
   terms and whether the transport extends to the physical cycle. Extend
   control to the remaining transverse domains and their endpoints, together
   with the auxiliary and base/frame faces. Specify the order of regulator
   removal and other cutoff limits and justify any interchange. Example 60
   gives a diagonal check for this construction.
2. **Establish the oriented global identity or exhibit its obstruction.**
   Put the geometric and Feynman cycles in a common resolved complement of
   the regulated polar divisor with explicit relative boundary data. Give a
   comparison chain whose oriented boundary is their difference modulo the
   agreed relative faces, and control the periods of those faces. If the
   identity fails, identify the surviving relative class and compute its
   contribution in the normalization of section 1. A matching scalar number
   in the absence of this boundary calculation does not close the gate.
