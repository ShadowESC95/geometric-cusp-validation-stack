# Gate 6 contour identity

**Status.** The four-point benchmark is reproduced through two loops. The
global identity between the independently constructed weighted-divisor chain
and the Lorentzian Feynman chain remains open. The scope, baseline, and commands
for the completed suite are recorded in [benchmark validation](BENCHMARK_VALIDATION.md).

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

## 4. Conditions for closing Gate 6

1. **Construct and control the regulated comparison through the rim
   collisions.** Fix a common prescription and cutoff chain for the capped
   lift and the Lorentzian integral. Section 3 gives both local rim detours,
   their signed meridian relations, and the bound on their period difference.
   Transport the physical cycle into these patches, choose the detours from
   its prescription, and compute or bound the interface terms. Extend that
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
