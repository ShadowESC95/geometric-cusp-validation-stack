# Mathematical model implemented by the stack

## 1. User-defined geometric pushforward

The user's proposed object is retained as a formal template:

\[
\Omega^{(k,n)}
\equiv
\int
\frac{d^{k\times n}C}{\mathrm{Vol}[GL(k)]}
\,\delta^{(4k)}(Y-CZ)\,I(C).
\]

The code separates the formal expression from the numerical realization. A positive chart for \(C\) is sampled and \(Y=CZ\) is computed. For k=1, the GL(1) quotient is fixed by \(\sum_i C_i=1\). For k=2, a positive polar-column representation or a gauge-fixed \(G_+(2,5)\) chart is used.

## 2. Four-point one-loop chart

Use

\[
AB=
\begin{pmatrix}
1 & x & 0 & -y\\
0 & z_0 & 1 & w
\end{pmatrix},
\qquad x,y,z_0,w>0.
\]

For \(Z_i=e_i\),

\[
\langle AB12\rangle=y,
\quad
\langle AB23\rangle=w,
\quad
\langle AB34\rangle=z_0,
\quad
\langle AB14\rangle=x,
\]

\[
\langle AB13\rangle=-(wx+yz_0),
\qquad
\langle AB24\rangle=-1.
\]

The canonical form is

\[
\Omega_4^{(1)}
=
\frac{dx}{x}\wedge\frac{dy}{y}\wedge\frac{dz_0}{z_0}\wedge\frac{dw}{w}.
\]

## 3. Cusp/collinear face and blow-up

For the cusp associated with the adjacent propagators 12 and 23, use

\[
y=0,
\qquad
w=0.
\]

Blow up the intersection by

\[
y=\rho\eta,
\qquad
w=\rho.
\]

Then

\[
\frac{dy}{y}=\frac{d\rho}{\rho}+\frac{d\eta}{\eta},
\qquad
\frac{dw}{w}=\frac{d\rho}{\rho},
\]

and the canonical form has a simple \(d\rho/\rho\) pole. The exceptional divisor carries the angular/projective normal variable \(\eta=y/w\).

The four-point cross-ratio is

\[
\chi=
\frac{\langle AB12\rangle\langle AB34\rangle}
{\langle AB23\rangle\langle AB14\rangle}
=
\frac{\eta z_0}{x}
\]

on the blow-up chart.

Consequently

\[
\boxed{
\frac{d\chi}{\chi}=\frac{d\eta}{\eta}
}
\]

when \(x,z_0\) are held fixed.

This is the precise mechanism tested by the linking-contour experiment: the standard \(\chi\)-contour is the image of a nontrivial cycle in the complexified projective normal direction.

## 4. Strict cusp line

The geometric collinear line may be written

\[
Z_A=Z_i,
\qquad
Z_B=Z_{i-1}+\alpha Z_{i+1}.
\]

A separate Grassmannian chart is used because the line can lie outside a convenient affine patch of the positive one-loop coordinates. This distinction is essential: \(\alpha\) is the tangential modulus of the strict cusp line, while \(\chi\) can become a coordinate on the projectivized normal/approach direction.

## 5. Two-loop logarithm integrand

The exact four-point two-loop logarithm integrand used by the research code is

\[
\widetilde I_4^{(2)}
=
-
\frac{\langle1234\rangle^3
\left(
\langle AB13\rangle\langle CD24\rangle
+
\langle AB24\rangle\langle CD13\rangle
\right)}
{D_{AB}D_{CD}\langle ABCD\rangle}.
\]

For the source chart,

\[
D_i=x_i y_i z_i w_i,
\]

\[
\langle ABCD\rangle
=
-(x_1-x_2)(w_1-w_2)-(y_1-y_2)(z_1-z_2),
\]

and

\[
n_{12}
=
 w_1x_1+w_2x_2+y_1z_1+y_2z_2.
\]

## 6. Collective radial scaling

Set

\[
y_i=\rho a_i,
\qquad
w_i=\rho b_i,
\qquad i=1,\ldots,L.
\]

For a tree graph with \(L\) vertices and \(L-1\) edges,

\[
D_{\rm ext}\sim\rho^{2L},
\qquad
D_{\rm mutual}\sim\rho^{L-1},
\qquad
N_{\rm tree}\sim\rho^{L-1}.
\]

Therefore the rational scalar integrand scales as

\[
I_{\rm tree}^{(L)}\sim\rho^{-2L}.
\]

The Jacobian from the \(2L\) small variables is

\[
d^{2L}q
\sim
\rho^{2L-1}d\rho\wedge d\Sigma,
\]

and hence

\[
\boxed{
\Omega_{\rm tree}^{(L)}
\sim
\frac{d\rho}{\rho}\wedge\omega_{\rm exc}^{(L)}.
}
\]

This power-counting result explains the mild logarithmic radial divergence of the logarithm integrand in the collective collinear limit.

## 7. Contour extraction benchmark

The four-point finite object is written

\[
F(g,z)=\sum_{L\ge0}g^{2L}F^{(L)}(z),
\]

with

\[
F^{(0)}(z)=-1,
\qquad
F^{(1)}(z)=\log^2z+\pi^2.
\]

The source contour prescription is implemented as a convention-controlled functional. The code deliberately does not merge the source's related \(F\)-normalizations into a single symbol because the Wilson-loop insertion object and the resummed geometric function differ by a coupling/sign factor.
