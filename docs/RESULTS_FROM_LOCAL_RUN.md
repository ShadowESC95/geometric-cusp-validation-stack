# Results from the build machine (18 Sep 2026)

This file records only computations actually run in the current build environment. Optional external packages such as FEniCSx, Firedrake, Meep, QuTiP, Kwant, VASP, FORM, OpenFOAM, COMSOL, deal.II, SageMath, Arb/python-flint and pySecDec were not installed in this environment; their adapters were syntax-checked only where possible and are intended for the user's own scientific-computing environment.

## Core test suite

`5 passed`

## Amplituhedron pushforward

User-defined k=2,n=6 pushforward: 8000 samples, positive-minor minimum \(1.18896\times10^{-5}\), mean margin 0.0698003.

Standard A-style k=2,n=6 pushforward: 8000 samples, positive-minor minimum \(2.82289\times10^{-6}\), 15 Plücker coordinates.

Figures:

- `figures/amplituhedron_user_k1_affine.png`
- `figures/amplituhedron_user_positivity.png`
- `figures/amplituhedron_user_plucker_pca.png`
- `figures/amplituhedron_standard_plucker_pca.png`

## One-loop exact algebra

Exceptional blow-up residue coefficient:

\[
\operatorname{Res}_{\rho=0}\Omega_4^{(1)}
=
\frac{1}{x\eta z_0}
\quad\text{(in the chosen orientation convention)}.
\]

Exceptional cross ratio:

\[
z_{\rm cross}=\frac{\eta z_0}{x}.
\]

Strict-cusp radial residue in the separate \(U_{12}\) chart:

\[
\operatorname{Res}_{\rho=0}\Omega_{\rm strict}=0.
\]

This distinction is important: the nonzero one-loop collinear logarithmic residue belongs to the blown-up adjacent-propagator intersection, while the strict one-dimensional cusp line is not itself a simple divisor of that affine canonical-form chart.

## Linking test

The numerical complex-contour image has winding approximately 1; the exact symbolic relation is

\[
\frac{dz_{\rm cross}}{z_{\rm cross}}
=
\frac{d\eta}{\eta}
\]

at fixed \(x,z_0\).

The code therefore supports a precise homological interpretation of the standard \(z\)-contour as the image of a loop in the projectivized normal direction, while it does not identify the tangential cusp modulus \(\alpha\) with \(z\).

## Two-loop collective scaling

Exact scalar leading coefficient after \(y_i=\rho a_i, w_i=\rho b_i\):

\[
\rho^4 I_4^{(2)}\to
\frac{a z_1+b x_1+c z_2+d x_2}
{a b c d x_1x_2z_1z_2
[a z_1-a z_2+b x_1-b x_2-c z_1+c z_2-d x_1+d x_2]}.
\]

For a tree graph with \(L\) vertices,

\[
I_{\rm scalar}\sim\rho^{-2L},
\qquad
J_{\rm radial}\sim\rho^{2L-1},
\qquad
\Omega\sim d\rho/\rho.
\]

This scaling was verified algebraically for \(L=1,2,3,4,5\).

## Contour benchmark

The numerical circle average of

\[
F^{(0)}(z)=-1
\]

is exactly \(-1\) up to floating-point representation, while

\[
F^{(1)}(z)=\log^2z+\pi^2
\]

on \(|z|=1\) gives

\[
\frac{1}{2\pi i}\oint\frac{dz}{z}F^{(1)}(z)
=
\frac{2\pi^2}{3}
\approx 6.5797361701131845.
\]

## JAX scaling check

The compiled JAX evaluation of the exact two-loop scalar integrand gives the fitted scaling

\[
I_4^{(2)}\propto\rho^{-4}
\]

with fitted exponent

\[
-3.9999999937.
\]

## Matrix substrate

A small Euclideanized Langevin run at \(N=8,D=4,g=0.7\) produced mean commutator density approximately 0.9816 and final value approximately 0.8432. This is a numerical probe of the matrix action, not a derivation of the continuum gauge dictionary.

## Continuum entropy proxy

The reduced 1D source-driven entropy transport sandbox produced maximum scalar proxy \(S\approx5.10\) and a final integrated value \(\approx36.33\). This is deliberately labelled a proxy and does not claim to establish the user's entropy map.
