# Research plan: adversarial validation stack

## Experiment A — user-defined amplituhedron pushforward

Test

\[
\Omega^{(k,n)}=
\int\frac{d^{k\times n}C}{\mathrm{Vol}[GL(k)]}
\delta^{(4k)}(Y-CZ)I(C).
\]

Compare:

1. direct positive-Grassmannian point-cloud pushforward,
2. Gaussian regularized delta/KDE,
3. standard \(A_{n,k,4}\)-style pushforward.

Primary plots: affine convex hull for k=1, Plücker-PCA cloud for k=2, positivity margins, and kernel-density slices.

## Experiment B — one-loop boundary blow-up

Use

\[
AB=\begin{pmatrix}1&x&0&-y\\0&z_0&1&w\end{pmatrix}
\]

and

\[
\Omega_4^{(1)}=d\log x\wedge d\log y\wedge d\log z_0\wedge d\log w.
\]

Blow up \(y=w=0\) with \(y=\rho\eta,\ w=\rho\). Verify symbolically that the total form has a simple radial pole and that

\[
z_{\rm cross}=\frac{\eta z_0}{x}.
\]

## Experiment C — strict cusp chart

Use a separate \(U_{12}\) Grassmannian chart near

\[
Z_A=Z_i,\qquad Z_B=Z_{i-1}+\alpha Z_{i+1}.
\]

Test whether \(\alpha\) is the same as the contour variable. The code is expected to show that it is not generically: \(\alpha\) is tangential while the contour variable can live in projectivized normal data.

## Experiment D — linking contour

Complexify \(\eta\in\mathbb C^\times\) and use

\[
C_B: \eta=Re^{i\varphi}.
\]

Because

\[
z_{\rm cross}=\frac{z_0}{x}\eta,
\]

we have

\[
\oint_{C_B}\frac{d\eta}{\eta}
=
\oint_{C_z}\frac{dz}{z}
=2\pi i
\]

for the corresponding homology class.

This is the strongest one-loop evidence in favour of interpreting \(C_B\) as a linking cycle in the exceptional divisor.

## Experiment E — two-loop collective scaling

Use

\[
\widetilde I_4^{(2)}
=
-\frac{n_{12}}{D_1D_2\langle ABCD\rangle}
\]

in the normalized four-point chart.

Scale

\[
y_i=\rho a_i,\qquad w_i=\rho b_i.
\]

The tree-graph power count predicts

\[
I^{(L)}\sim\rho^{-2L},
\qquad
d^{2L}q\sim\rho^{2L-1}d\rho,
\]

so

\[
\Omega^{(L)}\sim\frac{d\rho}{\rho}\wedge\omega_{\rm exc}^{(L)}.
\]

At L=2 the exact code expands the full integrand and reports the leading coefficient after factoring \(\rho^{-4}\).

## Experiment F — perturbative contour benchmark

Use the established finite functions

\[
F^{(0)}(z)=-1,
\qquad
F^{(1)}(z)=\log^2z+\pi^2.
\]

Evaluate the circle-contour averages numerically at increasing precision and compare against analytic moments.

Do not use the known \(\Gamma_{\rm cusp}\) coefficient as an input to the geometric period.

## Experiment G — matrix substrate

Sample the Euclideanized matrix action

\[
S_E=-\frac{N}{4g^2}\sum_{\mu,\nu}\mathrm{Tr}[X_\mu,X_\nu]^2
\]

using Langevin dynamics. Track commutator density and eigenvalue spectra.

This probes the proposed substrate but does not establish a continuum gauge dictionary by itself.

## Experiment H — downstream physical layers

Export geometric effective fields to:

- FEniCSx/Firedrake/deal.II: boundary PDEs,
- FiPy/OpenFOAM: entropy-source transport,
- Meep/COMSOL: gauge/EM propagation,
- Kwant: boundary quantum transport,
- QuTiP: open-system entropy dynamics,
- VASP: vacuum-field import,
- Geant4: effective-mass phenomenology.

Each backend is an independent consistency check, not a substitute for the exact geometry.
