# Current local status

The v3 stack is designed to expose, rather than hide, failures.

## Exact one-loop findings

The source chart canonical form is

```math
\Omega_4^{(1)}=d\log x\wedge d\log y\wedge d\log z_0\wedge d\log w.
```

For the adjacent-propagator blow-up

```math
y=\rho\eta,
\qquad
w=\rho,
```

the exceptional residue is nonzero and proportional to

```math
\frac{dx}{x}\wedge\frac{d\eta}{\eta}\wedge\frac{dz_0}{z_0}.
```

The cross-ratio is

```math
z=\frac{z_0}{x}\eta,
```

so the projective normal coordinate and the established contour variable are related by an exact multiplicative map.

For the homogeneous blow-up of the strict collinear line, the radial residue is zero at a generic point. The weighted blow-up with weights (1, 1, 2) in (b, c, d) yields a nonzero radial residue and preserves a variable cross-ratio. The one-loop CLI reports both checks.

## Two-loop finding

The collective blow-up gives the rational scalar scaling

```math
I_4^{(2)}\sim \rho^{-4},
```

with radial Jacobian $`\rho^3d\rho`$, giving a logarithmic differential-form scaling.

The mutual bracket couples the two loop directions, so independent one-loop residue factorization is not a valid two-loop prescription.

## Validation philosophy

These are not claims that the user's full theory is proven. They are computationally testable structural results from the specified model/benchmark.
