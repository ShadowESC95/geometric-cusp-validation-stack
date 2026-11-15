"""An explicit capped cycle on the weighted divisor and its scalar pushforward."""

from math import pi

import numpy as np
import sympy as sp

from geocusp.benchmark import cusp_coefficient_from_mean
from geocusp.residue import strict_cusp_weighted_blowup

from .complex_cross_ratio import integrated_projective
from .geometry import EXTERNAL_PAIRS, bracket


def certify_charts():
    a, b, c, d, r, u, v, z = sp.symbols("a b c d r u v z", nonzero=True)
    line = sp.Matrix([[1, 0, a, b], [0, 1, c, d]])
    form = 1/(a*d*(a*d-b*c))
    p = [bracket(line, *pair) for pair in EXTERNAL_PAIRS]
    assert sp.cancel(p[0]*p[2]/(p[1]*p[3])-(b*c-a*d)/(a*d)) == 0

    # The weight-two chart has a double cover with residual (u,v)->(-u,-v).
    lift = {b: r*u, c: r*v, d: r**2}
    jac = sp.Matrix([a, r*u, r*v, r**2]).jacobian([a, u, v, r]).det()
    residue = sp.cancel(r*jac*form.subs(lift))
    assert sp.cancel(residue-2/(a*(a-u*v))) == 0
    assert sp.cancel(residue.subs({u: -u, v: -v}, simultaneous=True)-residue) == 0
    # At the cap u=v=0 the coefficient is finite on the cover.
    assert residue.subs({u: 0, v: 0}) == 2/a**2

    # On the overlap b=-rho, the weight-two coordinate is D=1/u^2.
    # Check the orientation in radial-last convention.
    local_residue, _ = strict_cusp_weighted_blowup()
    local_a, local_z, local_D = (
        next(s for s in local_residue.free_symbols if s.name == name)
        for name in ("a", "z", "D")
    )
    overlap = sp.Matrix([a, u*v/a-1, 1/u**2])
    overlap_jac = overlap.jacobian([a, u, v]).det()
    assert sp.cancel(local_residue.subs(dict(zip(
        (local_a, local_z, local_D), overlap))) * overlap_jac-residue) == 0

    # Independently change gauge to columns (3,2). The invariant cross-ratio
    # remains z and the weight-two chart transition is regular for a,z != 0.
    gauge = line[:, [2, 1]].inv()*line
    other = [bracket(gauge, *pair) for pair in EXTERNAL_PAIRS]
    assert gauge == sp.Matrix([[1/a, 0, 1, b/a], [-c/a, 1, 0, d-b*c/a]])
    assert sp.cancel(other[0]*other[2]/(other[1]*other[3])-p[0]*p[2]/(p[1]*p[3])) == 0
    coords = sp.Matrix([gauge[0, 0], gauge[0, 3], gauge[1, 0], gauge[1, 3]])
    a_new, b_new, c_new, d_new = coords
    new_form = 1/(a_new*d_new*(a_new*d_new-b_new*c_new))
    assert sp.cancel(coords.jacobian([a, b, c, d]).det()*new_form-form) == 0
    assert sp.cancel(d_new.subs(lift).subs(v, a*(z+1)/u)+z*r**2) == 0

    # The orbit on the weight-two cover uses theta in [0, pi]. The quotient
    # identifies its two endpoints, so this is a single cycle downstairs.
    # The u*v=a(z+1) relation gives du^dv=(a/u)du^dz modulo da.
    parameterization = sp.Matrix([a, u, a*(z+1)/u])
    pullback = sp.cancel(residue.subs(v, parameterization[2])
                         *parameterization.jacobian([a, z, u]).det())
    assert sp.cancel(pullback-2/(a*z*u)) == 0

    # On the affine Lorentzian real slice c=-2*conj(b), a,d real.
    # In the b chart with b=-R*exp(i theta), d=R^2*D_real,
    # C=2*exp(-2i theta) and D=D_real*exp(-2i theta).
    R, D_real = sp.symbols("R D_real", positive=True)
    phase = sp.symbols("phase", nonzero=True)
    real_slice = {b: -R*phase, c: 2*R/phase, d: R**2*D_real}
    assert sp.cancel(((b*c-a*d)/(a*d)).subs(real_slice)+1+2/(a*D_real)) == 0


def lifted_cycle(a, z, theta):
    """Representative on the cover of the weight-two chart d=1."""
    product = a*(z+1)
    u = np.sqrt(abs(product))*np.exp(1j*theta)
    if product == 0:
        return 0j, 0j
    return u, product/u


def check_cycle_and_period():
    certify_charts()
    angles = np.linspace(-pi, pi, 65)[1:-1]
    for a_radius in (0.7, 1.3):
        for alpha in (0., pi/3, pi):
            a = a_radius*np.exp(1j*alpha)
            for phi in angles:
                z = np.exp(1j*phi)
                u, v = lifted_cycle(a, z, 0)
                end_u, end_v = lifted_cycle(a, z, pi)
                np.testing.assert_allclose(u*v, a*(z+1), atol=1e-13, rtol=0)
                np.testing.assert_allclose((end_u, end_v), (-u, -v), atol=1e-13, rtol=0)
                assert abs(a-u*v) > 0.69
            assert lifted_cycle(a, -1, 0) == (0j, 0j)

    # In (alpha,phi,theta) order the pullback is 2*i^3 times the
    # angular differentials. The real radius derivative of u drops out.
    # The a circle has length 2*pi and the quotient theta circle length pi.
    normalized_fiber = 2*(2j*pi)*(1j*pi)/(2j*pi)**2
    # This number is the factor multiplying the normalized phi integral.
    np.testing.assert_allclose(normalized_fiber, 1, atol=1e-14, rtol=0)

    nodes, weights = np.polynomial.legendre.leggauss(96)
    z_values = np.exp(1j*pi*nodes)
    np.testing.assert_allclose(weights @ z_values/2, 0, rtol=0, atol=1e-13)
    np.testing.assert_allclose(weights @ np.ones_like(z_values)/2, 1, rtol=0, atol=1e-13)
    # Evaluate the CD weight by the independently reduced integral, not by
    # imposing the claimed contour coefficient on the cycle.
    mean = normalized_fiber*(weights @ np.array([
        integrated_projective(z) for z in z_values
    ]))/2
    np.testing.assert_allclose(mean, 2*pi*pi/3, rtol=0, atol=1e-9)
    gamma1 = cusp_coefficient_from_mean(1, -normalized_fiber)
    gamma2 = cusp_coefficient_from_mean(2, mean)
    np.testing.assert_allclose([gamma1, gamma2], [4, -4*pi*pi/3], rtol=0, atol=1e-8)
    print(f"Oriented quotient fiber={normalized_fiber.real:+.0f}, "
          f"CD-weighted period={mean.real:.12f}")
    print(f"Benchmark conversion: Gamma[1]={gamma1.real:.9f}, Gamma[2]={gamma2.real:.9f}")
    print("PASS: capped weighted-divisor cycle and scalar pushforward")


if __name__ == "__main__":
    check_cycle_and_period()
