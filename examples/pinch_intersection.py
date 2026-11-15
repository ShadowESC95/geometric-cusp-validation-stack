"""Resolve two intersecting timelike Feynman pinch surfaces at fixed p_z."""

from math import pi

import numpy as np
import sympy as sp

from .energy_contour import finite_frame, finite_numerator, interval, upper_residues

HALF_WIDTH = 1/80
REGULAR_LOWER_POLES = (0, 1, 2, 4)


def patch_coordinates(u, v):
    """Use D2 and the insertion denominator as coordinates on D4=0."""
    energy, x = np.broadcast_arrays((1-u)/3, 7*u/8-v-0.5)
    y = np.sqrt((energy+0.5)**2-0.25-x*x)
    spatial = np.stack((x, y, np.full_like(x, -1)), axis=-1)
    return energy, spatial, (energy+0.5)/(3*y)


def lower_pole_terms(u, v, centers, eps):
    """Return lower-half-plane residues with the spatial Jacobian included."""
    _, spatial, jacobian = patch_coordinates(u, v)
    distances = np.sum((spatial[..., None, :] - centers[:, 1:])**2, axis=-1)
    omega = np.sqrt(distances-1j*eps)
    poles = centers[:, 0]+omega

    # Cancel the propagator whose energy pole is being taken.
    denominators = (poles[..., :, None]-centers[:, 0])**2-distances[..., None, :]+1j*eps
    index = np.arange(len(centers))
    denominators[..., index, index] = 1
    return (-4*spatial[..., 0]*jacobian)[..., None] / (
        2*omega*np.prod(denominators, axis=-1)
    )


def numerator_density(u, v, centers):
    """Real smooth factor multiplying (u+i0)^-1 (v+i0)^-1 at the D4 pole."""
    energy, spatial, jacobian = patch_coordinates(u, v)
    distances = np.sum((spatial[..., None, :] - centers[:, 1:])**2, axis=-1)
    denominators = (energy[..., None]-centers[:, 0])**2-distances
    d1, d3 = denominators[..., 0], denominators[..., 2]
    return -4*spatial[..., 0]*jacobian/(2*(energy+0.5)*d1*d3)


def boundary_value(centers, n):
    """Evaluate the two principal values by reflection on a symmetric box."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    q = HALF_WIDTH*(1+nodes)/2
    w = HALF_WIDTH*weights/2
    u = q[:, None]
    v = q

    def density(u, v):
        return numerator_density(u, v, centers)

    double_pv = w @ (
        (density(u, v)-density(-u, v)-density(u, -v)+density(-u, -v))/(u*v)
    ) @ w
    single_pvs = w @ ((density(0, q)-density(0, -q)+density(q, 0)-density(-q, 0))/q)
    intersection = -pi**2*density(0, 0)
    return double_pv - 1j*pi*single_pvs + intersection, intersection


def integrate_lower_poles(centers, eps, n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    # Resolve both narrow regulated poles on the same quadrature grid.
    extent = np.arcsinh(HALF_WIDTH/eps)
    q = eps*np.sinh(extent*nodes)
    w = weights*extent*eps*np.cosh(extent*nodes)
    residues = lower_pole_terms(q[:, None], q, centers, eps)
    return np.einsum("i,ijk,j->k", w, residues, w)


def regular_boundary_value(centers, n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    q, w = HALF_WIDTH*nodes, HALF_WIDTH*weights
    terms = lower_pole_terms(q[:, None], q, centers, eps=0)
    regular = terms[..., REGULAR_LOWER_POLES].sum(axis=-1)
    return w @ regular @ w


def check_pinch_intersection():
    R = sp.Rational
    exact_centers = finite_frame()
    E, x, y, z = sp.symbols("E x y z", real=True)
    loop = (E, x, y, z)
    d = [sp.expand(interval(loop, center)) for center in exact_centers]
    u, v = sp.symbols("u v", real=True)
    upper_2, lower_4, upper_insertion = d[1], d[3], d[4]
    affine_equations = [upper_2-lower_4-u, upper_insertion-lower_4-v]
    assert sp.solve(affine_equations, (E, x)) == {E: -(u+z)/3, x: 7*u/8-v-z/2-1}
    assert sp.factor(sp.Matrix([d[1], d[3], d[4]]).jacobian([E, x, y]).det()) == 6*y

    point = {E: R(1, 3), x: -R(1, 2), y: sp.sqrt(7)/6, z: -1}
    assert [propagator.subs(point) for propagator in d] == [-R(2, 3), 0, -R(4, 3), 0, 0]
    assert finite_numerator(loop, exact_centers).subs(point) == 2
    assert exact_centers[1][0] > point[E] > exact_centers[3][0]
    assert exact_centers[4][0] > point[E]

    # At the lower D4 pole the other propagators inherit the imaginary shift
    # eps*(1-(E-E_i)/(E-E_4)), positive for both upper-pole propagators.
    shifts = [sp.simplify(1-(point[E]-exact_centers[i][0])/(point[E]-exact_centers[3][0]))
              for i in (1, 4)]
    assert shifts == [R(9, 5), R(63, 40)]
    exact_contact = -9*sp.pi**2/(4*sp.sqrt(7))

    centers = np.array(exact_centers, dtype=float)
    target, contact = boundary_value(centers, n=32)
    np.testing.assert_allclose(contact, float(exact_contact), rtol=0, atol=1e-12)
    np.testing.assert_allclose(boundary_value(centers, n=24)[0], target, rtol=0, atol=1e-9)
    regular_target = regular_boundary_value(centers, n=32)
    np.testing.assert_allclose(
        regular_boundary_value(centers, n=24), regular_target, rtol=0, atol=1e-9,
    )

    # Confirm the closing sign by comparing upper and lower residue sums
    # off the intersection, where individual upper residues are well separated.
    pu, pv = HALF_WIDTH/2, -HALF_WIDTH/2
    _, spatial, jacobian = patch_coordinates(pu, pv)
    np.testing.assert_allclose(
        lower_pole_terms(pu, pv, centers, 1e-3).sum(),
        -jacobian*upper_residues(spatial, centers, 1e-3).sum(),
        rtol=1e-11, atol=1e-10,
    )

    print(f"Double-delta term={contact:.9f}, boundary value={target:.9f}")
    errors = []
    for eps in (1e-4, 1e-5, 1e-6):
        lower = integrate_lower_poles(centers, eps, n=96)
        value = lower[3]
        other_poles = lower[list(REGULAR_LOWER_POLES)].sum()
        errors.append(abs(value-target))
        print(f"eps={eps:.0e}: D4 residue integral={value:.9f}, "
              f"other poles Im={other_poles.imag:.3e}")
    np.testing.assert_allclose(
        integrate_lower_poles(centers, eps=1e-6, n=72), lower, rtol=0, atol=2e-6,
    )
    assert errors[2] < errors[1] < errors[0] and errors[2] < 0.0015
    np.testing.assert_allclose(other_poles, regular_target, rtol=0, atol=1e-6)
    raw_limit = -2j*pi*(target + regular_target)
    np.testing.assert_allclose(-2j*pi*lower.sum(), raw_limit, rtol=0, atol=0.01)
    print(f"PASS: intersecting pinches, raw patch energy integral limit={raw_limit:.9f}")


if __name__ == "__main__":
    check_pinch_intersection()
