"""Resolve the insertion crossing on the D3-D4 null edge at p_z=-1/4."""

from math import pi, sqrt

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .energy_contour import energy_poles, finite_frame, interval, upper_residues
from .null_insertion_crossing import (
    RADIUS,
    WIDTH,
    double_pole_integral,
    integrate_deformed_cylinder,
)

GROUP_CENTER = 1.8
GROUP_RADIUS = 0.19
GROUP_OFFSETS = GROUP_RADIUS*np.exp(2j*pi*np.arange(32)/32)
CHECK_REGULATOR = 1e-4


def paired_energy(u, r, x, eps, centers):
    """Normalized lower energy residue sum, averaged across p_x=0."""
    z = -(1+u)/4
    dz = z[..., None]-centers[:, 3]
    small = r*r-1j*eps
    transverse_ratio = x*x/small

    def axis_factors(energy):
        return (energy[..., None]-centers[:, 0])**2-dz**2-small[..., None]

    # At the lower D4 pole, D3=-delta. At the lower D3 pole, D2=-2*delta.
    # Cancel the common transverse factor before numerical integration.
    pair_sum = 0j
    for length, axial_energy, slope, spectators in (
        (z+0.5, z, -1, (0, 1)),
        (-z, -z, -2, (0, 3)),
    ):
        omega = np.sqrt(length**2+small)
        displacement = small/(omega+length)
        factors = axis_factors(axial_energy+displacement)
        insertion = factors[..., 4]-0.25
        pair_sum += (omega+length)/(2*omega*slope*factors[..., spectators[0]]
                                     *factors[..., spectators[1]]*(insertion**2-x*x))

    group = np.zeros(np.broadcast_shapes(u.shape, x.shape), dtype=complex)
    for offset in GROUP_OFFSETS:
        factors = axis_factors(np.asarray(GROUP_CENTER+offset))
        insertion = factors[..., 4]-0.25
        group += offset/(np.prod(factors[..., :4], axis=-1)*(insertion**2-x*x))
    return 8/pi*(transverse_ratio*pair_sum + x*x*group/len(GROUP_OFFSETS))


def certify_deformation(centers):
    """Bound the moving poles on the two independent contour deformations."""
    z_center = -0.25
    z_shift = WIDTH/4
    dz = z_center-centers[:, 3]
    lengths = np.sqrt(dz**2+centers[:, 1]**2)
    small_bound = RADIUS**2+CHECK_REGULATOR
    squared_shift = 2*abs(dz)*z_shift+z_shift**2+small_bound
    squared_shift[4] += RADIUS
    assert np.all(squared_shift < lengths**2)
    shifts = squared_shift/lengths
    upper, lower = energy_poles([0, 0, z_center], centers, 0)

    included = [0, 1, 4]
    assert np.all(abs(lower[included]-GROUP_CENTER)+shifts[included] < GROUP_RADIUS)
    assert np.all(abs(lower[[2, 3]]-GROUP_CENTER)-shifts[[2, 3]] > GROUP_RADIUS)
    assert np.all(abs(upper-GROUP_CENTER)-shifts > GROUP_RADIUS)
    for i, crossing in ((2, (1,)), (3, (2, 4))):
        for j in range(len(centers)):
            if j != i:
                assert abs(lower[j]-lower[i]) > shifts[j]+shifts[i]
            if j not in crossing:
                assert abs(upper[j]-lower[i]) > shifts[j]+shifts[i]

    # D4's insertion factors are H(u) +/- x. On the real interval
    # Im H >= 0, and Re H' > 0 throughout the upper half-disk.
    min_length = 0.25-z_shift
    min_omega = sqrt(min_length**2-z_shift**2-small_bound)
    assert (21/32)*small_bound/(min_omega*(min_omega+min_length)) < 1
    assert WIDTH-(21/8)*small_bound/(min_omega+min_length) > RADIUS


def cylinder_coefficient(centers):
    R = sp.Rational
    u, energy, x = sp.symbols("u E x", real=True)
    z = -(1+u)/4
    d = [interval((energy, x, 0, z), center) for center in centers]
    assert sp.expand(d[2]-d[3]) == z-energy
    assert sp.expand(d[1]-d[2]) == -2*(energy+z)
    assert sp.expand(d[4]-d[3]+x) == -21*energy/8-11*z/8-1
    low = [sp.expand(f.subs({energy: z, x: 0})) for f in d]
    high = [sp.expand(f.subs({energy: -z, x: 0})) for f in d]
    assert low == [u-1, u+1, 0, 0, u]
    assert high == [u/2-R(3, 2), 0, 0, u/2+R(1, 2), 3*u/16-R(13, 16)]
    # The radial and azimuth integrations give 4/slope. The absolute
    # axial Jacobian is 1/4, because p_z decreases with u.
    f = sp.cancel(1/(low[0]*low[1]))
    assert sp.diff(f, u).subs(u, 0) == 0
    ordinary = -R(1, 2)/(high[0]*high[3]*high[4]**2)
    smooth = quad(sp.lambdify(u, ordinary, "numpy"), -WIDTH, WIDTH)[0]
    return smooth-double_pole_integral(f, u)


def check_second_null_insertion_crossing():
    exact = finite_frame()
    centers = np.asarray(exact, dtype=float)
    certify_deformation(centers)
    u, r = np.asarray(0.0), np.asarray(RADIUS)
    point = np.array([r, 0, -0.25])
    reflected = point*np.array([-1, 1, 1])
    upper = sum(upper_residues(p, centers, CHECK_REGULATOR).sum()
                for p in (point, reflected))/pi
    np.testing.assert_allclose(paired_energy(u, r, r, CHECK_REGULATOR, centers),
                               upper, rtol=0, atol=1e-7)

    def cylinder(radius, eps, order):
        def density(u, r, x):
            return paired_energy(u, r, x, eps, centers)
        return integrate_deformed_cylinder(radius, eps, order, density, 0.25)/radius**2

    target = cylinder_coefficient(exact)
    print(f"D3-D4 shrinking-cylinder coefficient: {target:.9f}")
    errors = []
    for radius in (RADIUS, 1/128, 1/256, 1/512, 1/1024):
        value = cylinder(radius, 0, 16)
        errors.append(abs(value-target))
        print(f"R={radius:.6f}: integral/R^2={value:.9f}")
    np.testing.assert_allclose(cylinder(radius, 0, 24), value, rtol=0, atol=1e-8)
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 0.02

    boundary = cylinder(RADIUS, 0, 20)
    errors = []
    for eps in (1e-7, 1e-8, 1e-9):
        value = cylinder(RADIUS, eps, 20)
        errors.append(abs(value-boundary))
        print(f"eps={eps:.0e}: scaled boundary error={errors[-1]:.3e}")
    np.testing.assert_allclose(cylinder(RADIUS, eps, 28), value, rtol=0, atol=2e-4)
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 0.01
    print("PASS: second null-insertion crossing and the real area coefficient")


if __name__ == "__main__":
    check_second_null_insertion_crossing()
