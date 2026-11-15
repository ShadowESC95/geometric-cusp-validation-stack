"""Check local integrability across two projected external null edges."""

from math import pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from .energy_contour import direct_energy_integral, finite_frame, finite_numerator, interval

SEGMENT = (1/4, 3/8)
MAX_RADIUS = 1/32
MAX_REGULATOR = 1e-6
DIRECT_CHECK_REGULATOR = 1e-4
PINCHED_PAIRS = ((1, 2), (0, 3))
REGULAR_CONTOURS = ((-11/16, 1/8, (2, 3)), (1/10, 1/8, (4,)))


def paired_energy(lam, r, x, eps, centers):
    """Reflection average of the normalized energy integral at x and -x."""
    small = r*r-1j*eps
    dz = (lam-1)[..., None]-centers[:, 3]
    lengths = np.abs(dz[..., :4])
    omega = np.sqrt(lengths**2+small[..., None])
    displacement = small[..., None]/(omega+lengths)
    total = np.zeros(np.broadcast_shapes(lam.shape, x.shape), dtype=complex)
    for i, j in PINCHED_PAIRS:
        pole = centers[i, 0]-omega[..., i]
        # The two real branch energies agree at r=eps=0. Keeping the
        # difference of displacements avoids cancellation as r tends to zero.
        gap = -(displacement[..., i]+displacement[..., j])
        other = pole-centers[j, 0]+omega[..., j]
        spectators = [k for k in range(4) if k not in (i, j)]
        d = (pole[..., None]-centers[spectators, 0])**2-dz[..., spectators]**2-small[..., None]
        insertion = (pole-centers[4, 0])**2-dz[..., 4]**2-1/4-small
        total += 1/(-2*omega[..., i]*gap*other*np.prod(d, axis=-1)*(insertion**2-x*x))

    # The regular contours sum poles on the same side of the real axis.
    directions = np.exp(2j*pi*np.arange(32)/32)
    for center, radius, _ in REGULAR_CONTOURS:
        subtotal = np.zeros_like(total)
        for offset in radius*directions:
            d = (center+offset-centers[:, 0])**2-dz**2-small[..., None]
            insertion = d[..., 4]-1/4
            subtotal += offset/(np.prod(d[..., :4], axis=-1)*(insertion**2-x*x))
        total += subtotal/len(directions)

    # Reflection combines the two insertion factors and gives numerator
    # -4*x^2 in the real-energy integrand. Upper closure adds 2/pi.
    return -8/pi*x*x*total


def integrate_cylinder(radius, eps, centers, n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    lo, hi = SEGMENT
    lam = (lo+hi)/2+(hi-lo)*nodes[:, None, None]/2
    phi = pi*nodes[None, None, :]/2
    # Resolve the transverse crossover at r comparable to sqrt(eps).
    cutoff = np.sqrt(eps)
    knots = sorted({0, min(cutoff, radius), min(4*cutoff, radius), radius})
    total = 0j
    for a, b in zip(knots, knots[1:]):
        r = (a+b)/2+(b-a)*nodes[None, :, None]/2
        density = paired_energy(lam, r, r*np.cos(phi), eps, centers)
        # Half the azimuth suffices by reflection; the full disk doubles it.
        total += np.einsum("i,j,k,ijk->", weights, weights, weights, r*density)*(b-a)/2
    return pi*(hi-lo)/2*total


def certify_contours(centers):
    """Bound the pole locations using the endpoints of the axial segment."""
    assert 4*np.sqrt(MAX_REGULATOR) < MAX_RADIUS
    # All pole branches at zero transverse momentum are monotone in lambda.
    # Along this segment dz has fixed sign for every center, and each omega
    # is increasing in |dz|. The transverse shift is at most sqrt(|delta|).
    all_centers = [*centers, (centers[4][0], -centers[4][1], 0, centers[4][3])]
    dz = np.array(SEGMENT)[:, None]-1-np.array(all_centers)[:, 3]
    assert np.all(dz[0]*dz[1] > 0)
    times = np.array(all_centers)[:, 0]
    transverse_centers = np.array(all_centers)[:, 1]
    omega = np.sqrt(dz*dz+transverse_centers**2)
    poles = [times-omega, times+omega]
    # For the insertion, Re(omega+omega_axis) >= omega_axis gives a tighter bound.
    insertion_shift = (MAX_RADIUS**2+MAX_RADIUS+DIRECT_CHECK_REGULATOR)/omega[:, 4:].min()
    external_shift = np.sqrt(MAX_RADIUS**2+DIRECT_CHECK_REGULATOR)
    margins = np.array([external_shift]*4+[insertion_shift]*2)
    for center, radius, indices in REGULAR_CONTOURS:
        enclosed = set(indices)
        if 4 in enclosed:
            enclosed.add(5)
        for side, branches in enumerate(poles):
            for index, (left, right) in enumerate(zip(branches.min(axis=0), branches.max(axis=0))):
                if side == 0 and index in enclosed:
                    assert max(abs(left-center), abs(right-center))+margins[index] < radius
                else:
                    assert right+margins[index] < center-radius or left-margins[index] > center+radius
    for i, j in PINCHED_PAIRS:
        assert np.allclose(poles[0][:, i], poles[1][:, j])
        # The pinching energy is affine with slope of magnitude 1.
        # Every other axial branch has derivative in [-1, 1], so each
        # signed difference is monotone and can be bounded at the endpoints.
        pinching = poles[0][:, i]
        for side, branches in enumerate(poles):
            for k in range(len(all_centers)):
                if (side, k) in ((0, i), (1, j)):
                    continue
                assert np.min(np.abs(pinching-branches[:, k])) > external_shift+margins[k]
                assert np.prod(pinching-branches[:, k]) > 0


def check_null_pinch():
    exact = finite_frame()
    centers = np.array(exact, dtype=float)
    certify_contours(centers)
    lam = sp.symbols("lambda", real=True)
    energy, x, y, z = sp.symbols("E x y z", real=True)
    assert sp.expand(finite_numerator((energy, x, y, z), exact)) == -4*x
    assert all(point[1:3] == (0, 0) for point in exact[:4])
    assert exact[4][1:3] == (-sp.Rational(1, 2), 0)
    spectators = []
    for (i, j), energy in zip(PINCHED_PAIRS, (1-lam, -lam)):
        point = (energy, 0, 0, lam-1)
        d = [sp.factor(interval(point, center)) for center in exact]
        assert d[i] == d[j] == 0
        spectators.append(sp.prod(d[k] for k in range(4) if k not in (i, j))*d[4]**2)
    assert sp.factor(spectators[0]) == lam*(lam-1)*(3*lam+1)**2/4
    assert sp.factor(spectators[1]) == (2*lam-1)*(2*lam+1)*(10*lam+3)**2/64

    # -2*omega_i*gap*other tends to 2*L*(r^2-i*eps),
    # with L the axial separation of the two centers. The angular
    # integral of x^2/r^2 is pi, and the radial measure gives radius^2/2.
    coefficient = sum(
        -2/(abs(exact[i][3]-exact[j][3])*product)
        for (i, j), product in zip(PINCHED_PAIRS, spectators)
    )
    predicted = quad(sp.lambdify(lam, coefficient, "numpy"), *SEGMENT, epsabs=1e-12)[0]

    # Independent check from the two unsymmetrized real-axis integrals.
    p = np.array([MAX_RADIUS, 0, sum(SEGMENT)/2-1])
    reflected = p*np.array([-1, 1, 1])
    direct = sum(direct_energy_integral(q, centers, DIRECT_CHECK_REGULATOR)
                 for q in (p, reflected))/(2j*pi**2)
    closed = paired_energy(np.array(p[2]+1), np.array(p[0]), np.array(p[0]),
                           DIRECT_CHECK_REGULATOR, centers)
    np.testing.assert_allclose(closed, direct, rtol=0, atol=1e-7)

    print(f"Predicted shrinking-cylinder coefficient: {predicted:.12f}")
    errors = []
    for radius in (MAX_RADIUS, 1/64, 1/128, 1/256):
        boundary = integrate_cylinder(radius, 0, centers, n=12)/radius**2
        errors.append(abs(boundary-predicted))
        print(f"radius={radius:.6f}: normalized integral/radius^2={boundary.real:.12f}")
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 2e-5

    boundary = integrate_cylinder(MAX_RADIUS, 0, centers, n=16)
    errors = []
    for eps in (1e-6, 1e-7, 1e-8):
        regulated = integrate_cylinder(MAX_RADIUS, eps, centers, n=20)
        errors.append(abs(regulated-boundary)/MAX_RADIUS**2)
        print(f"eps={eps:.0e}: scaled difference from boundary={errors[-1]:.3e}")
    refined = integrate_cylinder(MAX_RADIUS, 1e-8, centers, n=28)
    assert abs(refined-regulated)/MAX_RADIUS**2 < 2e-6
    assert errors[2] < errors[1] < errors[0] and errors[2] < 4e-4
    print("PASS: external null-edge cylinder and its zero-regulator limit")


if __name__ == "__main__":
    check_null_pinch()
