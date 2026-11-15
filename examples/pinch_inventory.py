"""Classify opposite-side energy pinches and resolve the insertion center."""

from itertools import combinations
from math import pi

import numpy as np
import sympy as sp

from .energy_contour import energy_poles, finite_frame, finite_numerator, interval, upper_residues
from .null_endpoint import MAX_ETA, RADIUS, integrate_ball

REGULAR_CONTOURS = ((-1.32, 0.16, (2, 3)), (-0.09, 0.12, (0,)), (0.466, 0.065, (1,)))


def classify_pairs(centers):
    """Classify collisions of opposite-side poles at zero regulator."""
    null, timelike = [], []
    for i, j in combinations(range(len(centers)), 2):
        separation = interval(centers[i], centers[j])
        if separation == 0:
            null.append((i, j))
        elif separation > 0:
            timelike.append((i, j))
    # The meeting condition is |p-p_i|+|p-p_j|=|E_i-E_j|.
    # The triangle inequality excludes spacelike pairs, yields the segment
    # in the null case, and gives a prolate ellipsoid for timelike pairs.
    return null, timelike


def check_pair_geometry(centers):
    null, timelike = classify_pairs(centers)
    assert null == [(0, 1), (0, 3), (1, 2), (2, 3)]
    assert timelike == [(1, 3), (3, 4)]
    for i, j in (*null, *timelike):
        center_norms = sum(sp.sqrt(sum(x*x for x in centers[k][1:])) for k in (i, j))
        assert (abs(centers[i][0]-centers[j][0])+center_norms)/2 <= sp.Rational(9, 4)

    energy, x, y, z, t = sp.symbols("E x y z t", real=True)
    d = [interval((energy, x, y, z), center) for center in centers]
    ellipsoids = []
    for i, j in timelike:
        meeting, = sp.solve(d[i]-d[j], energy)
        ellipsoid = sp.expand(-d[i].subs(energy, meeting))
        ellipsoids.append(ellipsoid)
        assert sp.hessian(ellipsoid, (x, y, z)).is_positive_definite
    assert [sp.solve(p.subs({x: 0, y: 0}), z) for p in ellipsoids] == [
        [-sp.Rational(3, 2), 0], [-sp.Rational(13, 10), -sp.Rational(1, 4)],
    ]
    insertion_crossings = []
    for i, j in null:
        segment = tuple(a+t*(b-a) for a, b in zip(centers[i], centers[j]))
        for root in sp.solve(interval(segment, centers[4]), t):
            if 0 < root < 1:
                insertion_crossings.append(((i, j), root, segment[3].subs(t, root)))
    assert insertion_crossings == [
        ((0, 3), sp.Rational(1, 5), -sp.Rational(13, 10)),
        ((2, 3), sp.Rational(1, 2), -sp.Rational(1, 4)),
    ]
    print("Opposite-side pinches: 4 null segments and 2 timelike ellipsoids within |p| <= 9/4")
    print("D4-insertion ellipsoid crosses the null axis at p_z=-13/10 and -1/4")


def insertion_residues(point, eps, centers):
    """Return the soft upper insertion pole and the sum of the regular groups."""
    distances_sq = np.sum((point[..., None, :]-centers[:, 1:])**2, axis=-1)
    omega = np.sqrt(distances_sq[..., 4]-1j*eps)
    energy = centers[4, 0]-omega

    numerator = -4*point[..., 0]
    factors = (energy[..., None]-centers[:4, 0])**2-distances_sq[..., :4]+1j*eps
    soft = -numerator/(pi*omega*np.prod(factors, axis=-1))
    return soft, regular_residues(distances_sq, numerator, eps, centers)


def regular_residues(distances_sq, numerator, eps, centers):
    smooth = np.zeros_like(numerator, dtype=complex)
    for center, radius, _ in REGULAR_CONTOURS:
        for offset in radius*np.exp(2j*pi*np.arange(32)/32):
            factors = (center+offset-centers[:, 0])**2-distances_sq+1j*eps
            smooth += offset/np.prod(factors, axis=-1)
    return (2/pi)*numerator*smooth/32


def insertion_ball(radius, eta, centers, order):
    def density(r, mu, cos_phi, eps, _):
        sin_theta = np.sqrt(1-mu**2)
        x = r*sin_theta*cos_phi
        y = r*sin_theta*np.sqrt(1-cos_phi**2)
        z = np.broadcast_to(r*mu, x.shape)
        positive = np.stack((x, y, z), axis=-1)
        negative = positive*np.array([-1, 1, 1])
        # The ball is centered on the insertion. Pair the two x offsets
        # before using the quadrant rule in integrate_ball.
        first = insertion_residues(positive+centers[4, 1:], eps, centers)
        second = insertion_residues(negative+centers[4, 1:], eps, centers)
        return sum(first+second)/2

    return integrate_ball(radius, eta, centers, order, density)


def certify_insertion_poles(centers):
    upper, lower = energy_poles(centers[4, 1:], centers, 0)
    shift = RADIUS*(1+np.sqrt(MAX_ETA))
    for center, radius, indices in REGULAR_CONTOURS:
        for branch, is_upper in ((upper, True), (lower, False)):
            for i, pole in enumerate(branch):
                if is_upper and i in indices:
                    assert abs(pole-center)+shift < radius
                else:
                    assert abs(pole-center)-shift > radius
    for branch in (upper, lower):
        assert np.all(abs(branch[:4]-centers[4, 0]) > 2*shift)


def insertion_coefficients(centers):
    r, x, y, z = sp.symbols("r x y z", real=True)
    E0, x0, y0, z0 = centers[4]
    on_shell = E0-r, x0+r*x, y0+r*y, z0+r*z
    factors = [sp.expand(interval(on_shell, p)) for p in centers[:4]]
    constants = [f.subs(r, 0) for f in factors]
    assert constants == [-sp.Rational(1, 4), -sp.Rational(1, 4), -1, 1]
    assert sp.expand(finite_numerator(on_shell, centers)) == 2-4*r*x
    scaled = [sp.expand(f/c) for f, c in zip(factors, constants)]
    linear = [f.coeff(r, 1) for f in scaled]
    quadratic = [f.coeff(r, 2) for f in scaled]
    # Expand the smooth factor of the residue, with the 1/(pi*r) removed.
    coefficient = -2/sp.prod(constants)
    first = sp.expand(coefficient*(-sum(linear)-2*x))
    assert first.subs({x: 0, y: 0, z: 0}) == 0
    second = sp.expand(coefficient*((sum(linear)**2+sum(a*a for a in linear))/2
                                    - sum(quadratic)+2*x*sum(linear)))
    sphere_mean = second.subs({x: 0, y: 0, z: 0})+sum(second.coeff(a, 2) for a in (x, y, z))/3
    assert 2*coefficient == 64
    assert sp.simplify(sphere_mean-726) == 0
    return float(2*coefficient), float(sphere_mean)


def check_pinch_inventory():
    exact = finite_frame()
    check_pair_geometry(exact)
    leading, quartic = insertion_coefficients(exact)
    centers = np.asarray(exact, dtype=float)
    certify_insertion_poles(centers)
    eps = MAX_ETA*RADIUS**2
    probe = centers[4, 1:]+np.array([1, 2, 1])*RADIUS/3
    np.testing.assert_allclose(sum(insertion_residues(probe, eps, centers)),
                               2/pi*upper_residues(probe, centers, eps).sum(),
                               rtol=1e-10, atol=1e-8)

    # The regular residues are smooth on the ball. Their leading term is
    # the value at the center times its volume, while the linear term averages out.
    center = centers[4, 1:]
    distances_sq = np.sum((center-centers[:, 1:])**2, axis=-1)
    regular = regular_residues(distances_sq, -4*center[0], 0, centers)
    cubic = 4*pi/3*regular
    assert abs(cubic.imag) < 1e-9
    print(f"insertion ball cubic coefficient (regular residues): {cubic.real:.9f}")
    errors = []
    for radius in (RADIUS, 1/128, 1/512, 1/2048, 1/4096):
        value = insertion_ball(radius, 0, centers, order=16)
        correction = (value-leading*radius**2-cubic*radius**3)/radius**4
        errors.append(abs(correction-quartic))
        print(f"insertion ball R={radius:.6f}: integral/R^2={value.real/radius**2:.9f}, "
              f"quartic estimate={correction.real:.8f}")
    np.testing.assert_allclose(insertion_ball(radius, 0, centers, order=24), value,
                               rtol=0, atol=1e-11)
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 5

    eta = MAX_ETA
    # The leading scaled integral keeps only the 32/(pi*sqrt(r^2-i*eps)) term.
    nodes, weights = np.polynomial.legendre.leggauss(64)
    t = (1+nodes)/2
    target = 64*(weights @ (t*t/np.sqrt(t*t-1j*eta)))
    errors = []
    for radius in (RADIUS, 1/128, 1/512, 1/2048, 1/8192):
        value = insertion_ball(radius, eta, centers, order=28)/radius**2
        errors.append(abs(value-target))
    np.testing.assert_allclose(insertion_ball(radius, eta, centers, order=36)/radius**2,
                               value, rtol=0, atol=1e-9)
    print(f"eta={eta:g}: scaled limit={target:.9f}, smallest-ball result={value:.9f}")
    assert all(new < old for old, new in zip(errors, errors[1:])) and errors[-1] < 0.1
    print("PASS: pinch inventory and the insertion-center ball")


if __name__ == "__main__":
    check_pinch_inventory()
