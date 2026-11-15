"""Check the two null vertices whose incident edges point in opposite axial directions."""

from math import pi

import numpy as np
import sympy as sp

from .energy_contour import energy_poles, finite_frame, interval, upper_residues
from .null_endpoint import MAX_ETA, RADIUS
from .spatial_tail import reflected_centers

# Use E and p_z measured from the vertex, reversing both at D3.
VERTICES = ((0, 1), (2, -1))
VERTEX_SPECTATOR = (2, 0)
REGULAR_GROUPS = (((-2, 0.125, (2, 3)), (-0.277, 0.125, (4, 5))),
                  ((-2.05, 0.18, (0, 1, 4, 5)),))


def local_centers(centers, case):
    vertex, sign = VERTICES[case]
    return [(sign*(e-centers[vertex][0]), x, y, sign*(z-centers[vertex][3]))
            for e, x, y, z in centers]


def energy_density(r, mu, cos_phi, eps, centers, case):
    """Normalized upper residue sum for the reflected local integrand."""
    a = r*mu
    transverse_sq = r*r*(1-mu)*(1+mu)
    x_sq = transverse_sq*cos_phi**2
    small = transverse_sq-1j*eps
    omega = np.sqrt(r*r-1j*eps)
    ratio = cos_phi**2 if eps == 0 else x_sq/small
    dz = a[..., None]-centers[:, 3]

    def denominators(energy):
        return (energy[..., None]-centers[:, 0])**2-dz**2-small[..., None]

    def spectators(d):
        insertion = d[..., 4]-1/4
        return d[..., VERTEX_SPECTATOR[case]]*(insertion**2-x_sq)

    # The soft pole has null product -2*small. The upper pole of the
    # future neighbor has endpoint denominator -delta and past denominator
    # 4*a-3*delta, where delta is its displacement from the axial pole.
    soft = ratio/(4*omega*spectators(denominators(-omega)))
    length = 0.5-a
    neighbor_omega = np.sqrt(length*length+small)
    delta = small/(neighbor_omega+length)
    future = ratio*(neighbor_omega+length)/(
        2*neighbor_omega*(4*a-3*delta)*spectators(denominators(a-delta))
    )
    regular = np.zeros(np.broadcast_shapes(a.shape, x_sq.shape), dtype=complex)
    for center, radius, _ in REGULAR_GROUPS[case]:
        for offset in radius*np.exp(2j*pi*np.arange(32)/32):
            d = denominators(np.asarray(center+offset))
            insertion = d[..., 4]-0.25
            regular += offset/(np.prod(d[..., :4], axis=-1)*(insertion**2-x_sq))
    return -8/pi*(soft + future + x_sq*regular/32)


def angular_contour(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    theta = pi*(1-nodes)/2
    mu = np.exp(1j*theta)
    return mu, -1j*pi/2*mu*weights


def integrate_ball(radius, eta, centers, case, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (1+nodes[:, None, None])/2
    phi = pi*(1+nodes[None, None, :])/4
    mu, angular_weights = angular_contour(order)
    density = energy_density(radius*t, mu[None, :, None], np.cos(phi),
                             eta*radius**2, centers, case)
    return radius*pi/2*np.einsum("i,j,k,ijk->", weights, angular_weights, weights,
                                t*t*density)


def certify_deformation(centers, case):
    extended = np.asarray(reflected_centers(centers), dtype=float)
    upper, lower = energy_poles([0, 0, 0], extended, 0)
    # For |mu|<=1, bound the change in each squared spatial distance.
    # The same bounds apply to either reflected insertion.
    r = RADIUS
    small_bound = (2+MAX_ETA)*r*r
    distances = (lower-upper).real/2
    distance_shift = 2*abs(extended[:, 3])*r+r*r+small_bound
    distance_shift[4:] += np.sqrt(2)*r
    nonzero = distances > 0
    assert np.all(distance_shift[nonzero] < distances[nonzero]**2)
    shifts = np.full(len(extended), r*np.sqrt(1+MAX_ETA))
    shifts[nonzero] = distance_shift[nonzero]/distances[nonzero]
    for center, radius, included in REGULAR_GROUPS[case]:
        for branch, is_upper in ((upper, True), (lower, False)):
            for i, pole in enumerate(branch):
                if is_upper and i in included:
                    assert abs(pole-center)+shifts[i] < radius
                else:
                    assert abs(pole-center)-shifts[i] > radius
    soft_index, future_index = (0, 1) if case == 0 else (2, 3)
    for branch in (upper, lower):
        for index in (VERTEX_SPECTATOR[case], 4, 5):
            assert abs(branch[index]) > max(shifts[soft_index], shifts[future_index])+shifts[index]

    # For eps>0 the roots of small=0 have |mu|>1.
    # At eps=0 the transverse factor has already been canceled.
    # For the future pole's timelike denominator H(a)=4*a-3*delta,
    # Im H>=0 on the real axis. Bound |delta'| there and in the half-disk
    # using neighbor_omega=sqrt(1/4-a+r^2-i*eps).
    omega_min = np.sqrt(0.25-r-MAX_ETA*r*r)
    assert 3*(r+(1+MAX_ETA)*r*r)/(omega_min*(omega_min+0.5)) < 4
    # At eps=0, on the arc |a|=r, 3*|delta|<4*r uniformly down to r=0.
    assert 6*r/(omega_min+0.5-r) < 4


def endpoint_target(centers, case):
    vertex = VERTICES[case][0]
    future, past = (1, 3) if case == 0 else (3, 1)
    R = sp.Rational
    E, a, x, y, delta, omega = sp.symbols("E a x y delta omega", real=True)
    d = [interval((E, x, y, a), center) for center in centers]
    assert sp.expand(d[future]-d[vertex]) == a-E
    assert sp.expand(d[past]-d[vertex]) == 2*(E+a)
    assert sp.expand((d[vertex]-d[future]).subs(E, a-delta)) == -delta
    assert sp.expand((d[past]-d[future]).subs(E, a-delta)) == 4*a-3*delta
    pair = ((d[future]-d[vertex])*(d[past]-d[vertex])).subs(E, -omega)
    assert sp.expand(pair+2*(omega**2-a*a)) == 0
    spectator = VERTEX_SPECTATOR[case]
    at_origin = {E: 0, a: 0, x: 0, y: 0}
    product = d[spectator].subs(at_origin)*d[4].subs(at_origin)**2
    assert product == (-R(1, 8) if case == 0 else -2)
    # Both leading residues enter with -2/(pi*product). The angular
    # boundary value of 1/(mu+i0) over [-1,1] is -i*pi.
    return complex(-(2-sp.I*sp.pi)/product)


def check_mixed_null_endpoints():
    for case, (vertex, _) in enumerate(VERTICES):
        exact = local_centers(finite_frame(), case)
        centers = np.asarray(exact, dtype=float)
        certify_deformation(exact, case)
        predicted = endpoint_target(exact, case)
        mu, eps = np.array(-0.37), MAX_ETA*RADIUS**2
        point = np.array([RADIUS*np.sqrt(1-mu*mu), 0, RADIUS*mu])
        reflected = point*np.array([-1, 1, 1])
        independent = sum(upper_residues(q, centers, eps).sum()
                          for q in (point, reflected))/pi
        np.testing.assert_allclose(
            energy_density(np.array(RADIUS), mu, 1, eps, centers, case),
            independent, rtol=1e-9, atol=1e-8,
        )

        errors = []
        for radius in (RADIUS, 1/128, 1/256, 1/512, 1/1024):
            value = integrate_ball(radius, 0, centers, case, order=20)
            errors.append(abs(value-predicted))
            print(f"D{vertex+1}, R={radius:.6f}: integral/R^2={value:.9f}")
        np.testing.assert_allclose(
            integrate_ball(radius, 0, centers, case, order=28), value, rtol=0, atol=1e-7,
        )
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 0.01*abs(predicted)
        print(f"D{vertex+1}: predicted coefficient={predicted:.9f}")

        boundary = integrate_ball(RADIUS, 0, centers, case, order=28)
        errors = []
        for eta in (0.03, 0.01, 0.003):
            value = integrate_ball(RADIUS, eta, centers, case, order=40)
            errors.append(abs(value-boundary))
            print(f"D{vertex+1}, eta={eta:g}: regulator difference={errors[-1]:.3e}")
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 0.03*abs(predicted)
    print("PASS: both mixed null endpoints and their complex ball coefficients")


if __name__ == "__main__":
    check_mixed_null_endpoints()
