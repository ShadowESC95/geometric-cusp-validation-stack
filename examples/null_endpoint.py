"""Resolve the null-edge endpoint at the second external propagator center."""

from math import pi

import numpy as np
import sympy as sp

from .energy_contour import energy_poles, finite_frame, finite_numerator, interval, upper_residues
from .spatial_tail import reflected_centers

RADIUS = 1/64
MAX_ETA = 0.1
ENERGY_GROUPS = ((-1, (2, 3)), (0.28, (4, 5)))
CONTOUR_RADIUS = 1/8
CONTOUR_OFFSETS = CONTOUR_RADIUS*np.exp(2j*pi*np.arange(32)/32)


def endpoint_energy(r, mu, cos_phi, eps, centers):
    """Normalized energy integral averaged at opposite transverse x."""
    axial = r*mu
    angular_x_sq = cos_phi**2
    transverse_sq = r*r*(1-mu)*(1+mu)
    x_sq = transverse_sq*angular_x_sq
    omega = np.sqrt(r*r-1j*eps)
    # Evaluate the spectator factors on the upper D2 pole using D2=0.
    d4 = 2-3*omega+axial
    insertion = -0.25-3*(omega+axial)/8
    # The null pair product is 2*(p_perp^2-i*eps). At eps=0 the angular
    # cancellation is exact, including along either limiting axial direction.
    ratio = angular_x_sq if eps == 0 else x_sq/(transverse_sq-1j*eps)
    soft = 2/pi*ratio/(omega*d4*(insertion**2-x_sq))

    dz = axial[..., None]-1-centers[:, 3]
    # At a different energy the D1-D4 null pair also projects onto this ball.
    # On the upper D1 pole D4=-2*delta, with
    # delta=(p_perp^2-i*eps)/(omega1+length). Cancel that factor first.
    length = 0.5+axial
    small = transverse_sq-1j*eps
    omega1 = np.sqrt(length**2+small)
    pole1 = -axial-small/(omega1+length)
    d = (pole1[..., None]-centers[:, 0])**2-dz**2-small[..., None]
    insertion = d[..., 4]-0.25
    second_pair = -2/pi*ratio*(omega1+length)/(omega1*d[..., 1]*d[..., 2]
                                              *(insertion**2-x_sq))

    regular = np.zeros(np.broadcast_shapes(x_sq.shape, r.shape), dtype=complex)
    for center, _ in ENERGY_GROUPS:
        for offset in CONTOUR_OFFSETS:
            d = (center+offset-centers[:, 0])**2-dz**2-transverse_sq[..., None]+1j*eps
            insertion = d[..., 4]-0.25
            regular += offset/(np.prod(d[..., :4], axis=-1)*(insertion**2-x_sq))
    return soft + second_pair - 8/pi*x_sq*regular/len(CONTOUR_OFFSETS)


def integrate_ball(radius, eta, centers, order):
    """Integrate over a spatial ball with eps=eta*radius^2."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes[:, None, None]+1)/2
    mu = nodes[None, :, None]
    phi = (pi/4)*(nodes[None, None, :]+1)
    density = endpoint_energy(radius*t, mu, np.cos(phi), eta*radius**2, centers)
    return pi*radius**3/2*np.einsum("i,j,k,ijk->", weights, weights, weights, t*t*density)


def scaled_model(eta, order):
    """Integral of the leading endpoint density in the unit ball."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes[:, None]+1)/2
    transverse_sq = t*t*(1-nodes**2)
    ratio = transverse_sq/(transverse_sq-1j*eta)
    # The azimuth integral of cos(phi)^2 is pi.
    density = t*t*ratio/np.sqrt(t*t-1j*eta)
    return 8*(weights @ density @ weights)


def certify_poles(centers):
    extended = np.array(reflected_centers(centers), dtype=float)
    upper, lower = energy_poles([0, 0, -1], extended, 0)
    # Each pole moves by at most RADIUS+sqrt(eps) anywhere in the ball.
    shift = RADIUS*(1+np.sqrt(MAX_ETA))
    for center, included in ENERGY_GROUPS:
        for branches, is_upper in ((upper, True), (lower, False)):
            for index, pole in enumerate(branches):
                if is_upper and index in included:
                    assert abs(pole-center)+shift < CONTOUR_RADIUS
                else:
                    assert abs(pole-center)-shift > CONTOUR_RADIUS
    assert 2-RADIUS-3*shift > 0
    assert 0.25-3*(shift+RADIUS)/8-RADIUS > 0
    for branches, isolated in ((upper, 0), (lower, 3)):
        for index, pole in enumerate(branches):
            if index != isolated:
                assert abs(pole) > 2*shift
    assert 0.5 > shift
    # The remaining factors at the upper D2 pole have uniform margins,
    # and |omega| >= r gives an integrable O(1/r) bound in three dimensions.
    for branches in (upper, lower):
        assert np.all(abs(branches[[3, 4, 5]]-1) > 2*shift)


def check_null_endpoint():
    exact_centers = finite_frame()
    certify_poles(exact_centers)
    e, x, y, z = sp.symbols("E x y z", real=True)
    point = (e, x, y, z)
    assert sp.expand(finite_numerator(point, exact_centers)) == -4*x
    d = [sp.expand(interval(point, center)) for center in exact_centers]
    assert sp.expand(d[0]-d[1]) == e-z-2
    assert sp.expand(d[2]-d[1]) == 2*e+2*z
    assert sp.expand(d[3]-d[1]) == 3*e+z
    assert sp.expand(d[3]-d[0]) == 2*(e+z+1)
    assert sp.expand(d[4]-d[1]) == 3*e/8-x-3*z/8-1
    omega = sp.symbols("omega", real=True)
    pair = ((d[0]-d[1])*(d[2]-d[1])).subs(e, 1-omega)
    assert sp.expand(pair) == 2*(omega**2-(z+1)**2).expand()
    assert [f.subs({e: 1, x: 0, y: 0, z: -1}) for f in d] == [
        0, 0, 0, 2, -sp.Rational(1, 4),
    ]

    # At eps=0 the soft normalized energy density is
    # (16/pi)*cos(phi)^2/r to leading order. The second pair is O(1)
    # and the regular groups are O(r^2).
    s, phi = sp.symbols("s phi", real=True)
    leading = 16/sp.pi*sp.integrate(s, (s, 0, 1))*2*sp.integrate(
        sp.cos(phi)**2, (phi, 0, 2*sp.pi),
    )
    assert leading == 16
    centers = np.array(exact_centers, dtype=float)

    r, mu, eps = RADIUS, 0.37, MAX_ETA*RADIUS**2
    point = np.array([r*np.sqrt(1-mu**2), 0, -1+r*mu])
    reflected = point*np.array([-1, 1, 1])
    residue_sum = sum(upper_residues(p, centers, eps).sum() for p in (point, reflected))/pi
    np.testing.assert_allclose(
        endpoint_energy(np.array(r), np.array(mu), 1, eps, centers),
        residue_sum, rtol=1e-10, atol=1e-9,
    )

    for eta in (0, MAX_ETA):
        order = 24 if eta == 0 else 36
        target = 16 if eta == 0 else scaled_model(eta, 64)
        if eta:
            np.testing.assert_allclose(scaled_model(eta, 48), target, rtol=0, atol=1e-8)
        errors = []
        for radius in (RADIUS, 1/128, 1/256, 1/512):
            result = integrate_ball(radius, eta, centers, order)/radius**2
            errors.append(abs(result-target))
            print(f"eta={eta:g}, radius={radius:.6f}: integral/radius^2={result:.9f}")
        np.testing.assert_allclose(
            integrate_ball(radius, eta, centers, order+12)/radius**2, result,
            rtol=0, atol=5e-7,
        )
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 0.04
        print(f"eta={eta:g}: scaled prediction={target:.9f}")
    boundary = integrate_ball(RADIUS, 0, centers, 24)/RADIUS**2
    errors = []
    for eta in (0.03, 0.01, 0.003):
        value = integrate_ball(RADIUS, eta, centers, 48)/RADIUS**2
        errors.append(abs(value-boundary))
        print(f"fixed radius, eta={eta:g}: difference from zero regulator={errors[-1]:.3e}")
    assert errors[2] < errors[1] < errors[0] and errors[-1] < 0.8
    print("PASS: external null-edge endpoint and shrinking-ball limit")


if __name__ == "__main__":
    check_null_endpoint()
