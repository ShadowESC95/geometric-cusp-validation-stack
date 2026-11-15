"""Resolve the null-edge endpoint at the fourth external propagator center."""

from math import pi

import numpy as np
import sympy as sp

from .energy_contour import energy_poles, finite_frame, interval, upper_residues
from .null_endpoint import MAX_ETA, RADIUS, integrate_ball, scaled_model
from .spatial_tail import reflected_centers

GROUP_CENTER = 1.58
GROUP_RADIUS = 0.16
GROUP_OFFSETS = GROUP_RADIUS*np.exp(2j*pi*np.arange(48)/48)


def lower_endpoint_energy(r, mu, cos_phi, eps, centers):
    """Normalized energy integral closed below and averaged across x=0."""
    axial = r*mu
    transverse_sq = r*r*(1-mu)*(1+mu)
    x_sq = transverse_sq*cos_phi**2
    small = transverse_sq-1j*eps
    ratio = cos_phi**2 if eps == 0 else x_sq/small
    omega = np.sqrt(r*r-1j*eps)

    # At the lower D4 pole, D1*D3 = 2*(p_perp^2-i*eps).
    d2 = 2-3*omega-axial
    insertion = 1-(21*omega+11*axial)/8
    soft = 2/pi*ratio/(omega*d2*(insertion**2-x_sq))

    dz = axial[..., None]-0.5-centers[:, 3]
    # At the lower D3 pole, D2 = -2*displacement from the axial pole.
    length = 0.5-axial
    omega3 = np.sqrt(length**2+small)
    pole3 = length+small/(omega3+length)
    d = (pole3[..., None]-centers[:, 0])**2-dz**2-small[..., None]
    insertion = d[..., 4]-0.25
    second_pair = -2/pi*ratio*(omega3+length)/(omega3*d[..., 0]*d[..., 3]
                                              *(insertion**2-x_sq))

    regular = np.zeros_like(x_sq, dtype=complex)
    for offset in GROUP_OFFSETS:
        d = (GROUP_CENTER+offset-centers[:, 0])**2-dz**2-small[..., None]
        insertion = d[..., 4]-0.25
        regular += offset/(np.prod(d[..., :4], axis=-1)*(insertion**2-x_sq))
    # Clockwise closure reverses the sign of the reflected numerator.
    return soft + second_pair + 8/pi*x_sq*regular/len(GROUP_OFFSETS)


def certify_poles(centers):
    extended = np.array(reflected_centers(centers), dtype=float)
    upper, lower = energy_poles([0, 0, -0.5], extended, 0)
    shift = RADIUS*(1+np.sqrt(MAX_ETA))
    assert np.all(abs(lower[[0, 1, 4, 5]]-GROUP_CENTER)+shift < GROUP_RADIUS)
    assert np.all(abs(lower[[2, 3]]-GROUP_CENTER)-shift > GROUP_RADIUS)
    assert np.all(abs(upper-GROUP_CENTER)-shift > GROUP_RADIUS)

    # At the D4 endpoint, exclude spectator D2 and both insertions.
    for branch in (upper, lower):
        assert np.all(abs(branch[[1, 4, 5]]+0.5) > 2*shift)
    # At the other lower pole D3, only the upper D2 pole can meet it.
    assert np.all(abs(upper[[0, 2, 3, 4, 5]]-0.5) > 2*shift)
    assert np.all(abs(lower[[0, 1, 3, 4, 5]]-0.5) > 2*shift)


def endpoint_coefficients(centers):
    R = sp.Rational
    E, x, y, z, r, mu = sp.symbols("E x y z r mu", real=True)
    d = [sp.expand(interval((E, x, y, z), center)) for center in centers]
    on_shell = {E: -R(1, 2)+r, z: -R(1, 2)+r*mu}
    pair = ((d[0]-d[3])*(d[2]-d[3])).subs(on_shell)
    assert sp.expand(pair-2*r**2*(1-mu**2)) == 0
    assert sp.expand(d[1]-d[2]) == -2*(E+z)

    d2, insertion = [sp.expand((d[i]-d[3]).subs(on_shell)) for i in (1, 4)]
    assert d2 == 2-3*r-r*mu
    assert insertion == 1-(21*r+11*r*mu)/8-x
    # Retain x^2 in the reflected denominator by setting x^2=r^2*(1-mu^2)
    # times cos(phi)^2. It first matters one order beyond this expansion.
    radial_factor = 2/(d2*insertion.subs(x, 0)**2)
    leading_density = radial_factor.subs(r, 0)
    first_correction = sp.diff(radial_factor, r).subs(r, 0)
    second_pole = [q.subs({E: R(1, 2), x: 0, y: 0, z: -R(1, 2)}) for q in d]
    assert second_pole == [-1, 0, 0, 1, -R(5, 8)]
    second_density = -4/(second_pole[0]*second_pole[3]*second_pole[4]**2)
    leading = sp.integrate(leading_density, (mu, -1, 1))/2
    subleading = sp.integrate(first_correction+second_density, (mu, -1, 1))/3
    assert leading == 1 and subleading == R(1699, 150)
    return float(leading), float(subleading)


def check_lower_null_endpoint():
    exact_centers = finite_frame()
    certify_poles(exact_centers)
    leading, subleading = endpoint_coefficients(exact_centers)
    centers = np.array(exact_centers, dtype=float)

    r, mu, eps = RADIUS, -0.37, MAX_ETA*RADIUS**2
    point = np.array([r*np.sqrt(1-mu**2), 0, -0.5+r*mu])
    reflected = point*np.array([-1, 1, 1])
    upper = sum(upper_residues(p, centers, eps).sum() for p in (point, reflected))/pi
    np.testing.assert_allclose(
        lower_endpoint_energy(np.array(r), np.array(mu), 1, eps, centers),
        upper, rtol=1e-9, atol=1e-8,
    )

    def ball(radius, eta, order):
        return integrate_ball(radius, eta, centers, order, lower_endpoint_energy)/radius**2

    errors = []
    for radius in (RADIUS, 1/128, 1/256, 1/512, 1/1024):
        value = ball(radius, 0, 16)
        correction = (value.real-leading)/radius
        errors.append(abs(correction-subleading))
        print(f"radius={radius:.6f}: integral/radius^2={value.real:.9f}, "
              f"scaled correction={correction:.9f}")
    np.testing.assert_allclose(ball(radius, 0, 24), value, rtol=0, atol=1e-9)
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 0.03

    target = scaled_model(MAX_ETA, 64)/16
    errors = []
    for radius in (RADIUS, 1/128, 1/256, 1/512):
        value = ball(radius, MAX_ETA, 36)
        errors.append(abs(value-target))
        print(f"eta={MAX_ETA:g}, radius={radius:.6f}: integral/radius^2={value:.9f}")
    np.testing.assert_allclose(ball(radius, MAX_ETA, 48), value, rtol=0, atol=5e-7)
    assert all(new < old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < 0.03
    print(f"eta={MAX_ETA:g}: scaled prediction={target:.9f}")
    print("PASS: lower null-edge endpoint and the first two ball coefficients")


if __name__ == "__main__":
    check_lower_null_endpoint()
