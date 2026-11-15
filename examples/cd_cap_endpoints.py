"""Resolve the endpoint faces of the swept CD cap in two ratio charts."""

from math import pi

import numpy as np
import sympy as sp

from .cd_resolved_seam import normal_coordinates


def certify_endpoint_forms():
    x, y, z, a, b = sp.symbols("x y z a b")
    kernel = 2*(1+z)/(z*(1+x)*(1+y)*(z*x+y))
    reciprocal = kernel.subs({x: 1/x, y: 1/y, z: 1/z}, simultaneous=True)
    assert sp.cancel(reciprocal/(x*y*z)**2-kernel) == 0  # signed Jacobian is negative

    # First blow-up chart: y=x*a, with coordinate order x,a,z.
    ratio_form = sp.cancel(x*kernel.subs(y, x*a))
    expected = 2*(1+z)/(z*(1+x)*(1+x*a)*(z+a))
    assert sp.cancel(ratio_form-expected) == 0
    assert sp.residue(ratio_form, x, 0) == 0

    # Second chart: x=y*b, with coordinate order b,y,z.
    other_form = sp.cancel(y*kernel.subs(x, y*b))
    assert sp.cancel(other_form-2*(1+z)/(z*(1+y*b)*(1+y)*(1+z*b))) == 0
    assert sp.residue(other_form, y, 0) == 0
    assert (1+z*b).subs(b, 0) == 1  # the two legs join at y=b=0

    # At the upper endpoint use X=1/x, Y=1/y, keeping z projective.
    # The form in (X,Y,z) has moving factor z*Y+X instead of z*X+Y.
    upper = kernel.subs({x: 1/x, y: 1/y}, simultaneous=True)/(x*y)**2
    upper_ratio = sp.cancel(x*upper.subs(y, x*a))
    assert sp.cancel(upper_ratio-2*(1+z)/(z*(1+x)*(1+x*a)*(1+z*a))) == 0
    assert sp.residue(upper_ratio, x, 0) == 0

    # The transfer collar keeps z+a constant. In the reciprocal chart its
    # strict moving factor is (z+a)/(z*a).
    s, q, k, delta = sp.symbols("s q k delta")
    w = sp.symbols("w")
    # The join radius is |w|=x^(-1/2). In the second chart, b=x/y.
    y_of_w = (x*w+1)/(w-1)
    assert sp.cancel(y_of_w.subs(w, 1/k)-x-(1+x)*k/(1-k)) == 0
    # Along the log interpolation on the join circle both chart
    # coordinates tend to zero. The rational bounds are symmetric
    # in y and b after setting the endpoint distance to 1/x upstairs.
    assert sp.cancel((x/y_of_w).subs(w, 1/k)-(x*(1-k)/(x+k))) == 0

    ratio_y, base = 1+(1-s)*q, -1+s*q
    assert sp.expand(base+ratio_y) == q
    assert sp.cancel(1/base+1/ratio_y-q/(base*ratio_y)) == 0
    assert sp.expand((delta-k)+k) == delta  # a=1, homogeneous [k:delta-k]


def check_logarithmic_collar(x, normal, phase, delta):
    """Check the split collar, including the two charts on the joining circle."""
    y = normal-1
    q = delta*np.exp(-1j*phase)
    start = np.log(abs(normal)/(1+x))+1j*pi-np.log1p(-normal/(1+x))
    end = np.log((1+x)/(x*delta))+np.log1p(x*q/(1+x))
    split = -.5*np.log(x)
    # Polar radial segments at the start and end phases, with a phase
    # interpolation on the intermediate circle |w|=x^(-1/2).
    for alpha in np.linspace(0, 1, 5):
        logs = ((1-alpha)*start+alpha*(split+1j*start.imag),
                split+1j*((1-alpha)*start.imag+alpha*end.imag),
                (1-alpha)*(split+1j*end.imag)+alpha*end)
        for value in logs:
            w = np.exp(1j*phase+value)
            # Homogeneous inverse y=[w-1 : x*w+1] in the affine order
            # [denominator:numerator]. Both excluded y factors are nonzero.
            y0, y1 = w-1, x*w+1
            np.testing.assert_allclose((y1+y0)/(y1-x*y0), w, rtol=1e-9, atol=1e-12)

    def joins(argument):
        w = np.exp(split+1j*argument)
        if x <= 1:
            endpoint, other = x, (x*w+1)/(w-1)
        else:
            endpoint, other = 1/x, (w-1)/(x*w+1)
        return max(abs(other).max(), abs(endpoint/other).max())

    return joins(phase+start.imag), joins(phase+end.imag)


def check_transfer_and_disk(phase, delta):
    q = delta*np.exp(-1j*phase)
    for s in np.linspace(0, 1, 5):
        a = 1+(1-s)*q
        z = -1+s*q
        np.testing.assert_allclose(z+a, q, atol=1e-15, rtol=0)
        np.testing.assert_allclose(1/z+1/a, q/(z*a), atol=1e-15, rtol=0)
        assert np.min(abs(q/(z*a))) >= delta/(1+delta)**2

    # Verify the two homogeneous charts across the whole base disk via
    # positive normal loops around its two punctures, plus its outer rim.
    unit = np.exp(1j*phase)
    for k in (unit, delta*unit/4, delta+delta*unit/4):
        z0, z1 = k, delta-k
        np.testing.assert_allclose((z0+z1)/(abs(z0)+abs(z1)),
                                   delta/(abs(z0)+abs(z1)), atol=1e-15, rtol=0)
        assert np.min(delta/(abs(z0)+abs(z1))) >= delta/(2+delta)


def check_limiting_collar(phase, normal, delta):
    """Check the two radial legs and their phase path on the limiting faces."""
    y = normal-1
    q = delta*np.exp(-1j*phase)
    # With w=(1+y)/(y-x), use h=1/w at the lower end and h=1/(x*w)
    # at the upper end. Their limits at the seam are y/f and -1/f.
    # The proper endpoint legs are [1-h:h] in y, [1+h:-h] in 1/y.
    # On the exceptional legs use e=1/(a-1) in the relevant ratio a.
    start = (y/normal, -1/normal)
    end = (1/q, -(1+q)/q)
    scale = np.linspace(0, 1, 12)[:, None]
    for h0, e0, sign in zip(start, end, (1, -1)):
        h, e = scale*h0, scale*e0
        proper0, proper1 = 1-sign*h, sign*h
        exception0, exception1 = e, e+1
        np.testing.assert_allclose(proper0+proper1, 1, atol=1e-14, rtol=0)
        np.testing.assert_allclose(exception1-exception0, 1, atol=1e-14, rtol=0)
        # Homogeneous norms are bounded on both radial legs. The moved
        # section and the excluded y=-1 section stay nonzero as indicated
        # by the normalized constants, including at the common join.
        for norm in (abs(proper0)+abs(proper1), abs(exception0)+abs(exception1)):
            assert np.min(1/norm) >= delta/(2+3*delta)


def check_cd_cap_endpoints():
    certify_endpoint_forms()
    phase = 2*pi*np.arange(128)/128
    radius = .06
    normal, _, _ = normal_coordinates(radius, phase)
    # Hold the transferred meridian radius in z fixed on the full x sweep.
    # Within 0<delta<=R the collars and base disk are unchanged in type.
    check_transfer_and_disk(phase, radius)
    check_limiting_collar(phase, normal, radius)
    previous = np.inf
    for radial_distance in (1e-2, 1e-4, 1e-6):
        lower = max(check_logarithmic_collar(radial_distance, normal, phase, radius))
        upper = max(check_logarithmic_collar(1/radial_distance, normal, phase, radius))
        assert max(lower, upper) < previous
        previous = max(lower, upper)
        assert previous <= 1.25*np.sqrt(radial_distance)
        print(f"endpoint distance={radial_distance:g}: "
              f"chart join bounds={lower:.3e}, {upper:.3e}")

    print("PASS: both swept-cap endpoint faces lift off the strict moving divisor")


if __name__ == "__main__":
    check_cd_cap_endpoints()
