"""Signed intersection balance for projective CD contour families."""

from math import pi, sqrt

import numpy as np
import sympy as sp

from geocusp.contours import winding_number


def certify_divisor_balance():
    X0, X1, Y0, Y1, z = sp.symbols("X0 X1 Y0 Y1 z", nonzero=True)
    section = z*X1*(Y0-Y1)+Y1*(X0-X1)
    affine = section/(X0*Y0)
    X, Y = sp.symbols("X Y")
    assert sp.expand(section.subs({X0: 1, Y0: 1, X1: X, Y1: Y})) == (
        -(1+z)*X*Y+z*X+Y
    ).expand()
    assert affine.subs({X0: 1, X1: 1, Y0: 1, Y1: 0}) == z
    assert affine.subs({X0: 1, X1: 0, Y0: 1, Y1: 1}) == 1

    # The ratio extends to a single meromorphic function, independent
    # of homogeneous representatives, so the balance does not depend on charts.
    alpha, beta = sp.symbols("alpha beta", nonzero=True)
    scaled = affine.subs({X0: alpha*X0, X1: alpha*X1,
                          Y0: beta*Y0, Y1: beta*Y1}, simultaneous=True)
    assert sp.cancel(scaled-affine) == 0

    # In x,y coordinates the same function has a zero and two possible poles.
    x, y = sp.symbols("x y")
    converted = affine.subs({X0: 1+x, X1: x, Y0: 1+y, Y1: y})
    assert sp.cancel(converted-(z*x+y)/((1+x)*(1+y))) == 0


def trace_rotating_family(m, n):
    """Count local crossings of the three divisors on an explicit cylinder."""
    # Use unequal scales to keep the three sets of intersections distinct.
    x_scale, y_scale = 2., 3.
    v_x = x_scale/(1+x_scale)
    v_y = 1/(1+y_scale)
    v_section = sqrt(x_scale)/(sqrt(x_scale)+sqrt(y_scale))
    v_values = (v_section, v_x, v_y)
    assert len(set(v_values)) == 3
    exponents = (m+1-n, m, n)

    def factors(v, phase):
        z = np.exp(1j*phase)
        x = x_scale*z**m*(1-v)/v
        y = y_scale*z**n*v/(1-v)
        return np.array([z*x+y, 1+x, 1+y])

    indices = np.zeros(3, dtype=int)
    theta = 2*pi*np.arange(256)/256
    radius = min(0.025, min(abs(a-b) for a in v_values for b in v_values if a != b)/4)
    for factor_index, (v_center, exponent) in enumerate(zip(v_values, exponents)):
        # At the listed v there is one crossing per phase satisfying z^k=-1.
        for angle in (2*np.arange(abs(exponent))+1)*pi/max(1, abs(exponent)):
            at_center = factors(v_center, angle)
            assert abs(at_center[factor_index]) < 1e-12
            assert np.all(abs(np.delete(at_center, factor_index)) > .05)
            curves = factors(v_center+radius*np.cos(theta), angle+radius*np.sin(theta))
            local = np.array([winding_number(curve) for curve in curves])
            nonzero = np.flatnonzero(abs(local) > .5)
            np.testing.assert_array_equal(nonzero, [factor_index])
            np.testing.assert_allclose(local, np.rint(local), atol=1e-12, rtol=0)

            quotient = curves[0]/(curves[1]*curves[2])
            np.testing.assert_allclose(winding_number(quotient),
                                       local[0]-local[1]-local[2], atol=1e-12, rtol=0)
            indices += np.rint(local).astype(int)
    return indices


def check_cd_projective_intersections():
    certify_divisor_balance()
    theta = 2*pi*np.arange(256)/256
    corner_windings = (winding_number(np.exp(1j*theta)), winding_number(np.ones_like(theta)))
    assert corner_windings == (1, 0)

    # (0,1) rotates y alone. (-1,0) routes the crossing through X=infinity.
    # The other cases include either a moving-divisor crossing or extra
    # infinity crossings of cancelling signs.
    for m, n in ((0, 1), (-1, 0), (0, 0), (1, 2), (-2, -1), (2, 0)):
        actual = trace_rotating_family(m, n)
        expected = np.array([n-m-1, -m, n])
        np.testing.assert_array_equal(actual, expected)
        net = actual[0]-actual[1]-actual[2]
        assert net == corner_windings[1]-corner_windings[0]
        print(f"(m,n)=({m:+d},{n:+d}): intersections {tuple(int(i) for i in actual)}, "
              f"zero-minus-poles={net:+d}")

    print("PASS: compactified CD intersection balance for both infinity charts")


if __name__ == "__main__":
    check_cd_projective_intersections()
