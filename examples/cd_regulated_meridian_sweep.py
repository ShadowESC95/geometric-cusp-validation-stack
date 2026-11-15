"""Sweep the lower-corner regulated cap roots over the auxiliary triangle."""

from math import pi, sin

import numpy as np
from scipy.optimize import brentq

from geocusp.contours import winding_number

from .cd_regulated_base_cap import denominator_on_cap
from .cd_regulated_cap_meridians import RADIUS, limiting_crossing
from .full_finite_frame import ROTATION


CENTER = 5*RADIUS/2


def root_direction(s):
    """Choose the continuous root direction with negative imaginary part."""
    H = 1-(1+ROTATION)*s
    return -1j*RADIUS*np.sqrt(1j*H)


def swept_roots(scale, A, s, eps):
    """Use scale=x/(eps*A), extending the root formula to A=0."""
    x = eps*A*scale
    G = 1+(ROTATION-1)*s
    center = CENTER+RADIUS/2*(4*x+eps*scale*(G-A)/(1+x))
    offset = root_direction(s)*np.sqrt(scale)
    return np.stack((center+offset, center-offset))


def check_crossing_surfaces():
    s_equal = 1/(1+ROTATION.real)
    s = np.sort(np.r_[np.linspace(0., 1., 129), s_equal])
    direction = root_direction(s)
    np.testing.assert_allclose(direction**2, -1j*RADIUS**2*(1-(1+ROTATION)*s),
                               atol=2e-18, rtol=0)
    assert np.max(direction.imag) < 0
    thresholds = limiting_crossing(CENTER, np.stack((direction, -direction)))
    gap = sin(pi/8)
    assert np.min(thresholds) >= ((1-CENTER)/RADIUS)**2
    assert np.max(thresholds) <= ((1+CENTER)/RADIUS)**2/gap

    equality = np.flatnonzero(s == s_equal)[0]
    difference = thresholds[0]-thresholds[1]
    assert np.all(difference[:equality] < 0) and np.all(difference[equality+1:] > 0)
    equal_threshold = (1-CENTER**2)/(RADIUS**2*ROTATION.imag*s_equal)
    np.testing.assert_allclose(thresholds[:, equality], equal_threshold, atol=1e-12)
    roots = swept_roots(equal_threshold, .1, s_equal, 0.)
    np.testing.assert_allclose(roots, CENTER+np.array([-1j, 1j])*np.sqrt(1-CENTER**2),
                               atol=1e-14, rtol=0)
    # Differentiate the two threshold equations at their intersection.
    # The common length derivative cancels from their slope difference.
    at_equal = direction[equality]
    derivative = 1j*RADIUS**2*(1+ROTATION)/(2*at_equal)
    crossing_slope = -4*CENTER*derivative.real*np.sqrt(equal_threshold)/abs(at_equal)**2
    width = 1e-6
    near = root_direction(s_equal+np.array([-width, width]))
    observed = limiting_crossing(CENTER, near)-limiting_crossing(CENTER, -near)
    np.testing.assert_allclose((observed[1]-observed[0])/(2*width), crossing_slope,
                               rtol=1e-8, atol=0)
    assert crossing_slope > 0
    print(f"Rim crossings coincide in scale at s={s_equal:.9f}, "
          f"x/(eps*A)={equal_threshold:.9f}; the roots are distinct")

    # The analytic expression extends to the A=0 edge. Since x=eps*A*scale,
    # evaluate the original normalized denominator only for A>0.
    previous = np.inf
    for eps in (1e-5, 1e-6, 1e-7):
        errors = []
        for s_value in (0., .25, s_equal, .85, .98):
            slope = root_direction(s_value)
            limits = limiting_crossing(CENTER, np.array([slope, -slope]))
            for A in (0., .2*(1-s_value), .8*(1-s_value)):
                for branch, limit in enumerate(limits):
                    def rim_error(scale):
                        return abs(swept_roots(scale, A, s_value, eps)[branch])-1

                    scale = brentq(rim_error, .97*limit, 1.03*limit)
                    errors.append(abs(scale-limit))
                    if A > 0:
                        k = swept_roots(scale, A, s_value, eps)[branch]
                        denominator = denominator_on_cap(k, RADIUS, eps*A*scale, A, s_value, eps)
                        np.testing.assert_allclose(denominator, 0., atol=2e-14, rtol=0)
        error = max(errors)
        assert error < previous
        previous = error
        print(f"eps={eps:g}: maximum sampled crossing error={error:.3e}")
    assert previous < .1


def check_swept_meridians():
    unit = np.exp(2j*pi*np.arange(512)/512)
    punctures = np.array([0., RADIUS, CENTER])
    counts = set()
    for s in (.02, .25, .75, .95):
        slope = root_direction(s)
        limits = np.sort(limiting_crossing(CENTER, np.array([slope, -slope])))
        scales = (.5*limits[0], np.mean(limits), 1.2*limits[1])
        for fraction in (.05, .7):
            # Use a fraction of the available triangle width.
            auxiliary = fraction*(1-s)
            for scale in scales:
                roots = swept_roots(scale, auxiliary, s, 0.)
                inside = abs(roots) < 1.
                expected = np.count_nonzero(inside)
                counts.add(expected)

                def index(curve):
                    curve = np.asarray(curve)
                    polynomial = (curve-roots[0])*(curve-roots[1])
                    return winding_number(polynomial)

                holes = np.minimum(.01, np.min(abs(punctures[:, None]-roots), axis=1)/4)
                np.testing.assert_allclose([index(p+r*unit) for p, r in zip(punctures, holes)],
                                           0, atol=1e-12, rtol=0)
                indices = []
                for root_index in np.flatnonzero(inside):
                    root = roots[root_index]
                    clearance = min(1-abs(root), abs(root-roots[1-root_index]),
                                    np.min(abs(root-punctures)-holes))
                    assert clearance > 0
                    indices.append(index(root+clearance/4*unit))
                np.testing.assert_allclose(indices, 1, atol=1e-12, rtol=0)
                np.testing.assert_allclose(index(unit), expected, atol=1e-12, rtol=0)
                np.testing.assert_allclose(index(unit)-sum(indices), 0, atol=1e-12, rtol=0)
    assert counts == {0, 1, 2}
    print("All three winding counts balance on the sampled triangle slices")


def check_cd_regulated_meridian_sweep():
    check_crossing_surfaces()
    check_swept_meridians()
    print("PASS: lower-corner rim surfaces and the regulated meridian sweep")


if __name__ == "__main__":
    check_cd_regulated_meridian_sweep()
