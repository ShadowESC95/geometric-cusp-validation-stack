"""Track the regulated meridians through the rim of the diagonal base cap."""

from math import pi, sqrt

import numpy as np
from scipy.optimize import brentq

from geocusp.contours import winding_number

from .cd_regulated_base_cap import cap_root_coefficients, denominator_on_cap
from .full_finite_frame import ROTATION


RADIUS, AUXILIARY_A, AUXILIARY_S = .06, .5, .25


def limiting_crossing(center, slope):
    """Solve |center+slope*sqrt(chi)|=1 in the joint scaling limit."""
    projection = (center.conjugate()*slope).real
    length = abs(slope)**2
    distance = (-projection+sqrt(projection**2+length*(1-abs(center)**2)))/length
    assert distance > 0
    root = center+slope*distance
    np.testing.assert_allclose(abs(root), 1., atol=1e-14, rtol=0)
    assert (root.conjugate()*slope).real > 0
    return distance**2


def check_rim_crossings():
    center, slope = cap_root_coefficients(RADIUS, 0., AUXILIARY_A, AUXILIARY_S)
    slopes = (slope, -slope)
    limits = np.array([limiting_crossing(center, direction) for direction in slopes])
    np.testing.assert_allclose(limits, [189.052122729, 270.946552526], atol=1e-8, rtol=0)
    print(f"Limiting rim crossings x/eps={limits[0]:.9f}, {limits[1]:.9f}")

    previous = np.inf
    for eps in (1e-4, 1e-5, 1e-6):
        crossings = []
        for direction, limit in zip(slopes, limits):
            def root(chi):
                current, _ = cap_root_coefficients(RADIUS, eps*chi, AUXILIARY_A, AUXILIARY_S)
                return current+direction*sqrt(chi)

            chi = brentq(lambda chi: abs(root(chi))-1, .98*limit, 1.02*limit)
            crossings.append(chi)
            k = root(chi)
            denominator = denominator_on_cap(k, RADIUS, eps*chi, AUXILIARY_A, AUXILIARY_S, eps)
            np.testing.assert_allclose(denominator, 0., atol=3e-14, rtol=0)
            # The outward speed for increasing chi is positive at each crossing.
            G = 1+(ROTATION-1)*AUXILIARY_S
            center_derivative = eps*RADIUS/2*(4+(G-AUXILIARY_A)/(AUXILIARY_A*(1+eps*chi)**2))
            velocity = center_derivative+direction/(2*sqrt(chi))
            assert (k.conjugate()*velocity).real > 0
        error = np.max(abs(crossings-limits))
        assert error < previous
        previous = error
        print(f"eps={eps:g}: rim crossings={crossings[0]:.6f}, {crossings[1]:.6f}")
    assert previous < .01


def check_meridian_balance():
    unit = np.exp(2j*pi*np.arange(512)/512)
    holes = np.array([0., RADIUS, 2.5*RADIUS])
    eps = 1e-5
    for chi, count in ((.25, 2), (220., 1), (320., 0)):
        x = eps*chi
        center, slope = cap_root_coefficients(RADIUS, x, AUXILIARY_A, AUXILIARY_S)
        roots = center+np.array([slope, -slope])*sqrt(chi)
        interior = roots[abs(roots) < 1]
        assert len(interior) == count

        def index(curve):
            value = winding_number(denominator_on_cap(curve, RADIUS, x, AUXILIARY_A,
                                                       AUXILIARY_S, eps))
            assert abs(value-round(value)) < 1e-12
            return round(value)

        rim = index(unit)
        original = [index(hole+.02*unit) for hole in holes]
        additions = []
        for root in interior:
            # A smaller circle stays within the rim and away from all other holes.
            width = min(.005, (1-abs(root))/4)
            assert min(abs(root-holes)) > .02+width
            assert abs(roots[0]-roots[1]) > 2*width
            additions.append(index(root+width*unit))
        assert rim == count and original == [0, 0, 0] and additions == [1]*count
        # The clockwise boundary of each excision has the opposite index.
        assert rim-sum(original)-sum(additions) == 0
        print(f"x/eps={chi:g}: rim winding={rim}, regulated meridians={additions}")


def check_cd_regulated_cap_meridians():
    check_rim_crossings()
    check_meridian_balance()
    print("PASS: regulated rim crossings and relative meridian balance")


if __name__ == "__main__":
    check_cd_regulated_cap_meridians()
