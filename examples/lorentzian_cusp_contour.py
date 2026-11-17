"""The Lorentzian chain at a light-like edge and its contour in the z plane.

See docs/GATE_6_CLOSURE.md. Near the null edge through x2 and x3 the
transverse integral of the frozen-loop integrand is a one-dimensional
integral over x02^2. Tracking the sheet of z fixed by the Feynman
prescription, the coefficient of the collinear logarithm is the integral of
F(z) dz/z over a path that deforms to the unit circle on the edge segment and
to a point off it.
"""

import mpmath as mp
import numpy as np

from geocusp import lorentzian_cusp as lc


def check_lorentzian_cusp_contour():
    assert all(r == 0 for r in lc.local_identities().values())
    brackets = lc.weighted_chart_dictionary()
    print("Weighted chart: <AB12>=%s, <AB23>=%s, <AB14>=%s, 1+z=%s"
          % (brackets["12"], brackets["23"], brackets["14"], brackets["one_plus_z"]))

    expected = {(1, 1): [-1, -2, -1], (-1, -1): [1, 0], (-1, 1): [-1, 0], (1, -1): [-1]}
    for (s1, s3), seq in expected.items():
        _, phi = lc.sheet_path(s1 * 1e-3, s3 * 1e-7)
        assert lc.plateaus(phi) == seq
        print(f"sign(x01^2)={s1:+d}, sign(x03^2)={s3:+d}: arg z/pi visits {seq}")

    weights = {
        "z^0.3": lambda L: mp.exp(mp.mpf("0.3") * L),
        "log^2 z + pi^2": lambda L: L**2 + mp.pi**2,
    }
    ratio = mp.mpf("1e-12")
    with mp.workdps(20):
        for name, f in weights.items():
            mean = complex(lc.circle_mean(f))
            for beta in (0.25, 0.75):
                q1 = -mp.mpf(beta)
                jump = complex(lc.edge_jump(f, q1, abs(q1) * ratio))
                np.testing.assert_allclose(jump, -2j * np.pi * mean, rtol=0, atol=1e-7 * abs(mean))
            outside = complex(lc.edge_jump(f, mp.mpf("0.3"), mp.mpf("0.3") * ratio))
            assert abs(outside) < 1e-7 * abs(mean)
            one = lambda L: 1
            cusp = complex(lc.quadrant_sum(f, 1e-3, 1e-15) / lc.quadrant_sum(one, 1e-3, 1e-15))
            np.testing.assert_allclose(cusp, mean, rtol=0, atol=1e-7 * abs(mean))
            print(f"{name}: edge jump/(-2 pi i)={(jump/(-2j*np.pi)).real:.10f}, "
                  f"cusp quadrant ratio={cusp.real:.10f}, circle mean={mean.real:.10f}, "
                  f"off-segment jump={abs(outside):.1e}")
    print("PASS: Lorentzian edge chain deforms to the unit circle on the weighted divisor")


if __name__ == "__main__":
    check_lorentzian_cusp_contour()
