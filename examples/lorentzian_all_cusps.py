"""The edge statement at every cusp and in every real kinematic region.

See docs/GATE_6_BOUNDS.md. For each of the four cusps of the light-like
square, each of its two edges and each of the three sign patterns of
(s, t), the jump of the transverse kernel across the collinear invariant is
-2 pi i times the unit-circle mean of the weight on the edge segment, and
zero on the rest of the null line. The two-loop weight is regular at z = -1.
"""

import mpmath as mp
import numpy as np

from geocusp import lorentzian_edges as le

REGIONS = {"--": (-1.0, -1.7), "+-": (1.3, -1.0), "-+": (-1.0, 0.8)}


def check_lorentzian_all_cusps():
    weight = lambda L: mp.exp(mp.mpf("0.3") * L)
    with mp.workdps(18):
        mean = complex(le.circle_mean(weight))
        for region, (s, t) in REGIONS.items():
            p, s_vec, t_vec = le.light_like_square(region)
            assert (s_vec > 0) == (s > 0) and (t_vec > 0) == (t > 0)
            worst_on = worst_off = 0.0
            sheets = set()
            for index in (1, 2, 3, 4):
                exact, local, _ = le.exact_invariants(index, 0.07, -0.05, 0.02, 1.3, region)
                assert max(abs(exact[j] - local[j]) for j in exact) < 1e-12
                cusp = le.Cusp(index, s, t)
                for a in (1, -1):
                    for b in (1, -1):
                        sheets |= set(cusp.sheet_sequence(a * 1e-3 * abs(cusp.D), b * 1e-7 * abs(cusp.D)))
                for edge in ("next", "prev"):
                    on = complex(cusp.edge_jump(weight, edge, 0.3))
                    off = complex(cusp.edge_jump(weight, edge, -0.4))
                    worst_on = max(worst_on, abs(on + 2j * np.pi * mean))
                    worst_off = max(worst_off, abs(off))
            assert worst_on < 1e-8 and worst_off < 1e-8 and max(abs(k) for k in sheets) <= 2
            print(f"(s,t) signs {region}: 8 edges, max |jump + 2 pi i M| = {worst_on:.1e}, "
                  f"max off-segment jump = {worst_off:.1e}, sheets arg z/pi in {sorted(sheets)}")

    with mp.workdps(25):
        point = mp.mpc("0.3", "-4.0")
        gap = abs(le.f2_cover(point, "inner") - le.f2_cover(point, "outer"))
        assert gap < mp.mpf("1e-18")
        print(f"Two continuations of the two-loop weight differ by {float(gap):.1e}: regular at z = -1")
    print("PASS: edge statement at all four cusps in all three kinematic regions")


if __name__ == "__main__":
    check_lorentzian_all_cusps()
