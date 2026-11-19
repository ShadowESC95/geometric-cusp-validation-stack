"""Tests for the phase bound with a timelike diagonal (docs/GATE_6_PHASE.md)."""
import numpy as np
import pytest

from geocusp import feynman_contour as fc

from geocusp import phase_bound as pb


SQUARES = {
    "default +-": (lambda: fc.Contour(region="+-"), 0.359),
    "default -+": (lambda: fc.Contour(region="-+"), 0.359),
    "tuned +- 0.7": (lambda: fc.tuned_contour("+-", 0.7), 0.534),
    "tuned -+ 1.1": (lambda: fc.tuned_contour("-+", 1.1), 0.537),
}
CONDITIONS = ["certificate", "R4_centres", "R4_trivial", "G1", "G2", "R5_trivial",
              "B1", "B2", "B3", "B4", "B5", "R6_trivial"]


@pytest.mark.parametrize("name", list(SQUARES))
def test_conditions_hold_and_the_bound_is_below_the_corner(name):
    make, value = SQUARES[name]
    out = pb.phase_bound(make())
    failed = [k for k in CONDITIONS if not out[k]]
    assert not failed, failed
    assert out["all_conditions"]
    assert out["bound"] == pytest.approx(value, abs=5e-4)
    assert out["bound"] < np.pi - fc.PHASE_CORNER_MARGIN
    # the bound is the largest of the three zone bounds; the blend zone decides
    assert out["bound"] == max(out["bound_inner"], out["bound_blend"], out["bound_outer"])
    assert out["bound"] == pytest.approx(out["e"] + out["K_bl"])
    assert out["R7_outer_exact"] and out["bound_outer"] == 0.0
    assert out["bound_inner"] < 0.1
    assert out["K_in"] >= out["K_in_cases_3ab"] and out["K_near_cusp"] < 1e-5


def test_the_constants_of_the_default_square():
    out = pb.phase_bound(fc.Contour(region="+-"))
    assert out["e4"] == pytest.approx(0.02619, abs=1e-5)
    assert out["e_sp"] == pytest.approx(0.01261, abs=1e-5)
    assert out["c_3d"] == pytest.approx(0.00288, abs=1e-5)
    assert out["bound_inner"] == pytest.approx(0.0459, abs=2e-4)
    assert out["K_in"] == pytest.approx(0.1744, abs=1e-4)
    assert out["s"] == pytest.approx(out["s0"]) == pytest.approx(0.2090, abs=1e-4)
    assert out["K_bl"] == pytest.approx(0.2847, abs=1e-4)
    assert (out["s1"], out["s2"], out["s3"]) == pytest.approx((0.0675, 0.0300, 0.00057), abs=1e-4)


@pytest.mark.parametrize("region", ["+-", "-+"])
def test_the_bound_is_uniform_over_the_family_of_squares(region):
    values = []
    for angle in (0.02, 0.5, 1.5, 2.2, 2.5, 2.9, 3.087):
        out = pb.phase_bound(fc.tuned_contour(region, angle))
        assert out["all_conditions"], angle
        values.append(out["bound"])
    assert 0.53 < min(values) and max(values) < 0.57


def test_no_timelike_diagonal_is_rejected():
    with pytest.raises(ValueError):
        pb.phase_bound(fc.Contour(region="--"))


def _transverse_basis(C, j):
    pm, pp = C.p[(j - 1) % 4], C.p[j]
    D = 2 * fc.mdot(pm, pp)
    basis = []
    for seed in np.eye(4):
        v = seed - 2 * fc.mdot(seed, pp) / D * pm - 2 * fc.mdot(seed, pm) / D * pp
        for e in basis:
            v = v + fc.mdot(v, e) * e
        n = -fc.mdot(v, v)
        if n > 1e-12:
            basis.append(v / np.sqrt(n))
        if len(basis) == 2:
            return basis


def _points(C, n, seed):
    """Real points in all three zones: the stress sample of the contour, points near the null lines
    through adjacent cusps (on and beyond the segments, offsets 1e-9 ... 1), points near one cone at
    radii up to 100, and generic directions."""
    rng = np.random.default_rng(seed)
    yield from fc.sample_points(C, n, seed)
    bases = [_transverse_basis(C, j) for j in range(4)]
    for _ in range(2 * n):
        j = int(rng.integers(0, 4))
        big = rng.choice([-1, 1]) * 10 ** rng.uniform(-3, 1.5)
        small = rng.choice([-1, 1]) * 10 ** rng.uniform(-9, 0)
        a, b = (big, small) if rng.random() < 0.5 else (small, big)
        th = rng.uniform(0, 2 * np.pi)
        yield (C.x[j] + a * C.p[(j - 1) % 4] + b * C.p[j]
               + 10 ** rng.uniform(-9, 0) * (np.cos(th) * bases[j][0] + np.sin(th) * bases[j][1]))
    for _ in range(n):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        lam = 10 ** rng.uniform(0.5, 2)
        yield (C.x[int(rng.integers(0, 4))] + lam * np.array([rng.choice([-1.0, 1.0]), *v])
               + rng.normal(size=4) * 10 ** rng.uniform(-6, 0))
    for _ in range(n):
        v = rng.normal(size=4)
        yield C.centre + v / np.linalg.norm(v) * C.radius * rng.uniform(0.5, 6)
    for _ in range(n):                               # the blend zone: generic, and along the null lines
        v = rng.normal(size=4)
        yield C.centre + v / np.linalg.norm(v) * C.radius * rng.uniform(1.0, 2.0)
        j = int(rng.integers(0, 4))
        big = rng.choice([-1, 1]) * rng.uniform(4, 13)
        small = rng.choice([-1, 1]) * 10 ** rng.uniform(-9, 0)
        a, b = (big, small) if rng.random() < 0.5 else (small, big)
        yield C.x[j] + a * C.p[(j - 1) % 4] + b * C.p[j] + 10 ** rng.uniform(-9, -1) * bases[j][0]


@pytest.mark.parametrize("name", ["default +-", "tuned -+ 1.1"])
def test_K_and_the_phase_stay_below_the_proved_constants_on_samples(name):
    C = SQUARES[name][0]()
    out = pb.phase_bound(C)
    limits = (out["K_in_sharp"], out["K_bl"], out["K_out"])
    zone_bounds = (out["bound_inner"], out["bound_blend"], out["bound_outer"])
    largest, count, phase = [-np.inf] * 3, [0] * 3, [0.0] * 3
    for x0 in _points(C, 500, seed=11):
        if C.distance_to_edges(x0) < 1e-5:
            continue
        val = pb.k_values(C, x0)
        # float64: skip points where an invariant is within rounding error of zero
        if val["min_modulus"] < 1e-7 * (1 + float((x0 - C.centre) @ (x0 - C.centre))):
            continue
        z = val["zone"]
        count[z] += 1
        largest[z] = max(largest[z], max(val["K"].values()))
        phase[z] = max(phase[z], abs(val["arg_z"]))
    assert min(count) > 650, count
    for z in range(3):
        assert largest[z] <= limits[z], (z, largest[z], limits[z])
        assert phase[z] <= 2 * np.pi + zone_bounds[z] + 1e-12, (z, phase[z])
    # what is observed is far below what is proved
    assert max(largest) < 0.01 and max(phase) < 2 * np.pi + 1e-6


# ---- interval arithmetic (phase_bound_iv)
from geocusp import phase_bound_iv as pbi


@pytest.mark.parametrize("region,lo,hi,boxes", [("+-", 1.1, 1.103, 12), ("-+", 2.236, 2.239, 12),
                                                 ("+-", 0.02, 0.023, 12), ("-+", 3.084, 3.087, 12)])
def test_interval_conditions_hold_on_ranges_of_squares(region, lo, hi, boxes):
    result = pbi.cover(region, lo, hi, boxes)
    assert result["R7_outer_exact"]
    assert result["all"] and result["failed_count"] == 0
    assert result["bound_upper"] < 0.58
    assert result["Lambda_lower"] >= 2


def test_interval_bound_dominates_the_floating_point_bound():
    for region, angle in [("+-", 1.1), ("-+", 2.5)]:
        r = pbi.interval_phase_bound(region, angle - 1e-4, angle + 1e-4)
        o = pb.phase_bound(fc.tuned_contour(region, angle))
        assert r["all_conditions"]
        assert o["bound"] <= r["bound_upper"] <= o["bound"] + 0.02
        assert o["e"] <= r["e_upper"] and o["K_in"] <= r["K_in_upper"] and o["K_bl"] <= r["K_bl_upper"]


def test_interval_version_rejects_the_region_without_a_timelike_diagonal():
    with pytest.raises(ValueError):
        pbi.interval_phase_bound("--", 1.0, 1.1)
