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
    # the bound is e plus the largest of the three zone constants
    assert out["bound"] == pytest.approx(out["e"] + max(out["K_in"], out["K_bl"], out["K_out"]))
    assert out["K_in"] >= out["K_in_cases_3ab"] and out["K_near_cusp"] < 1e-5


def test_the_constants_of_the_default_square():
    out = pb.phase_bound(fc.Contour(region="+-"))
    assert out["e4"] == pytest.approx(0.02619, abs=1e-5)
    assert out["e_sp"] == pytest.approx(0.01261, abs=1e-5)
    assert out["c1_prime"] == pytest.approx(0.01042, abs=1e-5)
    assert out["K_in"] == pytest.approx(0.1744, abs=1e-4)
    assert out["s"] == pytest.approx(out["s0"]) == pytest.approx(0.2090, abs=1e-4)
    assert out["K_bl"] == pytest.approx(0.2847, abs=1e-4)
    assert (out["s1"], out["s2"], out["s3"]) == pytest.approx((0.0675, 0.0300, 0.00057), abs=1e-4)


def test_the_bound_is_not_uniform_in_the_square_but_stays_below_the_corner():
    out = pb.phase_bound(fc.tuned_contour("+-", 2.5))
    assert out["all_conditions"]
    assert out["bound"] == pytest.approx(1.578, abs=2e-3)
    assert out["K_in"] == out["K_in_case_3d"]            # c_1' grows when the other diagonal is small
    assert out["bound"] < np.pi - fc.PHASE_CORNER_MARGIN


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
    limits = (out["K_in"], out["K_bl"], out["K_out"])
    largest, count, phase = [-np.inf] * 3, [0] * 3, 0.0
    for x0 in _points(C, 1500, seed=11):
        if C.distance_to_edges(x0) < 1e-5:
            continue
        val = pb.k_values(C, x0)
        # float64: skip points where an invariant is within rounding error of zero
        if val["min_modulus"] < 1e-7 * (1 + float((x0 - C.centre) @ (x0 - C.centre))):
            continue
        z = val["zone"]
        count[z] += 1
        largest[z] = max(largest[z], max(val["K"].values()))
        phase = max(phase, abs(val["arg_z"]))
    assert min(count) > 2000, count
    for z in range(3):
        assert largest[z] <= limits[z], (z, largest[z], limits[z])
    assert phase <= 2 * np.pi + out["bound"]
    # what is observed is far below what is proved
    assert max(largest) < 0.01 and phase < 2 * np.pi + 1e-6
