"""The certificate in interval arithmetic (GATE_6_CONTOUR section 9)."""
import pytest

from geocusp import certificate_iv as ci
from geocusp import feynman_contour as fc


@pytest.mark.parametrize("region", ["--", "+-", "-+"])
def test_interval_certificate_holds_on_a_range_of_squares(region):
    result = ci.cover(region, 0.9, 1.3, 400)
    assert result["all"] and result["failed_count"] == 0
    assert result["excursion_upper"] < 0.2 and result["mu_Q_upper"] < 0.05


@pytest.mark.parametrize("region", ["--", "+-", "-+"])
def test_interval_certificate_agrees_with_the_floating_point_one(region):
    box = ci.interval_certificate(region, 1.0999, 1.1001)
    point = fc.certificate(fc.tuned_contour(region, 1.1))
    assert box["all"] and point["all"]
    assert point["mu_Q"] <= box["mu_Q_upper"] * 1.0001
    for key, value in point.items():
        if isinstance(value, bool) and key in box:
            assert box[key] == value


def test_interval_comparisons_are_strict():
    from mpmath import iv
    wide = iv.mpf([0.9, 1.1])
    assert not ci.lt(wide, 1) and not ci.lt(1, wide)
    assert ci.lt(wide, 1.2) and ci.le(wide, 1.2)


def test_interval_certificate_rejects_a_bad_contour(monkeypatch):
    """With the far field made too strong the excursion condition must fail."""
    real = fc.tuned_contour

    def strong(region="--", angle=1.1, **kw):
        contour = real(region, angle, **kw)
        contour.kappa_far = 0.2
        return contour

    monkeypatch.setattr(fc, "tuned_contour", strong)
    result = ci.interval_certificate("--", 1.0999, 1.1001)
    assert not result["all"] and not result["D_excursion_small"]


@pytest.mark.parametrize("region", ["--", "+-", "-+"])
def test_field_has_the_right_sign_next_to_the_edges_in_high_precision(region):
    result = ci.edge_sign_check(region, 1.1, dps=100, samples=30, seed=1)
    assert result["pairs"] >= 30 and result["all_negative"]


@pytest.mark.parametrize("region,lo,hi", [("--", 0.06, 0.08), ("+-", 0.02, 0.04), ("-+", 3.06, 3.08)])
def test_interval_certificate_close_to_the_degenerate_squares(region, lo, hi):
    result = ci.cover(region, lo, hi, 80)
    assert result["all"] and result["failed_count"] == 0


def test_planar_transverse_bound_dominates_the_floating_point_value():
    from mpmath import iv
    contour = fc.tuned_contour("--", 1.1)
    angle = iv.mpf([1.1, 1.1])
    p = ci._square("--", angle)
    for j in range(4):
        bound = ci._transverse_bound_planar(p[(j - 1) % 4], p[j])
        assert fc.transverse_norm(contour, j) <= float(bound.b) * (1 + 1e-12)


@pytest.mark.parametrize("extra", [1, 2, 3])
def test_contour_in_higher_dimensions_stays_on_the_feynman_side(extra):
    """Proposition 7: near the D-dimensional cones the imaginary part is negative."""
    import numpy as np
    contour = fc.tuned_contour("+-", 1.1)
    rng = np.random.default_rng(5)
    worst, count = -np.inf, 0
    for _ in range(1200):
        i = rng.integers(0, 4)
        d = np.array([rng.choice([-1.0, 1.0]), *(rng.normal(size=3) * 0.4)])
        if fc.mdot(d, d) <= 0:
            continue
        x0 = contour.x[i] + 10 ** rng.uniform(-1.5, 1.6) * d
        mu = fc.mdot(x0 - contour.x[i], x0 - contour.x[i]) * (1 + rng.normal() * 10 ** rng.uniform(-4, -1))
        if mu <= 0:
            continue
        e = rng.normal(size=extra)
        e *= np.sqrt(mu) / np.linalg.norm(e)
        real = contour.invariants(x0) - mu
        displaced, _ = fc.displaced_extended(contour, x0, e)
        r = np.sqrt((x0 - contour.centre) @ (x0 - contour.centre) + mu)
        tol = contour.delta0 / 2 if r < contour.radius else 0.1 * r
        for j in range(4):
            if abs(real[j]) <= tol:
                count += 1
                worst = max(worst, displaced[j].imag)
    assert count > 500 and worst < 0
