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
