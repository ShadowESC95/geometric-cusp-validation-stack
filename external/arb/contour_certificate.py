"""python-flint/Arb certification hook for contour moments."""
from __future__ import annotations


def unit_circle_z(k, n=4096):
    try:
        from flint import acb, pi
    except Exception as exc:
        raise RuntimeError('python-flint unavailable') from exc
    twopi=2*pi
    return [acb.exp(acb(0, twopi*k/n)) for k in range(n)]


def contour_mean_constant(value, n=4096):
    try:
        from flint import arb
    except Exception as exc:
        raise RuntimeError('python-flint unavailable') from exc
    return arb(value)
