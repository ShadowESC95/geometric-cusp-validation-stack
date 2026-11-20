"""Tests for geocusp.scaled_limit: the conditions of the certificate and of the phase
bound in scaled variables, uniformly up to the degenerate squares (GATE_6_PHASE 7.1)."""
import math

import numpy as np
import pytest
import sympy as sp
from mpmath import iv, mp

from geocusp import certificate_iv as ci
from geocusp import phase_bound_iv as pbi
from geocusp import scaled_limit as sl


@pytest.fixture
def dps50():
    old = iv.dps
    iv.dps = 50
    yield
    iv.dps = old


def _angle(end, d):
    return float(mp.pi - mp.mpf(d)) if end == "pi" else float(d)


# ---------------------------------------------------------------------------------------
# The two closed forms used in place of interval evaluation
# ---------------------------------------------------------------------------------------
def test_closed_forms_symbolic():
    """2 x.y and the planar transverse bound 1 + |v|^2/(x.y)^2 for every pair of edges."""
    c, sn, d = sp.symbols("c sn d", real=True)
    one, zero = sp.Integer(1), sp.Integer(0)
    vec = {"a": [one, zero, zero, one], "b": [one, zero, zero, -one],
           "m": [-one, -sn, zero, -c], "q": [-one, sn, zero, c]}
    dot = lambda u, v: u[0] * v[0] - u[1] * v[1] - u[2] * v[2] - u[3] * v[3]

    def bound(a, b):                       # certificate_iv._transverse_bound_planar, squared
        v = [a[1] * b[3] - a[3] * b[1], a[0] * b[3] - a[3] * b[0], 0, a[1] * b[0] - a[0] * b[1]]
        return 1 + sum(x * x for x in v) / dot(a, b) ** 2

    on_circle = {sn**2: 1 - c**2}
    expected = {("a", "b"): (4, 2), ("m", "q"): (4, 2),
                ("q", "a"): (-2 * (1 + c), 4 / (1 + c)), ("b", "m"): (-2 * (1 + c), 4 / (1 + c)),
                ("a", "m"): (-2 * (1 - c), 4 / (1 - c)), ("b", "q"): (-2 * (1 - c), 4 / (1 - c))}
    for (x, y), (D, B2) in expected.items():
        for u, v in ((x, y), (y, x)):
            assert sp.simplify(sp.expand(2 * dot(vec[u], vec[v])).subs(on_circle) - D) == 0
            assert sp.simplify(sp.expand(sp.numer(sp.together(bound(vec[u], vec[v]) - B2))).subs(on_circle)) == 0
    # the half-angle forms used for the factors f = 1 +- cos(angle)
    s = sp.sin(d / 2)
    assert sp.simplify(1 + sp.cos(sp.pi - d) - 2 * s**2) == 0
    assert sp.simplify(1 - sp.cos(sp.pi - d) - 2 * (1 - s**2)) == 0
    assert sp.simplify(sp.sin(d) - 2 * s * sp.cos(d / 2)) == 0


@pytest.mark.parametrize("region,end", sl.CASES)
def test_geometry_matches_interval_code(region, end, dps50):
    """The scaled primitives at sigma = sin(d/2), d > 0, enclose the same numbers as the
    interval code of the repository at the angle."""
    for d in (0.3, 1e-2, 1e-5):
        angle = _angle(end, d)
        p = ci._square(region, iv.mpf(angle))
        dd = (iv.pi - iv.mpf(angle)) if end == "pi" else iv.mpf(angle)
        sig = iv.sin(dd / 2)
        g = sl.geometry(region, end, sl.Box(sig.a, sig.b))
        for j in range(4):
            Dref = 2 * ci._mdot(p[(j - 1) % 4], p[j])
            Bref = ci._transverse_bound_planar(p[(j - 1) % 4], p[j])
            for mine, ref in ((sl.value(g["Dsigned"][j]), Dref), (sl.value(g["Bj"][j]), Bref)):
                assert mine.a <= ref.b and ref.a <= mine.b
                with mp.workdps(60):
                    assert abs(mp.mpf(mine.mid) / mp.mpf(ref.mid) - 1) < 1e-25


# ---------------------------------------------------------------------------------------
# The rewriting against the unscaled repository code
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("region,end", sl.CASES)
def test_rewriting_matches_repository_exact_constants(region, end, dps50):
    """Same constants (as intervals), same eps and omega: every named intermediate of the
    scaled evaluation equals the repository's local variable of that name, at 50 digits."""
    for d in (0.3, 1e-3, 1e-6):
        r = sl.compare_with_repository(region, end, d, exact=True)
        assert r["count"] >= (60 if region == "--" else 160)
        assert r["missing"] == [] and r["unmatched"] == ["angle"]
        assert r["worst_rel"] < 1e-30, r
        # kappa and nu are one interval in the repository here, so its (B2) kappa.b <= nu.a fails;
        # the scaled version has kappa = nu as a point at the exponent 4
        assert {k for _, k, _, _ in r["bool_mismatch"]} <= {"B2", "all_conditions"}
        # arcsin bounds are point intervals at the upper end in both codes, and differ below 1e-30
        assert {k for _, k in r["overlap_fail"]} <= {"eps_V", "asin_s", "b_bl", "bound"}


@pytest.mark.parametrize("region,end", sl.CASES)
def test_repository_passes_with_family_contour(region, end):
    """A real feynman_contour.Contour with the family's constants (float64) at a square:
    the repository's own interval certificate and phase bound pass, and agree with the scaled
    evaluation to rounding."""
    for d in (0.3, 1e-2, 1e-4, 1e-6):
        r = sl.compare_with_repository(region, end, d, exact=False)
        assert all(r["repository_all"]) and all(r["scaled_all"])
        assert r["bool_mismatch"] == []
        assert r["worst_rel"] < 1e-14


# ---------------------------------------------------------------------------------------
# The scaled arithmetic itself
# ---------------------------------------------------------------------------------------
def test_scaled_enclosures_contain_point_values():
    """For every named quantity, the scaled interval on a wide box contains its value at points
    of the box (same exponent k): the bookkeeping m sigma^k is an enclosure.  The arcsin bounds
    and what is built from them are upper bounds (a point at the upper end), not enclosures:
    for those only the upper ends are compared."""
    upper_only = {"eps_V", "K_in", "K_3d", "b_in", "asin_s", "b_bl", "bound"}
    wide = sl.Box(0, 0.05)
    big_c, big_p = {}, {}
    sl.scaled_certificate("+-", "pi", wide, vals=big_c)
    sl.scaled_phase_bound("+-", "pi", wide, vals=big_p)
    for sig in (0.0, 1e-9, 1e-3, 0.02, 0.05):
        pt_c, pt_p = {}, {}
        sl.scaled_certificate("+-", "pi", sl.Box(sig, sig), vals=pt_c)
        sl.scaled_phase_bound("+-", "pi", sl.Box(sig, sig), vals=pt_p)
        for big, pt in ((big_c, pt_c), (big_p, pt_p)):
            for name, x in pt.items():
                if not isinstance(x, sl.S):
                    continue
                y = big[name]
                if y.k != x.k:      # smin at sigma = 0 alone keeps only the leading term
                    continue
                assert x.m.b <= y.m.b, name
                if name not in upper_only:
                    assert y.m.a <= x.m.a, name


def test_smin_and_comparisons():
    box = sl.Box(1e-4, 1e-2)
    a, b = sl.S(4, 2, box), sl.S(3, 0, box)                     # 4 sigma^2 < 3
    m = sl.smin([a, b])
    assert m.k == 2 and m.m.a >= 4 - 1e-25 and m.m.b <= 4 + 1e-25
    assert sl.lt(a, b) and not sl.lt(b, a)
    assert sl.lt(0, a) and sl.le(a, sl.S(4, 2, box))
    assert not sl.lt(sl.S(1, 2, box), sl.S(1, 2, box))
    s = sl.ssqrt(sl.S(9, 3, box))
    assert s.k == sl.Fraction(3, 2)


# ---------------------------------------------------------------------------------------
# The result
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("region,end", sl.CASES)
def test_all_conditions_up_to_the_limit(region, end):
    """One box: every square with 0 < d <= 1e-6 (sigma in [0, sin(5e-7)])."""
    r = sl.scaled_cover(region, end, 1e-6, boxes=1)
    assert r["all"], r["failed"]
    assert r["sigma"][0] == 0.0 and r["sigma"][1] >= math.sin(5e-7)
    assert r["excursion_upper"] < 0.054 and r["jacobian_upper"] < 0.302 and r["mu_Q_upper"] < 1e-7
    if region != "--":
        assert r["bound_upper"] < 0.1875 and r["b_in_upper"] < 0.0985 and r["Lambda_lower"] > 5.08


@pytest.mark.parametrize("region,end", sl.CASES)
def test_overlap_with_the_uniform_cover(region, end):
    """The family alone covers 0 < d <= 0.1, which overlaps the uniform interval cover of the
    repository (angles up to 3.087, i.e. d >= 0.0546; in (-,-) also angles from 0.055)."""
    r = sl.scaled_cover(region, end, 0.1, boxes=100)
    assert r["all"], r["failed"]
    assert 0.1 > math.pi - 3.087 and 0.1 > 0.055
    if region != "--":
        assert r["bound_upper"] < 0.189


@pytest.mark.parametrize("region,end", sl.CASES)
def test_limit_values(region, end):
    """At sigma = 0 itself (the scaled conditions are continuous there)."""
    box = sl.Box(0, 0)
    vals = {}
    c = sl.scaled_certificate(region, end, box, vals=vals)
    assert c["all"]
    assert abs(sl.scaled_gap(region, end, box) - (1 - math.sqrt(3) / 2)) < 1e-8
    half_root5 = iv.sqrt(5) / 2
    assert vals["Delta"].m.a <= half_root5.b and half_root5.a <= vals["Delta"].m.b
    assert vals["Dmin"].k == 2 and vals["lam0"].k == 2 and vals["kf"].k == 1 and vals["rho1"].k == 6
    if region != "--":
        pv = {}
        ph = sl.scaled_phase_bound(region, end, box, vals=pv)
        assert ph["all_conditions"] and ph["R7_outer_exact"]
        assert abs(ph["Bn2_Dp_upper"] - 8) < 1e-12
        assert pv["mu_m"].k == sl.Fraction(3, 2) and pv["F"].k == sl.Fraction(1, 2)
        assert 0.187 < ph["bound_upper"] < 0.1875


def test_regions_are_mirror_images():
    """(-,+) is a cyclic relabelling of (+,-) (edges [q,a,b,m] against [a,b,m,q]): the same
    results up to interval widths, in the scaled code, and in the repository code up to the
    float64 rounding of tuned_contour's constants (computed from the edges in another order)."""
    def same(a, b, skip, tol):
        assert a.keys() == b.keys()
        for k in a:
            if k in skip:
                continue
            if isinstance(a[k], float):
                assert abs(a[k] - b[k]) <= tol * max(abs(a[k]), 1e-300), k
            else:
                assert a[k] == b[k], k

    box = sl.box_from_d(0, 1e-6)
    same(sl.scaled_phase_bound("+-", "pi", box), sl.scaled_phase_bound("-+", "pi", box), {"region"}, 1e-10)
    same(sl.scaled_certificate("+-", "pi", box), sl.scaled_certificate("-+", "pi", box), {"region"}, 1e-10)
    lo, hi = 3.1, 3.1001
    same(pbi.interval_phase_bound("+-", lo, hi), pbi.interval_phase_bound("-+", lo, hi),
         {"region", "ratio_s_over_t"}, 1e-5)


def test_conditions_are_not_vacuous():
    """Changing one constant of the family makes the expected condition fail on the same box."""
    box = sl.box_from_d(0, 1e-6)

    def fails(**change):
        fam = dict(sl.FAMILY, **change)
        out = sl.scaled_certificate("+-", "pi", box, fam)
        return {k for k, v in out.items() if v is False} - {"all"}

    assert fails(c2=3.5 * sl.FAMILY["c2"]) == {"E_local_diffeomorphism"}
    assert fails(c2=4 * sl.FAMILY["c2"]) == {"E_local_diffeomorphism", "D_excursion_small"}
    assert {"C0_centres", "I0_c_u_positive"} <= fails(tau0=3.0)
    assert fails(e_h1=0.14) == {"H1_gap_certified"}            # above 1 - sqrt3/2 = 0.13397
    assert fails(c1=0.07) == {"H1_frame_on"}                   # 2 delta0 = 0.14 |D'| > eps
    assert "D_excursion_small" in fails(c3=16 * sl.FAMILY["c3"])
    assert "I0_c_u_positive" in fails(omega=0.6)
    with pytest.raises(ArithmeticError):                       # (B4) etc. need v^2 < delta0
        sl.scaled_phase_bound("+-", "pi", box, dict(sl.FAMILY, c3=64 * sl.FAMILY["c3"]))
