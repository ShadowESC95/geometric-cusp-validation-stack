"""Exact expansion of harmonic polylogarithms at x = -1 (geocusp.hpl_exact)."""
import os
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from geocusp import hpl
from geocusp import hpl_exact as hx

PI2 = hx.PI2
needs_data = pytest.mark.skipif(not os.environ.get("HKM_ANC"),
                                reason="set HKM_ANC to the ancillary directory of arXiv:1911.10174")


def same(a, b):
    return (a - b).is_zero()


def test_double_shuffle_relations_determine_every_multiple_zeta_value_to_weight_six():
    log = []
    table = hx.solve_mzv(6, log)
    # (weight, unknowns, equations): full rank and consistency are asserted by the solver
    assert [(n, k) for n, k, _ in log] == [(2, 1), (3, 2), (4, 4), (5, 8), (6, 16)]
    z = lambda *c: table[hx.word_of(c)]
    assert same(z(2, 1), hx.Z3)                                   # Euler
    assert same(z(4), PI2 ** 2 / 90)
    assert same(z(3, 1), PI2 ** 2 / 360)
    assert same(z(2, 2), PI2 ** 2 / 120)
    assert same(z(2, 1, 1), z(4))                                 # duality
    assert same(z(6), PI2 ** 3 / 945)
    assert same(z(3, 2), hx.ZETA2 * hx.Z3 * 3 - hx.Z5 * Fraction(11, 2))
    assert same(z(4, 1), hx.Z5 * 2 - hx.ZETA2 * hx.Z3)
    assert same(z(3, 3), (hx.Z3 * hx.Z3 - z(6)) / 2)


def test_parser_refuses_what_it_cannot_represent_exactly():
    for text in ["1.5*HPL[{0}, x]", "1e3*HPL[{0}, x]", "Pi*HPL[{0}, x]", "HPL[{1}, x]",
                 "HPL[{0}, x]*HPL[{-1}, x]", "Log[x]", "HPL[{0}, x]/Zeta[3]^2"]:
        with pytest.raises(Exception):
            hx.parse_exact(text)


def test_table_agrees_with_direct_summation():
    # a coarse guard only: the double sums converge slowly
    mp.mp.dps = 20
    table = hx.solve_mzv(5)
    for a, b in [(2, 1), (3, 1), (2, 2), (3, 2), (2, 3), (4, 1)]:
        inner = (lambda k: mp.harmonic(k - 1)) if b == 1 else (lambda k: mp.zeta(b) - mp.zeta(b, k))
        direct = mp.nsum(lambda k: inner(k) / k ** a, [2, mp.inf])
        assert abs(table[hx.word_of((a, b))].numeric() - direct) < 2e-3   # the sum converges slowly


def test_shuffle_regularisation_of_divergent_words():
    # zeta(1) = 0, and (1) sh (0,1) = (1,0,1) + 2 (0,1,1): zeta_reg(1,0,1) = -2 zeta(2,1)
    assert hx.zeta_reg((1,)).is_zero() and hx.zeta_reg((0,)).is_zero()
    assert same(hx.zeta_reg((1, 0, 1)), hx.Z3 * -2)


def test_logarithms_are_detected():
    # log(1+x) and Li_2(-x) are not holomorphic at x = -1
    for text in ["HPL[{-1}, x]", "HPL[{0, -1}, x]", "HPL[{-1, 0}, x]"]:
        assert hx.holomorphy_report(hx.parse_exact(text))["logarithmic"]
    # H(0,0;x) = log^2(x)/2 is, with value (i pi)^2/2
    report = hx.holomorphy_report(hx.parse_exact("HPL[{0, 0}, x]"))
    assert not report["logarithmic"] and same(report["value"], hx.IPI * hx.IPI / 2)


def test_two_loop_weight_is_holomorphic_at_minus_one_exactly():
    """Second proof of GATE_6_BOUNDS section 2: no logarithm, value 11 pi^4/45."""
    coefficients = hx.parse_exact(hpl.F2_HKM)
    weights = hpl.HPLWeights({"F2": hpl.F2_HKM})
    mp.mp.dps = 20
    for sign in (+1, -1):
        report = hx.holomorphy_report(coefficients, sign)
        assert report["conditions"] == 3 and not report["logarithmic"]
        assert same(report["value"] * 64, PI2 ** 2 * Fraction(11, 45))
        y = 0.3 * np.exp(0.4j * sign)
        direct = weights.value("F2", -1 + y, weights.at(-1 + y))
        assert abs(complex(hx.series_at_minus_one(report["coefficients"], y)) - direct) < 1e-11


def test_a_wrong_coefficient_breaks_holomorphy():
    text = hpl.F2_HKM.replace("(Pi^2*HPL[{-1, 0}, x])/32", "(Pi^2*HPL[{-1, 0}, x])/31")
    assert hx.holomorphy_report(hx.parse_exact(text))["logarithmic"]


@needs_data
def test_three_loop_weight_is_holomorphic_at_minus_one_exactly():
    read = lambda name: open(os.path.join(os.environ["HKM_ANC"], name)).read()
    texts = {"P": read("HKM_F3_P_HPL.txt"), "NP1": read("HKM_F3_NP1_HPL.txt"), "NP2": read("HKM_F3_NP2_HPL.txt")}
    weights = hpl.HPLWeights(texts)
    mp.mp.dps = 20
    reports = {}
    for name, text in texts.items():
        coefficients = hx.parse_exact(text)
        for sign in (+1, -1):
            report = reports[name, sign] = hx.holomorphy_report(coefficients, sign)
            assert not report["logarithmic"], name
            y = 0.3 * np.exp(0.4j * sign)
            direct = weights.value(name, -1 + y, weights.at(-1 + y))
            assert abs(complex(hx.series_at_minus_one(report["coefficients"], y)) - direct) < 1e-10
    assert [reports[n, 1]["conditions"] for n in ("P", "NP1", "NP2")] == [15, 4, 7]
    z3sq, pi6 = hx.Z3 * hx.Z3, PI2 ** 3
    assert same(reports["P", 1]["value"], z3sq / 32 - pi6 * Fraction(67, 51840))
    assert reports["NP1", 1]["value"].is_zero() and reports["NP1", -1]["value"].is_zero()
    nonplanar = reports["NP1", 1]["slope"] + reports["NP2", 1]["value"]
    assert same(nonplanar, z3sq * Fraction(21, 8) + PI2 * hx.Z3 * Fraction(3, 8) - hx.Z5 * Fraction(15, 16)
                + pi6 * Fraction(49, 8640))
    assert same(nonplanar, reports["NP1", -1]["slope"] + reports["NP2", -1]["value"])
