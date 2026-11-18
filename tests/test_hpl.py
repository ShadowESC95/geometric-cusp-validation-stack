import os

import mpmath as mp
import numpy as np
import pytest

from geocusp import hpl


def test_index_notation_and_parser():
    assert hpl.expand_index((-3, 0, -1)) == (0, 0, -1, 0, -1)
    source, words = hpl.parse("Pi^2*HPL[{-2, 0}, x]/3 + Zeta[3]")
    assert words == {(0, -1, 0)} and "PI**2" in source and "Z3" in source
    with pytest.raises(ValueError):
        hpl.parse("HPL[{1}, x]")
    with pytest.raises(ValueError):
        hpl.parse("Log[x] + HPL[{0}, x]")


def test_series_against_classical_polylogarithms():
    words = hpl.suffix_closure([(0, 0, -1), (-1, -1), (0, 0)])
    v = hpl.series_values(words, 0.3)
    assert v[(-1,)] == pytest.approx(np.log(1.3), abs=1e-14)
    assert v[(0, -1)] == pytest.approx(-complex(mp.polylog(2, -0.3)), abs=1e-14)
    assert v[(0, 0, -1)] == pytest.approx(-complex(mp.polylog(3, -0.3)), abs=1e-14)
    assert v[(-1, -1)] == pytest.approx(np.log(1.3) ** 2 / 2, abs=1e-14)
    assert v[(0, 0)] == pytest.approx(np.log(0.3) ** 2 / 2, abs=1e-14)


def test_continuation_to_the_unit_circle():
    w = hpl.HPLWeights({"a": "HPL[{0, 0, 0, -1}, x]"})
    y = w.at(1.0)
    assert w.value("a", 1.0, y) == pytest.approx(-complex(mp.polylog(4, -1)), abs=1e-11)
    x = np.exp(2.5j)
    assert w.value("a", x, w.at(x)) == pytest.approx(-complex(mp.polylog(4, -x)), abs=1e-10)


def test_circle_means_give_the_two_and_three_loop_cusp():
    w = hpl.HPLWeights({"F1": hpl.F1_HKM, "F2": hpl.F2_HKM})
    means = w.circle_mean(nodes=24)
    assert -2 * means["F1"] == pytest.approx(-float(mp.zeta(2)) / 2, abs=1e-11)
    assert -4 / 3 * means["F2"] == pytest.approx(11 * float(mp.zeta(4)) / 8, abs=1e-10)
    # normalisation against the weights used elsewhere in the repository
    assert 64 * means["F2"] == pytest.approx(-11 * np.pi ** 4 / 15, abs=1e-8)


def test_two_loop_weight_is_single_valued_at_minus_one_on_the_principal_sheet():
    w = hpl.HPLWeights({"F2": hpl.F2_HKM})
    x = -1 + 0.3 * np.exp(0.4j)
    y = w.at(x)
    assert abs(w.value("F2", x, w.around(y, x, -1)) - w.value("F2", x, y)) < 1e-10


@pytest.mark.skipif(not os.environ.get("HKM_ANC"), reason="set HKM_ANC to the ancillary directory of arXiv:1911.10174")
def test_three_loop_weight_gives_the_four_loop_cusp():
    report = hpl.three_loop_report(os.environ["HKM_ANC"], nodes=40)
    for _name, got, expected in report["cusp"]:
        assert got == pytest.approx(expected, abs=2e-9)
    assert max(report["monodromy_at_minus_one"].values()) < 1e-9
    assert abs(report["values_near_minus_one"]["NP1"]) < 1e-5


def test_two_loop_weight_has_no_branch_point_at_minus_one_to_thirty_digits():
    result = hpl.monodromy_at_minus_one({"F2": hpl.F2_HKM}, dps=35)
    assert result["F2"] < 1e-30
    # a function that does have one
    assert hpl.monodromy_at_minus_one({"log": "HPL[{-1}, x]"}, dps=35)["log"] == pytest.approx(2 * np.pi, abs=1e-12)


@pytest.mark.skipif(not os.environ.get("HKM_ANC"), reason="set HKM_ANC to the ancillary directory of arXiv:1911.10174")
def test_three_loop_weight_has_no_branch_point_at_minus_one_to_forty_digits():
    read = lambda name: open(os.path.join(os.environ["HKM_ANC"], name)).read()
    result = hpl.monodromy_at_minus_one(
        {"P": read("HKM_F3_P_HPL.txt"), "NP1": read("HKM_F3_NP1_HPL.txt"), "NP2": read("HKM_F3_NP2_HPL.txt")},
        combine=lambda v, x: {"NP": v["NP1"] / (1 + x) + v["NP2"]})
    assert result["P"] < 1e-40 and result["NP"] < 1e-40
