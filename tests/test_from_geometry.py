"""The two-loop form from its inequalities, and an independent membership test."""
import random

import numpy as np
import pytest
import sympy as sp
from scipy.optimize import linprog

from geocusp import amplituhedron as am
from geocusp import negative_geometry as ng


def test_two_loop_form_from_inequalities_equals_the_integrand():
    """The sum of the d log forms of the five pieces is -n12 / (D1 D2 X12)."""
    _, pieces, total = ng.two_loop_form_from_inequalities()
    _, integrand, _, _ = ng.symbolic_two_loop_integrand()
    assert len(pieces) == 5
    assert sp.simplify(total - integrand) == 0


def test_the_pieces_tile_the_region():
    """Every point with Q > 0 lies in exactly one piece; points with Q < 0 in none."""
    variables, _, _ = ng.two_loop_form_from_inequalities()
    random.seed(1)
    inside = 0
    for _ in range(5000):
        pt = {str(s): random.uniform(0.05, 3) for s in variables}
        Q = (pt["x1"] - pt["x2"]) * (pt["w1"] - pt["w2"]) + (pt["y1"] - pt["y2"]) * (pt["z1"] - pt["z2"])
        count = sum(ng.region_piece(i, pt) for i in range(5))
        assert count == (1 if Q > 0 else 0)
        inside += Q > 0
    assert 2000 < inside < 3000


def test_each_piece_is_positive_on_itself():
    variables, pieces, _ = ng.two_loop_form_from_inequalities()
    functions = [sp.lambdify(variables, density) for _, density in pieces]
    random.seed(2)
    seen = [0] * 5
    for _ in range(4000):
        values = [random.uniform(0.05, 3) for _ in variables]
        pt = dict(zip(map(str, variables), values))
        for i in range(5):
            if ng.region_piece(i, pt):
                seen[i] += 1
                assert functions[i](*values) > 0
    assert min(seen) > 50


@pytest.mark.parametrize("k,n", [(1, 6), (1, 8), (2, 6), (2, 8)])
def test_sampler_passes_the_sign_flip_test_and_random_points_fail(k, n):
    result = am.standard_amplituhedron_pushforward(k, n, 200, seed=3)
    assert all(am.sign_flip_membership(Y, result.Z, k)["inside"] for Y in result.Y)
    rng = np.random.default_rng(0)
    rejected = sum(not am.sign_flip_membership(rng.normal(size=(k, n)) @ result.Z, result.Z, k)["inside"]
                   for _ in range(200))
    assert rejected > 150


def test_sign_flip_test_agrees_with_linear_programming_for_k1():
    """For k = 1 the amplituhedron is a cyclic polytope: Y is inside iff
    Y = c Z with c >= 0. Compare the two tests on mixed points."""
    n = 7
    result = am.standard_amplituhedron_pushforward(1, n, 10, seed=5)
    Z = result.Z
    rng = np.random.default_rng(4)
    agree = inside = 0
    trials = 400
    for _ in range(trials):
        c = rng.normal(size=n) + rng.choice([0.0, 0.8])       # some mostly positive
        Y = (c @ Z)[None, :]
        flips = am.sign_flip_membership(Y, Z, 1)["inside"]
        lp = linprog(np.zeros(n), A_eq=Z.T, b_eq=Y[0], bounds=[(0, None)] * n, method="highs")
        agree += flips == (lp.status == 0)
        inside += flips
    assert agree == trials
    assert 20 < inside < trials - 20


def test_residue_on_a_facet_is_the_form_of_the_facet():
    """One level of the residue recursion at w1 = 0."""
    residue, facet, difference = ng.two_loop_boundary_check()
    assert difference == 0
    assert residue != 0


def test_residue_on_the_mutual_facet_is_the_form_of_that_facet():
    """Q = 0, with the orientation reversed where A < 0."""
    residue, facet, ratio = ng.two_loop_mutual_facet_check()
    assert ratio == 1


def test_second_level_of_the_residue_recursion():
    differences = ng.two_loop_second_level_check()
    assert set(differences) == {"w1=w2=0", "w1=x1=0", "w1=x2=0"}
    assert all(d == 0 for d in differences.values())
