"""Checks for docs/CUSP_LOCAL_MODEL.md."""
import mpmath as mp
import pytest

from geocusp import cusp_local as cl


@pytest.fixture(autouse=True)
def precision():
    mp.mp.dps = 20


def test_measure_identity_is_exact():
    ratio, residual = cl.measure_identity()
    assert ratio in (1, -1) and residual == 0


@pytest.mark.parametrize("p", [0.3, -0.5, mp.mpc(0.2, 1.5)])
def test_triangle_double_pole_is_the_circle_mean(p):
    limit = mp.sin(mp.pi * p) / (mp.pi * p)
    errors = [abs(cl.triangle_closed_form(p, eta) - limit) for eta in (1e-2, 1e-3, 1e-4)]
    assert errors[2] < 1e-3 * abs(limit) * 5
    assert errors[1] / errors[2] == pytest.approx(10, rel=0.1)      # linear in eta


def test_triangle_is_symmetric_in_the_limit_under_p_to_minus_p():
    assert abs(cl.triangle_closed_form(0.4, 1e-6) - cl.triangle_closed_form(-0.4, 1e-6)) < 1e-4


def test_dirichlet_integral_of_the_triangle():
    numeric, closed = cl.dirichlet_check(0.3, 0.25)
    assert numeric == pytest.approx(closed, rel=1e-10)


@pytest.mark.parametrize("p", [0.3, -0.3, mp.mpc(0.2, 1.0)])
def test_vertex_integral_identity(p):
    numeric, closed, main = cl.vertex_main_term(p, 0.25, 0.1)
    assert abs(numeric - closed) < 1e-9 * abs(closed)
    assert abs(main) > 0


@pytest.mark.parametrize("p,ratio", [(0.3, 1), (-0.5, 2.5), (mp.mpc(0.2, 1.5), 1)])
def test_remainder_is_a_single_pole(p, ratio):
    """(N_eta - 4 sin(pi p)/(pi p)) / eta stays bounded as eta -> 0."""
    values = [abs(cl.remainder_ratio(p, eta, ratio)) for eta in (0.2, 0.1, 0.05, 0.025)]
    assert max(values) < 60
    assert abs(values[3] - values[2]) < abs(values[1] - values[0])     # settling


def test_four_cusps_add_to_the_box():
    """2 x triangle(p) + 2 x triangle(-p) and N_eta have the same limit."""
    p, eta = mp.mpf("0.3"), mp.mpf("1e-4")
    cusps = 2 * cl.triangle_closed_form(p, eta) + 2 * cl.triangle_closed_form(-p, eta)
    assert cusps == pytest.approx(float(4 * mp.sin(mp.pi * p) / (mp.pi * p)), rel=2e-3)
