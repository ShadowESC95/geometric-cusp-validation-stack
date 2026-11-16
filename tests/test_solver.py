import json
from math import pi

import numpy as np
import pytest

from geocusp import solver
from geocusp.solver import Weight


def test_double_pole_tends_to_sinc():
    for p in (0., .3, -.6, .85):
        fine, finer = solver.double_pole(p, .004), solver.double_pole(p, .002)
        np.testing.assert_allclose(2*finer-fine, 4*solver.sinc(p), rtol=0, atol=2e-4)
    # mpmath reference values from examples/mellin_functional.py
    np.testing.assert_allclose(solver.double_pole(0., .01), 4.02185407684, rtol=0, atol=1e-10)
    np.testing.assert_allclose(solver.double_pole(.85, .003), 0.686852694195, rtol=0, atol=1e-10)


def test_box_single_pole_in_ratio():
    slope = (solver.double_pole(0., .004, 3.)-solver.double_pole(0., .004))/(-.004)
    np.testing.assert_allclose(slope, -2*np.log(3), rtol=0, atol=.02)


@pytest.mark.parametrize("weight", [
    Weight((-1.,)),
    Weight((pi**2, 0., 1.)),
    Weight((.5, 1.2, -.7, .3, .1), 2., -.6),
    Weight((0.,), 1., .85),
])
def test_three_functionals_agree(weight):
    exact = solver.closed_mean(weight)
    np.testing.assert_allclose(solver.circle_mean(weight), exact, rtol=0, atol=1e-10)
    for radius in (.7, 1.3):
        np.testing.assert_allclose(solver.cycle_period(weight, radius), exact, rtol=0, atol=1e-6)
    extrapolated = 2*solver.lorentzian(weight, .002)-solver.lorentzian(weight, .004)
    np.testing.assert_allclose(extrapolated, exact, rtol=0, atol=2e-3*max(1, abs(exact)))


def test_cycle_period_detects_a_wrong_cycle(monkeypatch):
    real = solver.lifted_cycle

    def broken(alpha, phi, theta, radius):
        a, u, v = real(alpha, phi, theta, radius)
        return a, u, v*u/np.conj(u)

    monkeypatch.setattr(solver, "lifted_cycle", broken)
    assert abs(solver.cycle_period(Weight((1.,)))-1) > .1


def test_known_cusp_coefficients():
    first = solver.solve_functional(Weight((-1.,)), loop_order=1)
    second = solver.solve_functional(Weight((pi**2, 0., 1.)), loop_order=2)
    assert first["agree"] and second["agree"]
    assert first["matches_reference"] and second["matches_reference"]
    np.testing.assert_allclose([first["gamma"], second["gamma"]], [4, -4*pi**2/3])


def test_cd_weight_matches_closed_form():
    for z in (.35+.8j, -.65+.4j, 2+0j, 1+1e-5j, -1.5-.35j, 1j):
        np.testing.assert_allclose(solver.cd_weight(z), pi**2+np.log(z)**2, rtol=0, atol=1e-10)
    with pytest.raises(ValueError):
        solver.cd_weight(0)


def test_rim_collisions_and_grouped_period():
    upper, lower = solver.collision_points()
    for point, expected in ((upper, -6.052350324174+6.052350324174j),
                            (lower, 3.045968993180-3.045968993180j)):
        result = solver.rim_slice(point["t"], point["lam"], point["s"])
        assert result["regime"] == "collision" and result["single_periods"] is None
        np.testing.assert_allclose(complex(*result["roots"][0]), point["K"], rtol=0, atol=1e-5)
        np.testing.assert_allclose(complex(*result["paired_period"]), expected, rtol=0, atol=1e-8)

    band = solver.rim_slice(1., 470., .27)
    assert band["regime"] == "one-meridian" and band["inside"] == 1
    singles = [complex(*value) for value in band["single_periods"]]
    np.testing.assert_allclose(sum(singles), complex(*band["paired_period"]), rtol=1e-8, atol=1e-14)


def test_bad_input_is_rejected():
    with pytest.raises(ValueError):
        Weight((1.,), 1., 1.2)
    with pytest.raises(ValueError):
        solver.double_pole(.3, 0.)
    with pytest.raises(ValueError):
        solver.rim_slice(1., -5., .3)
    with pytest.raises(SystemExit):
        solver.main(["functional", "--exponent", "2"])


def test_command_line_prints_json(capsys):
    solver.main(["functional", "--logs", "-1", "--loop-order", "1"])
    first = json.loads(capsys.readouterr().out)
    assert first["agree"] and first["gamma"] == pytest.approx(4)
    solver.main(["cd", "0.35+0.8j"])
    assert json.loads(capsys.readouterr().out)["error"] < 1e-10
    solver.main(["rim", "--t", "1", "--lam", "470", "--s", "0.27"])
    assert json.loads(capsys.readouterr().out)["regime"] == "one-meridian"


def test_desktop_window_computes_and_redraws():
    tk = pytest.importorskip("tkinter")
    from geocusp.solver_gui import SolverApp
    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display")
    root.withdraw()
    try:
        app = SolverApp(root)
        root.update()
        assert app.results["functional"]["matches_reference"]
        app.preset(-1, 0, 0, 1)
        assert app.results["functional"]["gamma"] == pytest.approx(4)
        app.preset(0, 0, 1, 2)
        app.f["p"].set(.5)
        root.update()
        app.update_functional()
        assert app.results["functional"]["closed_mean"] == pytest.approx(2/pi)
        app.to_circle()
        app.update_cd()
        assert app.results["cd"]["modulus"] == pytest.approx(1) and app.results["cd"]["error"] < 1e-10
        app.across_cut()
        app.update_cd()
        assert app.results["cd"]["value"][1]*app.z["im"].get() < 0 or app.results["cd"]["error"] < 1e-8
        app.jump(solver.collision_points()[0])
        assert app.results["rim"]["regime"] == "collision"
        app.jump({"t": 1., "lam": 470., "s": .27})
        assert app.results["rim"]["regime"] == "one-meridian"
        app.r["lam"].set(1300.)
        app.update_rim()
        assert app.results["rim"]["regime"] == "outside"
        root.update()
    finally:
        root.destroy()
