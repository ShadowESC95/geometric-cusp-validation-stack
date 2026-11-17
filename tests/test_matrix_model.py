import numpy as np

from geocusp.matrix_model import action, gradient, hermitian_noise, langevin


def test_gradient_matches_finite_difference_of_action():
    rng = np.random.default_rng(7)
    X = 0.3 * hermitian_noise(4, 3, rng)
    V = hermitian_noise(4, 3, rng)
    g, h = 0.7, 1e-5
    numeric = (action(X + h * V, g) - action(X - h * V, g)) / (2 * h)
    analytic = np.einsum("mij,mji->", gradient(X, g), V).real
    np.testing.assert_allclose(numeric, analytic, rtol=1e-7)


def test_langevin_keeps_matrices_hermitian_and_traceless():
    run = langevin(N=4, D=3, thermalization=20, steps=40, measure_every=10)
    np.testing.assert_allclose(run.X, run.X.conj().transpose(0, 2, 1), atol=1e-13)
    np.testing.assert_allclose(np.trace(run.X, axis1=1, axis2=2), 0, atol=1e-12)
    assert np.all(run.commutator_density >= 0)
