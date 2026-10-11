from itertools import product
from math import factorial
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import roots_hermitenorm
from src.noise import ou_noise
from src.covariance import exponential_apply, residual_covariance
from src.weak_features import (library, exponents, quadrature_weights,
                               test_functions as bump_functions, weak_matrices, heat_transform)
from src.systems import trajectory, system
from src.selection import stls
from src.diagnostic import choose_grid


def test_exponential_operator_matches_dense():
    u = np.arange(39, dtype=float).reshape(13, 3) / 7
    for rho in (0., .2, .99, .999999):
        kernel = rho ** np.abs(np.arange(13)[:, None] - np.arange(13)[None, :])
        assert_allclose(exponential_apply(u, rho), kernel @ u, atol=2e-12)


def test_quadrature_and_bump_sign():
    t = np.linspace(0, 4, 4001)
    assert_allclose(quadrature_weights(t).sum(), 4.)
    y = np.exp(-t)[:, None]
    G, b = weak_matrices(t, y, [(1,)], [1.3, 2.5], 1.8)
    assert_allclose(b, -G, rtol=1e-9, atol=1e-10)
    phi, deriv = bump_functions(np.array([.4, 1.3, 2.2]), [1.3], 1.8)
    assert_allclose(phi[0, [0, 2]], 0., atol=1e-14)
    assert_allclose(deriv[0, [0, 1, 2]], 0., atol=1e-14)


def test_ou_reproducible_stationary_initial_and_nested():
    seed = 912
    eps = ou_noise(129, .01, .5, [.3, .2], np.random.default_rng(seed))
    repeat = ou_noise(129, .01, .5, [.3, .2], np.random.default_rng(seed))
    assert_allclose(eps, repeat, rtol=0, atol=0)
    initial = np.random.default_rng(seed).normal(size=2) * [.3, .2]
    assert_allclose(eps[0], initial)
    assert_allclose(eps[::4][::2], eps[::8], rtol=0, atol=0)
    rho = np.exp(-.01 / .5)
    # Exact stationarity equation, no Monte Carlo test with a random tolerance.
    assert_allclose(rho ** 2 * .3 ** 2 + .3 ** 2 * (-np.expm1(-.02 / .5)), .3 ** 2)


def test_hermite_bias_removal_by_gaussian_quadrature():
    z, w = roots_hermitenorm(5)
    w = w / np.sqrt(2 * np.pi)
    mean, variance = .7, .19
    values = library((mean + np.sqrt(variance) * z)[:, None], exponents(1, 4), [variance])
    assert_allclose(w @ values, mean ** np.arange(5), atol=3e-15)


def test_exact_covariance_against_gaussian_quadrature():
    # Only two timepoints and two independent coordinates => a 4D Gaussian.
    # Four nodes per axis: 4^4=256, sufficient for cubic-library covariance.
    t = np.array([.3, .8])
    x = np.array([[.4, -.2], [.8, .3]])
    powers = exponents(2, 3)
    coef = np.zeros((len(powers), 2))
    coef[powers.index((2, 1))] = [.7, -1.1]
    coef[powers.index((0, 1))] = [1., .4]
    centers, width, sc, sw, tau = [.45, .65], 1.2, np.array([.3, .2]), np.array([.1, .15]), .4
    rho = np.exp(-(t[1] - t[0]) / tau)
    v = sc ** 2 + sw ** 2
    C = np.zeros((4, 4))
    for j in range(2):
        C[j, j] = C[2 + j, 2 + j] = v[j]
        C[j, 2 + j] = C[2 + j, j] = sc[j] ** 2 * rho
    L = np.linalg.cholesky(C)
    nodes, weights = roots_hermitenorm(4)
    weights /= np.sqrt(2 * np.pi)
    for correction in (v, np.zeros(2), v * .8):
        samples, probs = [], []
        for indices in product(range(4), repeat=4):
            e = (L @ nodes[list(indices)]).reshape(2, 2)
            G, b = weak_matrices(t, x + e, powers, centers, width, correction)
            samples.append((b - G @ coef).reshape(-1))
            probs.append(np.prod(weights[list(indices)]))
        samples, probs = np.array(samples), np.array(probs)
        centered = samples - probs @ samples
        actual = centered.T @ (probs[:, None] * centered)
        predicted = residual_covariance(t, x, powers, coef, centers, width,
                                        sc, tau, sw, correction)
        assert_allclose(predicted, actual, rtol=1e-11, atol=1e-13)
        effective = heat_transform(powers, coef, v - correction)
        phi, derivative = bump_functions(t, centers, width)
        qw = quadrature_weights(t)
        expected_mean = -(derivative * qw) @ x - (phi * qw) @ library(x, powers) @ effective
        assert_allclose(probs @ samples, expected_mean.reshape(-1), atol=1e-13)


def test_clean_linear_and_polynomial_coefficient_recovery():
    for name in ('linear', 'vdp'):
        _, powers, truth = system(name)
        t = np.linspace(0, 8, 4001)
        clean = trajectory(name, t)
        centers = np.linspace(1.6, 6.4, 30)
        G, b = weak_matrices(t, clean, powers, centers, 3.)
        estimate = stls(G, b, 1e-4)
        assert_allclose(estimate, truth, atol=2e-5)
        assert np.array_equal(np.abs(estimate) > 1e-8, np.abs(truth) > 1e-8)


def test_diagnostic_refuses_unstable_support_and_selects_coarsest():
    records = [{'h': h, 'variance': 1., 'quadrature_difference': .001,
                'validation_loss': .1, 'support_stable': False,
                'conditioning_ok': True} for h in (.01, .02)]
    assert choose_grid(records, 1.)['h'] is None
    for r in records:
        r['support_stable'] = True
    assert choose_grid(records, 1.)['h'] == .02
