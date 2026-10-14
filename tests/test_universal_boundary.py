import numpy as np
from src.boundary_weak import covariance, quadrature_weights
from src.weak_efficiency import covariance_blocks
from src.universal_boundary import augmented_covariance, endpoint_blocks


def test_endpoint_block_covariance_matches_vector_moment_covariance():
    rng = np.random.default_rng(228); t = np.linspace(0, 1, 9)
    model = rng.normal(size=(9, 2)); phi = rng.normal(size=(3, 9)); der = rng.normal(size=phi.shape)
    sigma = np.array([.3, .4]); q = np.array([1., -1., .8, -.9]); tau = .3
    blocks, higher = covariance_blocks(t, model, phi, der, sigma, tau, quadrature_weights(t))
    result = augmented_covariance(q, blocks, higher, endpoint_blocks(t, model, phi, der, sigma, tau), sigma)
    vector = np.zeros((9, 2, 6)); derivative = np.zeros_like(vector)
    vector[:, 0, 0::2] = phi.T; vector[:, 1, 1::2] = phi.T
    derivative[:, 0, 0::2] = der.T; derivative[:, 1, 1::2] = der.T
    reference = covariance(t, model, vector, derivative, q, sigma, tau)
    np.testing.assert_allclose(result, reference, atol=2e-15, rtol=2e-13)
