import numpy as np
from src.polynomial_gls import covariance
from src.weak_efficiency import covariance_blocks, assemble_covariance
from src.boundary_weak import quadrature_weights
from src.support_information import POWERS


def test_full_polynomial_covariance_matches_exact_vdp_specialization():
    rng = np.random.default_rng(714); t = np.linspace(0, 1, 9)
    x = rng.normal(size=(9, 2)); phi = rng.normal(size=(3, 9)); der = rng.normal(size=phi.shape)
    sigma = np.array([.2, .3]); q = np.array([.9, -1., .8, -.7]); coeff = np.zeros((9, 2))
    coeff[1,0] = q[0]; coeff[0,1] = q[1]; coeff[1,1] = q[2]; coeff[6,1] = q[3]
    B, H = covariance_blocks(t, x, phi, der, sigma, .3, quadrature_weights(t))
    reference = assemble_covariance(q, B, H, sigma)
    actual = covariance(t, x, coeff, phi, der, POWERS, sigma, .3)
    np.testing.assert_allclose(actual, reference, atol=2e-15, rtol=2e-13)
