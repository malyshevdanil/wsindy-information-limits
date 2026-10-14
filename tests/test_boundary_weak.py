"""Algebra-only checks; no ODE research runs on the local machine."""
import numpy as np
from src.boundary_weak import covariance, moments


def test_known_initial_boundary_identity_for_constant_path():
    t = np.linspace(0, 1, 11); y = np.tile([2., -1.], (len(t), 1))
    slope = np.arange(8).reshape(2, 4)/10
    phi = np.ones((len(t), 2, 4))+t[:, None, None]*slope
    der = np.broadcast_to(slope, phi.shape)
    _, b = moments(t, y, phi, der, y[0], np.array([.2, .3]))
    np.testing.assert_allclose(b, 0, atol=2e-15)


def test_polynomial_covariance_and_endpoint_cross_terms_against_gaussian_draws():
    rng = np.random.default_rng(412)
    t = np.array([0., .5, 1.]); model = np.array([[2., 0.], [1.8, -.5], [1.3, -.9]])
    phi = rng.normal(size=(3, 2, 4)); der = rng.normal(size=phi.shape)
    q = np.array([1., -1., 1., -1.]); sigma = np.array([.3, .4]); tau = .7
    C = covariance(t, model, phi, der, q, sigma, tau)
    K = np.exp(-np.abs(t[:, None]-t[None, :])/tau)
    noise = np.einsum('ij,rjc->ric', np.linalg.cholesky(K), rng.normal(size=(120000, 3, 2)))*sigma
    y = model[None]+noise; w = np.array([1/6, 2/3, 1/6])
    z, v = y[:, :, 0], y[:, :, 1]
    f = np.stack((q[0]*v, q[1]*z+q[2]*v+q[3]*(z*z-sigma[0]**2)*v), axis=2)
    m = -np.einsum('nik,rni,n->rk', der, y, w)-np.einsum('nik,rni,n->rk', phi, f, w)
    m += y[:, -1]@phi[-1]-model[0]@phi[0]
    m = np.column_stack((m, noise[:, 0]))
    empirical = np.cov(m, rowvar=False)
    scale = np.sqrt(np.diag(C)[:, None]*np.diag(C)[None, :])
    np.testing.assert_allclose(empirical/scale, C/scale, atol=.012, rtol=0)
    Amean, bmean = moments(t, model, phi, der, model[0], np.zeros(2))
    mean = bmean-Amean@q
    # The expectation uses uncorrected clean features. Correcting the noiseless
    # model a second time would introduce an artificial sigma-squared defect.
    np.testing.assert_allclose((m.mean(axis=0)-mean)/np.sqrt(np.diag(C)/len(m)), 0,
                               atol=3.5, rtol=0)
