"""Exact corrected polynomial residual covariance on arbitrary supplied tests."""
from math import factorial
import numpy as np
from .weak_features import exponents, derivative_polynomial, library
from .boundary_weak import quadrature_weights
from .covariance import exponential_apply
from .weak_efficiency import gls_map


def moments(t, y, phi, der, powers, sigma):
    w = quadrature_weights(t); G = (phi*w)@library(y, [tuple(a) for a in powers], np.asarray(sigma)**2)
    d = y.shape[1]; n = len(powers); A = np.zeros((len(phi)*d, n*d))
    for j in range(d): A[j::d, j*n:(j+1)*n] = G
    return A, (-(der*w)@y).ravel()


def covariance(t, model, coefficients, phi, der, powers, sigma, tau, first_order=False):
    w = quadrature_weights(t); d = model.shape[1]; degree = max(map(sum, powers)); size = len(phi)*d
    C = np.zeros((size, size))
    for alpha in exponents(d, 1 if first_order else degree):
        q = sum(alpha)
        if not q: continue
        jac = derivative_polynomial(model, powers, coefficients, alpha)
        values = -phi.T[:, :, None]*jac[:, None, :]
        if q == 1: values[:, :, alpha.index(1)] -= der.T
        U = (values*w[:, None, None]).reshape(len(t), size)
        factor = float(np.prod([factorial(a) for a in alpha])*np.prod(np.asarray(sigma)**(2*np.array(alpha))))
        active = np.flatnonzero(np.any(U != 0., axis=0))
        if len(active):
            V = U[:, active]
            C[np.ix_(active, active)] += factor*(V.T@exponential_apply(V, np.exp(-q*(t[1]-t[0])/tau)))
    return (C+C.T)/2


def solve(A, b, C, relaxation=1e-10):
    scales = np.linalg.norm(A, axis=0)
    if np.any(scales == 0): raise ValueError('Zero polynomial design column')
    P, ridge = gls_map(A/scales, C, relaxation); P = P/scales[:, None]
    estimate = P@b; V = P@C@P.T
    if not np.isfinite(estimate).all() or np.any(np.diag(V) <= 0): raise ValueError('Invalid GLS estimate/variance')
    return estimate, (V+V.T)/2, {'covariance_ridge': ridge,
                               'scaled_design_condition': float(np.linalg.cond(A/scales))}
