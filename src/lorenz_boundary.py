"""Endpoint-enabled Lorenz residuals with exact degree-two covariance."""
import numpy as np
from .boundary_weak import quadrature_weights
from .weak_efficiency import gls_map
from scripts.run_lorenz_efficiency import moments, blocks, assemble, jacobians


def fit(t, y, model, qpilot, phi, der, initial, sigma, tau, controls=True,
        iterations=100, tolerance=1e-7, relaxation=1e-10, first_order=False, return_variance=False):
    w = quadrature_weights(t); A, b = moments(t, y, phi, der, w)
    b += (phi[:, -1, None]*y[-1]-phi[:, 0, None]*initial).ravel()
    size = len(b)
    if controls: A = np.vstack((A, np.zeros((3, 3)))); b = np.r_[b, y[0]-initial]
    B = blocks(t, model, phi, der, sigma, tau, first_order=first_order, weights=w); J = jacobians(model)
    q = np.array(qpilot).copy(); change = np.inf
    for iteration in range(iterations):
        C = assemble(q, B); full = np.zeros((size+3, size+3)); full[:size, :size] = C
        actual_J = J[0]+np.einsum('q,qnji->nji', q, J[1:])
        for j in range(3):
            E = np.zeros(size); E[j::3] = phi[:, -1]
            U = -phi.T[:, :, None]*actual_J[:, None, :, j]; U[:, :, j] -= der.T
            U = (U*w[:, None, None]).reshape(len(t), size)
            terminal = sigma[j]**2*(U.T@np.exp(-(t[-1]-t)/tau))
            initial_cross = sigma[j]**2*(U.T@np.exp(-t/tau))+sigma[j]**2*np.exp(-t[-1]/tau)*E
            full[:size,:size] += np.outer(terminal,E)+np.outer(E,terminal)+sigma[j]**2*np.outer(E,E)
            full[:size,size+j] = initial_cross; full[size+j,:size] = initial_cross
            full[size+j,size+j] = sigma[j]**2
        C = (full+full.T)/2 if controls else (full[:size,:size]+full[:size,:size].T)/2
        P, ridge = gls_map(A, C, relaxation); new = P@b
        change = float(np.linalg.norm(new-q)/max(np.linalg.norm(q),1e-30)); q = new
        if change < tolerance: break
    diag = {'iterations': iteration+1, 'converged': change < tolerance, 'last_change': change, 'ridge': ridge}
    if return_variance: diag['fixed_design_covariance'] = (P@C@P.T).tolist()
    return q, diag
