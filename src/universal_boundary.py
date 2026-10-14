"""Model-independent endpoint windows; direct exact covariance plug-in GLS.

Test construction uses time alone. The covariance pilot is an ODE trajectory,
as in the original comparison, so this implementation is not fully ODE-free.
"""
import numpy as np
from scipy.linalg import eigh
from .boundary_weak import quadrature_weights
from .weak_features import test_functions
from .weak_efficiency import nested_bumps, vdp_moments, covariance_blocks, assemble_covariance, gls_map


def universal_tests(t, widths, cutoff=1e-10, endpoints=('initial', 'terminal')):
    phi, der, original = nested_bumps(t, 480, cutoff)
    centers = [0. if side == 'initial' else t[-1] for side in endpoints]
    assert set(endpoints) <= {'initial', 'terminal'}
    additions = [test_functions(t, centers, width) for width in widths] if centers else []
    phi = np.vstack([phi]+[p for p, d in additions])
    der = np.vstack([der]+[d for p, d in additions])
    w = quadrature_weights(t)
    gram = (phi*w)@phi.T+(der*w)@der.T
    vals, vectors = eigh((gram+gram.T)/2)
    keep = vals > cutoff*vals.max()
    transform = (vectors[:, keep]/np.sqrt(vals[keep])).T
    return transform@phi, transform@der, {'original_basis': original,
        'half_windows': len(centers)*len(widths), 'retained_scalar_tests': int(keep.sum()),
        'basis_is_model_independent': True, 'relative_cutoff': cutoff}


def endpoint_blocks(t, model, phi, der, sigma, tau):
    w = quadrature_weights(t); x, v = model.T; h = t[1]-t[0]
    result = {}
    for key, kernel in [('initial', np.exp(-t/tau)), ('terminal', np.exp(-(t[-1]-t)/tau))]:
        result[key] = {'P': (phi*w)@kernel, 'D': (der*w)@kernel,
                       'XV': (phi*(w*x*v))@kernel, 'XX': (phi*(w*x*x))@kernel}
    result['trace'] = phi[:, -1]
    result['rho_full'] = float(np.exp(-(len(t)-1)*h/tau))
    return result


def augmented_covariance(q, blocks, higher, endpoints, sigma):
    C = assemble_covariance(q, blocks, higher, sigma)
    K = len(endpoints['trace']); size = 2*K; a, b, c, d = q
    result = np.zeros((size+2, size+2)); result[:size, :size] = C
    for j in range(2):
        E = np.zeros(size); E[j::2] = endpoints['trace']
        def cross(values):
            U = np.empty(size)
            if j == 0:
                U[0::2] = -values['D']; U[1::2] = -b*values['P']-2*d*values['XV']
            else:
                U[0::2] = -a*values['P']; U[1::2] = -values['D']-c*values['P']-d*values['XX']
            return sigma[j]**2*U
        terminal = cross(endpoints['terminal'])
        result[:size, :size] += np.outer(terminal, E)+np.outer(E, terminal)+sigma[j]**2*np.outer(E, E)
        initial = cross(endpoints['initial'])+sigma[j]**2*endpoints['rho_full']*E
        result[:size, size+j] = initial; result[size+j, :size] = initial
        result[size+j, size+j] = sigma[j]**2
    return (result+result.T)/2


def fit(t, y, model, qpilot, phi, der, initial, sigma, tau, iterations=100,
        tolerance=1e-7, relaxation=1e-10, mean_moments=False, controls=True,
        first_order=False, return_variance=False):
    weights = quadrature_weights(t)
    A, b = vdp_moments(t, y, phi, der, np.zeros(2) if mean_moments else sigma, weights)
    b = b+((phi[:, -1, None]*y[-1])-phi[:, 0, None]*initial).ravel()
    size = len(b)
    if controls:
        A = np.vstack((A, np.zeros((2, 4)))); b = np.r_[b, y[0]-initial]
    blocks, higher = covariance_blocks(t, model, phi, der, sigma, tau, weights)
    if first_order: higher = np.zeros_like(higher)
    ends = endpoint_blocks(t, model, phi, der, sigma, tau)
    q = np.array(qpilot).copy(); history = []
    for _ in range(iterations):
        C = augmented_covariance(q, blocks, higher, ends, sigma)
        if not controls: C = C[:size, :size]
        P, level = gls_map(A, C, relaxation); new = P@b
        history.append(float(np.linalg.norm(new-q)/max(np.linalg.norm(q), 1e-30))); q = new
        if history[-1] < tolerance: break
    diag = {'iterations': len(history), 'converged': history[-1] < tolerance,
                'last_change': history[-1], 'relaxation': level,
                'tests_built_from_model': False, 'covariance_pilot_requires_ode': True}
    if return_variance: diag['fixed_design_covariance'] = (P@C@P.T).tolist()
    return q, diag
