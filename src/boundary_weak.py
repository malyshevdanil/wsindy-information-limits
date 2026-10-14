"""Known-initial-state weak moments with endpoint controls, without caches.

Four vector adjoint tests and two initial-observation controls. All covariances
are exact for the specified discrete polynomial moments at a fixed model path.
Pilot plug-in and quadrature errors remain separate from that identity.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import cho_factor, cho_solve
from .covariance import exponential_apply
from .profile_likelihood import ensemble_sensitivities


def quadrature_weights(t):
    h = t[1]-t[0]
    if len(t) % 2 != 1 or not np.allclose(np.diff(t), h):
        raise ValueError('Composite Simpson requires an even number of uniform intervals')
    w = np.ones(len(t)); w[1:-1:2] = 4; w[2:-1:2] = 2
    return w*h/3


def optimal_tests(t, q, initial, sigma, tau):
    model, U = ensemble_sensitivities(t, q, np.asarray(initial)[None], 1e-10, 1e-12)
    x = model[0]; S = U[0, :, :, :4]; a, b, c, d = q
    z, v = x.T
    f = np.column_stack((a*v, b*z+c*v+d*z*z*v))
    J = np.zeros((len(t), 2, 2)); J[:, 0, 1] = a
    J[:, 1, 0] = b+2*d*z*v; J[:, 1, 1] = c+d*z*z
    Jp = np.zeros_like(J)
    Jp[:, 1, 0] = 2*d*(f[:, 0]*v+z*f[:, 1])
    Jp[:, 1, 1] = 2*d*z*f[:, 0]
    B = np.zeros_like(S); B[:, 0, 0] = v
    B[:, 1, 1] = z; B[:, 1, 2] = v; B[:, 1, 3] = z*z*v
    Bp = np.zeros_like(B); Bp[:, 0, 0] = f[:, 1]
    Bp[:, 1, 1] = f[:, 0]; Bp[:, 1, 2] = f[:, 1]
    Bp[:, 1, 3] = 2*z*f[:, 0]*v+z*z*f[:, 1]
    Sp = np.einsum('nij,njk->nik', J, S)+B
    Spp = np.einsum('nij,njk->nik', Jp, S)+np.einsum('nij,njk->nik', J, Sp)+Bp
    density = (S-tau*tau*Spp)/(2*tau*np.asarray(sigma)[None, :, None]**2)
    terminal = (S[-1]+tau*Sp[-1])/(2*np.asarray(sigma)[:, None]**2)
    atom0 = (S[0]-tau*Sp[0])/(2*np.asarray(sigma)[:, None]**2)
    js, ls = CubicSpline(t, J), CubicSpline(t, density)
    def rhs(s, p):
        return (-js(s).T@p.reshape(2, 4)-ls(s)).ravel()
    sol = solve_ivp(rhs, (float(t[-1]), float(t[0])), terminal.ravel(),
                    t_eval=t[::-1], method='DOP853', rtol=1e-10, atol=1e-10)
    if not sol.success: raise RuntimeError('Adjoint integration failed')
    phi = sol.y.T[::-1].reshape(len(t), 2, 4)
    derivative = -np.einsum('nji,njk->nik', J, phi)-density
    # Invertible scaling avoids large score units; leaves the moment span intact.
    scale = np.maximum(np.max(np.abs(phi), axis=(0, 1)), 1.)
    return x, phi/scale, derivative/scale, atom0/scale


def moments(t, y, phi, derivative, initial, sigma):
    w = quadrature_weights(t); z, v = y.T
    A = np.column_stack(((phi[:, 0]*w[:, None]).T@v,
                         (phi[:, 1]*w[:, None]).T@z,
                         (phi[:, 1]*w[:, None]).T@v,
                         (phi[:, 1]*w[:, None]).T@((z*z-sigma[0]**2)*v)))
    b = -np.einsum('nik,ni,n->k', derivative, y, w)
    b += phi[-1].T@y[-1]-phi[0].T@np.asarray(initial)
    return np.vstack((A, np.zeros((2, 4)))), np.r_[b, y[0]-initial]


def covariance(t, model, phi, derivative, q, sigma, tau):
    w = quadrature_weights(t); z, v = model.T; a, b, c, d = q
    K = phi.shape[-1]
    U = np.zeros((len(t), 2, K+2))
    U[:, 0, :K] = (-derivative[:, 0]-(b+2*d*z*v)[:, None]*phi[:, 1])*w[:, None]
    U[:, 1, :K] = (-derivative[:, 1]-a*phi[:, 0]-(c+d*z*z)[:, None]*phi[:, 1])*w[:, None]
    U[-1, :, :K] += phi[-1]
    U[0, :, K:] = np.eye(2)
    rho = np.exp(-(t[1]-t[0])/tau)
    C = sum(sigma[j]**2*(U[:, j].T@exponential_apply(U[:, j], rho)) for j in range(2))
    # Hermite multi-indices (2,0), (1,1), (2,1); endpoint controls are linear.
    for coefficient, factor, order in ((-d*v, 2*sigma[0]**4, 2),
                                       (-2*d*z, sigma[0]**2*sigma[1]**2, 2),
                                       (-d*np.ones(len(t)), 2*sigma[0]**4*sigma[1]**2, 3)):
        H = np.zeros((len(t), K+2)); H[:, :K] = phi[:, 1]*(w*coefficient)[:, None]
        C += factor*(H.T@exponential_apply(H, rho**order))
    return (C+C.T)/2


def fit(t, y, qpilot, initial, sigma, tau, moment_sigma=None):
    model, phi, derivative, atom0 = optimal_tests(t, qpilot, initial, sigma, tau)
    A, b = moments(t, y, phi, derivative, initial, sigma if moment_sigma is None else moment_sigma)
    C = covariance(t, model, phi, derivative, qpilot, sigma, tau)
    chol = cho_factor(C, lower=True)
    WA = cho_solve(chol, A)
    V = np.linalg.inv(A.T@WA)
    estimate = V@WA.T@b
    return estimate, V, {'moment_covariance_condition': float(np.linalg.cond(C)),
                         'initial_atom_norm': float(np.linalg.norm(atom0))}
