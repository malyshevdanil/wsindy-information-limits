"""Local Schur-complement information for an added x^3 term in Van der Pol."""
import numpy as np
from scipy.integrate import solve_ivp
from .path_scoring import ou_whiten

POWERS = np.array([[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]])
TARGET = 14  # Equation-major coefficients: 9 in x', 9 in v'; x^3 is index 5.


def augmented_sensitivities(t, initial):
    def rhs(_, state):
        z, v = state[:2]; S = state[2:].reshape(2, 20)
        J = np.array([[0., 1.], [-1.-2*z*v, 1.-z*z]])
        features = np.prod(np.array([z, v])[None]**POWERS, axis=1)
        B = np.zeros_like(S); B[0, :9] = features; B[1, 9:18] = features
        return np.r_[v, -z+v-z*z*v, (J@S+B).ravel()]
    initial_state = np.zeros(42); initial_state[:2] = initial
    initial_state[2:].reshape(2, 20)[:, 18:] = np.eye(2)
    sol = solve_ivp(rhs, (0., float(t[-1])), initial_state, t_eval=t,
                    method='DOP853', rtol=1e-10, atol=1e-12)
    if not sol.success: raise RuntimeError('Added-term sensitivity integration failed')
    return sol.y[:2].T, sol.y[2:].T.reshape(len(t), 2, 20)


def efficient_column_information(W, target, nuisance, cutoff):
    """Residualize a score column; avoid subtraction of nearly equal Gramians."""
    A = W[:, nuisance]
    scales = np.linalg.norm(A, axis=0)
    if np.any(scales == 0): raise ValueError('Zero nuisance score column')
    U, values, _ = np.linalg.svd(A/scales, full_matrices=False)
    keep = values > cutoff*values.max(); Q = U[:, keep]
    residual = W[:, target]-Q@(Q.T@W[:, target])
    information = float(residual@residual)
    if not information > 0: raise ValueError('Nonpositive effective information')
    return {'information': information, 'standard_deviation': information**-.5,
            'nuisance_rank': int(keep.sum()), 'nuisance_columns': len(nuisance),
            'nuisance_minimum_relative_singular_value': float(values[-1]/values[0])}


def references(t, S, sigma, tau, cutoff):
    W = ou_whiten(S, t[1]-t[0], tau, sigma).reshape(-1, 20)
    core = [1, 9, 10, 15]
    full = [j for j in range(18) if j != TARGET]
    return {label: efficient_column_information(W, TARGET, nuisance, cutoff)
            for label, nuisance in [('full_library_known_initial', full),
                                    ('full_library_boundary_free', full+[18, 19]),
                                    ('core5_known_initial', core),
                                    ('core5_boundary_free', core+[18, 19])]}
