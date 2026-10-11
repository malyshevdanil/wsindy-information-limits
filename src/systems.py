import numpy as np
from scipy.integrate import solve_ivp
from .weak_features import exponents, library


DEFINITIONS = {
    'linear': {'d': 1, 'degree': 1, 'initial': [2.],
               'terms': {(1,): [-1.]}},
    'damped': {'d': 2, 'degree': 2, 'initial': [2., 0.],
               'terms': {(1, 0): [0., -1.], (0, 1): [1., -.2]}},
    'vdp': {'d': 2, 'degree': 3, 'initial': [2., 0.],
            'terms': {(1, 0): [0., -1.], (0, 1): [1., 1.], (2, 1): [0., -1.]}},
    'duffing': {'d': 2, 'degree': 3, 'initial': [1.5, .5],
                'terms': {(1, 0): [0., -1.], (0, 1): [1., -.1], (3, 0): [0., -1.]}}
}


def system(name):
    info = DEFINITIONS[name]
    powers = exponents(info['d'], info['degree'])
    coef = np.zeros((len(powers), info['d']))
    for a, values in info['terms'].items():
        coef[powers.index(a)] = values
    return info, powers, coef


def trajectory(name, t, initial=None):
    info, powers, coef = system(name)
    initial = info['initial'] if initial is None else initial
    if name == 'linear':
        return np.asarray(initial)[None, :] * np.exp(-np.asarray(t)[:, None])
    def rhs(_, state):
        return (library(np.asarray(state)[None, :], powers) @ coef)[0]
    sol = solve_ivp(rhs, (float(t[0]), float(t[-1])), initial, t_eval=t,
                    method='DOP853', rtol=1e-10, atol=1e-12)
    if not sol.success or sol.y.shape[1] != len(t):
        raise RuntimeError(sol.message)
    return sol.y.T


def forecast_error(name, coefficients, horizon=2., points=201):
    info, powers, _ = system(name)
    initial = np.asarray(info['initial']) * .73 + .17
    t = np.linspace(0, horizon, points)
    truth = trajectory(name, t, initial)
    def rhs(_, state):
        return (library(np.asarray(state)[None, :], powers) @ coefficients)[0]
    def escape(_, state):
        return 100. - np.linalg.norm(state)
    escape.terminal = True
    sol = solve_ivp(rhs, (0., horizon), initial, t_eval=t, events=escape,
                    method='DOP853', rtol=1e-7, atol=1e-9)
    if not sol.success or sol.y.shape[1] != points:
        return float('inf'), 'failed_or_escaped'
    return float(np.linalg.norm(sol.y.T - truth) / np.linalg.norm(truth)), 'ok'
