import numpy as np
from scipy.signal import lfilter


def stationary_ar1(n, sigma, rho, rng):
    """Exact stationary Gaussian AR(1), independent state components."""
    sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
    if n < 2 or np.any(sigma < 0) or not 0 <= rho < 1:
        raise ValueError('Need n>=2, nonnegative sigma, and 0<=rho<1')
    out = np.empty((n, len(sigma)))
    out[0] = rng.normal(size=len(sigma)) * sigma
    innovation = sigma * np.sqrt((1 - rho) * (1 + rho))
    for j in range(len(sigma)):
        out[1:, j], _ = lfilter([1.], [1., -rho],
                               rng.normal(size=n - 1) * innovation[j],
                               zi=[rho * out[0, j]])
    return out


def ou_noise(n, h, tau, sigma, rng):
    if h <= 0 or tau <= 0:
        raise ValueError('h and tau must be positive')
    # expm1 avoids cancellation in the tiny-step innovation variance.
    sigma = np.atleast_1d(np.asarray(sigma, dtype=float))
    if n < 2 or np.any(sigma < 0):
        raise ValueError('Invalid n or sigma')
    rho = np.exp(-h / tau)
    out = np.empty((n, len(sigma)))
    out[0] = rng.normal(size=len(sigma)) * sigma
    innovation = sigma * np.sqrt(-np.expm1(-2 * h / tau))
    for j in range(len(sigma)):
        out[1:, j], _ = lfilter([1.], [1., -rho],
                               rng.normal(size=n - 1) * innovation[j],
                               zi=[rho * out[0, j]])
    return out


def calibrated_ar1(noise, h):
    """Independent zero-mean sensor calibration; not a residual-based fit.

    Moment estimates have finite-record error. No nugget separation is attempted.
    """
    a = np.asarray(noise, float)
    variance = np.mean(a * a, axis=0)
    rho_components = np.sum(a[:-1] * a[1:], axis=0) / np.maximum(
        np.sum(a[:-1] ** 2, axis=0), np.finfo(float).tiny)
    rho = float(np.clip(np.mean(rho_components), 1e-8, 1 - 1e-8))
    return {'variance': variance, 'rho': rho, 'tau': -h / np.log(rho),
            'assumption': 'zero_mean_independent_components_pure_AR1'}
