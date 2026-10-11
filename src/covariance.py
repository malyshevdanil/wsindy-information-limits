from math import factorial
import numpy as np
from scipy.signal import lfilter
from .weak_features import (exponents, test_functions, quadrature_weights,
                            derivative_polynomial, heat_transform)


def exponential_apply(values, rho):
    """Return Toeplitz(rho^|i-j|) @ values without a dense N by N matrix."""
    if not 0 <= rho < 1:
        raise ValueError('Invalid rho')
    a = np.asarray(values, float)
    forward = lfilter([1.], [1., -rho], a, axis=0)
    backward = lfilter([1.], [1., -rho], a[::-1], axis=0)[::-1]
    return forward + backward - a


def residual_covariance(t, x, powers, coefficients, centers, width,
                        sigma_c, tau, sigma_w=None, correction_variance=None,
                        first_order=False):
    """Exact covariance of polynomial weak residuals (window, equation order).

    x and coefficients may be oracle truth or a explicitly labelled pilot plug-in.
    Gaussian components are independent and share tau. Calibration uncertainty
    is not marginalized. first_order means Hermite q=1 only.
    """
    t, x = np.asarray(t), np.asarray(x)
    h = t[1] - t[0]
    if tau <= 0 or not np.allclose(np.diff(t), h):
        raise ValueError('Need uniform grid and positive tau')
    d, K = x.shape[1], len(centers)
    sc = np.broadcast_to(sigma_c, (d,))
    sw = np.zeros(d) if sigma_w is None else np.broadcast_to(sigma_w, (d,))
    true_variance = sc ** 2 + sw ** 2
    correction = (true_variance if correction_variance is None
                  else np.broadcast_to(correction_variance, (d,)))
    effective_coef = heat_transform(powers, coefficients, true_variance - correction)
    degree = max(1, max(map(sum, powers)))
    phi, derivative = test_functions(t, centers, width)
    w = quadrature_weights(t)
    result = np.zeros((K * d, K * d))
    for alpha in exponents(d, 1 if first_order else degree):
        q = sum(alpha)
        if not q:
            continue
        deriv = derivative_polynomial(x, powers, effective_coef, alpha)
        c = -phi[:, :, None] * deriv[None, :, :]
        if q == 1:
            c[:, :, alpha.index(1)] -= derivative
        U = (c * w[None, :, None]).transpose(1, 0, 2).reshape(len(t), K * d)
        factor = np.prod([factorial(a) for a in alpha])
        correlated = float(np.prod((sc ** 2) ** np.array(alpha)))
        total = float(np.prod(true_variance ** np.array(alpha)))
        result += factor * (correlated * (U.T @ exponential_apply(U, np.exp(-q * h / tau)))
                            + (total - correlated) * (U.T @ U))
    return (result + result.T) / 2


def coefficient_surrogate(G, residual_cov, equation, d, support):
    """Fixed-support clean/pilot-design surrogate; not actual estimator moments."""
    rows = np.arange(equation, residual_cov.shape[0], d)
    S = np.asarray(support, bool)
    if not np.any(S) or np.linalg.matrix_rank(G[:, S]) < S.sum():
        return None
    P = np.linalg.pinv(G[:, S])
    return P @ residual_cov[np.ix_(rows, rows)] @ P.T
