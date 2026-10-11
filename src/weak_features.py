from itertools import product
from math import factorial
import numpy as np


def exponents(d, degree):
    if d < 1 or degree < 0:
        raise ValueError('Invalid dimension or degree')
    return sorted((a for a in product(range(degree + 1), repeat=d)
                   if sum(a) <= degree), key=lambda a: (sum(a), tuple(-v for v in a)))


def hermite(u, variance, degree):
    values = [np.ones_like(u)]
    if degree:
        values.append(u)
    for q in range(1, degree):
        values.append(u * values[-1] - q * variance * values[-2])
    return values


def library(y, powers, correction_variance=None):
    y = np.asarray(y, float)
    if y.ndim != 2:
        raise ValueError('States must be an N by d array')
    degree = max(map(sum, powers))
    variance = (np.zeros(y.shape[1]) if correction_variance is None
                else np.broadcast_to(correction_variance, (y.shape[1],)))
    polynomials = [hermite(y[:, j], variance[j], degree)
                   for j in range(y.shape[1])]
    return np.column_stack([np.prod([polynomials[j][q] for j, q in enumerate(a)], axis=0)
                            for a in powers])


def quadrature_weights(t):
    t = np.asarray(t, float)
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError('Need increasing time grid')
    w = np.empty_like(t)
    w[0], w[-1] = (t[1] - t[0]) / 2, (t[-1] - t[-2]) / 2
    w[1:-1] = (t[2:] - t[:-2]) / 2
    return w


def test_functions(t, centers, width, order=4):
    """C^(order-1) compact polynomial bumps; width is the full support length."""
    if width <= 0 or order < 3:
        raise ValueError('Need positive width and order>=3')
    t = np.asarray(t)
    z = (t[None, :] - np.asarray(centers)[:, None]) / (width / 2)
    interior = np.abs(z) < 1
    base = np.maximum(1 - z * z, 0)
    phi = base ** order
    derivative = np.where(interior, -4 * order * z / width * base ** (order - 1), 0)
    return phi, derivative


def weak_matrices(t, y, powers, centers, width, correction_variance=None):
    phi, derivative = test_functions(t, centers, width)
    w = quadrature_weights(t)
    G = (phi * w) @ library(y, powers, correction_variance)
    b = -(derivative * w) @ np.asarray(y)
    return G, b


def derivative_polynomial(x, powers, coefficients, alpha):
    """Evaluate partial^alpha f / alpha!, equation columns in coefficients."""
    out = np.zeros((len(x), coefficients.shape[1]))
    for beta, coef in zip(powers, coefficients):
        if all(b >= a for b, a in zip(beta, alpha)):
            factor = np.prod([factorial(b) / (factorial(a) * factorial(b - a))
                              for b, a in zip(beta, alpha)])
            term = np.prod(np.asarray(x) ** np.array([b - a for b, a in zip(beta, alpha)]), axis=1)
            out += factor * term[:, None] * coef
    return out


def heat_transform(powers, coefficients, variance_difference):
    """Finite Gaussian heat operator: raw and imperfect Hermite corrections."""
    delta = np.asarray(variance_difference)
    lookup = {a: i for i, a in enumerate(powers)}
    out = np.zeros_like(coefficients)
    for beta, coef in zip(powers, coefficients):
        for m in product(*(range(b // 2 + 1) for b in beta)):
            alpha = tuple(b - 2 * q for b, q in zip(beta, m))
            factor = np.prod([factorial(b) / (factorial(b - 2 * q) * factorial(q))
                              * (delta[j] / 2) ** q
                              for j, (b, q) in enumerate(zip(beta, m))])
            out[lookup[alpha]] += factor * coef
    return out
