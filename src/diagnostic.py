import numpy as np


def choose_grid(records, reference_variance, eps_p=.05, eps_q=.01, eps_validation=.05,
                variance_floor=1e-12):
    """Frozen-support pilot diagnostic. Reference is finite refined-grid proxy.

    records contain h, variance, quadrature_difference, validation_loss,
    support_stable, conditioning_ok. Never interpret this as a sensor bandwidth.
    """
    valid = [r for r in records if np.isfinite(r['validation_loss'])]
    if not valid or not np.isfinite(reference_variance):
        return {'status': 'saturation_not_established', 'h': None}
    best_validation = min(r['validation_loss'] for r in valid)
    candidates = []
    for r in records:
        close = abs(r['variance'] - reference_variance) / max(reference_variance, variance_floor)
        if (r['support_stable'] and r['conditioning_ok'] and np.isfinite(close)
                and close <= eps_p and r['quadrature_difference'] <= eps_q
                and r['validation_loss'] <= (1 + eps_validation) * best_validation + 1e-12):
            candidates.append(r)
    if not candidates:
        return {'status': 'saturation_not_established', 'h': None}
    chosen = max(candidates, key=lambda r: r['h'])
    return {'status': 'candidate_selected', 'h': chosen['h'],
            'reference': 'finite_refined_grid_proxy',
            'relative_variance_gap': abs(chosen['variance'] - reference_variance) /
                                     max(reference_variance, variance_floor)}
