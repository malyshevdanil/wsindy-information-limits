import numpy as np


def stls(G, b, threshold, max_iterations=30):
    """Simple normalized sequential threshold least squares, not MSTLS."""
    if threshold < 0:
        raise ValueError('Negative threshold')
    b = np.asarray(b)
    if b.ndim == 1:
        b = b[:, None]
    norms = np.linalg.norm(G, axis=0)
    valid = norms > 1e-12 * max(float(norms.max()), 1.)
    safe_norms = np.where(valid, norms, 1.)
    A = G / safe_norms
    scale_b = np.maximum(np.linalg.norm(b, axis=0), 1e-14)
    result = np.zeros((G.shape[1], b.shape[1]))
    for i in range(b.shape[1]):
        active = valid.copy()
        z = b[:, i] / scale_b[i]
        gamma = np.zeros(G.shape[1])
        for _ in range(max_iterations):
            gamma[:] = 0
            if not active.any():
                break
            gamma[active] = np.linalg.lstsq(A[:, active], z, rcond=None)[0]
            new_active = active & (np.abs(gamma) >= threshold)
            if np.array_equal(active, new_active):
                break
            active = new_active
        # Refit even if iteration budget was exhausted; rejected terms stay zero.
        gamma[:] = 0
        if active.any():
            gamma[active] = np.linalg.lstsq(A[:, active], z, rcond=None)[0]
        result[:, i] = gamma * scale_b[i] / safe_norms
    return result


def select_threshold(G, b, validation, candidates, complexity_penalty=1e-3):
    records = []
    for threshold in candidates:
        coef = stls(G, b, threshold)
        loss = np.mean([np.linalg.norm(bv - Gv @ coef) ** 2 /
                        max(np.linalg.norm(bv) ** 2, 1e-20) for Gv, bv in validation])
        loss += complexity_penalty * np.count_nonzero(coef) / coef.size
        records.append((float(loss), float(threshold), coef))
    loss, threshold, coef = min(records, key=lambda r: (r[0], -r[1]))
    return threshold, coef, loss
