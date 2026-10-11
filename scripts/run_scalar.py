"""Separate exact one-window estimator experiment; execute on Kaggle."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.noise import ou_noise
from src.weak_features import test_functions, quadrature_weights


def run(output, repeats=1000):
    output.mkdir(parents=True, exist_ok=True)
    T, tau, sigma, a, x0 = 4., .4, .1, -1., 2.
    t = np.linspace(0, T, 4097)
    x = x0 * np.exp(a*t)
    phi, derivative = test_functions(t, [2.], 3.)
    records = []
    for seed in range(8000, 8000+repeats):
        e = ou_noise(len(t), t[1]-t[0], tau, [sigma], np.random.default_rng(seed))[:, 0]
        y = x+e
        for stride in (32,16,8,4,2,1):
            w = quadrature_weights(t[::stride])
            denominator = float(np.dot(phi[0, ::stride]*w, y[::stride]))
            numerator = -float(np.dot(derivative[0, ::stride]*w, y[::stride]))
            D = float(np.dot(phi[0, ::stride]*w, x[::stride]))
            estimate = numerator/denominator if denominator != 0 else None
            records.append({'seed': seed, 'h': float(t[1]*stride), 'denominator': denominator,
                            'clean_denominator': D, 'estimate': estimate,
                            'absolute_error': abs(estimate-a) if estimate is not None else None,
                            'denominator_guard_failed': abs(denominator-D) > .5*abs(D)})
    with (output/'scalar_ratio.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    summary = []
    for h in sorted({r['h'] for r in records}):
        selected = [r for r in records if r['h']==h]
        errors = [r['absolute_error'] for r in selected if r['absolute_error'] is not None]
        summary.append({'h': h, 'repeats': repeats,
                        'absolute_error_quantiles_50_90_99': np.quantile(errors, [.5,.9,.99]).tolist(),
                        'denominator_guard_failure_rate': sum(r['denominator_guard_failed'] for r in selected)/repeats})
    (output/'scalar_summary.json').write_text(json.dumps({'T': T, 'tau': tau, 'sigma': sigma,
        'a': a, 'initial': x0, 'seed_range': [8000,8000+repeats-1], 'summary': summary,
        'interpretation': 'quantile Monte Carlo experiment; no finite unconditional MSE is claimed'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=1000)
    args = parser.parse_args()
    run(args.output, args.repeats)
