"""Paired Monte Carlo uncertainty for the frozen 32-path boundary screen."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    stage = 'boundary_screen_v4'
    report = json.loads((ROOT/f'reports/{stage}_validation.json').read_text())
    assert report['paired_inputs_verified'] and report['completion']['paths'] == 32
    rows = json.loads((ROOT/f'results/kaggle_{stage}/rows.json').read_text())
    new = {r['seed']: r for r in rows if r['method'] == 'feasible_boundary_weak'}
    known = {r['seed']: r for r in json.loads((ROOT/'results/kaggle_efficiency_fisher_v2/rows.json').read_text())
             if r['method'] == 'fisher_known_initial'}
    weak = {}
    for shard in range(4):
        weak.update({r['seed']: r for r in json.loads((ROOT/f'results/kaggle_efficiency_power_{shard}/rows.json').read_text())
                     if r['weight_model'] == 'ou'})
    seeds = sorted(new)
    assert all(new[s]['input_sha256'] == known[s]['input_sha256'] == weak[s]['input_sha256'] for s in seeds)
    arrays = {label: np.array([source[s]['relative_error']**2 for s in seeds])
              for label, source in [('endpoint_weak', new), ('known_initial_scoring', known), ('boundary_free_480', weak)]}
    # Fixed postprocessing seed; bootstrap quantifies uncertainty, not a new
    # stopping/acceptance rule for the already frozen point-estimate threshold.
    rng = np.random.default_rng(8400000); indices = rng.integers(0, len(seeds), (20000, len(seeds)))
    groups = []
    boot = {}
    for label, errors in arrays.items():
        boot[label] = np.sqrt(errors[indices].mean(axis=1))
        groups.append({'method': label, 'paths': len(seeds), 'relative_rmse': float(np.sqrt(errors.mean())),
                       'bootstrap_rmse_ci95': np.quantile(boot[label], [.025, .975]).tolist()})
    differences = {label: np.quantile(boot['endpoint_weak']-boot[label], [.025, .975]).tolist()
                   for label in ('known_initial_scoring', 'boundary_free_480')}
    report['paired_prefix_comparison'] = {'groups': groups, 'difference_rmse_ci95': differences,
        'bootstrap_draws': 20000, 'bootstrap_seed': 8400000,
        'threshold_decision_unchanged': '32-path point RMSE exceeds frozen .00140 screen threshold; no full run triggered'}
    report['online_seconds'] = {'median_three_updates': float(np.median([r['seconds'] for r in new.values()])),
                                'scope': 'Three directly integrated adjoint updates; preliminary 30-window pilot excluded.'}
    (ROOT/f'reports/{stage}_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report['paired_prefix_comparison'], indent=2))


if __name__ == '__main__': main()
