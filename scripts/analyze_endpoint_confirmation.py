"""Full matched held-out comparison; refuse partial/selected Monte Carlo rows."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    from scripts.collect_stopped_runs import snapshot
    p = json.loads((ROOT/'configs/boundary_confirmation.json').read_text())
    rows = []; audits = []; bases = []
    refpath = ROOT/'configs'/p['references_file']; refsha = hashlib.sha256(refpath.read_bytes()).hexdigest()
    for shard in range(p['shards']):
        stage = f'boundary_confirm_{shard}'; folder = ROOT/f'results/kaggle_{stage}'
        manifest = json.loads((folder/'manifest.json').read_text())
        assert manifest['protocol'] == p and manifest['reference_file_sha256'] == refsha
        audits.append(snapshot(stage, manifest))
        assert json.loads((folder/'completion.json').read_text())['status'] == 'complete'
        assert json.loads((folder/'quadrature_gate.json').read_text())['relative_bias'] <= p['clean_mean_quadrature_bias_max']
        bases.append(json.loads((folder/'basis.json').read_text()))
        rows.extend(json.loads((folder/'rows.json').read_text()))
    seeds = list(range(p['first_seed_index'], p['last_seed_index']+1))
    methods = ['model_adjoint', 'universal_boundary', 'boundary_free_480', 'known_initial_scoring']
    indexed = {method: {r['seed']: r for r in rows if r['method'] == method} for method in methods}
    assert len(rows) == len(seeds)*len(methods)
    for method in methods: assert sorted(indexed[method]) == seeds
    assert all(len({indexed[m][seed]['input_sha256'] for m in methods}) == 1 for seed in seeds)
    rng = np.random.default_rng(p['bootstrap_seed'])
    bootstrap_indices = rng.integers(0, len(seeds), (p['bootstrap_draws'], len(seeds)))
    groups = []; boot = {}; rmses = {}
    truth = np.array([1., -1., 1., -1.])
    for method in methods:
        sample = [indexed[method][s] for s in seeds]
        errors = np.array([r['coefficients'] for r in sample])-truth
        assert np.isfinite(errors).all()
        loss = np.sum(errors*errors, axis=1)/4
        assert np.allclose(loss, np.array([r['relative_error'] for r in sample])**2, rtol=1e-12)
        rmses[method] = float(np.sqrt(loss.mean()))
        boot[method] = np.sqrt(loss[bootstrap_indices].mean(axis=1))
        group = {'method': method, 'paths': len(seeds), 'relative_rmse': rmses[method],
            'relative_bias': float(np.linalg.norm(errors.mean(axis=0))/2),
            'relative_dispersion': float(np.sqrt(np.mean(np.sum((errors-errors.mean(axis=0))**2, axis=1)))/2),
            'rmse_bootstrap_ci95': np.quantile(boot[method], [.025, .975]).tolist()}
        if 'seconds' in sample[0]:
            group['median_online_seconds_with_pilot'] = float(np.median([r['seconds']+r['pilot_seconds'] for r in sample]))
        groups.append(group)
    ratio = rmses['model_adjoint']/rmses['known_initial_scoring']
    gap = rmses['boundary_free_480']-rmses['known_initial_scoring']
    closed = None if gap <= 0 else (rmses['boundary_free_480']-rmses['universal_boundary'])/gap
    bootgap = boot['boundary_free_480']-boot['known_initial_scoring']; valid = bootgap > 0
    bootclosed = (boot['boundary_free_480'][valid]-boot['universal_boundary'][valid])/bootgap[valid]
    criteria = {'model_adjoint_ratio': ratio, 'model_adjoint_pass': ratio <= p['criteria']['model_adjoint_rmse_over_known_initial_likelihood_max'],
        'universal_fraction_gap_closed': closed, 'universal_pass': None if closed is None else closed >= p['criteria']['universal_fraction_of_rmse_gap_closed_min'],
        'model_adjoint_ratio_ci95': np.quantile(boot['model_adjoint']/boot['known_initial_scoring'], [.025, .975]).tolist(),
        'gap_closed_ci95_conditional_positive_denominator': np.quantile(bootclosed, [.025, .975]).tolist() if len(bootclosed) else None,
        'bootstrap_nonpositive_gap_count': int((~valid).sum())}
    report = {'completed': True, 'protocol': p, 'audits': audits, 'paired_input_fingerprints_verified': True,
              'reference_file_sha256': refsha, 'groups': groups, 'criteria': criteria, 'bases': bases,
              'limitations': ['Model-adjoint branch realizes efficient score directions in first order.',
                 'Universal windows depend on time alone; covariance pilot still uses ODE integration.',
                 'A failure of this finite family cannot rule out recovery by a richer universal family.']}
    (ROOT/'reports/endpoint_confirmation_results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({'groups': groups, 'criteria': criteria}, indent=2))


if __name__ == '__main__': main()
