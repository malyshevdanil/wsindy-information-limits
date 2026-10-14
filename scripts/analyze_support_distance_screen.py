"""Audit and visualize the deterministic support-distance filter."""
import json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    stage = 'support_distance'; folder = ROOT/f'results/kaggle_{stage}'; folder.mkdir(exist_ok=True)
    with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
        assert z.testzip() is None
        assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
        z.extractall(folder)
    from scripts.collect_stopped_runs import snapshot
    manifest = json.loads((folder/'manifest.json').read_text())
    audit = snapshot(stage, manifest)
    completion = json.loads((folder/'completion.json').read_text()); assert completion['status'] == 'complete'
    p = manifest['protocol']; rows = json.loads((folder/'rows.json').read_text())
    assert [r['delta'] for r in rows] == p['delta_grid']
    assert all(np.isclose(r['equal_prior_error'], norm.cdf(-r['path_distance']/2)) for r in rows)
    assert all(np.isfinite(r['path_distance']) and r['path_distance'] >= 0 for r in rows)
    successful = [i for i, r in enumerate(rows) if r['equal_prior_error'] <= p['target_equal_prior_error']]
    first = successful[0] if successful else None
    report = {'completed': True, 'snapshot': audit, 'protocol': p, 'completion': completion, 'rows': rows,
              'first_grid_crossing': None if first is None else
                  [rows[max(0, first-1)]['delta'], rows[first]['delta']],
              'caveats': ['Direct numerical Gaussian-path distances, with ODE and quadrature tolerances; not symbolic exact numbers.',
                          'Simple hypotheses with fixed other coefficients and fixed equal initial states.',
                          'Anti-damped null has a growing trajectory; no extrapolation from original limit-cycle sensitivities.',
                          'No support-recovery Monte Carlo and no composite weak-class detection theorem completed.']}
    (ROOT/'reports/support_distance_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    positive = [r for r in rows if r['delta'] > 0]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    axes[0].loglog([r['delta'] for r in positive], [r['path_distance'] for r in positive], 'o-')
    axes[0].axhline(completion['target_distance'], color='black', linestyle='--', label='5% simple-testing error')
    axes[0].set(xlabel='Magnitude of cubic coefficient', ylabel='OU Cameron–Martin path distance')
    axes[0].legend(fontsize=8)
    axes[1].semilogx([r['delta'] for r in positive], [r['equal_prior_error'] for r in positive], 'o-')
    axes[1].axhline(.05, color='black', linestyle='--')
    axes[1].set(xlabel='Magnitude of cubic coefficient', ylabel='Optimal equal-prior testing error', ylim=(-.01, .51))
    fig.suptitle('Fixed other coefficients; anti-damped linear null; T=20')
    fig.tight_layout(); dest = ROOT/'reports/figures'; dest.mkdir(exist_ok=True)
    for ext in ('png', 'pdf'): fig.savefig(dest/f'support_distance_screen.{ext}', dpi=200)
    print(json.dumps({'first_grid_crossing': report['first_grid_crossing'], 'completion': completion}, indent=2))


if __name__ == '__main__': main()
