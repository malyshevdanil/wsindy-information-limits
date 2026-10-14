"""Added-cubic information references and author WSINDy inclusion frequencies."""
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def wilson(successes, n):
    z = norm.ppf(.975); rate = successes/n
    center = (rate+z*z/(2*n))/(1+z*z/n)
    radius = z*np.sqrt(rate*(1-rate)/n+z*z/(4*n*n))/(1+z*z/n)
    return [float(center-radius), float(center+radius)]


def main():
    from scripts.collect_stopped_runs import snapshot
    stage = 'support_extra'; folder = ROOT/f'results/kaggle_{stage}'
    manifest = json.loads((folder/'manifest.json').read_text()); p = manifest['protocol']
    audit = snapshot(stage, manifest)
    assert json.loads((folder/'completion.json').read_text())['status'] == 'complete'
    information = json.loads((folder/'information.json').read_text())
    references = [r for r in information['rows'] if r['h'] == p['information_refinement_h'] and r['cutoff'] == p['primary_svd_cutoff']]
    rows = json.loads((folder/'rows.json').read_text()); grids = json.loads((folder/'delta_grids.json').read_text())
    assert len(rows) == len(p['levels'])*len(p['delta_multipliers_of_weak_reference'])*p['paths_per_delta_per_level']
    groups = []
    for level in p['levels']:
        for di, delta in enumerate(grids[str(level)]):
            sample = [r for r in rows if r['level'] == level and r['delta_index'] == di]
            assert sorted(r['seed'] for r in sample) == list(range(p['paths_per_delta_per_level']))
            assert all(r['delta'] == delta and r['status'] == 'ok' for r in sample)
            count = sum(r['included'] for r in sample)
            groups.append({'level': level, 'delta_index': di, 'delta': delta,
                'magnitude_over_local_weak_reference': abs(p['delta_multipliers_of_weak_reference'][di]),
                'paths': len(sample), 'included': count, 'inclusion_fraction': count/len(sample),
                'wilson_ci95': wilson(count, len(sample)),
                'median_fit_seconds': float(np.median([r['fit_seconds'] for r in sample]))})
    report = {'completed': True, 'protocol': p, 'audit': audit, 'information_references': references, 'groups': groups,
        'author_provenance': json.loads((folder/'author_provenance.json').read_text()),
        'limitations': ['Local Gaussian two-point reference; not an exact finite-perturbation support guarantee.',
            'Full author adaptive selection has no prescribed 5% false-inclusion calibration here; report the null rate.',
            'Information uses matched dictionary nuisance coefficients, not the unprofiled diagonal of a Fisher matrix.',
            'Negative restoring cubic alternatives; no conclusion about unstable large positive perturbations.']}
    (ROOT/'reports/support_extra_term_results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(p['levels']), figsize=(9, 3.5), sharey=True)
    for ax, level in zip(axes, p['levels']):
        sample = [g for g in groups if g['level'] == level]
        x = np.array([g['magnitude_over_local_weak_reference'] for g in sample])
        y = np.array([g['inclusion_fraction'] for g in sample]); ci = np.array([g['wilson_ci95'] for g in sample])
        ax.errorbar(x, y, yerr=np.maximum(0., np.vstack((y-ci[:, 0], ci[:, 1]-y))), fmt='o-', label='Author WSINDy; Wilson 95%')
        ax.axvline(1., color='black', linestyle='--', label='Local weak reference')
        known = next(r for r in references if r['level'] == level and r['label'] == 'full_library_known_initial')
        weak = next(r for r in references if r['level'] == level and r['label'] == 'full_library_boundary_free')
        ax.axvline(known['delta_min_local']/weak['delta_min_local'], color='gray', linestyle=':', label='Known-initial likelihood reference')
        ax.set(title=f'{level*100:g}% noise', xlabel='Magnitude / local weak detection reference', ylim=(-.02, 1.02))
    axes[0].set_ylabel('Fraction including x³ in v′'); axes[-1].legend(fontsize=7)
    fig.tight_layout(); destination = ROOT/'reports/figures'; destination.mkdir(exist_ok=True)
    for ext in ('pdf', 'png'): fig.savefig(destination/f'support_extra_inclusion.{ext}', dpi=200)
    print(json.dumps({'information_references': references, 'groups': groups}, indent=2))


if __name__ == '__main__': main()
