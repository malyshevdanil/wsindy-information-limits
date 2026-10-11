"""Create tables/figures only from an actual empirical run, retaining failures."""
import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def wilson(success, total, z=1.959963984540054):
    if not total:
        return None, None
    p = success / total
    center = (p + z*z/(2*total)) / (1 + z*z/total)
    radius = z*np.sqrt(p*(1-p)/total + z*z/(4*total*total)) / (1 + z*z/total)
    return float(center-radius), float(center+radius)


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(folder):
    manifest = json.loads((folder / 'manifest.json').read_text())
    if manifest['inference_status'] == 'installation_check_not_research':
        print('Smoke run: research figures and inferential summaries intentionally omitted.')
        return
    latest = {}
    for line in (folder / 'runs.jsonl').read_text().splitlines():
        row = json.loads(line)
        latest[row['key']] = row
    grouped = defaultdict(list)
    fields = ('system', 'regime', 'tau_over_width', 'noise_level', 'method', 'h')
    for row in latest.values():
        grouped[tuple(row[k] for k in fields)].append(row)
    summaries = []
    for key, rows in sorted(grouped.items()):
        good = [r for r in rows if r['status'] == 'ok']
        errors = [r['coefficient_error'] for r in good if r['coefficient_error'] is not None]
        n = len(rows)
        support_count = sum(r['exact_support'] for r in good)
        low, high = wilson(support_count, n)
        result = dict(zip(fields, key))
        result.update({'attempts': n, 'fit_failures': n-len(good),
                       'support_success_fraction_failure_counted': support_count/n,
                       'support_wilson_low': low, 'support_wilson_high': high,
                       'coefficient_error_median_successful_fits': float(np.median(errors)) if errors else None,
                       'coefficient_error_q25': float(np.quantile(errors, .25)) if errors else None,
                       'coefficient_error_q75': float(np.quantile(errors, .75)) if errors else None,
                       'forecast_failures': sum(r.get('forecast_status') != 'ok' for r in good),
                       'fit_seconds_median': float(np.median([r['fit_seconds'] for r in good])) if good else None})
        if len(good) >= 2:
            arrays = [np.load(folder / r['covariance_file']) for r in good]
            centered = np.array([a['residual'] - a['expected_mean'] for a in arrays])
            # Known-mean empirical second moment, independent seeds are the replicates.
            empirical = centered.T @ centered / len(good)
            predicted = np.mean([a['exact_covariance'] for a in arrays], axis=0)
            first = np.mean([a['first_order_covariance'] for a in arrays], axis=0)
            result['exact_covariance_relative_frobenius_error'] = float(
                np.linalg.norm(empirical-predicted) / max(np.linalg.norm(predicted), 1e-20))
            result['first_order_covariance_relative_frobenius_error'] = float(
                np.linalg.norm(empirical-first) / max(np.linalg.norm(first), 1e-20))
            result['residual_mean_error_norm'] = float(np.linalg.norm(centered.mean(axis=0)))
            for a in arrays:
                a.close()
        else:
            result.update({'exact_covariance_relative_frobenius_error': None,
                           'first_order_covariance_relative_frobenius_error': None,
                           'residual_mean_error_norm': None})
        summaries.append(result)
    write_csv(folder / 'summary.csv', summaries)
    pairs = []
    base_cells = defaultdict(dict)
    for key, rows in grouped.items():
        base_cells[key[:-1]][key[-1]] = rows
    rng = np.random.default_rng(98421)
    for cell, grids in sorted(base_cells.items()):
        steps = sorted(grids, reverse=True)
        for coarse, fine in zip(steps, steps[1:]):
            a = {r['seed']: r for r in grids[coarse] if r['status'] == 'ok'}
            b = {r['seed']: r for r in grids[fine] if r['status'] == 'ok'}
            seeds = [s for s in sorted(set(a)&set(b)) if a[s]['coefficient_error'] is not None
                     and b[s]['coefficient_error'] is not None]
            if len(seeds) < 2:
                continue
            ea = np.array([a[s]['coefficient_error'] for s in seeds])
            eb = np.array([b[s]['coefficient_error'] for s in seeds])
            def improvement(indices):
                return 1 - np.median(eb[indices], axis=-1) / np.maximum(np.median(ea[indices], axis=-1), 1e-15)
            bootstrap = improvement(rng.integers(0, len(seeds), size=(2000, len(seeds))))
            ci = np.quantile(bootstrap, [.025, .975])
            result = dict(zip(fields[:-1], cell))
            result.update({'h_coarse': coarse, 'h_fine': fine, 'paired_seeds': len(seeds),
                           'relative_median_improvement': float(improvement(np.arange(len(seeds)))),
                           'bootstrap_low': float(ci[0]), 'bootstrap_high': float(ci[1]),
                           'equivalence_within_5pct_on_this_pair': bool(ci[0] >= -.05 and ci[1] <= .05),
                           'excluded_unpaired_or_failed': len(grids[coarse])+len(grids[fine])-2*len(seeds),
                           'pairing_meaning': 'same_record' if cell[1] != 'fixed_index' else 'seed_block_different_process'})
            pairs.append(result)
    write_csv(folder / 'paired_grid_comparisons.csv', pairs)
    figures = folder / 'figures'
    figures.mkdir(exist_ok=True)
    for name in manifest['config']['systems']:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        for regime, marker in [('physical_ou', 'o'), ('iid', 's'), ('fixed_index', '^')]:
            # Only the declared representative pilot cell; tables preserve all cells.
            tau = manifest['config']['tau_over_width'][1 if len(manifest['config']['tau_over_width']) > 1 else 0]
            noise = manifest['config']['noise_levels'][-1]
            chosen = sorted([r for r in summaries if r['system']==name and r['regime']==regime
                             and r['method']=='known_debias' and r['tau_over_width']==tau
                             and r['noise_level']==noise], key=lambda r: r['h'])
            if not chosen:
                continue
            h = [r['h'] for r in chosen]
            error = [r['coefficient_error_median_successful_fits'] for r in chosen]
            axes[0].loglog(h, error, marker+'-', label=regime)
            axes[1].semilogx(h, [r['support_success_fraction_failure_counted'] for r in chosen], marker+'-', label=regime)
        axes[0].set_ylabel('Median relative coefficient error (successful fits)')
        axes[1].set_ylabel('Exact-support frequency (failures included)')
        axes[1].set_ylim(-.02, 1.02)
        for ax in axes:
            ax.set_xlabel('Sampling step h')
            ax.grid(alpha=.2)
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize=8)
        fig.suptitle(f'{name}: pilot only; independent seeds, known Gaussian correction')
        fig.savefig(figures / f'{name}_pilot.png', dpi=220)
        fig.savefig(figures / f'{name}_pilot.pdf')
        plt.close(fig)
    (folder / 'summary_status.json').write_text(json.dumps({
        'stage': manifest['inference_status'], 'unique_attempts': len(latest),
        'failed_attempts': sum(r['status']!='ok' for r in latest.values()),
        'inference_unit': 'independent noise seed', 'bootstrap_replicates': 2000,
        'limitations': ['pilot, not confirmatory', 'no multiple-comparison correction',
                       'coefficient quantiles conditional on successful fits; failures counted separately',
                       'covariance agreement with 10 seeds is imprecise',
                       'single-pair equivalence does not establish two-successive-step saturation']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    summarize(parser.parse_args().folder)
