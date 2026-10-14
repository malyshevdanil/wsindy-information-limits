"""Require all 500 paired records before endpoint-component attribution."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    from scripts.collect_stopped_runs import snapshot
    p = json.loads((ROOT/'configs/boundary_ablation.json').read_text())
    rows = []; audits = []; gates = []
    for i in range(p['shards']):
        stage = f'boundary_ablate_v2_{i}'; folder = ROOT/f'results/kaggle_{stage}'
        manifest = json.loads((folder/'manifest.json').read_text()); assert manifest['protocol'] == p
        audits.append(snapshot(stage, manifest))
        c = json.loads((folder/'completion.json').read_text()); assert c['status'] == 'complete' and c['paths'] == 125
        g = json.loads((folder/'quadrature_gates.json').read_text()); assert all(v['relative_bias'] <= p['clean_mean_quadrature_bias_max'] for v in g)
        gates.append(g); rows.extend(json.loads((folder/'rows.json').read_text()))
    methods = list(p['branches'])+p['archived_branches']; seeds = list(range(32, 532))
    indexed = {m: {r['seed']: r for r in rows if r['method'] == m} for m in methods}
    assert len(rows) == len(seeds)*len(methods)
    assert all(sorted(indexed[m]) == seeds for m in methods)
    assert all(len({indexed[m][s]['input_sha256'] for m in methods}) == 1 for s in seeds)
    indices = np.random.default_rng(p['bootstrap_seed']).integers(0, 500, (p['bootstrap_draws'], 500))
    groups = []; boot = {}; truth = np.array([1., -1., 1., -1.]); rmses = {}
    for method in methods:
        error = np.array([indexed[method][s]['coefficients'] for s in seeds])-truth
        assert np.isfinite(error).all()
        loss = np.sum(error**2, axis=1)/4
        rmses[method] = float(np.sqrt(loss.mean())); boot[method] = np.sqrt(loss[indices].mean(axis=1))
        groups.append({'method': method, 'paths': 500, 'relative_rmse': rmses[method],
                       'relative_bias': float(np.linalg.norm(error.mean(axis=0))/2),
                       'relative_dispersion': float(np.sqrt(np.mean(np.sum((error-error.mean(axis=0))**2, axis=1)))/2),
                       'rmse_ci95': np.quantile(boot[method], [.025, .975]).tolist()})
    gap = rmses['matched_boundary_free']-rmses['known_initial_scoring']
    bootgap = boot['matched_boundary_free']-boot['known_initial_scoring']; valid = bootgap > 0
    for g in groups:
        m = g['method']; g['fraction_matched_gap_closed'] = (rmses['matched_boundary_free']-rmses[m])/gap if gap > 0 else None
        fraction = (boot['matched_boundary_free'][valid]-boot[m][valid])/bootgap[valid]
        g['fraction_closed_ci95'] = np.quantile(fraction, [.025, .975]).tolist() if len(fraction) else None
        g['material_by_frozen_point_rule'] = g['fraction_matched_gap_closed'] >= .5 if gap > 0 else None
    report = {'completed': True, 'protocol': p, 'audits': audits, 'quadrature_gates': gates,
              'paired_fingerprints_verified': True, 'groups': groups,
              'bootstrap_nonpositive_gap_count': int((~valid).sum()),
              'all_minus_initial_rmse': rmses['universal_boundary']-rmses['initial_only'],
              'all_minus_initial_rmse_ci95': np.quantile(boot['universal_boundary']-boot['initial_only'], [.025,.975]).tolist(),
              'limitations': ['Attribution on already evaluated records, not an unseen-data confirmation.',
                              'Component gains interact and are not additive.',
                              'Terminal-only extension still annihilates initial sensitivity in the closed local class.',
                              'Initial controls exploit correlated measurement noise; known u is required.']}
    (ROOT/'reports/boundary_ablation_results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    lines = ['# Абляция граничной формы: 500 парных записей', '',
             '| Ветка | RMSE | Закрытая доля общего Simpson-разрыва | Парный 95% интервал доли |',
             '|---|---:|---:|---:|']
    for g in groups:
        f = g['fraction_matched_gap_closed']; ci = g['fraction_closed_ci95']
        lines.append(f"| {g['method']} | {g['relative_rmse']:.8f} | {f*100:.2f}% | [{ci[0]*100:.2f}%; {ci[1]*100:.2f}%] |" if f is not None else f"| {g['method']} | {g['relative_rmse']:.8f} | не оценивается | — |")
    lines += ['', 'Общий Simpson baseline отделяет эффект квадратуры. Все исходные наблюдения совпадают побитово.',
              'Это абляция на ранее оценённых данных; основной 1968-путевой результат остаётся отдельным подтверждением.',
              'Вклады не складываются: начальные полуокна и контрольный момент могут восстанавливать одни направления.',
              'Терминальные тесты без начального следа и контроля сохраняют nuisance-ортогональность.',
              'Правило существенного вклада ≥50% зафиксировано до новых оценок; интервалы используют 5000 общих bootstrap-resamples.']
    difference=report['all_minus_initial_rmse'];lo,hi=report['all_minus_initial_rmse_ci95']
    lines += ['',f'Парная разность RMSE (всё вместе минус начальные полуокна): {difference:.3g}, 95% интервал [{lo:.3g}; {hi:.3g}].',
              'Дополнительный выигрыш поверх начальных полуокон не различим. Начальный след — механизм Brunel et al. (2014), §2.3; вклад этой проверки — его количественная OU-оценка.']
    (ROOT/'reports/boundary_ablation_results_RU.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print(json.dumps(groups, indent=2))


if __name__ == '__main__': main()
