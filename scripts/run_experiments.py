"""Run explicitly on Kaggle; checkpoints are not manuscript-ready claims."""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_name, '1')

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time
import traceback
import numpy as np
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.noise import stationary_ar1, ou_noise, calibrated_ar1
from src.weak_features import weak_matrices
from src.covariance import residual_covariance, coefficient_surrogate
from src.systems import trajectory, system, forecast_error
from src.selection import select_threshold, stls
from src.diagnostic import choose_grid


def rng_for(*parts):
    key = json.dumps(parts, sort_keys=True).encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(key).digest()[:8], 'little'))


def source_hash():
    h = hashlib.sha256()
    for folder in ('src', 'scripts', 'tests'):
        for path in sorted((ROOT / folder).glob('*.py')):
            h.update(str(path.relative_to(ROOT)).replace('\\', '/').encode())
            h.update(path.read_bytes())
    return h.hexdigest()


def serializable(obj):
    if isinstance(obj, np.ndarray):
        return serializable(obj.tolist())
    if isinstance(obj, dict):
        return {str(k): serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [serializable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def make_noise(regime, n, h, tau, sigma, seed, cell, level=0, rho=.8):
    rng = rng_for(seed, cell, regime, level if regime == 'fixed_index' else 'nested')
    if regime == 'physical_ou':
        return ou_noise(n, h, tau, sigma, rng)
    if regime == 'iid':
        return rng.normal(size=(n, len(sigma))) * sigma
    return stationary_ar1(n, sigma, rho, rng)


def sum_coefficient_variance(G, cov, coef, mask):
    total = 0.
    d = coef.shape[1]
    for equation in range(d):
        v = coefficient_surrogate(G, cov, equation, d, mask[:, equation])
        if v is None:
            return float('nan')
        total += np.trace(v)
    return float(total)


def pilot_diagnostic(cfg, t, y, powers, coef, centers, width, sigma, tau,
                     correction, strides, validation_paths):
    """Pilot plug-in x is smoothed training data; never test truth."""
    h = t[1] - t[0]
    size = max(5, int(round(cfg['pilot_smoothing_width'] / h)) | 1)
    size = min(size, (len(t) - 1) | 1)
    smooth = savgol_filter(y, size, 3, axis=0)
    mask = np.abs(coef) > cfg['support_metric_threshold']
    Gref, bref = weak_matrices(t, smooth, powers, centers, width)
    covref = residual_covariance(t, smooth, powers, coef, centers, width,
                                 sigma, tau, correction_variance=correction, first_order=True)
    reference = sum_coefficient_variance(Gref, covref, coef, mask)
    # A reference at the finest recorded mesh is a proxy, never labelled P_infinity.
    records = []
    for stride in strides:
        tc, xs, ys = t[::stride], smooth[::stride], y[::stride]
        Gclean, bclean = weak_matrices(tc, xs, powers, centers, width)
        cov = residual_covariance(tc, xs, powers, coef, centers, width,
                                  sigma, tau, correction_variance=correction, first_order=True)
        variance = sum_coefficient_variance(Gclean, cov, coef, mask)
        G, b = weak_matrices(tc, ys, powers, centers, width, correction)
        fit = stls(G, b, pilot_diagnostic.threshold)
        losses = []
        for yv in validation_paths:
            Gv, bv = weak_matrices(tc, yv[::stride], powers, centers, width, correction)
            losses.append(np.linalg.norm(bv - Gv @ fit) ** 2 / max(np.linalg.norm(bv) ** 2, 1e-20))
        change = np.linalg.norm(Gclean - Gref) / max(np.linalg.norm(Gref), 1e-20)
        change = max(change, np.linalg.norm(bclean - bref) / max(np.linalg.norm(bref), 1e-20))
        well_conditioned = all(mask[:, i].any() and
                               np.linalg.cond(Gclean[:, mask[:, i]]) < 1e8
                               for i in range(coef.shape[1]))
        records.append({'h': float(h * stride), 'variance': variance,
                        'quadrature_difference': float(change),
                        'validation_loss': float(np.mean(losses)),
                        'support_stable': bool(np.array_equal(np.abs(fit) > cfg['support_metric_threshold'], mask)),
                        'conditioning_ok': well_conditioned})
    return {'selection': choose_grid(records, reference, **cfg['diagnostic']),
            'records': records, 'reference_variance': reference,
            'state': 'smoothed_train_trajectory', 'support': mask}


def run(cfg, output):
    for a, b in (('train_seeds', 'validation_seeds'), ('train_seeds', 'test_seeds'),
                 ('validation_seeds', 'test_seeds')):
        if set(cfg[a]) & set(cfg[b]):
            raise ValueError('Train/validation/test seeds must be disjoint')
    output.mkdir(parents=True, exist_ok=True)
    signature = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    manifest_path = output / 'manifest.json'
    fingerprint = source_hash()
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old['config_sha256'] != signature or old['source_sha256'] != fingerprint:
            raise ValueError('Output contains a different configuration/source; use a new output directory')
    manifest = {'config': cfg, 'config_sha256': signature, 'source_sha256': fingerprint,
                'python': platform.python_version(), 'platform': platform.platform(),
                'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'matplotlib')},
                'inference_status': cfg['stage'], 'method': 'custom_weak_STLS_not_MSTLS'}
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    rows_path = output / 'runs.jsonl'
    completed = set()
    if rows_path.exists():
        for line in rows_path.read_text().splitlines():
            row = json.loads(line)
            if row['status'] == 'ok':
                completed.add(row['key'])
    cells = list(__import__('itertools').product(cfg['systems'], cfg['regimes'],
                                                cfg['tau_over_width'], cfg['noise_levels']))
    for cell_index, (name, regime, tau_ratio, level) in enumerate(cells):
        info, powers, truth = system(name)
        width, T = cfg['window_width'][name], cfg['T'][name]
        finest_h = cfg['base_step'] / 2 ** (cfg['refinement_levels'] - 1)
        intervals = int(round(T / finest_h))
        if not np.isclose(intervals * finest_h, T):
            raise ValueError('T must be a multiple of the finest step')
        t = np.linspace(0, T, intervals + 1)
        clean = trajectory(name, t)
        sigma = level * np.std(clean, axis=0)
        tau = tau_ratio * width
        centers = np.linspace(width / 2 + finest_h, T - width / 2 - finest_h, cfg['window_count'])
        strides = [2 ** k for k in reversed(range(cfg['refinement_levels']))]
        cell = (name, regime, tau_ratio, level)
        paths = {s: clean + make_noise(regime, len(t), finest_h, tau, sigma, s, cell,
                                       rho=cfg['fixed_index_rho'])
                 for s in cfg['train_seeds'] + cfg['validation_seeds'] + cfg['test_seeds']}
        cal_n = max(200, int(round(cfg['calibration_duration'] / finest_h)) + 1)
        calibration = make_noise(regime, cal_n, finest_h, tau, sigma,
                                 cfg['calibration_seed'], cell, rho=cfg['fixed_index_rho'])
        calibration_fit = calibrated_ar1(calibration, finest_h)
        if regime == 'iid':
            calibration_fit['variance'] = np.mean(calibration ** 2, axis=0)
        for method in cfg['methods']:
            correction = (np.zeros(info['d']) if method == 'raw' else sigma ** 2
                          if method == 'known_debias' else calibration_fit['variance'])
            training = [weak_matrices(t, paths[s], powers, centers, width, correction)
                        for s in cfg['train_seeds']]
            validation_paths = [paths[s] for s in cfg['validation_seeds']]
            validation = [weak_matrices(t, yv, powers, centers, width, correction)
                          for yv in validation_paths]
            Gtrain = np.vstack([a[0] for a in training])
            btrain = np.vstack([a[1] for a in training])
            threshold, pilot_coef, loss = select_threshold(
                Gtrain, btrain, validation, cfg['threshold_candidates'], cfg['complexity_penalty'])
            diagnostic = {'selection': {'status': 'control_regime_no_physical_infill_claim', 'h': None}}
            if regime == 'physical_ou':
                predicted_sigma = np.sqrt(calibration_fit['variance']) if method == 'calibrated_debias' else sigma
                predicted_tau = calibration_fit['tau'] if method == 'calibrated_debias' else tau
                pilot_diagnostic.threshold = threshold
                diagnostic = pilot_diagnostic(cfg, t, paths[cfg['train_seeds'][0]], powers,
                                               pilot_coef, centers, width, predicted_sigma, predicted_tau,
                                               correction, strides, validation_paths)
            diagnostic.update({'cell': cell, 'method': method, 'threshold': threshold,
                               'validation_score': loss, 'calibration': calibration_fit})
            tag = hashlib.sha256(json.dumps([cell, method]).encode()).hexdigest()[:16]
            (output / f'diagnostic_{tag}.json').write_text(
                json.dumps(serializable(diagnostic), indent=2, allow_nan=False), encoding='utf-8')
            for seed in cfg['test_seeds']:
                for mesh_index, stride in enumerate(strides):
                    key = json.dumps([cell, method, seed, stride], separators=(',', ':'))
                    if key in completed:
                        continue
                    row = {'key': key, 'system': name, 'regime': regime, 'method': method,
                           'tau': tau, 'tau_over_width': tau_ratio, 'noise_level': level,
                           'seed': seed, 'h': finest_h * stride, 'T': T, 'K': len(centers),
                           'width': width, 'sigma': sigma, 'threshold': threshold,
                           'selected_h': diagnostic['selection']['h'], 'powers': powers,
                           'noise_pairing': 'new_process_per_h' if regime == 'fixed_index' else 'nested_record',
                           'calibration_seed': cfg['calibration_seed']}
                    try:
                        tc, xc = t[::stride], clean[::stride]
                        if regime == 'fixed_index':
                            yc = xc + make_noise(regime, len(tc), row['h'], tau, sigma,
                                                  seed, cell, mesh_index + 1, cfg['fixed_index_rho'])
                        else:
                            yc = paths[seed][::stride]
                        start = time.perf_counter()
                        G, b = weak_matrices(tc, yc, powers, centers, width, correction)
                        estimate = stls(G, b, threshold)
                        if not np.all(np.isfinite(estimate)):
                            raise FloatingPointError('Nonfinite coefficient estimate')
                        row['fit_seconds'] = time.perf_counter() - start
                        mask, target = np.abs(estimate) > cfg['support_metric_threshold'], np.abs(truth) > 0
                        row.update({'coefficients': estimate, 'true_coefficients': truth,
                                    'exact_support': bool(np.array_equal(mask, target)),
                                    'false_positive': int(np.count_nonzero(mask & ~target)),
                                    'false_negative': int(np.count_nonzero(~mask & target)),
                                    'coefficient_error': float(np.linalg.norm(estimate - truth) / np.linalg.norm(truth)),
                                    'design_condition': float(np.linalg.cond(G))})
                        residual = (b - G @ truth).reshape(-1)
                        noise_tau = tau if regime == 'physical_ou' else (-row['h'] / np.log(cfg['fixed_index_rho']))
                        sc, sw = (np.zeros_like(sigma), sigma) if regime == 'iid' else (sigma, None)
                        start = time.perf_counter()
                        cov = residual_covariance(tc, xc, powers, truth, centers, width,
                                                   sc, noise_tau, sw, correction)
                        cov_first = residual_covariance(tc, xc, powers, truth, centers, width,
                                                         sc, noise_tau, sw, correction, first_order=True)
                        Gclean, bclean = weak_matrices(tc, xc, powers, centers, width)
                        # For imperfect/raw correction exact mean differs from the clean residual.
                        from src.weak_features import heat_transform, library, test_functions, quadrature_weights
                        effective = heat_transform(powers, truth, sigma ** 2 - correction)
                        phi, _ = test_functions(tc, centers, width)
                        mean = (bclean - (phi * quadrature_weights(tc)) @ library(xc, powers) @ effective).reshape(-1)
                        row['oracle_linear_surrogate_variance'] = sum_coefficient_variance(
                            Gclean, cov_first, truth, target)
                        row['covariance_seconds'] = time.perf_counter() - start
                        record_name = hashlib.sha256(key.encode()).hexdigest()[:20] + '.npz'
                        np.savez_compressed(output / record_name, residual=residual, expected_mean=mean,
                                            exact_covariance=cov, first_order_covariance=cov_first,
                                            clean_weak_residual=bclean - Gclean @ truth)
                        row['covariance_file'] = record_name
                        row['forecast_error'], row['forecast_status'] = forecast_error(name, estimate)
                        row['status'] = 'ok'
                    except Exception as error:
                        row.update({'status': 'failed', 'error': str(error), 'traceback': traceback.format_exc()})
                    with rows_path.open('a', encoding='utf-8') as stream:
                        stream.write(json.dumps(serializable(row), allow_nan=False) + '\n')
                    print(json.dumps({'cell': cell_index + 1, 'total_cells': len(cells),
                                      'method': method, 'seed': seed, 'h': row['h'],
                                      'status': row['status']}), flush=True)
        print(f'Completed cell {cell_index + 1}/{len(cells)}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/pilot.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'results/kaggle_pilot')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    configuration = json.loads(args.config.read_text(encoding='utf-8'))
    if args.dry_run:
        number = np.prod([len(configuration[k]) for k in
                          ('systems', 'regimes', 'tau_over_width', 'noise_levels', 'methods', 'test_seeds')])
        print(json.dumps({'test_fits': int(number) * configuration['refinement_levels'],
                          'stage': configuration['stage'], 'run_location': 'Kaggle CPU'}))
    else:
        run(configuration, args.output)
