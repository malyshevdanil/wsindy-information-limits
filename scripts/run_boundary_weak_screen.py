"""Direct preregistered known-initial-state screen; research execution on Kaggle."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[key] = '1'
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from src.boundary_weak import fit
from src.weak_efficiency import nested_bumps, vdp_moments
from src.systems import trajectory
from src.noise import ou_noise
from scripts.run_efficiency_followup import input_digest
from scripts.run_experiments import source_hash


def write(path, data): path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')


def run(stage, out):
    if not str(ROOT).startswith('/kaggle/working/'):
        raise RuntimeError('Research calculations must run on Kaggle')
    assert stage in ('boundary_screen', 'boundary_screen_v2', 'boundary_screen_v3', 'boundary_screen_v4')
    cfg_path = ROOT/('configs/'+stage.replace('boundary_screen', 'boundary_weak_screen')+'.json')
    p = json.loads(cfg_path.read_text(encoding='utf-8'))
    out.mkdir(parents=True, exist_ok=True); start = time.perf_counter()
    write(out/'manifest.json', {'stage': stage, 'protocol': p, 'source_sha256': source_hash(),
                                'config_sha256': hashlib.sha256(cfg_path.read_bytes()).hexdigest()})
    t = np.linspace(0, p['T'], round(p['T']/p['h'])+1)
    x = trajectory('vdp', t)
    base = np.std(trajectory('vdp', np.linspace(0, 20, 2001)), axis=0)
    sigma = p['level']*base; truth = np.array([1., -1., 1., -1.])
    phi, derivative, _ = nested_bumps(t, p['pilot_windows'])
    low, high = np.array(p['pilot_box']).T
    # E[Hermite-corrected features(Y)] = features(x), not corrected features(x).
    # Keep the real noise covariance for GLS; disable correction only for the
    # deterministic mean-moment quadrature check.
    clean, _, diagnostic = fit(t, x, truth, p['initial'], sigma, p['tau'],
                               moment_sigma=np.zeros(2))
    defect = float(np.linalg.norm(clean-truth)/2)
    write(out/'quadrature_gate.json', {'clean_coefficients': clean.tolist(), 'relative_bias': defect,
                                      'threshold': p['clean_relative_quadrature_bias_max'], **diagnostic})
    if defect > p['clean_relative_quadrature_bias_max']:
        write(out/'completion.json', {'status': 'numerical_gate_failed', 'relative_bias': defect})
        raise RuntimeError('Prespecified clean quadrature gate failed; no Monte Carlo')
    rows = []
    for seed in range(p['screen_seed_indices'][0], p['screen_seed_indices'][1]+1):
        # Match the archived power protocol bit for bit: scale *after* OU
        # generation. Scaling the innovation standard deviation first changes
        # floating-point rounding and therefore the paired input fingerprint.
        y = x+p['level']*ou_noise(len(t), p['h'], p['tau'], base, np.random.default_rng(p['seed_base']+seed))
        A, b = vdp_moments(t, y, phi, derivative, sigma)
        q = np.clip(np.linalg.lstsq(A, b, rcond=None)[0], low, high)
        timer = time.perf_counter(); history = []
        for _ in range(p['updates']):
            new, V, diagnostic = fit(t, y, q, p['initial'], sigma, p['tau'])
            history.append(float(np.linalg.norm(new-q)/2)); q = np.clip(new, low, high)
        rows.append({'seed': seed, 'method': 'feasible_boundary_weak', 'coefficients': new.tolist(),
                     'relative_error': float(np.linalg.norm(new-truth)/2), 'input_sha256': input_digest(t, y),
                     'seconds': time.perf_counter()-timer, 'update_changes': history,
                     'clipped': bool(np.any(q != new)), 'plugin_covariance': V.tolist(), **diagnostic})
        if p['oracle_diagnostic']:
            oracle, _, _ = fit(t, y, truth, p['initial'], sigma, p['tau'])
            rows.append({'seed': seed, 'method': 'oracle_boundary_weak', 'coefficients': oracle.tolist(),
                         'relative_error': float(np.linalg.norm(oracle-truth)/2), 'input_sha256': input_digest(t, y)})
        write(out/'rows.json', rows)
        print(f'Direct boundary screen {seed+1}/{p["screen_paths"]}', flush=True)
    sample = [r for r in rows if r['method'] == 'feasible_boundary_weak']
    rmse = float(np.sqrt(np.mean([r['relative_error']**2 for r in sample])))
    write(out/'completion.json', {'status': 'complete', 'paths': len(sample), 'relative_rmse': rmse,
                                  'provisional_screen_pass': rmse <= p['success_relative_rmse_max'],
                                  'threshold': p['success_relative_rmse_max'],
                                  'elapsed_seconds': time.perf_counter()-start})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    run(args.stage, args.output)
