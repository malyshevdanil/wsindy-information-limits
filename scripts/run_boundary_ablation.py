"""Frozen endpoint-component ablation, on exactly paired original records."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[key] = '1'
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_efficiency_followup import generate, input_digest
from src.noise import ou_noise
from src.weak_efficiency import nested_bumps, vdp_moments, vdp_pilot_trajectory
from src.universal_boundary import universal_tests, fit


def write(path, data): path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')


def run(stage, out):
    if not str(ROOT).startswith('/kaggle/working/'): raise RuntimeError('Research execution must run on Kaggle')
    p = json.loads((ROOT/'configs/boundary_ablation.json').read_text())
    shard = int(stage.rsplit('_', 1)[-1]); assert 0 <= shard < p['shards']
    out.mkdir(parents=True, exist_ok=True); started = time.perf_counter()
    refpath = ROOT/'configs'/p['reference_file']
    refs = {r['seed']: r for r in json.loads(refpath.read_text())['rows']}
    write(out/'manifest.json', {'stage': stage, 'protocol': p, 'source_sha256': source_hash(),
                               'reference_sha256': hashlib.sha256(refpath.read_bytes()).hexdigest()})
    t, x, base = generate(p); sigma = p['level']*base; truth = np.array([1., -1., 1., -1.])
    phi0, der0, _ = nested_bumps(t, p['pilot_windows'])
    bases = {}; gates = []
    kwargs = dict(iterations=p['iterations'], tolerance=p['tolerance'], relaxation=p['covariance_relative_relaxation'])
    for method, spec in p['branches'].items():
        phi, der, meta = universal_tests(t, p['endpoint_half_window_full_widths'], p['basis_relative_cutoff'], spec['endpoints'])
        bases[method] = (phi, der)
        q, diag = fit(t, x, x, truth, phi, der, np.array(p['initial']), sigma, p['tau'],
                      controls=spec['controls'], mean_moments=True, **kwargs)
        defect = float(np.linalg.norm(q-truth)/2)
        gates.append({'method': method, 'relative_bias': defect, **meta, **diag})
        write(out/'quadrature_gates.json', gates)
        if defect > p['clean_mean_quadrature_bias_max']: raise RuntimeError('Frozen quadrature gate failed')
    rows = []
    for seed in range(p['first_seed_index'], p['last_seed_index']+1):
        if (seed-p['first_seed_index']) % p['shards'] != shard: continue
        y = x+p['level']*ou_noise(len(t), p['h'], p['tau'], base, np.random.default_rng(p['seed_base']+seed))
        digest = input_digest(t, y); assert digest == refs[seed]['input_sha256']
        A, b = vdp_moments(t, y, phi0, der0, sigma); pilot = np.linalg.lstsq(A, b, rcond=None)[0]
        model, constrained = vdp_pilot_trajectory(t, pilot, p['initial'], p['pilot_box'])
        for method, spec in p['branches'].items():
            timer = time.perf_counter(); phi, der = bases[method]
            q, diag = fit(t, y, model, constrained, phi, der, np.array(p['initial']), sigma, p['tau'], controls=spec['controls'], **kwargs)
            if not diag['converged']: raise RuntimeError('Ablation iteration failed')
            rows.append({'seed': seed, 'method': method, 'coefficients': q.tolist(),
                         'relative_error': float(np.linalg.norm(q-truth)/2), 'input_sha256': digest,
                         'fit_seconds': time.perf_counter()-timer, **diag})
        for method in p['archived_branches']:
            rows.append({**refs[seed][method], 'seed': seed, 'method': method, 'input_sha256': digest, 'archived_reference': True})
        write(out/'rows.json', rows)
        if len(rows) % 70 == 0: print(f'Ablation shard {shard}: {len(rows)//7}/125', flush=True)
    write(out/'completion.json', {'status': 'complete', 'paths': len(rows)//7, 'elapsed_seconds': time.perf_counter()-started})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run(args.stage, args.output)
