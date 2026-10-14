"""Held-out endpoint tests on paired archived observations, executed on Kaggle."""
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
from src.boundary_weak import fit as adjoint_fit
from src.universal_boundary import universal_tests, fit as universal_fit


def write(path, data): path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')


def run(stage, out):
    if not str(ROOT).startswith('/kaggle/working/'): raise RuntimeError('Research execution must be on Kaggle')
    p = json.loads((ROOT/'configs/boundary_confirmation.json').read_text())
    shard = int(stage.rsplit('_', 1)[-1]); assert 0 <= shard < p['shards']
    out.mkdir(exist_ok=True, parents=True); start = time.perf_counter()
    refpath = ROOT/'configs'/p['references_file']
    refs = {r['seed']: r for r in json.loads(refpath.read_text())['rows']}
    write(out/'manifest.json', {'stage': stage, 'protocol': p, 'source_sha256': source_hash(),
          'reference_file_sha256': hashlib.sha256(refpath.read_bytes()).hexdigest()})
    t, x, base = generate(p); sigma = p['level']*base; truth = np.array([1., -1., 1., -1.])
    phi0, der0, _ = nested_bumps(t, p['pilot_windows'])
    basis_start = time.perf_counter()
    phi, der, meta = universal_tests(t, p['endpoint_half_window_full_widths'], p['basis_relative_cutoff'])
    write(out/'basis.json', {**meta, 'seconds': time.perf_counter()-basis_start})
    kwargs = dict(iterations=p['iterations'], tolerance=p['tolerance'], relaxation=p['covariance_relative_relaxation'])
    mean, diag = universal_fit(t, x, x, truth, phi, der, np.array(p['initial']), sigma, p['tau'], mean_moments=True, **kwargs)
    defect = float(np.linalg.norm(mean-truth)/2)
    write(out/'quadrature_gate.json', {'relative_bias': defect, 'coefficients': mean.tolist(), **diag})
    if defect > p['clean_mean_quadrature_bias_max']:
        write(out/'completion.json', {'status': 'numerical_gate_failed'})
        raise RuntimeError('Unchanged clean mean quadrature gate failed')
    low, high = np.array(p['pilot_box']).T; rows = []
    for seed in range(p['first_seed_index'], p['last_seed_index']+1):
        if (seed-p['first_seed_index']) % p['shards'] != shard: continue
        y = x+p['level']*ou_noise(len(t), p['h'], p['tau'], base, np.random.default_rng(p['seed_base']+seed))
        digest = input_digest(t, y); assert digest == refs[seed]['input_sha256']
        timer = time.perf_counter(); A, b = vdp_moments(t, y, phi0, der0, sigma)
        pilot = np.linalg.lstsq(A, b, rcond=None)[0]; clipped = np.clip(pilot, low, high)
        pilot_seconds = time.perf_counter()-timer
        timer = time.perf_counter(); q = clipped.copy()
        for _ in range(p['model_adjoint_updates']):
            new, _, diag = adjoint_fit(t, y, q, p['initial'], sigma, p['tau'])
            q = np.clip(new, low, high)
        rows.append({'seed': seed, 'method': 'model_adjoint', 'coefficients': new.tolist(),
            'relative_error': float(np.linalg.norm(new-truth)/2), 'input_sha256': digest,
            'seconds': time.perf_counter()-timer, 'pilot_seconds': pilot_seconds, **diag})
        timer = time.perf_counter(); model, constrained = vdp_pilot_trajectory(t, pilot, p['initial'], p['pilot_box'])
        q, diag = universal_fit(t, y, model, constrained, phi, der, np.array(p['initial']), sigma, p['tau'], **kwargs)
        if not diag['converged']: raise RuntimeError('Universal covariance iteration failed to converge')
        rows.append({'seed': seed, 'method': 'universal_boundary', 'coefficients': q.tolist(),
            'relative_error': float(np.linalg.norm(q-truth)/2), 'input_sha256': digest,
            'seconds': time.perf_counter()-timer, 'pilot_seconds': pilot_seconds, **diag})
        for method in ('boundary_free_480', 'known_initial_scoring'):
            r = refs[seed][method]
            rows.append({'seed': seed, 'method': method, **r, 'input_sha256': digest, 'archived_reference': True})
        write(out/'rows.json', rows)
        if len(rows) % 40 == 0: print(f'Boundary confirmation shard {shard}: {len(rows)//4}/492', flush=True)
    write(out/'completion.json', {'status': 'complete', 'paths': len(rows)//4, 'elapsed_seconds': time.perf_counter()-start})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    run(args.stage, args.output)
