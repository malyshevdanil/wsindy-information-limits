"""Added-term local information and complete author WSINDy inclusion study."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[key] = '1'
import argparse, hashlib, io, json, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from src.support_information import augmented_sensitivities, references
from src.systems import trajectory
from src.noise import ou_noise
from scripts.run_experiments import source_hash
from scripts.run_full_adaptive import fetch, sha, run_octave
from scripts.run_efficiency_followup import input_digest


def write(path, data): path.write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')


def run(stage, out):
    if not str(ROOT).startswith('/kaggle/working/'): raise RuntimeError('Research calculations must run on Kaggle')
    assert stage == 'support_extra'
    p = json.loads((ROOT/'configs/support_extra_term.json').read_text())
    out.mkdir(parents=True, exist_ok=True); started = time.perf_counter()
    write(out/'manifest.json', {'stage': stage, 'protocol': p, 'source_sha256': source_hash()})
    base = np.std(trajectory('vdp', np.linspace(0, 20, 2001)), axis=0)
    constant = float(2*norm.ppf(1-p['target_equal_prior_error']))
    checks = []; info = {}; grids = {}; paths = {}
    for h in [p['h'], p['information_refinement_h']]:
        t = np.linspace(0, p['T'], round(p['T']/h)+1)
        _, S = augmented_sensitivities(t, p['initial'])
        for level in p['levels']:
            for cutoff in p['information_svd_cutoffs']:
                ref = references(t, S, level*base, p['tau'], cutoff)
                for label, r in ref.items():
                    checks.append({'h': h, 'level': level, 'cutoff': cutoff, 'label': label,
                                   'delta_min_local': constant*r['standard_deviation'], **r})
                    if h == p['information_refinement_h'] and cutoff == p['primary_svd_cutoff']:
                        info[level, label] = r
    write(out/'information.json', {'constant': constant, 'rows': checks})
    for level in p['levels']:
        for label in [r['label'] for r in checks if r['h'] == p['h'] and r['level'] == level and r['cutoff'] == p['primary_svd_cutoff']]:
            selected = [r for r in checks if r['level'] == level and r['label'] == label]
            fine = info[level, label]['information']
            if max(abs(r['information']/fine-1) for r in selected) > max(p['relative_information_refinement_tolerance'], p['relative_svd_sensitivity_tolerance']):
                write(out/'completion.json', {'status': 'information_gate_failed', 'level': level, 'label': label})
                raise RuntimeError('Information refinement/rank gate failed before Monte Carlo')
        delta_ref = constant*info[level, 'full_library_boundary_free']['standard_deviation']
        grids[level] = [m*delta_ref for m in p['delta_multipliers_of_weak_reference']]
    t = np.linspace(0, p['T'], round(p['T']/p['h'])+1)
    for level in p['levels']:
        for index, delta in enumerate(grids[level]):
            def rhs(_, state):
                z, v = state; return [v, -z+v-z*z*v+delta*z**3]
            sol = solve_ivp(rhs, (0., p['T']), p['initial'], t_eval=t, method='DOP853', rtol=1e-10, atol=1e-12)
            if not sol.success or np.max(np.abs(sol.y)) > p['maximum_clean_state_abs']:
                write(out/'completion.json', {'status': 'boundedness_gate_failed', 'delta': delta})
                raise RuntimeError('Prespecified boundedness gate failed')
            paths[level, index] = sol.y.T
    write(out/'delta_grids.json', {str(level): grid for level, grid in grids.items()})
    pin = json.loads((ROOT/'configs/author_pins.json').read_text())['WSINDy_ODE']
    data = fetch(f'https://codeload.github.com/{pin["repo"]}/zip/{pin["commit"]}')
    assert sha(data) == pin['zip_sha256']
    dest = ROOT/'author_sources/WSINDy_ODE'; dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for member in z.infolist():
            path = Path(member.filename)
            if len(path.parts) == 2 and path.suffix in ('.m', '.md'): (dest/path.name).write_bytes(z.read(member))
    write(out/'author_provenance.json', {**pin, 'files_sha256': {f.name: sha(f.read_bytes()) for f in dest.iterdir() if f.is_file()}})
    rows = []
    for li, level in enumerate(p['levels']):
        for seed in range(p['paths_per_delta_per_level']):
            noise = level*ou_noise(len(t), p['h'], p['tau'], base, np.random.default_rng(p['seed_base']+li*p['level_seed_stride']+seed))
            for di, delta in enumerate(grids[level]):
                y = paths[level, di]+noise; fitdir = out/f'fit_l{li}_d{di}_s{seed:03d}'; fitdir.mkdir(exist_ok=True)
                timer = time.perf_counter(); result = run_octave('wsindy', t, y, fitdir, {'WSINDy_ODE': dest})
                tags = np.asarray(result['tags']); match = np.flatnonzero(np.all(tags == [3, 0], axis=1))
                assert len(match) == 1
                coefficient = float(np.asarray(result['coefficients'])[match[0], 1])
                rows.append({'level': level, 'seed': seed, 'delta_index': di, 'delta': delta,
                    'input_sha256': input_digest(t, y), 'candidate_coefficient': coefficient,
                    'included': coefficient != 0., 'fit_seconds': time.perf_counter()-timer,
                    'method': 'complete_author_wsindy', 'status': 'ok'})
                write(out/'rows.json', rows)
            if seed % 10 == 0: print(f'Added-term author WSINDy level={level}: {seed+1}/100', flush=True)
    write(out/'completion.json', {'status': 'complete', 'rows': len(rows), 'elapsed_seconds': time.perf_counter()-started})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    run(args.stage, args.output)
