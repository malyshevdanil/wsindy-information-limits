"""Direct Gaussian path-distance filter; all numerical research on Kaggle."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[key] = '1'
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.stats import norm
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from src.systems import trajectory
from scripts.run_experiments import source_hash


def run(stage, out):
    if not str(ROOT).startswith('/kaggle/working/'): raise RuntimeError('Research calculations must run on Kaggle')
    assert stage == 'support_distance'
    cp = ROOT/'configs/support_detection_screen.json'; p = json.loads(cp.read_text())
    out.mkdir(parents=True, exist_ok=True); start = time.perf_counter()
    def write(name, data): (out/name).write_text(json.dumps(data, indent=2, allow_nan=False), encoding='utf-8')
    write('manifest.json', {'stage': stage, 'protocol': p, 'source_sha256': source_hash(),
                            'config_sha256': hashlib.sha256(cp.read_bytes()).hexdigest()})
    sigma = p['noise_level']*np.std(trajectory('vdp', np.linspace(0, 20, 2001)), axis=0)
    def rhs(delta, state):
        z, v = state; return np.array([v, -z+v-delta*z*z*v])
    def solve(delta):
        sol = solve_ivp(lambda t, y: rhs(delta, y), (0., p['T']), p['initial'],
                        method='DOP853', dense_output=True, rtol=p['ode_rtol'], atol=p['ode_atol'])
        if not sol.success: raise RuntimeError('Direct ODE integration failed')
        return sol.sol
    null = solve(0.); rows = []
    for delta in p['delta_grid']:
        alternative = solve(delta)
        terms = []; errors = []
        for j in range(2):
            def integrand(t):
                xn, xa = null(t), alternative(t)
                gp = rhs(delta, xa)-rhs(0., xn); g = xa-xn
                return (gp[j]+g[j]/p['tau'])**2
            value, error = quad(integrand, 0., p['T'], epsabs=1e-15,
                                 epsrel=p['distance_quadrature_relative_tolerance'], limit=1000)
            terms.append(p['tau']*value/(2*sigma[j]**2)); errors.append(error)
        distance = float(np.sqrt(sum(terms)))
        row = {'delta': delta, 'path_distance': distance, 'equal_prior_error': float(norm.cdf(-distance/2)),
               'quadrature_absolute_error_integrals': errors,
               'alternative_terminal_state': alternative(p['T']).tolist()}
        rows.append(row); write('rows.json', rows)
        print(f'Direct support-distance delta={delta:g}, distance={distance:.8g}', flush=True)
    write('completion.json', {'status': 'complete', 'rows': len(rows), 'sigma': sigma.tolist(),
                               'null_terminal_state': null(p['T']).tolist(),
                               'target_distance': float(2*norm.ppf(1-p['target_equal_prior_error'])),
                               'elapsed_seconds': time.perf_counter()-start,
                               'scope': p['scope']})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    run(args.stage, args.output)
