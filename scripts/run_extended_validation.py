"""Targeted synthetic checks; run on Kaggle, never infer sparse recovery from oracle fits."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.noise import ou_noise
from src.systems import system, trajectory
from src.weak_features import weak_matrices
from src.covariance import residual_covariance


def run(output, residual_repeats=4000, repetition_repeats=400):
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    info, powers, truth = system('vdp')
    t = np.linspace(0, 20, 2001)
    clean = trajectory('vdp', t)
    centers = np.linspace(1.505, 18.495, 30)
    unit_sigma = np.std(clean, axis=0)
    settings = {'stage': 'targeted_exploratory_validation', 'system': 'vdp', 'T': 20,
                'h': .01, 'width': 3., 'tau': .3, 'K': 30,
                'noise_levels': [.05, .2, .5], 'residual_repeats': residual_repeats,
                'residual_seeds': [20000, 20000+residual_repeats-1],
                'independent_repetition_blocks': repetition_repeats,
                'replication_seed_formula': '50000 + 8*block + replicate',
                'repetition_counts': [1, 2, 4, 8], 'repetition_noise_level': .05,
                'repetition_fit': 'oracle true support, stacked noisy corrected weak design',
                'selection': 'none; no sparse-support claim',
                'versions': {p: importlib.metadata.version(p) for p in ('numpy','scipy')},
                'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output/'manifest.json').write_text(json.dumps(settings, indent=2))
    residual_summaries = []
    for level in settings['noise_levels']:
        sigma = level*unit_sigma
        for method in ('raw', 'known_debias'):
            correction = np.zeros(2) if method == 'raw' else sigma**2
            covariance = residual_covariance(t, clean, powers, truth, centers, 3., sigma, .3,
                                             correction_variance=correction)
            first = residual_covariance(t, clean, powers, truth, centers, 3., sigma, .3,
                                        correction_variance=correction, first_order=True)
            # An independent pointwise expectation calculation supplies quadrature mean.
            from src.weak_features import heat_transform, library, test_functions, quadrature_weights
            Gclean, bclean = weak_matrices(t, clean, powers, centers, 3.)
            effective = heat_transform(powers, truth, sigma**2-correction)
            phi, _ = test_functions(t, centers, 3.)
            mean = (bclean-(phi*quadrature_weights(t))@library(clean,powers)@effective).reshape(-1)
            residuals = np.empty((residual_repeats, covariance.shape[0]))
            for k in range(residual_repeats):
                noise = ou_noise(len(t), .01, .3, sigma, np.random.default_rng(20000+k))
                G, b = weak_matrices(t, clean+noise, powers, centers, 3., correction)
                residuals[k] = (b-G@truth).reshape(-1)-mean
            empirical = residuals.T@residuals/residual_repeats
            squared_norm = np.sum(residuals**2, axis=1)
            trace_se = float(np.std(squared_norm,ddof=1)/np.sqrt(residual_repeats))
            tag = f'{method}_{level:g}'
            np.savez_compressed(output/f'residual_{tag}.npz', centered_residuals=residuals,
                                exact_covariance=covariance, first_order_covariance=first,
                                expected_mean=mean)
            residual_summaries.append({'method':method,'noise_level':level,'n':residual_repeats,
                'exact_relative_frobenius_error':float(np.linalg.norm(empirical-covariance)/np.linalg.norm(covariance)),
                'first_relative_frobenius_error':float(np.linalg.norm(empirical-first)/np.linalg.norm(covariance)),
                'higher_order_trace_fraction':float(np.trace(covariance-first)/np.trace(covariance)),
                'predicted_trace':float(np.trace(covariance)), 'empirical_trace':float(np.trace(empirical)),
                'empirical_trace_standard_error':trace_se,
                'trace_normal_95_low':float(np.trace(empirical)-1.96*trace_se),
                'trace_normal_95_high':float(np.trace(empirical)+1.96*trace_se)})
            print(json.dumps(residual_summaries[-1]),flush=True)
    (output/'residual_summary.json').write_text(json.dumps(residual_summaries,indent=2))
    # Increasing the number of independent recordings changes the information budget.
    # Each block supplies eight records, with paired prefixes for M=1,2,4,8.
    sigma = .05*unit_sigma
    target = truth != 0
    repetition_rows = []
    for block in range(repetition_repeats):
        matrices = []
        for replicate in range(8):
            noise = ou_noise(len(t), .01, .3, sigma,
                             np.random.default_rng(50000+8*block+replicate))
            matrices.append(weak_matrices(t, clean+noise, powers, centers, 3., sigma**2))
        for count in (1,2,4,8):
            G=np.vstack([pair[0] for pair in matrices[:count]])
            b=np.vstack([pair[1] for pair in matrices[:count]])
            estimate=np.zeros_like(truth)
            for equation in range(2):
                estimate[target[:,equation],equation]=np.linalg.lstsq(
                    G[:,target[:,equation]],b[:,equation],rcond=None)[0]
            repetition_rows.append({'block':block,'M':count,
                'coefficient_error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth))})
    (output/'repetitions.json').write_text(json.dumps(repetition_rows))
    summaries=[]
    for count in (1,2,4,8):
        errors=[r['coefficient_error'] for r in repetition_rows if r['M']==count]
        summaries.append({'M':count,'blocks':repetition_repeats,
                          'error_quantiles_25_50_75_90':np.quantile(errors,[.25,.5,.75,.9]).tolist()})
    (output/'repetition_summary.json').write_text(json.dumps(summaries,indent=2))
    # Refine deterministic state quadrature beyond the pilot reference.
    reference_rows=[]
    for tau in (.06,.3,1.5):
        predictions=[]
        for h in (.01,.005,.0025):
            tref=np.linspace(0,20,round(20/h)+1)
            xref=trajectory('vdp',tref)
            cov=residual_covariance(tref,xref,powers,truth,centers,3.,sigma,tau,
                                    correction_variance=sigma**2)
            predictions.append((h,cov))
        ref=predictions[-1][1]
        for h,cov in predictions:
            reference_rows.append({'tau':tau,'h':h,'reference_h':.0025,
                'relative_frobenius_gap':float(np.linalg.norm(cov-ref)/np.linalg.norm(ref)),
                'relative_trace_gap':float(abs(np.trace(cov)-np.trace(ref))/np.trace(ref))})
    (output/'reference_convergence.json').write_text(json.dumps(reference_rows,indent=2))
    (output/'completion.json').write_text(json.dumps({'status':'complete',
        'elapsed_seconds':time.perf_counter()-started,
        'residual_vectors':6*residual_repeats,'oracle_support_fits':4*repetition_repeats}))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.output)
