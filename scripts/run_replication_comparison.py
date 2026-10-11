"""Compare three estimating procedures on aligned independent noisy recordings."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.systems import system,trajectory
from src.noise import ou_noise
from src.weak_features import weak_matrices
from scripts.run_experiments import source_hash


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    _,powers,truth=system('vdp')
    target=truth!=0
    t=np.linspace(0,20,2001)
    clean=trajectory('vdp',t)
    sigma=.05*np.std(clean,axis=0)
    centers=np.linspace(1.505,18.495,30)
    counts=(1,2,4,8,16,32)
    methods=('stacked_design','mean_weak_matrices','mean_sensor_record')
    blocks=400
    (output/'manifest.json').write_text(json.dumps({'stage':'exploratory_estimator_comparison',
        'system':'vdp','T':20,'h':.01,'tau':.3,'width':3.,'K':30,'noise_level':.05,
        'blocks':blocks,'M':counts,'methods':methods,'seed_formula':'70000+32*block+replicate',
        'alignment':'same deterministic initial condition and time grid',
        'support':'oracle true support in every fit; not sparse discovery',
        'correction':'known variance; sigma^2/M for the averaged sensor record',
        'source_sha256':source_hash()},indent=2))
    rows=[]
    population_gram=[np.zeros((sum(target[:,i]),sum(target[:,i]))) for i in range(2)]
    population_rhs=[np.zeros(sum(target[:,i])) for i in range(2)]
    mean_G=np.zeros((30,len(powers)))
    mean_b=np.zeros((30,2))
    for block in range(blocks):
        matrices=[]
        sum_record=np.zeros_like(clean)
        for replicate in range(32):
            e=ou_noise(len(t),.01,.3,sigma,np.random.default_rng(70000+32*block+replicate))
            y=clean+e
            G,b=weak_matrices(t,y,powers,centers,3.,sigma**2)
            matrices.append((G,b))
            sum_record+=y
            mean_G+=G;mean_b+=b
            for equation in range(2):
                GS=G[:,target[:,equation]]
                population_gram[equation]+=GS.T@GS
                population_rhs[equation]+=GS.T@b[:,equation]
            M=replicate+1
            if M not in counts:
                continue
            stacked=(np.vstack([p[0] for p in matrices]),np.vstack([p[1] for p in matrices]))
            weakmean=(np.mean([p[0] for p in matrices],axis=0),np.mean([p[1] for p in matrices],axis=0))
            sensor=weak_matrices(t,sum_record/M,powers,centers,3.,sigma**2/M)
            for method,(design,response) in zip(methods,(stacked,weakmean,sensor)):
                estimate=np.zeros_like(truth)
                for equation in range(2):
                    estimate[target[:,equation],equation]=np.linalg.lstsq(
                        design[:,target[:,equation]],response[:,equation],rcond=None)[0]
                rows.append({'block':block,'M':M,'method':method,
                    'coefficient_error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth))})
        if block%50==0:
            print(f'Completed independent block {block+1}/{blocks}',flush=True)
    population=np.zeros_like(truth)
    averaged=np.zeros_like(truth)
    for equation in range(2):
        population[target[:,equation],equation]=np.linalg.solve(population_gram[equation],population_rhs[equation])
        averaged[target[:,equation],equation]=np.linalg.lstsq(mean_G[:,target[:,equation]],mean_b[:,equation],rcond=None)[0]
    summaries=[]
    for method in methods:
        for M in counts:
            errors=[r['coefficient_error'] for r in rows if r['method']==method and r['M']==M]
            summaries.append({'method':method,'M':M,'blocks':blocks,
                              'error_quantiles_25_50_75_90':np.quantile(errors,[.25,.5,.75,.9]).tolist()})
    (output/'rows.json').write_text(json.dumps(rows))
    (output/'summary.json').write_text(json.dumps(summaries,indent=2))
    (output/'population_check.json').write_text(json.dumps({'records':blocks*32,
        'interpretation':'Monte Carlo pooled designs, not an analytical limiting theorem',
        'stacked_design_relative_error':float(np.linalg.norm(population-truth)/np.linalg.norm(truth)),
        'mean_weak_design_relative_error':float(np.linalg.norm(averaged-truth)/np.linalg.norm(truth)),
        'stacked_coefficients':population.tolist(),'mean_weak_coefficients':averaged.tolist()},indent=2))
    (output/'completion.json').write_text(json.dumps({'status':'complete','elapsed_seconds':time.perf_counter()-start,
        'independent_records':blocks*32,'oracle_support_fits':len(rows)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
