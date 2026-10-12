"""Matched 5% path-likelihood fits using analytic sensitivities, Kaggle only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_information_followup import pilot_rows
from scripts.run_efficiency_followup import generate,input_digest
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash
from src.noise import ou_noise
from src.path_scoring import vdp_sensitivities,ou_whiten


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical runs execute on Kaggle')
    p=json.loads((ROOT/'configs/information_followup.json').read_text())['noise_sensitivity'];c=json.loads((ROOT/'configs/moderate_likelihood.json').read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();t,x,base=generate(p);pilots=pilot_rows();truth=np.array([1,-1,1,-1]);rows=[]
    write(out/'manifest.json',dict(stage=stage,protocol=c,observation_protocol=p,source_sha256=source_hash()))
    for seed in range(c['paths']):
        y=x+c['level']*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed));digest=input_digest(t,y)
        for unknown in (False,True):
            initial=y[0] if unknown else np.array(p['initial']);q0=np.array(pilots[seed,c['level']]['coefficients'])
            z0=np.r_[q0,initial] if unknown else q0.copy();lo=np.array(c['coefficient_lower']+c['initial_lower'] if unknown else c['coefficient_lower']);hi=np.array(c['coefficient_upper']+c['initial_upper'] if unknown else c['coefficient_upper'])
            z0=np.clip(z0,lo+1e-6,hi-1e-6);cache={}
            def evaluate(z):
                if 'z' not in cache or not np.array_equal(z,cache['z']):
                    model,U=vdp_sensitivities(t,z[:4],z[4:] if unknown else initial)
                    cache.update(z=z.copy(),residual=ou_whiten(model-y,p['h'],p['tau'],c['level']*base).ravel(),
                        jacobian=ou_whiten(U if unknown else U[:,:,:4],p['h'],p['tau'],c['level']*base).reshape(-1,len(z)))
                return cache
            fit=least_squares(lambda z:evaluate(z)['residual'],z0,jac=lambda z:evaluate(z)['jacobian'],bounds=(lo,hi),
                xtol=c['tolerance'],ftol=c['tolerance'],gtol=c['tolerance'],max_nfev=c['max_nfev'],x_scale='jac')
            if not fit.success:raise RuntimeError('Likelihood optimizer failed; no fits discarded')
            rows.append(dict(seed=seed,method='unknown_initial' if unknown else 'known_initial',coefficients=fit.x[:4].tolist(),
                relative_error=float(np.linalg.norm(fit.x[:4]-truth)/2),input_sha256=digest,status=int(fit.status),nfev=int(fit.nfev),
                optimality=float(fit.optimality),bound_active=fit.active_mask.tolist(),initial_estimate=fit.x[4:].tolist() if unknown else initial.tolist()))
        if seed%8==0:write(out/'rows.json',rows);print(f'Moderate likelihood {seed+1}/{c["paths"]}',flush=True)
    write(out/'rows.json',rows);write(out/'completion.json',dict(status='complete',rows=len(rows),elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);p=a.parse_args();run(p.stage,p.output)
