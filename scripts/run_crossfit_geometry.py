"""Kaggle-only instrument cache, with exact covariance accuracy checks."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_heterogeneous_records import instruments
from src.heterogeneous_gmm import interpolate_grid
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory,vdp_moments,covariance_blocks,assemble_covariance
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical calculations run on Kaggle')
    p=json.loads((ROOT/'configs/crossfit_records.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);box=np.array(p['pilot_box']);anchor=box.mean(axis=1)
    baseline,_=vdp_pilot_trajectory(t,p['truth'],p['initial_center'],p['pilot_box']);base=baseline.std(axis=0)
    phi,der,meta=nested_bumps(t,p['windows'],p['basis_relative_cutoff']);xs=p['instrument_grid']['x'];ys=p['instrument_grid']['v'];step=p['cache_difference_step'];arrays={};checks=[]
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash()))
    for level in p['levels']:
        sigma=np.sqrt(2)*level*base
        L,_=instruments(t,anchor,xs,ys,phi,der,sigma,p['tau'],p['covariance_relative_relaxation']);cache=[L]
        for j in range(4):
            offset=np.eye(4)[j]*step
            plus,_=instruments(t,anchor+offset,xs,ys,phi,der,sigma,p['tau'],p['covariance_relative_relaxation'])
            minus,_=instruments(t,anchor-offset,xs,ys,phi,der,sigma,p['tau'],p['covariance_relative_relaxation'])
            cache.append((plus-minus)/(2*step));print(f'Crossfit cache {level}: derivative {j+1}/4',flush=True)
        cache=np.array(cache);arrays['cache_'+str(level)]=cache
        for delta in p['cache_validation_offsets']:
            q=anchor+delta
            for u in p['cache_validation_initials']:
                i=xs.index(u[0]);j=ys.index(u[1]);approx=cache[0,i,j]+np.einsum('p,pmq->mq',delta,cache[1:,i,j])
                model,_=vdp_pilot_trajectory(t,q,u,p['pilot_box']);D,b=vdp_moments(t,model,phi,der,np.zeros(2))
                blocks,H=covariance_blocks(t,model,phi,der,sigma,p['tau']);C=assemble_covariance(q,blocks,H,sigma)
                C+=p['covariance_relative_relaxation']*np.trace(C)/len(C)*np.eye(len(C))
                exact=cho_solve(cho_factor(C,lower=True,check_finite=False),D,check_finite=False).T
                normalization=approx@D;inverse=np.linalg.inv(normalization)
                Va=inverse@(approx@C@approx.T)@inverse.T;Ve=np.linalg.inv(exact@D)
                loss=float(np.trace(Va)/np.trace(Ve)-1)
                checks.append(dict(level=level,offset=delta,initial=u,variance_loss=loss,weak_defect=float(np.linalg.norm(b-D@q))))
        print(f'Crossfit cache {level}: validation complete',flush=True)
    maximum=max(c['variance_loss'] for c in checks)
    write(out/'geometry.json',dict(anchor=anchor.tolist(),anchor_uses_training_records=False,sigma_base=base.tolist(),
        family=meta,checks=checks,maximum_variance_loss=maximum,threshold=p['maximum_cache_variance_loss'],
        cache_validated=maximum<=p['maximum_cache_variance_loss']))
    if maximum>p['maximum_cache_variance_loss']:raise RuntimeError('Prespecified instrument-cache accuracy check failed')
    np.savez_compressed(out/'crossfit_instruments.npz',anchor=anchor,base=base,**arrays)
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
