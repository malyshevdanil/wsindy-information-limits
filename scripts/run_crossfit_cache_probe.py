"""Kaggle-only diagnostic of normalized instrument interpolation, before MC."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory,vdp_moments,covariance_blocks,assemble_covariance
from scripts.run_efficiency_validation import write


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Calculations run on Kaggle')
    p=json.loads((ROOT/'configs/crossfit_records.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);anchor=np.mean(p['pilot_box'],axis=1)
    model,_=vdp_pilot_trajectory(t,anchor,p['initial_center'],p['pilot_box']);base=model.std(axis=0)
    phi,der,_=nested_bumps(t,p['windows'],p['basis_relative_cutoff']);checks=[]
    def exact(q,u,sigma):
        model,_=vdp_pilot_trajectory(t,q,u,p['pilot_box']);D,_=vdp_moments(t,model,phi,der,np.zeros(2))
        blocks,H=covariance_blocks(t,model,phi,der,sigma,p['tau']);C=assemble_covariance(q,blocks,H,sigma)
        C+=p['covariance_relative_relaxation']*np.trace(C)/len(C)*np.eye(len(C))
        L=cho_solve(cho_factor(C,lower=True,check_finite=False),D,check_finite=False).T
        return np.linalg.solve(L@D,L),D,C
    for level in p['levels']:
        sigma=np.sqrt(2)*level*base
        for u in p['cache_validation_initials']:
            M,_,_=exact(anchor,u,sigma);derivatives=[]
            for j in range(4):
                offset=np.eye(4)[j]*p['cache_difference_step']
                plus,_,_=exact(anchor+offset,u,sigma);minus,_,_=exact(anchor-offset,u,sigma)
                derivatives.append((plus-minus)/(2*p['cache_difference_step']))
            for delta in p['cache_validation_offsets']:
                delta=np.array(delta);Me,D,C=exact(anchor+delta,u,sigma)
                approx=M+np.einsum('p,pmq->mq',delta,derivatives);approx=np.linalg.solve(approx@D,approx)
                loss=float(np.trace(approx@C@approx.T)/np.trace(Me@C@Me.T)-1)
                checks.append(dict(level=level,initial=u,offset=delta.tolist(),variance_loss=loss))
            print(f'Normalized instrument probe {level}, {u}: complete',flush=True)
    write(out/'probe.json',dict(checks=checks,maximum_variance_loss=max(c['variance_loss'] for c in checks)))
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
