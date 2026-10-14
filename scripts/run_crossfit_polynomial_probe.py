"""Kaggle-only validation of a richer, bounded instrument approximation."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.crossfit_cache import powers,features
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory,vdp_moments,covariance_blocks,assemble_covariance
from scripts.run_efficiency_validation import write


def exact(t,q,u,sigma,phi,der,p):
    model,_=vdp_pilot_trajectory(t,q,u,p['pilot_box']);D,_=vdp_moments(t,model,phi,der,np.zeros(2))
    blocks,H=covariance_blocks(t,model,phi,der,sigma,p['tau']);C=assemble_covariance(q,blocks,H,sigma)
    C+=p['covariance_relative_relaxation']*np.trace(C)/len(C)*np.eye(len(C))
    L=cho_solve(cho_factor(C,lower=True,check_finite=False),D,check_finite=False).T
    return np.linalg.solve(L@D,L),D,C


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Calculations run on Kaggle')
    p=json.loads((ROOT/'configs/crossfit_records.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);anchor=np.mean(p['pilot_box'],axis=1);radius=.06
    model,_=vdp_pilot_trajectory(t,anchor,p['initial_center'],p['pilot_box']);base=model.std(axis=0)
    phi,der,_=nested_bumps(t,p['windows'],p['basis_relative_cutoff']);checks=[]
    nodes=anchor+radius*(2*qmc.Sobol(4,scramble=True,seed=3800000).random_base2(9)-1)
    exponents=powers(5);design=features(nodes,anchor,radius,exponents);inverse=np.linalg.pinv(design)
    # Include the most sensitive of the predeclared initial-state probes.
    u=p['cache_validation_initials'][-1]
    for level in p['levels']:
        sigma=np.sqrt(2)*level*base;maps=[]
        for i,q in enumerate(nodes):
            M,_,_=exact(t,q,u,sigma,phi,der,p);maps.append(M)
            if (i+1)%64==0:print(f'Polynomial instrument probe {level}: {i+1}/{len(nodes)}',flush=True)
        coefficients=np.einsum('pn,nim->pim',inverse,np.array(maps))
        for delta in p['cache_validation_offsets']:
            q=anchor+delta;Me,D,C=exact(t,q,u,sigma,phi,der,p)
            approx=np.einsum('p,pim->im',features(q,anchor,radius,exponents),coefficients);approx=np.linalg.solve(approx@D,approx)
            loss=float(np.trace(approx@C@approx.T)/np.trace(Me@C@Me.T)-1)
            checks.append(dict(level=level,initial=u,offset=delta,variance_loss=loss))
        print(f'Polynomial instrument probe {level}: maximum loss {max(c["variance_loss"] for c in checks):.6g}',flush=True)
    write(out/'probe.json',dict(degree=5,radius=radius,nodes=512,node_seed=3800000,checks=checks,maximum_variance_loss=max(c['variance_loss'] for c in checks)))
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
