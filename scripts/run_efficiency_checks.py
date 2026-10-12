"""Follow-up checks with preserved one-step results; execute only on Kaggle."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_efficiency_followup import generate,input_digest
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash
from src.path_scoring import one_step_vdp,ou_whiten
from src.noise import ou_noise
from src.weak_features import test_functions
from src.covariance import residual_covariance
from scripts.run_lorenz_efficiency import blocks,assemble

def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Kaggle only')
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    manifest=dict(stage=stage,source_sha256=source_hash())
    if stage=='efficiency_fixture':
        t=np.linspace(0,4,401);x=np.column_stack([np.sin(t)+1,np.cos(t),np.sin(2*t)+2]);q=np.array([10.,28.,8/3]);sig=np.array([.03,.04,.05]);tau=.1
        centers=[1.,2.,3.];phi,der=test_functions(t,centers,1.5)
        powers=[(1,0,0),(0,1,0),(0,0,1),(1,0,1),(1,1,0)]
        coef=np.array([[-q[0],q[1],0],[q[0],-1,0],[0,0,-q[2]],[0,-1,0],[0,0,1.]])
        exact=residual_covariance(t,x,powers,coef,centers,1.5,sig,tau)
        fast=assemble(q,blocks(t,x,phi,der,sig,tau));error=np.linalg.norm(exact-fast)/np.linalg.norm(exact)
        assert error<1e-12
        N=19;h=.01;rho=np.exp(-h/tau);C=sig[0]**2*rho**np.abs(np.arange(N)[:,None]-np.arange(N));V=np.arange(N*2,dtype=float).reshape(N,2)/N
        W=ou_whiten(V[:,None,:],h,tau,np.array([sig[0]])).reshape(N,2)
        whitening=np.linalg.norm(W.T@W-V.T@np.linalg.solve(C,V))/np.linalg.norm(W.T@W);assert whitening<1e-12
        write(out/'checks.json',dict(lorenz_exact_covariance_relative_error=float(error),ou_precision_whitening_relative_error=float(whitening)))
    else:
        p=json.loads((ROOT/'configs/information_followup.json').read_text(encoding='utf-8'))['fisher'];shard=int(stage[-1]);t,x,base=generate(p);rows=[]
        fixture=ROOT/'configs/refinement_pilots.npz';manifest['pilot_sha256']=hashlib.sha256(fixture.read_bytes()).hexdigest()
        manifest['protocol']=dict(paths=2000,shards=2,additional_steps=[2,3],initialization='recorded unknown-initial one-step output; follow-up after inspecting its bias',**{k:v for k,v in p.items() if k!='paths'})
        with np.load(fixture) as d:
            for seed in range(shard,2000,2):
                y=x+.005*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed));digest=input_digest(t,y);assert digest==d['low_hash'][seed].decode()
                q=d['joint_q'][seed].copy();u=d['joint_u'][seed].copy()
                for step in [2,3]:
                    q,diag=one_step_vdp(t,y,q,u,p['tau'],.005*base,True);u=np.array(diag['initial_estimate'])
                    rows.append(dict(seed=seed,step=step,coefficients=q.tolist(),relative_error=float(np.linalg.norm(q-np.array([1,-1,1,-1]))/2),input_sha256=digest,**diag))
                if seed%100<2:write(out/'rows.json',rows);print(f'Refinement shard {shard}: {seed}/2000',flush=True)
        write(out/'rows.json',rows)
    write(out/'manifest.json',manifest);write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--output',required=True,type=Path);args=p.parse_args();run(args.stage,args.output)
