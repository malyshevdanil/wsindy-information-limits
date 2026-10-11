"""Same frozen protocol, with batched deterministic VdP integration only.

No statistical settings or seeds are changed. Batching avoids thousands of
Python polynomial-library callbacks for the heterogeneous-recording endpoint.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import scripts.run_revision_validation as revision
from src.noise import ou_noise
from src.weak_features import weak_matrices


def batch_vdp(t,initials):
    def rhs(_,state):
        out=np.empty_like(state);x=state[::2];v=state[1::2]
        out[::2]=v;out[1::2]=v-x-x*x*v
        return out
    solution=solve_ivp(rhs,(float(t[0]),float(t[-1])),np.asarray(initials).ravel(),
        t_eval=t,method='DOP853',rtol=1e-10,atol=1e-12)
    if not solution.success:raise RuntimeError('Batched ODE integration failed')
    return solution.y.T.reshape(len(t),len(initials),2)


def heterogeneous(out,cfg):
    t,x,powers,truth,centers,L=revision.setup('vdp',.02);sigma=.05*np.std(x,axis=0)
    rows=[]
    for block in range(cfg['heterogeneous_blocks']):
        initials=[];noises=[]
        for r in range(max(cfg['heterogeneous_M'])):
            rng=np.random.default_rng(cfg['heterogeneous_seed_base']+16*block+r)
            initials.append(np.array([2.,0.])+rng.uniform(-.25,.25,size=2))
            noises.append(ou_noise(len(t),.02,.3,sigma,rng))
        clean=batch_vdp(t,initials);matrices=[];records=[]
        for r in range(max(cfg['heterogeneous_M'])):
            y=clean[:,r,:]+noises[r]
            matrices.append(weak_matrices(t,y,powers,centers,L,sigma**2));records.append(y);M=r+1
            if M not in cfg['heterogeneous_M']:continue
            pairs=[(np.vstack([p[0] for p in matrices]),np.vstack([p[1] for p in matrices])),
                (np.mean([p[0] for p in matrices],axis=0),np.mean([p[1] for p in matrices],axis=0)),
                weak_matrices(t,np.mean(records,axis=0),powers,centers,L,sigma**2/M)]
            for method,(G,b) in zip(('stacked','mean_weak','naive_mean_sensor'),pairs):
                estimate=revision.oracle(G,b,truth)
                rows.append({'block':block,'M':M,'method':method,'error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth))})
        if block%100==0:print(f'Batched heterogeneous block: {block+1}',flush=True)
    revision.write(out/'heterogeneous_rows.json',rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    revision.heterogeneous=heterogeneous
    revision.run(parser.parse_args().output)
