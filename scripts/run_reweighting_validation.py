import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_revision_validation import setup,oracle,write
from scripts.run_experiments import source_hash
from src.weak_features import weak_matrices
from src.weak_irls import fit_irls
from src.noise import ou_noise

def run(out):
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    cfg=json.loads((ROOT/'configs/reweighting.json').read_text())
    write(out/'manifest.json',{'stage':'reweighting','protocol':cfg,'source_sha256':source_hash()})
    t,x,powers,truth,centers,L=setup('vdp',.005);sigma=.05*np.std(x,axis=0)
    rows=[]
    for regime in cfg['regimes']:
        for seed in range(cfg['paths']):
            rng=np.random.default_rng(cfg['seed_base']+seed)
            e=ou_noise(len(t),.005,.3,sigma,rng) if regime=='physical_ou' else rng.normal(size=x.shape)*sigma
            for h in cfg['h']:
                stride=round(h/.005);ts=t[::stride];y=(x+e)[::stride]
                G,b=weak_matrices(ts,y,powers,centers,L)
                for method in cfg['methods']:
                    began=time.perf_counter()
                    if method=='oracle_raw_ols':estimate=oracle(G,b,truth);meta={}
                    else:estimate,meta=fit_irls(ts,y,powers,truth!=0,centers,L,sigma,
                        tau=.3 if method=='ou_irls_core' else None,
                        iterations=cfg['iterations'],tolerance=cfg['tolerance'])
                    rows.append({'regime':regime,'seed':seed,'h':h,'method':method,
                        'error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth)),
                        'coefficients':estimate.tolist(),'fit_seconds':time.perf_counter()-began,**meta})
            if seed%32==0:print(f'{regime} path {seed+1}',flush=True)
    write(out/'rows.json',rows)
    write(out/'completion.json',{'status':'complete','fits':len(rows),'elapsed_seconds':time.perf_counter()-start})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
