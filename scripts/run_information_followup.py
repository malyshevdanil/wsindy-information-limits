"""Paired likelihood scoring and unknown-covariance checks, Kaggle only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_efficiency_followup import generate,input_digest
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash
from src.noise import ou_noise
from src.path_scoring import vdp_sensitivities,one_step_vdp,ou_whiten,profile_ou_shape
from src.weak_efficiency import nested_bumps,vdp_moments,covariance_blocks,fit_covariance_iteration,vdp_pilot_trajectory


def pilot_rows():
    with np.load(ROOT/'configs/information_pilots.npz') as data:
        result={}
        for prefix,level in [('low',.005),('high',.05)]:
            for seed,q in enumerate(data[prefix+'_q']):
                result[seed,level]=dict(coefficients=q,pilot=data[prefix+'_pilot'][seed],
                    input_sha256=data['low_hash'][seed].decode() if prefix=='low' else None)
        return result



def fisher(out,p):
    t,x,base=generate(p);pilots=pilot_rows();truth=np.array([1.,-1.,1.,-1.]);rows=[]
    _,U=vdp_sensitivities(t,truth,p['initial']);bounds={}
    for unknown in [False,True]:
        V=U if unknown else U[:,:,:4];W=ou_whiten(V,p['h'],p['tau'],.005*base).reshape(-1,V.shape[-1])
        scale=np.linalg.norm(W,axis=0);C=np.linalg.inv((W/scale).T@(W/scale))/scale[:,None]/scale[None,:]
        bounds['unknown' if unknown else 'known']=dict(covariance=C[:4,:4].tolist(),relative_rmse=float(np.sqrt(np.trace(C[:4,:4]))/2))
    write(out/'discrete_information.json',bounds)
    for seed in range(p['paths']):
        y=x+.005*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed));digest=input_digest(t,y)
        pilot=pilots[seed,.005];assert digest==pilot['input_sha256'];q=np.array(pilot['coefficients'])
        for unknown in [False,True]:
            initial=y[0] if unknown else p['initial']
            estimate,diag=one_step_vdp(t,y,q,initial,p['tau'],.005*base,unknown)
            rows.append(dict(seed=seed,level=.005,method='fisher_unknown_initial' if unknown else 'fisher_known_initial',
                coefficients=estimate.tolist(),relative_error=float(np.linalg.norm(estimate-truth)/2),input_sha256=digest,**diag))
        if seed%50==0:write(out/'rows.json',rows);print(f'Fisher scoring {seed+1}/{p["paths"]}',flush=True)
    write(out/'rows.json',rows)


def noise_sensitivity(out,p):
    t,x,base=generate(p);pilots=pilot_rows();truth=np.array([1.,-1.,1.,-1.]);rows=[]
    phi,der,meta=nested_bumps(t,480,p['basis_relative_cutoff'])
    for seed in range(p['paths']):
        noise=ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed))
        for level in p['levels']:
            y=x+level*noise;pilot=pilots[seed,level];q=np.array(pilot['coefficients'])
            refined,_=one_step_vdp(t,y,q,p['initial'],p['tau'],level*base)
            fitted_mean,_=vdp_sensitivities(t,refined,p['initial'])
            fitted_tau,fitted_sigma=profile_ou_shape(y-fitted_mean,p['h'],p['noise_tau_fit_bounds'])
            model,pilot_box=vdp_pilot_trajectory(t,np.array(pilot['pilot']),p['initial'],p['pilot_box'])
            for name,tau,sigma in [('tau_half',p['tau']/2,level*base),('tau_double',2*p['tau'],level*base),('fitted_shape',fitted_tau,fitted_sigma)]:
                A,b=vdp_moments(t,y,phi,der,sigma);blocks,H=covariance_blocks(t,model,phi,der,sigma,tau)
                estimate,diag=fit_covariance_iteration(A,b,pilot_box,blocks,H,sigma,p['iterations'],p['tolerance'],p['covariance_relative_relaxation'])
                rows.append(dict(seed=seed,level=level,method=name,coefficients=estimate.tolist(),
                    relative_error=float(np.linalg.norm(estimate-truth)/2),estimated_tau=fitted_tau,
                    estimated_sigma=fitted_sigma.tolist(),used_tau=tau,used_sigma=sigma.tolist(),**meta,**diag))
        if seed%5==0:write(out/'rows.json',rows);print(f'Noise sensitivity {seed+1}/{p["paths"]}',flush=True)
    write(out/'rows.json',rows)


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Run numerical experiments on Kaggle')
    cfg=json.loads((ROOT/'configs/information_followup.json').read_text(encoding='utf-8'));out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    inputs=[ROOT/'configs/information_pilots.npz']
    write(out/'manifest.json',dict(stage=stage,protocol=cfg,source_sha256=source_hash(),
        pilot_files={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}))
    if stage.startswith('efficiency_fisher'):fisher(out,cfg['fisher'])
    else:noise_sensitivity(out,cfg['noise_sensitivity'])
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();run(args.stage,args.output)
