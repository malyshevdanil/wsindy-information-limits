"""Prespecified paired power study and matched author WENDy, on Kaggle."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,io,json,sys,time,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_efficiency_validation import write
from src.weak_efficiency import (nested_bumps,vdp_moments,vdp_pilot_trajectory,covariance_blocks,
    fit_covariance_iteration,scalar_generalized_test,scalar_optimal_test)
from src.systems import trajectory
from src.noise import ou_noise
from src.weak_features import test_functions


def generate(p):
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1)
    return t,trajectory('vdp',t),np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0)


def input_digest(t,y):return hashlib.sha256(t.tobytes()+y.tobytes()).hexdigest()


def curves(out):
    p=json.loads((ROOT/'configs/efficiency.json').read_text(encoding='utf-8'))['scalar']
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1)
    args=[t,p['a'],p['x0'],p['T'],p['tau'],.1]
    exact,der=scalar_generalized_test(*args)
    smooth,sd=scalar_optimal_test(*args,.1)
    bell,bd=test_functions(t,[2.],3.)
    # Display all three on the same clean-moment normalization.
    x=p['x0']*np.exp(p['a']*t)
    bell=bell[0]/np.trapezoid(bell[0]*x,t)
    smooth=smooth/np.trapezoid(smooth*x,t)
    np.savez_compressed(out/'optimal_test_curves.npz',t=t,exact=exact,smooth=smooth,bell=bell,
        endpoint_traces=exact[[0,-1]],smoothing_width=.1)


def power(out,p,shard):
    t,x,base=generate(p);truth=np.array([1.,-1.,1.,-1.]);sigma=.005*base
    phi0,der0,_=nested_bumps(t,30,p['basis_relative_cutoff'])
    phi,der,meta=nested_bumps(t,480,p['basis_relative_cutoff']);rows=[]
    if shard==0:curves(out)
    for seed in range(p['paths']):
        if seed%p['shards']!=shard:continue
        y=x+.005*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed))
        digest=input_digest(t,y)
        A0,b0=vdp_moments(t,y,phi0,der0,sigma);pilot=np.linalg.lstsq(A0,b0,rcond=None)[0]
        model,constrained=vdp_pilot_trajectory(t,pilot,p['initial'],p['pilot_box'])
        A,b=vdp_moments(t,y,phi,der,sigma)
        for weight,tau in [('ou',p['tau']),('iid',None)]:
            start=time.perf_counter();blocks,H=covariance_blocks(t,model,phi,der,sigma,tau)
            q,diag=fit_covariance_iteration(A,b,constrained,blocks,H,sigma,p['iterations'],p['tolerance'],p['covariance_relative_relaxation'])
            rows.append(dict(seed=seed,level=.005,windows=480,weight_model=weight,input_sha256=digest,
                coefficients=q.tolist(),relative_error=float(np.linalg.norm(q-truth)/2),pilot=pilot.tolist(),
                pilot_clipped=bool(np.any(pilot!=constrained)),seconds=time.perf_counter()-start,**meta,**diag))
        if len(rows)%20==0:write(out/'rows.json',rows);print(f'Power shard={shard}: {len(rows)//2}/500 paired paths',flush=True)
    write(out/'rows.json',rows)


def author(out,p,shard=None):
    from scripts.run_full_adaptive import fetch,sha,run_octave,canonical
    pins=json.loads((ROOT/'configs/author_pins.json').read_text(encoding='utf-8'));pin=pins['WENDy']
    data=fetch(f"https://codeload.github.com/{pin['repo']}/zip/{pin['commit']}")
    assert sha(data)==pin['zip_sha256']
    dest=ROOT/'author_sources/WENDy';dest.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        for entry in archive.infolist():
            path=Path(entry.filename)
            if len(path.parts)==2 and path.suffix in ('.m','.md'):(dest/path.name).write_bytes(archive.read(entry))
    write(out/'author_provenance.json',{'WENDy':{**pin,'files_verified':{f.name:sha(f.read_bytes()) for f in dest.iterdir() if f.is_file()}}})
    t,x,base=generate(p);rows=[];inputs=out/'inputs';inputs.mkdir(exist_ok=True)
    truth=np.array([1.,-1.,1.,-1.])
    for seed in range(p['paths']):
        if shard is not None and seed%4!=shard:continue
        y=x+.05*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed))
        np.savez_compressed(inputs/f'seed_{seed:04d}.npz',t=t,y=y)
        fit=out/'fits'/f'seed_{seed:04d}';fit.mkdir(parents=True,exist_ok=True)
        start=time.perf_counter();result=run_octave('wendy',t,y,fit,{'WENDy':dest})
        full=canonical(result['coefficients'],result['tags']);q=full[[1,0,1,6],[0,1,1,1]]
        assert np.all(np.isfinite(q))
        rows.append(dict(seed=seed,level=.05,method='author_wendy',status='ok',input_sha256=input_digest(t,y),
            coefficients=q.tolist(),relative_error=float(np.linalg.norm(q-truth)/2),seconds=time.perf_counter()-start,
            basis_counts=result['basis_counts'],mts=result['mts'],sigma_estimates=result['sigma_estimates']))
        write(out/'rows.json',rows);print(f'Matched author WENDy {seed+1}/{p["paths"]}',flush=True)


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical experiments must run on Kaggle')
    cfg=json.loads((ROOT/'configs/efficiency_followup.json').read_text(encoding='utf-8'))
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',dict(stage=stage,protocol=cfg,source_sha256=source_hash(),
        compatibility_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'src/author_compat').glob('*.m')}))
    if stage.startswith('efficiency_match') or stage=='efficiency_matched_wendy':author(out,cfg['matched'],int(stage[-1]) if stage.startswith('efficiency_match_') else None)
    else:power(out,cfg['power'],int(stage[-1]))
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();run(args.stage,args.output)
