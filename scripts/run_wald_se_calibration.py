"""Frozen same-path comparison of nominal, random-design sandwich and bootstrap SE."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_support_wald import integrate
from scripts.run_experiments import source_hash
from scripts.run_efficiency_followup import input_digest
from src.support_information import POWERS,TARGET
from src.polynomial_gls import moments,covariance,solve
from src.weak_efficiency import nested_bumps
from src.systems import trajectory
from src.noise import ou_noise
CORE=[1,9,10,15,TARGET]


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def estimator(t,y,phi,der,phi0,der0,sigma,tau,initial,parent):
    A0,b0=moments(t,y,phi0,der0,POWERS,sigma)
    raw=np.linalg.lstsq(A0[:,CORE],b0,rcond=None)[0];low,high=np.array(parent['pilot_coefficient_box']).T
    pilot=np.clip(raw,low,high);model=integrate(t,pilot,initial)
    coefficients=np.zeros(18);coefficients[CORE]=pilot
    A,b=moments(t,y,phi,der,POWERS,sigma)
    Am,bm=moments(t,model,phi,der,POWERS,np.zeros(2));b=b-(bm-Am@coefficients)
    C=covariance(t,model,coefficients.reshape(2,9).T,phi,der,POWERS,sigma,tau)
    q,V,diag=solve(A,b,C,parent['covariance_relative_relaxation'])
    return q,V,A,C,pilot,diag


def sandwich_target(t,mean,q0,A,C,phi,der,sigma,tau,p,seed,draws):
    level=p['parent']['covariance_relative_relaxation']*max(float(np.trace(C)/len(C)),np.finfo(float).tiny)
    factor=cho_factor(C+level*np.eye(len(C)),lower=True,check_finite=False)
    scales=np.linalg.norm(A,axis=0);As=A/scales
    H=As.T@cho_solve(factor,As,check_finite=False)
    e=np.zeros(18);e[TARGET]=1/scales[TARGET]
    a=np.linalg.solve(H,e)/scales
    Am,bm=moments(t,mean,phi,der,POWERS,np.zeros(2));defect=bm-Am@q0
    values=[]
    for draw in range(draws):
        rng=np.random.default_rng(p['score_seed_base']+p['resampling_seed_stride']*seed+draw)
        y=mean+ou_noise(len(t),t[1]-t[0],tau,sigma,rng)
        Ast,b=moments(t,y,phi,der,POWERS,sigma);r=b-defect-Ast@q0
        values.append(float((Ast@a)@cho_solve(factor,r,check_finite=False)))
    values=np.array(values)
    return float(np.std(values,ddof=1)),{'score_influence_mean':float(values.mean()),'score_draws':draws}


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Research simulations execute on Kaggle')
    p=json.loads((ROOT/'configs/wald_se_calibration.json').read_text());parent=json.loads((ROOT/'configs'/p['parent_protocol_file']).read_text())
    acquisition=json.loads((ROOT/'configs'/p['observation_protocol_file']).read_text())
    refs=json.loads((ROOT/'configs'/p['reference_file']).read_text());smoke=stage=='secal_smoke'
    shard=0 if smoke else int(stage.rsplit('_',1)[1]);assert 0<=shard<p['shards']
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':p,'source_sha256':source_hash(),'smoke':smoke})
    t=np.linspace(0,acquisition['T'],round(acquisition['T']/acquisition['h'])+1)
    base=np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0);sigma=p['level']*base
    truth=np.array([1.,-1.,1.,-1.,0.]);clean=integrate(t,truth,acquisition['initial'])
    phi,der,meta=nested_bumps(t,p['primary_windows'],parent['basis_relative_cutoff']);phi0,der0,_=nested_bumps(t,parent['pilot_windows'])
    qtruth=np.zeros(18);qtruth[CORE]=truth
    A,b=moments(t,clean,phi,der,POWERS,np.zeros(2));b-=b-A@qtruth
    C=covariance(t,clean,qtruth.reshape(2,9).T,phi,der,POWERS,sigma,acquisition['tau'])
    error=abs(solve(A,b,C,parent['covariance_relative_relaxation'])[0][TARGET]);assert error<=p['clean_mean_quadrature_gate']
    write(out/'quadrature_gate.json',{'candidate_absolute_error':float(error),'threshold':p['clean_mean_quadrature_gate']});write(out/'basis.json',meta)
    ref={r['seed']:r for r in refs['rows']};indices=[p['numerical_smoke']['seed']] if smoke else [s for s in range(p['paths']) if s%p['shards']==shard]
    rows=[];critical=float(norm.ppf(1-p['alpha']/2));bcount=p['numerical_smoke']['bootstrap_draws'] if smoke else p['bootstrap_draws'];scount=p['numerical_smoke']['score_draws'] if smoke else p['score_draws']
    for seed in indices:
        timer=time.perf_counter();rng=np.random.default_rng(acquisition['seed_base']+acquisition['level_seed_stride']+seed)
        y=clean+p['level']*ou_noise(len(t),acquisition['h'],acquisition['tau'],base,rng)
        digest=input_digest(t,y);assert digest==ref[seed]['input_sha256']
        estimate,V,A,C,pilot,diag=estimator(t,y,phi,der,phi0,der0,sigma,acquisition['tau'],acquisition['initial'],parent)
        nominal=float(np.sqrt(V[TARGET,TARGET]));assert np.max(np.abs(estimate-np.array(ref[seed]['coefficients'])))<=p['coefficient_replay_absolute_tolerance']
        assert abs(nominal/ref[seed]['standard_error']-1)<=p['se_replay_relative_tolerance']
        nullpilot=pilot.copy();nullpilot[-1]=0.;mean=integrate(t,nullpilot,acquisition['initial']);q0=np.zeros(18);q0[CORE]=nullpilot
        score_start=time.perf_counter();sandwich,details=sandwich_target(t,mean,q0,A,C,phi,der,sigma,acquisition['tau'],{**p,'parent':parent},seed,scount)
        score_seconds=time.perf_counter()-score_start
        boot=[];boot_start=time.perf_counter()
        for draw in range(bcount):
            rng=np.random.default_rng(p['bootstrap_seed_base']+p['resampling_seed_stride']*seed+draw)
            yy=mean+ou_noise(len(t),acquisition['h'],acquisition['tau'],sigma,rng)
            q,*_=estimator(t,yy,phi,der,phi0,der0,sigma,acquisition['tau'],acquisition['initial'],parent)
            boot.append(float(q[TARGET]))
            if (draw+1)%50==0:print(f'{stage}: seed {seed}, bootstrap {draw+1}/{bcount}',flush=True)
        boot=np.array(boot);bootstrap=float(np.std(boot,ddof=1))
        assert all(np.isfinite(v) and v>0 for v in [nominal,sandwich,bootstrap])
        for method,se in [('nominal',nominal),('design_noise_sandwich',sandwich),('parametric_bootstrap',bootstrap)]:
            rows.append({'seed':seed,'method':method,'candidate_coefficient':float(estimate[TARGET]),'standard_error':se,
                'z_score':float(estimate[TARGET]/se),'rejected':bool(abs(estimate[TARGET]/se)>critical),'input_sha256':digest})
        write(out/'rows.json',rows)
        write(out/f'diagnostics_s{seed:03d}.json',{'seed':seed,'pilot':pilot.tolist(),'null_generating_pilot':nullpilot.tolist(),
            'bootstrap_targets':boot.tolist(),'bootstrap_mean':float(boot.mean()),'bootstrap_draws':bcount,
            'score_seconds':score_seconds,'bootstrap_seconds':time.perf_counter()-boot_start,
            'total_seconds':time.perf_counter()-timer,**details,**diag})
        print(f'{stage}: completed seed {seed}; {len(rows)//3}/{len(indices)} paired records',flush=True)
    write(out/'completion.json',{'status':'complete','paths':len(indices),'rows':len(rows),'elapsed_seconds':time.perf_counter()-start,'smoke':smoke})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();run(a.stage,a.output)
