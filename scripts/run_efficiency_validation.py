"""Execute the two efficiency experiments on Kaggle only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from src.weak_efficiency import (scalar_information,scalar_optimal_test,nested_bumps,vdp_pilot_trajectory,
    vdp_moments,covariance_blocks,assemble_covariance,gls_map,fit_covariance_iteration)
from src.weak_features import test_functions,quadrature_weights
from src.covariance import exponential_apply,residual_covariance
from src.systems import trajectory,system
from src.noise import ou_noise


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')


def scalar(out,p):
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1)
    n=len(t)-1;assert n%2==0
    w=np.ones(len(t));w[1:-1:2]=4;w[2:-1:2]=2;w*=p['h']/3
    x=p['x0']*np.exp(p['a']*t);phi,der=test_functions(t,[2.],3.)
    pf,pd=test_functions(t,np.linspace(.805,3.195,30),1.6)
    rows=[];reference=[]
    for sigma in p['sigmas']:
        known,_,_,eff=scalar_information(p['a'],p['x0'],p['T'],p['tau'],sigma)
        tests={}
        for delta in p['smoothing_widths']:
            ph,dp=scalar_optimal_test(t,p['a'],p['x0'],p['T'],p['tau'],sigma,delta)
            D=float((ph*w)@x);ell=((-dp-p['a']*ph)*w)/D
            V=float(sigma*sigma*ell@exponential_apply(ell,np.exp(-p['h']/p['tau'])))
            annihilation=float(ell@(x/p['x0']))
            assert abs(annihilation)<1e-6, f'Annihilation delta={delta}: {annihilation}'
            assert V*eff>=1-1e-5
            reference.append({'sigma':sigma,'delta':delta,'variance':V,'weak_variance_bound':1/eff,
                              'known_variance_bound':1/known,'weak_efficiency':1/(V*eff),
                              'initial_annihilation':annihilation,'clean_quadrature_error':float(-(dp*w)@x/D-p['a'])})
            tests[delta]=(ph,dp,D)
        for seed in range(p['paths']):
            y=x+ou_noise(len(t),p['h'],p['tau'],[sigma],np.random.default_rng(p['seed_base']+seed))[:,0]
            G=(pf*w)@y;b=-(pd*w)@y;pilot=float(np.clip((G@b)/(G@G),-2.5,.5))
            estimate=float(-(der[0]*w)@y/((phi[0]*w)@y))
            rows.append({'sigma':sigma,'seed':seed,'method':'single_window','delta':None,'estimate':estimate,'bounded_estimate':float(np.clip(estimate,-2.5,.5)),'pilot':pilot})
            for delta,(ph,dp,D) in tests.items():
                for method in ('oracle_phi','pilot_phi'):
                    if method=='pilot_phi':ph,dp=scalar_optimal_test(t,pilot,p['x0'],p['T'],p['tau'],sigma,delta)
                    denominator=float((ph*w)@y)
                    estimate=float(-(dp*w)@y/denominator)
                    cleanD=float((ph*w)@x)
                    rows.append({'sigma':sigma,'seed':seed,'method':method,'delta':delta,
                        'estimate':estimate,'bounded_estimate':float(np.clip(estimate,-2.5,.5)),
                        'pilot':pilot,'guard_passed':abs(denominator)>=.05*abs(cleanD)})
            if seed%200==0:print(f'Scalar sigma={sigma} path={seed+1}',flush=True)
    write(out/'rows.json',rows);write(out/'scalar_reference.json',reference)


def geometry(out,p):
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);x=trajectory('vdp',t)
    calibration=trajectory('vdp',np.linspace(0,20,2001));sigma=.05*np.std(calibration,axis=0)
    truth=np.array([1.,-1.,1.,-1.]);rows=[]
    for K in p['window_counts']:
        phi,der,meta=nested_bumps(t,K,p['basis_relative_cutoff'])
        A,b=vdp_moments(t,x,phi,der,np.zeros(2))
        blocks,H=covariance_blocks(t,x,phi,der,sigma,p['tau']);OU=assemble_covariance(truth,blocks,H*0,sigma)
        blocks2,H2=covariance_blocks(t,x,phi,der,sigma,None);IID=assemble_covariance(truth,blocks2,H2*0,sigma)
        for model,C in [('ou',OU),('iid',IID)]:
            P,level=gls_map(A,C,p['covariance_relative_relaxation'])
            actual=P@OU@P.T
            row={**meta,'weight_model':model,'relative_rmse':float(np.sqrt(np.trace(actual))/2),
                 'actual_ou_covariance':actual.tolist(),'normalization_error':float(np.linalg.norm(P@A-np.eye(4))),'relaxation':level,
                 'clean_quadrature_relative_error':float(np.linalg.norm(P@b-truth)/2)}
            rows.append(row)
        print(f'Geometry K={K} retained={len(phi)}',flush=True)
    # Independent fixture compares precomputed polynomial covariance to Prop3 implementation.
    _,powers,coef=system('vdp');phi0,der0=test_functions(t,np.linspace(1.505,18.495,30),3.)
    blocks,H=covariance_blocks(t,x,phi0,der0,sigma,.3)
    C=assemble_covariance(truth,blocks,H,sigma)
    target=residual_covariance(t,x,powers,coef,np.linspace(1.505,18.495,30),3.,sigma,.3)
    discrepancy=float(np.linalg.norm(C-target)/np.linalg.norm(target));assert discrepancy<1e-10
    write(out/'geometry.json',{'results':rows,'exact_covariance_fixture_relative_error':discrepancy})


def vdp(out,p,shard):
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);x=trajectory('vdp',t)
    calibration=trajectory('vdp',np.linspace(0,20,2001));base_sigma=np.std(calibration,axis=0)
    bases={K:nested_bumps(t,K,p['basis_relative_cutoff']) for K in p['window_counts']}
    phi0,der0,_=bases[30];truth=np.array([1.,-1.,1.,-1.]);rows=[]
    for seed in range(p['paths']):
        if seed%2!=shard:continue
        # Same standardized OU path across amplitudes; algorithms are paired.
        noise=ou_noise(len(t),p['h'],p['tau'],base_sigma,np.random.default_rng(p['seed_base']+seed))
        for amplitude in p['levels']:
            sigma=amplitude*base_sigma;y=x+amplitude*noise
            A0,b0=vdp_moments(t,y,phi0,der0,sigma);pilot=np.linalg.lstsq(A0,b0,rcond=None)[0]
            model,pilot_box=vdp_pilot_trajectory(t,pilot,p['initial'],p['pilot_box'])
            for K,(phi,der,meta) in bases.items():
                A,b=vdp_moments(t,y,phi,der,sigma)
                for weight,tau in [('ou',p['tau']),('iid',None)]:
                    began=time.perf_counter();blocks,H=covariance_blocks(t,model,phi,der,sigma,tau)
                    q,diagnostics=fit_covariance_iteration(A,b,pilot_box,blocks,H,sigma,p['iterations'],p['tolerance'],p['covariance_relative_relaxation'])
                    rows.append({'seed':seed,'level':amplitude,'windows':K,'weight_model':weight,
                        'coefficients':q.tolist(),'relative_error':float(np.linalg.norm(q-truth)/2),
                        'pilot':pilot.tolist(),'pilot_clipped':bool(np.any(pilot!=pilot_box)),
                        'seconds':time.perf_counter()-began,**meta,**diagnostics})
        print(f'VdP shard={shard} seed={seed} complete',flush=True)
        write(out/'rows.json',rows)
    write(out/'rows.json',rows)


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):
        raise RuntimeError('Run numerical experiments on Kaggle')
    cfg=json.loads((ROOT/'configs/efficiency.json').read_text());start=time.perf_counter()
    out.mkdir(parents=True,exist_ok=True)
    write(out/'manifest.json',{'stage':stage,'protocol':cfg,'source_sha256':source_hash()})
    if stage.startswith('efficiency_scalar'):scalar(out,cfg['scalar'])
    elif stage in ('efficiency_geometry','efficiency_geometry_v2'):geometry(out,cfg['vdp'])
    else:vdp(out,cfg['vdp'],int(stage[-1]))
    write(out/'completion.json',{'status':'complete','elapsed_seconds':time.perf_counter()-start})
    print(f'{stage} complete',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True,choices=['efficiency_scalar','efficiency_scalar_v2','efficiency_scalar_v3','efficiency_geometry','efficiency_geometry_v2','efficiency_vdp_0','efficiency_vdp_1','efficiency_vdp_polish_0','efficiency_vdp_polish_1','efficiency_vdp_final_0','efficiency_vdp_final_1']);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.stage,args.output)
