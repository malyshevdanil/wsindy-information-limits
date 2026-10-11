"""Frozen fresh-seed validation. Run numerical experiments on Kaggle only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.noise import ou_noise
from src.systems import system,trajectory
from src.weak_features import weak_matrices,test_functions,quadrature_weights,derivative_polynomial
from src.covariance import residual_covariance
from src.mstls import relative_mstls
from scripts.run_experiments import source_hash


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def oracle(G,b,truth):
    estimate=np.zeros_like(truth)
    for eq in range(truth.shape[1]):
        active=truth[:,eq]!=0
        estimate[active,eq]=np.linalg.lstsq(G[:,active],b[:,eq],rcond=None)[0]
    return estimate


def setup(name,h):
    T,L=(4.,1.6) if name=='linear' else (20.,3.)
    t=np.linspace(0,T,round(T/h)+1)
    centers=np.linspace(L/2+.005,T-L/2-.005,30)
    _,powers,truth=system(name)
    return t,trajectory(name,t),powers,truth,centers,L


def benchmark(out,cfg):
    rows=[]
    for name in cfg['systems']:
        tf,xf,powers,truth,centers,L=setup(name,min(cfg['h']))
        for regime in cfg['regimes']:
            for level in cfg['noise_levels']:
                sigma=level*np.std(xf,axis=0)
                for seed in range(cfg['benchmark_repeats']):
                    # Equal seeds across conditions give deliberate common random numbers.
                    rng=np.random.default_rng(cfg['benchmark_seed_base']+seed)
                    e=(ou_noise(len(tf),tf[1]-tf[0],.1*L,sigma,rng) if regime=='physical_ou'
                       else rng.normal(size=xf.shape)*sigma)
                    for h in cfg['h']:
                        stride=round(h/min(cfg['h']))
                        t,x,y=tf[::stride],xf[::stride],(xf+e)[::stride]
                        Gr,b=weak_matrices(t,y,powers,centers,L)
                        Gc,_=weak_matrices(t,y,powers,centers,L,sigma**2)
                        for method,G in [('relative_mstls_raw',Gr),('relative_mstls_hermite',Gc),('oracle_hermite',Gc)]:
                            began=time.perf_counter()
                            if method=='oracle_hermite':
                                estimate=oracle(G,b,truth);choices=[]
                            else:
                                estimate,choices=relative_mstls(G,b)
                            selected=estimate!=0;target=truth!=0
                            rows.append({'system':name,'regime':regime,'level':level,'seed':seed,'h':h,
                                'method':method,'error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth)),
                                'support_correct':bool(np.array_equal(selected,target)),
                                'false_positive':int(np.sum(selected & ~target)),
                                'false_negative':int(np.sum(~selected & target)),
                                'coefficients':estimate.tolist(),'thresholds':choices,
                                'fit_seconds':time.perf_counter()-began})
                print(f'Benchmark complete: {name}, {regime}, {level}',flush=True)
    write(out/'benchmark_rows.json',rows)
    return len(rows)


def fisher(T,a=-1.,tau=.4,sigma=.1,x0=2.):
    return tau*x0*x0/(2*sigma*sigma)*quad(lambda t:np.exp(2*a*t)*(1+(a+1/tau)*t)**2,0,T,epsabs=1e-10)[0]


def scalar(out,cfg):
    a,x0,sigma,tau,T=-1.,2.,.1,.4,4.
    tf=np.linspace(0,T,4097)
    steps=[32,16,8,4,2,1]
    weak=np.zeros((cfg['scalar_repeats'],len(steps)))
    mle=np.zeros_like(weak)
    info=[]
    grids=[]
    for stride in steps:
        t=tf[::stride];h=t[1]-t[0];rho=np.exp(-h/tau)
        s=x0*t*np.exp(a*t)
        I=s[0]**2/sigma**2+np.sum((s[1:]-rho*s[:-1])**2)/(sigma**2*(-np.expm1(-2*h/tau)))
        phi,der=test_functions(t,[2.],3.)
        weights=quadrature_weights(t)
        theta_grid=np.linspace(-2.5,.5,61)
        means=x0*np.exp(theta_grid[:,None]*t[None,:])
        info.append({'h':h,'fisher_information':float(I),'cr_variance':1/float(I)})
        grids.append((t,rho,phi[0]*weights,der[0]*weights,theta_grid,means))
    for seed in range(cfg['scalar_repeats']):
        y=x0*np.exp(a*tf)+ou_noise(len(tf),tf[1]-tf[0],tau,[sigma],np.random.default_rng(cfg['scalar_seed_base']+seed))[:,0]
        for j,(stride,data) in enumerate(zip(steps,grids)):
            t,rho,p,pd,theta_grid,means=data;ys=y[::stride]
            weak[seed,j]=-pd@ys/(p@ys)
            innovation=ys[1:]-rho*ys[:-1]
            residuals=innovation[None,:]-(means[:,1:]-rho*means[:,:-1])
            scores=np.sum(residuals**2,axis=1)
            candidates=[(-2.5,float(scores[0])),(.5,float(scores[-1]))]
            for index in range(1,len(theta_grid)-1):
                if scores[index]<=scores[index-1] and scores[index]<=scores[index+1]:
                    def objective(theta):
                        m=x0*np.exp(theta*t)
                        return float(np.sum((innovation-(m[1:]-rho*m[:-1]))**2))
                    opt=minimize_scalar(objective,bounds=(theta_grid[index-1],theta_grid[index+1]),method='bounded',options={'xatol':1e-10})
                    candidates.append((float(opt.x),float(opt.fun)))
            mle[seed,j]=min(candidates,key=lambda pair:pair[1])[0]
        if seed%200==0:print(f'Scalar likelihood: {seed+1}',flush=True)
    # Linearized weak variance and information efficiency, not ratio moments.
    phi,der=test_functions(tf,[2.],3.);w=quadrature_weights(tf)
    D=float((phi[0]*w)@(x0*np.exp(a*tf)))
    q=(-der[0]-a*phi[0])*w/D
    from src.covariance import exponential_apply
    variance=float(sigma*sigma*q@exponential_apply(q,np.exp(-(tf[1]-tf[0])/tau)))
    np.savez_compressed(out/'scalar_likelihood.npz',weak=weak,mle=mle,h=np.array([row['h'] for row in info]))
    write(out/'information.json',{'scalar_parameters':{'a':a,'x0':x0,'sigma':sigma,'tau':tau,'T':T},
        'likelihood_bounds':[-2.5,.5],'grid_information':info,'continuous_information':fisher(T),
        'weak_linearized_variance':variance,'weak_linearized_efficiency':1/(fisher(T)*variance),
        'duration': [{'T':v,'decay_information':fisher(v),'constant_mean_effective_n':1+v/(2*tau)} for v in [1.,2.,4.,8.,16.,32.,64.]],
        'independent_recordings':[{'M':M,'information':M*fisher(T)} for M in [1,2,4,8,16,32]],
        'comparison_scope':'Bounded Gaussian MLE uses known noise covariance and initial condition; neither unbiasedness nor finite-sample CR attainment is asserted. Weak ratio moments are undefined.'})


def coefficients(out,cfg):
    t,x,powers,truth,centers,L=setup('vdp',.01)
    G0,b0=weak_matrices(t,x,powers,centers,L)
    d=truth.shape[1];K=len(centers)
    indexes=[np.flatnonzero(truth[:,i]) for i in range(d)]
    P=np.zeros((sum(map(len,indexes)),K*d));position=0
    clean_min=[]
    for i,S in enumerate(indexes):
        P[position:position+len(S),np.arange(i,K*d,d)]=np.linalg.pinv(G0[:,S])
        clean_min.append(np.linalg.svd(G0[:,S],compute_uv=False)[-1]);position+=len(S)
    clean_est=oracle(G0,b0,truth)
    def compact(est):return np.concatenate([est[S,i] for i,S in enumerate(indexes)])
    phi,der=test_functions(t,centers,L);w=quadrature_weights(t)
    summaries=[]
    for level in cfg['coefficient_noise_levels']:
        sigma=level*np.std(x,axis=0)
        C1=residual_covariance(t,x,powers,truth,centers,L,sigma,.3,first_order=True)
        predicted=P@C1@P.T
        actual=[];linear=[];guards=[]
        for seed in range(cfg['coefficient_repeats']):
            e=ou_noise(len(t),.01,.3,sigma,np.random.default_rng(cfg['coefficient_seed_base']+seed))
            G,b=weak_matrices(t,x+e,powers,centers,L,sigma**2)
            est=oracle(G,b,truth)
            actual.append(compact(est-clean_est))
            residual=np.zeros((K,d))
            for j in range(d):
                alpha=tuple(int(k==j) for k in range(d))
                jac=derivative_polynomial(x,powers,truth,alpha)
                residual-= (phi*w)@(jac*e[:,j,None])
                residual[:,j]-=(der*w)@e[:,j]
            linear.append(P@residual.ravel())
            guards.append(all(np.linalg.svd(G[:,S],compute_uv=False)[-1]>=.5*clean_min[i] for i,S in enumerate(indexes)))
        actual=np.array(actual);linear=np.array(linear);guards=np.array(guards)
        empirical_linear=linear.T@linear/len(linear)
        guarded=actual[guards]
        empirical_actual=np.cov(guarded,rowvar=False,bias=True)
        np.savez_compressed(out/f'coefficient_{level:g}.npz',actual=actual,linear=linear,guards=guards,predicted_covariance=predicted)
        summaries.append({'level':level,'paths':len(actual),'guarded_paths':int(guards.sum()),
            'linear_relative_frobenius':float(np.linalg.norm(empirical_linear-predicted)/np.linalg.norm(predicted)),
            'guarded_actual_relative_frobenius':float(np.linalg.norm(empirical_actual-predicted)/np.linalg.norm(predicted)),
            'guarded_bias_norm':float(np.linalg.norm(guarded.mean(axis=0))),
            'paired_linearization_relative_rms':float(np.linalg.norm(guarded-linear[guards])/np.linalg.norm(linear[guards])),
            'gaussian_known_mean_covariance_relative_rms_floor':float(np.sqrt((np.trace(predicted)**2+np.sum(predicted**2))/len(linear))/np.linalg.norm(predicted)),
            'scope':'Actual covariance conditional on singular-value guard; no unconditional moment claim. Errors centred at discrete clean fit.'})
        print(f'Coefficient surrogate: {level}',flush=True)
    write(out/'coefficient_summary.json',summaries)


def covariance_null(out,cfg):
    summaries=[]
    source=ROOT/'inputs/covariance'
    for path in sorted(source.glob('residual_*.npz')):
        data=np.load(path)
        R=data['centered_residuals'];C=data['exact_covariance'];n,p=R.shape
        S=R.T@R/n
        fourth=np.mean(np.sum(R*R,axis=1)**2)
        rms=np.sqrt(max(0.,fourth-np.sum(C*C))/n)/np.linalg.norm(C)
        gaussian_rms=np.sqrt((np.trace(C)**2+np.sum(C*C))/n)/np.linalg.norm(C)
        # Affine null bootstrap: force population covariance C, retain approximate shape.
        vals,U=np.linalg.eigh(S);target_vals,V=np.linalg.eigh(C)
        floor=max(vals.max()*1e-12,1e-30)
        whiten=U@np.diag(1/np.sqrt(np.maximum(vals,floor)))@U.T
        root=V@np.diag(np.sqrt(np.maximum(target_vals,0)))@V.T
        Z=(R-R.mean(axis=0))@whiten@root
        # Rewhiten centred covariance to avoid a tiny mean-centering mismatch.
        zS=Z.T@Z/n; zv,zU=np.linalg.eigh(zS)
        Z=Z@(zU@np.diag(1/np.sqrt(np.maximum(zv,max(zv.max()*1e-12,1e-30))))@zU.T)@root
        rng=np.random.default_rng(cfg['null_seed_base'])
        discrepancy=[]
        for _ in range(cfg['null_covariance_repeats']):
            counts=rng.multinomial(n,np.full(n,1/n))
            simulated=Z.T@(counts[:,None]*Z)/n
            discrepancy.append(float(np.linalg.norm(simulated-C)/np.linalg.norm(C)))
        observed=float(np.linalg.norm(S-C)/np.linalg.norm(C))
        summaries.append({'file':path.name,'n':n,'dimension':p,'observed_relative_frobenius':observed,
            'estimated_fourth_moment_relative_rms_floor':float(rms),'gaussian_reference_relative_rms_floor':float(gaussian_rms),
            'affine_null_bootstrap_quantiles':np.quantile(discrepancy,[.025,.5,.975]).tolist(),
            'bootstrap_exceedance_fraction':float(np.mean(np.array(discrepancy)>=observed)),
            'interpretation':'Approximate covariance-null calibration using empirical non-Gaussian shape; Gaussian floor is exact only for Gaussian residuals; exceedance fraction is not an exact test p-value.'})
        print(f'Covariance null: {path.name}',flush=True)
    if len(summaries)!=6:raise ValueError('Expected six archived covariance input arrays')
    write(out/'covariance_null.json',summaries)


def heterogeneous(out,cfg):
    t,x,powers,truth,centers,L=setup('vdp',.02)
    sigma=.05*np.std(x,axis=0)
    rows=[]
    for block in range(cfg['heterogeneous_blocks']):
        matrices=[];records=[]
        for r in range(max(cfg['heterogeneous_M'])):
            rng=np.random.default_rng(cfg['heterogeneous_seed_base']+16*block+r)
            initial=np.array([2.,0.])+rng.uniform(-.25,.25,size=2)
            clean=trajectory('vdp',t,initial)
            y=clean+ou_noise(len(t),.02,.3,sigma,rng)
            matrices.append(weak_matrices(t,y,powers,centers,L,sigma**2));records.append(y)
            M=r+1
            if M not in cfg['heterogeneous_M']:continue
            pairs=[(np.vstack([p[0] for p in matrices]),np.vstack([p[1] for p in matrices])),
                   (np.mean([p[0] for p in matrices],axis=0),np.mean([p[1] for p in matrices],axis=0)),
                   weak_matrices(t,np.mean(records,axis=0),powers,centers,L,sigma**2/M)]
            for method,(G,b) in zip(('stacked','mean_weak','naive_mean_sensor'),pairs):
                estimate=oracle(G,b,truth)
                rows.append({'block':block,'M':M,'method':method,'error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth))})
        if block%100==0:print(f'Heterogeneous block: {block+1}',flush=True)
    write(out/'heterogeneous_rows.json',rows)


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/'configs/revision.json').read_text())
    start=time.perf_counter()
    write(out/'manifest.json',{'stage':'revision','protocol':cfg,'source_sha256':source_hash(),
        'input_note':'Residual arrays from earlier exploratory study reused only to calibrate its covariance discrepancy.'})
    fits=benchmark(out,cfg)
    scalar(out,cfg)
    coefficients(out,cfg)
    covariance_null(out,cfg)
    heterogeneous(out,cfg)
    write(out/'completion.json',{'status':'complete','benchmark_fits':fits,'elapsed_seconds':time.perf_counter()-start})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
