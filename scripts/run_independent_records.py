"""Fixed-noise repeated-record GMM and coverage; simulations run only on Kaggle."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.record_gmm import corrected_design,fit_mean
from src.weak_features import test_functions,quadrature_weights
from src.weak_efficiency import vdp_pilot_trajectory,vdp_moments,covariance_blocks,assemble_covariance,gls_map
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def summaries(y,P,D):
    R,N,d=y.shape;K=len(P);u,v=y[:,:,0],y[:,:,1]
    A=np.zeros((R,2*K,4));A[:,0::2,0]=v@P.T
    A[:,1::2,1]=u@P.T;A[:,1::2,2]=v@P.T;A[:,1::2,3]=(u*u*v)@P.T
    b=np.empty((R,2*K));b[:,0::2]=-u@D.T;b[:,1::2]=-v@D.T
    difference=(y[0::2]-y[1::2])/np.sqrt(2)
    vs=(difference*difference).mean(axis=1)
    lag=np.sum(difference[:,:-1]*difference[:,1:],axis=(1,2))
    previous=np.sum(difference[:,:-1]**2,axis=(1,2))
    return A,b,vs,lag,previous


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical runs must execute on Kaggle')
    p=json.loads((ROOT/'configs/independent_records.json').read_text(encoding='utf-8'));out.mkdir(parents=True,exist_ok=True)
    shard=int(stage.rsplit('_',1)[1]);start=time.perf_counter()
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash()))
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);qtrue=np.array(p['truth'])
    x,_=vdp_pilot_trajectory(t,qtrue,p['initial'],[[-3,3]]*4)
    sigma=p['noise_level']*np.std(x,axis=0);variance=sigma*sigma;rho=np.exp(-p['h']/p['tau'])
    phi,der=test_functions(t,p['test_centers'],p['test_width'],p['test_order']);w=quadrature_weights(t);P=phi*w;D=der*w
    cleanA,cleanb=vdp_moments(t,x,phi,der,np.zeros(2))
    blocks,higher=covariance_blocks(t,x,phi,der,sigma,p['tau']);C=assemble_covariance(qtrue,blocks,higher,sigma)
    C1=assemble_covariance(qtrue,blocks,np.zeros_like(higher),sigma)
    mapping,_=gls_map(cleanA,C,0);pseudo=mapping@cleanb
    geometry=dict(sigma=sigma.tolist(),clean_weak_defect=(cleanb-cleanA@qtrue).tolist(),
        clean_gls_relative_error=float(np.linalg.norm(pseudo-qtrue)/np.linalg.norm(qtrue)),
        full_covariance=C.tolist(),first_order_covariance=C1.tolist(),
        higher_order_trace_fraction=float(np.trace(C-C1)/np.trace(C)),
        known_variance_gmm_sqrt_n_relative_sd=float(np.sqrt(np.trace(mapping@C@mapping.T))/np.linalg.norm(qtrue)))
    write(out/'geometry.json',geometry)
    rows=[];hashes=[];total=0;sumr=np.zeros(len(C));sumrr=np.zeros_like(C)
    repetitions=list(range(shard,p['monte_carlo_repetitions'],p['shards']))
    for rep in repetitions:
        rng=np.random.default_rng(p['seed_base']+rep);R=max(p['record_counts']);N=len(t)
        innovations=rng.normal(size=(R,N,2))*sigma*np.sqrt(1-rho*rho)
        innovations[:,0]/=np.sqrt(1-rho*rho)
        y=x+lfilter([1],[1,-rho],innovations,axis=1)
        raw,b,vs,lag,previous=summaries(y,P,D)
        digest=hashlib.sha256()
        for a in (raw,b,vs,lag,previous):digest.update(np.ascontiguousarray(a).tobytes())
        hashes.append(dict(repetition=rep,seed=p['seed_base']+rep,sufficient_statistics_sha256=digest.hexdigest()))
        trueA=corrected_design(raw,variance[0]);r=b-np.einsum('rmp,p->rm',trueA,qtrue)
        total+=R;sumr+=r.sum(axis=0);sumrr+=r.T@r
        rawpair=(raw[0::2]+raw[1::2])/2;bpair=(b[0::2]+b[1::2])/2
        for n in p['record_counts']:
            B=n//2;vhat=vs[:B].mean(axis=0);rhohat=float(lag[:B].sum()/previous[:B].sum())
            if not 0<rhohat<1:raise RuntimeError('Unplanned covariance boundary event')
            tauhat=-p['h']/np.log(rhohat)
            for method,noise_variance,estimated in [('gmm_known',variance,False),('gmm_estimated',vhat,True)]:
                q,V,diagnostic=fit_mean(rawpair[:B],bpair[:B],vs[:B],noise_variance,estimated)
                sd=np.sqrt(np.diag(V));error=q-qtrue
                joint=float(error@np.linalg.solve(V,error));z=np.abs(error/sd)
                assert np.all(np.isfinite(q)) and np.linalg.eigvalsh(V).min()>0
                rows.append(dict(repetition=rep,n=n,method=method,coefficients=q.tolist(),standard_errors=sd.tolist(),
                    covariance=V.tolist(),marginal_covered=(z<=1.959963984540054).tolist(),
                    joint_covered=bool(joint<=chi2.ppf(.95,4)),joint_statistic=joint,
                    estimated_variances=vhat.tolist(),estimated_tau=tauhat,**diagnostic))
            As=trueA[:n];bs=b[:n]
            stacked=np.linalg.solve(np.einsum('rmp,rmq->pq',As,As),np.einsum('rmp,rm->p',As,bs))
            rows.append(dict(repetition=rep,n=n,method='stacked_known',coefficients=stacked.tolist()))
        if rep==repetitions[0]:np.savez_compressed(out/'sufficient_statistics_fixture.npz',raw=raw,b=b,vs=vs,lag=lag,previous=previous)
        if (len(hashes)%10)==0:
            write(out/'rows.json',rows);print(f'Independent records shard {shard}: {len(hashes)}/{len(repetitions)} repetitions',flush=True)
    write(out/'rows.json',rows);write(out/'input_hashes.json',hashes)
    write(out/'pooled_residual_moments.json',dict(records=total,sum=sumr.tolist(),sum_outer=sumrr.tolist()))
    write(out/'completion.json',dict(status='complete',stage=stage,repetitions=len(hashes),rows=len(rows),
        simulated_independent_records=total,elapsed_seconds=time.perf_counter()-start,numerical_failures=0))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);p=a.parse_args();run(p.stage,p.output)
