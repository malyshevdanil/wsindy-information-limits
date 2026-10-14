"""Kaggle-only technical and 50-repetition incidental-state likelihood pilot."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.profile_likelihood import fit_profile,ensemble_sensitivities,whiten_records
from src.path_scoring import vdp_sensitivities
from src.crossfit_gmm import pilot_from_means
from src.weak_efficiency import vdp_pilot_trajectory
from src.weak_features import test_functions,quadrature_weights
from scripts.run_heterogeneous_records import trajectories,noise
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical simulations run on Kaggle')
    config='incidental_pilot_v2.json' if '_v2_' in stage else 'incidental_pilot.json'
    p=json.loads((ROOT/'configs'/config).read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    smoke=stage=='incidental_smoke';shard=0 if smoke else int(stage.rsplit('_',1)[1]);reps=[0] if smoke else list(range(shard,p['repetitions'],p['shards']))
    ns=[8,p['records']] if smoke else [p['records']];truth=np.array(p['truth']);rows=[];checks=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);baseline,_=vdp_pilot_trajectory(t,truth,p['initial_center'],[[-3,3]]*4);base=baseline.std(axis=0)
    phi,der=test_functions(t,[2.1,5.2,8.3,11.4,14.5,17.6],4,4);w=quadrature_weights(t)
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash(),repetitions=reps,record_counts=ns,smoke=smoke))
    for rep in reps:
        rng=np.random.default_rng(p['pilot_seed_base']+rep);u=np.array(p['initial_center'])+rng.uniform(-p['initial_half_width'],p['initial_half_width'],(p['records'],2))
        x=trajectories(t,truth,u);eta=noise(rng,x.shape,base,np.exp(-p['h']/p['tau']));xi=noise(rng,x.shape,base,np.exp(-p['h']/p['tau']));y=x+p['level']*eta
        digest=hashlib.sha256()
        for array in (u,x,eta,xi):digest.update(np.ascontiguousarray(array).tobytes())
        for n in ns:
            from scripts.run_heterogeneous_records import moments
            A,b=moments(y[:n],phi*w,der*w,(p['level']*base)**2)
            q0=pilot_from_means(A.mean(axis=0),b.mean(axis=0),np.column_stack((p['coefficient_lower'],p['coefficient_upper'])))
            for known in (False,True):
                q,uh,V,diagnostic=fit_profile(t,y[:n],q0,u[:n] if known else y[:n,0],p['tau'],p['level']*base,p,known_initial=known)
                error=q-truth;stat=float(error@np.linalg.solve(V,error));sd=np.sqrt(np.diag(V))
                rows.append(dict(repetition=rep,n=n,level=p['level'],method='oracle_likelihood' if known else 'profile_likelihood',coefficients=q.tolist(),
                    initial_estimates=uh.tolist(),covariance=V.tolist(),relative_error=float(np.linalg.norm(error)/2),relative_information_sd=float(np.sqrt(np.trace(V))/2),
                    wald_joint_covered=bool(stat<=chi2.ppf(.95,4)),wald_marginal_covered=(np.abs(error/sd)<=1.959963984540054).tolist(),
                    input_sha256=digest.hexdigest(),diagnostic=diagnostic))
                write(out/'rows.json',rows)
                if not diagnostic['success']:raise RuntimeError('Likelihood convergence check failed; no repetition excluded')
                if smoke and n==8 and not known:
                    def evaluate(z):
                        model,U=ensemble_sensitivities(t,z[:4],z[4:].reshape(n,2),p['ode_rtol'],p['ode_atol'])
                        r=whiten_records(model-y[:n],p['h'],p['tau'],p['level']*base).ravel();W=whiten_records(U,p['h'],p['tau'],p['level']*base).reshape(n,-1,6)
                        J=np.zeros((len(r),4+2*n));size=W.shape[1]
                        for k in range(n):J[k*size:(k+1)*size,:4]=W[k,:,:4];J[k*size:(k+1)*size,4+2*k:6+2*k]=W[k,:,4:]
                        return r,J
                    z=np.r_[q0,y[:n,0].ravel()]
                    fit=least_squares(lambda z:evaluate(z)[0],z,jac=lambda z:evaluate(z)[1],
                        bounds=(np.r_[p['coefficient_lower'],np.tile(p['initial_lower'],n)],np.r_[p['coefficient_upper'],np.tile(p['initial_upper'],n)]),
                        ftol=1e-12,xtol=1e-10,gtol=1e-8,max_nfev=100,x_scale='jac')
                    discrepancy=float(np.max(np.abs(fit.x[:4]-q)));checks.append(dict(check='dense_scipy_joint_fit',success=bool(fit.success),maximum_coefficient_difference=discrepancy))
                    if not fit.success or discrepancy>1e-5:raise RuntimeError('Dense optimizer validation failed')
        print(f'Incidental likelihood pilot shard {shard}: {reps.index(rep)+1}/{len(reps)} repetitions; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    write(out/'technical_checks.json',checks);write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start,repetitions=len(reps),rows=len(rows),discarded_repetitions=0))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
