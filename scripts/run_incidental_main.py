"""Matched likelihood components of experiment A, including full refit bootstrap.

Physical inputs are reusable by the independent cross-fitted weak component.
No weak estimate is replaced by a likelihood estimate while its cache is pending.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.profile_likelihood import fit_profile,ensemble_sensitivities,whiten_records,normal_blocks,block_step
from src.crossfit_gmm import pilot_from_means
from src.weak_efficiency import vdp_pilot_trajectory
from src.weak_features import test_functions,quadrature_weights
from scripts.run_heterogeneous_records import trajectories,noise,moments
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Simulations run on Kaggle')
    config='incidental_compact.json' if stage.startswith('inc_short_') else 'incidental_trimmed.json' if stage.startswith('inc_trim_') else 'incidental_main.json'
    p=json.loads((ROOT/'configs'/config).read_text());out.mkdir(parents=True,exist_ok=True)
    smoke=stage=='inc_main_smoke'
    if smoke:li=0;n=8;reps=[0];draws=2
    else:
        # Separate n/amplitude blocks keep 50-refit jobs below notebook limits.
        _,_,ls,ns,bs=stage.split('_');li=int(ls[1:]);n=int(ns[1:]);block=int(bs[1:])
        first=block*p['repetitions_per_job'];reps=list(range(first,min(first+p['repetitions_per_job'],p['monte_carlo_repetitions'])))
        assert n in p['record_counts'] and reps
        draws=p['parametric_bootstrap_draws']
        if config=='incidental_trimmed.json' and not any(c['level']==p['levels'][li] and c['n']==n for c in p['bootstrap_correction_cells']):draws=0
    level=p['levels'][li];truth=np.array(p['truth']);start=time.perf_counter();rows=[];bootrows=[];hashes=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);rho=np.exp(-p['h']/p['tau'])
    baseline,_=vdp_pilot_trajectory(t,truth,p['initial_center'],[[-3,3]]*4);base=baseline.std(axis=0);sigma=level*base
    phi,der=test_functions(t,[2.1,5.2,8.3,11.4,14.5,17.6],4,4);w=quadrature_weights(t)
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash(),smoke=smoke,repetitions=reps,n=n,level=level,bootstrap_draws=draws))
    for rep in reps:
        rng=np.random.default_rng(p['seed_base']+rep)
        u=np.array(p['initial_center'])+rng.uniform(-p['initial_half_width'],p['initial_half_width'],(800,2))
        x=trajectories(t,truth,u);eta=noise(rng,x.shape,base,rho);xi=noise(rng,x.shape,base,rho)
        digest=hashlib.sha256()
        for array in (u,x,eta,xi):digest.update(np.ascontiguousarray(array).tobytes())
        hashes.append(dict(repetition=rep,input_sha256=digest.hexdigest()))
        y=(x+level*eta)[:n];ur=u[:n]
        weak_start_timer=time.perf_counter()
        A,b=moments(y,phi*w,der*w,sigma**2)
        q0=pilot_from_means(A.mean(axis=0),b.mean(axis=0),np.column_stack((p['coefficient_lower'],p['coefficient_upper'])))
        weak_start_seconds=time.perf_counter()-weak_start_timer
        _,Ut=ensemble_sensitivities(t,truth,ur,p['ode_rtol'],p['ode_atol'])
        Wt=whiten_records(Ut,p['h'],p['tau'],sigma)
        It,gt=normal_blocks(Wt,np.zeros_like(y));_,_,Ie,_=block_step(It,gt)
        local_sd=float(np.sqrt(np.trace(np.linalg.inv(Ie)))/2)
        def save(method,q,V,diagnostic,extra=None):
            e=q-truth;stat=float(e@np.linalg.solve(V,e));sd=np.sqrt(np.diag(V))
            row=dict(repetition=rep,n=n,level=level,method=method,coefficients=q.tolist(),covariance=V.tolist(),
                relative_error=float(np.linalg.norm(e)/2),joint_statistic=stat,
                wald_joint_covered=bool(stat<=chi2.ppf(.95,4)),wald_marginal_covered=(np.abs(e/sd)<=1.959963984540054).tolist(),
                local_information_relative_sd=local_sd,input_sha256=digest.hexdigest(),diagnostic=diagnostic,
                shared_weak_start_seconds=weak_start_seconds,
                runtime_excludes='simulated truth/data generation and information-reference diagnostics')
            if extra:row.update(extra)
            rows.append(row);write(out/'rows.json',rows)
            if not diagnostic['success']:raise RuntimeError('Required fit failed; no repetition discarded')
        q,uh,V,diagnostic=fit_profile(t,y,q0,y[:,0],p['tau'],sigma,p)
        save('profile_likelihood',q,V,diagnostic,dict(initial_estimates=uh.tolist()))
        qo,uo,Vo,do=fit_profile(t,y,q0,ur,p['tau'],sigma,p,known_initial=True)
        save('oracle_likelihood',qo,Vo,do)
        write(out/'input_hashes.json',hashes)
        if not draws:
            write(out/'bootstrap_fits.json',bootrows)
            print(f'Incidental main {stage}: {reps.index(rep)+1}/{len(reps)}; no correction scheduled',flush=True)
            continue
        bootstrap_start=time.perf_counter();mean=trajectories(t,q,uh);estimates=[]
        # The decimal blocks separate amplitude, prefix n, repetition and draw.
        boot_seed=p['parametric_bootstrap_seed_base']+li*10000000+n*10000+rep*100
        for draw in range(draws):
            brng=np.random.default_rng(boot_seed+draw)
            yb=mean+noise(brng,mean.shape,sigma,rho)
            qb,ub,Vb,db=fit_profile(t,yb,q,uh,p['tau'],sigma,p)
            bootrows.append(dict(repetition=rep,n=n,level=level,draw=draw,seed=boot_seed+draw,
                coefficients=qb.tolist(),diagnostic=db))
            write(out/'bootstrap_fits.json',bootrows)
            if not db['success']:raise RuntimeError('Bootstrap refit failed; no draw discarded')
            estimates.append(qb)
        estimates=np.array(estimates);bias=estimates.mean(axis=0)-q;corrected=q-bias
        elapsed=time.perf_counter()-bootstrap_start
        dc=dict(success=True,elapsed_seconds=diagnostic['elapsed_seconds']+elapsed,bootstrap_refits=draws,
                bound_active_coefficients=diagnostic['bound_active_coefficients'],
                bound_active_initial_coordinates=diagnostic['bound_active_initial_coordinates'])
        save('bias_corrected_likelihood',corrected,V,dc,dict(bootstrap_bias=bias.tolist(),
            bootstrap_bias_mc_standard_error=(estimates.std(axis=0,ddof=1)/np.sqrt(draws)).tolist(),
            covariance_convention='uncorrected profile Fisher Wald covariance; finite-bootstrap correction uncertainty assessed by Monte Carlo coverage'))
        write(out/'input_hashes.json',hashes)
        print(f'Incidental main {stage}: {reps.index(rep)+1}/{len(reps)}; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start,
        repetitions=len(reps),rows=len(rows),bootstrap_refits=len(bootrows),discarded_records_or_fits=0,
        weak_component='separate matching-input crossfit job required'))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True)
    args=a.parse_args();run(args.stage,args.output)
