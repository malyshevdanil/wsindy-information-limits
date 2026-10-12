"""Prespecified calibration on original independent-record seeds, Kaggle only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.record_gmm import fit_mean
from src.record_bootstrap import batch_fit,bootstrap_t
from src.weak_features import test_functions,quadrature_weights
from src.weak_efficiency import vdp_pilot_trajectory
from scripts.run_independent_records import summaries
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical runs must execute on Kaggle')
    p=json.loads((ROOT/'configs/independent_records.json').read_text());c=json.loads((ROOT/'configs/record_bootstrap.json').read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();shard=int(stage.rsplit('_',1)[1])
    write(out/'manifest.json',dict(stage=stage,protocol=c,source_protocol=p,source_sha256=source_hash()))
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);truth=np.array(p['truth'])
    x,_=vdp_pilot_trajectory(t,truth,p['initial'],[[-3,3]]*4)
    sigma=p['noise_level']*np.std(x,axis=0);rho=np.exp(-p['h']/p['tau'])
    phi,der=test_functions(t,p['test_centers'],p['test_width'],p['test_order']);w=quadrature_weights(t);P=phi*w;D=der*w
    rows=[];hashes=[];max_batch_discrepancy=0.;reps=list(range(shard,c['monte_carlo_repetitions'],c['shards']))
    for rep in reps:
        rng=np.random.default_rng(p['seed_base']+rep);R=max(p['record_counts']);N=len(t)
        innovations=rng.normal(size=(R,N,2))*sigma*np.sqrt(1-rho*rho);innovations[:,0]/=np.sqrt(1-rho*rho)
        y=x+lfilter([1],[1,-rho],innovations,axis=1)
        raw,b,vs,lag,previous=summaries(y,P,D);digest=hashlib.sha256()
        for a in (raw,b,vs,lag,previous):digest.update(np.ascontiguousarray(a).tobytes())
        hashes.append(dict(repetition=rep,sufficient_statistics_sha256=digest.hexdigest()))
        rawpair=(raw[0::2]+raw[1::2])/2;bpair=(b[0::2]+b[1::2])/2
        for n in c['record_counts']:
            B=n//2;vhat=vs[:B].mean(axis=0)
            q,V,_=fit_mean(rawpair[:B],bpair[:B],vs[:B],vhat,True)
            qb,Vb=batch_fit(rawpair[None,:B],bpair[None,:B],vs[None,:B])
            discrepancy=max(float(np.max(np.abs(qb[0]-q))),float(np.max(np.abs(Vb[0]-V))))
            max_batch_discrepancy=max(max_batch_discrepancy,discrepancy)
            if discrepancy>1e-8:raise RuntimeError('Batched estimator differs from original')
            seed=c['bootstrap_seed_base']+rep*1000+n
            critical=bootstrap_t(rawpair[:B],bpair[:B],vs[:B],q,V,np.random.default_rng(seed),c['bootstrap_draws'],c['batch_size'])
            sd=np.sqrt(np.diag(V));limits=np.array(critical['marginal_t_quantiles'])
            lower=q-limits[1]*sd;upper=q-limits[0]*sd;error=q-truth
            statistic=float(error@np.linalg.solve(V,error))
            rows.append(dict(repetition=rep,n=n,coefficients=q.tolist(),covariance=V.tolist(),bootstrap_seed=seed,
                marginal_lower=lower.tolist(),marginal_upper=upper.tolist(),
                marginal_covered=((truth>=lower)&(truth<=upper)).tolist(),
                joint_covered=bool(statistic<=critical['joint_critical_value']),joint_statistic=statistic,**critical))
        if len(hashes)%5==0:
            write(out/'rows.json',rows);print(f'Record bootstrap shard {shard}: {len(hashes)}/{len(reps)} repetitions',flush=True)
    write(out/'rows.json',rows);write(out/'input_hashes.json',hashes)
    write(out/'completion.json',dict(status='complete',repetitions=len(reps),rows=len(rows),
        bootstrap_fits=len(rows)*c['bootstrap_draws'],maximum_batch_fit_discrepancy=max_batch_discrepancy,
        failed_or_discarded_draws=0,elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);p=a.parse_args();run(p.stage,p.output)
