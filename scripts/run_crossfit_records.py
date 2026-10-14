"""Kaggle-only matched heterogeneous records with five-fold pilot refitting."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_heterogeneous_records import trajectories,noise,moments
from src.crossfit_gmm import crossfit_equations,crossfit_bootstrap,pilot_from_means
from src.crossfit_cache import features
from src.heterogeneous_gmm import interpolate_grid,record_fit,record_bootstrap
from src.weak_efficiency import nested_bumps
from src.weak_features import test_functions,quadrature_weights
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def evaluate(Bcache,acache,pilots,anchor,folds,exponents=None,radius=None):
    affine=(np.column_stack((np.ones(len(folds)),pilots[folds]-anchor)) if exponents is None
            else features(pilots[folds],anchor,radius,exponents))
    return np.einsum('rp,prij->rij',affine,Bcache),np.einsum('rp,pri->ri',affine,acache)


def summarize(q,V,truth):
    error=q-truth;sd=np.sqrt(np.diag(V));stat=float(error@np.linalg.solve(V,error))
    if not np.isfinite(q).all() or np.linalg.eigvalsh(V).min()<=0:raise RuntimeError('Invalid fit; no results discarded')
    return dict(coefficients=q.tolist(),covariance=V.tolist(),relative_error=float(np.linalg.norm(error)/2),
        joint_statistic=stat,wald_joint_covered=bool(stat<=chi2.ppf(.95,4)),
        wald_marginal_covered=(np.abs(error/sd)<=1.959963984540054).tolist())


def add_intervals(row,q,V,truth,limits,critical,prefix):
    sd=np.sqrt(np.diag(V));lower=q-limits[1]*sd;upper=q-limits[0]*sd
    row.update({prefix+'_marginal_lower':lower.tolist(),prefix+'_marginal_upper':upper.tolist(),
        prefix+'_t_quantiles':limits.tolist(),prefix+'_joint_critical':float(critical),
        prefix+'_marginal_covered':((truth>=lower)&(truth<=upper)).tolist(),
        prefix+'_joint_covered':bool(row['joint_statistic']<=critical)})


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical experiments run on Kaggle')
    small=stage.startswith(('crossfit120_','repair120_'))
    main_a=stage.startswith(('crossfitA_','repairA'))
    repair=stage.split('_')[0] if stage.startswith('repair') else None
    config=repair+'_records.json' if repair else 'incidental_weak.json' if main_a else 'crossfit120_records.json' if small else 'crossfit_records.json'
    prefix=repair if repair else 'crossfitA' if main_a else 'crossfit120' if small else 'crossfit'
    p=json.loads((ROOT/'configs'/config).read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    smoke=stage.endswith('_smoke');shard=0 if smoke else int(stage.rsplit('_',1)[1])
    reps=[0,1] if smoke else list(range(shard,p['monte_carlo_repetitions'],p['shards']))
    counts=p.get('smoke_record_counts',[50,800,1050]) if smoke else p['record_counts'];truth=np.array(p['truth'])
    cache_path=ROOT/'configs'/(prefix+'_instruments.npz')
    with np.load(cache_path) as z:
        anchor=z['anchor'];base=z['base'];libraries={level:z['cache_'+str(level)] for level in p['levels']}
        exponents=z['exponents'] if 'exponents' in z else None;radius=float(z['radius']) if 'radius' in z else None
        oracles={key:z[key] for key in z.files if key.startswith('oracle_')}
    geometry=json.loads((ROOT/'configs'/(prefix+'_cache_geometry.json')).read_text());assert geometry['cache_validated']
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash(),
        instrument_cache_sha256=hashlib.sha256(cache_path.read_bytes()).hexdigest(),smoke=smoke,
        repetitions=reps,record_counts=counts,cache_geometry=geometry))
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);rho=np.exp(-p['h']/p['tau']);w=quadrature_weights(t)
    phi,der,meta=nested_bumps(t,p['windows'],p['basis_relative_cutoff']);P=phi*w;D=der*w
    pilot_phi,pilot_der=test_functions(t,[2.1,5.2,8.3,11.4,14.5,17.6],4,4)
    xs=p['instrument_grid']['x'];ys=p['instrument_grid']['v'];rows=[];hashes=[];pilot_rows=[]
    ref=json.loads((ROOT/'configs/crossfit_information.json').read_text())
    for rep in reps:
        rng=np.random.default_rng(p['seed_base']+rep)
        u=np.array(p['initial_center'])+rng.uniform(-p['initial_half_width'],p['initial_half_width'],(800,2))
        x=trajectories(t,truth,u);eta=noise(rng,x.shape,base,rho);xi=noise(rng,x.shape,base,rho)
        old_hash=hashlib.sha256()
        for array in (u,x,eta,xi):old_hash.update(np.ascontiguousarray(array).tobytes())
        if p.get('include_split_diagnostic',True):
            extra=np.random.default_rng(p['extra_record_seed_base']+rep)
            ue=np.array(p['initial_center'])+extra.uniform(-p['initial_half_width'],p['initial_half_width'],(250,2))
            xe=trajectories(t,truth,ue);ee=noise(extra,xe.shape,base,rho);ze=noise(extra,xe.shape,base,rho)
            u=np.concatenate((u,ue));x=np.concatenate((x,xe));eta=np.concatenate((eta,ee));xi=np.concatenate((xi,ze))
        digest=hashlib.sha256()
        for array in (u,x,eta,xi):digest.update(np.ascontiguousarray(array).tobytes())
        hashrow=dict(repetition=rep,first_800_input_sha256=old_hash.hexdigest(),all_input_sha256=digest.hexdigest(),physical_records=len(u))
        if len(u)==1050:hashrow['all_1050_input_sha256']=digest.hexdigest()
        hashes.append(hashrow)
        Ir=interpolate_grid(np.array(ref['information_grid']),ref['information_x'],ref['information_v'],u)
        for level in p['levels']:
            preprocessing_start=time.perf_counter()
            y=x+level*eta;ya=y+level*xi;yb=y-level*xi
            Ap,bp=moments(y,pilot_phi*w,pilot_der*w,(level*base)**2)
            Aa,ba=moments(ya,P,D,2*(level*base)**2);Ab,bb=moments(yb,P,D,2*(level*base)**2)
            Bcache=[];acache=[]
            for cache in libraries[level]:
                La=interpolate_grid(cache,xs,ys,ya[:,0]);Lb=interpolate_grid(cache,xs,ys,yb[:,0])
                Bcache.append(.5*(np.einsum('rpm,rmq->rpq',La,Ab)+np.einsum('rpm,rmq->rpq',Lb,Aa)))
                acache.append(.5*(np.einsum('rpm,rm->rp',La,bb)+np.einsum('rpm,rm->rp',Lb,ba)))
            Bcache=np.array(Bcache);acache=np.array(acache)
            preprocessing_seconds=time.perf_counter()-preprocessing_start
            for n in counts:
                fit_start=time.perf_counter()
                folds=np.arange(n)%p['folds'];q,V,pilots=crossfit_equations(Ap[:n],bp[:n],Bcache[:,:n],acache[:,:n],anchor,p['pilot_box'],folds,cache_exponents=exponents,cache_radius=radius)
                q=q[0];V=V[0];pilots=pilots[0];B,a=evaluate(Bcache[:,:n],acache[:,:n],pilots,anchor,folds,exponents,radius)
                qe,Ve=record_fit(B,a)
                if not (np.allclose(q,qe,rtol=1e-10,atol=1e-11) and np.allclose(V,Ve,rtol=1e-9,atol=1e-15)):
                    raise RuntimeError('Expanded crossfit equations disagree')
                row=dict(repetition=rep,n=n,level=level,method='crossfit'+str(p['windows']),inference_records=n,training_records_extra=0,**summarize(q,V,truth))
                row['local_information_relative_sd']=float(np.sqrt(np.trace(np.linalg.inv(Ir[:n].sum(axis=0)*(.05/level)**2)))/2)
                row['estimate_seconds']=time.perf_counter()-fit_start
                row['shared_preprocessing_seconds']=preprocessing_seconds
                row['shared_preprocessing_records']=len(u)
                seed=p['bootstrap_seed_base']+rep*10000+n+round(level*1000)
                if 'bootstrap_amplitude_index' in p:
                    seed=p['bootstrap_seed_base']+rep*1000000+p['bootstrap_amplitude_index']*100000+n
                boot_start=time.perf_counter()
                limits,critical=crossfit_bootstrap(Ap[:n],bp[:n],Bcache[:,:n],acache[:,:n],anchor,p['pilot_box'],folds,q,np.random.default_rng(seed),p['bootstrap_draws'],exponents,radius)
                row['bootstrap_refit_seconds']=time.perf_counter()-boot_start
                add_intervals(row,q,V,truth,limits,critical,'refit')
                if p.get('compare_frozen_bootstrap',True):
                    boot_start=time.perf_counter()
                    limits,critical=crossfit_bootstrap(Ap[:n],bp[:n],Bcache[:,:n],acache[:,:n],anchor,p['pilot_box'],folds,q,np.random.default_rng(seed),p['bootstrap_draws'],exponents,radius,refit_pilots=False)
                    add_intervals(row,q,V,truth,limits,critical,'frozen')
                    row['bootstrap_frozen_seconds']=time.perf_counter()-boot_start
                row['bootstrap_seed']=seed;rows.append(row)
                pilot_rows.append(dict(repetition=rep,n=n,level=level,coefficients=pilots.tolist(),
                    maximum_relative_error=float(np.max(np.linalg.norm(pilots-truth,axis=1)/2)),
                    clipped_coordinates=int(np.sum((pilots==np.array(p['pilot_box'])[:,0])|(pilots==np.array(p['pilot_box'])[:,1]))),
                    cache_clipped_coordinates=int(np.sum(np.abs(pilots-anchor)>radius)) if radius is not None else 0))
            if oracles:
                Ao,bo=moments(y,P,D,(level*base)**2)
                for label,Ao_,bo_ in (('original',Ao,bo),('fission',.5*(Aa+Ab),.5*(ba+bb))):
                    L=interpolate_grid(oracles['oracle_'+label+'_'+str(level)],xs,ys,u)
                    Bo=np.einsum('rpm,rmq->rpq',L,Ao_);ao=np.einsum('rpm,rm->rp',L,bo_)
                    for n in counts:
                        qo,Vo=record_fit(Bo[:n],ao[:n])
                        row=dict(repetition=rep,n=n,level=level,method='oracle120_'+label,
                            inference_records=n,training_records_extra=0,**summarize(qo,Vo,truth))
                        row['local_information_relative_sd']=float(np.sqrt(np.trace(np.linalg.inv(Ir[:n].sum(axis=0)*(.05/level)**2)))/2)
                        rows.append(row)
            if p.get('include_split_diagnostic',True):
                pilot=pilot_from_means(Ap[50:].mean(axis=0),bp[50:].mean(axis=0),p['pilot_box'])
                B,a=evaluate(Bcache[:,:50],acache[:,:50],pilot[None],anchor,np.zeros(50,dtype=int),exponents,radius)
                q,V=record_fit(B,a)
                rows.append(dict(repetition=rep,n=1050,level=level,method='split1000_50',inference_records=50,training_records_extra=1000,
                    local_information_relative_sd=float(np.sqrt(np.trace(np.linalg.inv(Ir.sum(axis=0)*(.05/level)**2)))/2),**summarize(q,V,truth)))
            if rep==reps[0] and level==p['levels'][0]:
                fixture=dict(A=Ap,b=bp,Bcache=Bcache,acache=acache,anchor=anchor)
                if exponents is not None:fixture.update(exponents=exponents,radius=radius)
                np.savez_compressed(out/'crossfit_fixture.npz',**fixture)
        write(out/'rows.json',rows);write(out/'pilot_rows.json',pilot_rows);write(out/'input_hashes.json',hashes)
        print(f'Crossfit records shard {shard}: {reps.index(rep)+1}/{len(reps)} repetitions; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    write(out/'completion.json',dict(status='complete',repetitions=len(reps),rows=len(rows),elapsed_seconds=time.perf_counter()-start,
        discarded_records_or_bootstrap_fits=0,bootstrap_pilots_refitted=True))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
