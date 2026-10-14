"""Frozen independent calibration and small-noise guarded-method checks."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_support_wald import integrate
from scripts.run_wald_se_calibration import estimator,CORE
from scripts.run_efficiency_followup import input_digest
from src.support_information import POWERS,TARGET
from src.noise import ou_noise
from src.systems import trajectory
from src.weak_efficiency import nested_bumps
from src.polynomial_gls import moments,covariance
from src.guarded_gls import fit_guarded_polynomial
from src.self_instrumented_gls import fit as self_fit,dyadic_cubic_tests


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Research simulations execute on Kaggle')
    p=json.loads((ROOT/'configs/independent_validation_v2.json').read_text());parent=json.loads((ROOT/'configs'/p['parent_protocol_file']).read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':p,'source_sha256':source_hash()})
    base=np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0);smoke=stage.startswith('iq_smoke')
    if stage.startswith('iq_geo_'):
        geometry(stage,out,p,parent,base,start);return
    moderate=stage.startswith('iq_boot') or smoke
    level=p['moderate_level'] if moderate else p['small_level'];sigma=level*base
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1)
    phi,der,basis=nested_bumps(t,p['windows'],parent['basis_relative_cutoff']);phi0,der0,_=nested_bumps(t,parent['pilot_windows'])
    write(out/'basis.json',basis)
    if smoke:cells=[(0,0.,[p['smoke']['seed']])]
    elif stage.startswith('iq_nom_'):
        part=int(stage.rsplit('_',1)[1]);cells=[]
        cells.append((0,0.,[s for s in range(p['small_null_paths']) if s%p['nominal_shards']==part]))
        for index,multiplier in enumerate(p['small_power_multipliers'],1):
            cells.append((index,-multiplier*p['delta_reference_small'],[s for s in range(p['small_power_paths']) if s%p['nominal_shards']==part]))
    else:
        part=int(stage.rsplit('_',1)[1]);power=stage.startswith('iq_bootp_')
        n=p['moderate_power_paths'] if power else p['moderate_null_paths'];size=p['bootstrap_paths_per_shard']
        cells=[(1 if power else 0,-p['moderate_power_multiplier']*p['delta_reference_moderate'] if power else 0.,list(range(part*size,min((part+1)*size,n))))]
    rows=[];gates=[]
    for cell,delta,seeds in cells:
        truth=np.array([1.,-1.,1.,-1.,delta]);clean=integrate(t,truth,p['initial'])
        qtruth=np.zeros(18);qtruth[CORE]=truth
        A,b=moments(t,clean,phi,der,POWERS,np.zeros(2));gates.append({'cell':cell,'candidate_discrete_defect':float((b-A@qtruth)@ (b-A@qtruth))})
        write(out/'clean_checks.json',gates)
        for seed in seeds:
            timer=time.perf_counter();noise_seed=seed if smoke else (p['observation_seed_moderate'] if moderate else p['observation_seed_small'])+seed
            y=clean+ou_noise(len(t),p['h'],p['tau'],sigma,np.random.default_rng(noise_seed))
            q,V,_,_,pilot,diag=estimator(t,y,phi,der,phi0,der0,sigma,p['tau'],p['initial'],parent)
            nominal=float(np.sqrt(V[TARGET,TARGET]));critical=norm.ppf(1-p['alpha']/2)
            common={'seed':seed,'cell':cell,'level':level,'delta':delta,'input_sha256':input_digest(t,y),'candidate_coefficient':float(q[TARGET]),'coefficients':q.tolist()}
            rows.append({**common,'method':'nominal','se':nominal,'rejected':bool(abs(q[TARGET]/nominal)>critical),**diag})
            if moderate:
                nullpilot=pilot.copy();nullpilot[-1]=0.;mean=integrate(t,nullpilot,p['initial']);targets=[]
                count=p['smoke']['bootstrap_draws'] if smoke else p['bootstrap_draws']
                for draw in range(count):
                    bseed=p['smoke']['bootstrap_seed']+draw if smoke else p['bootstrap_seed_base']+p['bootstrap_seed_stride']*seed+p['bootstrap_cell_offset']*cell+draw
                    yy=mean+ou_noise(len(t),p['h'],p['tau'],sigma,np.random.default_rng(bseed))
                    fitted,*_=estimator(t,yy,phi,der,phi0,der0,sigma,p['tau'],p['initial'],parent)
                    targets.append(float(fitted[TARGET]))
                    if (draw+1)%100==0:print(f'{stage}: seed {seed}, bootstrap {draw+1}/{count}',flush=True)
                se=float(np.std(targets,ddof=1));assert np.isfinite(se) and se>0
                rows.append({**common,'method':'parametric_bootstrap','se':se,'rejected':bool(abs(q[TARGET]/se)>critical)})
                write(out/f'bootstrap_c{cell}_s{seed:04d}.json',{'targets':targets,'draws':count,'pilot':pilot.tolist(),'null_pilot':nullpilot.tolist()})
            write(out/'rows.json',rows);print(f'{stage}: cell {cell}, seed {seed}, seconds {time.perf_counter()-timer:.1f}',flush=True)
    if smoke:
        geometry(stage,out,p,parent,base,start,smoke=True)
    write(out/'completion.json',{'status':'complete','rows':len(rows),'elapsed_seconds':time.perf_counter()-start,'smoke':smoke})


def geometry(stage,out,p,parent,base,start,smoke=False):
    g=p['geometry'];index=1 if smoke else int(stage.rsplit('_',1)[1]);epsilon=g['levels'][index]
    refinement=max(0,int(np.ceil(.75*np.log2(.005/epsilon))));h=.005/2**refinement
    t=np.linspace(0,p['T'],round(p['T']/h)+1);truth=np.zeros(18);truth[CORE]=[1.,-1.,1.,-1.,0.]
    clean=integrate(t,truth[CORE],p['initial']);sigma=epsilon*base
    phi,der,basis=nested_bumps(t,g['windows'],parent['basis_relative_cutoff']);phi0,der0,_=nested_bumps(t,parent['pilot_windows'])
    self_phi,self_der=dyadic_cubic_tests(t,g['self_family_level'])
    A,b=moments(t,clean,phi,der,POWERS,np.zeros(2));first=covariance(t,clean,truth.reshape(2,9).T,phi,der,POWERS,base,p['tau'],first_order=True)
    check={'level':epsilon,'h':h,'basis':basis,'scaled_clean_design_singular_values':np.linalg.svd(A/np.linalg.norm(A,axis=0),compute_uv=False).tolist(),
        'unit_covariance_eigenvalues':np.linalg.eigvalsh(first).tolist(),'clean_residual_norm':float(np.linalg.norm(b-A@truth))}
    write(out/'geometry_clean_checks.json',check)
    bounds=np.tile(g['coefficient_box'],(18,1));records=[]
    n=p['smoke']['geometry_paths'] if smoke else g['paths_per_level']
    for seed in range(n):
        noise_seed=p['smoke']['seed']+1000000 if smoke else g['seed_base']+g['level_stride']*index+seed
        y=clean+ou_noise(len(t),h,p['tau'],sigma,np.random.default_rng(noise_seed))
        old,old_diag=fit_guarded_polynomial(t,y,phi,der,phi0,der0,POWERS,epsilon,base,p['tau'],CORE,parent['pilot_coefficient_box'],bounds,lambda t,q:integrate(t,q,p['initial']))
        new,new_diag=self_fit(t,y,epsilon,base,p['tau'],POWERS,self_phi,self_der,bounds)
        for method,estimate,diag in [('original_guarded',old,old_diag),('self_instrumented_guarded',new,new_diag)]:
            assert np.isfinite(estimate).all() and np.max(np.abs(estimate))<=3
            records.append({'seed':seed,'level':epsilon,'h':h,'method':method,'coefficients':estimate.tolist(),'input_sha256':input_digest(t,y),'diagnostics':diag})
        write(out/'geometry_rows.json',records);print(f'{stage}: geometry seed {seed+1}/{n}',flush=True)
    if not smoke:write(out/'completion.json',{'status':'complete','rows':len(records),'elapsed_seconds':time.perf_counter()-start,'smoke':False})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();run(a.stage,a.output)
