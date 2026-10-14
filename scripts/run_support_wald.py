"""Paired full-library weak-GLS Wald test and untuned author MSTLS refinement."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import argparse, json, sys, time, io, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.support_information import POWERS, TARGET
from src.polynomial_gls import moments, covariance, solve
from src.weak_efficiency import nested_bumps
from src.systems import trajectory
from src.noise import ou_noise
from scripts.run_efficiency_followup import input_digest
from scripts.run_experiments import source_hash
from scripts.run_full_adaptive import fetch, sha, run_octave


def write(path,data): path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def integrate(t,q,initial):
    a,b,c,d,delta=q
    def rhs(_,state):
        z,v=state;return [a*v,b*z+c*v+d*z*z*v+delta*z**3]
    sol=solve_ivp(rhs,(0.,float(t[-1])),initial,t_eval=t,method='DOP853',rtol=1e-10,atol=1e-12)
    if not sol.success: raise RuntimeError('Core-five pilot integration failed')
    return sol.y.T


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Research runs execute on Kaggle')
    filename='support_wald_v2.json' if stage.startswith('support_wald_v2_') else 'support_wald.json'
    p=json.loads((ROOT/'configs'/filename).read_text()); parent=json.loads((ROOT/'configs'/p['parent_protocol_file']).read_text())
    ref=json.loads((ROOT/'configs'/p['reference_file']).read_text()); shard=int(stage.rsplit('_',1)[-1]); assert 0<=shard<4
    li,part=divmod(shard,2);level=parent['levels'][li]; out.mkdir(parents=True,exist_ok=True); started=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':p,'parent_protocol':parent,'source_sha256':source_hash()})
    t=np.linspace(0,parent['T'],round(parent['T']/parent['h'])+1)
    base=np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0);sigma=level*base
    delta_ref=ref['delta_references'][str(level)]; deltas=[m*delta_ref for m in p['multipliers']]
    paths={}; phi,der,meta=nested_bumps(t,p['primary_windows'],p['basis_relative_cutoff'])
    phi0,der0,_=nested_bumps(t,p['pilot_windows']);core=[1,9,10,15,TARGET];gates=[]
    for di,delta in enumerate(deltas):
        truth=np.zeros(18);truth[core]=[1.,-1.,1.,-1.,delta]
        x=integrate(t,truth[core],parent['initial']); assert np.max(np.abs(x)) <= parent['maximum_clean_state_abs']
        paths[di]=x
        A,b=moments(t,x,phi,der,POWERS,np.zeros(2))
        uncorrected_b=b.copy()
        if p.get('pilot_quadrature_defect_correction'):b=b-(b-A@truth)
        C=covariance(t,x,truth.reshape(2,9).T,phi,der,POWERS,sigma,parent['tau'])
        estimate,V,diag=solve(A,b,C,p['covariance_relative_relaxation'])
        defect=float(abs(estimate[TARGET]-delta))
        gates.append({'delta_index':di,'delta':delta,'candidate_absolute_bias':defect,
                      'uncorrected_candidate_bias':float(abs(solve(A,uncorrected_b,C,p['covariance_relative_relaxation'])[0][TARGET]-delta)),
                      'all_coefficients_maximum_absolute_bias':float(np.max(np.abs(estimate-truth))),**diag})
        write(out/'quadrature_gates.json',gates)
        if defect > p['clean_mean_quadrature_coefficient_absolute_error_max']:raise RuntimeError('Candidate quadrature gate failed before Monte Carlo')
    write(out/'basis.json',meta)
    pin=json.loads((ROOT/'configs/author_pins.json').read_text())['WSINDy_ODE']
    data=fetch(f'https://codeload.github.com/{pin["repo"]}/zip/{pin["commit"]}');assert sha(data)==pin['zip_sha256']
    dest=ROOT/'author_sources/WSINDy_ODE';dest.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for member in z.infolist():
            f=Path(member.filename)
            if len(f.parts)==2 and f.suffix in ('.m','.md'):(dest/f.name).write_bytes(z.read(member))
    write(out/'author_provenance.json',pin)
    refs={(r['level'],r['seed'],r['multiplier']):r for r in ref['rows']};rows=[]
    low,high=np.array(p['pilot_coefficient_box']).T; critical=float(norm.ppf(1-p['alpha']/2))
    for seed in range(p['paths_per_level']):
        if seed%2 != part:continue
        noise=level*ou_noise(len(t),parent['h'],parent['tau'],base,np.random.default_rng(parent['seed_base']+li*parent['level_seed_stride']+seed))
        for di,(delta,multiplier) in enumerate(zip(deltas,p['multipliers'])):
            y=paths[di]+noise;digest=input_digest(t,y);key=(level,seed,multiplier)
            if key in refs: assert digest==refs[key]['input_sha256']
            timer=time.perf_counter();A0,b0=moments(t,y,phi0,der0,POWERS,sigma)
            pilot=np.linalg.lstsq(A0[:,core],b0,rcond=None)[0];bounded=np.clip(pilot,low,high)
            model=integrate(t,bounded,parent['initial']);coef=np.zeros(18);coef[core]=bounded
            A,b=moments(t,y,phi,der,POWERS,sigma)
            if p.get('pilot_quadrature_defect_correction'):
                model_A,model_b=moments(t,model,phi,der,POWERS,np.zeros(2))
                b=b-(model_b-model_A@coef)
            C=covariance(t,model,coef.reshape(2,9).T,phi,der,POWERS,sigma,parent['tau'])
            estimate,V,diag=solve(A,b,C,p['covariance_relative_relaxation']);se=float(np.sqrt(V[TARGET,TARGET])); zscore=float(estimate[TARGET]/se)
            common={'level':level,'seed':seed,'delta_index':di,'delta':delta,'multiplier':multiplier,'input_sha256':digest}
            rows.append({**common,'method':'weak_gls_wald','included':abs(zscore)>critical,'candidate_coefficient':float(estimate[TARGET]),
                         'standard_error':se,'z_score':zscore,'coefficients':estimate.tolist(),'fit_seconds':time.perf_counter()-timer,
                         'pilot_clipped':bool(np.any(pilot!=bounded)),**diag})
            if key in refs:
                rows.append({**common,'method':'author_mstls','included':refs[key]['included'],'candidate_coefficient':refs[key]['candidate_coefficient'],'archived_reference':True})
            else:
                timer=time.perf_counter();fitdir=out/f'fit_d{di}_s{seed:03d}';fitdir.mkdir(exist_ok=True)
                result=run_octave('wsindy',t,y,fitdir,{'WSINDy_ODE':dest});tags=np.asarray(result['tags']);match=np.flatnonzero(np.all(tags==[3,0],axis=1));assert len(match)==1
                candidate=float(np.asarray(result['coefficients'])[match[0],1])
                rows.append({**common,'method':'author_mstls','included':candidate!=0.,'candidate_coefficient':candidate,'fit_seconds':time.perf_counter()-timer})
            write(out/'rows.json',rows)
        if len(rows)%160==0:print(f'Wald/MSTLS shard {shard}: {len(rows)//16}/50 paths',flush=True)
    write(out/'completion.json',{'status':'complete','rows':len(rows),'elapsed_seconds':time.perf_counter()-started})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.stage,args.output)
