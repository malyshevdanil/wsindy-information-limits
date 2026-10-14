"""Known-initial endpoint checks on moderate-noise VdP and three-state Lorenz."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_efficiency_followup import generate,input_digest
from scripts.run_lorenz_efficiency import integrate,moments
from src.weak_efficiency import nested_bumps,vdp_moments,vdp_pilot_trajectory
from src.universal_boundary import universal_tests,fit as vdp_fit
from src.lorenz_boundary import fit as lorenz_fit
from src.boundary_weak import quadrature_weights
from src.path_scoring import ou_whiten
from src.noise import ou_noise


def write(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Research calculations run on Kaggle')
    protocol=json.loads((ROOT/'configs/boundary_generalization.json').read_text()); system=stage.rsplit('_',1)[-1]
    p=protocol[system];out.mkdir(parents=True,exist_ok=True);started=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':protocol,'system':system,'source_sha256':source_hash()})
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);initial=np.array(p['initial']);kwargs=dict(iterations=protocol['iterations'],tolerance=protocol['tolerance'],relaxation=protocol['covariance_relative_relaxation'])
    if system=='vdp5':
        t,x,base=generate(p);truth=np.array([1.,-1.,1.,-1.]);sigma=p['level']*base
        refs={r['seed']:r for r in json.loads((ROOT/'configs'/p['reference_file']).read_text())['rows']}
        fitting=vdp_fit
    else:
        truth=np.array(p['truth']);x,U=integrate(t,truth,initial,True);sigma=p['level']*np.std(x,axis=0);fitting=lorenz_fit
        information={}
        for unknown in [False,True]:
            V=U if unknown else U[:,:,:3];W=ou_whiten(V,p['h'],p['tau'],sigma).reshape(-1,V.shape[-1]);scale=np.linalg.norm(W,axis=0)
            C=np.linalg.inv((W/scale).T@(W/scale))/scale[:,None]/scale[None,:]
            information['unknown' if unknown else 'known']={'covariance':C[:3,:3].tolist(),'relative_rmse':float(np.sqrt(np.trace(C[:3,:3]))/np.linalg.norm(truth))}
        write(out/'information.json',information)
    phi0,der0,_=nested_bumps(t,p['pilot_windows']);phi,der,meta=universal_tests(t,p['widths'],protocol['basis_relative_cutoff'])
    free,fd,fmeta=nested_bumps(t,protocol['primary_interior_windows'],protocol['basis_relative_cutoff']);gates=[]
    geometry=[]
    for method,ph,de,controls in [('universal_boundary',phi,der,True),('matched_boundary_free',free,fd,False)]:
        extra={'mean_moments':True} if system=='vdp5' else {}
        q,diag=fitting(t,x,x,truth,ph,de,initial,sigma,p['tau'],controls=controls,**extra,**kwargs)
        defect=float(np.linalg.norm(q-truth)/np.linalg.norm(truth));gates.append({'method':method,'relative_bias':defect,**diag})
        write(out/'quadrature_gates.json',gates)
        if defect>protocol['clean_mean_quadrature_bias_max']:raise RuntimeError('Generalization quadrature gate failed')
        for first in (True,False):
            _,gd=fitting(t,x,x,truth,ph,de,initial,sigma,p['tau'],controls=controls,**extra,
                         first_order=first,return_variance=True,iterations=1,tolerance=protocol['tolerance'],relaxation=protocol['covariance_relative_relaxation'])
            V=np.array(gd['fixed_design_covariance'])
            geometry.append({'method':method,'first_order':first,'fixed_clean_design_covariance':V.tolist(),
                             'relative_rmse_surrogate':float(np.sqrt(np.trace(V))/np.linalg.norm(truth))})
    write(out/'oracle_geometry.json',geometry)
    write(out/'basis.json',{'universal':meta,'free':fmeta});rows=[]
    for seed in range(p['paths']):
        noise=ou_noise(len(t),p['h'],p['tau'],sigma,np.random.default_rng(p['seed_base']+seed))
        # Preserve the exact historical VdP scaling arithmetic for fingerprint checks.
        if system=='vdp5':noise=p['level']*ou_noise(len(t),p['h'],p['tau'],base,np.random.default_rng(p['seed_base']+seed))
        y=x+noise;digest=input_digest(t,y)
        if system=='vdp5':assert digest==refs[seed]['input_sha256']
        A0,b0=vdp_moments(t,y,phi0,der0,sigma) if system=='vdp5' else moments(t,y,phi0,der0)
        pilot=np.linalg.lstsq(A0,b0,rcond=None)[0]
        if system=='vdp5':model,q0=vdp_pilot_trajectory(t,pilot,initial,p['pilot_box'])
        else:q0=np.clip(pilot,np.array(p['pilot_box'])[:,0],np.array(p['pilot_box'])[:,1]);model=integrate(t,q0,initial)
        estimates={}
        for method,ph,de,controls in [('universal_boundary',phi,der,True),('matched_boundary_free',free,fd,False)]:
            timer=time.perf_counter();q,diag=fitting(t,y,model,q0,ph,de,initial,sigma,p['tau'],controls=controls,**kwargs)
            if not diag['converged']:raise RuntimeError('Generalization covariance iteration failed')
            estimates[method]=q
            rows.append({'seed':seed,'method':method,'coefficients':q.tolist(),'relative_error':float(np.linalg.norm(q-truth)/np.linalg.norm(truth)),
                         'input_sha256':digest,'fit_seconds':time.perf_counter()-timer,**diag})
        if system=='vdp5':
            for method in ('known_initial_mle','unknown_initial_mle','archived_boundary_free'):
                rows.append({**refs[seed][method],'seed':seed,'method':method,'input_sha256':digest,'archived_reference':True})
        else:
            cache={};timer=time.perf_counter()
            def evaluate(z):
                if 'q' not in cache or not np.array_equal(z,cache['q']):
                    path,V=integrate(t,z,initial,True)
                    cache.update(q=z.copy(),residual=ou_whiten(path-y,p['h'],p['tau'],sigma).ravel(),jacobian=ou_whiten(V[:,:,:3],p['h'],p['tau'],sigma).reshape(-1,3))
                return cache
            low,high=np.array(p['pilot_box']).T;z0=np.clip(estimates['matched_boundary_free'],low+1e-6,high-1e-6)
            fitted=least_squares(lambda z:evaluate(z)['residual'],z0,jac=lambda z:evaluate(z)['jacobian'],bounds=(low,high),
                                 xtol=p['likelihood_tolerance'],ftol=p['likelihood_tolerance'],gtol=p['likelihood_tolerance'],
                                 max_nfev=p['likelihood_max_nfev'],x_scale='jac')
            if not fitted.success:raise RuntimeError('Lorenz known-initial likelihood failed; no dropping')
            rows.append({'seed':seed,'method':'known_initial_mle','coefficients':fitted.x.tolist(),'relative_error':float(np.linalg.norm(fitted.x-truth)/np.linalg.norm(truth)),
                         'input_sha256':digest,'fit_seconds':time.perf_counter()-timer,'nfev':int(fitted.nfev),'active_bounds':fitted.active_mask.tolist()})
        write(out/'rows.json',rows)
        if seed%10==0:print(f'Endpoint generalization {system}: {seed+1}/{p["paths"]}',flush=True)
    write(out/'completion.json',{'status':'complete','paths':p['paths'],'elapsed_seconds':time.perf_counter()-started})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.stage,args.output)
