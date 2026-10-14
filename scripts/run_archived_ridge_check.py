"""Deterministic ridge sensitivity on the actual archived test spans."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_support_wald import integrate as vdp_integrate
from scripts.run_lorenz_efficiency import integrate as lorenz_integrate
from src.weak_efficiency import nested_bumps,vdp_moments,covariance_blocks,assemble_covariance,gls_map
from src.systems import trajectory
import src.universal_boundary as universal
import src.lorenz_boundary as lorenz

def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False),encoding='utf-8')

def check(A,C,name,rho):
    P0,_=gls_map(A,C,0.);Pr,level=gls_map(A,C,rho)
    V0=P0@C@P0.T;Vr=Pr@C@Pr.T;delta=Pr-P0;increment=delta@C@delta.T
    relative=float(np.trace(increment)/np.trace(V0));assert relative>=0
    return {'family':name,'moments':len(C),'parameters':A.shape[1],'V0':V0.tolist(),'Vrho':Vr.tolist(),
        'psd_increment':increment.tolist(),'relative_trace_increment':relative,
        'relative_rmse_increment':float(relative/(np.sqrt(1+relative)+1)),
        'relative_map_change':float(np.linalg.norm(delta)/np.linalg.norm(P0)),
        'normalized_identity_residual':float(np.linalg.norm(Vr-V0-increment)/np.linalg.norm(V0)),
        'ridge':float(level),'normalization_error':float(np.linalg.norm(Pr@A-np.eye(A.shape[1])))}

def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('ODE work runs on Kaggle')
    p=json.loads((ROOT/'configs/archived_ridge_check.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':p,'source_sha256':source_hash()});rows=[];captured=[]
    def spy(A,C,relaxation=1e-10):
        captured[:]=[(A.copy(),C.copy())];return gls_map(A,C,relaxation)
    universal.gls_map=spy;lorenz.gls_map=spy
    for system in ['vdp','lorenz']:
        q=p[system];t=np.linspace(0,q['T'],round(q['T']/q['h'])+1)
        truth=np.array([1.,-1.,1.,-1.]) if system=='vdp' else np.array(q['truth'])
        x=vdp_integrate(t,np.r_[truth,0.],q['initial']) if system=='vdp' else lorenz_integrate(t,truth,np.array(q['initial']))
        sigma=np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0) if system=='vdp' else np.std(x,axis=0)
        free,fd,_=nested_bumps(t,480,p['basis_relative_cutoff'])
        if system=='vdp':
            A,_=vdp_moments(t,x,free,fd,np.zeros(2));B,H=covariance_blocks(t,x,free,fd,sigma,q['tau']);C=assemble_covariance(truth,B,H*0,sigma)
            rows.append(check(A,C,'vdp_compact_trapezoid_480',p['relative_ridge']))
            widths=[3.,1.5,.75,.375,.1875,.09375]
            families=[('compact_simpson',free,fd,False),('initial',*universal.universal_tests(t,widths,p['basis_relative_cutoff'],('initial',))[:2],False),
                ('terminal',*universal.universal_tests(t,widths,p['basis_relative_cutoff'],('terminal',))[:2],False),
                ('initial_control',free,fd,True),('full',*universal.universal_tests(t,widths,p['basis_relative_cutoff'])[:2],True)]
        else:
            widths=[.6,.3,.15,.075,.0375,.01875]
            families=[('compact_simpson',free,fd,False),('full',*universal.universal_tests(t,widths,p['basis_relative_cutoff'])[:2],True)]
        for name,phi,der,control in families:
            extra={'mean_moments':True} if system=='vdp' else {}
            fitting=universal.fit if system=='vdp' else lorenz.fit
            fitting(t,x,x,truth,phi,der,np.array(q['initial']),sigma,q['tau'],controls=control,first_order=True,iterations=1,**extra)
            A,C=captured[0];rows.append(check(A,C,system+'_'+name,p['relative_ridge']))
            write(out/'ridge_check.json',{'rows':rows});print('Completed ridge geometry '+system+' '+name,flush=True)
    write(out/'completion.json',{'status':'complete','families':len(rows),'elapsed_seconds':time.perf_counter()-start})

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);v=a.parse_args();run(v.stage,v.output)
