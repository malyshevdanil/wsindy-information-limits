"""Deterministic unregularized covariance reference for the S28 test family."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,cho_factor,cho_solve
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_support_wald import integrate
from scripts.run_wald_se_calibration import CORE
from src.support_information import POWERS
from src.systems import trajectory
from src.weak_features import library,derivative_polynomial
from src.self_instrumented_gls import dyadic_cubic_tests
from src.boundary_weak import quadrature_weights
from src.polynomial_gls import moments,covariance

def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')

def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('ODE integrations execute on Kaggle')
    p=json.loads((ROOT/'configs/guarded_geometry_reference.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',{'stage':stage,'protocol':p,'source_sha256':source_hash()})
    truth=np.array(p['coefficients'],dtype=float);base=np.std(trajectory('vdp',np.linspace(0,20,2001)),axis=0);rows=[]
    for h in p['reference_steps']:
        t=np.linspace(0,p['T'],round(p['T']/h)+1);clean=integrate(t,truth[CORE],p['initial'])
        rhs=np.column_stack([clean[:,1],-clean[:,0]+clean[:,1]-clean[:,0]**2*clean[:,1]])
        z=(t-10)/9.5;v=np.maximum(1-z*z,0);w=v**4;wp=np.where(np.abs(z)<1,-8*z*v**3/9.5,0.)
        g=library(clean,POWERS,np.zeros(2));gd=np.zeros_like(g)
        for j in range(2):
            alpha=[0,0];alpha[j]=1
            gd+=derivative_polynomial(clean,POWERS,np.eye(9),tuple(alpha))*rhs[:,j,None]
        cubic,cd=dyadic_cubic_tests(t,p['spline_level']);phi=np.concatenate([(w[:,None]*g).T,cubic]);der=np.concatenate([(wp[:,None]*g+w[:,None]*gd).T,cd])
        q=quadrature_weights(t);gram=(phi*q)@phi.T+(der*q)@der.T;values,U=eigh((gram+gram.T)/2)
        assert values[0]>0;transform=(U/np.sqrt(values)).T;phi=transform@phi;der=transform@der
        D,b=moments(t,clean,phi,der,POWERS,np.zeros(2));C=covariance(t,clean,truth.reshape(2,9).T,phi,der,POWERS,base,p['tau'],first_order=True)
        chol=cho_factor(C,lower=True);scaled=D/np.linalg.norm(D,axis=0);weighted=cho_solve(chol,scaled);normal=scaled.T@weighted
        normal=(normal+normal.T)/2;cholN=cho_factor(normal,lower=True);scales=np.linalg.norm(D,axis=0)
        V=cho_solve(cholN,np.eye(18))/scales[:,None]/scales[None,:];P=cho_solve(cholN,weighted.T)/scales[:,None]
        identity_error=float(np.linalg.norm(P@C@P.T-V)/np.linalg.norm(V));assert identity_error<1e-6
        row={'h':h,'tests':len(phi),'covariance':V.tolist(),'sqrt_trace_V':float(np.sqrt(np.trace(V))),
            'raw_gram_minimum_eigenvalue':float(values[0]),'residual_covariance_minimum_eigenvalue':float(eigh(C,eigvals_only=True)[0]),
            'scaled_normal_condition':float(np.linalg.cond(normal)),'sandwich_identity_relative_error':identity_error,
            'clean_quadrature_residual_norm':float(np.linalg.norm(b-D@truth))}
        rows.append(row);write(out/'reference.json',{'sigma_shape':base.tolist(),'rows':rows});print('Reference completed at h='+str(h),flush=True)
    relative=abs(rows[-1]['sqrt_trace_V']/rows[-2]['sqrt_trace_V']-1)
    assert relative<p['refinement_check_relative_tolerance']
    write(out/'completion.json',{'status':'complete','rows':len(rows),'refinement_relative_change':relative,'elapsed_seconds':time.perf_counter()-start})

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
