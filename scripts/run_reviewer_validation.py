"""Reviewer calculations on Kaggle; no algorithm training is repeated."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from scripts.run_revision_fast import batch_vdp
from scripts.run_revision_validation import setup, oracle
from src.covariance import exponential_apply
from src.noise import ou_noise
from src.systems import trajectory
from src.weak_features import test_functions, quadrature_weights, weak_matrices


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def augmented_rhs(_,state,theta):
    x,v=state[:2];a,b,c,d=theta
    f=np.array([a*v,b*x+c*v+d*x*x*v])
    J=np.array([[0.,a],[b+2*d*x*v,c+d*x*x]])
    B=np.array([[v,0.,0.,0.],[0.,x,v,x*x*v]])
    S=state[2:10].reshape(2,4);Z=state[10:].reshape(2,2)
    return np.r_[f,(J@S+B).ravel(),(J@Z).ravel()]


def vdp_information(out,cfg):
    p=cfg['vdp'];T=p['T'];theta=np.array(p['truth']);tau=p['tau']
    fine=min(p['information_h']);t=np.linspace(0,T,round(T/fine)+1)
    initial=np.r_[p['initial'],np.zeros(8),np.eye(2).ravel()]
    sol=solve_ivp(lambda z,u:augmented_rhs(z,u,theta),(0,T),initial,t_eval=t,
                  method='DOP853',rtol=p['rtol'],atol=p['atol'])
    assert sol.success
    state=sol.y.T;S=state[:,2:10].reshape(-1,2,4);Z=state[:,10:].reshape(-1,2,2)
    A=np.concatenate([S,Z],axis=2)
    derivatives=np.array([augmented_rhs(ti,u,theta) for ti,u in zip(t,state)])
    Adot=np.concatenate([derivatives[:,2:10].reshape(-1,2,4),derivatives[:,10:].reshape(-1,2,2)],axis=2)
    tc=np.linspace(0,T,round(T/p['calibration_h'])+1)
    sigma=p['level']*np.std(trajectory('vdp',tc),axis=0)
    # Calibration is held fixed under parameter perturbation: covariance is not a parameter.
    outputs=[]
    for h in p['information_h']:
        stride=round(h/fine);ts=t[::stride];sens=A[::stride];sd=Adot[::stride]
        innovations=sens[1:]-np.exp(-h/tau)*sens[:-1]
        discrete=np.einsum('di,dj,d->ij',sens[0],sens[0],1/sigma**2)
        discrete+=np.einsum('ndi,ndj,d->ij',innovations,innovations,
                            1/(sigma**2*(-np.expm1(-2*h/tau))))
        Q=sd+sens/tau;w=quadrature_weights(ts)
        continuous=np.einsum('di,dj,d->ij',sens[0],sens[0],1/sigma**2)
        continuous+=np.einsum('ndi,ndj,n,d->ij',Q,Q,w,tau/(2*sigma**2))
        row={'h':h}
        for name,I in [('discrete',discrete),('continuous_quadrature',continuous)]:
            I=(I+I.T)/2;known=I[:4,:4]
            effective=known-I[:4,4:]@np.linalg.solve(I[4:,4:],I[4:,:4])
            C=np.linalg.inv(known);CW=np.linalg.inv(effective)
            assert np.linalg.eigvalsh(CW-C).min()>-1e-10
            direction=np.array([0.,0.,1.,-1.])
            row[name]={'joint_information':I.tolist(),'active_covariance_bound':C.tolist(),
                'initial_nuisance_covariance_bound':CW.tolist(),
                'relative_rmse_bound':float(np.sqrt(np.trace(C))/np.linalg.norm(theta)),
                'weak_initial_annihilation_relative_rmse_bound':float(np.sqrt(np.trace(CW))/np.linalg.norm(theta)),
                'mu_only_standard_deviation_bound':float(1/np.sqrt(direction@known@direction)),
                'condition_number':float(np.linalg.cond(known))}
        outputs.append(row)
    # Independent central-difference checks of all four coefficient sensitivities.
    checks=[];eps=p['finite_difference_step'];indices=np.linspace(0,len(t)-1,101,dtype=int)
    for j in range(4):
        solutions=[]
        for sign in (-1,1):
            th=theta.copy();th[j]+=sign*eps
            r=solve_ivp(lambda z,u:augmented_rhs(z,np.r_[u,np.zeros(8),np.zeros(4)],th)[:2],
                        (0,T),p['initial'],t_eval=t[indices],method='DOP853',rtol=p['rtol'],atol=p['atol'])
            assert r.success;solutions.append(r.y.T)
        fd=(solutions[1]-solutions[0])/(2*eps)
        discrepancy=float(np.linalg.norm(fd-S[indices,:,j])/np.linalg.norm(S[indices,:,j]))
        assert discrepancy<1e-5;checks.append(discrepancy)
    assert abs(outputs[-1]['continuous_quadrature']['relative_rmse_bound']/outputs[-2]['continuous_quadrature']['relative_rmse_bound']-1)<1e-4
    write(out/'vdp_information.json',{'parameters':p,'sigma_fixed':sigma.tolist(),'grids':outputs,
        'finite_difference_relative_errors':checks,
        'scope':'Local covariance bound for regular locally unbiased estimation of four independently free active coefficients with fixed known covariance; not a bound on median errors or biased sparse discovery.'})
    np.savez_compressed(out/'vdp_sensitivities.npz',t=t,clean=state[:,:2],S=S,Z=Z,Sdot=Adot[:,:,:4],Zdot=Adot[:,:,4:],sigma=sigma)
    print('VDP information complete',flush=True)


def window_families(out,cfg):
    p=cfg['scalar'];a,x0,sigma,tau,T=[p[k] for k in ('a','x0','sigma','tau','T')]
    outputs=[];arrays={}
    for n in p['grid_intervals']:
        t=np.linspace(0,T,n+1);w=quadrature_weights(t);x=x0*np.exp(a*t)
        for family in p['families']:
            phis=[];ders=[]
            for L in family['widths']:
                centers=([2.] if family['count']==1 else np.linspace(L/2+.005,T-L/2-.005,family['count']))
                phi,der=test_functions(t,centers,L);phis.extend(phi);ders.extend(der)
            phi=np.array(phis);der=np.array(ders);D=(phi*w)@x
            U=((-der-a*phi)*w).T
            C=sigma**2*(U.T@exponential_apply(U,np.exp(-(t[1]-t[0])/tau)))
            C=(C+C.T)/2;vals,V=np.linalg.eigh(C)
            for cutoff in p['relative_eigenvalue_cutoffs']:
                keep=vals>cutoff*vals.max();Q=V[:,keep]
                invD=Q@((Q.T@D)/vals[keep]);information=float(D@invD)
                coef=invD/information;variance=float(coef@C@coef)
                ell=((-der-a*phi).T@coef);weighted=ell*w
                norm=float(weighted@(t*x));annihilation=float(weighted@(x/x0))
                assert abs(norm-1)<2e-6 and abs(annihilation)<2e-6
                outputs.append({'family':family['name'],'grid_intervals':n,'windows':len(D),
                    'eigenvalue_cutoff':cutoff,'rank':int(keep.sum()),'variance':variance,
                    'optimal_linear_information':information,'normalization':norm,'initial_annihilation':annihilation})
                if n==max(p['grid_intervals']) and cutoff==min(p['relative_eigenvalue_cutoffs']):
                    arrays[family['name']+'_weights']=coef
                    arrays[family['name']+'_covariance']=C
                    arrays[family['name']+'_normalizing_moments']=D
    write(out/'window_families.json',{'parameters':p,'results':outputs,
        'scope':'Oracle covariance-weighted optimal linear combination within each declared finite bump span, using true a and covariance; not global optimization over every compactly supported function and not variance of the unguarded ratio.'})
    np.savez_compressed(out/'window_family_arrays.npz',**arrays)
    print('Window optimization complete',flush=True)


def replay_heterogeneous(out,cfg):
    p=cfg['heterogeneous'];t,x,powers,truth,centers,L=setup('vdp',.02)
    sigma=.05*np.std(x,axis=0);rows=[]
    for block in range(p['blocks']):
        initials=[];noises=[]
        for r in range(max(p['recordings'])):
            rng=np.random.default_rng(p['seed_base']+16*block+r)
            initials.append(np.array([2.,0.])+rng.uniform(-.25,.25,size=2))
            noises.append(ou_noise(len(t),.02,.3,sigma,rng))
        clean=batch_vdp(t,initials);matrices=[];records=[]
        for r in range(max(p['recordings'])):
            y=clean[:,r,:]+noises[r];matrices.append(weak_matrices(t,y,powers,centers,L,sigma**2));records.append(y)
            M=r+1
            if M not in p['recordings']:continue
            pairs=[(np.vstack([v[0] for v in matrices]),np.vstack([v[1] for v in matrices])),
                   (np.mean([v[0] for v in matrices],axis=0),np.mean([v[1] for v in matrices],axis=0)),
                   weak_matrices(t,np.mean(records,axis=0),powers,centers,L,sigma**2/M)]
            for method,(G,b) in zip(('stacked','mean_weak','naive_mean_sensor'),pairs):
                estimate=oracle(G,b,truth)
                rows.append({'block':block,'M':M,'method':method,
                    'error':float(np.linalg.norm(estimate-truth)/np.linalg.norm(truth)),
                    'coefficients':estimate.tolist()})
        if block%50==0:print(f'Coefficient replay block {block+1}/400',flush=True)
    write(out/'heterogeneous_coefficients.json',{'truth':truth.tolist(),'rows':rows})


def run(out):
    if not str(ROOT).startswith('/kaggle/working/'):
        raise RuntimeError('Numerical experiments must run on Kaggle')
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    cfg=json.loads((ROOT/'configs/reviewer.json').read_text())
    write(out/'manifest.json',{'stage':'reviewer','protocol':cfg,'source_sha256':source_hash()})
    vdp_information(out,cfg);window_families(out,cfg);replay_heterogeneous(out,cfg)
    write(out/'completion.json',{'status':'complete','elapsed_seconds':time.perf_counter()-start})
    print('Reviewer validation complete',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
