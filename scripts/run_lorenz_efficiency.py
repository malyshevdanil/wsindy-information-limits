"""Three-state polynomial check with exact OU weak-residual covariance."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash
from src.weak_efficiency import nested_bumps,gls_map
from src.weak_features import quadrature_weights
from src.covariance import exponential_apply
from src.noise import ou_noise
from src.path_scoring import ou_whiten


def jacobians(x):
    N=len(x);J=np.zeros((4,N,3,3));u,v,w=x.T
    J[0,:,1,0]=-w;J[0,:,1,1]=-1;J[0,:,1,2]=-u
    J[0,:,2,0]=v;J[0,:,2,1]=u
    J[1,:,0,0]=-1;J[1,:,0,1]=1;J[2,:,1,0]=1;J[3,:,2,2]=-1
    return J


def integrate(t,q,initial,sensitivities=False):
    def rhs(_,state):
        u,v,w=state[:3];s,r,b=q;dx=np.array([s*(v-u),r*u-v-u*w,u*v-b*w])
        if not sensitivities:return dx
        J=jacobians(state[:3][None,:])[:,0];J=J[0]+np.einsum('q,qij->ij',q,J[1:])
        U=state[3:].reshape(3,6);dU=J@U
        dU[:,:3]+=np.diag([v-u,u,-w])
        return np.r_[dx,dU.ravel()]
    start=np.asarray(initial)
    if sensitivities:
        U=np.zeros((3,6));U[:,3:]=np.eye(3);start=np.r_[start,U.ravel()]
    sol=solve_ivp(rhs,(0,float(t[-1])),start,t_eval=t,method='DOP853',rtol=1e-10,atol=1e-12)
    if not sol.success:raise RuntimeError('Lorenz integration failed')
    return (sol.y[:3].T,sol.y[3:].T.reshape(len(t),3,6)) if sensitivities else sol.y.T


def moments(t,y,phi,der,weights=None):
    w=quadrature_weights(t) if weights is None else weights;u,v,z=y.T;K=len(phi)
    A=np.zeros((3*K,3));A[0::3,0]=(phi*w)@(v-u);A[1::3,1]=(phi*w)@u;A[2::3,2]=-(phi*w)@z
    b=-(der*w)@y;b[:,1]+=(phi*w)@(v+u*z);b[:,2]-=(phi*w)@(u*v)
    return A,b.ravel()


def blocks(t,model,phi,der,sigma,tau,first_order=False,weights=None):
    J=jacobians(model);N=len(t);K=len(phi);w=quadrature_weights(t) if weights is None else weights;B={}
    rho=np.exp(-(t[1]-t[0])/tau)
    for coordinate in range(3):
        arrays={}
        for k in range(4):
            U=-phi.T[:,:,None]*J[k,:,None,:,coordinate]
            if k==0:U[:,:,coordinate]-=der.T
            U=(U*w[:,None,None]).reshape(N,3*K)
            if np.any(U):arrays[k]=U
        applied={k:exponential_apply(U,rho) for k,U in arrays.items()}
        for k,U in arrays.items():
            for l in arrays:
                if l<k:continue
                key=(k,l);value=sigma[coordinate]**2*(U.T@applied[l])
                B[key]=B.get(key,0)+value
    if not first_order:
        P=(phi*w).T;H=P.T@exponential_apply(P,rho*rho)
        B[0,0][1::3,1::3]+=sigma[0]**2*sigma[2]**2*H
        B[0,0][2::3,2::3]+=sigma[0]**2*sigma[1]**2*H
    return B


def assemble(q,B):
    weights=np.r_[1,q];C=np.zeros_like(B[0,0])
    for (k,l),value in B.items():C+=weights[k]*weights[l]*(value if k==l else value+value.T)
    return (C+C.T)/2


def fit(A,b,q0,B,p):
    q=q0.copy()
    for it in range(p['iterations']):
        P,_=gls_map(A,assemble(q,B),p['covariance_relative_relaxation']);new=P@b
        change=float(np.linalg.norm(new-q)/np.linalg.norm(q));q=new
        if change<p['tolerance']:break
    return q,dict(converged=change<p['tolerance'],iterations=it+1,last_change=change)


def run(out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical runs must execute on Kaggle')
    p=json.loads((ROOT/'configs/lorenz_efficiency.json').read_text(encoding='utf-8'));out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    write(out/'manifest.json',dict(stage='efficiency_lorenz',protocol=p,source_sha256=source_hash()))
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);truth=np.array(p['truth']);norm=np.linalg.norm(truth)
    x,U=integrate(t,truth,p['initial'],True);sigma=p['level']*np.std(x,axis=0)
    information={}
    for unknown in [False,True]:
        V=U if unknown else U[:,:,:3];W=ou_whiten(V,p['h'],p['tau'],sigma).reshape(-1,V.shape[-1])
        scale=np.linalg.norm(W,axis=0);C=np.linalg.inv((W/scale).T@(W/scale))/scale[:,None]/scale[None,:]
        information['unknown' if unknown else 'known']=dict(covariance=C[:3,:3].tolist(),relative_rmse=float(np.sqrt(np.trace(C[:3,:3]))/norm),
            scaled_information_condition=float(np.linalg.cond((W/scale).T@(W/scale))))
    write(out/'information.json',dict(**information,sigma=sigma.tolist()))
    geometry=[];bases={}
    for K in p['geometry_windows']:
        phi,der,meta=nested_bumps(t,K,p['basis_relative_cutoff']);A,b=moments(t,x,phi,der)
        B=blocks(t,x,phi,der,sigma,p['tau'],True);P,_=gls_map(A,assemble(truth,B),p['covariance_relative_relaxation'])
        C=P@assemble(truth,B)@P.T
        geometry.append(dict(windows=K,relative_rmse=float(np.sqrt(np.trace(C))/norm),covariance=C.tolist(),
            clean_quadrature_relative_error=float(np.linalg.norm(P@b-truth)/norm),normalization_error=float(np.linalg.norm(P@A-np.eye(3))),**meta))
        if K in (30,p['fit_windows']):bases[K]=(phi,der,meta)
        print(f'Lorenz geometry K={K}',flush=True)
    write(out/'geometry.json',geometry)
    rows=[];phi0,der0,_=bases[30];phi,der,meta=bases[p['fit_windows']]
    for seed in range(p['paths']):
        y=x+ou_noise(len(t),p['h'],p['tau'],sigma,np.random.default_rng(p['seed_base']+seed))
        A0,b0=moments(t,y,phi0,der0);pilot=np.linalg.lstsq(A0,b0,rcond=None)[0]
        q0=np.clip(pilot,np.array(p['pilot_box'])[:,0],np.array(p['pilot_box'])[:,1]);model=integrate(t,q0,p['initial'])
        A,b=moments(t,y,phi,der);B=blocks(t,model,phi,der,sigma,p['tau']);q,diag=fit(A,b,q0,B,p)
        rows.append(dict(seed=seed,coefficients=q.tolist(),relative_error=float(np.linalg.norm(q-truth)/norm),
            pilot_clipped=bool(np.any(q0!=pilot)),**diag,**meta))
        if seed%10==0:write(out/'rows.json',rows);print(f'Lorenz {seed+1}/{p["paths"]}',flush=True)
    write(out/'rows.json',rows);write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();run(args.output)
