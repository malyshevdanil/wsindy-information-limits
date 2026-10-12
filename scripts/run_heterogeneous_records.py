"""Model-weighted heterogeneous-record benchmark with independent instruments."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor,cho_solve
from scipy.signal import lfilter
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.heterogeneous_gmm import interpolate_grid,record_fit,record_bootstrap
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory,vdp_moments,covariance_blocks,assemble_covariance
from src.weak_features import test_functions,quadrature_weights
from src.path_scoring import vdp_sensitivities,ou_whiten
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def trajectories(t,q,u):
    def rhs(_,state):
        x,v=state.reshape(-1,2).T;a,b,c,d=q
        return np.column_stack((a*v,b*x+c*v+d*x*x*v)).ravel()
    sol=solve_ivp(rhs,(0.,t[-1]),np.asarray(u).ravel(),t_eval=t,method='DOP853',rtol=1e-10,atol=1e-12)
    if not sol.success:raise RuntimeError('Ensemble trajectory integration failed')
    return sol.y.T.reshape(len(t),len(u),2).transpose(1,0,2)


def noise(rng,shape,sigma,rho):
    innovations=rng.normal(size=shape)*sigma*np.sqrt(1-rho*rho);innovations[:,0]/=np.sqrt(1-rho*rho)
    return lfilter([1],[1,-rho],innovations,axis=1)


def moments(y,P,D,variance):
    x,v=y[:,:,0],y[:,:,1];R=len(y);K=len(P);A=np.zeros((R,2*K,4));b=np.empty((R,2*K))
    A[:,0::2,0]=v@P.T;A[:,1::2,1]=x@P.T;A[:,1::2,2]=v@P.T;A[:,1::2,3]=((x*x-variance[0])*v)@P.T
    b[:,0::2]=-x@D.T;b[:,1::2]=-v@D.T
    return A,b


def instruments(t,q,xs,ys,phi,der,sigma,tau,relaxation):
    result=[];max_defect=0.
    for ux in xs:
        row=[]
        for uv in ys:
            model,_=vdp_pilot_trajectory(t,q,[ux,uv],[[-3,3]]*4)
            D,b=vdp_moments(t,model,phi,der,np.zeros(2));max_defect=max(max_defect,float(np.linalg.norm(b-D@q)))
            blocks,H=covariance_blocks(t,model,phi,der,sigma,tau);C=assemble_covariance(q,blocks,H,sigma)
            C+=relaxation*np.trace(C)/len(C)*np.eye(len(C))
            L=cho_solve(cho_factor(C,lower=True,check_finite=False),D,check_finite=False).T
            row.append(L)
        result.append(row)
    return np.array(result),max_defect


def information(t,q,u,sigma,tau):
    _,U=vdp_sensitivities(t,q,u);W=ou_whiten(U,t[1]-t[0],tau,sigma).reshape(-1,6);I=W.T@W
    return I[:4,:4]-I[:4,4:]@np.linalg.solve(I[4:,4:],I[4:,:4])


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Run numerical experiments on Kaggle')
    p=json.loads((ROOT/'configs/heterogeneous_records.json').read_text());out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();shard=int(stage.rsplit('_',1)[1]);truth=np.array(p['truth']);t=np.linspace(0,p['T'],round(p['T']/p['h'])+1)
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash()))
    baseline,_=vdp_pilot_trajectory(t,truth,p['initial_center'],[[-3,3]]*4);base=np.std(baseline,axis=0);rho=np.exp(-p['h']/p['tau']);w=quadrature_weights(t)
    rng=np.random.default_rng(p['training_seed']);u=np.array(p['initial_center'])+rng.uniform(-p['initial_half_width'],p['initial_half_width'],(p['training_records'],2))
    x=trajectories(t,truth,u);y=x+noise(rng,x.shape,p['training_level']*base,rho)
    phi,der=test_functions(t,[2.1,5.2,8.3,11.4,14.5,17.6],4,4)
    A,b=moments(y,phi*w,der*w,(p['training_level']*base)**2);pilot=np.linalg.lstsq(A.mean(axis=0),b.mean(axis=0),rcond=None)[0]
    pilot=np.clip(pilot,np.array(p['pilot_box'])[:,0],np.array(p['pilot_box'])[:,1]);del x,y,A,b
    xs=p['instrument_grid']['x'];ys=p['instrument_grid']['v'];families={};library={};oracle={};defects={}
    for K in p['windows']:
        phi,der,meta=nested_bumps(t,K,p['basis_relative_cutoff']);families[K]=(phi*w,der*w,meta)
        for level in p['levels']:
            library[K,level],defects[f'{K}_{level}']=instruments(t,pilot,xs,ys,phi,der,np.sqrt(2)*level*base,p['tau'],p['covariance_relative_relaxation'])
            if K==120:oracle[level],_=instruments(t,truth,xs,ys,phi,der,level*base,p['tau'],p['covariance_relative_relaxation'])
    m=p['information_grid_points'];ix=np.linspace(1.75,2.25,m);iv=np.linspace(-.25,.25,m)
    info=np.array([[information(t,truth,[ux,uv],.05*base,p['tau']) for uv in iv] for ux in ix])
    checks=[];probe=np.random.default_rng(3300000).uniform([1.75,-.25],[2.25,.25],(12,2))
    for u in probe:
        exact=information(t,truth,u,.05*base,p['tau']);approx=interpolate_grid(info,ix,iv,u[None])[0]
        checks.append(float(np.linalg.norm(exact-approx)/np.linalg.norm(exact)))
    if max(checks)>.003:raise RuntimeError('Information interpolation error exceeds prespecified 0.3%')
    write(out/'geometry.json',dict(training_coefficients=pilot.tolist(),training_relative_error=float(np.linalg.norm(pilot-truth)/2),
        sigma_base=base.tolist(),families={str(K):families[K][2] for K in families},maximum_model_weak_defects=defects,
        information_interpolation_relative_errors=checks,information_grid=info.tolist(),information_x=ix.tolist(),information_v=iv.tolist(),
        covariance_model='full Prop3 including higher Hermite orders, at independent training coefficient pilot',
        covariance_known=True,initial_states_supplied_to_feasible_instruments=False))
    np.savez_compressed(out/'instrument_fixture.npz',pilot=pilot,**{f'L{K}_{level}':L for (K,level),L in library.items()})
    rows=[];hashes=[];reps=list(range(shard,p['monte_carlo_repetitions'],p['shards']))
    for rep in reps:
        rng=np.random.default_rng(p['seed_base']+rep);R=max(p['record_counts'])
        u=np.array(p['initial_center'])+rng.uniform(-p['initial_half_width'],p['initial_half_width'],(R,2));x=trajectories(t,truth,u)
        eta=noise(rng,x.shape,base,rho);xi=noise(rng,x.shape,base,rho)
        digest=hashlib.sha256()
        for array in (u,x,eta,xi):digest.update(np.ascontiguousarray(array).tobytes())
        hashes.append(dict(repetition=rep,input_sha256=digest.hexdigest()))
        Ir=interpolate_grid(info,ix,iv,u)
        for level in p['levels']:
            y=x+level*eta;ya=y+level*xi;yb=y-level*xi;equations={}
            for K,(P,D,meta) in families.items():
                Aa,ba=moments(ya,P,D,2*(level*base)**2);Ab,bb=moments(yb,P,D,2*(level*base)**2)
                La=interpolate_grid(library[K,level],xs,ys,ya[:,0]);Lb=interpolate_grid(library[K,level],xs,ys,yb[:,0])
                B=.5*(np.einsum('rpm,rmq->rpq',La,Ab)+np.einsum('rpm,rmq->rpq',Lb,Aa))
                a=.5*(np.einsum('rpm,rm->rp',La,bb)+np.einsum('rpm,rm->rp',Lb,ba));equations[f'fission{K}']=(B,a)
                if K==120:
                    A,b=moments(y,P,D,(level*base)**2);L=interpolate_grid(oracle[level],xs,ys,u)
                    equations['oracle120']=(np.einsum('rpm,rmq->rpq',L,A),np.einsum('rpm,rm->rp',L,b))
            for method,(B,a) in equations.items():
                for n in p['record_counts']:
                    q,V=record_fit(B[:n],a[:n]);error=q-truth;sd=np.sqrt(np.diag(V));stat=float(error@np.linalg.solve(V,error))
                    if not np.isfinite(q).all() or np.linalg.eigvalsh(V).min()<=0:raise RuntimeError('Invalid fit; no records discarded')
                    bound=np.sqrt(np.trace(np.linalg.inv(Ir[:n].sum(axis=0)*(0.05/level)**2)))/2
                    row=dict(repetition=rep,n=n,level=level,method=method,coefficients=q.tolist(),covariance=V.tolist(),
                        relative_error=float(np.linalg.norm(error)/2),local_information_relative_sd=float(bound),
                        wald_marginal_covered=(np.abs(error/sd)<=1.959963984540054).tolist(),
                        wald_joint_covered=bool(stat<=chi2.ppf(.95,4)),joint_statistic=stat)
                    if method.startswith('fission'):
                        # Method-specific K, independent of dictionary iteration order.
                        seed=p['bootstrap_seed_base']+rep*10000+n+int(method[7:])*10+round(level*1000)
                        limits,critical=record_bootstrap(B[:n],a[:n],q,np.random.default_rng(seed),p['bootstrap_draws'])
                        lower=q-limits[1]*sd;upper=q-limits[0]*sd
                        row.update(bootstrap_seed=seed,bootstrap_t_quantiles=limits.tolist(),bootstrap_joint_critical=float(critical),
                            bootstrap_marginal_lower=lower.tolist(),bootstrap_marginal_upper=upper.tolist(),
                            bootstrap_marginal_covered=((truth>=lower)&(truth<=upper)).tolist(),bootstrap_joint_covered=bool(stat<=critical))
                    rows.append(row)
            if rep==reps[0] and level==.05:
                np.savez_compressed(out/'equation_fixture.npz',**{f'{method}_{key}':array for method,arrays in equations.items() for key,array in zip(('B','a'),arrays)})
        if len(hashes)%5==0:
            write(out/'rows.json',rows);print(f'Heterogeneous records shard {shard}: {len(hashes)}/{len(reps)} repetitions',flush=True)
    write(out/'rows.json',rows);write(out/'input_hashes.json',hashes)
    write(out/'completion.json',dict(status='complete',repetitions=len(reps),rows=len(rows),elapsed_seconds=time.perf_counter()-start,
        records_simulated=len(reps)*R,levels_paired=True,discarded_records_or_bootstrap_fits=0))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);p=a.parse_args();run(p.stage,p.output)
