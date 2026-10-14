"""Kaggle-only deterministic first incidental-bias coefficient and checks."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.profile_likelihood import ensemble_sensitivities,initial_curvature,whiten_records
from src.weak_efficiency import vdp_pilot_trajectory
from scripts.run_heterogeneous_records import trajectories
from scripts.run_experiments import source_hash
from scripts.run_efficiency_validation import write


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Numerical calculations run on Kaggle')
    p=json.loads((ROOT/'configs/incidental_curvature.json').read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);q=np.array(p['truth']);x,_=vdp_pilot_trajectory(t,q,p['initial_center'],[[-3,3]]*4);base=x.std(axis=0)
    write(out/'manifest.json',dict(stage=stage,protocol=p,source_sha256=source_hash()))
    for order in p['initial_gauss_legendre_orders']:
        nodes,w=np.polynomial.legendre.leggauss(order);u=np.array([[p['initial_center'][0]+p['initial_half_width']*a,p['initial_center'][1]+p['initial_half_width']*b] for a in nodes for b in nodes]);weights=np.outer(w,w).ravel()/4
        _,U=ensemble_sensitivities(t,q,u,p['ode_rtol'],p['ode_atol']);W=whiten_records(U,p['h'],p['tau'],base).reshape(len(u),-1,6)
        H=whiten_records(initial_curvature(t,q,u,p['ode_rtol'],p['ode_atol']),p['h'],p['tau'],base).reshape(len(u),-1,2,2)
        S=W[:,:,:4];Z=W[:,:,4:];N=np.einsum('rmi,rmj->rij',Z,Z);Ni=np.linalg.inv(N);cross=np.einsum('rmi,rmj->rij',Z,S)
        projected=S-np.einsum('rmi,rij,rjk->rmk',Z,Ni,cross);I=np.einsum('rmi,rmj->rij',projected,projected)
        curvature=np.einsum('rmij,rij->rm',H,Ni);term=np.einsum('rmi,rm->ri',projected,curvature)
        meanI=np.einsum('r,rij->ij',weights,I);coefficient=-.5*np.linalg.solve(meanI,np.einsum('r,ri->i',weights,term))
        rows.append(dict(quadrature_order=order,bias_coefficient=coefficient.tolist(),mean_information=meanI.tolist(),predictions=[dict(level=level,
            coefficient_bias=(level**2*coefficient).tolist(),relative_bias=float(level**2*np.linalg.norm(coefficient)/2),
            relative_information_sd=float(level*np.sqrt(np.trace(np.linalg.inv(meanI))/p['reference_records'])/2)) for level in p['noise_levels']]))
        print(f'Incidental curvature order {order}: complete',flush=True)
    probe=u[:2];H=initial_curvature(t,q,probe,p['ode_rtol'],p['ode_atol']);step=1e-3;checks=[]
    center=trajectories(t,q,probe)
    for j in range(2):
        delta=np.eye(2)[j]*step;finite=(trajectories(t,q,probe+delta)-2*center+trajectories(t,q,probe-delta))/step**2
        checks.append(float(np.linalg.norm(finite-H[:,:,:,j,j])/np.linalg.norm(H[:,:,:,j,j])))
    if max(checks)>1e-3:raise RuntimeError('Curvature finite-difference validation failed')
    write(out/'curvature.json',dict(results=rows,second_initial_sensitivity_relative_errors=checks,
        quadrature_relative_difference=float(np.linalg.norm(np.array(rows[1]['bias_coefficient'])-rows[0]['bias_coefficient'])/max(np.linalg.norm(rows[1]['bias_coefficient']),1e-15))))
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
