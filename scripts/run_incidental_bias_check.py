"""Kaggle-only finite-n shared-parameter curvature diagnostic."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.profile_likelihood import ensemble_hessian,ensemble_sensitivities,whiten_records
from src.weak_efficiency import vdp_pilot_trajectory
from scripts.run_experiments import source_hash
from scripts.run_efficiency_validation import write


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('ODE calculations run on Kaggle')
    config=json.loads((ROOT/'configs/incidental_bias_check.json').read_text());p=json.loads((ROOT/'configs/incidental_curvature.json').read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);q=np.array(p['truth']);model,_=vdp_pilot_trajectory(t,q,p['initial_center'],[[-3,3]]*4);base=model.std(axis=0)
    write(out/'manifest.json',dict(stage=stage,protocol=config,source_sha256=source_hash()))
    for order in config['initial_quadrature_orders']:
        nodes,w=np.polynomial.legendre.leggauss(order);weights=np.outer(w,w).ravel()/4
        u=np.array([np.array(p['initial_center'])+p['initial_half_width']*np.array([a,b]) for a in nodes for b in nodes])
        _,U,H=ensemble_hessian(t,q,u,p['ode_rtol'],p['ode_atol'])
        W=whiten_records(U,p['h'],p['tau'],base).reshape(len(u),-1,6)
        H=whiten_records(H,p['h'],p['tau'],base).reshape(len(u),-1,6,6)
        S=W[:,:,:4];Z=W[:,:,4:];Ni=np.linalg.inv(np.einsum('rmi,rmj->rij',Z,Z))
        NiC=Ni@np.einsum('rmi,rmj->rij',Z,S);projected=S-Z@NiC
        individualI=np.einsum('rmi,rmj->rij',projected,projected);meanI=np.einsum('r,rij->ij',weights,individualI);Ii=np.linalg.inv(meanI)
        h=np.einsum('rmij,rij->rm',H[:,:,4:,4:],Ni)
        terms=np.einsum('rmi,rm->ri',projected,h);meanterm=np.einsum('r,ri->i',weights,terms)
        incidental=-.5*Ii@meanterm
        tangent=np.zeros((len(u),6,4));tangent[:,:4]=np.eye(4);tangent[:,4:]=-NiC
        covariance=np.einsum('rai,ij,rbj->rab',tangent,Ii,tangent)
        shared_h=np.einsum('rmab,rab->rm',H,covariance)
        shared=-.5*Ii@np.einsum('r,rmi,rm->i',weights,projected,shared_h)
        # Delta-method correction for randomness of the average information and
        # curvature term: E[Ibar^-1 tbar] differs from I^-1 E[t].
        deltaI=individualI-meanI;deltat=terms-meanterm
        adjustment=np.einsum('r,rij,jk,rkl,l->i',weights,deltaI,Ii,deltaI,Ii@meanterm)
        adjustment-=np.einsum('r,rij,jk,rk->i',weights,deltaI,Ii,deltat)
        random_design=-.5*Ii@adjustment
        rows.append(dict(order=order,incidental_coefficient=incidental.tolist(),shared_parameter_finite_n_coefficient=shared.tolist(),
            random_initial_design_finite_n_coefficient=random_design.tolist(),
            predictions=[dict(n=n,level=level,coefficient_bias=(level**2*(incidental+(shared+random_design)/n)).tolist()) for n in config['record_counts'] for level in config['levels']]))
        write(out/'bias_check.json',dict(rows=rows))
        print(f'Full-joint curvature order {order}: complete',flush=True)
    errors=[];probe=u[:2];_,_,analytic=ensemble_hessian(t,q,probe,p['ode_rtol'],p['ode_atol']);step=1e-4
    for coordinate in (0,3,4):
        qp=q.copy();qm=q.copy();up=probe.copy();um=probe.copy()
        if coordinate<4:qp[coordinate]+=step;qm[coordinate]-=step
        else:up[:,coordinate-4]+=step;um[:,coordinate-4]-=step
        _,plus=ensemble_sensitivities(t,qp,up,p['ode_rtol'],p['ode_atol']);_,minus=ensemble_sensitivities(t,qm,um,p['ode_rtol'],p['ode_atol'])
        finite=(plus-minus)/(2*step);expected=analytic[...,coordinate]
        errors.append(dict(coordinate=coordinate,relative_error=float(np.linalg.norm(finite-expected)/np.linalg.norm(expected))))
    write(out/'finite_difference_checks.json',errors)
    assert max(c['relative_error'] for c in errors)<=config['relative_finite_difference_tolerance']
    for key in ('incidental_coefficient','shared_parameter_finite_n_coefficient','random_initial_design_finite_n_coefficient'):
        assert np.linalg.norm(np.array(rows[0][key])-rows[1][key])/max(np.linalg.norm(rows[1][key]),1e-15)<=config['relative_quadrature_tolerance']
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',required=True,type=Path)
    args=a.parse_args();run(args.stage,args.output)
