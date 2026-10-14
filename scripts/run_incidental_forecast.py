"""Kaggle-only deterministic weak-variance quadrature for crossing forecasts."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory
from scripts.run_crossfit_polynomial_probe import exact
from scripts.run_experiments import source_hash
from scripts.run_efficiency_validation import write


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Deterministic ODE computations run on Kaggle')
    config=json.loads((ROOT/'configs/incidental_forecast.json').read_text());p=json.loads((ROOT/'configs'/config['reference_protocol']).read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);q=np.array(p['truth'])
    model,_=vdp_pilot_trajectory(t,q,p['initial_center'],p['pilot_box']);base=model.std(axis=0)
    write(out/'manifest.json',dict(stage=stage,protocol=config,source_sha256=source_hash()))
    for windows in config['families']:
        phi,der,_=nested_bumps(t,windows,p['basis_relative_cutoff'])
        for order in config['initial_quadrature_orders']:
            nodes,weights=np.polynomial.legendre.leggauss(order)
            initials=[np.array(p['initial_center'])+p['initial_half_width']*np.array([a,b]) for a in nodes for b in nodes]
            weights=np.outer(weights,weights).ravel()/4
            for level in config['levels']:
                for method,scale,variance_factor in [('oracle_original',1.,1.),('oracle_fission',np.sqrt(2),.5)]:
                    informations=[]
                    for u in initials:
                        M,D,C=exact(t,q,u,scale*level*base,phi,der,p)
                        informations.append(np.linalg.inv(M@C@M.T))
                    meanI=np.einsum('r,rij->ij',weights,np.array(informations));V=variance_factor*np.linalg.inv(meanI)
                    rows.append(dict(windows=windows,order=order,level=level,method=method,
                        covariance_of_sqrt_n_estimate=V.tolist(),relative_variance=float(np.trace(V)/4)))
                    write(out/'forecast_variance.json',dict(rows=rows))
                    print(f'Deterministic variance {windows} tests, order {order}, noise {level}, {method}',flush=True)
    checks=[]
    for row in rows:
        if row['order']!=config['initial_quadrature_orders'][-1]:continue
        earlier=next(r for r in rows if r['windows']==row['windows'] and r['level']==row['level'] and r['method']==row['method'] and r['order']==config['initial_quadrature_orders'][0])
        loss=abs(earlier['relative_variance']/row['relative_variance']-1)
        checks.append(dict(windows=row['windows'],level=row['level'],method=row['method'],relative_difference=loss))
    write(out/'quadrature_checks.json',checks)
    assert max(c['relative_difference'] for c in checks)<=config['quadrature_relative_tolerance']
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',required=True,type=Path)
    args=a.parse_args();run(args.stage,args.output)
