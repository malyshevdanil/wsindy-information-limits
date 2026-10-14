"""Frozen bounded polynomial cache, checked against exact model GLS before MC."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.crossfit_cache import powers,features,design_nodes
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory,vdp_moments
from scripts.run_crossfit_polynomial_probe import exact
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Calculations run on Kaggle')
    config=(stage.split('_')[0]+'_records.json' if stage.startswith('repair') else
            'incidental_weak.json' if stage=='crossfitA_geometry' else
            'crossfit120_records.json' if stage=='crossfit120_geometry' else 'crossfit_records.json')
    p=json.loads((ROOT/'configs'/config).read_text());c=json.loads((ROOT/'configs'/p.get('cache_protocol_file','crossfit_cache.json')).read_text());out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);anchor=np.mean(p['pilot_box'],axis=1);radius=c['radius']
    model,_=vdp_pilot_trajectory(t,anchor,p['initial_center'],p['pilot_box']);base=model.std(axis=0)
    phi,der,meta=nested_bumps(t,p['windows'],p['basis_relative_cutoff']);checks=[];arrays={}
    exponents=powers(c['degree']);nodes=design_nodes(c,anchor)
    inverse=np.linalg.pinv(features(nodes,anchor,radius,exponents));xs=p['instrument_grid']['x'];ys=p['instrument_grid']['v']
    write(out/'manifest.json',dict(stage=stage,protocol=p,cache_protocol=c,source_sha256=source_hash()))
    for level in p['levels']:
        sigma=np.sqrt(2)*level*base;cache=np.empty((len(exponents),len(xs),len(ys),4,2*len(phi)))
        for i,ux in enumerate(xs):
            for j,uv in enumerate(ys):
                maps=[exact(t,q,[ux,uv],sigma,phi,der,p)[0] for q in nodes]
                cache[:,i,j]=np.einsum('pn,nim->pim',inverse,np.array(maps))
                # Restore anchor information so pooling preserves the relative
                # precision of different trajectories at the declared anchor.
                Ma,Da,Ca=exact(t,anchor,[ux,uv],sigma,phi,der,p)
                info=np.linalg.inv(Ma@Ca@Ma.T)
                cache[:,i,j]=np.einsum('ij,pjm->pim',info,cache[:,i,j])
                if stage.startswith('repair'):
                    np.savez_compressed(out/f'checkpoint_{level}_{i}_{j}.npz',coefficients=cache[:,i,j])
                print(f'Crossfit polynomial cache {level}: grid {i*len(ys)+j+1}/{len(xs)*len(ys)}; elapsed {time.perf_counter()-start:.1f}s',flush=True)
        arrays['cache_'+str(level)]=cache
        if p.get('include_oracle_fission_comparison'):
            for label,scale in (('original',1.),('fission',np.sqrt(2))):
                oracle=np.empty((len(xs),len(ys),4,2*len(phi)))
                for i,ux in enumerate(xs):
                    for j,uv in enumerate(ys):
                        Mo,Do,Co=exact(t,np.asarray(p['truth']),[ux,uv],scale*level*base,phi,der,p)
                        oracle[i,j]=np.linalg.solve(Mo@Co@Mo.T,Mo)
                arrays['oracle_'+label+'_'+str(level)]=oracle
        validation_offsets=[[0.,0.,0.,0.]]+p['cache_validation_offsets']
        validation_offsets+=list(radius*np.eye(4))+list(-radius*np.eye(4))
        for delta in validation_offsets:
            delta=np.asarray(delta)
            q=anchor+delta
            for u in p['cache_validation_initials']:
                Me,D,C=exact(t,q,u,sigma,phi,der,p)
                approx=np.einsum('p,pim->im',features(q,anchor,radius,exponents),cache[:,xs.index(u[0]),ys.index(u[1])]);approx=np.linalg.solve(approx@D,approx)
                loss=float(np.trace(approx@C@approx.T)/np.trace(Me@C@Me.T)-1)
                checks.append(dict(level=level,initial=u,offset=delta.tolist(),variance_loss=loss))
        write(out/'checks.json',checks)
    maximum=max(v['variance_loss'] for v in checks);validated=maximum<=p['maximum_cache_variance_loss']
    write(out/'geometry.json',dict(anchor=anchor.tolist(),anchor_uses_training_records=False,sigma_base=base.tolist(),family=meta,
        cache_protocol=c,aggregation_scale='anchor model moment information, fixed before MC',checks=checks,maximum_variance_loss=maximum,threshold=p['maximum_cache_variance_loss'],cache_validated=validated))
    np.savez_compressed(out/'crossfit_instruments.npz',anchor=anchor,base=base,exponents=exponents,radius=radius,**arrays)
    if not validated:raise RuntimeError('Prespecified cache accuracy check failed; diagnostic cache retained')
    write(out/'completion.json',dict(status='complete',elapsed_seconds=time.perf_counter()-start))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args();run(args.stage,args.output)
