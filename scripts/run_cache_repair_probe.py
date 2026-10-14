"""Kaggle-only endpoint-complete probe before expensive cache rebuilds."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.crossfit_cache import powers,features,design_nodes
from src.weak_efficiency import nested_bumps,vdp_pilot_trajectory
from scripts.run_crossfit_polynomial_probe import exact
from scripts.run_efficiency_validation import write
from scripts.run_experiments import source_hash


def run(stage,out):
    if not str(ROOT).startswith('/kaggle/working/'):raise RuntimeError('Calculations run on Kaggle')
    p=json.loads((ROOT/'configs/crossfit_records.json').read_text())
    design_file=('cache120_main_probe.json' if stage=='cache120a' else 'cache_repair_probe_v4.json' if stage=='cache_fix_v4' else 'cache_repair_probe_A.json' if stage=='cache_fix_a' else 'cache_repair_probe.json')
    design=json.loads((ROOT/'configs'/design_file).read_text())
    c=json.loads((ROOT/'configs'/design.get('cache_protocol_file','crossfit_cache_v3.json')).read_text())
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();checks=[]
    t=np.linspace(0,p['T'],round(p['T']/p['h'])+1);anchor=np.mean(p['pilot_box'],axis=1);radius=c['radius']
    model,_=vdp_pilot_trajectory(t,anchor,p['initial_center'],p['pilot_box']);base=model.std(axis=0)
    exponents=powers(c['degree']);nodes=design_nodes(c,anchor)
    inverse=np.linalg.pinv(features(nodes,anchor,radius,exponents))
    offsets=[[0.]*4]+p['cache_validation_offsets']+list(radius*np.eye(4))+list(-radius*np.eye(4))
    offsets+=design.get('extra_validation_offsets',[])
    write(out/'manifest.json',dict(stage=stage,protocol=design,cache_protocol=c,source_sha256=source_hash()))
    for windows in design['families']:
        phi,der,_=nested_bumps(t,windows,p['basis_relative_cutoff'])
        for level in design['levels']:
            sigma=np.sqrt(2)*level*base;maps=[]
            for i,q in enumerate(nodes):
                maps.append(exact(t,q,design['initial'],sigma,phi,der,p)[0])
                if (i+1)%128==0:print(f'Cache repair probe {windows} {level}: {i+1}/{len(nodes)}',flush=True)
            coefficients=np.einsum('pn,nim->pim',inverse,np.array(maps))
            np.savez_compressed(out/f'probe_coefficients_{windows}_{level}.npz',coefficients=coefficients)
            for delta in offsets:
                q=anchor+delta;Me,D,C=exact(t,q,design['initial'],sigma,phi,der,p)
                approximate=np.einsum('p,pim->im',features(q,anchor,radius,exponents),coefficients)
                approximate=np.linalg.solve(approximate@D,approximate)
                loss=float(np.trace(approximate@C@approximate.T)/np.trace(Me@C@Me.T)-1)
                checks.append(dict(windows=windows,level=level,offset=np.asarray(delta).tolist(),variance_loss=loss))
            write(out/'probe.json',dict(checks=checks,maximum_variance_loss=max(v['variance_loss'] for v in checks),threshold=design['threshold']))
    passed=max(v['variance_loss'] for v in checks)<=design['threshold']
    write(out/'completion.json',dict(status='complete' if passed else 'accuracy_failure',elapsed_seconds=time.perf_counter()-start))
    if not passed:raise RuntimeError('Repair probe failed unchanged 5% threshold; full rebuild forbidden')


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',required=True);a.add_argument('--output',type=Path,required=True)
    args=a.parse_args();run(args.stage,args.output)
