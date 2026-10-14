"""Compare archived estimation errors to the matching deterministic covariance."""
import json,hashlib,shutil,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    folder=ROOT/'results/guarded_geometry_reference';reference=json.loads((folder/'reference.json').read_text())
    from scripts.collect_stopped_runs import snapshot
    audit=snapshot('guarded_geometry_reference',json.loads((folder/'manifest.json').read_text()))
    completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete'
    V=np.array(reference['rows'][-1]['covariance']);scale=np.sqrt(np.trace(V));truth=np.array(json.loads((ROOT/'configs/guarded_geometry_reference.json').read_text())['coefficients'])
    levels=json.loads((ROOT/'configs/independent_validation_v2.json').read_text())['geometry']['levels'];rows=[]
    for i,epsilon in enumerate(levels):
        records=json.loads((ROOT/f'data/independent_validation/iq_geo_{i}/geometry_rows.json').read_text())
        selected=sorted((r for r in records if r['method']=='self_instrumented_guarded'),key=lambda r:r['seed'])
        assert [r['seed'] for r in selected]==list(range(100))
        e=(np.array([r['coefficients'] for r in selected])-truth)/epsilon;loss=np.sum(e*e,axis=1)
        resamples=np.random.default_rng(271828+i).integers(0,100,(10000,100))
        ratio=np.sqrt(loss.mean())/scale;ci=np.quantile(np.sqrt(loss[resamples].mean(axis=1))/scale,[.025,.975])
        cov=np.cov(e,rowvar=False);bias=e.mean(axis=0)
        rows.append({'epsilon':epsilon,'paths':100,'scaled_rmse':float(np.sqrt(loss.mean())),
            'reference_sqrt_trace_V':float(scale),'rmse_ratio':float(ratio),'pointwise_bootstrap95_ratio':ci.tolist(),
            'scaled_bias_norm':float(np.linalg.norm(bias)),'empirical_covariance':cov.tolist(),
            'covariance_relative_frobenius_error':float(np.linalg.norm(cov-V)/np.linalg.norm(V)),
            'mean_reference_mahalanobis_per_parameter':float(np.mean(np.sum(e*np.linalg.solve(V,e.T).T,axis=1))/18),
            'fallbacks':sum(r['diagnostics']['fallback'] for r in selected),'clipping':sum(r['diagnostics'].get('clipped',False) for r in selected)})
    result={'source_audit':audit,'deterministic_completion':completion,'tests':134,'library_coefficients':18,
        'sqrt_trace_V':float(scale),'reference_covariance':V.tolist(),'rows':rows,
        'resampling_scope':'postprocessing archived coefficient errors; 10000 fixed-seed record resamples, pointwise intervals, no new observation paths',
        'interpretation':'finite-family reference only; ratio near one does not establish a joint growing-family limit or arbitrary-model efficiency'}
    (ROOT/'reports/guarded_geometry_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    bundle=ROOT/'data/guarded_geometry_reference';bundle.mkdir(exist_ok=True);inventory=[]
    for original in sorted(folder.glob('*.json')):
        target=bundle/original.name;shutil.copyfile(original,target)
        inventory.append({'path':target.name,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
    (bundle/'sha256_manifest.json').write_text(json.dumps({'original_bytes_preserved':True,'files':inventory},indent=2),encoding='utf-8')
    print(json.dumps([{k:v for k,v in r.items() if k!='empirical_covariance'} for r in rows],indent=2))

if __name__=='__main__':main()
