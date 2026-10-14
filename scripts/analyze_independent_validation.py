"""Analyze only complete frozen independent cells; never inspect partial rates."""
import json,sys,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.stats import beta,norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.collect_stopped_runs import snapshot

def interval(k,n):
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]

def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(center-half),float(center+half)]

def analyze():
    protocol=json.loads((ROOT/'configs/independent_validation_v2.json').read_text())
    result={'protocol':protocol['protocol'],'complete_cells':{},'geometry':{},'partial_rates_reported':False}
    families={'small':[f'iq_nom_{i}' for i in range(4)],'moderate_null':[f'iq_bootn_{i:02d}' for i in range(50)],'moderate_power':[f'iq_bootp_{i:02d}' for i in range(10)]}
    for family,stages in families.items():
        if not all((ROOT/'results'/s/'completion.json').exists() for s in stages):continue
        rows=[]
        for stage in stages:
            folder=ROOT/'results'/stage;manifest=json.loads((folder/'manifest.json').read_text())
            assert manifest['protocol']==protocol;snapshot(stage,manifest)
            assert json.loads((folder/'completion.json').read_text())['status']=='complete'
            rows.extend(json.loads((folder/'rows.json').read_text()))
        for cell in sorted({r['cell'] for r in rows}):
            for method in sorted({r['method'] for r in rows}):
                selected=sorted((r for r in rows if r['cell']==cell and r['method']==method),key=lambda r:r['seed'])
                n=(protocol['small_null_paths'] if cell==0 else protocol['small_power_paths']) if family=='small' else protocol['moderate_null_paths'] if cell==0 else protocol['moderate_power_paths']
                assert [r['seed'] for r in selected]==list(range(n))
                for r in selected:
                    assert np.isfinite(r['se']) and r['se']>0
                    assert r['rejected']==(abs(r['candidate_coefficient']/r['se'])>norm.ppf(.975))
                k=sum(r['rejected'] for r in selected);ci=interval(k,n)
                wc=wilson(k,n)
                row={'paths':n,'rejections':k,'rate':k/n,'clopper_pearson95':ci,'wilson95':wc,'wilson_compatible_with_5_percent':wc[0]<=.05<=wc[1]}
                if cell==0:row['within_frozen_size_band']=ci[0]>=.03 and ci[1]<=.07
                elif family=='small':
                    multiplier=protocol['small_power_multipliers'][cell-1];d=2*norm.ppf(.95)*multiplier
                    prediction=float(norm.cdf(-norm.ppf(.975)-d)+norm.sf(norm.ppf(.975)-d))
                    row.update(local_prediction=prediction,within_frozen_power_band=ci[0]>=prediction-.05 and ci[1]<=prediction+.05,size_validation_required=True)
                result['complete_cells'][f'{family}_c{cell}_{method}']=row
    for i,epsilon in enumerate(protocol['geometry']['levels']):
        folder=ROOT/'results'/f'iq_geo_{i}'
        if not (folder/'completion.json').exists():continue
        manifest=json.loads((folder/'manifest.json').read_text());snapshot(f'iq_geo_{i}',manifest);assert manifest['protocol']==protocol
        records=json.loads((folder/'geometry_rows.json').read_text());truth=np.zeros(18);truth[[1,9,10,15]]=[1,-1,1,-1]
        for method in protocol['geometry']['methods']:
            rows=sorted((r for r in records if r['method']==method),key=lambda r:r['seed'])
            assert [r['seed'] for r in rows]==list(range(100))
            error=np.array([r['coefficients'] for r in rows])-truth
            result['geometry'][f'{epsilon}_{method}']={'paths':100,'scaled_bias':(error.mean(axis=0)/epsilon).tolist(),
                'scaled_rmse':float(np.sqrt(np.mean(np.sum(error**2,axis=1)))/epsilon),
                'scaled_fourth_error_moment':float(np.mean(np.sum(error**2,axis=1)**2)/epsilon**4),
                'fallbacks':sum(r['diagnostics']['fallback'] for r in rows),
                'clipped':sum(r['diagnostics'].get('clipped',False) for r in rows)}
    result['all_studies_complete']=len(result['complete_cells'])==7 and len(result['geometry'])==10
    # Three small cells, two methods in each of the two moderate cells.
    completed_stages=[s for group in families.values() for s in group if (ROOT/'results'/s/'completion.json').exists()]
    completed_stages += [f'iq_geo_{i}' for i in range(5) if (ROOT/'results'/f'iq_geo_{i}'/'completion.json').exists()]
    # Preserve complete shards of the budget-limited study without examining
    # their rejection rates or presenting them as a completed primary cell.
    for stage in completed_stages:
        folder=ROOT/'results'/stage
        manifest=json.loads((folder/'manifest.json').read_text());assert manifest['protocol']==protocol;snapshot(stage,manifest)
        completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete'
        if stage.startswith('iq_boot'):
            draws=list(folder.glob('bootstrap_c*.json'));assert len(draws)==20
            assert completion['rows']==40 and all(len(json.loads(p.read_text())['targets'])==399 for p in draws)
    result['completed_bootstrap_shards']=[s for s in completed_stages if s.startswith('iq_boot')]
    if completed_stages:
        bundle=ROOT/'data/independent_validation';inventory=[]
        for stage in completed_stages:
            for original in sorted((ROOT/'results'/stage).glob('*.json')):
                target=bundle/stage/original.name;target.parent.mkdir(parents=True,exist_ok=True)
                if target.exists():assert target.read_bytes()==original.read_bytes()
                else:shutil.copyfile(original,target)
                inventory.append({'path':target.relative_to(bundle).as_posix(),'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
        (bundle/'manifest.json').write_text(json.dumps({'original_bytes_preserved':True,'all_studies_complete':result['all_studies_complete'],'completed_stages':completed_stages,'files':inventory},indent=2))
    path=ROOT/'reports/independent_validation_results.json';path.write_text(json.dumps(result,indent=2,allow_nan=False))
    return result

if __name__=='__main__':print(json.dumps(analyze(),indent=2))
