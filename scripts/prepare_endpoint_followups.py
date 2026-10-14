"""Prepare frozen references from complete archives, without new simulations."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def write(path,data):path.write_text(json.dumps(data,indent=2),encoding='utf-8')


def main():
    parent=json.loads((ROOT/'configs/support_extra_term.json').read_text())
    result=json.loads((ROOT/'data/endpoint_confirmation/support_extra_term_results.json').read_text())
    rows=json.loads((ROOT/'data/endpoint_confirmation/support_extra/rows.json').read_text())
    assert result['completed'] and len(rows)==1000
    reference={str(level):next(r['delta_min_local'] for r in result['information_references'] if r['level']==level and r['label']=='full_library_boundary_free') for level in parent['levels']}
    write(ROOT/'configs/support_wald_references.json',{'delta_references':reference,'information_references':result['information_references'],
          'rows':[{**r,'multiplier':parent['delta_multipliers_of_weak_reference'][r['delta_index']]} for r in rows]})
    mlepath=ROOT/'results/kaggle_moderate_likelihood/rows.json';mle=json.loads(mlepath.read_text())
    weakpaths=[ROOT/f'results/kaggle_efficiency_vdp_final_{i}/rows.json' for i in range(2)]
    weak=[r for path in weakpaths for r in json.loads(path.read_text()) if r['level']==.05 and r['windows']==480 and r['weight_model']=='ou']
    indexed={r['seed']:{'seed':r['seed'],'archived_boundary_free':r} for r in weak};assert len(indexed)==128
    for r in mle:
        record=indexed[r['seed']];record[r['method']+'_mle']=r
        if 'input_sha256' in record:assert record['input_sha256']==r['input_sha256']
        record['input_sha256']=r['input_sha256']
    assert all(all(k in r for k in ('known_initial_mle','unknown_initial_mle','input_sha256')) for r in indexed.values())
    write(ROOT/'configs/boundary_vdp5_references.json',{'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [mlepath]+weakpaths},
          'rows':[indexed[s] for s in sorted(indexed)]})
    from scripts.build_efficiency_kaggle import build
    stages=['boundary_general_vdp5','boundary_general_lorenz']+[f'support_wald_{i}' for i in range(4)]
    for stage in stages:
        build(stage)
        assert (ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').stat().st_size < 1_000_000


if __name__=='__main__':main()
