"""Validate frozen cross-fit runs and replay their archived record equations."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scripts.validate_nuisance_runs import snapshot
from src.crossfit_gmm import crossfit_equations
ROOT=Path(__file__).resolve().parents[1]


def validate(stages,full=False):
    checks=[];seen=[];rows_all=[]
    repair=stages[0].split('_')[0] if stages[0].startswith('repair') else None
    small=stages[0].startswith(('crossfit120_','repair120_'));main_a=stages[0].startswith(('crossfitA_','repairA'))
    config=repair+'_records.json' if repair else 'incidental_weak.json' if main_a else 'crossfit120_records.json' if small else 'crossfit_records.json'
    planned=json.loads((ROOT/'configs'/config).read_text());main_a=planned.get('matched_incidental_main',main_a)
    small=planned['windows']==120;method='crossfit'+str(planned['windows'])
    geometry_stage=planned.get('geometry_source_stage',repair+'_geometry' if repair else 'crossfitA_geometry' if main_a else 'crossfit120_geometry' if small else 'crossfit_geometry_v2')
    if full:
        checks.append(snapshot(geometry_stage,planned.get('geometry_protocol_file',config)))
        actual_cache_sha256=hashlib.sha256((ROOT/f'results/kaggle_{geometry_stage}/crossfit_instruments.npz').read_bytes()).hexdigest()
    old={r['repetition']:r['input_sha256'] for i in range(4) for r in json.loads((ROOT/f'results/kaggle_heterogeneous_records_{i}/input_hashes.json').read_text())}
    for stage in stages:
        checks.append(snapshot(stage,config));folder=ROOT/f'results/kaggle_{stage}'
        m=json.loads((folder/'manifest.json').read_text());p=m['protocol'];truth=np.array(p['truth'])
        completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete' and completion['discarded_records_or_bootstrap_fits']==0
        assert m['instrument_cache_sha256']==m['cache_geometry']['cache_sha256']
        if full:assert m['instrument_cache_sha256']==actual_cache_sha256
        assert m['cache_geometry']['cache_validated'] and p['external_training_records']==0
        rows=json.loads((folder/'rows.json').read_text());hashes=json.loads((folder/'input_hashes.json').read_text());rows_all+=rows
        assert [r['repetition'] for r in hashes]==m['repetitions'];seen+=m['repetitions']
        for h in hashes:
            if not main_a and p['seed_base']==3000000:assert h['first_800_input_sha256']==old[h['repetition']]
            elif main_a:
                matches=[]
                pattern='kaggle_inc_short_l*_n800_b*/input_hashes.json' if p['seed_base']==7100000 else 'kaggle_inc_main_l*_n800_b*/input_hashes.json'
                for f in (ROOT/'results').glob(pattern):
                    matches.extend(r['input_sha256'] for r in json.loads(f.read_text()) if r['repetition']==h['repetition'])
                assert all(value==h['first_800_input_sha256'] for value in matches)
        for rep in m['repetitions']:
            for level in p['levels']:
                for n in m['record_counts']:
                    assert sum(r['repetition']==rep and r['level']==level and r['n']==n and r['method']==method for r in rows)==1
                    if p.get('include_oracle_fission_comparison',small):
                        for label in ('original','fission'):
                            assert sum(r['repetition']==rep and r['level']==level and r['n']==n and r['method']=='oracle120_'+label for r in rows)==1
                assert sum(r['repetition']==rep and r['level']==level and r['method']=='split1000_50' for r in rows)==int(p.get('include_split_diagnostic',True))
        for r in rows:
            q=np.array(r['coefficients']);V=np.array(r['covariance']);assert np.isfinite(q).all() and np.linalg.eigvalsh(V).min()>0
            error=q-truth;assert np.isclose(error@np.linalg.solve(V,error),r['joint_statistic'])
            if r['method']==method:
                assert r['training_records_extra']==0 and r['inference_records']==r['n']
                for prefix in ('refit','frozen'):
                    if prefix=='frozen' and not p.get('compare_frozen_bootstrap',True):continue
                    lower=np.array(r[prefix+'_marginal_lower']);upper=np.array(r[prefix+'_marginal_upper'])
                    assert ((truth>=lower)&(truth<=upper)).tolist()==r[prefix+'_marginal_covered']
                    assert (r['joint_statistic']<=r[prefix+'_joint_critical'])==r[prefix+'_joint_covered']
        with np.load(folder/'crossfit_fixture.npz') as f:
            for n in m['record_counts']:
                q,V,_=crossfit_equations(f['A'][:n],f['b'][:n],f['Bcache'][:,:n],f['acache'][:,:n],f['anchor'],p['pilot_box'],np.arange(n)%p['folds'],cache_exponents=f['exponents'],cache_radius=float(f['radius']))
                r=next(r for r in rows if r['repetition']==m['repetitions'][0] and r['n']==n and r['level']==p['levels'][0] and r['method']==method)
                np.testing.assert_allclose(q[0],r['coefficients'],rtol=1e-10,atol=1e-11);np.testing.assert_allclose(V[0],r['covariance'],rtol=1e-9,atol=1e-15)
    if full:
        assert sorted(seen)==list(range(p['monte_carlo_repetitions']))
        assert len(rows_all)==p['monte_carlo_repetitions']*len(p['levels'])*(len(p['record_counts'])*(3 if p.get('include_oracle_fission_comparison',small) else 1)+int(p.get('include_split_diagnostic',True)))
    report=dict(completed=True,full_inferential_runs=full,stages=checks,input_hashes_match_previous_records=p['seed_base']==3000000,physical_seed_base=p['seed_base'],fold_equations_replayed=True,coverage_indicators_verified=True,rows=len(rows_all))
    prefix=repair if repair else 'crossfitA' if main_a else 'crossfit120' if small else 'crossfit'
    path=ROOT/('reports/'+prefix+('_run_validation.json' if full else '_smoke_validation.json'));path.write_text(json.dumps(report,indent=2));return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--smoke',action='store_true');args=parser.parse_args()
    print(json.dumps(validate(['crossfit_records_smoke'] if args.smoke else [f'crossfit_records_{i}' for i in range(4)],full=not args.smoke)))
