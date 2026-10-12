"""Audit new execution snapshots, hashes, record equations and interval coverage."""
import ast,base64,hashlib,io,json,sys,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.heterogeneous_gmm import record_fit


def load(p):return json.loads(p.read_text(encoding='utf-8'))


def snapshot(stage,config):
    folder=ROOT/f'results/kaggle_{stage}';manifest=load(folder/'manifest.json');nb=load(ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb');payload=None
    for cell in nb['cells']:
        for statement in ast.parse(''.join(cell['source'])).body:
            if isinstance(statement,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in statement.targets):payload=ast.literal_eval(statement.value)
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as bundle:
        assert bundle.testzip() is None;fingerprint=hashlib.sha256()
        for directory in ('src','scripts','tests'):
            for name in sorted(n for n in bundle.namelist() if n.startswith(directory+'/') and n.endswith('.py')):
                fingerprint.update(name.encode());fingerprint.update(bundle.read(name))
        assert fingerprint.hexdigest()==manifest['source_sha256']
        assert json.loads(bundle.read('configs/'+config))==manifest['protocol']
    with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as archive:
        assert archive.testzip() is None
        for name in archive.namelist():
            if not name.endswith('/'):assert archive.read(name)==(folder/name).read_bytes()
    return dict(stage=stage,execution_source_and_protocol_verified=True,archive_crc_and_bytes_verified=True)


def main():
    stages=[]
    for i in range(4):stages.append(snapshot(f'record_bootstrap_{i}','record_bootstrap.json'))
    report=load(ROOT/'reports/record_bootstrap_results.json');assert report['observations_identical_to_original_sha256'] and report['point_estimates_and_covariances_identical']
    assert report['bootstrap_fits']==4995000 and report['failed_or_discarded_draws']==0
    if (ROOT/'reports/heterogeneous_records_results.json').exists():
        report=load(ROOT/'reports/heterogeneous_records_results.json');p=report['protocol'];truth=np.array(p['truth']);count=0
        for i in range(4):
            stage=f'heterogeneous_records_{i}';stages.append(snapshot(stage,'heterogeneous_records.json'));folder=ROOT/f'results/kaggle_{stage}'
            rows=load(folder/'rows.json');expected=list(range(i,p['monte_carlo_repetitions'],p['shards']));hashes=load(folder/'input_hashes.json')
            assert [r['repetition'] for r in hashes]==expected
            for n in p['record_counts']:
                for level in p['levels']:
                    for method in p['methods']:assert [r['repetition'] for r in rows if r['n']==n and r['level']==level and r['method']==method]==expected
            for r in rows:
                q=np.array(r['coefficients']);V=np.array(r['covariance']);error=q-truth;assert np.isfinite(q).all() and np.linalg.eigvalsh(V).min()>0
                assert np.isclose(error@np.linalg.solve(V,error),r['joint_statistic'])
                if r['method'].startswith('fission'):
                    lower=np.array(r['bootstrap_marginal_lower']);upper=np.array(r['bootstrap_marginal_upper'])
                    assert ((truth>=lower)&(truth<=upper)).tolist()==r['bootstrap_marginal_covered']
                    assert (r['joint_statistic']<=r['bootstrap_joint_critical'])==r['bootstrap_joint_covered']
                count+=1
            with np.load(folder/'equation_fixture.npz') as fixture:
                for method in p['methods']:
                    for n in p['record_counts']:
                        q,V=record_fit(fixture[method+'_B'][:n],fixture[method+'_a'][:n]);r=next(r for r in rows if r['repetition']==expected[0] and r['level']==.05 and r['method']==method and r['n']==n)
                        assert np.allclose(q,r['coefficients'],rtol=1e-11,atol=1e-12) and np.allclose(V,r['covariance'],rtol=1e-10,atol=1e-15)
        assert count==p['monte_carlo_repetitions']*len(p['record_counts'])*len(p['levels'])*len(p['methods'])
        for g in report['groups']:assert np.isclose(g['relative_rmse']**2,g['relative_bias']**2+g['relative_dispersion']**2)
    if (ROOT/'reports/moderate_likelihood_results.json').exists():
        stages.append(snapshot('moderate_likelihood','moderate_likelihood.json'))
        likelihood=load(ROOT/'reports/moderate_likelihood_results.json')
        assert likelihood['observation_hashes_match_author_128']
        assert all(g['paths']==128 and g['bound_active_paths']==0 for g in likelihood['groups'])
        assert hashlib.sha256((ROOT/'results/kaggle_moderate_likelihood_download/wsindy_moderate_likelihood_outputs.zip').read_bytes()).hexdigest()==likelihood['output_sha256']
    result=dict(completed=True,stages=stages,bootstrap_data_identity_verified=True,coverage_and_equation_fixtures_verified=True)
    (ROOT/'reports/nuisance_run_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))


if __name__=='__main__':main()
