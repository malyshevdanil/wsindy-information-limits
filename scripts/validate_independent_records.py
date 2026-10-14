"""Verify replicate source snapshots, seed counts, fixtures and interval bookkeeping."""
import ast,base64,hashlib,io,json,zipfile,sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.record_gmm import fit_mean


def load(p):return json.loads(p.read_text(encoding='utf-8'))


def main():
    report=load(ROOT/'reports/independent_records_results.json');p=report['protocol'];checked=[];all_rows=[]
    for shard in range(p['shards']):
        stage=f'independent_records_{shard}';folder=ROOT/f'results/kaggle_{stage}';manifest=load(folder/'manifest.json')
        nb=load(ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb');payload=None
        for cell in nb['cells']:
            for statement in ast.parse(''.join(cell['source'])).body:
                if isinstance(statement,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in statement.targets):payload=ast.literal_eval(statement.value)
        with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:
            assert archive.testzip() is None;digest=hashlib.sha256()
            for directory in ['src','scripts','tests']:
                for name in sorted(n for n in archive.namelist() if n.startswith(directory+'/') and n.endswith('.py')):digest.update(name.encode());digest.update(archive.read(name))
            if digest.hexdigest()!=manifest['source_sha256']:
                from scripts.source_archive_audit import verify_redacted_source
                verify_redacted_source(archive,stage,manifest['source_sha256'])
            assert json.loads(archive.read('configs/independent_records.json'))==p==manifest['protocol']
        with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as archive:
            assert archive.testzip() is None
            for entry in archive.infolist():
                if not entry.is_dir():assert archive.read(entry.filename)==(folder/entry.filename).read_bytes()
        expected=list(range(shard,p['monte_carlo_repetitions'],p['shards']));hashes=load(folder/'input_hashes.json')
        assert [r['repetition'] for r in hashes]==expected
        assert [r['seed'] for r in hashes]==[p['seed_base']+r for r in expected]
        rows=load(folder/'rows.json');all_rows+=rows
        for n in p['record_counts']:
            for method in ['gmm_known','gmm_estimated','stacked_known']:
                assert [r['repetition'] for r in rows if r['n']==n and r['method']==method]==expected
        with np.load(folder/'sufficient_statistics_fixture.npz') as fixture:
            digest=hashlib.sha256()
            for name in ['raw','b','vs','lag','previous']:digest.update(np.ascontiguousarray(fixture[name]).tobytes())
            assert digest.hexdigest()==hashes[0]['sufficient_statistics_sha256']
            raw=(fixture['raw'][0::2]+fixture['raw'][1::2])/2;b=(fixture['b'][0::2]+fixture['b'][1::2])/2
            variance=np.array(load(folder/'geometry.json')['sigma'])**2
            for n in p['record_counts']:
                B=n//2
                for method,estimated in [('gmm_known',False),('gmm_estimated',True)]:
                    v=fixture['vs'][:B].mean(axis=0) if estimated else variance
                    q,V,_=fit_mean(raw[:B],b[:B],fixture['vs'][:B],v,estimated)
                    original=next(r for r in rows if r['n']==n and r['method']==method and r['repetition']==expected[0])
                    assert np.allclose(q,original['coefficients'],rtol=1e-11,atol=1e-12)
                    assert np.allclose(V,original['covariance'],rtol=1e-10,atol=1e-15)
        checked.append(dict(stage=stage,execution_fingerprint_verified=True,source_protocol_verified=True,archive_crc_and_bytes_verified=True,fixture_recomputed=True))
    qtrue=np.array(p['truth'])
    for r in all_rows:
        q=np.array(r['coefficients']);assert np.isfinite(q).all()
        if r['method'].startswith('gmm'):
            V=np.array(r['covariance']);assert np.linalg.eigvalsh(V).min()>0
            error=q-qtrue;stat=float(error@np.linalg.solve(V,error));assert np.isclose(stat,r['joint_statistic'])
            assert (stat<=chi2.ppf(.95,4))==r['joint_covered']
            assert (np.abs(error/np.sqrt(np.diag(V)))<=1.959963984540054).tolist()==r['marginal_covered']
    assert len(all_rows)==p['monte_carlo_repetitions']*len(p['record_counts'])*3
    for r in report['groups']:assert np.isclose(r['relative_rmse']**2,r['relative_bias']**2+r['relative_dispersion']**2)
    result=dict(completed=True,stages=checked,all_repetitions_retained=True,coverage_and_sandwich_verified=True,
        independent_records=report['independent_records_simulated'],planned_estimates=len(all_rows))
    (ROOT/'reports/independent_records_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':main()
