"""Verify frozen sources, downloaded bytes and all follow-up seed counts."""
import ast,base64,hashlib,io,json,zipfile,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
    follow=load(ROOT/'reports/efficiency_followup_results.json');info=load(ROOT/'reports/information_followup_results.json')
    assert 'matched' in follow and 'refinement' in info,'Required experiments unfinished'
    stages=list(follow['execution_sources'])+list(info['execution_sources'])+['efficiency_fixture'];checked=[]
    for stage in stages:
        folder=ROOT/f'results/kaggle_{stage}';manifest=load(folder/'manifest.json');assert load(folder/'completion.json')['status']=='complete'
        nb=load(ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb');payload=None
        for cell in nb['cells']:
            for statement in ast.parse(''.join(cell['source'])).body:
                if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='payload' for t in statement.targets):payload=ast.literal_eval(statement.value)
        with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as bundle:
            assert bundle.testzip() is None
            fingerprint=hashlib.sha256()
            for directory in ['src','scripts','tests']:
                for name in sorted(n for n in bundle.namelist() if n.startswith(directory+'/') and n.endswith('.py')):fingerprint.update(name.encode());fingerprint.update(bundle.read(name))
            if fingerprint.hexdigest()!=manifest['source_sha256']:
                from scripts.source_archive_audit import verify_redacted_source
                verify_redacted_source(bundle,stage,manifest['source_sha256'])
            if 'pilot_files' in manifest:
                for name,digest in manifest['pilot_files'].items():assert hashlib.sha256(bundle.read(name)).hexdigest()==digest
            if 'pilot_sha256' in manifest:assert hashlib.sha256(bundle.read('configs/refinement_pilots.npz')).hexdigest()==manifest['pilot_sha256']
            if 'compatibility_hashes' in manifest:
                for name,digest in manifest['compatibility_hashes'].items():assert hashlib.sha256(bundle.read('src/author_compat/'+name)).hexdigest()==digest
        archive_path=ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip'
        with zipfile.ZipFile(archive_path) as archive:
            assert archive.testzip() is None
            for entry in archive.infolist():
                if not entry.is_dir():assert archive.read(entry.filename)==(folder/entry.filename).read_bytes()
        checked.append(dict(stage=stage,source_sha256=manifest['source_sha256'],output_sha256=hashlib.sha256(archive_path.read_bytes()).hexdigest()))
    assert follow['original_128_coefficient_max_discrepancy']<1e-12
    assert follow['all_power_converged'] and all(r['converged']==128 for r in info['noise'])
    for r in list(follow['power'].values())+list(info['fisher'].values())+list(info['refinement'].values())+list(info['noise'])+[info['lorenz']]:
        assert abs(r['relative_rmse']**2-r['relative_bias']**2-r['relative_dispersion']**2)<1e-12
    fixture=load(ROOT/'results/kaggle_efficiency_fixture/checks.json');assert max(fixture.values())<1e-12
    for r in info['lorenz']['geometry']:
        assert np.linalg.eigvalsh(r['covariance']).min()>0 and r['normalization_error']<1e-10
    result=dict(completed=True,stages=checked,unique_paired_seeds_verified=True,
        recorded_one_step_and_followup_refinement_separate=True,all_downloaded_bytes_verified=True,
        exact_covariance_and_ou_score_fixture=fixture)
    (ROOT/'reports/efficiency_followup_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
if __name__=='__main__':main()
