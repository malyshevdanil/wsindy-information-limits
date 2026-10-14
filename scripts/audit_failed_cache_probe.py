"""Verify a failed probe's frozen snapshot and audit each completed regime.

The raw outputs have no finished ZIP; their SHA256 is recorded separately.
Passing regimes may enter full geometry validation, never direct MC inference.
"""
import ast,base64,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    stage='cache_fix_probe';folder=ROOT/'results/failed_cache_fix_probe/results_cache_fix_probe'
    manifest=json.loads((folder/'manifest.json').read_text())
    notebook=json.loads((ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').read_text())
    payload=None
    for cell in notebook['cells']:
        for statement in ast.parse(''.join(cell['source'])).body:
            if isinstance(statement,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in statement.targets):
                payload=ast.literal_eval(statement.value)
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as bundle:
        assert bundle.testzip() is None
        fingerprint=hashlib.sha256()
        for directory in ('src','scripts','tests'):
            for name in sorted(n for n in bundle.namelist() if n.startswith(directory+'/') and n.endswith('.py')):
                fingerprint.update(name.encode());fingerprint.update(bundle.read(name))
        assert fingerprint.hexdigest()==manifest['source_sha256']
        assert json.loads(bundle.read('configs/cache_repair_probe.json'))==manifest['protocol']
        assert json.loads(bundle.read('configs/crossfit_cache_v3.json'))==manifest['cache_protocol']
    data=json.loads((folder/'probe.json').read_text());assert data['threshold']==.05
    families={}
    for family,windows,levels in [('repair120',120,[.005,.05]),('repair480',480,[.005,.05]),('repairA0',480,[.05])]:
        checks=[r for r in data['checks'] if r['windows']==windows and r['level'] in levels]
        assert len(checks)==13*len(levels)
        assert all(len({tuple(r['offset']) for r in checks if r['level']==level})==13 for level in levels)
        loss=max(r['variance_loss'] for r in checks);assert np.isfinite(loss)
        families[family]=dict(passed=loss<=.05,maximum_variance_loss=loss,checks=checks)
    report=dict(stage=stage,execution_source_and_protocol_verified=True,raw_output_sha256={
        name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ('manifest.json','probe.json')},
        archive_available=False,threshold=.05,families=families,
        scope='Geometry promotion only; each full geometry must independently pass all initial-state checks')
    (ROOT/'reports/cache_repair_probe_groups.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:{'passed':v['passed'],'maximum_loss':v['maximum_variance_loss']} for k,v in families.items()}))


if __name__=='__main__':main()
