"""Validate recorded reviewer calculations, identities and publication wording."""
import argparse
import ast
import base64
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main(publication=False):
    out=ROOT/'results/kaggle_reviewer'
    assert json.loads((out/'completion.json').read_text())['status']=='complete'
    manifest=json.loads((out/'manifest.json').read_text())
    notebook=json.loads((ROOT/'kaggle/reviewer/wsindy_reviewer.ipynb').read_text())
    payload=None
    for cell in notebook['cells']:
        code=ast.parse(''.join(cell['source']))
        for stmt in code.body:
            if isinstance(stmt,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='payload' for t in stmt.targets):
                payload=ast.literal_eval(stmt.value)
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as bundle:
        assert bundle.testzip() is None
        fingerprint=hashlib.sha256()
        for folder in ('src','scripts','tests'):
            for name in sorted(n for n in bundle.namelist() if n.startswith(folder+'/') and n.endswith('.py')):
                fingerprint.update(name.encode());fingerprint.update(bundle.read(name))
        assert fingerprint.hexdigest()==manifest['source_sha256']
        assert json.loads(bundle.read('configs/reviewer.json'))==manifest['protocol']
    with zipfile.ZipFile(ROOT/'results/kaggle_reviewer_download/wsindy_reviewer_outputs.zip') as archive:
        assert archive.testzip() is None
    report=json.loads((ROOT/'reports/reviewer_results.json').read_text())
    assert report['replay_maximum_error_difference']<1e-8
    assert len(report['heterogeneous_bias_variance'])==9
    for r in report['heterogeneous_bias_variance']+report['adaptive_ou_finest_bias_variance']:
        assert abs(r['relative_rmse']**2-r['relative_bias']**2-r['relative_dispersion']**2)<1e-10
    info=report['vdp_information']
    assert max(info['finite_difference_relative_errors'])<1e-5
    for r in info['grids']:
        for I in (r['discrete'],r['continuous_quadrature']):
            C=np.array(I['active_covariance_bound']);W=np.array(I['initial_nuisance_covariance_bound'])
            assert np.linalg.eigvalsh(C).min()>0 and np.linalg.eigvalsh(W-C).min()>-1e-10
    for r in report['optimal_window_families']:
        assert 0<r['known_initial_efficiency']<=r['weak_class_efficiency']<=1+1e-6
        assert abs(r['normalization']-1)<2e-6 and abs(r['initial_annihilation'])<2e-6
    result={'completed':True,'recorded_execution_fingerprint_verified':True,'output_crc_verified':True,
            'sensitivities_checked':True,'information_psd_checked':True,'window_normalization_checked':True,
            'replay_3600_fits_matches_original':True,'bias_variance_identity_checked':True,'publication_checked':publication}
    if publication:
        from docx import Document
        import pymupdf
        article=(ROOT/'manuscript/article_revised.md').read_text(encoding='utf-8')
        supplement=(ROOT/'manuscript/supplement.md').read_text(encoding='utf-8')
        assert '### 4.3.' in article and '### 4.4.' in article
        assert '## Declaration of generative AI and AI-assisted technologies' in article
        assert '## Data availability' in article
        archive=json.loads((ROOT/'configs/public_archive.json').read_text())
        assert archive['version_url'] in article
        assert not re.search('abandoned grid diagnostic|initial custom-STLS pilot|historical|earlier 4,000-path experiment|fresh validation',article,re.I)
        assert 'malyshevdanil' not in supplement and 'credentials' not in supplement
        assert len(Document(ROOT/'manuscript/article_revised.docx').tables)==9
        with pymupdf.open(ROOT/'manuscript/article_revised.pdf') as pdf:
            text='\n'.join(p.get_text() for p in pdf)
            assert all(f'Table {i}.' in text for i in range(1,10))
            assert 'Data availability' in text and 'Van der Pol information reference' in text
        result['nine_tables_and_main_disclosures_verified']=True
    (ROOT/'reports/reviewer_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--publication',action='store_true')
    main(p.parse_args().publication)
