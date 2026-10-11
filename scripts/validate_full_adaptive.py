"""Verify archived full-pipeline executions and optionally the final release."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import zipfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main(publication=False,smoke_only=False):
    cfg=json.loads((ROOT/'configs/full_adaptive.json').read_text())
    pins=json.loads((ROOT/'configs/author_pins.json').read_text())
    stages=['adaptive_smoke'] if smoke_only else ['adaptive_smoke','adaptive_classical',*[f'adaptive_wm_{i}' for i in range(4)]]
    records=[];hashes={};fit_count=0
    for stage in stages:
        folder=ROOT/f'results/kaggle_{stage}'
        manifest=json.loads((folder/'manifest.json').read_text())
        assert manifest['completed'] and manifest['location']=='Kaggle'
        if stage!='adaptive_smoke':assert manifest['protocol']==cfg
        provenance=json.loads((folder/'author_provenance.json').read_text())
        for name,pin in pins.items():
            assert provenance[name]['commit']==pin['commit']
            if 'zip_sha256' in pin:assert provenance[name]['zip_sha256']==pin['zip_sha256']
            else:
                for file,expected in pin['files'].items():assert provenance[name]['files_verified'][file]==expected
        rows=json.loads((folder/'rows.json').read_text())
        for row in rows:
            with np.load(folder/'inputs'/(row['case']+'.npz')) as data:
                actual=hashlib.sha256(data['t'].tobytes()+data['y'].tobytes()).hexdigest()
                assert actual==row['input_sha256']
                if stage!='adaptive_smoke':hashes.setdefault(row['case'],set()).add(actual)
            if row['method']=='wmsindy' and row['status']=='ok':
                fit=json.loads((folder/'fits/wmsindy'/row['case']/'fit.json').read_text())
                assert len(fit['history'])==8
                assert all(h['training_steps']==5000 for h in fit['history'])
                assert [h['optimizer_iteration'] for h in fit['history']]==list(range(5000,40001,5000))
                assert all(len(h['mts_after'])==2 for h in fit['history'])
                with np.load(folder/'fits/wmsindy'/row['case']/'training_arrays.npz') as arrays:
                    assert arrays['iterates'].shape==(8,9,2)
                    assert np.all(np.isfinite(arrays['noise_estimate']))
        if stage=='adaptive_smoke':assert len(rows)==3 and all(r['status']=='ok' for r in rows)
        else:fit_count+=len(rows)
        with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as archive:
            assert archive.testzip() is None
            assert not any('author_sources/' in n for n in archive.namelist())
        records.append({'stage':stage,'fits':len(rows),'failed':sum(r['status']!='ok' for r in rows),
                        'manifest_verified':True,'source_pins_verified':True,'input_hashes_verified':True})
    if not smoke_only:
        assert len(hashes)==cfg['repeats']*4 and all(len(h)==1 for h in hashes.values())
        assert fit_count==7*4*cfg['repeats']
    result={'completed':True,'comparative_fits':fit_count,'paired_cases':len(hashes),
            'full_wm_budget_verified':True,'stages':records,'publication_checked':publication}
    if publication:
        from docx import Document
        import pymupdf
        doc=Document(ROOT/'manuscript/article_revised.docx')
        source=(ROOT/'manuscript/article_revised.md').read_text(encoding='utf-8')
        table_count=len(re.findall(r'^Table \d+\.',source,flags=re.M))
        assert len(doc.inline_shapes)==7 and len(doc.tables)==table_count
        with pymupdf.open(ROOT/'manuscript/article_revised.pdf') as pdf:
            text='\n'.join(p.get_text() for p in pdf)
            assert 'Complete adaptive author pipelines' in text
            assert all(f'Table {i}.' in text for i in range(1,table_count+1))
        for item in json.loads((ROOT/'dist/manifest.json').read_text()):
            path=ROOT/'dist'/item['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256']
            with zipfile.ZipFile(path) as archive:
                assert archive.testzip() is None
                for entry in archive.infolist():
                    assert not any(part in ('author_repos','author_sources','.kaggle') for part in Path(entry.filename).parts)
                    assert Path(entry.filename).name not in ('access_token','kaggle.json','.env')
                    if entry.file_size<30_000_000:
                        assert not re.search(rb'KGAT_[A-Za-z0-9]{20,}',archive.read(entry))
        result['release_crc_sha_and_secret_scan']=True
    (ROOT/('reports/full_adaptive_smoke_validation.json' if smoke_only else 'reports/full_adaptive_validation.json')).write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--publication',action='store_true');parser.add_argument('--smoke-only',action='store_true')
    args=parser.parse_args();main(args.publication,args.smoke_only)
