"""Create allowlisted submission and reproducibility archives; exclude secrets."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript';DIST=ROOT/'dist';DIST.mkdir(exist_ok=True)
for suffix in ('.md','.docx','.tex','.pdf'):
    shutil.copy2(BASE/('article_revised'+suffix),BASE/('article_en'+suffix))
shutil.copy2(BASE/'publication_references.bib',BASE/'references.bib')
sources=[]
for folder,pattern in [('src','*.py'),('scripts','*.py'),('tests','*.py'),('configs','*.json')]:sources+=list((ROOT/folder).glob(pattern))
sources+=list((ROOT/'src/author_compat').glob('*.m'))
sources += [ROOT/'requirements.txt',ROOT/'requirements-publication.txt',ROOT/'README.md']
sources += list((ROOT/'kaggle').glob('*/*.ipynb'))
sources += list((ROOT/'kaggle').glob('*/*metadata.json'))
sources += [ROOT/'kaggle/README_RU.md',ROOT/'literature/q1_revision_notes_RU.md',
            ROOT/'literature/publication_metadata_audit.json',ROOT/'literature/author_code_audit.json',
            ROOT/'reports/mstls_parity.json']
sources += [ROOT/'literature/efficiency_notes_RU.md']
sources += [BASE/name for name in ('article_revised.md','supplement.md','information_theory.md',
                                 'publication_references.bib','numeric.csl')]
sources += [BASE/name for name in ('weak_efficiency_theory.md','efficiency_article_template.md',
            'efficiency_supplement_template.md','supporting_validation.md')]
sources += [BASE/'archive/pre_efficiency'/name for name in ('article_revised.md','supplement.md')]
sources += [ROOT/'reports/full_adaptive_protocol_RU.md']
if (ROOT/'reports/reviewer_results.json').exists():
    sources += [ROOT/'reports/reviewer_protocol_RU.md']
sources = [p for p in sources if p.name != 'inspect_revision_live_log.py' and p.parent.name not in ('efficiency_vdp_polish_0','efficiency_vdp_polish_1')]
manuscripts=[BASE/('article_revised'+suffix) for suffix in ('.md','.docx','.tex','.pdf')]
manuscripts += [BASE/('supplement'+suffix) for suffix in ('.md','.docx','.tex','.pdf')]
manuscripts += [BASE/('weak_efficiency_theory'+suffix) for suffix in ('.md','.docx','.tex','.pdf')]
manuscripts += [BASE/'publication_references.bib',BASE/'numeric.csl',BASE/'publication_build.json']
manuscripts += list((BASE/'figures/publication').glob('*'))+list((BASE/'submission').glob('*'))
data=[]
for stage in ('revision_fast','reweighting','geometry'):
    data += list((ROOT/f'results/kaggle_{stage}').glob('*.json'))
    data += list((ROOT/f'results/kaggle_{stage}').glob('*.npz'))
    data.append(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip')
data += list((ROOT/'results/kaggle_extended').glob('*.npz'))
data += list((ROOT/'results/kaggle_extended').glob('*.json'))
data += [ROOT/'kaggle/revision_inputs/input_manifest.json',
         ROOT/'reports/revision_results.json',ROOT/'reports/artifact_validation.json',
         ROOT/'reports/q1_revision_results_RU.md',ROOT/'reports/q1_revision_results_RU.docx',
         ROOT/'reports/research_progress_RU.md',ROOT/'reports/research_progress_RU.docx']
if (ROOT/'reports/full_adaptive_results.json').exists():
    for stage in ('adaptive_smoke','adaptive_classical',*[f'adaptive_wm_{i}' for i in range(4)]):
        folder=ROOT/f'results/kaggle_{stage}'
        data += [p for p in folder.rglob('*') if p.suffix in ('.json','.npz','.log','.txt')]
        data.append(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip')
    data += [ROOT/'reports/full_adaptive_results.json',ROOT/'reports/full_adaptive_validation.json',
             ROOT/'reports/full_adaptive_results_RU.md',ROOT/'reports/full_adaptive_results_RU.docx']
if (ROOT/'reports/reviewer_results.json').exists():
    data += list((ROOT/'results/kaggle_reviewer').glob('*.json'))
    data += list((ROOT/'results/kaggle_reviewer').glob('*.npz'))
    data += [ROOT/'results/kaggle_reviewer_download/wsindy_reviewer_outputs.zip',
             ROOT/'reports/reviewer_results.json',ROOT/'reports/reviewer_results_RU.md',
             ROOT/'reports/reviewer_results_RU.docx',ROOT/'reports/reviewer_validation.json']
archives=[]
if (ROOT/'reports/efficiency_results.json').exists():
    for stage in ('efficiency_scalar_v3','efficiency_geometry','efficiency_geometry_v2','efficiency_vdp_0','efficiency_vdp_1',
                  'efficiency_vdp_final_0','efficiency_vdp_final_1'):
        data += list((ROOT/f'results/kaggle_{stage}').glob('*.json'))
        data.append(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip')
    data += [ROOT/'reports'/name for name in ('efficiency_results.json','efficiency_results_RU.md',
               'efficiency_results_RU.docx','efficiency_validation.json')]
for name,files in [('submission_package',manuscripts),('reproducibility_code',sources),('reproducibility_data',data)]:
    path=DIST/(name+'.zip')
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for file in sorted(set(files)):
            if not file.is_file():raise FileNotFoundError(file)
            archive.write(file,file.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(path) as check:
        if check.testzip() is not None:raise ValueError('CRC failed')
    archives.append({'file':path.name,'bytes':path.stat().st_size,'files':len(set(files)),
                     'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(DIST/'manifest.json').write_text(json.dumps(archives,indent=2))
print(json.dumps(archives))
