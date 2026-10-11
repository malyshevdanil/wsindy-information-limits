"""Prepare an allowlisted public reproducibility tree from verified archives."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'reports/public_release_v1'


def main():
    TARGET.mkdir(exist_ok=False)
    assert TARGET.resolve().is_relative_to(ROOT.resolve())
    for name in ('reproducibility_code.zip','reproducibility_data.zip'):
        with zipfile.ZipFile(ROOT/'dist'/name) as archive:
            assert archive.testzip() is None
            for entry in archive.infolist():
                assert (TARGET/entry.filename).resolve().is_relative_to(TARGET.resolve())
            archive.extractall(TARGET)
    # The public tree omits journal administration and personal author fields.
    public_readme='''# Information limits in weak-form identification

Code and synthetic data supporting the study of fixed-horizon OU measurement
error, polynomial weak-residual covariance, complete adaptive author pipelines,
Van der Pol Fisher information, finite-family window efficiency and aggregation
bias. This repository is a scientific reproducibility snapshot.

## Recorded results

- `reports/revision_results.json`: fixed-window, likelihood and covariance comparisons.
- `reports/full_adaptive_results.json`: all 448 full-pipeline comparative fits.
- `reports/reviewer_results.json`: information bounds, window combinations and bias/dispersion.
- `results/kaggle_*/`: frozen numerical outputs, inputs and original output archives.
- `kaggle/*/*.ipynb`: immutable execution-source snapshots.
- `configs/`: protocols and pinned author-source identifiers.

All empirical and deterministic numerical experiments execute on Kaggle.
Local commands below only postprocess recorded outputs and build documents.
Install `requirements.txt` and `requirements-publication.txt` first.

```console
python scripts/analyze_revision.py --stage revision_fast
python scripts/analyze_full_adaptive.py
python scripts/finalize_full_adaptive.py
python scripts/prepare_reviewer_manuscript.py
python scripts/analyze_reviewer_validation.py
python scripts/validate_reviewer.py
```

To rebuild Word/PDF on Windows, install the publication tools using
`python scripts/install_publication_tools.py`, then run
`python scripts/build_publication.py`. See the script for the Tectonic binary
path on another operating system. Coefficient medians are not CR efficiencies;
the information comparison uses relative RMS errors and states its assumptions.
The coefficient replay uses the same recorded seeds and is not new evidence.

## Reproducibility and provenance

Original third-party author code is not redistributed. `configs/author_pins.json`
provides repository revisions and SHA256 values; runners fetch those exact
sources. Our compatibility adapters are in `src/author_compat/`. Statistical
core functions, numerical tolerances and changed host operations are documented
in the manuscript supplement. Source fingerprints match frozen notebooks,
not subsequent local manuscript edits.

`release_manifest.json` provides SHA256 for every tracked artifact. The tag
`v1.0.1` identifies the manuscript's public code/data version. No DOI has been
assigned. Authorship and journal submission are separate from this archive.
'''
    (TARGET/'README.md').write_text(public_readme,encoding='utf-8')
    (TARGET/'.gitignore').write_text('__pycache__/\n*.pyc\n.pytest_cache/\ntools/\n.env\n.kaggle/\nreports/release_check*/\nreports/public_release*/\n')
    # Keep the scientific artifact bytes identical on every Git platform.
    (TARGET/'.gitattributes').write_text('* -text\n')
    # Figures are analysis outputs; include them to make the snapshot usable immediately.
    for p in (ROOT/'manuscript/figures/publication').glob('*'):
        dest=TARGET/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    # Never add the original Word plan or unpublished author declarations.
    entries=[]
    for p in sorted(TARGET.rglob('*')):
        if not p.is_file():continue
        assert not any(part in ('author_sources','author_repos','.kaggle') for part in p.parts)
        assert p.name not in ('access_token','kaggle.json','.env')
        data=p.read_bytes()
        assert not re.search(rb'(?:KGAT_[A-Za-z0-9]{20,}|ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,})',data)
        entries.append({'path':p.relative_to(TARGET).as_posix(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    (TARGET/'release_manifest.json').write_text(json.dumps(entries,indent=2))
    print(json.dumps({'prepared':True,'path':str(TARGET),'files':len(entries),'bytes':sum(r['bytes'] for r in entries),'credential_scan_passed':True}))


if __name__=='__main__':main()
