# Information limits in weak-form identification

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
`v1.0.0` identifies the manuscript's public code/data version. No DOI has been
assigned. Authorship and journal submission are separate from this archive.
