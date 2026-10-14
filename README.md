# The information cost of boundary-free weak estimation under correlated noise

Research manuscript and reproducibility archive prepared in Electronic Journal
of Statistics format. The contributions are exact weak adjoint-range geometry,
constructive full-library small-noise efficiency under a checkable jet condition,
and paired correlated-noise benchmarks with independently checked Wald size.
Known-initial trace adaptation is attributed to Brunel and coauthors.

OpenAI Codex assisted with literature search, mathematical derivations and proof
development, research code, and manuscript drafting and editing. The author's
confirmed verification and responsibility are recorded in the manuscript's
generative-AI declaration. AI tools are not authors. Prior commits and published
tags remain unchanged; historical redaction commands are disabled by AGENTS.md.

- [EJS manuscript](manuscript/article_ejs.pdf), [IMS LaTeX](manuscript/article_ejs.tex),
  [Word](manuscript/article_revised.docx), [editable source](manuscript/article_revised.md).
- [Submission supplement](manuscript/supplement.pdf): core proofs and numerical audits.
- [Separate research archive](manuscript/research_extensions_archive.pdf): earlier
  independent-record GMM, fission, cross-fitting, incidental-state studies and
  historical comparison protocols. Their original section labels are preserved.
- [Constructive proof](manuscript/self_instrumented_gls_theory.md) and
  [matched finite-span covariance](reports/guarded_geometry_comparison.json).
- Original-byte outputs and SHA-256 manifests are under `data/endpoint_followups`,
  `data/wald_se_calibration`, `data/independent_validation` and
  `data/guarded_geometry_reference`. Immutable notebooks are under `kaggle`.

Theorems concern known positive OU covariance, deterministic polynomial ODEs
and sequential small-noise/dense-span limits. The full-library pilot uses
smoothed observations. The older model-pilot GLS and the targeted Wald
procedure have different theoretical scopes; see supplement Table S8.
The independent small-noise size tolerance is met, while the full power
tolerance check fails at half a reference unit. Moderate-noise bootstrap
validation remains incomplete after a recorded compute-budget decision.
Unsubmitted frozen notebooks are protocol artifacts, not completed experiments.

## Reproduce and rebuild

Install `requirements-publication.txt`. Research ODE integrations and new
observation simulations run on Kaggle. Local analysis verifies archived
coefficient errors, frozen source and byte manifests. The portable Tectonic
binary is installed with `scripts/install_publication_tools.py` on Windows.
Document equations are editable Word OMML; EJS uses the official IMS class.

```powershell
python scripts/analyze_independent_validation.py
python scripts/analyze_guarded_geometry_reference.py
python scripts/verify_vdp_jets.py
python -m unittest discover -s tests -p test_calibration_repairs.py -v
python -m unittest discover -s tests -p test_self_instrumented_algebra.py -v
python scripts/build_publication.py
python scripts/build_ejs.py
python scripts/validate_endpoint_followups.py
```

`build_publication.py` builds the current main paper, core supplement and
separate archive; historical preparation scripts should not overwrite the
current manuscript. Input trajectories are regenerated from frozen seed
protocols. Pinned download specifications identify original author sources.
The preceding v1.2.0 archive remains at
[its Zenodo DOI](https://doi.org/10.5281/zenodo.23313935); the current version
has separate metadata in `configs/public_archive.json`.
