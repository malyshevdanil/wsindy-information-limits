# Boundary-free weak estimating equations under correlated noise: efficiency with nuisance initial states

Statistical theory and a reproducible benchmark for boundary-free weak ODE estimation when initial states are nuisance parameters. The target format is **Electronic Journal of Statistics**. The main results are the exact adjoint influence class, its sharp initial-state Schur-complement information reference, the explicit optimal OU test and a conditional heterogeneous-record estimating-equation extension.

- [EJS PDF](manuscript/article_ejs.pdf), [IMS LaTeX](manuscript/article_ejs.tex), [Word](manuscript/article_revised.docx), [editable manuscript](manuscript/article_revised.md).
- [Supplement](manuscript/supplement.pdf), [weak-class proof](manuscript/weak_efficiency_theory.pdf), [heterogeneous-record proof and assumptions](manuscript/heterogeneous_records_theory.pdf).
- [Heterogeneous-record results](reports/heterogeneous_records_results.json), [record-bootstrap audit](reports/record_bootstrap_results.json), [matched 5% likelihood](reports/moderate_likelihood_results.json).
- [Full author-pipeline reproduction](reports/full_adaptive_results.json), [paired temporal weighting and WENDy](reports/efficiency_followup_results.json), [literature metadata](literature/publication_metadata_audit.json).

Rich model-weighted heterogeneous records use known OU covariance, independent instruments and a frozen independent training pilot. Their accuracy and coverage are conditional on that training sample; n counts inference records, not the full training-plus-inference budget. The local summed-information reference and fixed-noise sandwich are distinguished. Consistency of covariance estimated from heterogeneous incidental-state residuals requires additional calibration assumptions. Common-trajectory GMM is a supporting lemma and diagnostic: Gaussian path averaging is its natural full-record reference.

## Rebuild from recorded outputs

```powershell
python scripts/analyze_revision.py --stage revision_fast
python scripts/analyze_full_adaptive.py
python scripts/analyze_reviewer_validation.py
python scripts/prepare_efficiency_manuscript.py
python scripts/analyze_efficiency.py
python scripts/analyze_efficiency_followup.py
python -m scripts.analyze_information_followup
python scripts/prepare_efficiency_followup_manuscript.py
python scripts/analyze_independent_records.py
python scripts/prepare_independent_records_manuscript.py
python -m scripts.analyze_record_bootstrap
python scripts/prepare_nuisance_manuscript.py
python -m scripts.analyze_heterogeneous_records
python scripts/prepare_heterogeneous_results.py
python -m scripts.analyze_moderate_likelihood
python scripts/prepare_moderate_likelihood_results.py
python scripts/build_publication.py
python scripts/build_ejs.py
python scripts/validate_nuisance_runs.py
python scripts/validate_efficiency.py --publication
python scripts/validate_independent_records.py
python scripts/package_publication.py
python scripts/reproduce_efficiency_release.py
```

Numerical simulations run on Kaggle only. Local scripts analyze existing outputs, check source snapshots and checksums, and rebuild documents. Immutable notebooks contain frozen execution source and configs. Archives in `dist/` separate manuscripts, code and synthetic-data outputs with SHA256 manifests. Full time series can be regenerated from seeds; sufficient-statistic fixtures and all coefficient estimates are archived. Install `requirements-publication.txt`; the document builder uses the portable Tectonic binary from `scripts/install_publication_tools.py` on Windows. Other systems require the matching binary path. Word equations use native editable OMML. The IMS class and style are unmodified official support files; the explicitly requested `[ejs]` option is still supported, while the upstream package recommends `ejsv2` for new manuscripts.

Author: **Danil Romanovich Malyshev**, Independent researcher, Almaty, Kazakhstan; malyshevdanil04@gmail.com; [ORCID 0009-0007-0947-3866](https://orcid.org/0009-0007-0947-3866).

The paper has not been submitted. Human review, funding/conflict declarations and exclusive-submission confirmation remain administrative requirements before actual journal submission. DOI metadata is recorded only after a real Zenodo release is verified.
