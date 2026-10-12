# Nuisance initial-state revision

The common-trajectory Gaussian path average is sufficient at known covariance. The old averaged polynomial moment experiment is retained as a calibration diagnostic (Lemma 2), not an efficiency contribution.

Heterogeneous initial states introduce record-specific nuisance parameters. The summed Schur-complement reference is an exact local block-information calculation. Fixed-noise mean identities hold conditional on each initial state, but same-record random instruments need not preserve them. A symmetric known-covariance Gaussian fission construction supplies independent instruments; this is credited to Leiner et al., JASA 120(549), 135?146 (2025), DOI 10.1080/01621459.2023.2270748. The full author preprint https://arxiv.org/html/2112.11079v9 was inspected at Sections 2.2 and B.1.1. Online publication was 2023; the final volume is 2025.

Recentered studentized GMM bootstrap is credited to Hall and Horowitz (1996), Econometrica 64(4), 891?916, DOI 10.2307/2171849. Metadata and publisher entry were verified; the implementation rationale was checked in Horowitz?s author-written Handbook chapter, https://users.ssc.wisc.edu/~behansen/718/Horowitz%20Handbook.pdf, Section 3.7. The original Hall?Horowitz full proof was not claimed as read.

The common-trajectory bootstrap overcovers severely at n=50 (99.9% joint, median volume ratio 248) and approaches nominal coverage at n>=200. This is reported honestly; no critical-value tuning to Monte Carlo coverage is used. Heterogeneous model weights use Prop3, avoiding empirical high-dimensional covariance weighting; a four-dimensional record-score sandwich and record bootstrap remain feasible at small n.

Known-covariance fission does not prove consistency of incidental-state residual covariance estimation. An external consistent calibration and its influence, or a demonstrated residual bias correction, are required for an estimated-covariance extension. No claim of finite-noise global likelihood efficiency or joint n/span/noise rates is made.

Official EJS manuscript support was downloaded from https://github.com/vtex-soft/texsupport.ims_cosponsored-ejs. Class files are unmodified. The user requested option ejs, which remains supported; the upstream 2025 update recommends ejsv2.
