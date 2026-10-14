# Literature filter: endpoint information and support detection

Checked 2026-10-13. Primary sources only; distinctions below are provisional
novelty assessments, not claims that a comprehensive search proves priority.

## Known initial values: closest prior work

Brunel, Clairon and d'Alche-Buc (2014), *Parametric estimation of ordinary
differential equations with orthogonality conditions*, JASA 109, 173–185,
DOI [10.1080/01621459.2013.841583](https://doi.org/10.1080/01621459.2013.841583).
Author manuscript: https://arxiv.org/pdf/1410.7566.
Section 2.3, manuscript pp. 12–13 (PDF indices 11–12), explicitly identifies
information loss from endpoint-vanishing tests. For known initial values it
allows a nonzero initial test trace while retaining a zero terminal trace.
Section 5.2.1 also constrains the spline proxy to the known initial value.
**Decision:** endpoint-enabled weak moments and their motivation are established
prior work. Neither should be claimed as an invention here. The proposed
contribution is a sharp correlated-noise information comparison and an explicit
OU realization, conditional on successful numerical validation.

Tran and Bortz, *Weak Form Scientific Machine Learning: Test Function
Construction for System Identification*, SISC 48(4), C890–C915, published online
2026-07-28, DOI [10.1137/25M1776020](https://doi.org/10.1137/25M1776020).
Author manuscript: https://arxiv.org/pdf/2507.03206.
Sections 1.2–3 construct compactly supported tests and select local support by
estimating quadrature error. This addresses implementable test selection,
including errors-in-variables weighting, rather than the known-initial OU
endpoint-score recovery investigated here. Its quadrature analysis is directly
relevant: nonzero endpoints require an independent clean-defect check.

Messenger and Bortz (2021), *Weak SINDy: Galerkin-Based Data-Driven Model
Selection*, DOI [10.1137/20M1343166](https://doi.org/10.1137/20M1343166),
https://pmc.ncbi.nlm.nih.gov/articles/PMC10795802/;
Bortz et al. (2023), *Direct Estimation of Parameters in ODE Models Using WENDy*,
https://pmc.ncbi.nlm.nih.gov/articles/PMC10556152/.
These establish weak regression, noise-aware covariance weighting and adaptive
compact tests. The candidate is a known-support information benchmark, not a
replacement reproduction of either complete method.

## Support-detection direction: filter, not an experiment result

Messenger and Bortz, *Asymptotic consistency of the WSINDy algorithm in the
limit of continuum data*, IMA J. Numerical Analysis 45(6), 3264–3312 (2025),
DOI [10.1093/imanum/drae086](https://doi.org/10.1093/imanum/drae086),
author manuscript https://arxiv.org/abs/2211.16000.
Theorems 3.1 and 4.2 address continuum moment behavior and recovery conditions.
Their continuum-noise experiment must be distinguished from an OU path with
fixed physical correlation time; denser observations of the latter retain a
finite information limit.

Chen, Fasel and Bizyaeva, *Fourier Weak SINDy: Spectral Test Function Selection
for Robust Model Identification*, https://arxiv.org/abs/2604.20141 (2026 preprint),
is another recent test-selection comparator. Its existence prevents treating
spectral test selection as new. A full comparison remains pending.

A September 2026 author preprint, Tala, Shah and Patel,
*A recovery frontier for weak-form SINDy: exact support recovery and automatic
threshold selection with PACE*, https://doi.org/10.2139/ssrn.7447203,
already studies empirical recovery frontiers and threshold selection. Read its
full methods before claiming a new phase-transition experiment. A Gaussian
fixed-horizon impossibility reference could complement such recovery studies.

For two Gaussian path means with a common covariance, the optimal equal-prior
testing error is Phi(-||x1-x0||_H/2). Discrete observation cannot improve on the
full-path experiment. This is an established Gaussian testing identity; the
new application would require actual model-specific distances. If coefficients
other than the candidate term are unknown, the local reference profiles both
those coefficients and initial states: it is NOT simply the diagonal entry of
I_e. The relevant scalar information is 1/(I_e^{-1})[j,j].
Local delta_min = 2 Phi^{-1}(1-alpha) sqrt((I_e^{-1})[j,j]) is a local benchmark,
not an exact finite-amplitude guarantee about support recovery.
Keep the measurement covariance fixed across the alternatives; recalibrating
noise by each alternative's trajectory variance changes the testing model.

For the proposed Van der Pol x^2 v experiment the null coefficient d=0 gives
an anti-damped linear oscillator when the remaining coefficients are held at
(1,-1,1). Its trajectory over T=20 differs qualitatively from the d=-1 limit
cycle. Extrapolating a sensitivity at d=-1 all the way to a zero-coefficient
support test is therefore inappropriate. The screening calculation must first
evaluate exact null/alternative path distances on a prespecified coefficient
grid. Only then can a useful recovery simulation and its threshold criterion be
frozen. A finite-noise testing bound for all nonlinear weak statistics also
requires more than the first-order weak influence-space theorem.

## Second paper

Leiner et al., *Data fission: splitting a single data point*,
https://arxiv.org/abs/2112.11079, provides the established selection/inference
split mechanism. An ODE application also needs the selected-model target,
misspecification treatment and known-covariance scope to be stated. Independent
fission pieces alone do not establish finite-sample confidence coverage for a
nonlinear fitted weak estimator. Data thinning and recent conditional-inference
reviews remain to be examined before that separate project is frozen.
