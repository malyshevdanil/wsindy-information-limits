# Paired null calibration of Wald standard errors at 5% noise

This diagnostic compares three standard errors on the same 100 previously
evaluated Van der Pol null paths. The original 18-coefficient estimate,
480 time windows, nominal level 0.05, known initial state and known stationary
OU covariance remain fixed. The original nominal procedure rejected 23 of
100 null paths. The preregistered decision is compatibility of each rejection
count with 0.05 under its two-sided 95% Wilson interval. These previously
examined paths provide a diagnostic, rather than independent confirmation.

The protocol and reference observations are frozen in
`configs/wald_se_calibration.json` and `configs/wald_se_references.json`.
`reports/wald_se_protocol_decision.json` records the freeze time and hashes
before any new SE results. Every observed path is replayed against its
original input hash, full coefficient vector and nominal standard error.

## Random-design sandwich

Write the corrected weak residual as m(theta)=b-A theta. For this branch
the realized pilot covariance is held fixed and
W=(C+lambda I)^(-1), with the original ridge rule. Set H=A'WA. The target
influence in a model draw is

    psi* = e_j' H^(-1) A*' W (b* - d0 - A* theta0),

where theta0 is the null-imposed bounded core-five pilot and d0 is the
deterministic quadrature defect of that generating trajectory. A* is random:
the same OU noise enters its Hermite-corrected polynomial regressors and b*.
The standard error is the centered sample standard deviation of 999 psi*
draws. Thus the score meat includes design noise and covariance between the
design and residual. It is a model Monte Carlo estimate of the score meat,
not an exact analytic finite-sample variance formula. With a deterministic
design, the expression reduces to the usual fixed-weight GLS variance.

This is a local sandwich diagnostic at the observed bread. Variation in
fitted covariance weights is held fixed in this branch. A nonzero score
mean is recorded separately; changing the standard error does not remove
errors-in-variables bias or guarantee Gaussian studentization.

## Full-refit parametric bootstrap

Generate 399 independent stationary OU paths around the same null-imposed
pilot trajectory, with the original known covariance and initial state.
For each replicate, refit the 30-window bounded core-five pilot, integrate
its trajectory, rebuild the exact polynomial Hermite residual covariance,
apply the original pilot quadrature correction and refit all 18 final
coefficients. The bootstrap SE is the sample standard deviation of the
target estimates, using denominator B-1. The observed coefficient and
normal critical value remain unchanged. Neither bootstrap centering nor
an empirical rejection cutoff is used.

The generating model retains the four fitted background coefficients and
sets the target coefficient to zero. It shares the original core-five pilot
structure; robustness to an incorrectly restricted pilot is outside this
diagnostic. All bootstrap targets and timings are saved. Any integration or
numerical failure stops the shard; failed draws are not dropped.

## Computation and interpretation

The score and bootstrap seeds occupy distinct frozen ranges, disjoint from
the observed data. Five shards each contain 20 observed paths. A separate
one-path numerical smoke uses eight score draws and four bootstrap draws
solely to check execution and replay; it cannot assess calibration. Only
the complete 100-path sample is analyzed, without outcome-dependent changes
to counts, thresholds, pilot bounds or standard error scaling.

Wilson compatibility at 100 paths is a preliminary screening criterion.
If either revised SE passes, an independent null study and power study
remain necessary before claiming a calibrated inferential procedure.
If neither passes, the report retains the negative result and separates
score variability, pilot variability and possible centering errors without
adjusting this protocol retrospectively.
