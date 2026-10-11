---
title: "Supplementary material: information bounds, covariance validation and reproducibility"
lang: en-US
---

## S1. Continuous and discrete likelihood

All parameters share known positive stationary OU covariance. The innovation transformation separates an independent normal initial observation from a Brownian drift experiment. The translation theorem [@cameron1944] suffices for this common-covariance case; a changing-covariance criterion [@feldman1958] is unnecessary. The family-wide consistency obstruction is classical [@zhang2004]. Unknown initial conditions would enter the sensitivity matrix as nuisance parameters.

For the exponential model, set $q_a(t)=x_0(a+\tau^{-1})e^{at}$ and $\nu^2=2\sigma^2/\tau$. Continuous likelihood is

$$
\log L(a;Y)=\nu^{-2}\left[\int_0^Tq_a(t)(dY(t)+Y(t)dt/\tau)-\frac12\int_0^Tq_a(t)^2dt\right]+\text{constant}.
$$

The initial term is independent of $a$ when $x_0$ is fixed. The exact discrete negative log likelihood, with $e_n(a)=y_n-x_0e^{at_n}$, is

$$
-\log L_h(a)=\frac{e_0(a)^2}{2\sigma^2}+
\frac{\sum_{n=1}^N[e_n(a)-\rho e_{n-1}(a)]^2}{2\sigma^2(1-\rho^2)}+\text{constant}.
$$

The independent validation constrains optimization to $[-2.5,0.5]$ before simulation. Sixty-one parameter-grid points identify candidate local minima, each refined by bounded scalar optimization; endpoints are checked. This is a numerical bounded-likelihood estimate, without a proof that all minima between coarse points are found. It need not be unbiased or attain the Cramér–Rao bound. Its squared error is finite by boundedness. The weak ratio is summarized by quantiles rather than a population MSE.

The separate geometry execution performs deterministic sensitivity quadrature for the scalar example. It computes the joint information for $(a,x_0)$, its Schur complement and the weak influence's normalization and initial-state annihilation. Both identities are verified to $10^{-9}$ on the finest reference grid. Its config, source fingerprint and output archive are retained; this computation adds no stochastic paths and is not an additional Monte Carlo replication.

## S2. Mixed-noise covariance

Let each sample error contain OU variance $v_{c,j}$ and independent nugget variance $v_{w,j}$, with total $v_j=v_{c,j}+v_{w,j}$. For assumed correction $\widehat v$, use $\bar f=\exp[\frac12\sum_j(v_j-\widehat v_j)\partial_j^2]f$ in the residual coefficients. The exact discrete covariance is

$$
\Sigma_r=\sum_{\alpha\ne0}\alpha!\left[
v_c^\alpha U_\alpha^TR_{|\alpha|,h}U_\alpha+
(v^\alpha-v_c^\alpha)U_\alpha^TU_\alpha\right].
$$

Distinct times share only OU covariance; equal times contain total marginal variance. The response–library covariance is included in $U_\alpha$. Correctness fixtures compare this expression with Gaussian quadrature, including raw and miscalibrated features.

## S3. Covariance-estimation Monte Carlo error

For independent centered vectors $r_m$ with covariance $\Sigma$, use $S_n=n^{-1}\sum_mr_mr_m^T$ and theoretical centering. Independence gives

$$
E\|S_n-\Sigma\|_F^2=\frac{E\|r\|^4-\|\Sigma\|_F^2}{n}.
$$

For Gaussian vectors $E\|r\|^4=(\operatorname{tr}\Sigma)^2+2\|\Sigma\|_F^2$, making the relative RMS floor

$$
\frac{\sqrt{(\operatorname{tr}\Sigma)^2+\|\Sigma\|_F^2}}{\sqrt n\,\|\Sigma\|_F}.
$$

Polynomial residuals above order one are non-Gaussian. The protocol estimates their fourth moment and labels the Gaussian expression as a reference. An affine-null bootstrap centers archived residuals, whitens empirical covariance, recolors to the analytical target and resamples 250 multinomial count vectors with the original path count. Covariance eigenvalues are floored at $10^{-12}$ of the largest during whitening. This approximate null retains empirical shape; its exceedance fraction is not an exact test p-value.

Coefficient validation uses 2,000 paths at amplitudes $0.001,0.01,0.05,0.2$. Errors subtract the clean discrete fit. The residual surrogate uses the true coefficients; the exact finite-grid derivative about the clean fit would instead use its slightly quadrature-shifted coefficients. This remaining quadrature-scale distinction is retained rather than described as an exact finite-grid coefficient identity. A guard requires every support-restricted design's smallest singular value to exceed half the clean value. Actual covariance is conditional on this event, whereas the linearized Gaussian covariance is unguarded. Counts, bias and paired linearization errors accompany discrepancy values. The protocol does not presume finite unconditional nonlinear estimator moments.

## S4. Protocol, implementation and provenance

`configs/revision.json` was fixed before the corresponding numerical outputs; it is not a public preregistration. Numerical experiments execute in CPU Kaggle notebooks. Batched deterministic Van der Pol integration uses the recorded solver and tolerances; a correctness check compares batched and individual trajectories. Duplicate executions are not counted as independent evidence. Six synthetic residual arrays are used for covariance-error calibration. Notebook source fingerprints identify the executed implementation; local manuscript edits do not regenerate those execution snapshots.

Seed ranges start at 1,000,000 for sparse regression, 1,100,000 for scalar likelihood, 1,200,000 for coefficient validation, 1,300,000 for covariance bootstrap and 1,500,000 for heterogeneous records. Other experiment ranges are disjoint. Seed reuse within a condition pairs grids and methods. Conditions use common random numbers and cannot be pooled as independent experiments.

Four systems, two regimes, two amplitudes, four grids, three methods and 128 paths produce 24,576 fits. Thirty fourth-power bump functions have full support width 1.6 for decay and three for oscillators. All centers remain fixed. The correlation ratio is $\tau/L=0.1$. The independent protocol has no forecasting endpoint. Scalar likelihood has 1,000 paths on six nested grids. Heterogeneous acquisition uses 400 blocks with up to 16 recordings, initial offsets uniform on $[-0.25,0.25]^2$ around $(2,0)$ and fixed sensor noise amplitudes derived from the baseline path.

Relative MSTLS uses zero ridge penalty, unit rescaling, 50 logarithmic thresholds $10^{-4}$ to one and loss weight 0.8. Seven equation fixtures matched the pinned author regression code exactly; the repository commit and SHA-256 are recorded in the parity report. This component does not reproduce adaptive WSINDy or simultaneous noise reconstruction in WmSINDy [@messenger2021; @lopez2026]. Author source is inspected separately and excluded from the distribution package because no redistribution license was found in the inspected tree.

Every output retains source fingerprint, execution environment and completion status. Completed notebooks remain immutable snapshots. Seven correctness tests check exponential covariance application, integration by parts, OU sampling, Hermite correction, exact raw/miscalibrated/mixed covariance, clean recovery and conservative diagnostic refusal. Additional checks verify common-variance scaling of IRLS and batched versus individual ODE integration. Local postprocessing creates tables and figures from downloaded results. Median intervals resample independent paths; support intervals are Wilson intervals. No global multiplicity-adjusted comparison is claimed.

## S5. Aggregation and approximation diagnostics

Coefficient-level aggregation diagnostics use the same 400 blocks, random seeds, initial conditions and sensor errors as Table S3. A separate execution preserves the coefficient matrices without changing the estimators. It replays the existing experiment and is not an additional set of independent records. Error norms are checked against the original output before combining these results.

For coefficient error vectors $e_m$ normalized by $\|\xi^*\|_F$, define $B=\|\bar e\|$, $S^2=n^{-1}\sum_m\|e_m-\bar e\|^2$ and $R^2=n^{-1}\sum_m\|e_m\|^2$. The empirical identity $R^2=B^2+S^2$ separates squared bias of the sample mean from dispersion; finite-sample $B^2$ includes sampling error of that mean. An adjusted estimate subtracts $\operatorname{tr}(\widehat{\operatorname{Cov}}e)/n$ using the unbiased sample covariance. Medians alone do not identify this decomposition.

The relative RMS linearization diagnostic in Table S2 divides the RMS difference between actual and linearized coefficient displacements by the RMS linear displacement. It can exceed one. All 2000 paths pass the singular-value guard at each amplitude, so guarded and complete empirical samples coincide; this is not proof of finite population moments.

## S6. Complete adaptive author-pipeline reproduction

A separate protocol, `configs/full_adaptive.json`, was frozen before comparative outputs. It uses Van der Pol with parameter one, initial state (2,0), duration 20, noise amplitude five percent, OU correlation time 0.3 or independent errors, and meshes 0.04 and 0.01. Seeds 910000--910015 give 16 independent records per cell. All algorithms receive the same fine-grid error realization and its coarse restriction; SHA256 checks verify paired inputs. This seed range is separate from the earlier studies. Discovery libraries have nine nonconstant monomials of total degree one through three, matching the released WmSINDy example. Fixed-window discovery and supplied-support comparators are rerun on that library, rather than reusing ten-feature results.

Pinned author repositories are WSINDy_ODE commit `72d8d6e5be2b`, WENDy commit `65c52a893b1f`, and the WmSINDy Van der Pol example at commit `f1a781009a5c` [@messenger2021; @bortz2023; @lopez2026]. Files or repository ZIPs are fetched from those revisions and verified with the full revision identifiers and SHA256 values in `configs/author_pins.json`. Core author functions are executed without rewriting the statistical algorithms. Author sources are excluded from distributed archives; download instructions and our adapters are provided.

WSINDy uses FFT-corner widths and polynomial degrees, gradient-adapted placement with relax_AG=1, automatic SMAF smoothing initialized as in the released driver, scaling parameter two, 100 thresholds from 1e-4 to one, zero ridge penalty and loss weight 0.8. The adaptive placement is explicitly enabled relative to the released driver's uniform-placement default. GLS remains disabled, as in that driver. Selected widths, degrees, centers, smoothing-filter lengths and losses are saved.

WENDy uses the released set-and-forget settings: estimated noise, the exp(-9/(1-x²)) bump, automatic minimum radius and four dyadic multiples, SVD corner truncation, maximum 100 iterations, relative tolerance 1e-10, diagonal relaxation 1e-10, and normality checks after iteration ten with cutoff 1e-4. Its four-term structure is prespecified. Actual coefficients are withheld: the diagnostic true_vec is replaced by a unit vector and use_true=0. The returned normality value and iterates are archived. These are the released iid-noise equations; under physical OU error they are evaluated as a misspecified observation model, without claiming a full adaptive OU-WENDy extension.

The MATLAB functions execute under GNU Octave 8.4.0 with statistics and signal packages. Separate compatibility functions provide exact analytic polynomial Jacobians in the declared feature order, exact values and first derivatives of the specified bump with author normalization, and the two-input pagemtimes operation used in covariance construction. These replace Symbolic Toolbox and missing host operations. They do not change radius selection, SVD, covariance updates, regression or stopping logic. Bitwise agreement with a proprietary MATLAB runtime is not asserted.

WmSINDy executes the author Python initializer, RK4 predictions, weak-residual loss, gradient updates, adaptive window function, pruning and least-squares refits. Every completed fit performs eight outer cycles of 5000 Adam steps (40,000 total), with threshold 0.15, learning rate 0.001, epsilon 1e-7, q=1 and decay 0.9. Noise is initialized to zero. Legacy Adam is used with TensorFlow 2.20.0 / tf-keras 2.20.1. Slots for new coefficient variables are created outside traced functions, and the author loss functions are retraced after each window update so changed variable-container shapes are actually used. Adam is initialized independently for each recording; noise-variable slots and iteration count continue through its eight cycles. The released all-valid-convolution-row calculation is retained. Loss checkpoints, eight coefficient iterates, window widths, spectral corners, optimizer counts and the reconstructed noise signal are archived.

The compatibility control uses the author's Van der Pol example: parameter 0.5, initial state (-2,1), duration ten, 1000 linspace points, iid Gaussian amplitude 25 percent and seed one. Its nominal WmSINDy prediction step 0.01 is retained although the linspace spacing is 10/999. The benchmark instead uses its exact grid spacing. Clean trajectories are integrated with DOP853 at relative and absolute tolerances 1e-12 for this control. The single control checks execution and recovery and does not reproduce the published 50-run aggregate figure. Reported comparison coefficient errors use the Frobenius norm.

The comparison contains 448 planned fits: four discovery algorithms and three supplied-support parameter estimators, four observation cells and 16 records. Discovery success is exact support equality; supplied-support estimators receive no structure-discovery score. Failure counts remain in success denominators. Error medians and bootstrap intervals condition on completed fits; all failures and their logs are retained. Mesh differences are paired by recording. Wilson limits accompany support counts in `reports/full_adaptive_results.json`. Computational runtimes are execution diagnostics, not cross-language speed claims. With one oscillator and 16 paths per cell, these results do not justify universal rankings or fine precision of support probabilities.

## S7. Information, window optimization and coefficient replay

The configuration `configs/reviewer.json` specifies supplementary calculations for Sections S9.7–S9.9. The Van der Pol calibration uses the identical clean h=0.01 trajectory and coordinate standard deviations as the complete-pipeline comparison. Noise variances are held fixed under every parameter perturbation. Four independently free coefficients are ordered as velocity in equation one, then position, velocity and position-squared times velocity in equation two. The augmented DOP853 system integrates the two states, eight coefficient sensitivities and four initial-state sensitivities, with relative tolerance 1e-11 and absolute tolerance 1e-13. Stationary OU innovations give exact sampled-trajectory information from these numerical sensitivities; trapezoidal integration evaluates the continuous norm on grids from h=0.04 down to 0.00125. Central coefficient perturbations of size 1e-5 check sensitivities at 101 points. The covariance matrices, joint information, sensitivities and their derivatives are archived.

The scalar window optimum uses the complete cross-window residual covariance and the normalization constraint in Section 4.6. Covariance eigenvalues below 1e-10 or 1e-12 times the largest are discarded. The normalization vector is projected onto the retained eigenspace, where the constrained quadratic optimum is computed. The resulting influences are checked for unit coefficient sensitivity and zero initial-state sensitivity. Quadrature uses 4096 and 8192 intervals. The spectral-cutoff dependence of the enriched 120-window family is retained in S9.8 and numerical report; no global shape optimum is claimed. The calculation uses generating parameters and noise covariance, so the numbers concern an oracle local Gaussian reference.

The heterogeneous coefficient replay exactly reproduces all 3600 method/block/record-count errors of the 400-block aggregation experiment. It saves coefficients in addition to norms and does not supply new independent evidence. Bias and dispersion obey the empirical squared-error identity in S5. The bias-fraction adjustment subtracts the estimated sampling variance of the coefficient sample mean; neither adjusted nor unadjusted fractions are claimed as exact population quantities.

For the 16-path adaptive comparison, relative RMSE uses the complete coefficient matrices, including falsely selected terms, and the same Frobenius denominator as the error medians. Its 95% intervals resample independent records 5000 times. Comparing that RMSE with a local unbiased information reference indicates the scale of remaining error; it does not establish unbiasedness, regularity, calibrated uncertainty or attained efficiency.

