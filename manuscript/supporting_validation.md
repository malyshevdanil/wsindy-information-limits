### S9.1. Full-record information and scalar likelihood

For $a=-1$, $x_0=2$, $T=4$, $\tau=0.4$ and $\sigma=0.1$, continuous Fisher information is 144.187. The chosen scalar weak function has linearized information efficiency 0.063; Section S9.8 evaluates optimal combinations within declared window families. On 1,000 paired paths, the finest-grid median absolute errors were 0.2271 (95% interval [0.2132, 0.2454]) for the weak ratio and 0.0561 ([0.0533, 0.0603]) for bounded Gaussian likelihood. The latter uses known covariance and initial condition. Its finite-sample bias and boundedness prevent interpreting its MSE as an unbiased Cramér–Rao efficiency measure.

![Figure S1. Scalar coefficient errors with median bootstrap intervals, and continuous Fisher information as duration increases.](figures/publication/information_and_likelihood.png)

The coarsest scalar grid already captures 99.937% of continuous-record information. Extending the decay record from $T=4$ to $8$ adds only 0.564%, whereas a second independent record doubles information. This is a quantitative acquisition comparison under the stated known-covariance model, rather than a universal cost-optimal allocation.

Initial-state annihilation raises the weak influence variance lower bound from 0.006935 to 0.008213, an increase of 18.4%. The chosen window's efficiency relative to the nuisance-adjusted information reference is 0.075, compared with 0.063 relative to known-initial-state information. This separates the weak class restriction from additional window inefficiency. Deterministic quadrature verifies both influence identities to $10^{-9}$; it adds no Monte Carlo paths.

### S9.2. Fixed-window sparse regression

The fixed-window comparison contains 24,576 fits: four systems, two noise regimes, two amplitudes, four grids, three methods and 128 paths per condition. Thirty fourth-power polynomial bump windows remain fixed as the grid changes. For decay, $T=4$ and $L=1.6$; for oscillators, $T=20$ and $L=3$. Correlation ratio is $\tau/L=0.1$. Relative MSTLS uses 50 thresholds between $10^{-4}$ and one, loss weight 0.8 and no ridge penalty; raw and known-variance corrected features share these settings. Settings and seed ranges were fixed before outputs; this is not a public preregistration. Initial states and libraries are specified in the supplement.

Table S1. Corrected relative MSTLS at 5% noise and $h=0.005$. Error and exact support are separate endpoints; 128 paths enter each row.

| System | Error model | Median relative error | 95% median interval | Exact support | 95% support interval |
|---|---|---:|---|---:|---|
| Exponential decay | OU | 0.0213 | [0.0160, 0.0248] | 128/128 | [0.971, 1.000] |
| Exponential decay | Independent | 0.0026 | [0.0021, 0.0034] | 128/128 | [0.971, 1.000] |
| Damped oscillator | OU | 0.0131 | [0.0114, 0.0152] | 128/128 | [0.971, 1.000] |
| Damped oscillator | Independent | 0.0013 | [0.0012, 0.0016] | 128/128 | [0.971, 1.000] |
| Van der Pol | OU | 0.0536 | [0.0402, 0.0658] | 81/128 | [0.547, 0.711] |
| Van der Pol | Independent | 0.0031 | [0.0028, 0.0034] | 126/128 | [0.945, 0.996] |
| Duffing | OU | 0.0804 | [0.0698, 0.0990] | 63/128 | [0.407, 0.578] |
| Duffing | Independent | 0.0061 | [0.0053, 0.0077] | 128/128 | [0.971, 1.000] |

In the Van der Pol OU condition, correct-support paths had median error 0.0340, whereas incorrect-support paths had median 0.7189. The aggregate median therefore combines parameter-estimation error with structure-selection error and cannot be read as an intrinsic information limit.

![Figure S2. Median coefficient errors with 95% path bootstrap intervals. Known-support fits separate parameter estimation from sparse selection.](figures/publication/sparse_refinement.png)

![Figure S3. Error distributions conditional on correct or incorrect selected support in the finest-grid Van der Pol OU condition. Counts accompany each distribution; bimodality is not presumed.](figures/publication/support_conditional.png)

Intervals resample independent paths and retain paired mesh comparisons. Conditions share random-number indices and are not pooled as independent experiments. Exact-support intervals use Wilson limits. No multiplicity-adjusted global support claim is made. No forecast endpoint is included in this validation.

For corrected MSTLS at 5% OU noise, 3 of four system-specific paired median-improvement intervals for $h=0.01$ to $0.005$ lie wholly inside the $\pm5\%$ band. This is a descriptive equivalence check of coefficient medians and does not guarantee support recovery.

### S9.3. Coefficient covariance and Monte Carlo calibration

Joint coefficient covariance was tested on 2,000 Van der Pol paths at each amplitude, on known support. Actual errors subtract the clean discrete fit; a guard retains designs whose smallest support singular value is at least half its clean value. Every simulated path passed: the guard did not filter these samples, so guarded summaries equal the full empirical summaries. This observation does not establish finite unconditional estimator moments or control unobserved near-singular designs. The guard defines the intended validation region; it is not evidence that the linearization is accurate.

Table S2. Relative Frobenius covariance discrepancies and paired relative RMS linearization errors on known support.

| Noise amplitude | Guarded paths | Linear covariance discrepancy | Guarded actual discrepancy | Paired linearization RMS |
|---|---:|---:|---:|---:|
| 0.1% | 2000/2000 | 0.040 | 0.040 | 0.025 |
| 1.0% | 2000/2000 | 0.040 | 0.041 | 0.245 |
| 5.0% | 2000/2000 | 0.040 | 0.294 | 1.156 |
| 20.0% | 2000/2000 | 0.040 | 0.481 | 2.591 |

The linear covariance discrepancy is close to its Monte Carlo RMS floor at every amplitude because the standardized linear experiment is paired across amplitudes. Actual covariance agrees at 0.1% and 1% noise but has substantially larger discrepancy at 5% and 20%. Paired pathwise linearization already has appreciable error at 1%. Thus the coefficient surrogate is a local small-noise description, and the residual covariance theorem alone cannot validate it at the moderate-noise benchmark setting.

The last column is $R=\{\sum_m\|e_m-\ell_m\|^2/\sum_m\|\ell_m\|^2\}^{1/2}$, where $e_m$ is the actual coefficient displacement from the clean fit and $\ell_m$ its first-order surrogate. It is not a relative error of the estimated coefficients and is not constrained to $[0,1]$. At 5% and 20% noise, omitted nonlinear terms and random-design coupling exceed the RMS magnitude of the first-order surrogate. A stable singular value alone does not suppress those terms.

![Figure S4. Joint coefficient covariance discrepancies, bootstrap intervals for guarded actual covariance, Gaussian Monte Carlo RMS references, and paired pathwise linearization error.](figures/publication/coefficient_covariance.png)

In the 4,000-path nonlinear residual experiment, higher orders contributed 20.3% of corrected trace at 50% noise. The covariance calibration places the observed Frobenius discrepancy against the fourth-moment Monte Carlo scale and an affine-null bootstrap. The null procedure retains empirical non-Gaussian shape and forces the analytical covariance; its range is approximate, and its exceedance fraction is not an exact p-value.

For the corrected 50% condition, observed matrix discrepancy was 0.064; estimated fourth-moment relative RMS floor was 0.063, and the approximate null 95% range was [0.053, 0.072]. The 6.4% discrepancy therefore must be read relative to this uncertainty, rather than as an isolated accuracy percentage.

![Figure S5. Observed covariance discrepancies and approximate covariance-null bootstrap ranges for three discrete amplitudes and two correction choices.](figures/publication/covariance_null.png)

### S9.4. Different initial conditions

A separate frozen experiment used 400 blocks with up to 16 independent Van der Pol recordings. Initial offsets were uniform on $[-0.25,0.25]^2$ about $(2,0)$; sensor noise covariance remained fixed. All methods supplied true support.

Table S3. Known-support aggregation with different initial conditions; 400 independent blocks per row.

| Method | M | Median relative error | 95% median interval |
|---|---:|---:|---|
| Stacked design | 1 | 0.03577 | [0.03305, 0.03888] |
| Stacked design | 16 | 0.03785 | [0.03684, 0.03857] |
| Mean weak matrices | 1 | 0.03577 | [0.03305, 0.03888] |
| Mean weak matrices | 16 | 0.00735 | [0.00693, 0.00773] |
| Mean sensor path | 1 | 0.03577 | [0.03305, 0.03888] |
| Mean sensor path | 16 | 0.10847 | [0.10523, 0.11193] |

Mean weak matrices improve with independent acquisition in this heterogeneous setting, whereas stacking retains a stable error and sensor-path averaging deteriorates. The latter estimates the dynamics of an averaged path that generally fails the original nonlinear equation. The evidence supports a distinction between aggregation operations, under the supplied support and excitation used here. Coefficient-level bias and dispersion are reported in Section S9.9; the near-equal single-record and stacked median errors alone do not identify the bias fraction.

![Figure S6. Independent recordings with different initial states. Mean weak matrices and sensor-path averaging define different estimating equations.](figures/publication/heterogeneous_recordings.png)

### S9.5. Covariance reweighting component

An additional comparison supplies known support and compares ordinary raw-feature least squares, the iid covariance IRLS step from WENDy Algorithm 2, and its known-OU covariance extension. It retains our fixed windows and uses ten iterations, relative stopping tolerance $10^{-6}$, known coordinate variances and scale-aware covariance relaxation $10^{-10}$. These choices make it a component experiment rather than a reproduction of complete WENDy.

Table S4. Van der Pol covariance-reweighting components: parameter one, initial state (2,0), T=20, 5% state-scaled noise, OU correlation time 0.3, 30 fixed windows of full width three, and 128 paths per condition. True support is supplied.

| Error model | Component | h | Median relative error | 95% median interval | Converged |
|---|---|---:|---:|---|---:|
| OU | OLS, known support | 0.04 | 0.0322 | [0.0289, 0.0378] | n/a |
| OU | OLS, known support | 0.005 | 0.0325 | [0.0289, 0.0363] | n/a |
| OU | WENDy iid IRLS core | 0.04 | 0.0185 | [0.0164, 0.0205] | 126/128 |
| OU | WENDy iid IRLS core | 0.005 | 0.0186 | [0.0165, 0.0203] | 127/128 |
| OU | OU IRLS core | 0.04 | 0.0184 | [0.0163, 0.0202] | 124/128 |
| OU | OU IRLS core | 0.005 | 0.0181 | [0.0161, 0.0198] | 121/128 |
| Independent | OLS, known support | 0.04 | 0.0080 | [0.0071, 0.0092] | n/a |
| Independent | OLS, known support | 0.005 | 0.0034 | [0.0031, 0.0037] | n/a |
| Independent | WENDy iid IRLS core | 0.04 | 0.0057 | [0.0047, 0.0065] | 128/128 |
| Independent | WENDy iid IRLS core | 0.005 | 0.0032 | [0.0030, 0.0034] | 128/128 |
| Independent | OU IRLS core | 0.04 | 0.0058 | [0.0051, 0.0065] | 128/128 |
| Independent | OU IRLS core | 0.005 | 0.0032 | [0.0029, 0.0035] | 128/128 |

At h=0.005 under OU error, the median decreases from 0.0325 for ordinary least squares to 0.0186 for iid IRLS and 0.0181 for OU IRLS. Reweighting therefore improves this component comparison, but its two IRLS medians are close and their displayed intervals overlap; the table alone does not establish superiority of the OU extension. OU errors remain nearly unchanged across the two meshes, whereas independent-noise errors decrease. Known covariance and supplied support distinguish this comparison from the complete adaptive pipelines in Section S9.6. All iteration-capped estimates enter the summaries. The convergence column reports the prescribed tolerance; it does not assess model correctness.





### S9.6. Complete adaptive author pipelines

A separate frozen comparison executes the complete pinned author WSINDy, WENDy and WmSINDy algorithms [@messenger2021; @bortz2023; @lopez2026]. WSINDy selects its test-function widths spectrally, places centers by the gradient, applies automatic SMAF smoothing and library scaling, and selects the MSTLS threshold. WENDy estimates noise, selects four adaptive radii and an SVD basis, and runs its full covariance iteration and stopping rules. WmSINDy jointly learns coefficients and noise with forward/backward RK4 and weak-residual losses, updates windows from the denoised data, and prunes and refits coefficients after each of eight cycles of 5000 Adam steps. Versions and compatibility adapters are detailed in Supplement S6.

The benchmark uses Van der Pol, five-percent state-scaled noise, fixed duration 20, correlation time 0.3, and 16 independent paths per condition. Meshes 0.04 and 0.01 are nested. All methods receive identical observations verified by SHA256; discovery methods share the nine nonconstant monomials of degrees one through three. Fixed-window MSTLS is rerun on this library. WENDy receives the prespecified correct four-term structure and is assessed as parameter estimation. Hermite MSTLS uses the known calibration variance, whereas author pipelines estimate or reconstruct noise from observations.

Table S5. Discovery pipelines at h=0.01; 16 independent records per row. Median intervals are path bootstraps. Support counts include all planned fits.

| Error model | Method | Median error [95% interval] | Exact support | Failed fits |
|---|---|---:|---:|---:|
| OU | Fixed raw MSTLS | 0.0402 [0.0300, 0.5657] | 11/16 | 0/16 |
| OU | Fixed Hermite MSTLS | 0.0420 [0.0311, 0.5656] | 11/16 | 0/16 |
| OU | Adaptive author WSINDy | 0.0430 [0.0354, 0.0536] | 14/16 | 0/16 |
| OU | Adaptive author WmSINDy | 0.0153 [0.0110, 0.0180] | 16/16 | 0/16 |
| Independent | Fixed raw MSTLS | 0.0043 [0.0036, 0.0059] | 16/16 | 0/16 |
| Independent | Fixed Hermite MSTLS | 0.0029 [0.0023, 0.0046] | 16/16 | 0/16 |
| Independent | Adaptive author WSINDy | 0.0086 [0.0070, 0.0118] | 16/16 | 0/16 |
| Independent | Adaptive author WmSINDy | 0.0059 [0.0036, 0.5396] | 10/16 | 0/16 |

Table S6. Parameter estimation with supplied support at h=0.01; 16 independent records per row.

| Error model | Method | Median error [95% interval] | Failed fits |
|---|---|---:|---:|
| OU | Adaptive author WENDy (known support) | 0.0151 [0.0120, 0.0212] | 0/16 |
| OU | Fixed raw, known support | 0.0300 [0.0224, 0.0367] | 0/16 |
| OU | Fixed Hermite, known support | 0.0312 [0.0235, 0.0384] | 0/16 |
| Independent | Adaptive author WENDy (known support) | 0.0029 [0.0023, 0.0035] | 0/16 |
| Independent | Fixed raw, known support | 0.0043 [0.0036, 0.0058] | 0/16 |
| Independent | Fixed Hermite, known support | 0.0029 [0.0023, 0.0042] | 0/16 |

On the finer OU grid, WmSINDy recovered support on 16/16 paths with median error 0.0153, whereas adaptive WSINDy recovered 14/16 with median 0.0430. The paired median error difference, WmSINDy minus WSINDy, was -0.0279 with a descriptive 95% bootstrap interval [-0.0378, -0.0204]. Under independent errors, WmSINDy recovered support on 10/16 paths, compared with 16/16 for fixed Hermite MSTLS. Its wide median interval reflects selection failures despite the full training budget. Thus the favorable OU result does not establish uniform superiority across observation regimes. These unadjusted comparisons are descriptive and limited to the frozen example.

The supplied-support comparison likewise concerns finite-sample accuracy: adaptive WENDy improves the OU coefficient median relative to these fixed-window fits, although its released noise model assumes independent observations. This empirical behavior does not establish calibrated uncertainty under OU error. Better window selection or reconstruction can improve accuracy within the finite information budget; it does not contradict the full-path infill bound.

All 448 planned fits are archived; 0 computational failures occurred. Error medians condition on completed fits. Wilson intervals, coarse-grid values, and paired mesh differences are retained in the reproducibility report. With 16 paths, support fractions have substantial binomial uncertainty. These conditions compare complete algorithms on one oscillator and do not establish a general ranking or an asymptotic plateau.

![Figure S7. Complete adaptive author pipelines and matched fixed-window comparisons. Error bars are 95% path bootstrap intervals with n=16 per cell. WENDy uses known support; the discovery methods do not.](figures/publication/full_adaptive.png)





### S9.7. Van der Pol information reference

The OU experiment in Section S9.6 has $x'=a v$ and $v'=b x+c v+d x^2v$, with four independently free active coefficients $(a,b,c,d)=(1,-1,1,-1)$. The initial state (2,0) is known. We hold the calibrated coordinate noise variances fixed when differentiating the mean, rather than re-estimating state-scaled variance at perturbed parameters. The calibration is exactly that of the benchmark fine grid. Coefficient and initial-state sensitivities solve $S'=JS+\partial_\theta f$ and $Z'=JZ$, with $S(0)=0$ and $Z(0)=I$. OU innovations give discrete information; the continuous norm in Section 4.1 gives its limit.

For four active coefficients, the relative RMS reference is $\{\operatorname{tr}(\mathcal I_{\theta\theta}^{-1})\}^{1/2}/\|\theta\|_2$. The weak initial-state-annihilation reference replaces information by its Schur complement. Both are local covariance bounds, rather than bounds on error medians or on every biased estimator.

Table S7. Van der Pol OU relative RMS information references for the exact Section S9.6 calibration. Four active coefficients are independently free. The second column assumes known initial state; the third also imposes weak initial-sensitivity annihilation.

| Observation | Known-initial reference | Weak influence reference |
|---|---:|---:|
| Grid h=0.04 | 0.013115 | 0.016929 |
| Grid h=0.01 | 0.013094 | 0.016906 |
| Continuous record | 0.013093 | 0.016904 |

For the same 16 OU records, WmSINDy has relative RMSE 0.01926 (descriptive 95% bootstrap interval [0.01324, 0.02576]) and WENDy 0.02154 ([0.01446, 0.02927]). Their RMSEs are respectively 1.47 and 1.65 times the known-initial local reference. This scale comparison is more informative than comparing their medians near 0.015 with a variance bound. It is not an efficiency certificate: sparse selection and nonlinear bias are present, the sample is small, and the regular locally unbiased assumption is unverified. WmSINDy also includes trajectory-prediction terms, so it is not restricted to boundary-free weak influences.

A one-parameter Van der Pol model fixes a=1 and b=-1 and ties c=mu, d=-mu. Its continuous standard-deviation reference for mu is 0.002492. That substantially stronger constraint changes the estimation problem and must not be used to judge four-coefficient discovery. Central finite differences check the four sensitivities with relative discrepancies below 1.6e-07; refining information quadrature from h=0.0025 to 0.00125 checks the reported reference.

### S9.8. Optimal efficiency within declared window spans

The scalar parameters are those of Section S9.1. We optimize the covariance-weighted linear combination in Section 5.4, including every cross-window covariance. Family F30 contains 30 fourth-power bump windows of full width 1.6; centers range uniformly from 0.805 to 3.195. Family F120 contains 60 centers for each of widths 1.6 and 3.0, with centers 0.005 inside the admissible support range. These are declared finite families, rather than an exhaustive search over all shapes and widths. The reference single window has width three and center two.

Table S8. Oracle optimal local influence variance and efficiency within each scalar window span. Information references are 144.187 for known initial state and 121.752 after initial-state projection. The table uses 8192 quadrature intervals and relative covariance-eigenvalue cutoff 1e-12.

| Family | Retained rank | Variance | Known-initial efficiency | Weak-class efficiency |
|---|---:|---:|---:|---:|
| Single window | 1 | 0.109776 | 0.0632 | 0.0748 |
| F30 (30 windows) | 30 | 0.010656 | 0.6508 | 0.7708 |
| F120 (120 windows) | 68 | 0.009339 | 0.7426 | 0.8795 |

Thus 0.063 describes one window, not the optimal efficiency of the weak family. F120 attains 0.743 of known-initial information in the retained span. Tightening the cutoff from 1e-10 to 1e-12 raises retained rank from 41 to 68; at 1e-10 its known-initial efficiency is 0.718. This sensitivity is disclosed rather than treated as an exact infinite-dimensional optimum. Doubling quadrature intervals changes every reported variance by less than 0.0004%. All optimal influences satisfy coefficient normalization and initial-state annihilation numerically. Their efficiencies remain below the initial-state-projected ceiling.

### S9.9. Bias and dispersion of aggregation

The coefficient-preserving execution replays exactly the 400 independent blocks in Table S3; it adds no independent samples. For normalized coefficient errors, the empirical identity states that squared RMSE equals squared bias plus squared dispersion. Table S9 reports the norm of the sample mean error as bias and the square root of the second term as dispersion. Bias here is an empirical estimate, not a proved population limit. We also subtract sample-mean sampling variance when estimating the squared-bias fraction; the adjustment is small for 400 blocks.

Table S9. Relative coefficient bias, dispersion and RMS error for heterogeneous Van der Pol aggregation; the same 400 blocks as Table S3. Values are dimensionless. Bias fraction uses the unadjusted empirical squared mean error.

| Method | M | Bias | Dispersion | RMSE | Squared-bias fraction |
|---|---:|---:|---:|---:|---:|
| Stacked design | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Stacked design | 16 | 0.03738 | 0.00835 | 0.03830 | 95.3% |
| Mean weak matrices | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Mean weak matrices | 16 | 0.00239 | 0.00893 | 0.00924 | 6.7% |
| Mean sensor path | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Mean sensor path | 16 | 0.11044 | 0.02914 | 0.11422 | 93.5% |

For stacking, dispersion decreases from 0.03229 at M=1 to 0.00835 at M=16, while bias changes from 0.02967 to 0.03738. The estimated bias contribution is 45.8% of single-record MSE and 95.3% at M=16; the sample-mean-adjusted fractions are 45.6% and 95.2%. The stacking floor is therefore dominated by systematic random-design error as acquisition increases. Whether the single-record error is mainly bias must be judged by this decomposition, not by the two median norms. Mean weak matrices reduce both errors under the stated excitation, whereas sensor-path averaging produces nonlinear averaging bias.




