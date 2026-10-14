# Independent-record extension retained outside the main narrative

## 7. Independent records with nuisance initial states

### 7.1. Conditional moment identities and independent instruments

Let independent records have means $x_{\theta,u_r}$, with a common coefficient vector and a separate unknown initial state for each record. Hermite-corrected weak moments obey $E\{m_r(\theta_0,\nu_0)\mid u_r\}=0$. Consequently the initial-state distribution can remain unspecified. At finite noise, a same-record pilot used in $D_r^\top\Omega_r^{-1}$ can correlate with the residual and spoil that identity. We use an explicit independent-instrument construction to state a feasible fixed-noise result.

For known OU covariance $C$, draw independent $\zeta_r\sim N(0,C)$ and set $Y_r^\pm=Y_r\pm\zeta_r$. These independent Gaussian views have covariance $2C$ [@leiner2025]. Instruments $L_r^\pm$ use their own view and independent training data; moments use the opposite view with its doubled correction variance. Define $B_r=(L_r^+G_r^-+L_r^-G_r^+)/2$ and $a_r=(L_r^+b_r^-+L_r^-b_r^+)/2$.

**Assumptions (A1)–(A5).** Trajectories and sensitivities are regular on a compact initial-state set; records and randomizations are independent with known positive covariance; the instruments obey the stated view separation and uniform moment bounds; the average derivative has a nonsingular limit and the centered record scores satisfy a Lindeberg CLT; the average instrumented quadrature defect is $o(n^{-1/2})$. The precise numbered conditions are in S16.1.

**Theorem S16.1 (heterogeneous-record weak estimating equations).** With known sensor covariance $C_{\nu_0}$ and under (A1)–(A5), conditional on the independent training data,

$$
\widehat\theta_n=\left(\sum_r B_r\right)^{-1}\sum_r a_r
\xrightarrow{P}\theta_0,\qquad
\sqrt n(\widehat\theta_n-\theta_0)
\Rightarrow N(0,B^{-1}V_\psi B^{-\top}),
$$

where $B=\lim n^{-1}\sum E B_r$ and $V_\psi=\lim n^{-1}\sum\operatorname{Var}(a_r-B_r\theta_0)$. The empirical covariance of the coefficient-dimensional record scores consistently estimates its sandwich. The proof uses $E(L_r^+m_r^-\mid u_r)=0$ and the symmetric identity, followed by the standard estimating-equation CLT [@hansen1982; @newey1994]; see S16.1.

**Corollary S16.2 (record-wise cross-fitting; assumptions and proof in S16.4).** Split independent records into a fixed number of balanced folds before observing the data. For each held-out fold, fit a coefficient pilot on its complement and construct the independent-view instruments from that pilot. If the pilots converge to a deterministic limit, their score maps satisfy the uniform continuity and envelope conditions (A6), and the quadrature defects satisfy (A5), the pooled estimator has the same first-order Gaussian limit and consistent record-score sandwich as the deterministic limiting instrument map. Every record enters the final estimating equation, with no external training sample.

The key mean identity holds for every instrument pilot: conditional on a fold complement, each held-out score is centered. Consistent complement pilots therefore change the normalized score sum by $o_P(1)$. Independence is used separately within each fold. Cross-fitting follows established sample-reuse methodology [@chernozhukov2018]; its role here is to preserve the corrected weak moment identity with heterogeneous initial states. The known covariance assumption remains in force. A record bootstrap can refit all complement pilots and retain every physical record's copies within its original fold; its finite-sample calibration requires numerical assessment.

### 7.2. Model weights and the product information reference

For deterministic clean-trajectory instruments, $D_r=E G_r$ and Proposition 3's full covariance $\Omega_r$ yield $L_r=D_r^\top\Omega_r^{-1}$. The oracle fixed-family covariance is $(\sum_rD_r^\top\Omega_r^{-1}D_r)^{-1}$. The feasible experiment estimates these instrument maps on independent training trajectories. Its sandwich estimates only a four-dimensional record-score covariance, even when hundreds of weak moments supply the within-record model weights.

The joint local information for $\theta,u_1,\ldots,u_n$ has separate nuisance blocks. Taking their Schur complement gives

$$
\mathcal I_{e,n}=\sum_r\mathcal I_e(u_r).
$$

Theorem 2 gives this sharp boundary-free reference record by record. For fixed regular records, dense test spans recover it in the sequential small-noise limit. Symmetric view separation preserves the leading original-record influence when pilot instruments converge regularly, because the auxiliary noises cancel at first order (S16.2). At finite amplitude the full Hermite covariance and random instruments determine a different sandwich. A simultaneous record-count, span-growth and noise-decay rate theorem remains open.

### 7.3. Covariance calibration and the common-trajectory case

The known-covariance construction isolates the initial-state nuisance question. Estimated covariance requires a consistent calibration and its influence contribution. Separate fixed-horizon residual fits estimate a new initial state for every record; their residual quadratic forms can retain incidental-parameter bias. Thus consistency of a residual covariance estimate needs a bias correction or an independently validated calibration scheme (S16.3). The pair-difference estimator in S14 applies to equal means only.

If all records share one trajectory and covariance is known, their Gaussian path average is sufficient for the mean parameter and has covariance $C/n$. Path averaging therefore provides the natural route to Corollary 1. For nonlinear polynomial moments it also suppresses higher Hermite orders faster than averaging separate moments. The common-trajectory fixed-family GMM result is retained as Lemma S14.1 and a calibration diagnostic in S14–S17.

For pure OU measurement error and a continuously differentiable deterministic trajectory, quadratic variation removes the finite-variation signal: $[Y]_T=2\sigma^2T/\tau$ coordinatewise. Thus infill realized quadratic variation divided by $2T$ consistently estimates $\sigma^2/\tau$ from one record, without fitting an initial state. This identifies the microergodic ratio; full covariance calibration for heterogeneous records remains open.


### 8.6. Heterogeneous records with rich model-weighted tests

Four hundred independent repetitions each contain 800 Van der Pol records, with initial states uniform on $[1.75,2.25]\times[-0.25,0.25]$, matching the acquisition range in S9.4. Counts 50, 100, 200, 400 and 800 are prefixes; 0.5% and 5% amplitudes share standard OU realizations. Marginal sensor variances use the fixed baseline trajectory calibration and $\tau=0.3$, with $T=20$ and $h=0.01$. Known covariance isolates the initial-state nuisance question. The initial states are hidden from feasible instruments.

A separate 1000-record training sample supplies one frozen coefficient pilot. Model trajectories at a 5-by-5 instrument grid (25 initial-state points) compute $D^\top\Omega^{-1}$ from the full Proposition 3 covariance. Bilinear interpolation at each fission view's noisy first observation supplies the initial-state instrument pilot; the weak moments use the independent opposite view with doubled correction variance, and the two directions are averaged. Thus finite-noise validity follows from instrument separation even with imperfect trajectory pilots. Inference is conditional on the frozen training sample: n counts inference records only, and the plotted information reference concerns those records rather than the total training-plus-inference budget.

The 120- and 480-window families are those in Section 8.2. Within-record weights are model based, while a four-dimensional empirical record-score covariance supplies the sandwich. A prespecified studentized bootstrap resamples complete record equations 999 times, with marginal bootstrap-t intervals and a joint quadratic critical value. An oracle 120-window reference receives generating coefficients and clean initial states; it is distinguished from both feasible constructions. No stacking comparator enters this experiment.

The training pilot is [0.99966, -0.98571, 0.97544, -0.96713], with relative error 0.02172. A separate 9-by-9 sensitivity-information grid (81 benchmark points) has maximum relative interpolation discrepancy 0.0007 on 12 separately seeded checks. The record-specific information matrices enter the summed Schur-complement reference; they are used only for benchmarking.

Table 6. Heterogeneous initial states at 5% noise, 400 repetitions per count. All RMSE columns are scaled by $\sqrt n$. The local reference uses the sum of initial-state-adjusted information; finite-noise estimates also have nonlinear and random-instrument terms. Joint coverage is for nominal 95% studentized record-bootstrap regions. Monte Carlo Wilson intervals and marginal coverage are in S18.

| n | Oracle 120 | Feasible 120 | Feasible 480 | Local reference | Joint 120 | Joint 480 |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 0.01762 | 0.02247 | 0.02227 | 0.01690 | 95.5% | 95.8% |
| 100 | 0.01707 | 0.02183 | 0.02179 | 0.01690 | 94.5% | 94.5% |
| 200 | 0.01735 | 0.02311 | 0.02293 | 0.01690 | 93.8% | 94.5% |
| 400 | 0.01697 | 0.02236 | 0.02219 | 0.01690 | 94.0% | 94.8% |
| 800 | 0.01726 | 0.02195 | 0.02170 | 0.01690 | 93.8% | 94.8% |

For 480 tests at 0.5% noise, the scaled RMSE is 0.001818 at 50 records and 0.001802 at 800, compared with local references 0.001690 and 0.001690. At 5% noise the corresponding scaled errors are 0.02227 and 0.02170; the finite-amplitude comparison is reported separately from local efficiency. At 50 records the 480-test joint bootstrap coverage is 95.8% with Wilson interval [93.3%, 97.3%]. The model-weighted construction avoids estimating a moment-dimensional empirical weight covariance from those 50 records.

![Figure 5. Different initial states: square-root-record-count-scaled RMSE against the summed nuisance-information reference, and nominal 95% joint record-bootstrap coverage. The oracle receives clean initial states; feasible instruments use independent views and a frozen training pilot. Left-panel whiskers resample independent Monte Carlo repetitions; right-panel whiskers are binomial Wilson intervals. The local Gaussian reference is a small-noise covariance scale. Training data are additional and fixed across repetitions; the reference counts inference records only.](figures/publication/heterogeneous_gmm.png)

The feasible 480-test construction remains about 7% above the local RMSE reference at 0.5% noise and 28% above it at 5% for 800 records. The near-reference oracle uses clean initial states and generating coefficients; its performance separates span approximation from feasible random-instrument accuracy. The coefficient training pilot is fixed at 5% noise, so this experiment also retains a pilot discrepancy as the inference amplitude decreases. It checks finite-noise consistency and accuracy of the declared construction, rather than the regular-pilot small-noise efficiency limit. At 0.5% noise and 800 records, joint bootstrap coverage is 92.25%, and its Wilson interval excludes 95%. Consequently small-sample calibration at 5% improves substantially, while uniform calibration over noise amplitudes and record counts remains unresolved.

### 8.7. Incidental initial states: local forecast and pilot

For Gaussian records with heterogeneous unknown initial states, the profile score can acquire a curvature term of order $\sigma^2$ [@neyman1948; @lancaster2000; @hahn2004; @arellano2007]. The local calculation in the companion curvature derivation gives $b(\sigma)=\sigma^2 a+o(\sigma^2)$. A full joint-sensitivity calculation adds the shared-parameter and random-initial-design terms $\sigma^2 c/n$. At 5% noise and $n=800$, these terms change the predicted relative bias from 0.00004795 to 0.00004811, an increase of 0.34%.

The independently seeded 50-repetition pilot has relative bias norm 0.0001031 and relative dispersion 0.0005823. The estimated bias norm is noisy and upward biased. A centered Monte Carlo bootstrap comparison of its four-dimensional mean with the local prediction nevertheless gives an exploratory discrepancy probability of approximately 0.008. The computed $1/n$ correction is too small to explain that discrepancy. Higher-noise remainder terms and Monte Carlo variation therefore need assessment before the numerical coefficient is interpreted as a validated finite-noise bias prediction.

Writing $A_W(\sigma)$ and $A_L(\sigma)$ for the relative variances of the square-root-record-count-scaled weak and likelihood estimates, a leading-order crossing forecast is

$$
n_*(\sigma)\simeq\frac{[A_W(\sigma)-A_L(\sigma)]_+}{\|b(\sigma)\|^2/\|\theta_0\|^2}.
$$

Deterministic integration of the full Hermite covariance over the declared initial-state distribution gives crossing counts of approximately 4580 and 2660 for oracle 120-test fission moments at 5% and 10% noise. These weights receive the generating coefficients and clean initial states; the likelihood variance and bias use their local expansions. The figures are order-of-magnitude forecasts for that oracle benchmark. They predict a crossing above 800 records, while the matched finite-noise comparison is required to establish empirical absence of a crossing over the checked counts. Conditional centering of weak moments supplies structural protection against incidental initial states; finite-sample accuracy depends on instruments, quadrature and estimator nonlinearity.
