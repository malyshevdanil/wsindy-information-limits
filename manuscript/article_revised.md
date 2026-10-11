---
title: "Fixed-horizon information limits and polynomial weak-residual covariance under correlated measurement error"
lang: en-US
---

## Abstract

Increasing the sampling frequency of a noisy trajectory improves numerical resolution but may leave parameter uncertainty essentially unchanged. We study deterministic ordinary differential equations observed through stationary Ornstein–Uhlenbeck measurement error with fixed physical correlation time. An explicit Cameron–Martin norm shows that smooth parameter-dependent trajectories induce equivalent continuous-path laws on a fixed observation interval. A quantitative two-point bound therefore excludes simultaneous consistency throughout a nontrivial parameter family, for every estimator based on increasingly dense samples. Dynamical sensitivities give finite continuous-record Fisher information. Boundary-free weak influences annihilate initial-state sensitivities, yielding a stricter local covariance bound even when the initial state is known. For polynomial vector fields, a finite Hermite expansion yields exact weak-residual covariance, including higher Gaussian orders and library–response cross-covariance; exponential recursions evaluate it without a dense time covariance matrix. Numerical comparisons examine full Gaussian likelihood, sparse weak regression, complete adaptive author pipelines, coefficient covariance and independent recordings. They separate support-selection failures from parameter uncertainty and locate the small-noise range of coefficient linearization. The analysis distinguishes denser sampling from longer observation and repeated acquisition: independent recordings add information, whereas longer observation requires persistent parameter sensitivity. These results provide an explicit statistical reference for interpreting sampling saturation in weak identification.

**Keywords:** weak-form identification; correlated measurement error; fixed-domain asymptotics; Fisher information; polynomial covariance; errors in variables.

## 1. Introduction

Data acquisition for dynamical identification involves three choices: how often to sample, how long to observe, and how many independent recordings to acquire. These choices have different statistical effects. Sampling more frequently improves numerical integration while repeatedly observing the same correlated error path. A longer recording may encounter new noise but little additional parameter sensitivity. Independent recordings can add information even when each individual trajectory remains statistically unresolved.

Sparse identification of nonlinear dynamics represents a vector field through a library and selects active coefficients [@brunton2016]. Integral formulations reduce sensitivity to numerical differentiation [@schaeffer2017], and weak formulations transfer differentiation to test functions [@messenger2021; @reinbold2020]. Their regression matrices still contain noisy states. Numerical resolution, feature bias, residual covariance and model selection interact, so reducing quadrature error does not by itself establish statistical consistency.

The distinction between fixed-domain and increasing-domain observations is classical in Gaussian process inference [@stein1999]. Ying established fixed-domain results for exponential covariance parameters [@ying1991], and Zhang used equivalence of Gaussian measures to explain inconsistent covariance estimation in the Matérn family [@zhang2004]. Those results concern covariance parameters. Here the covariance is common across parameters, and the parameter enters the deterministic mean through an ODE. The same measure-equivalence mechanism supplies an estimator-independent limit, with an explicit norm determined by dynamical sensitivities.

Continuous-time regression with correlated errors has an established optimal-design theory. Dette, Pepelyshev and Zhigljavsky derive continuous-record best linear unbiased estimators and efficient discrete designs for triangular covariance kernels [@dette2016]. Dynamical sensitivities define a local regression in this framework. Finite limiting precision and endpoint information are therefore precedents for our weak-estimation analysis.

The continuum-data theory of WSINDy studies uncorrelated measurement errors and library-dependent consistency [@messenger2025]. Fixed physical correlation time defines a different experiment: adjacent correlations approach one as the mesh is refined. Correlation induced by filtering independently sampled errors also requires distinguishing a physical filter scale from a fixed number of sample indices. Our analysis concerns a fixed covariance process observed at increasingly dense times.

We develop an explicit OU Gaussian-shift representation, a two-point lower bound, finite Fisher information and a weak influence bound that identifies loss of initial-state information. We then specialize Hermite identities to the coupled polynomial weak residual and derive its exact covariance and efficient evaluation. Finally, we contrast grid refinement, observation duration and independent acquisition. Classical Gaussian-measure, nuisance projection and polynomial identities are used with attribution; the contribution is their explicit dynamical formulation and numerical examination in this sampling regime.

## 2. Related estimation methods

WENDy treats weak parameter estimation as an errors-in-variables problem and iteratively updates a linearized residual covariance [@bortz2023]. It retains perturbations of both the feature matrix and the response. We use the same coupled residual and compute every Gaussian order for polynomial features under the specified OU model. First-order truncation recovers the structure motivating weak reweighting, with temporal covariance inserted explicitly. Residual covariance is distinct from the covariance of the nonlinear coefficient estimator.

Errors-in-variables effects have a substantial statistical literature [@fuller1987]. Gaussian correction of nonlinear functions of a noisy mean also precedes weak identification [@stefanski1989]. Hermite correction removes pointwise polynomial bias when marginal noise variance is known; it does not remove random design–response coupling. Generalized smoothing and profiling estimate ODE parameters through a different compromise between model fidelity and observation error [@ramsay2007]. The information bound applies to these alternatives when they observe the same path experiment.

WmSINDy combines weak regression with simultaneous dynamics and noise reconstruction [@lopez2026]. We distinguish the large fixed-window regression-component study from a separate paired execution of the complete adaptive author WSINDy, WENDy and WmSINDy pipelines. The latter retains noise estimation or reconstruction and adaptive test-function selection; WENDy supplies structure for parameter estimation. Ensemble-SINDy uses bootstrap aggregation to quantify feature inclusion [@fasel2022]. Resampling one recorded path and acquiring independent physical recordings are different operations: the latter changes the number of independent path experiments in the information calculation.

## 3. Observation model and weak equations

Let $x_\theta:[0,T]\to\mathbb R^d$ be a deterministic ODE trajectory. In the polynomial setting,

$$
x'_i=f_i(x)=\sum_{j=1}^{p}\xi^*_{ji}\vartheta_j(x),\qquad i=1,\ldots,d.
$$

The library contains all true terms. Observations are $Y(t)=x_\theta(t)+\varepsilon(t)$, with independent stationary OU coordinates [@uhlenbeck1930] and common-across-parameter covariance

$$
C_j(s,t)=\sigma_j^2\exp(-|s-t|/\tau_j),\qquad \sigma_j>0,\quad\tau_j>0.
$$

The continuous record is a statistical reference experiment; actual data are evaluations $y_n=Y(t_n)$ at mesh $h$. Initial conditions are known unless stated otherwise; unknown initial conditions would introduce nuisance sensitivities.

Fix $K$ test functions $\phi_k\in C_c^2(0,T)$. Their positions and full support widths $L$ remain fixed as $h$ decreases. Let $\mathcal Q_h$ be composite trapezoidal quadrature with weights $w_n$. Clean weak quantities are

$$
G^0_{kj}=\int_0^T\phi_k\vartheta_j(x)\,dt,\qquad
b^0_{ki}=-\int_0^T\phi'_kx_i\,dt,\qquad b^0=G^0\xi^*.
$$

Noisy matrices use corrected polynomial features $\widetilde\vartheta_j$, defined in Section 5:

$$
(G_h)_{kj}=\mathcal Q_h(\phi_k\widetilde\vartheta_j(y)),\qquad
(b_h)_{ki}=-\mathcal Q_h(\phi'_ky_i).
$$

Physical OU infill sets $\rho(h)=\exp(-h/\tau)$ with fixed $T,\tau$. An iid discrete control has independent error of fixed per-sample variance. Fixed-index AR(1) instead keeps $\rho<1$ constant, making the physical correlation time $-h/\log\rho\to0$. These models have different continuum limits.

## 4. An estimator-independent fixed-horizon limit

### 4.1. OU translation geometry

The stationary OU process satisfies $d\varepsilon_j=-\varepsilon_jdt/\tau_j+\sqrt{2\sigma_j^2/\tau_j}\,dW_j$, with an independent initial value $\varepsilon_j(0)\sim N(0,\sigma_j^2)$. Its Cameron–Martin space is $H^1([0,T];\mathbb R^d)$ as a set, with norm

$$
\|g\|_{\mathcal H}^2=\sum_j\left\{\frac{g_j(0)^2}{\sigma_j^2}+
\frac{\tau_j}{2\sigma_j^2}\int_0^T(g'_j+g_j/\tau_j)^2\,dt\right\}.
$$

Transforming a path into its initial value and innovation path $z(t)-z(0)+\tau_j^{-1}\int_0^tz(s)ds$ separates the normal initial observation from scaled Brownian motion. Translation by $g$ gives an initial shift and Brownian drift $g'+g/\tau_j$. The classical translation theorem [@cameron1944] gives the norm. Equivalently,

$$
\|g\|_{\mathcal H}^2=\sum_j\left\{\frac{\tau_j}{2\sigma_j^2}\int_0^T(g'_j)^2dt+
\frac1{2\sigma_j^2\tau_j}\int_0^Tg_j^2dt+
\frac{g_j(0)^2+g_j(T)^2}{2\sigma_j^2}\right\}.
$$

The endpoint terms encode stationary initial noise. Zero initial error would define a different experiment.

**Theorem 1 (fixed-horizon two-point bound).** Suppose $x_\theta\in H^1$ and all parameters share the stated positive OU covariance. For distinct $\theta,\eta$, put $\delta=\|\theta-\eta\|$ and $d_{\theta\eta}=\|x_\theta-x_\eta\|_{\mathcal H}$. Every estimator from any finite grid satisfies

$$
\max_{v\in\{\theta,\eta\}}P_v\{\|\widehat\theta_h-v\|\geq\delta/2\}
\geq\Phi(-d_{\theta\eta}/2)>0.
$$

No denser-grid estimator sequence is consistent in probability at every parameter of a nontrivial family. Moreover,

$$
\max_{v\in\{\theta,\eta\}}E_v\|\widehat\theta_h-v\|^2
\geq\frac{\delta^2}{4}\Phi(-d_{\theta\eta}/2).
$$

**Proof.** The mean difference belongs to the common Cameron–Martin space, making the path laws equivalent; this is the common-covariance case of Gaussian measure equivalence [@feldman1958]. Under the first law the log likelihood ratio of the second is normal with mean $-d_{\theta\eta}^2/2$ and variance $d_{\theta\eta}^2$. Their total variation is $2\Phi(d_{\theta\eta}/2)-1$. Sampling cannot increase total variation. A nearest-parameter test derived from an estimator can err only if that estimator is at least $\delta/2$ from the generating parameter. The sum of testing errors is at least one minus total variation, giving the probability bound. Multiplying the estimation probability by $\delta^2/4$ gives the risk bound. Both bounds hold at every mesh. $\square$

A constant estimator can be consistent at one parameter. The theorem excludes consistency throughout the family. The equivalence-based obstruction is established fixed-domain theory [@zhang2004]; its norm and quantitative specialization here concern ODE mean trajectories. The bound includes adaptive weak functions and full-path inference.

### 4.2. Fisher information and weak efficiency

If $\theta\mapsto x_\theta$ is differentiable in $\mathcal H$ and $s_r=\partial_{\theta_r}x_\theta$, continuous-record Fisher information is

$$
\mathcal I_{rs}(\theta)=\langle s_r,s_s\rangle_{\mathcal H}.
$$

Its inverse is the Cramér–Rao covariance bound for regular unbiased estimators, when nonsingular. Biased estimators require the corresponding bias derivative; finite-horizon maximum likelihood need not attain the unbiased bound.

For $x_a(t)=x_0e^{at}$ with known $x_0$, $s(t)=x_0te^{at}$ and

$$
\mathcal I_T(a)=\frac{\tau x_0^2}{2\sigma^2}
\int_0^Te^{2at}[1+(a+\tau^{-1})t]^2\,dt.
$$

On a uniform grid, exact innovations give

$$
\mathcal I_h(a)=\frac{s(0)^2}{\sigma^2}+
\frac{\sum_{n=1}^N[s(t_n)-\rho s(t_{n-1})]^2}{\sigma^2(1-\rho^2)}
\longrightarrow\mathcal I_T(a),\qquad\rho=e^{-h/\tau}.
$$

Nested grids give increasing information bounded by $\mathcal I_T$. The limit follows from a Riemann sum of squared sensitivity innovations. Statistical saturation therefore concerns the full experiment.

For $D_\phi=\int\phi x_a\ne0$, a scalar weak ratio has limit

$$
\widehat a_h=\frac{-\mathcal Q_h(\phi'Y)}{\mathcal Q_h(\phi Y)}
\longrightarrow a+\frac{A}{D_\phi+B},\qquad
A=\int(-\phi'-a\phi)\varepsilon,\quad B=\int\phi\varepsilon.
$$

OU path continuity gives almost sure convergence of the quadratures. The weights $\phi$ and $-\phi'-a\phi$ are independent: dependence would make a nonzero compactly supported function solve a homogeneous first-order equation. The OU kernel is strictly positive on nonzero weights, so $(A,B)$ is nondegenerate Gaussian. Positive conditional variance of $A$ given $B$ and positive denominator density at zero imply infinite absolute and squared error moments of the unguarded ratio. Its empirical summaries should use quantiles or explicit guards.

The influence $\ell=(-\phi'-a\phi)/D_\phi$ has Gaussian variance $V_\phi=\iint\ell(u)C(u,v)\ell(v)dudv$. Differentiating the weak identity gives $\int\ell s=1$. Cauchy–Schwarz in the covariance RKHS gives

$$
1\leq V_\phi\mathcal I_T(a),\qquad
\operatorname{Eff}_\phi=[\mathcal I_T(a)V_\phi]^{-1}\leq1.
$$

This efficiency concerns the Gaussian linearization; it is not the unguarded ratio's variance.

### 4.3. Optimal combination within a finite window family

For scalar exponential dynamics, let $w_k=-\phi'_k-a\phi_k$. Define the normalizing moments and residual covariance by

$$
D_k=\int\phi_k x_a\,dt,\qquad
C_{k\ell}=\iint w_k(s)C(s,t)w_\ell(t)\,ds\,dt.
$$

A combined Gaussian influence $\ell=\sum_k\alpha_kw_k$ is normalized by $\alpha^TD=1$. If the residual covariance matrix $C$ is positive definite, minimizing $\alpha^TC\alpha$ under this constraint gives

$$
\alpha_* = \frac{C^{-1}D}{D^TC^{-1}D},\qquad
V_*=(D^TC^{-1}D)^{-1},\qquad
\operatorname{Eff}_{\mathcal F}=\frac{D^TC^{-1}D}{\mathcal I_T(a)}.
$$

This follows either from a Lagrange multiplier or covariance Cauchy–Schwarz. With a singular matrix, the inverse is restricted to its retained positive eigenspace and the constraint must remain feasible. The weights use the generating coefficient and known covariance: this is an oracle local optimum within the specified span, rather than a demonstrated adaptive estimator or a global optimum over all compactly supported functions. Every member still annihilates initial-state sensitivity, so $V_*\geq\mathcal I_{a\mid u}^{-1}$. Section 7.8 reports declared families, spectral cutoffs and quadrature checks.

### 4.4. Initial-state information lost by weak influences

Let $S=\partial_\theta x$ and $Z=\partial_u x$ be trajectory sensitivities to active coefficients and initial state $u=x(0)$, respectively. Write $J=\partial_x f(x)$, so $Z'=JZ$. For fixed compactly supported test functions and a full-rank clean restricted design, the coefficient influence has the form $L_\varepsilon=\int_0^T\ell(t)\varepsilon(t)dt$. Integration by parts in the linearized weak equations gives

$$
\int_0^T\ell S\,dt=I,\qquad \int_0^T\ell Z\,dt=0.
$$

Indeed each unnormalized residual influence is $-\phi'_k I-\phi_kJ$ and therefore annihilates $Z$; multiplication by the clean regression inverse preserves this property. Coefficient sensitivities obey $S'=JS+\partial_\theta f$, whose weak integral is the clean design, giving the first identity.

Form the joint sensitivity information blocks $\mathcal I_{\theta\theta}$, $\mathcal I_{\theta u}$ and $\mathcal I_{uu}$ using the inner product of Section 4.2. Define

$$
\mathcal I_{\theta\mid u}=\mathcal I_{\theta\theta}
-\mathcal I_{\theta u}\mathcal I_{uu}^{-1}\mathcal I_{u\theta}.
$$

**Proposition 2 (weak influence covariance bound).** If these matrices are nonsingular, the Gaussian influence covariance satisfies

$$
\operatorname{Cov}(L_\varepsilon)\succeq\mathcal I_{\theta\mid u}^{-1}
\succeq\mathcal I_{\theta\theta}^{-1}.
$$

**Proof.** Project the coefficient score orthogonally away from initial-state scores in the Gaussian observation space. Its covariance is $\mathcal I_{\theta\mid u}$. The two influence identities give cross-covariance $I$ with the projected score. Subtracting its information-normalized version from $L_\varepsilon$ leaves covariance $\operatorname{Cov}(L_\varepsilon)-\mathcal I_{\theta\mid u}^{-1}\succeq0$. The second inequality follows from the Schur complement. $\square$

The projection is classical efficient-information geometry, consistent with continuous-regression precision theory [@dette2016]. The weak-specific restriction is annihilation of initial-state sensitivities. It holds even when the actual experiment knows $u$: these weak equations do not exploit that knowledge. For scalar decay, the nuisance direction is $q(t)=e^{at}$ and the stronger bound is $V_\phi\geq[\mathcal I_{aa}-\langle s,q\rangle_{\mathcal H}^2/\|q\|_{\mathcal H}^2]^{-1}$. Additional window inefficiency is distinct from this class restriction. The bound describes local influence covariance, without asserting efficiency of a nonlinear ratio or attainability within a finite window family.

## 5. Exact polynomial weak-residual covariance

### 5.1. Finite Hermite expansion

Assume a common $\tau$ for independent coordinates, and set $v_j=\sigma_j^2$. Define $H_0=1$, $H_1=u$, $H_{q+1}=uH_q-qvH_{q-1}$ and corrected monomials $\widetilde\vartheta_\beta(y)=\prod_jH_{\beta_j}(y_j;v_j)$. Their expectation given $x$ is $x^\beta$. This Gaussian debiasing construction [@stefanski1989] centers the polynomial residual at the true coefficients.

For continuum residual $r_{ki}=b_{ki}-(G\xi^*)_{ki}$, polynomial expansion gives

$$
r_{ki}=\sum_{1\leq|\alpha|\leq p_*}\int_0^Tc_{ki,\alpha}(t)H_\alpha(\varepsilon(t);v)dt,
$$

where $p_*\geq1$ is the relevant maximum degree, $H_\alpha=\prod_jH_{\alpha_j}$ and

$$
c_{ki,\alpha}=-\mathbf1_{\{\alpha=e_i\}}\phi'_k
-\phi_k\frac{\partial^\alpha f_i(x)}{\alpha!}.
$$

The response and feature perturbations share this coefficient and retain their cross-covariance.

**Proposition 3 (exact covariance).** Under these assumptions,

$$
\operatorname{Cov}(r_{ki},r_{\ell m})=
\sum_{1\leq|\alpha|\leq p_*}\alpha!\prod_jv_j^{\alpha_j}
\int_0^T\int_0^Tc_{ki,\alpha}(s)c_{\ell m,\alpha}(t)
e^{-|\alpha||s-t|/\tau}\,ds\,dt.
$$

**Proof.** The generating function $\exp(zu-vz^2/2)$ yields $E[H_q(U;v)H_r(V;v)]=\mathbf1_{q=r}q!v^q\rho^q$ for jointly Gaussian variables of variance $v$ and correlation $\rho$. This is the Hermite form of classical Gaussian moment identities [@isserlis1918]. Coordinate independence gives orthogonality between different multi-indices and factor $\alpha!\prod_jv_j^{\alpha_j}$. Substitute $\rho=e^{-|s-t|/\tau}$ and integrate the finite expansion. $\square$

Order $q$ has correlation time $\tau/q$. Higher orders need not be negligible at moderate amplitude. Different windows and equations remain correlated despite independent sensor coordinates. The covariance is exact for a supplied clean trajectory and true polynomial; plugging in estimates introduces another approximation.

### 5.2. Raw features and efficient evaluation

For correction variance $\widehat v$, define

$$
\bar f=\exp\left[\frac12\sum_j(v_j-\widehat v_j)\partial_j^2\right]f.
$$

The exponential terminates for polynomials. Using $\bar f$ in $c$ gives covariance, while mean residual is $\int\phi_k[f_i(x)-\bar f_i(x)]dt$. Raw monomials use $\widehat v=0$. Mean bias and covariance are separate outputs.

On a uniform grid, collect weighted coefficients into an $N\times Kd$ matrix $U_\alpha$. Each term is $\alpha!v^\alpha U_\alpha^TR_{q,h}U_\alpha$, where $(R_{q,h})_{nm}=\rho_q^{|n-m|}$ and $\rho_q=e^{-qh/\tau}$. For vector $u$,

$$
a_n=u_n+\rho_qa_{n-1},\qquad z_n=u_n+\rho_qz_{n+1},\qquad
R_{q,h}u=a+z-u,
$$

with zero exterior initialization. Each multi-index costs $O(N(Kd)^2)$ time and $O(NKd+(Kd)^2)$ memory, avoiding an $N\times N$ temporal matrix. An independent fixed-variance nugget has additional diagonal-in-time covariance; its linear weak contribution is $O(h)$. The exact mixed-noise formula is in the supplement.

### 5.3. Local coefficient uncertainty

Fix a supplied support with full-rank clean design. In a small-noise coupling $\varepsilon=\sigma Z$, polynomial smoothness and pseudoinverse continuity give, for each fixed path near the clean design,

$$
\widehat\xi-\xi^*=\sigma P_0r_1+O_Z(\sigma^2),\qquad
\Sigma_{\mathrm{lin}}=P_0\Sigma_{r,1}P_0^T.
$$

$P_0$ maps the joint residual ordered by window and equation to the supplied coefficient coordinates, retaining cross-equation covariance. $\Sigma_{r,1}$ contains Gaussian order one. Discrete validation subtracts the noise-free discrete fit to isolate quadrature error. This pathwise expansion defines a Gaussian surrogate, rather than an unconditional moment expansion. A singular-value guard specifies the conditional region for actual coefficient-covariance validation.

## 6. Sampling, duration and independent acquisition

For a constant unknown mean in OU error, continuous information is $\sigma^{-2}[1+T/(2\tau)]$. Thus $1+T/(2\tau)$ is its information-equivalent number of independent observations of variance $\sigma^2$. This gives a precise interpretation of the familiar $T/\tau$ scale in a specified model.

ODE coefficients depend on excitation. For exponential decay $a<0$, $\mathcal I_T(a)$ approaches a finite limit as $T\to\infty$, because late observations have vanishing sensitivity. At $a=0$,

$$
\mathcal I_T(0)=\frac{\tau x_0^2}{2\sigma^2}
\left[T+T^2/\tau+T^3/(3\tau^2)\right].
$$

A universal duration prescription based only on $T/\tau$ misses this dependence. Known covariance and specified sensitivities allow $\mathcal I_h/\mathcal I_T$ to quantify sampling information capture. It is a model-based reference, not a data-only sparse support diagnostic.

Independent recordings add information: identical experiments give $M\mathcal I_T$, and differing known initial conditions give the sum of their sensitivity information matrices. Shared calibration errors or drift would violate the independence assumption.

Weak aggregation depends on the estimating equation. For identical clean trajectories, corrected matrices satisfy $E[G_m]=G_h^0$ and $E[b_m]=b_h^0$. Averaging before regression converges to the clean finite-grid fit, followed by the true coefficient as quadrature error vanishes. Stacking instead has limit

$$
\xi_{\mathrm{stack}}=[E(G_m^TG_m)]^{-1}E(G_m^Tb_m),
$$

when invertible. Random design–response moments generally survive, as in errors-in-variables regression [@fuller1987]. With differing initial conditions and shared coefficients, weak-matrix averaging remains valid for independent identically distributed records if moments are finite and $E[G_m^0]$ has full rank on the chosen support, with quadrature controlled. Opposing excitation can cancel and invalidate this rank condition. Averaging sensor paths before evaluating nonlinear features is generally invalid because $f(E[x_m])\ne E[f(x_m)]$.

## 7. Numerical evidence

The numerical experiments use exponential decay, damped, Van der Pol and Duffing oscillators. Fixed windows and nested grids distinguish numerical resolution from independent-noise averaging. They evaluate Gaussian likelihood, sparse weak regression, polynomial residual and coefficient covariance, covariance reweighting, heterogeneous-record aggregation and complete adaptive pipelines. Known marginal variance is used where explicitly stated. Settings and provenance are detailed in the supplement.

<!-- REVISION_RESULTS_START -->

### 7.1. Full-record information and scalar likelihood

For $a=-1$, $x_0=2$, $T=4$, $\tau=0.4$ and $\sigma=0.1$, continuous Fisher information is 144.187. The chosen scalar weak function has linearized information efficiency 0.063; Section 7.8 evaluates optimal combinations within declared window families. On 1,000 paired paths, the finest-grid median absolute errors were 0.2271 (95% interval [0.2132, 0.2454]) for the weak ratio and 0.0561 ([0.0533, 0.0603]) for bounded Gaussian likelihood. The latter uses known covariance and initial condition. Its finite-sample bias and boundedness prevent interpreting its MSE as an unbiased Cramér–Rao efficiency measure.

![Figure 1. Scalar coefficient errors with median bootstrap intervals, and continuous Fisher information as duration increases.](figures/publication/information_and_likelihood.png)

The coarsest scalar grid already captures 99.937% of continuous-record information. Extending the decay record from $T=4$ to $8$ adds only 0.564%, whereas a second independent record doubles information. This is a quantitative acquisition comparison under the stated known-covariance model, rather than a universal cost-optimal allocation.

Initial-state annihilation raises the weak influence variance lower bound from 0.006935 to 0.008213, an increase of 18.4%. The chosen window's efficiency relative to the nuisance-adjusted information reference is 0.075, compared with 0.063 relative to known-initial-state information. This separates the weak class restriction from additional window inefficiency. Deterministic quadrature verifies both influence identities to $10^{-9}$; it adds no Monte Carlo paths.

### 7.2. Fixed-window sparse regression

The fixed-window comparison contains 24,576 fits: four systems, two noise regimes, two amplitudes, four grids, three methods and 128 paths per condition. Thirty fourth-power polynomial bump windows remain fixed as the grid changes. For decay, $T=4$ and $L=1.6$; for oscillators, $T=20$ and $L=3$. Correlation ratio is $\tau/L=0.1$. Relative MSTLS uses 50 thresholds between $10^{-4}$ and one, loss weight 0.8 and no ridge penalty; raw and known-variance corrected features share these settings. Settings and seed ranges were fixed before outputs; this is not a public preregistration. Initial states and libraries are specified in the supplement.

Table 1. Corrected relative MSTLS at 5% noise and $h=0.005$. Error and exact support are separate endpoints; 128 paths enter each row.

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

![Figure 2. Median coefficient errors with 95% path bootstrap intervals. Known-support fits separate parameter estimation from sparse selection.](figures/publication/sparse_refinement.png)

![Figure 3. Error distributions conditional on correct or incorrect selected support in the finest-grid Van der Pol OU condition. Counts accompany each distribution; bimodality is not presumed.](figures/publication/support_conditional.png)

Intervals resample independent paths and retain paired mesh comparisons. Conditions share random-number indices and are not pooled as independent experiments. Exact-support intervals use Wilson limits. No multiplicity-adjusted global support claim is made. No forecast endpoint is included in this validation.

For corrected MSTLS at 5% OU noise, 3 of four system-specific paired median-improvement intervals for $h=0.01$ to $0.005$ lie wholly inside the $\pm5\%$ band. This is a descriptive equivalence check of coefficient medians and does not guarantee support recovery.

### 7.3. Coefficient covariance and Monte Carlo calibration

Joint coefficient covariance was tested on 2,000 Van der Pol paths at each amplitude, on known support. Actual errors subtract the clean discrete fit; a guard retains designs whose smallest support singular value is at least half its clean value. Every simulated path passed: the guard did not filter these samples, so guarded summaries equal the full empirical summaries. This observation does not establish finite unconditional estimator moments or control unobserved near-singular designs. The guard defines the intended validation region; it is not evidence that the linearization is accurate.

Table 2. Relative Frobenius covariance discrepancies and paired relative RMS linearization errors on known support.

| Noise amplitude | Guarded paths | Linear covariance discrepancy | Guarded actual discrepancy | Paired linearization RMS |
|---|---:|---:|---:|---:|
| 0.1% | 2000/2000 | 0.040 | 0.040 | 0.025 |
| 1.0% | 2000/2000 | 0.040 | 0.041 | 0.245 |
| 5.0% | 2000/2000 | 0.040 | 0.294 | 1.156 |
| 20.0% | 2000/2000 | 0.040 | 0.481 | 2.591 |

The linear covariance discrepancy is close to its Monte Carlo RMS floor at every amplitude because the standardized linear experiment is paired across amplitudes. Actual covariance agrees at 0.1% and 1% noise but has substantially larger discrepancy at 5% and 20%. Paired pathwise linearization already has appreciable error at 1%. Thus the coefficient surrogate is a local small-noise description, and the residual covariance theorem alone cannot validate it at the moderate-noise benchmark setting.

The last column is $R=\{\sum_m\|e_m-\ell_m\|^2/\sum_m\|\ell_m\|^2\}^{1/2}$, where $e_m$ is the actual coefficient displacement from the clean fit and $\ell_m$ its first-order surrogate. It is not a relative error of the estimated coefficients and is not constrained to $[0,1]$. At 5% and 20% noise, omitted nonlinear terms and random-design coupling exceed the RMS magnitude of the first-order surrogate. A stable singular value alone does not suppress those terms.

![Figure 4. Joint coefficient covariance discrepancies, bootstrap intervals for guarded actual covariance, Gaussian Monte Carlo RMS references, and paired pathwise linearization error.](figures/publication/coefficient_covariance.png)

In the 4,000-path nonlinear residual experiment, higher orders contributed 20.3% of corrected trace at 50% noise. The covariance calibration places the observed Frobenius discrepancy against the fourth-moment Monte Carlo scale and an affine-null bootstrap. The null procedure retains empirical non-Gaussian shape and forces the analytical covariance; its range is approximate, and its exceedance fraction is not an exact p-value.

For the corrected 50% condition, observed matrix discrepancy was 0.064; estimated fourth-moment relative RMS floor was 0.063, and the approximate null 95% range was [0.053, 0.072]. The 6.4% discrepancy therefore must be read relative to this uncertainty, rather than as an isolated accuracy percentage.

![Figure 5. Observed covariance discrepancies and approximate covariance-null bootstrap ranges for three discrete amplitudes and two correction choices.](figures/publication/covariance_null.png)

### 7.4. Different initial conditions

A separate frozen experiment used 400 blocks with up to 16 independent Van der Pol recordings. Initial offsets were uniform on $[-0.25,0.25]^2$ about $(2,0)$; sensor noise covariance remained fixed. All methods supplied true support.

Table 3. Known-support aggregation with different initial conditions; 400 independent blocks per row.

| Method | M | Median relative error | 95% median interval |
|---|---:|---:|---|
| Stacked design | 1 | 0.03577 | [0.03305, 0.03888] |
| Stacked design | 16 | 0.03785 | [0.03684, 0.03857] |
| Mean weak matrices | 1 | 0.03577 | [0.03305, 0.03888] |
| Mean weak matrices | 16 | 0.00735 | [0.00693, 0.00773] |
| Mean sensor path | 1 | 0.03577 | [0.03305, 0.03888] |
| Mean sensor path | 16 | 0.10847 | [0.10523, 0.11193] |

Mean weak matrices improve with independent acquisition in this heterogeneous setting, whereas stacking retains a stable error and sensor-path averaging deteriorates. The latter estimates the dynamics of an averaged path that generally fails the original nonlinear equation. The evidence supports a distinction between aggregation operations, under the supplied support and excitation used here. Coefficient-level bias and dispersion are reported in Section 7.9; the near-equal single-record and stacked median errors alone do not identify the bias fraction.

![Figure 6. Independent recordings with different initial states. Mean weak matrices and sensor-path averaging define different estimating equations.](figures/publication/heterogeneous_recordings.png)

### 7.5. Covariance reweighting component

An additional comparison supplies known support and compares ordinary raw-feature least squares, the iid covariance IRLS step from WENDy Algorithm 2, and its known-OU covariance extension. It retains our fixed windows and uses ten iterations, relative stopping tolerance $10^{-6}$, known coordinate variances and scale-aware covariance relaxation $10^{-10}$. These choices make it a component experiment rather than a reproduction of complete WENDy.

Table 4. Van der Pol covariance-reweighting components: parameter one, initial state (2,0), T=20, 5% state-scaled noise, OU correlation time 0.3, 30 fixed windows of full width three, and 128 paths per condition. True support is supplied.

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

At h=0.005 under OU error, the median decreases from 0.0325 for ordinary least squares to 0.0186 for iid IRLS and 0.0181 for OU IRLS. Reweighting therefore improves this component comparison, but its two IRLS medians are close and their displayed intervals overlap; the table alone does not establish superiority of the OU extension. OU errors remain nearly unchanged across the two meshes, whereas independent-noise errors decrease. Known covariance and supplied support distinguish this comparison from the complete adaptive pipelines in Section 7.6. All iteration-capped estimates enter the summaries. The convergence column reports the prescribed tolerance; it does not assess model correctness.

<!-- REVISION_RESULTS_END -->

<!-- FULL_ADAPTIVE_START -->

### 7.6. Complete adaptive author pipelines

A separate frozen comparison executes the complete pinned author WSINDy, WENDy and WmSINDy algorithms [@messenger2021; @bortz2023; @lopez2026]. WSINDy selects its test-function widths spectrally, places centers by the gradient, applies automatic SMAF smoothing and library scaling, and selects the MSTLS threshold. WENDy estimates noise, selects four adaptive radii and an SVD basis, and runs its full covariance iteration and stopping rules. WmSINDy jointly learns coefficients and noise with forward/backward RK4 and weak-residual losses, updates windows from the denoised data, and prunes and refits coefficients after each of eight cycles of 5000 Adam steps. Versions and compatibility adapters are detailed in Supplement S6.

The benchmark uses Van der Pol, five-percent state-scaled noise, fixed duration 20, correlation time 0.3, and 16 independent paths per condition. Meshes 0.04 and 0.01 are nested. All methods receive identical observations verified by SHA256; discovery methods share the nine nonconstant monomials of degrees one through three. Fixed-window MSTLS is rerun on this library. WENDy receives the prespecified correct four-term structure and is assessed as parameter estimation. Hermite MSTLS uses the known calibration variance, whereas author pipelines estimate or reconstruct noise from observations.

Table 5. Discovery pipelines at h=0.01; 16 independent records per row. Median intervals are path bootstraps. Support counts include all planned fits.

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

Table 6. Parameter estimation with supplied support at h=0.01; 16 independent records per row.

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

![Figure 7. Complete adaptive author pipelines and matched fixed-window comparisons. Error bars are 95% path bootstrap intervals with n=16 per cell. WENDy uses known support; the discovery methods do not.](figures/publication/full_adaptive.png)

<!-- FULL_ADAPTIVE_END -->

<!-- REVIEWER_RESULTS_START -->

### 7.7. Van der Pol information reference

The OU experiment in Section 7.6 has $x'=a v$ and $v'=b x+c v+d x^2v$, with four independently free active coefficients $(a,b,c,d)=(1,-1,1,-1)$. The initial state (2,0) is known. We hold the calibrated coordinate noise variances fixed when differentiating the mean, rather than re-estimating state-scaled variance at perturbed parameters. The calibration is exactly that of the benchmark fine grid. Coefficient and initial-state sensitivities solve $S'=JS+\partial_\theta f$ and $Z'=JZ$, with $S(0)=0$ and $Z(0)=I$. OU innovations give discrete information; the continuous norm in Section 4.1 gives its limit.

For four active coefficients, the relative RMS reference is $\{\operatorname{tr}(\mathcal I_{\theta\theta}^{-1})\}^{1/2}/\|\theta\|_2$. The weak initial-state-annihilation reference replaces information by its Schur complement. Both are local covariance bounds, rather than bounds on error medians or on every biased estimator.

Table 7. Van der Pol OU relative RMS information references for the exact Section 7.6 calibration. Four active coefficients are independently free. The second column assumes known initial state; the third also imposes weak initial-sensitivity annihilation.

| Observation | Known-initial reference | Weak influence reference |
|---|---:|---:|
| Grid h=0.04 | 0.013115 | 0.016929 |
| Grid h=0.01 | 0.013094 | 0.016906 |
| Continuous record | 0.013093 | 0.016904 |

For the same 16 OU records, WmSINDy has relative RMSE 0.01926 (descriptive 95% bootstrap interval [0.01324, 0.02576]) and WENDy 0.02154 ([0.01446, 0.02927]). Their RMSEs are respectively 1.47 and 1.65 times the known-initial local reference. This scale comparison is more informative than comparing their medians near 0.015 with a variance bound. It is not an efficiency certificate: sparse selection and nonlinear bias are present, the sample is small, and the regular locally unbiased assumption is unverified. WmSINDy also includes trajectory-prediction terms, so it is not restricted to boundary-free weak influences.

A one-parameter Van der Pol model fixes a=1 and b=-1 and ties c=mu, d=-mu. Its continuous standard-deviation reference for mu is 0.002492. That substantially stronger constraint changes the estimation problem and must not be used to judge four-coefficient discovery. Central finite differences check the four sensitivities with relative discrepancies below 1.6e-07; refining information quadrature from h=0.0025 to 0.00125 checks the reported reference.

### 7.8. Optimal efficiency within declared window spans

The scalar parameters are those of Section 7.1. We optimize the covariance-weighted linear combination in Section 4.3, including every cross-window covariance. Family F30 contains 30 fourth-power bump windows of full width 1.6; centers range uniformly from 0.805 to 3.195. Family F120 contains 60 centers for each of widths 1.6 and 3.0, with centers 0.005 inside the admissible support range. These are declared finite families, rather than an exhaustive search over all shapes and widths. The reference single window has width three and center two.

Table 8. Oracle optimal local influence variance and efficiency within each scalar window span. Information references are 144.187 for known initial state and 121.752 after initial-state projection. The table uses 8192 quadrature intervals and relative covariance-eigenvalue cutoff 1e-12.

| Family | Retained rank | Variance | Known-initial efficiency | Weak-class efficiency |
|---|---:|---:|---:|---:|
| Single window | 1 | 0.109776 | 0.0632 | 0.0748 |
| F30 (30 windows) | 30 | 0.010656 | 0.6508 | 0.7708 |
| F120 (120 windows) | 68 | 0.009339 | 0.7426 | 0.8795 |

Thus 0.063 describes one window, not the optimal efficiency of the weak family. F120 attains 0.743 of known-initial information in the retained span. Tightening the cutoff from 1e-10 to 1e-12 raises retained rank from 41 to 68; at 1e-10 its known-initial efficiency is 0.718. This sensitivity is disclosed rather than treated as an exact infinite-dimensional optimum. Doubling quadrature intervals changes every reported variance by less than 0.0004%. All optimal influences satisfy coefficient normalization and initial-state annihilation numerically. Their efficiencies remain below the initial-state-projected ceiling.

### 7.9. Bias and dispersion of aggregation

The coefficient-preserving execution replays exactly the 400 independent blocks in Table 3; it adds no independent samples. For normalized coefficient errors, the empirical identity states that squared RMSE equals squared bias plus squared dispersion. Table 9 reports the norm of the sample mean error as bias and the square root of the second term as dispersion. Bias here is an empirical estimate, not a proved population limit. We also subtract sample-mean sampling variance when estimating the squared-bias fraction; the adjustment is small for 400 blocks.

Table 9. Relative coefficient bias, dispersion and RMS error for heterogeneous Van der Pol aggregation; the same 400 blocks as Table 3. Values are dimensionless. Bias fraction uses the unadjusted empirical squared mean error.

| Method | M | Bias | Dispersion | RMSE | Squared-bias fraction |
|---|---:|---:|---:|---:|---:|
| Stacked design | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Stacked design | 16 | 0.03738 | 0.00835 | 0.03830 | 95.3% |
| Mean weak matrices | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Mean weak matrices | 16 | 0.00239 | 0.00893 | 0.00924 | 6.7% |
| Mean sensor path | 1 | 0.02967 | 0.03229 | 0.04385 | 45.8% |
| Mean sensor path | 16 | 0.11044 | 0.02914 | 0.11422 | 93.5% |

For stacking, dispersion decreases from 0.03229 at M=1 to 0.00835 at M=16, while bias changes from 0.02967 to 0.03738. The estimated bias contribution is 45.8% of single-record MSE and 95.3% at M=16; the sample-mean-adjusted fractions are 45.6% and 95.2%. The stacking floor is therefore dominated by systematic random-design error as acquisition increases. Whether the single-record error is mainly bias must be judged by this decomposition, not by the two median norms. Mean weak matrices reduce both errors under the stated excitation, whereas sensor-path averaging produces nonlinear averaging bias.


<!-- REVIEWER_RESULTS_END -->

## 8. Discussion and limitations

Proposition 2 identifies a cost of eliminating boundary terms: the resulting influences discard initial-state sensitivity information. When the initial state is known, this restriction raises the local covariance reference before any particular window inefficiency is considered. Boundary-augmented moments or an additional trajectory-likelihood term could retain that information; their endpoint noise and its correlation with integrated noise would require explicit treatment. The present experiments quantify the restriction without testing those extensions.

The full-path bound changes the interpretation of sampling saturation. A stable error curve of one weak estimator describes that estimator; Theorem 1 establishes that no procedure can obtain simultaneous consistency from the specified single-path infill experiment. Algorithms can nevertheless have very different finite-sample accuracy. Likelihood, weighting and test-function design can improve estimation within the finite information available.

Higher Gaussian orders affect residual uncertainty after exact marginal correction. Their covariance is computable independently of sparse selection. Propagating order-one covariance through a clean restricted design is a separate local approximation: the independent validation finds moderate-noise coefficient discrepancies well above the Monte Carlo floor. Small-noise covariance agreement and exact nonlinear residual covariance therefore do not imply calibrated uncertainty for the full sparse estimator.

The analysis assumes additive stationary Gaussian error, positive covariance common across parameters and fully observed deterministic ODEs. Unknown covariance, dependent coordinates, partial observation and process noise require further work. General information theory allows unequal coordinate correlation times; the implemented polynomial covariance uses a common time scale. Acquisition costs are not estimated and no universal budget optimizer is asserted.

The large sparse comparison concerns a published regression component with fixed windows. The additional full adaptive author-pipeline benchmark uses one oscillator, two meshes and 16 paths per condition. It establishes executable comparisons under those conditions, with substantial support uncertainty, but cannot establish universal method rankings or an asymptotic plateau. WENDy and its fixed-window parameter comparators receive supplied structure; discovery results are reported separately. Oracle states and supplied support are verification devices. 

Acquisition decisions should compare discretization, missed sampling information, sensitivity over additional time and independence of additional records. For decay, repeating the experiment may be more informative than observing an increasingly long tail. Persistent excitation can instead make longer recordings useful. Application-specific covariance and sensitivities determine which regime applies.

## Declaration of generative AI and AI-assisted technologies

OpenAI Codex assisted literature searches, mathematical drafting, source development, experiment orchestration, numerical interpretation and manuscript preparation. Numerical results are derived from archived execution outputs, and figures are produced by plotting code. The tool is not an author. The human authors must review and edit the manuscript and accept responsibility for its content before submission; this statement does not assert that their review has already occurred.

## Data availability

Code, synthetic observations, numerical outputs and frozen execution notebooks are available in the public [reproducibility repository, version v1.0.1](https://github.com/malyshevdanil/wsindy-information-limits/tree/v1.0.1). The repository includes protocol configurations, source fingerprints, pinned author-code download specifications, compatibility adapters, postprocessing scripts and SHA256 manifests. Original author sources are downloaded from their pinned repositories rather than redistributed. No external experimental dataset is used. The tagged snapshot identifies the version supporting this article; a DOI has not been assigned.

## References
