---
title: "Efficient weak-form estimation with correlated noise: the cost of boundary-free tests"
author: "Danil Romanovich Malyshev"
lang: en-US
---

Independent researcher, Almaty, Kazakhstan

Corresponding author: malyshevdanil04@gmail.com

## Abstract

Weak-form dynamical identification avoids differentiating noisy measurements, but removing boundary terms also discards information. For deterministic ODEs observed in stationary Ornstein-Uhlenbeck measurement error, we characterize the closure of boundary-free adjoint influences: it is exactly the covariance-dual subspace orthogonal to initial-state sensitivities. This gives a sharp local covariance bound equal to the inverse initial-state Schur complement, even when the initial state is known. We construct its optimal influence and an explicit generalized test function with endpoint jumps. Compact tests attain the bound as an infimum. Correct-covariance weak generalized least squares approaches this bound for dense test spans, in a sequential small-noise limit with controlled quadrature and regular pilots. A finite Hermite expansion supplies exact coupled polynomial residual weights. A fixed-horizon two-point bound separates this efficiency question from consistency under denser sampling. Experiments test a smoothed scalar optimal function with a fitted pilot and nested 30-, 120- and 480-window Van der Pol families, comparing physical OU and iid weights. Small-noise accuracy is compared with the information reference; finite-noise bias and dispersion are reported separately. The contribution identifies which information weak equations can retain and how covariance weighting and test-space approximation govern its recovery.

**Keywords:** weak-form identification; correlated measurement error; fixed-domain asymptotics; Fisher information; polynomial covariance; errors in variables.

## 1. Introduction

Boundary-vanishing weak equations can lose information even when an initial state is known, as recognized by orthogonality-condition estimation [@brunel2014] and recent weak identifiability work [@heitzman2026]. We quantify this effect: the closed class of boundary-free weak influences has exactly the local information available when the initial state is treated as a nuisance parameter in the full path likelihood. For an unknown initial state, this restriction causes no additional first-order information loss; for a known initial state, the lost information is the initial-sensitivity projection, and the retained information is its Schur complement. This is a local information statement, not an identity between finite-noise estimators.

Data acquisition for dynamical identification involves three choices: how often to sample, how long to observe, and how many independent recordings to acquire. These choices have different statistical effects. Sampling more frequently improves numerical integration while repeatedly observing the same correlated error path. A longer recording may encounter new noise but little additional parameter sensitivity. Independent recordings can add information even when each individual trajectory remains statistically unresolved.

Sparse identification of nonlinear dynamics represents a vector field through a library and selects active coefficients [@brunton2016]. Integral formulations reduce sensitivity to numerical differentiation [@schaeffer2017], and weak formulations transfer differentiation to test functions [@messenger2021; @reinbold2020]. Their regression matrices still contain noisy states. Numerical resolution, feature bias, residual covariance and model selection interact, so reducing quadrature error does not by itself establish statistical consistency.

The distinction between fixed-domain and increasing-domain observations is classical in Gaussian process inference [@stein1999]. Ying established fixed-domain results for exponential covariance parameters [@ying1991], and Zhang used equivalence of Gaussian measures to explain inconsistent covariance estimation in the Matérn family [@zhang2004]. Those results concern covariance parameters. Here the covariance is common across parameters, and the parameter enters the deterministic mean through an ODE. The same measure-equivalence mechanism supplies an estimator-independent limit, with an explicit norm determined by dynamical sensitivities.

Continuous-time regression with correlated errors has an established optimal-design theory, including optimal signed measures with endpoint masses. Dette, Pepelyshev and Zhigljavsky derive continuous-record best linear unbiased estimators and efficient discrete designs for triangular covariance kernels [@dette2016]. Dynamical sensitivities define a local regression in this framework. Finite limiting precision and endpoint information are therefore precedents for our weak-estimation analysis.

The continuum-data theory of WSINDy studies uncorrelated measurement errors and library-dependent consistency [@messenger2025]. Fixed physical correlation time defines a different experiment: adjacent correlations approach one as the mesh is refined. Correlation induced by filtering independently sampled errors also requires distinguishing a physical filter scale from a fixed number of sample indices. Our analysis concerns a fixed covariance process observed at increasingly dense times.

Our main result identifies the complete boundary-free adjoint influence class and constructs the influence attaining its sharp information reference in covariance-norm closure. Correct temporal covariance and dense test spans then yield a conditional, sequential small-noise efficiency result for WENDy-like estimation. An explicit OU adjoint construction and exact polynomial residual weights make the reference computable. Continuous-regression BLUE, Gaussian projection and Hermite identities are established results [@dette2016; @isserlis1918]; the contribution is their precise connection through the ODE adjoint range, including the boundary restriction and the order of limits. Specifically, we prove (i) the exact closed weak class, (ii) its sharp information reference and an explicit OU optimal test, (iii) conditional dense-span efficiency of correct-covariance GLS with a concrete spline density lemma, and (iv) an estimator-independent fixed-horizon limit. Three dynamical systems and paired likelihood references test the resulting predictions.

## 2. Related estimation methods

Statistical ODE estimation has long used smoothing to reduce repeated forward solves: spline least squares [@varah1982], nonparametric gradient matching [@brunel2008], measurement-error regression [@liang2008] and root-$n$ smooth-and-match estimation [@gugushvili2012]. Direct integral estimators have root-$n$ consistency for suitable linear-in-function parameterizations [@dattner2015]. Orthogonality-condition estimation already supplies weak moments, GMM weighting and boundary-information adaptations [@brunel2014]. Thus neither the qualitative boundary-loss idea nor optimally weighted moments originates here.

Generalized profiling can be asymptotically efficient without explicitly solving the ODE, under its approximation and i.i.d. sampling assumptions [@ramsay2007; @qi2010]. One-step likelihood refinement of an ODE pilot is also established [@dattner2018]; Dattner and Klaassen [@dattner2015] discuss such refinements as a direction beyond their direct-estimator result. We do not propose replacing these efficient full-observation methods. The question is the optimum available to boundary-free weak equations, which also provide tractable moments for library-based structural discovery. The present theorem assumes supplied identifiable structure; the one-step computation below is a likelihood reference on the same correlated-noise records. These small supplied-structure examples remain tractable for likelihood refinement; computational superiority of weak estimation is not inferred from them.

Gaussian regression geometry and small-noise efficiency belong to classical statistical theory [@parzen1961; @ibragimov1981; @kutoyants1994]. Here those principles are connected to the dynamical adjoint range. The distinction from prior qualitative boundary observations is the exact covariance-dual class, its sharp Schur-complement reference and the explicit OU representative.

WENDy treats weak parameter estimation as an errors-in-variables problem and iteratively updates a linearized residual covariance [@bortz2023]. It retains perturbations of both the feature matrix and the response. We use the same coupled residual and compute every Gaussian order for polynomial features under the specified OU model. First-order truncation recovers the structure motivating weak reweighting, with temporal covariance inserted explicitly. Residual covariance is distinct from the covariance of the nonlinear coefficient estimator.

WENDy-MLE optimizes an approximate Gaussian weak-residual likelihood with parameter-dependent covariance and analytic derivatives [@rummel2025]. Its likelihood contains both a log determinant and a covariance-weighted residual term; it is not identical to a frozen-weight GLS update. Tran and Bortz select local test radii by estimating quadrature error [@tran2026], while Fourier Weak SINDy selects informative sinusoidal frequencies spectrally [@chen2026]. These are relevant precedents for weak likelihoods and test design. The complementary question here is the complete adjoint influence class, its sharp initial-state information loss and the explicit OU optimal representative. We do not claim weak GLS or optimized test selection themselves as new methods.

Finite-noise bias and uncertainty are also studied directly for WENDy-IRLS [@chawla2026]. Our bias/dispersion reporting concerns the specified correlated-noise experiment and does not replace these broader coverage studies.

Errors-in-variables effects have a substantial statistical literature [@fuller1987]. Gaussian correction of nonlinear functions of a noisy mean also precedes weak identification [@stefanski1989]. Hermite correction removes pointwise polynomial bias when marginal noise variance is known; it does not remove random design–response coupling. Generalized smoothing and profiling estimate ODE parameters through a different compromise between model fidelity and observation error [@ramsay2007]. The information bound applies to these alternatives when they observe the same path experiment.

WmSINDy combines weak regression with simultaneous dynamics and noise reconstruction [@lopez2026]. We distinguish the large fixed-window regression-component study from a separate paired execution of the complete adaptive author WSINDy, WENDy and WmSINDy pipelines. The latter retains noise estimation or reconstruction and adaptive test-function selection; WENDy supplies structure for parameter estimation. Ensemble-SINDy uses bootstrap aggregation to quantify feature inclusion [@fasel2022]. Resampling one recorded path and acquiring independent physical recordings are different operations: the latter changes the number of independent path experiments in the information calculation.

## 3. Observation model and weak equations

Let $x_\theta:[0,T]\to\mathbb R^d$ be a deterministic ODE trajectory. In the polynomial setting,

$$
x'_i=f_i(x)=\sum_{j=1}^{p}\xi^*_{ji}\vartheta_j(x),\qquad i=1,\ldots,d.
$$

The library contains all true terms. Observations are $Y(t)=x_\theta(t)+\eta_C(t)$, with independent stationary OU coordinates [@uhlenbeck1930] and common-across-parameter covariance

$$
C_j(s,t)=\sigma_j^2\exp(-|s-t|/\tau_j),\qquad \sigma_j>0,\quad\tau_j>0.
$$

The continuous record is a statistical reference experiment; actual data are evaluations $y_n=Y(t_n)$ at mesh $h$. The initial state $u=x(0)$ is known in the baseline acquisition experiment; Sections 5 and 8 compare known and unknown $u$ explicitly. The symbol $Z$ denotes its fundamental sensitivity matrix, $h$ the observation mesh, $\epsilon$ a scalar noise amplitude, and $C$ the continuous covariance operator. Finite weak covariance matrices are denoted by $\Sigma$. In Corollary 1, $C$ is the unit-amplitude covariance shape and the physical covariance is $\epsilon^2 C$.

Fix $K$ test functions $\phi_k\in C_c^2(0,T)$. Their positions and full support widths $L$ remain fixed as $h$ decreases. Let $\mathcal Q_h$ be composite trapezoidal quadrature with weights $w_n$. Clean weak quantities are

$$
G^0_{kj}=\int_0^T\phi_k\vartheta_j(x)\,dt,\qquad
b^0_{ki}=-\int_0^T\phi'_kx_i\,dt,\qquad b^0=G^0\xi^*.
$$

Noisy matrices use corrected polynomial features $\widetilde\vartheta_j$, defined in Section 6:

$$
(G_h)_{kj}=\mathcal Q_h(\phi_k\widetilde\vartheta_j(y)),\qquad
(b_h)_{ki}=-\mathcal Q_h(\phi'_ky_i).
$$

Physical OU infill sets $\rho(h)=\exp(-h/\tau)$ with fixed $T,\tau$. An iid discrete control has independent error of fixed per-sample variance. Fixed-index AR(1) instead keeps $\rho<1$ constant, making the physical correlation time $-h/\log\rho\to0$. These models have different continuum limits.

## 4. Fixed-horizon information limits

### 4.1. OU translation geometry

The stationary OU process satisfies $d\eta_{C,j}=-\eta_{C,j}dt/\tau_j+\sqrt{2\sigma_j^2/\tau_j}\,dW_j$, with an independent initial value $\eta_{C,j}(0)\sim N(0,\sigma_j^2)$. Its Cameron–Martin space is $H^1([0,T];\mathbb R^d)$ as a set, with norm

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

### 4.2. Continuous-record information and scalar weak ratios

**Proposition 1 (continuous-record information).** If $\theta\mapsto x_\theta$ is differentiable in $\mathcal H$ and $s_r=\partial_{\theta_r}x_\theta$, continuous-record Fisher information is

$$
\mathcal I_{rs}(\theta)=\langle s_r,s_s\rangle_{\mathcal H}.
$$

The Gaussian mean-score identity proves the assertion: the directional score has variance equal to the squared Cameron-Martin norm of its sensitivity. Its inverse is the Cramér–Rao covariance bound for regular unbiased estimators, when nonsingular. Biased estimators require the corresponding bias derivative; finite-horizon maximum likelihood need not attain the unbiased bound.

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
A=\int(-\phi'-a\phi)\eta_C,\quad B=\int\phi\eta_C.
$$

OU path continuity gives almost sure convergence of the quadratures. The weights $\phi$ and $-\phi'-a\phi$ are independent: dependence would make a nonzero compactly supported function solve a homogeneous first-order equation. The OU kernel is strictly positive on nonzero weights, so $(A,B)$ is nondegenerate Gaussian. Positive conditional variance of $A$ given $B$ and positive denominator density at zero imply infinite absolute and squared error moments of the unguarded ratio. Its empirical summaries should use quantiles or explicit guards.

The influence $\ell=(-\phi'-a\phi)/D_\phi$ has Gaussian variance $V_\phi=\iint\ell(u)C(u,v)\ell(v)dudv$. Differentiating the weak identity gives $\int\ell s=1$. Cauchy–Schwarz in the covariance RKHS gives

$$
1\leq V_\phi\mathcal I_T(a),\qquad
\operatorname{Eff}_\phi=[\mathcal I_T(a)V_\phi]^{-1}\leq1.
$$

This efficiency concerns the Gaussian linearization; it is not the unguarded ratio's variance.

## 5. Boundary-free weak class and efficiency

### 5.1. Exact boundary-free influence class

Let $J\in C^1([0,T];\mathbb R^{d\times d})$, and let the invertible fundamental matrix satisfy $Z'=JZ$, $Z(0)=I_d$. Put $\mathcal D v=v'-Jv$ and $\mathcal D^*\phi=-\phi'-J^T\phi$. Sensor coordinates have independent stationary OU covariance with fixed positive scales. Its Cameron–Martin space $\mathcal H$ is $H^1$ with the OU inner product. Write $\mathcal E=\mathcal H'$ for its continuous dual, equipped with the covariance norm

$$
\|\lambda\|_C^2=\langle\lambda,C\lambda\rangle,
\qquad
\langle C^{-1}g,h\rangle=\langle g,h\rangle_{\mathcal H}.
$$

This definition includes ordinary integrated influences and endpoint evaluation measures. General dual influences act through the $L^2(P)$ closure of integrated Gaussian observations, rather than a pathwise pairing for every distribution. Their variance norm is $\|\lambda\|_C$; the noise path itself need not belong to $\mathcal H$.

**Proposition 2 (exact range and closure).** In $L^2$,

$$
\mathcal D^*H_0^1
=\left\{\ell\in L^2:\int_0^T Z(t)^T\ell(t)dt=0\right\}.
$$

The range of compactly supported smooth tests is dense in this space in $L^2$. Its closure in $\mathcal E$ is exactly

$$
\mathcal M=\overline{\mathcal D^*C_c^\infty}^{\|\cdot\|_C}
=\{\lambda\in\mathcal E:\langle\lambda,Z_j\rangle=0, j=1,\ldots,d\}.
$$

Consequently $\mathcal M$ has codimension $d$.

**Proof.** Since $Z^T\mathcal D^*\phi=-(Z^T\phi)'$, zero endpoint values imply the moment constraint. Conversely, for an $L^2$ influence with zero moment,

$$
\phi(t)=-Z(t)^{-T}\int_0^t Z(s)^T\ell(s)ds
$$

belongs to $H_0^1$ and satisfies $\mathcal D^*\phi=\ell$. Approximation in $H_0^1$ by $C_c^\infty$ and continuity of $\mathcal D^*:H_0^1\to L^2$ give $L^2$ density. For the stronger closure statement, a function $v\in\mathcal H$ annihilates every $\mathcal D^*\phi$ if and only if $\mathcal D v=0$ distributionally, hence $v=Zc$. The Hahn–Banach annihilator characterization of a closed subspace then gives the stated closure in $\mathcal H'$. The columns of $Z$ are linearly independent elements of $\mathcal H$, so their dual constraints have rank $d$. $\square$

The exact range uses zero-endpoint Sobolev tests; literal compact tests require closure. The covariance norm is equivalent to the dual Sobolev norm for this positive OU model.

**Lemma 1 (a concrete dense test family).** Let $V_m$ be vector-valued continuous piecewise-linear splines on the nested uniform dyadic partitions of $[0,T]$, with zero values at both endpoints. Then

$$
\overline{\bigcup_m\mathcal D^*V_m}^{\|\cdot\|_C}=\mathcal M.
$$

**Proof.** Given $v\in H_0^1$, choose $\psi\in C_c^\infty$ close to $v$ in $H^1$. Its nodal interpolant $I_m\psi\in V_m$ satisfies $\|I_m\psi-\psi\|_{H^1}\leq c\Delta_m\|\psi''\|_{L^2}$, where $\Delta_m$ is the mesh width. The union of $V_m$ is therefore dense in $H_0^1$. Continuity of $\mathcal D^*:H_0^1\to L^2$ and of the embedding $L^2\hookrightarrow\mathcal H'$ transfers this approximation to the covariance norm. Proposition 2 identifies the closure of $\mathcal D^*H_0^1$ with $\mathcal M$. $\square$

These splines are admissible zero-endpoint Sobolev tests; they need not vanish in an endpoint neighborhood. If literal smooth compact tests are required, approximate each spline in $H_0^1$ by such tests. Thus Lemma 1 provides an explicit family for the density condition in Corollary 1. It does not assert density of the finite multiscale bump families used in the experiments.

### 5.2. Sharp information and optimal influence

Let $S=[S_1,\ldots,S_p]\in\mathcal H^p$ be coefficient sensitivities, and form the joint information blocks for $S$ and $Z$. Assume $\mathcal I_{uu}$ and the Schur complement are positive definite. Set

$$
R=S-Z\mathcal I_{uu}^{-1}\mathcal I_{u\theta},
\qquad
\mathcal I_e=\mathcal I_{\theta\theta}
-\mathcal I_{\theta u}\mathcal I_{uu}^{-1}\mathcal I_{u\theta}.
$$

Projection and best linear unbiased estimation in correlated continuous regression are classical [@dette2016]. Here the adjoint-range characterization identifies exactly which dynamical weak influences obey that projection.

**Theorem 2 (sharp boundary-free influence covariance).** Let $\Lambda=[\lambda_1,\ldots,\lambda_p]\in\mathcal M^p$ satisfy $\langle\lambda_r,S_s\rangle=\delta_{rs}$. Its Gaussian influence covariance obeys

$$
V(\Lambda)\succeq\mathcal I_e^{-1}.
$$

The bound is attained in $\mathcal M$ by the matrix-valued dual influence

$$
\Lambda_*=C^{-1}R\,\mathcal I_e^{-1}.
$$

It is the infimum for normalized compact-test influences. It is a minimum among literal compact tests only when $\Lambda_*$ belongs to their adjoint range.

**Proof.** The $\mathcal H$ projection makes $\langle C^{-1}R,Z\rangle=0$ and $\langle C^{-1}R,S\rangle=\mathcal I_e$. Thus $\Lambda_*$ meets both constraints and has covariance $\mathcal I_e^{-1}$. The difference $\Delta=\Lambda-\Lambda_*$ pairs to zero with both $S$ and $Z$, hence with $R$. Covariance cross terms vanish and

$$
V(\Lambda)=\mathcal I_e^{-1}+V(\Delta).
$$

For attainability as an infimum, approximate each column of $\Lambda_*$ by compact-test influences in covariance norm. Their normalization matrices converge to $I_p$ by continuous pairing with $S$; right-multiplication by their inverse transpose enforces exact normalization while preserving the range and covariance convergence. $\square$

For an unknown initial state, $\mathcal I_e$ is the full likelihood nuisance-adjusted information, so this weak restriction has no further first-order cost. When the initial state is known, the information loss is exactly $\mathcal I_{\theta u}\mathcal I_{uu}^{-1}\mathcal I_{u\theta}$; covariance inflation is $\mathcal I_e^{-1}-\mathcal I_{\theta\theta}^{-1}$. There is no direction-independent scalar loss percentage in general. When $u$ is known, $Y(0)-u$ also reveals the initial stationary measurement error, which remains correlated with later observations. The unrestricted score can exploit that observation while the weak influence must retain its initial-sensitivity annihilation constraint.

### 5.3. Explicit optimal OU test

For a scalar coordinate with $g\in H^2$, integration by parts in the stationary OU norm gives the dual measure

$$
C^{-1}g=
\frac{g-\tau^2g''}{2\sigma^2\tau}\,dt
+\frac{g(0)-\tau g'(0)}{2\sigma^2}\delta_0
+\frac{g(T)+\tau g'(T)}{2\sigma^2}\delta_T.
$$

For independent coordinates the formula applies coordinatewise, using their respective scales. If only $g\in H^1$, the Riesz representative remains well-defined in $\mathcal H'$ but need not be a density plus finite endpoint measures.

For a smooth optimal representative write $\Lambda_*=m(t)dt+m_0\delta_0+m_T\delta_T$. On the open interval an adjoint test matrix is

$$
\Phi_*(t)=-Z(t)^{-T}
\left[Z(0)^Tm_0+\int_0^t Z(s)^Tm(s)ds\right].
$$

Extend it by zero outside the observation interval. Its endpoint jumps supply $m_0=-\Phi_*(0+)$ and $m_T=\Phi_*(T-)$, and distributionally $\mathcal D^*\Phi_*=\Lambda_*$. These formulas include the initial atom; omitting it in the cumulative integral is incorrect. The annihilation condition supplies the terminal jump. Cutoff and smoothing approximate this generalized test in covariance norm. These measures describe a statistical reference; the implemented estimators keep compact tests and add no explicit boundary observation terms.

For scalar exponential dynamics, $z=e^{at}$, $s=x_0te^{at}$, $k=\mathcal I_{au}/\mathcal I_{uu}$ and $g=e^{at}(x_0t-k)$. Define

$$
A_0=\frac{(1-\tau^2a^2)x_0}{2\sigma^2\tau},\qquad
B_0=\frac{-(1-\tau^2a^2)k-2\tau^2ax_0}{2\sigma^2\tau},
$$

$$
q_0=\frac{-k-\tau(x_0-ak)}{2\sigma^2},\qquad
F_j(t)=\int_0^t s^je^{2as}ds.
$$

Then

$$
\phi_*(t)=-\frac{e^{-at}}{\mathcal I_e}
[q_0+A_0F_1(t)+B_0F_0(t)].
$$

For $a\ne0$, $F_0=\operatorname{expm1}(2at)/(2a)$ and $F_1=[e^{2at}(2at-1)+1]/(4a^2)$; their continuous limits at $a=0$ are $t$ and $t^2/2$. A pilot coefficient substitutes into this expression before a compact cutoff is applied. The explicit multivariate fundamental-matrix expression is an analytical identity, not a guarantee of stable numerical inversion of $Z(t)$; strongly contracting dynamics can make that inversion ill-conditioned.

### 5.4. Correct-covariance GLS and dense spans

For tests $\phi_k$, put $\lambda_k=\mathcal D^*\phi_k$, $D_{kr}=\langle\lambda_k,S_r\rangle$ and $\Sigma_{k\ell}=\langle\lambda_k,C\lambda_\ell\rangle$. Remove exact linear dependencies. If $D$ has full column rank, covariance-weighted weak regression has local influence and covariance

$$
P=(D^T\Sigma^{-1}D)^{-1}D^T\Sigma^{-1},\qquad
V_K=(D^T\Sigma^{-1}D)^{-1}.
$$

WENDy's covariance iteration motivates this construction [@bortz2023], but the released iid covariance is not the correct physical-OU covariance. Under polynomial features Proposition 3 supplies the complete coupled residual weights. Its order-one term determines the small-noise efficient influence; higher orders correct residual uncertainty at finite amplitude without by themselves removing coefficient-estimator bias.

The same leading map arises from a parameter-dependent Gaussian weak likelihood of the form studied by WENDy-MLE [@rummel2025], provided the covariance shape is correct and the local solution is regular. For fixed dimension write $r(\theta)=\epsilon r_1+o_p(\epsilon)$, $\partial_\theta r=-D+o_p(1)$ and $\Sigma_\epsilon(\theta)=\epsilon^2\Sigma_0(\theta)$, with smooth positive definite $\Sigma_0$. The objective $\log\det\Sigma_\epsilon+r^T\Sigma_\epsilon^{-1}r$ has covariance-derivative and log-determinant score terms of order $O_p(1)$ in an $O_p(\epsilon)$ neighborhood; its mean-residual score is of order $O_p(\epsilon^{-1})$. The local score equation therefore gives

$$
\widehat\theta-\theta=\epsilon
(D^T\Sigma_0^{-1}D)^{-1}D^T\Sigma_0^{-1}r_1+o_p(\epsilon).
$$

Thus the log determinant and covariance derivatives distinguish finite-noise weak likelihood from frozen GLS, while both share the correct-covariance leading influence under these assumptions. This connects Corollary 1 to a regular local weak-likelihood estimator; it does not supply global optimization or covariance-calibration consistency for a released implementation.

**Corollary 1 (small-noise, dense-span efficiency).** Let the noise be $\epsilon\eta_C$ with known covariance shape $C$. Assume a smooth identifiable parameterization, a generating parameter interior to a compact parameter set, consistent pilots, consistent covariance-shape weights and quadrature error $o(\epsilon)$. For each fixed finite test family assume a regular local weak solution with the GLS linearization. For the covariance conclusion additionally assume uniform integrability of $\|\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)\|^2$. A sufficient implementation uses a uniformly stable smooth local solve, Gaussian moment control on the local Taylor remainder, and bounded fallback outside the regular region on a bad event whose probability is $o(\epsilon^2)$. Then

$$
\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)
\Rightarrow N(0,V_K),\qquad
\epsilon^{-2}\operatorname{Cov}(\widehat\theta_{K,\epsilon})\to V_K.
$$

If nested test-influence spans are dense in $\mathcal M$, then

$$
V_K\longrightarrow\mathcal I_{e,C}^{-1}.
$$

**Proof.** For fixed $K$, a first-order expansion of the weak equations gives $D(\widehat\theta-\theta)=\epsilon r_1+o_p(\epsilon)$ with the GLS map above. Consistent pilots and weights change the leading map by $o_p(1)$; derivatives of pilot-dependent tests multiply the zero clean residual. Gaussianity gives the distributional assertion; the stated risk condition justifies covariance convergence. In the Gaussian dual space, finite-family information is the Gram matrix of the orthogonal projection of $C^{-1}R$ onto its test span. Density gives projection convergence, and nested spans give increasing information; continuity of inversion at positive definite $\mathcal I_e$ gives covariance convergence. $\square$

The conclusion is sequential: first small noise and controlled quadrature for fixed $K$, then dense spans. A diagonal joint limit requires uniform remainder, conditioning and covariance-weight control that are not supplied here. This theorem concerns regular WENDy-like weak estimation with correct covariance and supplied identifiable structure. The conditional result covers the leading GLS map, rather than an algorithm name or finite-noise support-selection procedure. A fixed basis cutoff must be compatible with the density assumption; discarding an informative direction permanently can invalidate it.


## 6. Exact polynomial weak-residual covariance

### 6.1. Finite Hermite expansion

Assume a common $\tau$ for independent coordinates, and set $v_j=\sigma_j^2$. Define $H_0=1$, $H_1=u$, $H_{q+1}=uH_q-qvH_{q-1}$ and corrected monomials $\widetilde\vartheta_\beta(y)=\prod_jH_{\beta_j}(y_j;v_j)$. Their expectation given $x$ is $x^\beta$. This Gaussian debiasing construction [@stefanski1989] centers the polynomial residual at the true coefficients.

For continuum residual $r_{ki}=b_{ki}-(G\xi^*)_{ki}$, polynomial expansion gives

$$
r_{ki}=\sum_{1\leq|\alpha|\leq p_*}\int_0^Tc_{ki,\alpha}(t)H_\alpha(\eta_C(t);v)dt,
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

The order-one covariance supplies the GLS weights required by Corollary 1. Higher orders refine residual weights at finite amplitude; they do not alone remove coefficient bias. Order $q$ has correlation time $\tau/q$. Higher orders need not be negligible at moderate amplitude. Different windows and equations remain correlated despite independent sensor coordinates. The covariance is exact for a supplied clean trajectory and true polynomial; plugging in estimates introduces another approximation.

### 6.2. Raw features and efficient evaluation

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

## 7. Sampling and the scope of efficiency

Theorem 1 applies at fixed physical noise covariance and fixed horizon. Corollary 1 concerns a different limit: decreasing noise for each fixed test family, followed by dense test spans. Increasing the grid density improves the quadrature approximation; it does not replace either limit. Increasing the number of tests improves approximation of an admissible influence but cannot recover the initial-state information excluded by Proposition 2.

Independent records add their sensitivity information. Longer records add useful information only while parameter sensitivities remain appreciable. For scalar decay, information approaches a finite limit as duration increases. These distinctions and numerical acquisition and aggregation checks are detailed in S9. Averaging corrected weak matrices and stacking noisy regressions are different estimating equations; the latter can retain errors-in-variables bias [@fuller1987].

## 8. Tests of the efficiency prediction

The principal comparisons use supplied identifiable structure and known covariance shape; an additional sensitivity check fits the noise shape from data. They test feasible weak-estimation components, separately from the complete author pipelines archived in S6 and S9.6. No true parameter is supplied to a feasible pilot. Deterministic geometry calculations use true parameters only as local references. Protocols, seeds, source fingerprints and output archives are preserved in the public reproducibility archive.

@@NUMERICAL@@

## 9. Discussion and limitations

The adjoint-range characterization turns the initial-sensitivity constraint into an exact description of the closed influence class. Its sharp optimum follows from classical Gaussian projection, and the OU inverse covariance makes an optimal generalized test explicit. Endpoint atoms explain why compact tests generally approach the optimum without attaining it. The constructed scalar cutoff uses compact tests throughout; this article does not test an estimator with additional boundary terms.

Correct covariance weighting is optimal within a fixed span. Density in the covariance norm is required to approach the class bound. A well-conditioned representation changes computation, not the underlying span. An explicit fundamental-matrix inverse may itself be ill-conditioned for a contracting dynamical direction. The multiscale calculation therefore uses a fixed, data-independent basis transformation and reports retained dimensions and numerical relaxation.

The efficiency conclusion applies to the local GLS map under the regularity and risk assumptions of Corollary 1. Global covariance-iteration behavior, sparse support selection and simultaneous test-growth/noise-decay limits require additional analysis. Lemma 1 supplies a dense spline family, while the finite bump experiment tests approximation at three declared dimensions. Its iid comparison concerns this specific misspecification.

Exact residual covariance and marginal Hermite correction address distinct parts of the weak errors-in-variables problem. At finite noise, random design and pilot errors can still produce bias and excess dispersion. The numerical comparison therefore reports RMSE, bias, dispersion, stopping diagnostics and uncertainty, rather than equating a median error with a covariance bound. The separate 16-path complete-pipeline benchmark is retained as context; the matched 128-path WENDy comparison and the 2000-path small-noise weighting study address distinct practical questions.

The analysis assumes fully observed deterministic ODEs and additive stationary Gaussian measurement error with known positive covariance common across parameters. Unknown covariance, partial observation and process noise are outside the theorem proved here. The feasible vector implementation uses one fitted-trajectory solve to stabilize covariance weights. It therefore tests a statistical efficiency mechanism and does not establish a forward-solve-free or runtime-optimal algorithm. The result provides a local information reference for boundary-free weak estimation. The matched finite-noise comparison and the full-pipeline benchmark illustrate its practical scope; they do not establish a general algorithm ranking.

## Declaration of generative AI and AI-assisted technologies

AI-assisted tools supported literature searches, mathematical drafting, source development, experiment orchestration, numerical interpretation and manuscript preparation. Numerical results are derived from archived execution outputs, and figures are produced by plotting code. The tool is not an author. The human authors must review and edit the manuscript and accept responsibility for its content before submission; this statement does not assert that their review has already occurred.

## Data availability

Code, synthetic observations, numerical outputs and frozen execution notebooks are available in the public [reproducibility repository, version v1.0.1](https://github.com/malyshevdanil/wsindy-information-limits/tree/v1.0.1). The repository includes protocol configurations, source fingerprints, pinned author-code download specifications, compatibility adapters, postprocessing scripts and SHA256 manifests. Original author sources are downloaded from their pinned repositories rather than redistributed. No external experimental dataset is used. The tagged snapshot identifies the version supporting this article; a DOI has not been assigned.

## References
