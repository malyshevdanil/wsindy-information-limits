---
title: "The information cost of boundary-free weak estimation under correlated noise"
author: "Danil Romanovich Malyshev"
lang: en-US
---

Independent researcher, Almaty, Kazakhstan

Corresponding author: malyshevdanil04@gmail.com

ORCID: https://orcid.org/0009-0007-0947-3866

**MSC 2020:** Primary 62F12; secondary 62M10, 34A55.

## Abstract

Boundary-free weak equations discard local information when initial states are known. For deterministic ODEs observed with stationary Ornstein–Uhlenbeck measurement error and supplied covariance, we characterize their closed adjoint influence class and its sharp nuisance-adjusted information bound. The restriction raises the local RMSE reference by 29% for Van der Pol and almost doubles it for Lorenz. Allowing known-initial traces recovers the full local influence class. This boundary adaptation is due to Brunel and coauthors; universal time-only windows recover about 97% of the empirical gap in Van der Pol and 95% in Lorenz at small noise. An ablation attributes the recovery to initial traces. For a prespecified cubic candidate, independent small-noise Wald validation passes the frozen size tolerance; power agreement depends on perturbation magnitude. MSTLS inclusion concerns a different structure-selection task. At moderate noise, endpoint recovery is incomplete and nominal Wald errors require estimator-level calibration. The contribution is an exact weak-class characterization and quantitative correlated-noise benchmarks, rather than a new boundary remedy or general selection procedure.

**Keywords:** weak-form identification; correlated measurement error; fixed-domain asymptotics; Fisher information; polynomial covariance; errors in variables.

## 1. Introduction

Boundary-vanishing weak equations can lose information even when an initial state is known, as recognized by orthogonality-condition estimation [@brunel2014] and recent weak identifiability work [@heitzman2026]. We quantify this effect: the closed class of boundary-free weak influences has exactly the local information available when the initial state is treated as a nuisance parameter in the full path likelihood. For an unknown initial state, this restriction causes no additional first-order information loss; for a known initial state, the lost information is the initial-sensitivity projection, and the retained information is its Schur complement. The equivalence concerns local information.

Known initial states also permit a remedy. Initial-trace conditions are established in orthogonality estimation [@brunel2014, Section 2.3]. We characterize their complete local influence class, quantify finite-family recovery under physical temporal correlation, and distinguish their contribution from an initial-noise control. The control reveals a correlated error whose covariance with interior moments can reduce variance. The ablation determines its practical contribution rather than presuming it drives recovery.

Sparse identification of nonlinear dynamics represents a vector field through a library and selects active coefficients [@brunton2016]. Integral formulations reduce sensitivity to numerical differentiation [@schaeffer2017], and weak formulations transfer differentiation to test functions [@messenger2021; @reinbold2020]. Their regression matrices still contain noisy states. Statistical consistency depends jointly on feature correction, residual covariance and the acquisition limit.

The distinction between fixed-domain and increasing-domain observations is classical in Gaussian process inference [@stein1999]. Ying established fixed-domain results for exponential covariance parameters [@ying1991], and Zhang used equivalence of Gaussian measures to explain inconsistent covariance estimation in the Matérn family [@zhang2004]. Those results concern covariance parameters. Here the covariance is common across parameters, and the parameter enters the deterministic mean through an ODE. The same measure-equivalence mechanism supplies an estimator-independent limit, with an explicit norm determined by dynamical sensitivities.

Optimal estimating equations and their sandwich information are classical [@godambe1960]. Efficient scores remove nuisance-score directions by orthogonal projection [@bickel1993; @vandervaart1998, Chapter 25]; the initial-state Schur complement is the finite-dimensional instance here. Continuous-time BLUE and endpoint-weighted estimation likewise precede weak ODE analysis [@dette2016; @dette2017]. Our contribution identifies the adjoint range of the restricted weak class and quantifies its cost and recovery in this dynamical experiment.

The continuum-data theory of WSINDy studies uncorrelated measurement errors and library-dependent consistency [@messenger2025]. Fixed physical correlation time defines a different experiment: adjacent correlations approach one as the mesh is refined. Correlation induced by filtering independently sampled errors also requires distinguishing a physical filter scale from a fixed number of sample indices. Our analysis concerns a fixed covariance process observed at increasingly dense times.

The paper has three contributions: exact boundary-free and known-initial adjoint-range characterizations; a constructive efficient full-library estimator under a checkable jet condition, verified exactly for Van der Pol (Theorem 3); and quantitative correlated-noise benchmarks separating information cost, established boundary recovery and prespecified testing. Polynomial Hermite covariance supplies the feasible weights. Independent-record GMM, fission and cross-fitting are preserved in a separate research archive.

## 2. Related estimation methods

ODE estimation has long used trajectory smoothing and integral matching: spline least squares [@varah1982], nonparametric gradient matching [@brunel2008], measurement-error regression [@liang2008], smooth-and-match estimation [@gugushvili2012], and direct integral estimation [@dattner2015]. Generalized profiling can be efficient under its approximation and independent-sampling assumptions [@ramsay2007; @qi2010]; one-step likelihood refinement of an ODE pilot is also established [@dattner2018]. These provide full-observation references. Our question is the information attainable by restricted weak equations under physical temporal correlation.

WENDy treats weak parameter estimation as an errors-in-variables problem, retaining coupled feature and response perturbations and iteratively updating their linearized covariance [@bortz2023]. WENDy-MLE uses an approximate Gaussian weak-residual likelihood with parameter-dependent covariance and analytic derivatives [@rummel2025]. Its log determinant and covariance-derivative terms distinguish its finite-noise objective from frozen-weight GLS. Polynomial Hermite correction and exact residual covariance here extend that coupled calculation to all Gaussian orders for supplied OU noise. They do not alone correct the nonlinear coefficient estimator. Gaussian correction of noisy polynomial features has classical precedents [@fuller1987; @stefanski1989], and WENDy finite-noise bias and uncertainty remain active topics [@chawla2026].

Test design also has substantial precedents. Tran and Bortz choose local radii through quadrature-error estimates [@tran2026], while Fourier Weak SINDy selects sinusoidal frequencies spectrally [@chen2026]. WmSINDy combines weak regression with dynamics and noise reconstruction [@lopez2026]. Ensemble-SINDy quantifies feature inclusion by bootstrap aggregation [@fasel2022]. Our author-method checks preserve their distinct structural inputs and are documented in the research archive. The new class characterization and constructive full-library result complement these methods; the supplied-candidate Wald comparison establishes no general superiority for structure selection.

### 2.1. Complementarity with orthogonality-condition estimation

Brunel et al. [@brunel2014] combine a nonparametric proxy, orthogonality conditions and GMM in an independent-observation, increasing-point experiment, including delay equations. They recognize boundary-information loss and adapt moments to known boundary values. Our deterministic polynomial setting is narrower. It instead treats fixed-horizon physical OU error, identifies the complete adjoint influence classes, gives their sharp initial-state Schur-complement relation and an explicit OU optimal test, and constructs an attainable full-library estimator under a checkable jet condition. The small-noise and sequential dense-span limits differ from increasing-point asymptotics. The results complement each other.

## 3. Observation model and weak equations

Let $x_\theta:[0,T]\to\mathbb R^d$ be a deterministic ODE trajectory. In the polynomial setting,

$$
x'_i=f_i(x)=\sum_{j=1}^{p}\xi^*_{ji}\vartheta_j(x),\qquad i=1,\ldots,d.
$$

The library contains all true terms. Observations are $Y(t)=x_\theta(t)+\eta_C(t)$, with independent stationary OU coordinates [@uhlenbeck1930] and common-across-parameter covariance

$$
C_j(s,t)=\sigma_j^2\exp(-|s-t|/\tau_j),\qquad \sigma_j>0,\quad\tau_j>0.
$$

The continuous record is a statistical reference experiment; actual data are evaluations $y_n=Y(t_n)$ at mesh $h$. The initial state $u=x(0)$ is known in the baseline acquisition experiment; Sections 5 and 7 compare known and unknown $u$ explicitly. The symbol $Z$ denotes its fundamental sensitivity matrix, $h$ the observation mesh, $\epsilon$ a scalar noise amplitude, and $C$ the continuous covariance operator. Finite weak covariance matrices are denoted by $\Sigma$. In Proposition 3, $C$ is the unit-amplitude covariance shape and the physical covariance is $\epsilon^2 C$.

Fix $K$ test functions $\phi_k\in C_c^2(0,T)$. Their positions and full support widths $L$ remain fixed as $h$ decreases. Let $\mathcal Q_h$ denote composite trapezoidal quadrature with weights $w_n$ in the formulas; matched endpoint confirmations use composite Simpson quadrature as specified in their frozen protocols. Clean weak quantities are

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

**Lemma 1 (fixed-horizon two-point bound; proof in S1.1).** Suppose the mean trajectories belong to $\mathcal H$ and share the stated positive OU covariance. For $\delta=\|\theta-\eta\|>0$ and $d=\|x_\theta-x_\eta\|_{\mathcal H}$, every finite-grid estimator satisfies

$$
\max_{v\in\{\theta,\eta\}}P_v\{\|\widehat\theta_h-v\|\geq\delta/2\}
\geq\Phi(-d/2)>0.
$$

Its maximal two-point squared risk is at least $(\delta^2/4)\Phi(-d/2)$. Thus dense-grid consistency throughout a nontrivial fixed-horizon family is excluded by Gaussian path equivalence [@feldman1958; @zhang2004].

### 4.2. Continuous-record information

**Proposition 1 (continuous-record information).** If $\theta\mapsto x_\theta$ is differentiable in $\mathcal H$ and $s_r=\partial_{\theta_r}x_\theta$, continuous-record Fisher information is

$$
\mathcal I_{rs}(\theta)=\langle s_r,s_s\rangle_{\mathcal H}.
$$

The Gaussian mean-score identity proves the assertion: the directional score has variance equal to the squared Cameron-Martin norm of its sensitivity. Its inverse is the Cramér–Rao covariance bound for regular unbiased estimators, when nonsingular. Biased estimators require the corresponding bias derivative; finite-horizon maximum likelihood need not attain the unbiased bound.

Exact discrete OU innovations converge to this finite path information under grid refinement. S26.1 gives the scalar calculation and shows why unguarded weak ratios can have infinite error moments despite a finite Gaussian influence variance. Information bounds concern local regular influences; finite-noise estimator risk requires separate control.

## 5. Boundary-free weak class and efficiency

### 5.1. Exact boundary-free influence class

Along $x_{\theta,u}$, put $J(t)=D_xf(x_{\theta,u}(t),\theta)$ and $Z(t)=\partial_u x_{\theta,u}(t)$. Assume $J\in C^1([0,T];\mathbb R^{d\times d})$, and let the invertible fundamental matrix satisfy $Z'=JZ$, $Z(0)=I_d$. Put $\mathcal D v=v'-Jv$ and $\mathcal D^*\phi=-\phi'-J^T\phi$. Sensor coordinates have independent stationary OU covariance with fixed positive scales. Its Cameron–Martin space $\mathcal H$ is $H^1$ with the OU inner product. Write $\mathcal E=\mathcal H'$ for its continuous dual, equipped with the covariance norm

$$
\|\lambda\|_C^2=\langle\lambda,C\lambda\rangle,
\qquad
\langle C^{-1}g,v\rangle=\langle g,v\rangle_{\mathcal H}.
$$

This definition includes ordinary integrated influences and endpoint evaluation measures. General dual influences are defined by $L^2(P)$ limits of integrated Gaussian observations. Their variance norm is $\|\lambda\|_C$; the noise path itself need not belong to $\mathcal H$.

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

**Lemma 2 (a concrete dense test family).** Let $V_m$ be vector-valued continuous piecewise-linear splines on the nested uniform dyadic partitions of $[0,T]$, with zero values at both endpoints. Then

$$
\overline{\bigcup_m\mathcal D^*V_m}^{\|\cdot\|_C}=\mathcal M.
$$

**Proof.** Given $v\in H_0^1$, choose $\psi\in C_c^\infty$ close to $v$ in $H^1$. Its nodal interpolant $I_m\psi\in V_m$ satisfies $\|I_m\psi-\psi\|_{H^1}\leq c\Delta_m\|\psi''\|_{L^2}$, where $\Delta_m$ is the mesh width. The union of $V_m$ is therefore dense in $H_0^1$. Continuity of $\mathcal D^*:H_0^1\to L^2$ and of the embedding $L^2\hookrightarrow\mathcal H'$ transfers this approximation to the covariance norm. Proposition 2 identifies the closure of $\mathcal D^*H_0^1$ with $\mathcal M$. $\square$

These splines are admissible zero-endpoint Sobolev tests; they need not vanish in an endpoint neighborhood. If literal smooth compact tests are required, approximate each spline in $H_0^1$ by such tests. Thus Lemma 2 provides an explicit family for the density condition in Proposition 3. The experiments use finite multiscale approximations.

### 5.2. Sharp information and optimal influence

Let $S=[S_1,\ldots,S_p]\in\mathcal H^p$ be coefficient sensitivities, and form the joint information blocks for $S$ and $Z$. Assume $\mathcal I_{uu}$ and the Schur complement are positive definite. Set

$$
R=S-Z\mathcal I_{uu}^{-1}\mathcal I_{u\theta},
\qquad
\mathcal I_e=\mathcal I_{\theta\theta}
-\mathcal I_{\theta u}\mathcal I_{uu}^{-1}\mathcal I_{u\theta}.
$$

The nuisance-score projection and efficient influence are classical [@bickel1993; @vandervaart1998, Chapter 25], as is continuous-record BLUE [@dette2016; @dette2017]. Here the adjoint-range characterization identifies exactly which dynamical weak influences obey that projection.

**Theorem 1 (sharp boundary-free influence covariance).** Let $\Lambda=[\lambda_1,\ldots,\lambda_p]\in\mathcal M^p$ satisfy $\langle\lambda_r,S_s\rangle=\delta_{rs}$. Its Gaussian influence covariance obeys

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

Known-boundary adaptations have a direct precedent: Brunel et al. [@brunel2014, Section 2.3] allow a nonzero initial test trace when the initial state is known and retain a zero terminal trace. Their construction motivates endpoint-enabled moments. The information comparison above concerns the endpoint-vanishing class; restoring endpoint information changes that class. Recent test construction work also treats numerical integration error as a central design criterion [@tran2026]. An implementation of endpoint-enabled tests therefore needs a separate clean quadrature check before its empirical accuracy can be compared with the known-initial likelihood reference.

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

Extend it by zero outside the observation interval. Its endpoint jumps supply $m_0=-\Phi_*(0+)$ and $m_T=\Phi_*(T-)$, and distributionally $\mathcal D^*\Phi_*=\Lambda_*$. These formulas include the initial atom; omitting it in the cumulative integral is incorrect. The annihilation condition supplies the terminal jump. Cutoff and smoothing approximate this generalized test in covariance norm. These measures describe the boundary-free statistical reference. Section 7.2 separately evaluates an endpoint-enabled extension for a known initial state.

The scalar exponential specialization and fitted-pilot formula are in S26.5. Figure 1 shows how the generalized optimum and its compact cutoff differ near the endpoints.

![Figure 1. Scalar optimal generalized test for $a=-1$, $x_0=2$, $T=4$ and $\tau=0.4$, alongside its compact cutoff with width 0.1 and an ordinary width-three bell. All tests use the same clean-moment normalization. Endpoint panels show nonzero one-sided traces and jumps of the zero extension; their adjoint derivatives supply endpoint atoms. The curves describe influences annihilating the initial sensitivity, under the unknown-initial information reference.](figures/publication/optimal_test_function.png)

### 5.4. Correct-covariance GLS and dense spans

For tests $\phi_k$, put $\lambda_k=\mathcal D^*\phi_k$, $D_{kr}=\langle\lambda_k,S_r\rangle$ and $\Sigma_{k\ell}=\langle\lambda_k,C\lambda_\ell\rangle$. Remove exact linear dependencies. If $D$ has full column rank, covariance-weighted weak regression has local influence and covariance

$$
P=(D^T\Sigma^{-1}D)^{-1}D^T\Sigma^{-1},\qquad
V_K=(D^T\Sigma^{-1}D)^{-1}.
$$

WENDy's covariance iteration motivates this construction [@bortz2023], and physical-OU covariance changes its temporal weights. Under polynomial features Proposition 4 supplies the complete coupled residual weights. Its order-one term determines the small-noise efficient influence; higher orders correct residual uncertainty at finite amplitude without by themselves removing coefficient-estimator bias.

This is the Godambe-information construction for optimally weighted estimating equations [@godambe1960]: the local information is $D^\top\Sigma^{-1}D$. WENDy-MLE has the same leading small-noise map under regularity and correctly specified covariance; finite-noise covariance-derivative terms distinguish its objective (S26.3).

**Proposition 3 (conditional small-noise, dense-span efficiency).** For known unit-amplitude covariance shape $C$, assume the fixed-family GLS expansion, consistent feasible weights, quadrature error $o(\epsilon)$, and uniform integrability of squared scaled errors stated as (P1)–(P3) in S26.2. Then

$$
\epsilon^{-1}(\widehat\theta_{K,\epsilon}-\theta)\Rightarrow N(0,V_K),\qquad
\epsilon^{-2}\operatorname{Cov}(\widehat\theta_{K,\epsilon})\to V_K.
$$

If the nested influence spans are dense in $\mathcal M$ and $\mathcal I_{e,C}$ is positive definite, then $V_K\to\mathcal I_{e,C}^{-1}$. The limits are sequential: small noise for fixed $K$, then dense spans. S26.2 proves both assertions under the explicit assumptions. A joint limit or risk claim for unguarded ratios requires additional analysis.

The archived implementation has a fixed relative covariance ridge that survives normalization by $\epsilon^2$ and an unbounded final fit; S27 derives its regularized limiting covariance. Theorem 3 gives a separate constructive realization whose assumptions are verified from a jet condition. S29 independently validates the supplied-candidate Wald procedure within frozen tolerance bands. These three procedures have distinct theoretical and empirical scopes (Table S27).

### 5.5. Extended class and recovery with a known initial state

Throughout this subsection, $u$ is known exactly, covariance is supplied and independent of $\theta$, $J$ is continuous, and $S$ and $Z$ belong to the OU Cameron–Martin space $\mathcal H$. Assume $\mathcal I_{\theta\theta}$ is positive definite. The coordinate OU scales and correlation times are positive. Consequently $\mathcal E$ is the covariance-dual space of vector-valued $H^1([0,T])$, with continuous endpoint evaluation. The initial sensitivity is normalized by $Z(0)=I_d$ and the parameter sensitivity satisfies $S(0)=0$.

For $\phi\in H^1([0,T];\mathbb R^d)$, use the known-initial residual

$$
m_\phi(\theta)=-\int_0^T[\phi'^TY+\phi^Tf_\theta^H(Y)]dt
+\phi(T)^TY(T)-\phi(0)^Tu.
$$

Its first-order influence is $L\phi=(-\phi'-J^T\phi)dt+\phi(T)\delta_T$. The initial boundary term contains the exact state rather than its noisy observation. This adaptation of the weak identity is established [@brunel2014, Section 2.3]. The result below specifies its covariance-dual geometry and its relation to the correlated-noise control.

**Theorem 2 (known-initial extended influence class).** Under these assumptions,

$$
\overline{L(H^1)}^{\mathcal E}
=\overline{L\{\phi\in H^1:\phi(T)=0\}}^{\mathcal E}
=\mathcal E.
$$

In contrast, allowing only a terminal trace preserves the original class:

$$
\overline{L\{\phi\in H^1:\phi(0)=0\}}^{\mathcal E}=\mathcal M.
$$

Alternatively, adding the initial noise controls to the closed boundary-free class gives the exact equality

$$
\mathcal M+\{a\delta_0:a\in\mathbb R^d\}=\mathcal E.
$$

The sharp normalized influence covariance in either full extended closure is $\mathcal I_{\theta\theta}^{-1}$, with representative $\Lambda_{\mathrm{full}}=C^{-1}S\mathcal I_{\theta\theta}^{-1}$. For ordinary tests it is an infimum; an exact smooth-test minimum additionally requires the representative to belong to the specified range. Explicit initial observation controls can realize the initial atom.

*Proof.* For $g\in\mathcal H$, the continuous dual pairing obeys

$$
\langle g,L\phi\rangle
=\int_0^T\phi^T(g'-Jg)dt+\phi(0)^Tg(0).
$$

If this vanishes for all permitted tests, compact interior tests first imply $g'=Jg$ in the weak sense. Such solutions are $g=Zc$. For all $H^1$ tests, or for tests with $\phi(T)=0$, the initial trace is arbitrary, so $g(0)=c=0$. Their annihilator is trivial and their closure is $\mathcal E$. For tests with $\phi(0)=0$, every $Zc$ annihilates the range and there are no further annihilators. Its closure is therefore $\mathcal M$, as characterized in Theorem 1.

For any $\lambda\in\mathcal E$, define $a_j=\langle Z_{\cdot j},\lambda\rangle$. Continuity of this pairing follows from $Z\in\mathcal H$. Since $Z(0)=I_d$, $\lambda-a\delta_0$ annihilates every initial-sensitivity column and belongs to $\mathcal M$. This proves the exact sum identity. For covariance optimality, $\langle\Lambda_{\mathrm{full}},S\rangle=I_p$ and $\operatorname{Cov}(\Lambda_{\mathrm{full}})=\mathcal I_{\theta\theta}^{-1}$. Any other normalized influence has a covariance-orthogonal remainder, giving a positive-semidefinite covariance increment. Density, continuous pairing with $S$ and normalization by the convergent sensitivity matrix give attainability as an infimum, by the same approximation argument as Theorem 1. $\square$

Correct-covariance GLS on dense universal extended spans approaches this bound under the sequential small-noise conditions of Proposition 3, with the full rather than nuisance-adjusted information. The theorem concerns the closed local influence class; performance of finite families and finite-noise covariance pilots requires numerical assessment. In particular, finite-span gains of different endpoint additions can overlap.

For compact tests, let $\Omega$ be the exact Hermite residual covariance, $\Sigma_0=\operatorname{diag}(\sigma_j^2)$ and $B=\operatorname{Cov}(m_\phi,Y(0)-u)$. The optimal linear control replaces the residual covariance by

$$
\Omega-B\Sigma_0^{-1}B^T.
$$

Only the first Hermite order contributes to $B$. In coordinate $j$, its cross-covariance is the integral of the linear residual weight against $\sigma_j^2e^{-t/\tau_j}$. In a sampled iid model the compact-test weights and derivatives vanish at the initial node; hence this cross-covariance is exactly zero. This distinction explains the correlated-noise control mechanism. It does not establish an additional finite-family gain when initial-trace tests already approximate the full efficient directions.


### 5.6. A constructive efficient full-library estimator

**Theorem 3 (constructive attainability under a jet condition).** Consider an analytic polynomial ODE trajectory $x_0$ on $[0,T]$, $T>1$, with vector field linear in a polynomial library $g_1,\ldots,g_m$ in each of $d$ equations. Let the $dm$ true coefficients lie strictly inside a supplied compact box. Observe $Y=x_0+\epsilon\eta$ with known positive independent-coordinate stationary OU covariance. Assume that, at some $t_*$ in the observation interval,

$$
\det[\partial_t^r g_j(x_0(t_*))]_{r=0,\ldots,m-1;\ j=1,\ldots,m}\ne0.
$$

Use the full-library smoothed-data pilot, Hermite-corrected moments, cubic-spline tests, vanishing relative ridge and bounded final estimate specified in S28.0, with $h^2=o(\epsilon)$. For each fixed spline level $k$, the pilot is root-$\epsilon$ consistent, the normalized covariance weights are consistent, and

$$
\epsilon^{-1}(\widehat\theta_{k,\epsilon}-\theta_0)\Rightarrow N(0,V_k),\qquad
\epsilon^{-2}E[(\widehat\theta_{k,\epsilon}-\theta_0)(\widehat\theta_{k,\epsilon}-\theta_0)^\top]\to V_k.
$$

The squared scaled errors are uniformly integrable, and $V_k\to\mathcal I_{\theta\mid u}^{-1}$ as $k\to\infty$ after the fixed-level small-noise limit. The coefficient pilot uses all library terms; no generating support is supplied. The formulas and proof are in exact arithmetic. S28.0 proves the general result, including rank, weight consistency and moment control; S28.1 verifies it for Van der Pol with the full 18-coefficient cubic library and $u=(2,0)$, where the jet determinant is $-65265167726360395776$. This theorem concerns the specified guarded estimator and sequential limits, not finite-noise Wald calibration or the archived fixed-ridge fits.

## 6. Exact polynomial weak-residual covariance

### 6.1. Finite Hermite expansion

Assume a common $\tau$ for independent coordinates, and set $v_j=\sigma_j^2$. Define $H_0(z;v)=1$, $H_1(z;v)=z$, $H_{q+1}(z;v)=zH_q(z;v)-qvH_{q-1}(z;v)$ and corrected monomials $\widetilde\vartheta_\beta(y)=\prod_jH_{\beta_j}(y_j;v_j)$. Their expectation given $x$ is $x^\beta$. This Gaussian debiasing construction [@stefanski1989] centers the polynomial residual at the true coefficients.

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

**Proposition 4 (exact covariance).** Under these assumptions,

$$
\operatorname{Cov}(r_{ki},r_{\ell m})=
\sum_{1\leq|\alpha|\leq p_*}\alpha!\prod_jv_j^{\alpha_j}
\int_0^T\int_0^Tc_{ki,\alpha}(s)c_{\ell m,\alpha}(t)
e^{-|\alpha||s-t|/\tau}\,ds\,dt.
$$

**Proof.** The generating function $\exp(az-va^2/2)$ yields $E[H_q(U;v)H_r(V;v)]=\mathbf1_{q=r}q!v^q\rho^q$ for jointly Gaussian variables of variance $v$ and correlation $\rho$. This is the Hermite form of classical Gaussian moment identities [@isserlis1918]. Coordinate independence gives orthogonality between different multi-indices and factor $\alpha!\prod_jv_j^{\alpha_j}$. Substitute $\rho=e^{-|s-t|/\tau}$ and integrate the finite expansion. $\square$

The order-one covariance supplies the GLS weights required by Proposition 3. Higher orders refine residual weights at finite amplitude; coefficient bias also depends on random-design terms. Order $q$ has correlation time $\tau/q$. Higher orders need not be negligible at moderate amplitude. Different windows and equations remain correlated despite independent sensor coordinates. The covariance is exact for a supplied clean trajectory and true polynomial; plugging in estimates introduces another approximation.

### 6.2. Implementation and covariance misspecification

The finite Hermite sum can be evaluated on a grid without storing a dense temporal matrix. S26.4 gives the recursion and the extension to incorrectly corrected raw features. Under temporal weights $W$, leading covariance is $P_W\Sigma_{\mathrm{true},1}P_W^\top$; changing correction variances can also change the moment mean. The present comparisons supply the sensor covariance. Sensitivity to supplied and fitted correlation times is reported in S25.4.

## 7. Numerical benchmarks: cost, recovery and targeted selection

The benchmarks use fully observed polynomial systems, supplied structure and sensor covariance, and archived paired observations. Feasible covariance pilots receive no generating coefficients. S22–S25 preserve protocols, seed ranges, fingerprints and additional checks. Local references and finite-noise performance are reported separately.

### 7.1. Information cost and finite-span efficiency

For Van der Pol, we fit four independent coefficients in $x'=av$, $v'=bx+cv+dx^2v$, with truth $(1,-1,1,-1)$, initial state $(2,0)$, $T=20$, $h=0.005$ and OU range 0.3. Relative errors are divided by the true coefficient norm, two. The known-initial local RMSE reference at 5% noise is 0.013093, versus the boundary-free reference 0.016904: a 29% cost. Both scale linearly with noise amplitude.

Table 1. Oracle Gaussian influence RMSE under actual OU noise at 5%. Both methods use the same test span and clean design; they differ only in the covariance used for GLS weights. Actual variance uses the OU covariance for both.

| Prescribed windows | Retained tests | OU-weight RMSE | iid-weight RMSE | OU RMSE / class reference |
|---:|---:|---:|---:|---:|
| 30 | 30 | 0.019605 | 0.020296 | 1.1597 |
| 120 | 120 | 0.017043 | 0.017839 | 1.0082 |
| 480 | 479 | 0.016938 | 0.017734 | 1.0020 |

Rich spans produce most of the improvement: 480-window OU GLS is within 0.2% of the boundary-free local reference, while iid weights retain a roughly 5% gap. Exact residual-covariance calculations verify the stated local geometry. The original 128-path finite-span check is in S25.2; the larger paired comparison follows.

Table 2. Van der Pol on the same 2000 paths at 0.5% noise. Error norms are divided by the true four-coefficient norm, two. The reference concerns the leading local Gaussian influence. All paths are retained.

| Estimator | RMSE (95% interval) | Bias | Dispersion | Local reference |
|---|---:|---:|---:|---:|
| One step, known initial | 0.001301 [0.001269, 0.001332] | 0.000085 | 0.001298 | 0.001309 |
| One step, unknown initial | 0.001838 [0.001790, 0.001887] | 0.000469 | 0.001777 | 0.001690 |
| Three steps, unknown initial | 0.001669 [0.001630, 0.001708] | 0.000068 | 0.001668 | 0.001690 |
| Weak GLS, OU | 0.001685 [0.001645, 0.001723] | 0.000083 | 0.001683 | 0.001690 |
| Weak GLS, iid | 0.001759 [0.001718, 0.001801] | 0.000103 | 0.001756 | 0.001690 |

All 2000 prespecified paths at 0.5% noise are retained, including the original 128-path prefix once. Feasible weak GLS uses a bounded 30-window covariance pilot and 480 prescribed windows. Likelihood refinements use exact OU innovations and analytical ODE sensitivities. Paired iid-minus-OU RMSE is 0.00007426, with bootstrap interval [0.00005576, 0.00009327], a 4.22% benefit from temporal weighting. The initial unknown-state scoring step has an appreciable nonlinear initialization bias; three steps approach the nuisance-adjusted reference. The known-initial estimate has lower error because it observes the information discarded by compact tests.

![Figure 2. Van der Pol on 2000 shared paths at 0.5% noise: one scoring step with known initial state, one or three joint steps with unknown initial state, and weak GLS with OU or iid weights. Whiskers give 95% bootstrap intervals for RMSE over independent paths, not parameter confidence intervals or error standard deviations. The known-initial information reference is lower than the unknown-initial reference. The finite-step unknown-initial estimates are identified separately; reference agreement concerns local small-noise behavior.](figures/publication/likelihood_comparison.png)

### 7.2. Recovery with universal endpoint tests

For known $u$, allow nonzero test traces and use the corrected residual

$$
m_\phi(\theta)=-\int_0^T[\phi'(t)^TY(t)+\phi(t)^Tf_\theta^H(Y(t))]dt
+\phi(T)^TY(T)-\phi(0)^Tu.
$$

The universal family adds twelve time-only endpoint half-windows and $Y(0)-u$ controls to 480 interior windows. Covariance weights still require an integrated 30-window pilot. Theorem 2 identifies the full local influence class; the known-initial trace adaptation itself is established by Brunel et al. [@brunel2014, Section 2.3].

Table 3. Known-initial endpoint confirmation on 1968 paired Van der Pol paths at 0.5% noise. The likelihood reference is one Fisher-scoring step. All errors use the four-coefficient norm, two.

| Method | RMSE | Bias norm | Dispersion |
|---|---:|---:|---:|
| Boundary-free, 480 | 0.00168754 | 0.00008679 | 0.00168531 |
| Universal endpoints | 0.00130936 | 0.00007049 | 0.00130746 |
| Model-adjoint endpoints | 0.00129741 | 0.00006771 | 0.00129564 |
| Known-initial scoring | 0.00129800 | 0.00008753 | 0.00129505 |

These 1968 records were unused for endpoint fitting, while their compact and scoring references were archived. The frozen confirmation criterion requires at least 50% gap recovery. Universal endpoints recover 97.1%, with paired interval 95.5–98.6%, and lie 0.88% above scoring. Model-adjoint tests realize score directions and provide an implementation check, rather than a separate methodological gain. S22 and S25.5 record the pre-confirmation criterion amendment and quadrature checks.

![Figure 3. Known-initial information recovery on 1968 paired Van der Pol records at 0.5% OU noise. Whiskers are Monte Carlo 95% RMSE intervals from shared record bootstrap resamples. Universal time-only endpoint windows recover 97.1% of the original weak-to-scoring gap; model-adjoint tests serve as a score-realization check. The reference is a single known-initial Fisher-scoring step.](figures/publication/endpoint_confirmation.png)

### 7.3. Attribution to initial traces

The 500-path ablation holds Simpson quadrature and the covariance pilot fixed. Component gains are measured against the matched compact branch and known-initial scoring.

Table 4. Endpoint-component ablation on 500 paired paths at 0.5% noise. Recovery is measured against the matched Simpson boundary-free branch and known-initial scoring. Intervals resample common record indices.

| Added component | RMSE | Gap closed | Paired 95% interval |
|---|---:|---:|---:|
| Initial half-windows | 0.00132445 | 98.9% | [95.1%; 102.7%] |
| Terminal half-windows | 0.00164493 | -0.3% | [-1.8%; 1.2%] |
| Initial controls only | 0.00159795 | 14.3% | [7.0%; 21.4%] |
| All components | 0.00132453 | 98.9% | [95.1%; 102.6%] |
| Matched boundary-free | 0.00164395 | 0% | — |
| Known-initial scoring | 0.00132100 | 100% | — |

Initial traces alone recover essentially the entire gap, quantifying the Brunel mechanism under OU noise. Initial controls give a smaller standalone gain and terminal traces give none resolved. Combining all components changes RMSE relative to initial traces by $7.46\times10^{-8}$, with paired interval $[-7.18\times10^{-7},8.55\times10^{-7}]$. Their contributions overlap and cannot be added. The control's standalone improvement is specific to correlation between initial and interior error, rather than evidence for an additional gain after trace adaptation.

![Figure 4. Paired endpoint-component attribution on 500 Van der Pol records. Whiskers show 95% bootstrap intervals for the fraction of the matched weak-to-scoring RMSE gap recovered. The dashed line marks the prespecified 50% material-contribution rule. Initial traces recover the gap; initial-noise controls give a smaller standalone gain, and terminal traces give no resolved improvement. Contributions overlap.](figures/publication/endpoint_ablation.png)

### 7.4. Moderate-noise limits and Lorenz confirmation

The 5% Van der Pol check reuses all 128 matched observations, with known initial state and the same feasible pilot. The likelihood reference optimizes exact OU innovations.

Table 5. Moderate-noise endpoint confirmation on 128 paired Van der Pol paths at 5%. The known-initial reference uses full nonlinear likelihood fitting.

| Method | RMSE | Bias norm | Dispersion |
|---|---:|---:|---:|
| Matched boundary-free | 0.0255965 | 0.0081646 | 0.0242594 |
| Universal endpoints | 0.0218790 | 0.0056954 | 0.0211247 |
| Known-initial MLE | 0.0142088 | 0.0007255 | 0.0141902 |
| Unknown-initial MLE | 0.0171952 | 0.0006504 | 0.0171829 |

Universal endpoints recover 32.6% of the gap (18.3–45.8%), failing the frozen 50% criterion. Clean-design first-order and full-Hermite references are 0.013094 and 0.013636, below empirical dispersion 0.021125. Restored local geometry therefore leaves finite-noise design, pilot and fitting effects. The matched author-WENDy comparison is in S25.3.

For Lorenz, we fit $(s,r,b)=(10,28,8/3)$ in its supplied physical structure, with $u=(1,1,1)$, $T=4$, $h=0.001$, OU range 0.1 and 0.5% coordinate noise. The unknown-initial local relative RMSE reference is 0.000411, versus 0.000209 when the initial state is known: almost a factor two.

Table 6. Lorenz endpoint confirmation on 256 paired paths. Errors are relative to the norm of $(10,28,8/3)$. The common horizon is four and the sampling interval is 0.001.

| Method | RMSE | Bias norm | Dispersion |
|---|---:|---:|---:|
| Matched boundary-free | 0.00040994 | 0.00002170 | 0.00040936 |
| Universal endpoints | 0.00021686 | 0.00002153 | 0.00021579 |
| Known-initial MLE | 0.00020784 | 0.00002355 | 0.00020650 |

All 256 paths and all converged fits are retained. Universal endpoint windows recover 95.5% of the compact-to-likelihood gap (92.8–98.0%), passing the same 50% rule. The full-Hermite clean-design reference is 0.00020926, close to the observed endpoint dispersion. These two supplied structures support small-noise recovery; moderate-noise universality remains unestablished.

### 7.5. Prespecified testing versus sparse inclusion

A negative $\delta x^3$ perturbation to the second Van der Pol equation leaves the orbit bounded. Let $S_{-\delta}$ collect the other coefficient sensitivities. The local composite likelihood information is $\|s_\delta-\Pi_{S_{-\delta}}s_\delta\|_{\mathcal H}^2$; compact tests additionally project out $Z$. Profiling all 17 other library coefficients gives two-point reference magnitudes 0.01928 and 0.01941 at 0.5% noise. The initial-state cost for this direction is only 0.68%, despite the larger aggregate estimation cost. These references use equal-prior testing error 0.05, not the Wald significance level.

The targeted Wald procedure fits all eighteen degree-one-through-three coefficients with Hermite-corrected weak GLS. A bounded core-five pilot supplies covariance and the quadrature-defect correction. Its structural family matches the generating model, whereas a general library selector lacks that supplied structure. The author's MSTLS and observation paths are unchanged; Research archive S20 and supplement S23 give the pinned pipeline and intermediate-grid protocol.

Table 7. Targeted Wald rejection and author MSTLS inclusion, each based on 100 paired paths per cell. Magnitudes are multiples of the null local full-library boundary-free two-point reference, 0.01941 at 0.5% noise and 0.19410 at 5%. Alternatives have negative signs. The 5% Wald column fails the null-calibration gate and is descriptive.

| Magnitude | Wald (0.5%) | MSTLS (0.5%) | Wald (5%) | MSTLS (5%) |
|---:|---:|---:|---:|---:|
| 0 | 7% | 0% | 23% | 2% |
| 0.5 | 44% | 0% | 62% | 4% |
| 1 | 91% | 0% | 90% | 24% |
| 2 | 100% | 0% | 99% | 58% |
| 2.5 | 100% | 0% | 100% | 68% |
| 3 | 100% | 46% | 99% | 69% |
| 3.5 | 100% | 97% | 98% | 78% |
| 4 | 100% | 100% | 95% | 70% |

At small noise, nominal null rejection is 7/100 (Wilson 3.4–13.7%), compatible with 5%. At one reference unit the targeted test rejects on 91/100 paths, versus no MSTLS inclusions, consistent with the local Gaussian Wald prediction. MSTLS reaches its transition at larger magnitudes. This contrasts supplied-candidate testing with a particular sparse selector; it neither attributes the difference exclusively to its threshold nor establishes library-wide or post-selection validity.

An independently frozen confirmation uses 2000 new null paths and 1000 new paths at each of two alternatives, with the estimator unchanged. Null rejection is 5.25%, with exact 95% interval 4.31–6.32%, wholly inside the prespecified 3–7% band. At one reference unit, 92.9% rejection (91.13–94.41%) passes the tolerance of five percentage points around the 90.82% Gaussian prediction. At half a unit, 41.0% (37.93–44.12%) fails that interval-containment criterion around 37.65%, although the point difference is only 3.35 points. Thus size passes the frozen validation, but the complete local power curve does not; exact agreement with the Gaussian approximation is not established. S29.1 gives the protocol and intervals.

![Figure 5. Paired rejection/inclusion frequencies for the added negative cubic term, with Wilson 95% intervals and 100 paths per point. The dashed curve is the two-sided Gaussian Wald reference computed from null local efficient information; the horizontal scale uses a distinct equal-prior two-point error reference. At one unit its Wald prediction is about 90.8%. The small-noise null-calibration gate passes; the 5% gate fails, so the right-panel Wald frequencies are descriptive. Finite alternatives need not follow the null-local curve.](figures/publication/support_wald_mstls.png)

### 7.6. Moderate-noise standard-error diagnostic

At 5% noise, the nominal Wald null rate is 23/100 (Wilson 15.8–32.2%), so the corresponding alternative rejection frequencies remain descriptive. A separately frozen diagnostic on these same 100 null paths holds coefficient estimates and critical values fixed. A random-design score sandwich gives 2/100 rejections (0.6–7.0%), and a full-refit parametric bootstrap gives 5/100 (2.2–11.2%); both pass the Wilson-compatibility gate. Median paired SE multipliers are 2.81 and 1.79. S24 reports the complete counts, 399 bootstrap refits per path, 999 model score draws and numerical audit. These previously examined paths motivate independent size confirmation; power with the revised SEs remains untested.

### 7.7. Matching the constructive estimator to its finite-span reference

The full-library estimator of Theorem 3 uses nine smoothed-data feature instruments and 125 time-only cubic splines. The deterministic limiting covariance for this 134-test span has $\sqrt{\operatorname{tr}V_k}=9.8144$. The calculation retains all rows and uses neither ridge nor rank truncation. Refining the finest grid from $h=0.00015625$ to $0.000078125$ changes this reference by $4.7\times10^{-8}$ relatively. It concerns all eighteen coefficients, unlike the earlier supplied four-coefficient Van der Pol benchmarks.

Table 8. Constructive full-library estimator on 100 independent paths at each amplitude. Intervals resample archived coefficient errors 10,000 times and are pointwise; no new observations are generated for this comparison.

| Noise | RMSE/$\epsilon$ | RMSE/$(\epsilon\sqrt{\operatorname{tr}V_k})$ | 95% interval for ratio |
|---:|---:|---:|---:|
| 5% | 11.6265 | 1.1846 | 1.1140–1.2563 |
| 0.5% | 10.3505 | 1.0546 | 0.9883–1.1202 |
| 0.1% | 9.9850 | 1.0174 | 0.9452–1.0885 |
| 0.05% | 9.5899 | 0.9771 | 0.8994–1.0528 |
| 0.01% | 9.3516 | 0.9528 | 0.8891–1.0147 |

All four lower-noise ratio intervals include one; the 5% interval lies above it. There are no guard fallbacks or final clipping on these 500 paths. This supports the fixed-family small-noise covariance prediction and identifies finite-noise excess risk. It neither tests the sequential dense-span limit nor proves exact Wald calibration. Supplement S29.2 records matrix checks and covariance diagnostics; Table S27 distinguishes the three estimators used in this paper.

## 8. Discussion and limitations

Boundary-free equations retain exactly the local information available when the initial state is a nuisance parameter. Their restriction has a quantified cost when the state is known: 29% in the Van der Pol aggregate RMSE reference and almost a factor two for Lorenz. Rich spans and correct temporal weights approximate the sharp class reference. Theorem 2 supplies the complementary full-information class for known-initial boundary adaptations. The exact adjoint-range equalities, rather than classical projection formulas or the Brunel remedy, are the theoretical contribution.

The recovery experiment and ablation distinguish class theory from method novelty. Universal endpoint windows recover 97.1% of the weak-to-scoring gap on 1968 Van der Pol records and 95.5% of the weak-to-likelihood gap on 256 Lorenz records. On the 500-record Van der Pol ablation, initial traces alone explain essentially all that gain, reproducing the established Brunel adaptation. The correlated-noise initial control is useful separately, recovering 14.3% of the gap, but contributes no resolved improvement after initial traces. Thus the remedy is an established boundary mechanism with a quantitative OU benchmark. At 5% Van der Pol noise only 32.6% of the likelihood gap is recovered. Covariance weights still require an integrated pilot model; computational superiority has not been demonstrated.

Information loss is direction specific. For the added cubic coefficient, initial-state projection changes the local composite detection reference by only 0.68%, despite its larger aggregate coefficient-estimation cost. On the refined small-noise grid, the supplied MSTLS includes the term on 0%, 46% and 97% of paths at magnitudes 2.5, 3 and 3.5 reference units. Independent small-noise confirmation gives a 5.25% null rate and passes the frozen size band. Power passes its tolerance criterion at one reference unit and fails it at half a unit, delimiting agreement with the local approximation (S29.1). Its pilot information and target differ from the structure selector; the frequencies alone establish no general statistical superiority. This comparison identifies a practical gap between supplied-candidate testing and the particular sparse-selection procedure. It does not prove a threshold-only cause, identify a variance-loss ratio or provide inference after unrestricted library selection. At 5% noise the nominal Wald null rate is 23%. A separately frozen same-path diagnostic reduces it to 5% using a full-refit parametric bootstrap SE and to 2% using a random-design sandwich; both Wilson intervals include nominal 5%. This identifies estimator-level uncertainty as a practical calibration issue, while independent size confirmation and moderate-noise power remain unresolved. Large finite perturbations also leave the null local Gaussian approximation.

The scope is fully observed deterministic polynomial ODEs with known positive sensor covariance. Small-noise influence efficiency, finite-noise feasible accuracy and support selection are separate claims. Exact residual covariance does not by itself give exact finite-noise estimator variance when weights and designs depend on the record. On the matched moderate-noise benchmark, likelihood is more accurate than the original boundary-free weak estimator. Independent-record GMM, fission and cross-fitting remain supplementary extensions. Joint span/noise rates, unknown-covariance calibration with heterogeneous states, partial observation and process noise require further analysis.


## Data availability

The preceding benchmark is publicly archived in the [reproducibility repository, version v1.2.0](https://github.com/malyshevdanil/wsindy-information-limits/tree/v1.2.0) and its [Zenodo snapshot](https://doi.org/10.5281/zenodo.23313935). The present endpoint, selection and SE-calibration follow-ups are preserved separately in the accompanying reproducibility bundle; the v1.2.0 DOI identifies the preceding version. Configurations, seed ranges, source fingerprints, coefficient estimates, covariance summaries and SHA256 manifests permit exact input reconstruction and paired audit. Full synthetic time series are regenerated from the frozen seed protocol. Pinned download specifications retrieve original author sources. No external experimental dataset is used. A public versioned identifier for the current bundle must be added before submission.

## References
