---
title: "Boundary-free weak influences: range closure, sharp information and feasible GLS"
lang: en-US
---

## 1. Spaces and adjoint range

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

## 2. Sharp covariance reference

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

## 3. Explicit OU representative

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

## 4. GLS and the order of limits

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

## References
