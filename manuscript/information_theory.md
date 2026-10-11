# Fixed-horizon information bounds for deterministic dynamics in OU measurement error

The results below specialize classical Gaussian shift theory to the observation model. Equivalence-based impossibility is established methodology, including Zhang (2004, Section 3); the contribution here is an explicit dynamical sensitivity norm and its link to weak residuals and acquisition choices.

Continuous-record precision and efficient discrete designs also have precedents in correlated regression [@dette2016]. Local ODE sensitivities serve as regressors in that setting.

## Observation experiment and Cameron–Martin norm

Observe $Y(t)=\mu_\theta(t)+\varepsilon(t)$ on $[0,T]$, where the deterministic mean is in $H^1([0,T];\mathbb R^d)$. Every parameter uses the same known stationary OU covariance. Coordinates have positive variances $\sigma_j^2$ and correlation times $\tau_j>0$, and are independent. The RKHS is $H^1$ as a set, with norm

$$
\|g\|_{\mathcal H}^2=\sum_{j=1}^d\left[\frac{g_j(0)^2}{\sigma_j^2}+\frac{\tau_j}{2\sigma_j^2}\int_0^T\left(g'_j+\frac{g_j}{\tau_j}\right)^2dt\right].
$$

To verify this directly, write $d\varepsilon_j=-\varepsilon_jdt/\tau_j+\sqrt{2\sigma_j^2/\tau_j}\,dW_j$, with an independent stationary Gaussian initial value. The invertible path transformation $z\mapsto(z(0),z(t)-z(0)+\tau_j^{-1}\int_0^tz(s)ds)$ maps the experiment to its initial normal variable and scaled Brownian path. Translation by $g$ costs the initial normal shift plus Brownian drift $g'+g/\tau_j$. The Cameron–Martin theorem gives the stated norm. Expanding its square gives the symmetric form

$$
\|g\|_{\mathcal H}^2=\sum_j\left[\frac{\tau_j}{2\sigma_j^2}\int_0^T(g'_j)^2dt+\frac1{2\sigma_j^2\tau_j}\int_0^Tg_j^2dt+\frac{g_j(0)^2+g_j(T)^2}{2\sigma_j^2}\right].
$$

The endpoint terms are essential for stationary initial noise. Assuming zero initial noise would define a different experiment.

## Theorem: impossibility for an entire parameter family

For any two distinct parameters $\theta,\eta$ with means in $H^1$, their continuous-path laws are equivalent. Let $\delta=\|\theta-\eta\|$ and $d_{\theta\eta}=\|\mu_\theta-\mu_\eta\|_{\mathcal H}$. Every estimator $\widehat\theta_h$ based on any finite collection of evaluations of the same path satisfies

$$
\max_{v\in\{\theta,\eta\}}P_v\left(\|\widehat\theta_h-v\|\geq\delta/2\right)\geq\Phi(-d_{\theta\eta}/2)>0.
$$

Consequently no sequence of denser-grid estimators is consistent in probability at every parameter of a nontrivial family. This assertion permits an estimator that is consistent at an isolated parameter, for example a constant estimator.

**Proof.** A shift between the means belongs to the common RKHS, so the Cameron–Martin translation theorem gives mutual absolute continuity. Under the first law the log likelihood ratio of the second law is normal with mean $-d_{\theta\eta}^2/2$ and variance $d_{\theta\eta}^2$. Thus the total variation distance of the path laws is $2\Phi(d_{\theta\eta}/2)-1$. Evaluation of a path cannot increase total variation. Any estimator defines a nearest-parameter test. A test error can occur only when the estimator's error relative to the generating parameter is at least $\delta/2$. The sum of errors of every test is at least one minus total variation, or $2\Phi(-d_{\theta\eta}/2)$. At least one of the two estimation error probabilities is at least half this sum. The lower bound is independent of the grid, which excludes simultaneous consistency. The same reasoning and $E[Z^2]\geq r^2P(|Z|\geq r)$ yield

$$
\max_{v\in\{\theta,\eta\}}E_v\|\widehat\theta_h-v\|^2\geq\frac{\delta^2}{4}\Phi(-d_{\theta\eta}/2).
$$

This risk statement allows infinite moments and needs neither unbiasedness nor a particular weak estimator. It also covers randomized estimators by adjoining independent randomization to both experiments. For $d_{\theta\eta}=0$, both means define the same law and the bound remains valid. $\square$

## Fisher information and weak efficiency

Suppose $\theta\mapsto\mu_\theta$ is differentiable in the RKHS norm. For sensitivities $s_r=\partial_{\theta_r}\mu_\theta$, the continuous-record Fisher information is

$$
\mathcal I_{rs}(\theta)=\langle s_r,s_s\rangle_{\mathcal H}.
$$

For a regular unbiased estimator of a scalar parameter, $\operatorname{Var}(\widehat\theta)\geq\mathcal I(\theta)^{-1}$. With nuisance parameters the appropriate bound uses the corresponding element of the inverse information matrix, when nonsingular. Finite information does not imply that a finite-horizon maximum likelihood estimator attains this bound.

For $\mu_a(t)=x_0e^{at}$ with known $x_0$, sensitivity $s(t)=x_0te^{at}$ and $s(0)=0$, one obtains

$$
\mathcal I_T(a)=\frac{\tau x_0^2}{2\sigma^2}\int_0^T e^{2at}\left[1+(a+\tau^{-1})t\right]^2dt.
$$

At a uniform grid, $\rho=e^{-h/\tau}$, the exact information is

$$
\mathcal I_h(a)=\frac{s(0)^2}{\sigma^2}+\frac{\sum_{n=1}^N[s(t_n)-\rho s(t_{n-1})]^2}{\sigma^2(1-\rho^2)}\ \longrightarrow\ \mathcal I_T(a).
$$

For nested grids it increases monotonically by data processing. The convergence also follows directly from $s(t_n)-\rho s(t_{n-1})=h(s'+s/\tau)+o(h)$ and $1-\rho^2=2h/\tau+o(h)$.

For the scalar weak ratio, define $D_\phi=\int\phi\mu_a\ne0$. Its first-order noise influence is $\ell(t)=(-\phi'(t)-a\phi(t))/D_\phi$. Differentiating the weak identity gives $\int\ell s=1$. Let $V_\phi=\iint\ell(u)C(u,v)\ell(v)dudv$. The RKHS representer of this linear functional is $C\ell$, so Cauchy–Schwarz gives

$$
1=\langle C\ell,s\rangle_{\mathcal H}^2\leq V_\phi\mathcal I_T(a),\qquad \operatorname{Eff}_\phi=\frac1{\mathcal I_T(a)V_\phi}\leq1.
$$

Here the squared inner product equals one; the preceding display uses that normalization. $V_\phi$ is the variance of a Gaussian linearization. The unguarded ratio itself can have infinite second moment, so its empirical squared errors are not estimators of $V_\phi$ and must not be described as an attained Cramér–Rao variance.

## Initial-state annihilation and a stricter weak influence bound

Let $S=\partial_\theta x$ and $Z=\partial_u x$ for active vector-field coefficients and initial state $u$. Their variational equations are $S'=JS+F$ and $Z'=JZ$, where $J=\partial_xf$ and $F=\partial_\theta f$. For compactly supported test functions, the first-order joint residual influence is $A_k(t)=-\phi'_k(t)I-\phi_k(t)J(t)$. Integration by parts yields $\int A_kZ=0$ and $\int A_kS=\int\phi_kF$. With a full-rank clean weak design and a left inverse $P_0$, the normalized coefficient influence $\ell=P_0A$ therefore satisfies

$$
\int\ell Z=0,\qquad\int\ell S=I.
$$

Write the Gaussian sensitivity score vectors as $U_\theta,U_u$. They have joint information matrix given by the RKHS Gram matrix of $S,Z$. When its blocks and Schur complement are positive definite, the efficient score is

$$
U_{\mathrm{eff}}=U_\theta-\mathcal I_{\theta u}\mathcal I_{uu}^{-1}U_u,
\qquad
\mathcal I_{\mathrm{eff}}=\mathcal I_{\theta\theta}-\mathcal I_{\theta u}\mathcal I_{uu}^{-1}\mathcal I_{u\theta}.
$$

For $L_\varepsilon=\int\ell\varepsilon$, the score covariance identity gives $\operatorname{Cov}(L_\varepsilon,U_{\mathrm{eff}})=I$. Thus

$$
\operatorname{Cov}\big(L_\varepsilon-\mathcal I_{\mathrm{eff}}^{-1}U_{\mathrm{eff}}\big)
=\operatorname{Cov}(L_\varepsilon)-\mathcal I_{\mathrm{eff}}^{-1}\succeq0.
$$

The Schur-complement covariance bound is classical Gaussian regression geometry [@dette2016]. The weak restriction is its mandatory annihilation of initial-state sensitivities, even if the actual initial state is supplied. This bound concerns the Gaussian influence, not the unguarded nonlinear estimator's moments. It is not an attainment theorem for a chosen finite window family.

For $x=x_0e^{at}$, initial-state sensitivity is $q=e^{at}$, while coefficient sensitivity is $s=x_0te^{at}$. The projected information is $\mathcal I_{aa}-\langle s,q\rangle_{\mathcal H}^2/\|q\|_{\mathcal H}^2$. Independent deterministic numerical verification is recorded in the geometry execution manifest.

## Acquisition regimes and heterogeneous recordings

For a constant mean $\mu_\theta(t)=\theta$, the exact continuous information is $\sigma^{-2}[1+T/(2\tau)]$. The quantity $1+T/(2\tau)$ is its information-equivalent number of independent observations of variance $\sigma^2$. This identifies the often-used $T/\tau$ scale, including its factor of two, in a specific model.

For a decaying exponential $a<0$, $\mathcal I_T(a)$ converges to a finite limit as $T\to\infty$; a longer near-zero tail does not maintain parameter excitation. At $a=0$, information instead grows as $\tau x_0^2[ T+T^2/\tau+T^3/(3\tau^2)]/(2\sigma^2)$. Therefore a universal duration rule based only on $T/\tau$ is inappropriate for dynamical coefficients.

Independent recordings add information: identical records yield $M\mathcal I_T$, while records with known differing initial conditions yield $\sum_{m=1}^M\mathcal I_{T,m}$. An acquisition decision must therefore compare model sensitivities, physical covariance, and cost. A large $\mathcal I_h/\mathcal I_T$ diagnoses how much information a grid captures in the known Gaussian model; it is neither a support-recovery guarantee nor an operational estimate when sensitivities are unknown.

For weak estimation with different initial conditions and a shared polynomial coefficient vector, every clean pair satisfies $b_m^0=G_m^0\xi$. Corrected weak matrices have conditional expectations $E[G_m\mid x_m]=G_m^0$ and $E[b_m\mid x_m]=b_m^0$ up to deterministic quadrature error. Under independent identically distributed recordings, integrability and full rank of $E[G_m^0]$ on the supplied support, the law of large numbers implies consistency of least squares applied to $\overline G_M,\overline b_M$, first for the population finite-grid identity and subsequently as quadrature error vanishes. The full-rank condition is substantive: opposing excitation can cancel when matrices are averaged. Stacking designs instead involves $E[G_m^TG_m]$ and $E[G_m^Tb_m]$, which generally retain errors-in-variables terms. Averaging sensor paths before evaluating nonlinear features is invalid for differing initial conditions because $f(E[x_m])$ need not equal $E[f(x_m)]$.

## Scope

All equivalence and information conclusions use a fixed finite horizon and common positive OU covariance across parameters. They need not hold for noiseless coordinates, parameter-dependent singular covariances, repeated independent recordings with $M\to\infty$, or increasing-domain experiments with persistent excitation. The classical measure-equivalence mechanism is attributed to its predecessors; its specialization should not be advertised as a new Gaussian-measure theorem.
