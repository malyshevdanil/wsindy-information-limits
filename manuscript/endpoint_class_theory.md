### 5.5. Extended class and recovery with a known initial state

Throughout this subsection, $u$ is known exactly, covariance is supplied and independent of $\theta$, $J$ is continuous, and $S$ and $Z$ belong to the OU Cameron–Martin space $\mathcal H$. Assume $\mathcal I_{\theta\theta}$ is positive definite. The coordinate OU scales and correlation times are positive. Consequently $\mathcal E$ is the covariance-dual space of vector-valued $H^1([0,T])$, with continuous endpoint evaluation. The initial sensitivity is normalized by $Z(0)=I_d$ and the parameter sensitivity satisfies $S(0)=0$.

For $\phi\in H^1([0,T];\mathbb R^d)$, use the known-initial residual

$$
m_\phi(\theta)=-\int_0^T[\phi'^TY+\phi^Tf_\theta^H(Y)]dt
+\phi(T)^TY(T)-\phi(0)^Tu.
$$

Its first-order influence is $L\phi=(-\phi'-J^T\phi)dt+\phi(T)\delta_T$. The initial boundary term contains the exact state rather than its noisy observation. This adaptation of the weak identity is established [@brunel2014, Section 2.3]. The result below specifies its covariance-dual geometry and its relation to the correlated-noise control.

**Theorem 3 (known-initial extended influence class).** Under these assumptions,

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

If this vanishes for all permitted tests, compact interior tests first imply $g'=Jg$ in the weak sense. Such solutions are $g=Zc$. For all $H^1$ tests, or for tests with $\phi(T)=0$, the initial trace is arbitrary, so $g(0)=c=0$. Their annihilator is trivial and their closure is $\mathcal E$. For tests with $\phi(0)=0$, every $Zc$ annihilates the range and there are no further annihilators. Its closure is therefore $\mathcal M$, as characterized in Theorem 2.

For any $\lambda\in\mathcal E$, define $a_j=\langle Z_{\cdot j},\lambda\rangle$. Continuity of this pairing follows from $Z\in\mathcal H$. Since $Z(0)=I_d$, $\lambda-a\delta_0$ annihilates every initial-sensitivity column and belongs to $\mathcal M$. This proves the exact sum identity. For covariance optimality, $\langle\Lambda_{\mathrm{full}},S\rangle=I_p$ and $\operatorname{Cov}(\Lambda_{\mathrm{full}})=\mathcal I_{\theta\theta}^{-1}$. Any other normalized influence has a covariance-orthogonal remainder, giving a positive-semidefinite covariance increment. Density, continuous pairing with $S$ and normalization by the convergent sensitivity matrix give attainability as an infimum, by the same approximation argument as Theorem 2. $\square$

Correct-covariance GLS on dense universal extended spans approaches this bound under the sequential small-noise conditions of Corollary 1, with the full rather than nuisance-adjusted information. The theorem concerns the closed local influence class; performance of finite families and finite-noise covariance pilots requires numerical assessment. In particular, finite-span gains of different endpoint additions can overlap.

For compact tests, let $\Omega$ be the exact Hermite residual covariance, $\Sigma_0=\operatorname{diag}(\sigma_j^2)$ and $B=\operatorname{Cov}(m_\phi,Y(0)-u)$. The optimal linear control replaces the residual covariance by

$$
\Omega-B\Sigma_0^{-1}B^T.
$$

Only the first Hermite order contributes to $B$. In coordinate $j$, its cross-covariance is the integral of the linear residual weight against $\sigma_j^2e^{-t/\tau_j}$. In a sampled iid model the compact-test weights and derivatives vanish at the initial node; hence this cross-covariance is exactly zero. This distinction explains the correlated-noise control mechanism. It does not establish an additional finite-family gain when initial-trace tests already approximate the full efficient directions.
