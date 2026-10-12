# Fixed-noise weak GMM across independent records

## S14. Common-trajectory GMM reference and proof

The replicate asymptotic is distinct from single-path infill and from Corollary 1. Let $Y_r=x_{\theta_0,u}+\eta_r$, $r=1,\ldots,n$, be independent copies on the same fixed interval. The deterministic trajectory and initial state are shared, the covariance is stationary Gaussian, and its amplitude is any fixed positive finite value. Assume a supplied polynomial vector field linear in $\theta\in\mathbb R^p$, a fixed collection of boundary-vanishing tests, and exact integrals. Write the Hermite-corrected weak moment as

$$
m_r(\theta,\nu)=b_r-G_r(\nu)\theta\in\mathbb R^M,
\qquad D=E G_r(\nu_0),\qquad E m_r(\theta_0,\nu_0)=0,
$$

where $\nu$ contains the marginal variances used in the correction and any covariance-shape parameters. The expectation identity follows by the Gaussian heat/Hermite identity at every time, Fubini, and the deterministic weak ODE identity. It is exact at fixed noise; it does not assert that the ratio or least-squares estimator from one record is unbiased.

**Lemma 2 (fixed-noise independent-record GMM).** Suppose $D$ has column rank $p$, the covariance $\Omega=\operatorname{Cov}(m_r(\theta_0,\nu_0))$ is positive definite, and moments have finite second moments. With known $\nu_0$ and random positive-definite $W_n\to_p W$, define

$$
\widehat\theta_n=(\bar G_n^\top W_n\bar G_n)^{-1}\bar G_n^\top W_n\bar b_n,
\qquad P_W=(D^\top W D)^{-1}D^\top W.
$$

Then $\widehat\theta_n\to_p\theta_0$ and

$$
\sqrt n(\widehat\theta_n-\theta_0)\Rightarrow N(0,V_W),
\qquad V_W=P_W\Omega P_W^\top.
$$

Within this fixed moment family and the corresponding regular linear GMM influence class, $W=\Omega^{-1}$ minimizes covariance, giving $V_*=(D^\top\Omega^{-1}D)^{-1}$. For polynomial moments, $\Omega$ is the full covariance from Proposition 3, including every nonzero Hermite order. Empirical residual covariance at a consistent preliminary estimate, or a consistent full model-based covariance, supplies feasible two-step weights. First-order covariance alone generally gives a different, nonoptimal $V_W$ at fixed noise. The theorem does not assert full-path likelihood efficiency at finite amplitude.

For estimated $\nu$, additionally assume local differentiability and an asymptotically linear nuisance estimator

$$
\sqrt n(\widehat\nu-\nu_0)=n^{-1/2}\sum_{r=1}^n\xi_r+o_p(1),
\quad E\xi_r=0,\quad E\|\xi_r\|^2<\infty,
\quad H=E\partial_\nu m_r(\theta_0,\nu_0).
$$

Assume the corresponding sample derivative converges locally, with a uniform second-order Taylor remainder $o_p(n^{-1/2})$. Then the same closed-form estimator, evaluated at $\widehat\nu$, is consistent and asymptotically normal with

$$
\psi_r=m_r(\theta_0,\nu_0)+H\xi_r,
\quad \widetilde\Omega=\operatorname{Cov}(\psi_r),
\quad \widetilde V_W=P_W\widetilde\Omega P_W^\top.
$$

If $\widetilde\Omega$ is positive definite, estimated inverse $\widetilde\Omega$ is optimal **for this specified nuisance estimator and moment influence family**. Its optimum is not generally the known-covariance optimum or the globally optimal joint GMM estimator. Estimating only a shape parameter that changes weights but not mean correction has zero $H$ for that parameter; estimating a variance used in Hermite correction need not. A consistent empirical influence covariance, including cross-covariances between $m_r$ and $\xi_r$, yields marginal Wald intervals and joint ellipsoids with asymptotic nominal coverage. These are standard GMM and two-step-estimation arguments specialized to the weak Hermite moments [@hansen1982; @newey1994].

*Proof.* Polynomial Gaussian functionals on a bounded interval with continuous deterministic mean have moments of all finite orders; bounded test functions allow the required exchanges of integration and expectation. For known covariance, the law of large numbers gives $\bar G_n\to D$, $\bar b_n\to D\theta_0$, and convergence of the inverse normal matrix. The exact identity

$$
\widehat\theta_n-\theta_0=P_n\bar m_n(\theta_0,\nu_0)
$$

and the iid multivariate CLT give the stated limit, without a small-noise expansion. If $P D=I_p$, set $P_*= (D^\top\Omega^{-1}D)^{-1}D^\top\Omega^{-1}$. Then

$$
P\Omega P^\top-P_*\Omega P_*^\top
=(P-P_*)\Omega(P-P_*)^\top\succeq0,
$$

because the cross terms vanish. A consistent preliminary least-squares estimate makes the empirical covariance converge to $\Omega$ by the law of large numbers and polynomial continuity. In the nuisance case, the stipulated Taylor expansion gives

$$
\sqrt n\bar m_n(\theta_0,\widehat\nu)
=n^{-1/2}\sum_r[m_r(\theta_0,\nu_0)+H\xi_r]+o_p(1).
$$

The random design and random weights multiply a $O_p(n^{-1/2})$ centered mean, so their consistency suffices for the leading map. Apply the joint CLT to $(m_r,\xi_r)$ and the same covariance decomposition using $\widetilde\Omega$. Consistency of the sandwich and coverage follow by Slutsky. This proves the assertions.

**Numerical quadrature condition.** On a fixed discrete grid, Hermite correction is still unbiased for the corresponding *quadrature moment*. Its deterministic expectation is the weak quadrature defect, which is generally nonzero at $\theta_0$. Consistency for the continuous ODE therefore requires exact integration or vanishing defect. In a triangular array, require $\sqrt n\|E m_{r,h_n}(\theta_0,\nu_0)\|\to0$, convergence of designs and covariance, and uniform integrability sufficient for the triangular CLT and nuisance remainder. Mean-square convergence of quadrature moments to their continuous versions gives the stochastic part for this fixed test family. This condition concerns grid refinement with independent records; it does not prove a joint growing-window/small-noise rate. The numerical experiment uses a fixed fine grid and discloses its clean defect rather than claiming literal fixed-grid consistency for the ODE parameter.

### S14.1. Consistent OU shape from independent copies

For even $n=2B$, define independent pairs $d_r=(Y_{2r-1}-Y_{2r})/\sqrt2$. Their common deterministic mean cancels and $d_r$ has the original stationary OU covariance. At a fixed nondegenerate grid $t_0,\ldots,t_L$ with mesh $h>0$, let

$$
s_{rj}=\frac1{L+1}\sum_{i=0}^L d_{rj}(t_i)^2,
\quad a_r=\sum_{j,i<L}d_{rj}(t_i)d_{rj}(t_{i+1}),
\quad c_r=\sum_{j,i<L}d_{rj}(t_i)^2.
$$

Then $\widehat v_j=B^{-1}\sum_r s_{rj}\to v_j$, $\widehat\rho=\bar a/\bar c\to e^{-h/\tau}$, and $\widehat\tau=-h/\log\widehat\rho\to\tau$. Their joint $\sqrt B$ normality follows from the independent-pair CLT and the delta method; dependence of time points is retained inside each pair. Positivity of variances and $0<\rho<1$ make the transformation regular. The influence for $v_j$ is $s_{rj}-v_j$, and for $\rho$ it is $(a_r-\rho c_r)/E c_r$; the $\tau$ influence is $h\xi_{\rho,r}/[\rho(\log\rho)^2]$. If finite-sample clipping is needed, any fixed bounds containing the true parameters in their interior are asymptotically inactive and must be reported.

Use pair-averaged weak moments for inference, so the theorem applies with sample size $B$, independent pair influences, and variance divided by $B=n/2$. In Van der Pol the only correction derivative is $H_{v_1,2k}=\theta_4 E\int\phi_k Y_2$; the other rows and variance/shape derivatives are zero. Thus the empirical covariance of pair residuals plus this derivative times $(s_{r1}-\widehat v_1)$ supplies the nuisance-adjusted weights and sandwich. Averaging pair differences for shape estimation requires a common deterministic mean; heterogeneous initial conditions would require a different construction. Replication permits separate OU variance and range consistency that is unavailable under fixed-horizon infill of one observed path.

### S14.2. Why stacking can retain bias

Stacked corrected least squares converges, when its normal matrix is nonsingular, to

$$
\theta_{\mathrm{stack}}=\theta_0+
\{E(G_r^\top G_r)\}^{-1}E[G_r^\top m_r(\theta_0,\nu_0)].
$$

Exact $E m_r=0$ does not imply $E[G_r^\top m_r]=0$, since the design and residual share measurement error. Stacking is therefore generally biased and inconsistent in this errors-in-variables setting; special orthogonality exceptions are possible. Averaged-moment GMM instead uses $\bar G_n^\top W_n\bar m_n$, whose probability limit is zero at the true parameter. The stacking behavior in Table S11 is an instance of this distinction, not a claim that every stacked regression must be inconsistent.

With known covariance, the same argument extends to iid random initial states sharing the coefficient vector, provided the conditional weak identity, full-rank mean design and finite moments hold. This covers the mechanism in the heterogeneous-initial-state aggregation example. The pair-difference nuisance estimator above still requires a common deterministic mean and is not applied to that heterogeneous example.

### S14.3. Misspecified temporal weighting

For known marginal correction and any nonsingular assumed temporal shape, the fixed-family variance is the sandwich $P_W\Omega_{\mathrm{true}}P_W^\top/n$. In the small-noise single-record setting replace $\Omega_{\mathrm{true}}$ with the order-one covariance and the factor $1/n$ with $\epsilon^2$. Correct marginal correction preserves averaged-moment consistency even if $\tau$ is wrong; a wrong variance used for Hermite correction can destroy the moment identity. The finite-noise fitted-shape rows in S12 change both correction and weights and cannot be read as a pure temporal-weighting experiment. The $\tau/2$ and $2\tau$ rows retain the true marginal variances.

## References
