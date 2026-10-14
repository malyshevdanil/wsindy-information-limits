---
title: "Local curvature of the nuisance-profile score"
lang: en-US
---

## Incidental initial states: a local bias calculation

Increasing the number of record-specific initial states raises an incidental-parameter issue [@neyman1948; @lancaster2000]. Its presence alone does not establish a nonzero bias for a particular Gaussian mean model. Bias reduction for nonlinear panels has an extensive literature [@hahn2004; @arellano2007]; the calculation below concerns a fixed-grid Gaussian ODE mean with known covariance, and uses a small-noise expansion rather than a large within-record sample limit.

For one record whiten the observations by the known base covariance: $y=\mu(\theta_0,u)+\epsilon z$, $z\sim N(0,I)$. At the true parameters write $S=\partial_\theta\mu$, $Z=\partial_u\mu$, $N=Z^\top Z$ and $P=I-ZN^{-1}Z^\top$. Assume the mean map is three times continuously differentiable, $N$ is nonsingular, and the interior local nuisance minimizer admits a second-order expansion with integrable, uniform remainders over the allowed initial states. Define

$$
h(u)=\sum_{i,j}(N^{-1})_{ij}\partial_{u_i u_j}\mu.
$$

**Local profile-score expansion.** For the unscaled least-squares profile equation $\psi(\theta)=\partial_\theta\mu(\theta,\widehat u(\theta))^\top\{y-\mu(\theta,\widehat u(\theta))\}$,

$$
E\psi(\theta_0)=-\tfrac12\epsilon^2S^\top P h(u)+o(\epsilon^2).
$$

The log-likelihood score divides this expression by $\epsilon^2$; the order statement above concerns the unscaled equation. Suppose the average local information $\overline I_e=\lim_n n^{-1}\sum_r S_r^\top P_r S_r$ is positive definite, the population profile equation has a locally unique root, and its uniform expansion is valid. That root has

$$
\theta_\epsilon^\dagger-\theta_0
=-\tfrac12\epsilon^2\overline I_e^{-1}
\lim_n n^{-1}\sum_r S_r^\top P_r h(u_r)+o(\epsilon^2).
$$

**Proof.** Write $\widehat u=u+\epsilon a+\epsilon^2b+o_{L^2}(\epsilon^2)$. The nuisance normal equation gives $a=N^{-1}Z^\top z$. The whitened residual is $\epsilon Pz-\epsilon^2\{Zb+\tfrac12\mu_{uu}[a,a]\}+o_{L^1}(\epsilon^2)$. Gaussian orthogonality makes $a$ independent of $Pz$, and $E(aa^\top)=N^{-1}$. Expanding the nuisance normal equation to second order gives $E b=-\tfrac12N^{-1}Z^\top h$. The mixed coefficient-sensitivity term $E\{(\mu_{\theta u}a)^\top Pz\}$ is zero. Substitution into the coefficient equation therefore gives its stated mean. Its population derivative at zero noise is $-S^\top P S$, by differentiating the nuisance normal equation. The local implicit-function expansion after averaging records gives the displayed root displacement. $\square$

A nonzero curvature coefficient predicts a fixed-noise displacement that remains when the number of records increases. A consistency claim for the computed local profile estimator additionally requires a uniform law of large numbers and identification of this population root. Other nuisance branches, boundary optima and rare large-noise fits fall outside this local expansion. Linear mean models can have zero curvature term. The experiment tests the size and direction of the displacement rather than assuming a likelihood failure.

For the prescribed Van der Pol initial-state distribution, the deterministic Kaggle calculation gives the coefficient of $\epsilon^2$ as approximately $(0.0024214,-0.0037071,0.0349238,-0.0152314)$. Eight- and twelve-node Gauss--Legendre rules per initial-state coordinate agree to relative $1.8\times10^{-9}$; finite-difference checks of diagonal second initial sensitivities have relative errors below $3.3\times10^{-6}$. At $n=800$, the predicted relative bias divided by the nuisance-adjusted norm standard error is 0.0803, 0.1605 and 0.3210 at 5%, 10% and 20% noise. These are local predictions, separate from the Monte Carlo pilot and the main comparison.

## References
