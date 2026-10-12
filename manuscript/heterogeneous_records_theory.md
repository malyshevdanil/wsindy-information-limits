---
title: "Heterogeneous initial states: weak instruments and independent-record inference"
lang: en-US
---

## S16. Heterogeneous initial states and independent weak instruments

Consider independent records $Y_r=x_{\theta_0,u_r}+\eta_r$ of a deterministic polynomial ODE. The initial states may form a deterministic triangular array; alternatively they may be iid draws from a compact set, with the following conditions imposed conditionally on them. The sensor covariance is $C_\nu$. Each Hermite-corrected, exactly integrated weak moment obeys

$$
E\{m_r(\theta_0,\nu_0)\mid u_r\}=0,
\qquad m_r(\theta,\nu)=b_r-G_r(\nu)\theta.
$$

This identity follows by taking expectations of the corrected polynomial and applying the ODE weak identity separately to each trajectory. Its validity is uniform over allowed initial states; the distribution of $u_r$ need not be estimated. The numerical version needs the quadrature condition below.

### S16.1. Numbered assumptions and fixed-noise estimating equations

**(A1)** The coefficient vector has fixed dimension $p$ and belongs to an open identifiable neighborhood. Trajectories and their first coefficient/initial-state sensitivities exist on $[0,T]$. Initial states lie in a compact set on which the selected moments have uniformly bounded moments of order $2+\delta$ for some $\delta>0$.

**(A2)** Records and any auxiliary randomizations are independent across $r$. Noise is Gaussian with known positive covariance $C_{\nu_0}$. Weak tests and quadrature rules are fixed independently of each record's noise. The coefficient structure is supplied.

**(A3)** For every record there are two conditionally independent views $Y_r^+$ and $Y_r^-$ with mean $x_{\theta_0,u_r}$ and covariance $2C_{\nu_0}$. An instrument $L_r^\pm$ is measurable with respect to its own view and independent auxiliary training data. Corrected moments $m_r^\pm$ use the variance of their own view. The instrument map is fixed conditional on that training data and has the moment bounds required in (A1).

**(A4)** Put $B_r=(L_r^+G_r^-+L_r^-G_r^+)/2$, $a_r=(L_r^+b_r^-+L_r^-b_r^+)/2$ and $\psi_r=a_r-B_r\theta_0$. Conditional on the training data, $n^{-1}\sum E B_r\to B$ with $B$ nonsingular, and $n^{-1}\sum\operatorname{Var}(\psi_r)\to V_\psi$ positive definite. A law of large numbers holds for $B_r$ and the Lindeberg condition holds for $\psi_r$.

**(A5)** Exact integration is used, or the average instrumented quadrature defect is $o(n^{-1/2})$. This condition concerns centered continuous-parameter inference. A vanishing average defect alone suffices for consistency.

**Theorem 3 (heterogeneous-record weak estimating equations).** Under (A1)–(A5), conditional on the independent training data,

$$
\widehat\theta_n=\left(\sum_{r=1}^n B_r\right)^{-1}\sum_{r=1}^n a_r
\ \xrightarrow{P}\ \theta_0,
\qquad
\sqrt n(\widehat\theta_n-\theta_0)
\ \Rightarrow\ N(0,B^{-1}V_\psi B^{-\top}).
$$

The empirical covariance of the $p$-dimensional record scores, combined with the mean derivative $\bar B_n$, consistently estimates this sandwich. Thus a rich weak family requires a $p\times p$ empirical score covariance, while its within-record weights can come from the full model covariance in Proposition 3. Fixed-noise consistency applies to the stated independent-instrument construction. It does not follow merely by inserting a same-record pilot into $D_r^\top\Omega_r^{-1}$.

**Proof.** Conditional on $u_r$ and training data, $L_r^+$ is independent of $m_r^-(\theta_0,\nu_0)$, and conversely. Hence

$$
E(\psi_r\mid u_r,\text{training})
=\tfrac12 E L_r^+\,E m_r^-+
\tfrac12 E L_r^-\,E m_r^+=0.
$$

The exact identity $\widehat\theta_n-\theta_0=\bar B_n^{-1}\bar\psi_n$ gives consistency by the law of large numbers, followed by the displayed CLT by Lindeberg and Slutsky. The empirical score covariance converges to $V_\psi$: each conditional mean is zero, so variation of the deterministic trajectories creates no between-record mean term. Replacing $\theta_0$ by the consistent estimate changes this covariance by $o_P(1)$ under the moment conditions. The quadrature assumption makes its centered contribution $o_P(1)$ after multiplication by $\sqrt n$. This is the standard estimating-equation/GMM argument [@hansen1982; @newey1994], with the independent-instrument identity addressing the incidental-initial-state setting. $\square$

### S16.2. Gaussian construction and local information reference

For known $C_{\nu_0}$, draw $\zeta_r\sim N(0,C_{\nu_0})$ independently and set $Y_r^\pm=Y_r\pm\zeta_r$. Their cross-covariance is $C_{\nu_0}-C_{\nu_0}=0$; jointly Gaussian views are independent, have identical mean and covariance $2C_{\nu_0}$. Together they recover the original record. This Gaussian randomization is an established data-fission construction [@leiner2025]. It is used here to check finite-noise instruments, not proposed as a new general splitting method.

For a supplied clean trajectory, write $D_r=E G_r$ and $\Omega_r=\operatorname{Cov}(m_r(\theta_0,\nu_0))$. The deterministic optimal instrument within the prescribed linear moments is $L_r=D_r^\top\Omega_r^{-1}$. Direct calculation gives the oracle covariance

$$
\left(\sum_r D_r^\top\Omega_r^{-1}D_r\right)^{-1}.
$$

This is a fixed-family moment reference at finite noise. Proposition 3 evaluates $\Omega_r$ including all Hermite orders. A noisy $G_r^\top$ in place of $D_r^\top$ introduces $E(G_r^\top m_r)$ and can invalidate the zero-mean score. Independent instruments preserve the mean identity while generally changing the finite-noise covariance from this oracle expression.

In the local Gaussian experiment, each record has coefficient sensitivity $S_r$ and initial sensitivity $Z_r$. The joint information for $\theta,u_1,\ldots,u_n$ has separate nuisance blocks. Its coefficient Schur complement is exactly

$$
\mathcal I_{e,n}=\sum_{r=1}^n
\left(\mathcal I_{\theta\theta}(u_r)-
\mathcal I_{\theta u}(u_r)\mathcal I_{uu}(u_r)^{-1}
\mathcal I_{u\theta}(u_r)\right)
=\sum_r\mathcal I_e(u_r).
$$

Theorem 2 applies record by record. With known shape and noise covariance $\epsilon^2 C$, a fixed test span gives the order-one moment information $D_r^\top\Omega_{r,1}^{-1}D_r$. Dense spans approach $\mathcal I_e(u_r)$ for any fixed finite collection of regular trajectories. Thus $(\sum_r\mathcal I_e(u_r))^{-1}$ is the sharp local reference for boundary-free influences in the product experiment. This exact block calculation alone supplies neither a fixed-noise incidental-parameter MLE theorem nor a growing-record, growing-span uniform efficiency result.

The symmetric independent-view score also preserves the leading influence under the following additional local conditions: its instruments converge in $L^{2+\delta}$ to the deterministic correct-covariance instruments; their perturbation is $O_P(\epsilon)$; weak moments have an $L^2$ expansion with $O(\epsilon^2)$ remainder; and covariance solves remain regular for the fixed span. Then

$$
\tfrac12\{L_r^+m_r^-+L_r^-m_r^+\}
=L_{r,0}\,m_{r,1}(\eta_r)+O_{L^2}(\epsilon^2).
$$

Indeed the first-order moments are linear in $\eta_r\pm\zeta_r$, their average cancels $\zeta_r$, and instrument perturbations multiply moments of order $\epsilon$. The factor two in each view's covariance scales $L_{r,0}$ uniformly and cancels from the normalized estimating equation. Combined with Theorem 2, this proves the conditional local product-experiment efficiency statement for fixed $n$ and sequential span enrichment. The fixed-noise CLT above and the local information reference describe complementary properties of the same weak scores; joint rates remain open.

### S16.3. Covariance estimation and common-trajectory sufficiency

For an independently calibrated covariance estimate $\widehat\nu$, asymptotic linearity of the averaged score adds $H\xi$ to its influence, where $H=E\partial_\nu\psi$ and $\xi$ is the calibration influence with the appropriate relative sample-size scaling. The Newey–McFadden nuisance argument then gives its sandwich [@newey1994]. For estimated-covariance Gaussian fission, calibration error also perturbs cross-view independence; its contribution must be included or be negligible at the inference scale. The known-covariance construction proves the stated result without presuming those extra conditions.

Uncorrected residual-profile estimates from separate fixed-horizon trajectory fits need not be consistent as records accumulate: each fit estimates a new initial state with nonvanishing error. Even in a linear Gaussian model, fitting a mean projects part of the error out of the residual and changes its expected quadratic form. A residual estimator therefore requires a demonstrated bias correction or a separate consistent calibration scheme. Pair differences eliminate the mean only for equal trajectories. These facts prevent transferring the common-trajectory pair estimator to heterogeneous initial states.

If all records instead share one initial state and known covariance, the Gaussian path average is sufficient for the mean parameter: the sample contrasts are parameter-free independent Gaussian observations. Applying a single-record procedure to the path average reduces covariance to $C/n$. For a nonlinear polynomial, a Hermite covariance term of order $q$ then scales as $n^{-q}$, whereas averaging separate corrected moments retains its fixed-amplitude contribution divided by $n$. Path averaging is the natural reference in this common-trajectory experiment. The old fixed-family GMM calculation is retained as Lemma 2 and a calibration diagnostic, rather than as an efficiency contribution.

## References
