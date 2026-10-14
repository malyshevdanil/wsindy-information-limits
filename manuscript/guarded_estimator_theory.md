## S27. Archived regularization and a constructive guarded variant

### S27.1. What the archived solver actually approaches

In `src/weak_efficiency.py`, `gls_map` adds
$\lambda_\epsilon I=\rho\,\operatorname{tr}(\widehat\Sigma_\epsilon)I/m$
to the physical residual covariance, with fixed relative $\rho=10^{-10}$
in the primary protocols. If
$\epsilon^{-2}\widehat\Sigma_\epsilon\to\Sigma_1>0$, its normalized limit is

$$
W_\rho=\left[\Sigma_1+\rho\operatorname{tr}(\Sigma_1)I/m\right]^{-1}.
$$

Thus fixed relative ridge is not removed by small-noise scaling. Whenever a
regular local expansion exists, the corresponding leading map and covariance
are

$$
P_\rho=(D^\top W_\rho D)^{-1}D^\top W_\rho,\qquad
V_\rho=P_\rho\Sigma_1P_\rho^\top.
$$

For the ideal map $P_0$, both maps satisfy $P_\rho D=P_0D=I$ and

$$
V_\rho-V_0=(P_\rho-P_0)\Sigma_1(P_\rho-P_0)^\top\succeq0.
$$

Indeed, $P_0\Sigma_1=V_0D^\top$ annihilates the cross terms because
$(P_\rho-P_0)D=0$. Equality holds only when the two maps coincide; a small
numerical ridge alone gives no theorem-level equality. Its impact also
depends on conditioning. The archived final fits are unbounded and their
bounded covariance pilots do not control inverse random-design tails.
Neither uniform integrability nor exact ideal-GLS efficiency of these
archived algorithms is claimed. Their finite mesh and fixed spectral cutoff
also differ from the sequential limit in Proposition 3.

### S27.2. A specified algorithm for which the fixed-family assumptions follow

The separate routine `src/guarded_gls.py` specifies a bounded vanishing-ridge
variant. It retains fixed, deterministic test functions and the known sensor
covariance model. It divides physical residual covariance by $\epsilon^2$,
uses relative ridge $\rho_\epsilon=\epsilon^2$, and rejects a numerical solve
when normalized covariance or column-scaled normal-matrix relative minimum
eigenvalues fall below $10^{-8}\epsilon^2$. A rejected solve returns the
fixed coefficient-box midpoint. A successful solve is projected onto that
same compact box. These safeguards leave a regular local solution unchanged
asymptotically. This routine is a new analytically specified variant; none
of the reported Monte Carlo outcomes uses it.

**Proposition S27.1 (guarded fixed-family weak GLS).** Consider a correctly
specified polynomial vector field linear in its $p$ coefficients, known
initial state and supplied independent-coordinate OU covariance shape
with a common correlation time. Fix a finite deterministic family of
boundary-free tests. Assume:

- (G1) The true coefficient lies inside the final compact box. Its clean
  trajectory and nearby pilot-model trajectories exist on $[0,T]$ and
  depend smoothly on the pilot coefficients on a common bounded region.
- (G2) The clean corrected design $D$ and the clean design of a fixed
  finite-window pilot both have full column rank. The pilot's declared
  structural family contains the truth. Clipping the pilot to a compact
  convex box containing the truth is allowed, including truth on its boundary.
- (G3) The first-order residual covariance $\Sigma_1$ of this finite family
  is positive definite after exact dependencies are removed. The model
  covariance is evaluated from the full Hermite formula at the pilot,
  with known marginal correction variances.
- (G4) Deterministic quadrature weights have uniformly bounded total
  absolute mass for each test and converge to their continuum integrals.
  Their clean residual and pilot-design errors are $o(\epsilon)$.
  Any subtracted pilot-model quadrature defect has value $o(\epsilon)$
  and derivative $o(1)$ uniformly near the truth. Tests and their first
  derivatives are bounded and have the required quadrature regularity.

Build the pilot by corrected least squares on its stated family, using a
fixed bounded fallback if its design is singular. Construct the full
Hermite covariance at the integrated pilot, and apply the guarded routine
above to the final corrected moments. Then, for fixed family size,

$$
\epsilon^{-1}(\widehat\theta_\epsilon^{\mathrm{guard}}-\theta)
\Rightarrow N(0,V_K),\qquad
\epsilon^{-2}\operatorname{Cov}(\widehat\theta_\epsilon^{\mathrm{guard}})
\longrightarrow V_K=(D^\top\Sigma_1^{-1}D)^{-1}.
$$

The squared scaled errors are uniformly integrable. If this construction is
applied to nested dense admissible families satisfying (G1)–(G4) separately,
the subsequent family limit gives the boundary-free covariance reference
of Theorem 1. This is a sequential result; it supplies no uniform rate in
family size and no conclusion for an unchanged fixed mesh or ridge.

**Proof.** Write $Y=x+\epsilon\eta$ and $M=\|\eta\|_\infty$. For each OU
coordinate,

$$
\eta_j(t)=e^{-t/\tau_j}
\left[\eta_j(0)+\sqrt{2\sigma_{0,j}^2/\tau_j}
\int_0^t e^{s/\tau_j}dW_j(s)\right].
$$

The stochastic integral is a martingale with a Gaussian terminal value of
finite variance. Doob's $L^q$ inequality for $q>1$ bounds its supremum moment
by $(q/(q-1))^q$ times its terminal moment. Hence $E M^q<\infty$ for every
$q>1$, uniformly over evaluation grids, and
$P(M>c/\epsilon)=O(\epsilon^q)$ for every such $q$.

Polynomial corrected moments are locally Lipschitz in the path supremum norm;
their polynomial correction terms start at order $\epsilon^2$. By (G2),
pilot least squares is locally regular at the clean path. On
$\epsilon M<c$ for a sufficiently small fixed $c$, it gives
$\widetilde\theta-\theta_{\mathrm{pilot},0}=O(\epsilon M+o(\epsilon))$.
Projection onto the pilot box is nonexpansive and fixes the truth, so the
same bound survives clipping, including clipping at a true boundary value.
By (G1) and (G3), local continuity of the integrated trajectory and covariance
then gives

$$
\epsilon^{-2}\widehat\Sigma_\epsilon
=\Sigma_1+O_p(\epsilon)+o(1).
$$

Higher Hermite orders contribute $O(\epsilon^2)$ or smaller to this
normalized covariance. Consequently the model weights are consistent; this
is derived from pilot regularity and the covariance formula, rather than
assumed as a separate conclusion. By (G4), the final moment residual at
truth is $\epsilon r_1+o_p(\epsilon)$, with centered Gaussian $r_1$ of
covariance $\Sigma_1$, while the design converges to $D$.

The positive eigenvalues in (G2)–(G3) bound local inverses uniformly. Since
both guard thresholds and normalized ridge vanish, the guards are inactive
on a sufficiently small neighborhood of the clean path for all sufficiently
small $\epsilon$. In that neighborhood the final solve yields

$$
\widehat\theta_\epsilon-\theta
=\epsilon(D^\top\Sigma_1^{-1}D)^{-1}D^\top\Sigma_1^{-1}r_1
+o_p(\epsilon).
$$

The final box projection is inactive there because the truth is interior.
This proves the local distributional assertion for the specified algorithm.

For moment control, on the same good event local inverse bounds and (G4)
give $\|\epsilon^{-1}(\widehat\theta_\epsilon-\theta)\|
\leq C(1+M)$ with a constant depending on the fixed family. Its fourth
moment is bounded. Outside that event both estimate and truth are in the
fixed compact box, so the fourth scaled-error contribution is bounded by
$C\epsilon^{-4}P(M>c/\epsilon)=O(\epsilon^4)$, using $q=8$.
Thus scaled errors have bounded fourth moments, proving uniform
integrability of their squares and convergence of covariance. Finally,
S26.2's projection argument proves the separate dense-family limit. $\square$

This result strengthens the existence of a realizable efficient procedure
under verifiable clean-rank, trajectory and quadrature conditions. It neither
certifies the archived unguarded algorithm nor validates the Wald selector,
misspecified pilots or finite-noise SEs. Independent size and power studies
remain necessary for the reported testing procedure.
