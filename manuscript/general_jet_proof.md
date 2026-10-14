### S8.1. General proof of Theorem 3

Use the estimator defined in main §5.6, including even uniform grids,
known Hermite correction variances, the specified guards and projection
to a fixed coefficient box. Write $C$ for the unit-amplitude covariance
of $\eta$. The true trajectory is analytic on the full compact interval,
so it and all derivatives needed below are bounded. Constants in this
proof may depend on the fixed spline level and trajectory, but not on
$\epsilon$ or the observation grid. No uniform growing-level assertion
is made.

**Step 1: the jet condition gives clean rank.** Put
$G(t)=(g_1(x_0(t)),\ldots,g_m(x_0(t)))^\top$. If $c^\top G(t)$ vanishes
on any open subinterval, analyticity extends that identity to the whole
interval. Differentiation at $t_*$ through order $m-1$ and the nonzero
jet determinant give $c=0$. Since $w>0$ on $(1/2,T-1/2)$,

$$
c^\top\left(\int wGG^\top dt\right)c
=\int w(c^\top G)^2dt>0\quad(c\ne0).
$$

The clean feature-instrument pilot design contains $d$ copies of this
positive matrix. Hence it is invertible without a separate rank
assumption. The final design contains this block and has full column
rank as well.

**Step 2: the combined Gram and residual covariance are positive.**
Suppose a linear combination $w c^\top G$ equals a cubic-spline
combination. On the support interior the former is analytic. Across
every interior spline knot, analyticity forces adjacent cubic pieces
to have equal derivatives of all orders, and hence to be one cubic
polynomial. At $t=1/2$ the analytic trajectory is bounded and $w$ has
an order-four zero. The polynomial and its first three derivatives
therefore vanish at that edge, forcing it to be zero. Step 1 gives
$c=0$, and independence of the cardinal splines gives zero for their
coefficients. Thus the full clean $H^1$ Gram has a positive minimum
eigenvalue at every fixed level.

For vector tests the adjoint $D^*\phi=-\phi'-J^\top\phi$ is injective
on $H_0^1$: a homogeneous solution with zero initial trace is zero
by uniqueness for the linear equation. Each OU coordinate covariance
is strictly positive on nonzero $L^2$ weights, as its spectral density
$2\sigma_{0,j}^2\tau/(1+\tau^2\omega^2)$ is positive. Consequently the
first-order covariance $\Sigma_1$ of any independent finite test family
is positive definite. Both covariance and design inverse bounds are
available in a neighbourhood of the clean data.

**Step 3: proxy and raw instruments converge.** Couple the observations
on successive grids to one continuous unit-amplitude OU path and write
$M=\|\eta\|_\infty$. The OU representation as a Gaussian initial value
plus a time-changed Brownian martingale and the maximal inequality give
$E M^q<\infty$ for every finite $q$. Linear interpolation has clean
uniform error $O(h^2)$ and noise bounded by $M$. Kernel symmetry and
its unit mass give, on the support of $w$,

$$
\|z-x_0\|_{C^1}
\leq C(b^2+h^2/b+\epsilon M/b),\qquad
\|z\|_{C^3}\leq C(1+h^2/b^3+\epsilon M/b^3).
$$

The kernel never meets reflected boundaries on this support.
On the entire interval reflection introduces only $O(b)$ clean bias,
so the proxy remains uniformly consistent there. With
$b=\epsilon^{1/10}/4$ and $h^2=o(\epsilon)$, the first bound tends to
zero in every fixed finite moment. Polynomial composition and its
chain rule imply convergence of the feature instruments to $wG$ in
$H^1$. The cubic tests are deterministic. The combined raw Gram
therefore converges to its positive clean Gram.

**Step 4: adaptive tests contribute no clean first-order residual.**
For every boundary-free test, including the random raw instruments,
the continuous clean identity is exactly
$b_\phi(x_0)-G_\phi(x_0)\theta_0=0$. Its derivative with respect to
$\phi$ is zero. Polynomial Hermite expansion of the data residual gives

$$
m_\epsilon(\theta_0)
=\epsilon r_{\Phi_0}(\eta)+o_{L^2}(\epsilon),\qquad
r_\phi(\eta)=-\int(\phi'^\top\eta+\phi^\top J\eta)dt.
$$

Before whitening, replacing $\Phi_\epsilon$ by $\Phi_0$ in the linear
term has scaled error bounded by
$C\|\Phi_\epsilon-\Phi_0\|_{H^1}M$, whose second moment tends to zero
by Step 3 and the Gaussian maximal moments. Higher Hermite orders are
finite polynomials multiplied by at least $\epsilon^2$; their scaled
remainders vanish in $L^2$. For quadrature, the clean residual defect
is $O(h^2)$ on bounded-$C^3$ instrument neighbourhoods. Simpson weights
have bounded total mass. OU path continuity makes noise quadrature
errors converge to zero after scaling by $\epsilon$, with domination
by a fixed polynomial in $M$. One may first restrict to such
neighbourhoods and then use the Gaussian moment bounds to remove the
restriction. Thus no rate of differentiability of the OU path is
assumed. The clean pilot inverse from Step 1 now gives
$\widetilde\theta-\theta_0=O_P(\epsilon)$. The pilot projection is
locally inactive because the truth is strictly interior to its box.

**Step 5: weights, cutoff and guards converge.** The normalized exact
Hermite covariance at the proxy and pilot converges to $\Sigma_1$.
This follows from their consistency and continuity of the polynomial
coefficient weights and OU integrals; the order-$q\geq2$ terms carry
$\epsilon^{2q-2}$ after normalization and vanish. By Step 2 the
spectral cutoff eventually retains every row at fixed level.
For analysis, whitening can use the continuous symmetric inverse
square root of the Gram. The implemented eigenvector whitening differs
by an orthogonal row transformation, to which GLS, trace ridge and
the guards are invariant. Repeated eigenvalues therefore create no
extra eigenvector-consistency hypothesis.

The relative ridge and guard thresholds tend to zero, while the clean
matrix gaps stay positive. The final map is consequently

$$
\epsilon^{-1}(\widehat\theta-\theta_0)
=P_0 r_{\Phi_0}(\eta)+o_P(1),\qquad
P_0=(D^\top\Sigma_1^{-1}D)^{-1}D^\top\Sigma_1^{-1}.
$$

Here $D$ and $\Phi_0$ are expressed in the same clean whitened coordinates.
The Gaussian leading term has covariance
$V_k=(D^\top\Sigma_1^{-1}D)^{-1}$. This derives the expansion and weight
consistency, rather than assuming them.

**Step 6: moment convergence follows from the bounded estimator.**
For a fixed sufficiently small $c$, let the good event be
$\epsilon b^{-3}M\leq c$. On it, proxy derivatives, clean-neighbourhood
inverses and test transformations are uniformly bounded for sufficiently
small $\epsilon$. Since $\epsilon M\leq c b^3$, the finite polynomial
remainders and the quadrature bounds yield
$\epsilon^{-1}\|\widehat\theta-\theta_0\|\leq C(1+M)$.
The complement implies $M>C\epsilon^{-7/10}$. The box bounds the
unscaled output error on that event, including all guard failures.
Using the twelfth maximal moment,

$$
E[\epsilon^{-4}\|\widehat\theta-\theta_0\|^4;\mathrm{bad}]
\leq C\epsilon^{-4}\epsilon^{(7/10)12}
=C\epsilon^{22/5}\to0.
$$

The good-event fourth moment is bounded by $C E(1+M)^4$. The squared
scaled errors are therefore uniformly integrable. Together with the
Gaussian distributional limit, this proves the asserted second-moment
and covariance limits.

**Step 7: dense spans give the weak-class bound.** Interior cardinal
cubic spaces are nested under their refinement relation. Approximate
an $H_0^1$ function by a smooth compactly supported function, then
apply cubic quasi-interpolation on a sufficiently fine mesh; its
support remains interior. This proves density in $H_0^1$. The adjoint
is continuous into the covariance-dual norm, so Proposition 2 gives
density of the corresponding weak influences. Adding the fixed clean
feature instruments preserves density. The efficient information is
positive: if $S a=Z c$, evaluating at zero gives $c=0$, and the
sensitivity equation and Step 1 give $a=0$. The finite-span projection
argument of Theorem 1 therefore gives
$V_k\to\mathcal I_{e,C}^{-1}$ after the fixed-level small-noise limit.
All constants may depend on $k$, so this establishes a sequential
limit, not a joint growth rate. The concrete verification below checks
the jet hypothesis for Van der Pol; it is not used as a premise in
this general proof.
