### S8.1. General proof of Theorem 3

Use the estimator defined in main §5.6, including even uniform grids,
known Hermite correction variances, the specified guards and projection
to a fixed coefficient box. Write $C$ for the unit-amplitude covariance
of $\eta$, with known common OU correlation time as required by the covariance routine. The true trajectory is analytic on the full compact interval,
so it and all derivatives needed below are bounded. Constants in this
proof may depend on the fixed spline level and trajectory, but not on
$\epsilon$ or the observation grid. No uniform growing-level assertion
is made.

#### Auxiliary lemmas and the localization event

Here $m$ counts library functions, $r_k$ counts retained scalar tests,
and $q_k=d r_k=\dim(\widehat\Sigma_1)$ counts stacked residual moments.
The ridge is always
$\epsilon^2\operatorname{tr}(\widehat\Sigma_1)I_{q_k}/q_k$.
At a fixed level the cutoff eventually retains all raw tests, so $q_k$
is eventually deterministic. The library size $m$ is never the ridge divisor.

**Lemma S8.A (smoothed trajectory and instruments).** Let $I_h$ denote
linear interpolation and let $z=K_b*I_hY$, with the reflected extension
specified in the algorithm. Set $M=\|\eta\|_\infty$ and let
$\omega_\eta(h)=\sup_{|s-t|\leq h}\|\eta(s)-\eta(t)\|$.
For any finite $s\geq1$, $E M^s<\infty$ and
$\|\omega_\eta(h)\|_{L^s}\to0$. On the support of $w$,

$$
\|z-x_0\|_{C^1}\leq C(b^2+h^2/b+\epsilon M/b),\qquad
\|z\|_{C^3}\leq C(1+h^2/b^3+\epsilon M/b^3).
$$

On the whole observation interval,
$\|z-x_0\|_\infty\leq C(b+h^2+\epsilon M)$.
Consequently the raw feature tests and their derivatives converge in
$H^1$, in every fixed finite moment, to $w g_j(x_0)$.
The interior deterministic smoothing bias is
$b^2x_0''/22+O(b^4)$; its derivative has the corresponding expansion.
With $b=\epsilon^{1/10}/4$, its leading order is $\epsilon^{1/5}$.

**Proof.** The OU Gaussian initial value and time-changed Brownian
representation give all maximal moments. Almost-sure continuity and
$\omega_\eta(h)\leq2M$ give the modulus limit by dominated convergence.
Linear interpolation has clean error $O(h^2)$ and noise supremum at most
$M$. Convolution multiplies these errors by bounded kernel mass;
its $j$th derivative has $L^1$ norm $O(b^{-j})$.
On the support of $w$, $b<1/2$ excludes reflected boundaries.
Kernel symmetry cancels the first Taylor term; direct integration gives
$\int K_1=1$ and $\int s^2K_1(s)ds=1/11$. Taylor expansion proves
the stated bias and derivative bounds. Reflected clean trajectories are
Lipschitz, which gives the whole-interval $O(b)$ bound. Polynomial
composition, the chain rule and the maximal moments finish the proof.
$\square$

**Lemma S8.B (quadrature without differentiating OU noise).** For any
boundary-free test $\phi\in W^{3,\infty}$, including the piecewise cubic
tests, write $F_\phi=-\phi'^\top x_0-\phi^\top f_{\theta_0}(x_0)$.
On an even uniform grid, composite Simpson quadrature satisfies

$$
|\mathcal Q_hF_\phi-\textstyle\int F_\phi|
\leq C h^2\|\phi\|_{W^{3,\infty}}.
$$

For the random raw tests of Lemma S8.A, the clean discrete residual
is $o_{L^s}(\epsilon)$ for every fixed finite $s$ when $h^2=o(\epsilon)$.
For the fixed clean first-order weights
$\ell=-\phi'-J^\top\phi$, the stochastic quadrature error obeys

$$
|\mathcal Q_h(\ell^\top\eta)-\textstyle\int\ell^\top\eta|
\leq C[hM+\omega_\eta(h)]=o_{L^s}(1).
$$

**Proof.** The clean identity gives $\int F_\phi=0$.
$F_\phi'$ is absolutely continuous with bounded second derivative
almost everywhere; its bound uses the test's first three derivatives
and the bounded clean trajectory derivatives. The trapezoidal error
is $O(h^2\|F_\phi''\|_\infty)$, and Simpson's identity
$\mathcal Q_h=(4\mathcal T_h-\mathcal T_{2h})/3$ gives the claimed bound.
No fourth derivative across spline knots is needed. Lemma S8.A implies
that the random $W^{3,\infty}$ norms have uniformly bounded finite
moments: $h^2/b^3\to0$ and $\epsilon/b^3\to0$.
Thus division of the clean error by $\epsilon$ gives convergence to zero.
Simpson weights are nonnegative and have total mass $T$.
The modulus of $\ell^\top\eta$ is bounded by
$C[hM+\omega_\eta(h)]$ because the clean weight is Lipschitz.
Lemma S8.A proves its moment convergence. $\square$

**Lemma S8.C (self-instrumentation and the pilot rate).** Before
whitening, let $\Phi_\epsilon$ contain the random feature tests and the
deterministic splines, and let $\Phi_0$ be its clean limit. Then

$$
m_{\epsilon,h}(\theta_0)
=\epsilon r_{\Phi_0}(\eta)+o_{L^2}(\epsilon),\qquad
r_\phi(\eta)=-\int(\phi'^\top\eta+\phi^\top J\eta)dt.
$$

The full-library feature pilot satisfies
$\widetilde\theta-\theta_0=O_P(\epsilon)$ whenever its clean design
has full column rank. Data-dependent tests are not asserted to give
exactly unbiased finite-noise moments.

**Proof.** The clean continuous identity holds for every admissible
$\phi$, so perturbing the test produces exactly zero clean residual.
Its potentially $O(\epsilon^{1/5})$ drift therefore creates no
standalone coefficient-error term. In the polynomial Hermite expansion,
the first-order residual is the displayed Gaussian linear functional.
Replacing $\Phi_\epsilon$ by $\Phi_0$ contributes at most
$C\epsilon\|\Phi_\epsilon-\Phi_0\|_{H^1}M$, which is
$o_{L^2}(\epsilon)$ by Lemma S8.A and Holder's inequality.
Every higher-order term contains at least $\epsilon^2$ and a finite
polynomial in $M$, with uniformly bounded test moments, hence is
$o_{L^2}(\epsilon)$. Lemma S8.B removes clean and first-order
quadrature errors; bounded quadrature mass controls the higher orders.
The pilot design converges to its nonsingular clean design. Its local
least-squares inverse is continuous and the truth is interior to the
pilot box, giving the claimed $O_P(\epsilon)$ rate. $\square$

**Lemma S8.D (good-event conditions and moment control).** Fix $k$.
Let $\Gamma_0$ be the clean raw $H^1$ Gram matrix, $\Sigma_{1,0}$ the
clean whitened residual covariance, and $N_0$ the clean column-scaled
final GLS normal matrix. Assume the positive gaps established in Steps
1-2. Choose a fixed sufficiently small $c_k>0$ and put
$\mathcal E_\epsilon=\{\epsilon b^{-3}M\leq c_k\}$.
For all sufficiently small $\epsilon$ along $h^2=o(\epsilon)$,
the following conditions hold on this event:

- (E1) Proxy and test derivatives through order three are bounded;
  $\epsilon M\leq c_kb^3$, $h^2/b^3\to0$ and $b\to0$.
- (E2) The clean residual divided by $\epsilon$ is bounded;
  $h^2/\epsilon\to0$.
- (E3) The raw Gram gap is at least half its clean gap, the pilot
  design remains full rank, and the pilot lies inside its box.
- (E4) Normalized covariance and final normal-matrix gaps are at least
  half their clean gaps. Column scales stay bounded away from zero.
- (E5) The relative Gram cutoff, covariance ridge and relative
  covariance/design guards lie below those positive local gaps.
  Every raw test is retained, and all numerical-failure guards are inactive.
- (E6) The unprojected final estimate lies inside its box and
  $\epsilon^{-1}\|\widehat\theta-\theta_0\|\leq C_k(1+M)$.

The squared scaled errors are uniformly integrable.

**Proof.** Lemmas S8.A-B establish (E1)-(E2) uniformly on the event.
The proxy and test errors tend uniformly to zero there because
$\epsilon M/b\leq c_kb^2$; polynomial expansions give
pilot and final residuals bounded by $C_k(\epsilon M+h^2)$.
Positive clean gaps and continuity give (E3)-(E4), with Lemma S8.C
providing the local pilot inverse. Covariance continuity uses the
whole-interval proxy bound and the factor $\epsilon^{2q-2}$ for each
higher Hermite order $q\geq2$. The fixed positive gaps dominate the
vanishing cutoff, trace ridge and guards, proving (E5).
Local inverse bounds applied to the residual yield (E6); here
$\epsilon M+h^2\leq c_k b^3+h^2\to0$, so both projections are inactive.
All constants depend on fixed $k$ and its clean gaps.
On the complement $M>C_k'\epsilon^{-7/10}$, the box bounds every output,
including fallbacks. Markov's inequality with $E M^{12}<\infty$ gives

$$
E[\epsilon^{-4}\|\widehat\theta-\theta_0\|^4;
\mathcal E_\epsilon^c]\leq C_k\epsilon^{22/5}\to0.
$$

On $\mathcal E_\epsilon$, the fourth scaled moment is at most
$C_kE(1+M)^4$. This proves uniform integrability. $\square$

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

**Step 3: proxy and raw instruments converge.** Lemma S8.A gives
finite-moment convergence of the raw tests in $H^1$ and consistency of
the proxy on the full interval. Their Gram therefore converges to the
positive clean Gram proved in Step 2.

**Step 4: self-instrumentation and the pilot.** Lemmas S8.B-C give
$m_{\epsilon,h}(\theta_0)=\epsilon r_{\Phi_0}(\eta)+o_{L^2}(\epsilon)$,
including the numerical integration errors. Step 1 supplies the clean
pilot inverse, so $\widetilde\theta-\theta_0=O_P(\epsilon)$.
Its box projection is locally inactive. The $O(\epsilon^{1/5})$
instrument drift is absorbed in the vanishing scaled remainder because
the clean identity holds for every admissible test.

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

**Step 6: moment convergence.** The positive matrix gaps proved in
Steps 1-2 and the limits in Steps 3-5 verify Lemma S8.D's premises.
Conditions (E1)-(E6) localize the solve; its good-event fourth scaled
moment is bounded, and the complementary contribution is
$O(\epsilon^{22/5})$. Hence the squared scaled errors are uniformly
integrable. The Gaussian distributional limit in Step 5 then yields
both the second-moment matrix limit and covariance convergence.

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
