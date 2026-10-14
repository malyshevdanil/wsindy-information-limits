## S8. Constructive full-library efficiency under a jet condition

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

### S8.2. Verified Van der Pol instance

This section verifies the small-noise assumptions for a specified guarded
procedure, rather than placing them among the hypotheses. It concerns
`src/self_instrumented_gls.py`, not the archived fixed-relative-ridge solver.
The independent numerical study checks the specified implementation.

**Proposition S8.2.** Let the true trajectory solve
$x'=v$, $v'=-x+v-x^2v$, $(x(0),v(0))=(2,0)$ on $[0,20]$.
Observe $Y=x_0+\epsilon\eta$, where the two coordinate OU processes are
independent, have known positive scales and known range $\tau>0$.
Fit all eighteen coefficients in the two-equation library
$g=(x,v,x^2,xv,v^2,x^3,x^2v,xv^2,v^3)$, in the box $[-3,3]^{18}$.
For each fixed dyadic spline level $k\geq3$, use the procedure below with
$h^2=o(\epsilon)$ and $h\to0$. Then its pilot satisfies $\widetilde\theta-\theta_0=O_P(\epsilon)$, its normalized covariance weights are consistent, and

$$
\epsilon^{-1}(\widehat\theta_{k,\epsilon}-\theta_0)
\ \Rightarrow\ N(0,V_k),\qquad
\epsilon^{-2}E[(\widehat\theta_{k,\epsilon}-\theta_0)
(\widehat\theta_{k,\epsilon}-\theta_0)^\top]\to V_k.
$$

Here $V_k$ is the optimal covariance in the stated finite span of
boundary-free tests. In the sequential limit,
$V_k\to\mathcal I_{e,C}^{-1}$ as $k\to\infty$.
This is a small-noise statement at this identifiable trajectory. It gives
neither an exact finite-noise Wald size nor a simultaneous growing-$k$ limit.

#### Algorithm in the verified instance

Interpolate the observations linearly and convolve with
$K_b(s)=b^{-1}(315/256)(1-(s/b)^2)_+^4$, $b=\epsilon^{1/10}/4$.
Reflect the interpolant at the endpoints. Denote the proxy by $z$.
Let $w(t)=[1-((t-10)/9.5)^2]_+^4$ and form nine tests
$w g_j(z)$ and their analytic derivatives. Fit the eighteen coefficients
by ordinary least squares in these eighteen weak equations, using the
known-variance Hermite corrections. Clip the pilot to the coefficient box.
No generating support or ODE solution is supplied to this pilot.

Append the cardinal cubic splines
$B_3(t/a-j)$, $a=20/2^k$, $j=2,\ldots,2^k-2$.
Whiten their combined $H^1$ Gram matrix, retaining eigenvalues above
$10^{-10}\epsilon^{1/10}$ times its largest eigenvalue.
Compute the exact polynomial Hermite covariance at the proxy and clipped
pilot. Normalize this covariance by $\epsilon^2$, add relative ridge
$\epsilon^2\operatorname{tr}(\widehat\Sigma_1)I_{q_k}/q_k$, where $q_k=d r_k$ counts stacked residual moments, and perform GLS.
The pilot and final normal matrices use relative minimum-eigenvalue guards
$10^{-8}\epsilon^2$; the normalized covariance has the corresponding guard.
Any failure returns the box midpoint. Clip the final estimate to the box.
These are the choices in the implementation, including its bounded fallback.

### S8.3. Exact identifiability and existence

Put $L=v\partial_x+(-x+v-x^2v)\partial_v$. Exact integer arithmetic gives

$$
\det\big[(L^r g_j)(2,0)\big]_{r=0,\ldots,8;\ j=1,\ldots,9}
=-65265167726360395776\ne0.
$$

The script `scripts/verify_vdp_jets.py` constructs every entry and verifies
the determinant without integrating an ODE. The true solution is global:
for $E=(x^2+v^2)/2$, $E'=v^2(1-x^2)\leq2E$, and Gronwall prevents finite-time
escape. Polynomial ODE solutions are analytic on their interval of
existence. If a linear combination of $g_j(x_0(t))$ vanishes on an open
interval, analytic continuation and the displayed determinant imply that
all its coefficients vanish. Thus
$\int w(t)g(x_0(t))g(x_0(t))^\top dt$ is positive definite. The clean
full-library pilot design consists of two copies of this matrix and is
invertible. This proves rank instead of assuming it.

The combined clean test Gram is also positive definite. To see this,
suppose $w\sum c_jg_j(x_0)$ equals a linear combination of the cubic
splines. On $(0.5,19.5)$ the left side is analytic. A piecewise cubic
function that is analytic across its knots is one cubic polynomial there.
At $0.5$ the left side has a zero of order at least four, so this polynomial
is zero. Analytic independence gives $c=0$; independence of the cardinal
splines then gives zero for their coefficients. This argument works at
every fixed level, including knots at the support edges.

For vector tests, the adjoint $D^*\phi=-\phi'-J^\top\phi$ is injective
on $H_0^1$: $D^*\phi=0$ with $\phi(0)=0$ has the unique zero solution.
The OU covariance is strictly positive on nonzero $L^2$ functions.
Consequently the finite-family first-order residual covariance $\Sigma_1$
is positive definite. The final design has full column rank because it
contains the invertible pilot block. The full efficient information is
positive definite as well: a parameter sensitivity in the initial-state
span must have zero initial value, hence zero initial-state coefficient;
its differential equation and library independence then force its parameter
coefficient to vanish.

### S8.4. Pilot and covariance expansions

Write $M=\|\eta\|_\infty$. All moments of $M$ are finite, by the OU
representation as a Gaussian initial value plus a time-changed Brownian
martingale and the maximal inequality. On the support of $w$ the kernel
does not encounter a reflected boundary. Linear interpolation and kernel
moment cancellation give

$$
\|z-x_0\|_{C^1(\operatorname{supp}w)}
\leq C(b^2+h^2/b+\epsilon b^{-1}M).
$$

The proxy is uniformly consistent on the full interval too, with a
boundary bias of order $b$. Its first three derivatives on the support
are bounded by $C(1+h^2/b^3+\epsilon b^{-3}M)$.
Since $h^2=o(\epsilon)$ and $b=\epsilon^{1/10}/4$, all these statements
hold with errors tending to zero on events of probability tending to one.
Thus the random instruments converge in $H^1$ to $w g_j(x_0)$.

Crucially, for every zero-trace test $\phi$, including a data-dependent
one, the clean continuous weak residual is identically zero. Its variation
with respect to the test is therefore zero. Expanding the Hermite-corrected
residual gives $\epsilon r_{\Phi_0}(\eta)+o_{L^2}(\epsilon)$.
The random-instrument perturbation contributes only
$\epsilon\,O(\|\Phi-\Phi_0\|_{H^1}M)$, and the Hermite terms of order
at least two contribute $O(\epsilon^2(1+M^3))$ locally.
The clean composite-quadrature defect is $O(h^2)$ because the test derivatives
are piecewise smooth with bounded relevant derivatives. The noise
Riemann sums converge to their continuous integrals by OU path
continuity, with domination by $C M$. Hence their scaled errors vanish
in $L^2$. Rank stability gives the pilot expansion and $\widetilde\theta-\theta_0=O_P(\epsilon)$. Clipping is inactive with probability tending to one, because
the true coefficients are strictly inside the box.

Continuity of the exact Hermite covariance, proxy consistency and pilot
consistency yield $\epsilon^{-2}\widehat\Omega\to\Sigma_1$.
The higher Hermite orders vanish after this normalization. For a fixed
level the clean Gram gap is positive, so the vanishing spectral cutoff
eventually retains every test. The whitening is an invertible row change;
GLS is invariant to it before ridge is added. Its limiting ridge is zero.
All the guards are consequently inactive with probability tending to one.
It follows directly that

$$
\epsilon^{-1}(\widehat\theta-\theta_0)
= (D^\top\Sigma_1^{-1}D)^{-1}D^\top\Sigma_1^{-1}
r_{\Phi_0}(\eta)+o_P(1).
$$

The sign of $r$ follows the convention $b-G\theta$. This linear Gaussian
map has covariance $V_k=(D^\top\Sigma_1^{-1}D)^{-1}$.

### S8.5. Moment control and dense-span limit

For fixed $k$, choose a fixed small $c>0$ and the good event
$\epsilon b^{-3}M\leq c$. For sufficiently small $\epsilon$, all inverse
matrices there are uniformly bounded, all required test derivatives are
bounded, and the residual polynomial bounds give
$\epsilon^{-1}\|\widehat\theta-\theta_0\|\leq C(1+M)$.
On its complement, clipping bounds the unscaled error by the box diameter.
Since the complement implies $M>C\epsilon^{-7/10}$, Markov's inequality
with the twelfth moment gives

$$
E[\epsilon^{-4}\|\widehat\theta-\theta_0\|^4;
\text{bad event}]\leq C\epsilon^{-4}\epsilon^{(7/10)12}
=C\epsilon^{22/5}\to0.
$$

The good-event fourth moment is bounded by $C E(1+M)^4$.
Thus the squared scaled errors are uniformly integrable; covariance
convergence follows from the distributional limit. This supplies the
missing moment argument for the bounded estimator.

The retained interior cardinal cubic spline spaces are nested under the
cubic refinement relation. Their union is dense in $H_0^1$: first
approximate by smooth compactly supported functions, then apply cubic
spline quasi-interpolation on a sufficiently fine grid. The adjoint map
is continuous from $H_0^1$ to the OU covariance norm, so their weak
influences are dense in the boundary-free class identified in Proposition 2.
Adding the nine feature instruments preserves density. The finite-span
optimal covariances therefore converge to the sharp boundary-free bound.
The order of limits is essential: first $\epsilon\to0$ for each fixed
level, then the level tends to infinity. No statement about the old
fixed-ridge estimates follows from this proposition.

## S9. Interpretation of the independent validation rule

The new calibration study uses fresh observations, unchanged estimators,
fixed sample counts and frozen seed ranges. Its primary size rule requires
a two-sided 95% Clopper--Pearson interval wholly within $[0.03,0.07]$.
The small-noise power rules require their intervals wholly within five
percentage points of the specified local Gaussian predictions; those
rules are interpreted only after size validation. These are tolerance-band
checks at specified models, not proofs of exact size 5% or uniform validity.

For a conjunctive declaration that all such cells pass, the probability
of falsely declaring that all their true probabilities belong to their
bands is at most 5%. Indeed, if any probability lies outside its band,
the declaration implies that this cell's exact 95% interval excludes its
true probability. That event has probability at most 5%. This is an
intersection--union argument; dependence between paired power cells
does not invalidate it. The individual intervals remain pointwise, and
are not claimed to form a simultaneous 95% confidence set.
