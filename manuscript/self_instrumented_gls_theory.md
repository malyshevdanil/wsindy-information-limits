## S28. Constructive full-library efficiency under a jet condition

### S28.0. General statement and proof of Theorem 3

Let $x_0:[0,T]\to\mathbb R^d$ be an analytic polynomial ODE solution on
the full compact interval, $T>1$. The field is linear in $m$ polynomial
library functions in each equation; the true $dm$ coefficients are
interior to a fixed finite box. The observation is $x_0+\epsilon\eta$,
with known positive independent-coordinate stationary OU scales and
range. Suppose the $m$-by-$m$ derivative matrix of the library along the
trajectory, through order $m-1$, has nonzero determinant at some $t_*$.

The specified algorithm uses the kernel and bandwidth of S28.1 below,
and $w(t)=[1-((t-T/2)/(T/2-1/2))^2]_+^4$. Its $m$ feature instruments
are $w g_j(z)$, evaluated on the smoothed linear interpolant $z$, giving
$dm$ full-library pilot equations. Add the interior dyadic cubic spline
family at each fixed level. Whiten its full $H^1$ Gram matrix with
relative retention cutoff $10^{-10}\epsilon^{1/10}$. Use the exact
polynomial Hermite covariance at the proxy and box-clipped pilot,
relative covariance ridge $\epsilon^2$, normal-matrix and covariance
guards $10^{-8}\epsilon^2$, and final projection to the coefficient box.
Failed guards return its midpoint. Take $h^2=o(\epsilon)$.

We prove the conclusions of Theorem 3; no pilot-rank or weight-consistency
assumption is added. Analyticity and the jet condition imply that a
linear combination of library values vanishing on any open interval
has all coefficients zero. Therefore $\int w g g^\top$ is positive
definite and the clean pilot design has $d$ invertible copies of it.
The combined clean test Gram is positive definite too: a relation
between $w c^\top g$ and cubic splines makes the analytic left side a
single cubic polynomial on the support interior. Its order-four zero
at the left edge forces that cubic to vanish. The jet condition and
spline independence then eliminate both sets of coefficients. The
adjoint on zero-trace tests is injective, and positive OU covariance
makes the first-order finite-family residual covariance positive
definite. The final design includes the invertible pilot block.

With $M=\|\eta\|_\infty$, kernel moment cancellation gives proxy $C^1$
error at most $C(b^2+h^2/b+\epsilon M/b)$ on the instrument support,
and a $C^3$ bound $C(1+h^2/b^3+\epsilon M/b^3)$, where
$b=\epsilon^{1/10}/4$. The proxy is uniformly consistent on the full
interval; reflection adds only vanishing boundary bias. The clean
weak identity holds for every such instrument, so its derivative with
respect to the test is zero. Hence instrument estimation creates no
standalone first-order clean residual. Taylor expansion of the finite
polynomial library, exact Hermite correction, and OU path continuity
give the same linear Gaussian residual and root-$\epsilon$ pilot
expansion as in S28.3. The clean quadrature error is $o(\epsilon)$;
noise Riemann sums converge with domination by $CM$. For any fixed
polynomial degree, higher-order residual terms are bounded by $C\epsilon M$
on a neighbourhood where $\epsilon M$ is small, and their first-order
remainders vanish in $L^2$ by the finite Gaussian maximal moments.

Proxy and pilot consistency imply convergence of the normalized exact
Hermite covariance to its positive first-order limit. All clean Gram
gaps are positive at fixed level. Thus the vanishing cutoff eventually
retains every row, the guards become inactive, and the vanishing ridge
gives the GLS linear map in S28.3. Clipping is locally inactive because
the true coefficients are interior. For the good event
$\epsilon b^{-3}M\leq c$, local inverse and derivative bounds give
$\epsilon^{-1}\|\widehat\theta-\theta_0\|\leq C(1+M)$.
On its complement the box bounds the error, and the twelfth Gaussian
maximal moment gives $C\epsilon^{-4}\epsilon^{(7/10)12}\to0$ for the
scaled fourth moment there. This proves uniform integrability and
second-moment convergence, not just the distributional limit.

Finally the interior cubic spaces are nested and dense in $H_0^1$.
Continuity of the adjoint and Proposition 2 imply density in the weak
influence class. The full efficient information is positive definite:
if a parameter sensitivity equals $Zc$, its zero initial value forces
$c=0$; its differential equation and library independence force every
parameter coefficient to be zero. Theorem 1 then gives the stated
sequential covariance limit. The specific calculations below verify
the general jet hypothesis and spell out all bounds for Van der Pol.

### S28.1. Verified Van der Pol instance

This section verifies the small-noise assumptions for a specified guarded
procedure, rather than placing them among the hypotheses. It concerns
`src/self_instrumented_gls.py`, not the archived fixed-relative-ridge solver.
The mathematical statement treats its formulas in exact arithmetic; the
independent numerical study checks their floating-point implementation.

**Proposition S28.1.** Let the true trajectory solve
$x'=v$, $v'=-x+v-x^2v$, $(x(0),v(0))=(2,0)$ on $[0,20]$.
Observe $Y=x_0+\epsilon\eta$, where the two coordinate OU processes are
independent, have known positive scales and known range $\tau>0$.
Fit all eighteen coefficients in the two-equation library
$g=(x,v,x^2,xv,v^2,x^3,x^2v,xv^2,v^3)$, in the box $[-3,3]^{18}$.
For each fixed dyadic spline level $k\geq3$, use the procedure below with
$h^2=o(\epsilon)$ and $h\to0$. Then its pilot is root-$\epsilon$
consistent, its normalized covariance weights are consistent, and

$$
\epsilon^{-1}(\widehat\theta_{k,\epsilon}-\theta_0)
\ \Rightarrow\ N(0,V_k),\qquad
\epsilon^{-2}E[(\widehat\theta_{k,\epsilon}-\theta_0)
(\widehat\theta_{k,\epsilon}-\theta_0)^\top]\to V_k.
$$

Here $V_k$ is the optimal covariance in the stated finite span of
boundary-free tests. In the sequential limit,
$V_k\to\mathcal I_{\theta\mid u}^{-1}$ as $k\to\infty$.
This is a small-noise statement at this identifiable trajectory. It gives
neither an exact finite-noise Wald size nor a simultaneous growing-$k$ limit.

### S28.1. Specified algorithm

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
$\epsilon^2\operatorname{tr}(\widehat\Sigma_1)/m$, and perform GLS.
The pilot and final normal matrices use relative minimum-eigenvalue guards
$10^{-8}\epsilon^2$; the normalized covariance has the corresponding guard.
Any failure returns the box midpoint. Clip the final estimate to the box.
These are the choices in the implementation, including its bounded fallback.

### S28.2. Exact identifiability and existence

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

### S28.3. Pilot and covariance expansions

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
The clean trapezoidal defect is $O(h^2)$ because the test derivatives
are piecewise smooth with bounded relevant derivatives. The noise
Riemann sums converge to their continuous integrals by OU path
continuity, with domination by $C M$. Hence their scaled errors vanish
in $L^2$. Rank stability now gives the pilot expansion and root-$\epsilon$
consistency. Clipping is inactive with probability tending to one, because
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

### S28.4. Moment control and dense-span limit

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

## S29. Interpretation of the independent validation rule

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
