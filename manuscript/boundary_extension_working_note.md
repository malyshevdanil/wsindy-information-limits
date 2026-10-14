# Known-initial endpoint extension: working derivation

Status: internal mathematical note. The 32-path direct screen is complete;
the revised paired confirmation uses the remaining 1968 paths, including a
model-independent test family. Both paired criteria passed: universal tests recover 97.1% of the RMSE gap (95% interval 95.5–98.6%). See reports/endpoint_confirmation_results_RU.md.
The endpoint adaptation itself is due to Brunel et al. (2014, section 2.3).
The paired confirmation quantifies recovery by a finite universal family;
the endpoint adaptation is credited rather than presented as a new principle.

Let J(t)=D_x f_theta(x(t)), let u=x(0) be known, and let H be the OU
Cameron–Martin space (equivalent to H^1 for fixed positive sigma and tau).
For a vector test phi in H^1, the endpoint-enabled residual is

    m_phi = -integral(phi' Y + phi f_theta^H(Y))
            + phi(T)^T Y(T) - phi(0)^T u.

Its linear noise functional is

    L phi = (-phi' - J^T phi) dt + phi(T) delta_T.

There is no negative phi(0) delta_0 atom: u is exact and nonrandom.
For g in H, integration by parts gives

    <g,L phi> = integral phi^T(g' - Jg) + phi(0)^T g(0).

Consequently an annihilator satisfies g'=Jg and g(0)=0, hence g=0.
The range is dense in H*, equipped with the noise-variance norm. This establishes
the linearized full-class closure; an efficient feasible estimator additionally
requires controlled pilot, quadrature and small-noise remainders. Density by
itself proves none of those implementation conditions.

The OU inverse covariance applied to a smooth sensitivity S has interior
density (S-tau^2 S'')/(2 sigma^2 tau), initial atom
(S(0)-tau S'(0))/(2 sigma^2), and terminal atom
(S(T)+tau S'(T))/(2 sigma^2). With unequal coordinate noise scales this is
applied separately to each coordinate. Solve the adjoint ODE backwards with
terminal trace equal to the terminal atom. Its residual influence has the
correct interior density and terminal atom; add the initial atom times Y(0)-u
to obtain the full score influence. Four vector tests plus two initial controls
give a concrete finite moment span containing these four score influences at
the pilot trajectory. GLS uses the complete polynomial moment covariance.
This is score-equivalent in first order. Agreement with likelihood confirms
that realization; it is not independent evidence for universally constructed
weak tests. At finite noise the Hermite weak equations are not literally the
Gaussian likelihood score. The endpoint controls also make this a six-moment,
four-parameter system, rather than a four-moment exactly identified fit.

An explicit initial control is within the endpoint-enabled closure. For a smooth
cutoff chi equal to one at zero and zero beyond one, take phi_epsilon(t)=a chi(t/
epsilon). Then -phi_epsilon' dt converges in noise-variance norm to a delta_0;
the J term vanishes, by covariance continuity. This explains the numerical
implementation without silently identifying ordinary H^1 tests with atoms.
Even phi(T)=0 tests have dense range by the same annihilator argument; Brunel's
initial-trace adaptation already supplies the essential geometry.

The screen reports an endpoint-plus-initial-control implementation, not an
ordinary smooth-test minimum. It evaluates exact covariance for the discretized
moments on directly integrated pilot trajectories; it uses no model cache or
polynomial interpolation surrogate. Sampling the directly computed adjoint and
numerically integrating its moments still introduces quadrature error, checked
before Monte Carlo. Cubic interpolation inside the adjoint solver is numerical
ODE evaluation, not an information-exact continuous computation.

When u is unknown these residuals depend explicitly on u. The nuisance-immune
boundary-free theorem and its Schur-complement benchmark remain separate.

A finite universal family can approximate this enlarged class without solving
the model to build its tests. The confirming implementation uses the archived
480 interior bump windows, twelve clipped half-bumps and the initial controls.
Covariance weights still use a directly integrated pilot model, as before.
Consequently passing this test supports recovery by model-independent windows;
it does not establish a completely ODE-free estimator or a computation-time
advantage. Failing it only limits the tested family and weight implementation.
The dense-class argument prevents interpreting such a failure as proof that
only model-built test functions could ever recover the information.
