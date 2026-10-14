# Interpretation planned before the main comparison

This text is a conditional editorial plan, not a claim that the main experiment
has established a crossing, likelihood bias, or a computational speed advantage.

If the matched Monte Carlo comparison finds no crossing up to 800 records, report
that result explicitly. The deterministic curvature calculation predicts a
nonzero leading small-noise displacement of the local population profile root.
Its magnitude and direction require comparison with the full Monte Carlo bias
vector, including uncertainty of that vector. A ratio of estimated bias norms
alone can exaggerate disagreement when the Monte Carlo mean is noisy.

Boundary-free Hermite-corrected weak moments have exact conditional centering
under the specified known covariance and exact continuous integration. This
property protects estimating equations against heterogeneous initial states.
Finite-sample coefficient estimates still incur ratio bias, and fixed-grid
quadrature defects require separate control. Their centering is a structural
advantage; it does not establish an RMSE advantage over likelihood.

For the existing matched single-record 5% benchmark, local likelihood with
unknown initial state has relative RMSE 0.0172, compared with 0.0241 for author
WENDy and 0.0256 for weak GLS. Those numbers concern the earlier benchmark, not
the new heterogeneous-record experiment. Where likelihood is computationally
accessible and its finite-noise bias is small, its accuracy makes it the natural
reference. Weak equations may become preferable at larger sample sizes, higher
noise, or under computational constraints; each proposed advantage needs the
matched RMSE and timing evidence.

The timing comparison will distinguish online point estimation, interval
construction and offline numerical instrument preparation. Weak point
estimation reuses a precomputed instrument cache and avoids per-record ODE
optimization in the online phase. Building that cache requires model ODE
calculations. An unconditional description as estimation "without solving ODEs"
would conceal this cost. Bootstrap-based weak uncertainty calibration must also
be timed separately from the point estimate.

If a statistically resolved crossing is found, give the first checked sample
size and its paired Monte Carlo interval. Treat the analytic n* as a
leading-order extrapolation rather than an observed continuous crossing. If the
bootstrap-corrected likelihood is more accurate, report its accuracy and cost;
the result illustrates a trade-off between structural moment centering and
numerical bias correction.
