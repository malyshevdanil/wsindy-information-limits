# Added-term detection: local nuisance-adjusted reference

The null is the original Van der Pol trajectory. The added coefficient delta
multiplies x^3 in the second equation. It is not the removal of the stabilizing
x^2 v term. The covariance is fixed across alternatives and calibrated on the
original trajectory; recalibrating noise from each alternative would change
the testing experiment.

Let s_delta denote the added-term trajectory sensitivity. Let A collect other
coefficient sensitivities. The local known-initial likelihood information is

    I_delta|A = ||s_delta - Proj_H,A s_delta||_H^2.

For boundary-free weak influences, add the initial sensitivities Z to the
nuisance projection:

    I_delta|A,Z = ||s_delta - Proj_H,[A,Z] s_delta||_H^2.

This is also the local unknown-initial likelihood reference. The latter is
equivalently the Schur complement for delta after projecting other coefficients
out of I_e. The raw diagonal I_e[delta,delta] would treat those coefficients as
known and give the wrong composite-model reference.

For full-rank finite matrices, the variance is (I^{-1})[delta,delta], using
the appropriate complete information matrix. The numerical implementation
residualizes whitened score columns by SVD instead of subtracting close
Gram matrices. It checks the cutoff and grid refinement before Monte Carlo.

In the local Gaussian shift experiment the equal-prior error between the
least-favorable nearby hypotheses is approximately
Phi(-|delta| sqrt(I_delta|nuisance)/2). The frozen benchmark uses error 0.05,
so c=2 Phi^{-1}(0.95) and delta_min_local=c/sqrt(I_delta|nuisance).
This is not the same criterion as a size-0.05 test with specified power, and
it is not an exact finite-amplitude guarantee for a sparse selection algorithm.

The primary nuisance dictionary matches the complete author WSINDy polynomial
library of degrees 1 through 3 in both equations: eighteen coefficients,
seventeen other than the candidate, plus the two initial coordinates for the
boundary-free comparison. A separate five-coefficient reference profiles only
the original four coefficients. Neither reference assumes known values for
the other coefficients in its stated dictionary.

The full-dictionary projection is a regular local composite-model reference.
The actual perturbation experiment keeps the sparse background fixed, and
WSINDy selects a smaller support. Its inclusion curve can therefore beat that
conservative reference without contradicting an information bound: a uniform
composite testing bound would also need least-favorable nuisance perturbations.
The core5 reference helps separate this dictionary effect. Neither curve is
to be labeled a universal support-recovery impossibility frontier.

The author's adaptive support-selection algorithm is run without threshold
tuning on the new study. Its null inclusion rate and Wilson intervals must be
reported alongside inclusion under alternatives. A visually sharp inclusion
curve with large null inclusion does not establish calibrated detection.
The numerical alternative grid uses negative restoring cubic coefficients;
claims about large positive perturbations require a separate boundedness study.

Completed: 1000 author fits; frozen source and configuration verified.
At 0.5% noise full-library local references are 0.01927944 (known u) and 0.01940965 (weak).
Inclusion counts at 0, 0.5, 1, 2, 4 weak-reference units: 0, 0, 0, 0, 100;
at 5% noise: 2, 4, 24, 58, 70, each out of 100. Full results and Wilson intervals
are in reports/support_extra_term_results.json and supplement S20.
