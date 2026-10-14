# Record-wise cross-fitting

Chernozhukov, Chetverikov, Demirer, Duflo, Hansen, Newey and Robins (2018), *Double/debiased machine learning for treatment and structural parameters*, The Econometrics Journal 21(1), C1–C68, DOI https://doi.org/10.1111/ectj.12097. Crossref metadata verified. Author full text inspected at https://arxiv.org/html/1608.00060v7, Section 3.1 and the sample-splitting construction. Cross-fitting and orthogonal-score sample reuse are established statistical tools, not a novelty claim of this paper.

For this weak-score construction, the population moment is exactly zero for every coefficient instrument pilot at known covariance. This permits a direct conditional second-moment argument fold by fold; it does not justify conditioning on all complement pilots together. Assumption (A6) makes the required convergence, continuity and envelope conditions explicit. Covariance nuisance parameters entering Hermite corrections are excluded from this exact identity.

The finite-sample bootstrap resamples within the originally assigned folds, retaining all copies of a physical record in one fold and refitting each complement pilot. This accounts for pilot variability. Nominal coverage remains an empirical question; the cross-fitting theorem alone does not resolve the previous 92% result.

Computational validation preceded Monte Carlo inference. An affine cache failed the prespecified 5% trace-variance-loss threshold, and normalized affine weights also failed. Their notebooks and diagnostics are preserved. No inferential results from those approximations are used. A richer bounded approximation is being assessed against exact model instruments before promotion.
