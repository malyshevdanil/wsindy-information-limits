"""Prepare stable manuscript templates; never change execution snapshots."""
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'


def prepare():
    old=(BASE/'archive/pre_efficiency/article_revised.md').read_text(encoding='utf-8')
    theory=(BASE/'weak_efficiency_theory.md').read_text(encoding='utf-8')
    front=old[:old.index('### 4.3.')]
    front=front.replace('Fixed-horizon information limits and polynomial weak-residual covariance under correlated measurement error',
        'Efficient weak-form estimation under correlated measurement error: fixed-horizon information limits and the cost of boundary-free test functions')
    a=front.index('## Abstract');b=front.index('**Keywords:**')
    front=front[:a]+'''## Abstract

Weak-form dynamical identification avoids differentiating noisy measurements, but removing boundary terms also discards information. For deterministic ODEs observed in stationary Ornstein-Uhlenbeck measurement error, we characterize the closure of boundary-free adjoint influences: it is exactly the covariance-dual subspace orthogonal to initial-state sensitivities. This gives a sharp local covariance bound equal to the inverse initial-state Schur complement, even when the initial state is known. We construct its optimal influence and an explicit generalized test function with endpoint jumps. Compact tests attain the bound as an infimum. Correct-covariance weak generalized least squares approaches this bound for dense test spans, in a sequential small-noise limit with controlled quadrature and regular pilots. A finite Hermite expansion supplies exact coupled polynomial residual weights. A fixed-horizon two-point bound separates this efficiency question from consistency under denser sampling. Experiments test a smoothed scalar optimal function with a fitted pilot and nested 30-, 120- and 480-window Van der Pol families, comparing physical OU and iid weights. Small-noise accuracy is compared with the information reference; finite-noise bias and dispersion are reported separately. The contribution identifies which information weak equations can retain and how covariance weighting and test-space approximation govern its recovery.

'''+front[b:]
    a=front.index('We develop an explicit');b=front.index('## 2.',a)
    front=front[:a]+'''Our main result identifies the complete boundary-free adjoint influence class and constructs the influence attaining its sharp information reference in covariance-norm closure. Correct temporal covariance and dense test spans then yield a conditional, sequential small-noise efficiency result for WENDy-like estimation. An explicit OU adjoint construction and exact polynomial residual weights make the reference computable. Continuous-regression BLUE, Gaussian projection and Hermite identities are established results [@dette2016; @isserlis1918]; the contribution is their precise connection through the ODE adjoint range, including the boundary restriction and the order of limits. Two experiments test the resulting predictions.

'''+front[b:]
    front=front.replace('Continuous-time regression with correlated errors has an established optimal-design theory.',
        'Continuous-time regression with correlated errors has an established optimal-design theory, including optimal signed measures with endpoint masses.')
    front=front.replace('### 4.2. Fisher information and weak efficiency','### 4.2. Continuous-record information and scalar weak ratios')
    front=front.replace('If $\\theta\\mapsto x_\\theta$ is differentiable',
        '**Proposition 1 (continuous-record information).** If $\\theta\\mapsto x_\\theta$ is differentiable',1)
    front=front.replace('Its inverse is the Cram',
        'The Gaussian mean-score identity proves the assertion: the directional score has variance equal to the squared Cameron-Martin norm of its sensitivity. Its inverse is the Cram',1)
    pieces=[]
    headings=[('## 1. Spaces and adjoint range','### 4.3. Exact boundary-free influence class'),
        ('## 2. Sharp covariance reference','### 4.4. Sharp information and optimal influence'),
        ('## 3. Explicit OU representative','### 4.5. Explicit optimal OU test'),
        ('## 4. GLS and the order of limits','### 4.6. Correct-covariance GLS and dense spans')]
    body=theory[theory.index('## 1.'):theory.index('## References')]
    for source,target in headings:body=body.replace(source,target)
    covariance=old[old.index('## 5.'):old.index('## 6.')]
    covariance=covariance.replace('Order $q$ has correlation time',
        'The order-one covariance supplies the GLS weights required by Corollary 1. Higher orders refine residual weights at finite amplitude; they do not alone remove coefficient bias. Order $q$ has correlation time')
    design='''## 6. Sampling and the scope of efficiency

Theorem 1 applies at fixed physical noise covariance and fixed horizon. Corollary 1 concerns a different limit: decreasing noise for each fixed test family, followed by dense test spans. Increasing the grid density improves the quadrature approximation; it does not replace either limit. Increasing the number of tests improves approximation of an admissible influence but cannot recover the initial-state information excluded by Proposition 2.

Independent records add their sensitivity information. Longer records add useful information only while parameter sensitivities remain appreciable. For scalar decay, information approaches a finite limit as duration increases. These distinctions and numerical acquisition and aggregation checks are detailed in S9. Averaging corrected weak matrices and stacking noisy regressions are different estimating equations; the latter can retain errors-in-variables bias [@fuller1987].

## 7. Two tests of the efficiency prediction

Both experiments use known covariance shape and supplied identifiable structure. They test feasible weak-estimation components, separately from the complete author pipelines archived in S6 and S9.6. No true parameter is supplied to a feasible pilot. Deterministic geometry calculations use true parameters only as local references. All computations execute on Kaggle; protocols, seeds, source fingerprints and output archives are preserved.

@@NUMERICAL@@

## 8. Discussion and limitations

The adjoint-range characterization turns the initial-sensitivity constraint into an exact description of the closed influence class. Its sharp optimum follows from classical Gaussian projection, and the OU inverse covariance makes an optimal generalized test explicit. Endpoint atoms explain why compact tests generally approach the optimum without attaining it. The constructed scalar cutoff uses compact tests throughout; this article does not test an estimator with additional boundary terms.

Correct covariance weighting is optimal within a fixed span. Density in the covariance norm is required to approach the class bound. A well-conditioned representation changes computation, not the underlying span. An explicit fundamental-matrix inverse may itself be ill-conditioned for a contracting dynamical direction. The multiscale calculation therefore uses a fixed, data-independent basis transformation and reports retained dimensions and numerical relaxation.

Corollary 1 requires consistent pilots, correct covariance shape, a regular local estimator, controlled quadrature and risk conditions. Its limits are sequential. It does not prove a joint window-growth/noise-decay limit, global convergence of covariance iterations, efficiency after sparse support selection, or efficiency of the released iid WENDy implementation under OU noise. A finite nested family does not establish density of an infinite family. The iid comparison tests a specified misspecification and does not establish a universal limiting plateau.

Exact residual covariance and marginal Hermite correction address distinct parts of the weak errors-in-variables problem. At finite noise, random design and pilot errors can still produce bias and excess dispersion. The numerical comparison therefore reports RMSE, bias, dispersion, stopping diagnostics and uncertainty, rather than equating a median error with a covariance bound. Supplementary complete-pipeline results use a different 16-path protocol and are not pooled with the new 128-path experiment.

The analysis assumes fully observed deterministic ODEs and additive stationary Gaussian measurement error with known positive covariance common across parameters. Unknown covariance, partial observation and process noise are outside the theorem proved here. The scope is an information reference and its weak approximation, rather than a general ranking of identification algorithms.

'''
    disclosures=old[old.index('## Declaration of'):]
    (BASE/'efficiency_article_template.md').write_text(front+body+'\n'+covariance+design+disclosures,encoding='utf-8')
    # Supporting evidence remains in a stable source instead of being discarded.
    supporting=old[old.index('### 7.1.'):old.index('## 8.')]
    supporting=re.sub(r'<!--.*?-->','',supporting)
    supporting=re.sub(r'### 7\.(\d+)\.',r'### S9.\1.',supporting)
    supporting=re.sub(r'\bSection 7\.(\d+)',r'Section S9.\1',supporting)
    supporting=re.sub(r'\bSections 7\.(\d+)',r'Sections S9.\1',supporting)
    supporting=re.sub(r'\bTable (\d+)',r'Table S\1',supporting)
    supporting=re.sub(r'\bFigure (\d+)',r'Figure S\1',supporting)
    supporting=supporting.replace('Section 4.3','Section 4.6')
    supporting=supporting.replace('no global shape optimum is claimed',
        'the compact-test infimum is now characterized by Theorem 2; these values concern only the listed finite spans')
    (BASE/'supporting_validation.md').write_text(supporting,encoding='utf-8')
    supplement=(BASE/'archive/pre_efficiency/supplement.md').read_text(encoding='utf-8')
    supplement=supplement[:supplement.index('## References')]
    supplement=re.sub(r'\bTable (\d+)',r'Table S\1',supplement)
    supplement=supplement.replace('Sections 7.7–7.9','Sections S9.7–S9.9').replace('Section 4.3','Section 4.6')
    supplement=supplement.replace('retained in the main text','retained in S9.8')
    (BASE/'efficiency_supplement_template.md').write_text(supplement,encoding='utf-8')
    print('Efficiency templates and supporting evidence prepared')


if __name__=='__main__':prepare()
