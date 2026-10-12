"""Prepare stable manuscript templates; never change execution snapshots."""
from pathlib import Path
import re,json
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'


def prepare():
    old=(BASE/'archive/pre_efficiency/article_revised.md').read_text(encoding='utf-8')
    theory=(BASE/'weak_efficiency_theory.md').read_text(encoding='utf-8')
    front=old[:old.index('### 4.3.')]
    front=front.replace('Fixed-horizon information limits and polynomial weak-residual covariance under correlated measurement error',
        'Efficient weak-form estimation with correlated noise: the cost of boundary-free tests')
    author=json.loads((ROOT/'configs/authors.json').read_text(encoding='utf-8'))
    front=front.replace('lang: en-US', 'author: "'+author['name']+'"\nlang: en-US')
    front=front.replace('## Abstract', author['affiliation']+'\n\nCorresponding author: '+author['email']+'\n\n## Abstract',1)
    a=front.index('## Abstract');b=front.index('**Keywords:**')
    front=front[:a]+'''## Abstract

Weak-form dynamical identification avoids differentiating noisy measurements, but removing boundary terms also discards information. For deterministic ODEs observed in stationary Ornstein-Uhlenbeck measurement error, we characterize the closure of boundary-free adjoint influences: it is exactly the covariance-dual subspace orthogonal to initial-state sensitivities. This gives a sharp local covariance bound equal to the inverse initial-state Schur complement, even when the initial state is known. We construct its optimal influence and an explicit generalized test function with endpoint jumps. Compact tests attain the bound as an infimum. Correct-covariance weak generalized least squares approaches this bound for dense test spans, in a sequential small-noise limit with controlled quadrature and regular pilots. A finite Hermite expansion supplies exact coupled polynomial residual weights. A fixed-horizon two-point bound separates this efficiency question from consistency under denser sampling. Experiments test a smoothed scalar optimal function with a fitted pilot and nested 30-, 120- and 480-window Van der Pol families, comparing physical OU and iid weights. Small-noise accuracy is compared with the information reference; finite-noise bias and dispersion are reported separately. The contribution identifies which information weak equations can retain and how covariance weighting and test-space approximation govern its recovery.

'''+front[b:]
    a=front.index('We develop an explicit');b=front.index('## 2.',a)
    front=front[:a]+'''Our main result identifies the complete boundary-free adjoint influence class and constructs the influence attaining its sharp information reference in covariance-norm closure. Correct temporal covariance and dense test spans then yield a conditional, sequential small-noise efficiency result for WENDy-like estimation. An explicit OU adjoint construction and exact polynomial residual weights make the reference computable. Continuous-regression BLUE, Gaussian projection and Hermite identities are established results [@dette2016; @isserlis1918]; the contribution is their precise connection through the ODE adjoint range, including the boundary restriction and the order of limits. Specifically, we prove (i) the exact closed weak class, (ii) its sharp information reference and an explicit OU optimal test, (iii) conditional dense-span efficiency of correct-covariance GLS with a concrete spline density lemma, and (iv) an estimator-independent fixed-horizon limit. Three dynamical systems and paired likelihood references test the resulting predictions.

'''+front[b:]
    recent='''WENDy-MLE optimizes an approximate Gaussian weak-residual likelihood with parameter-dependent covariance and analytic derivatives [@rummel2025]. Its likelihood contains both a log determinant and a covariance-weighted residual term; it is not identical to a frozen-weight GLS update. Tran and Bortz select local test radii by estimating quadrature error [@tran2026], while Fourier Weak SINDy selects informative sinusoidal frequencies spectrally [@chen2026]. These are relevant precedents for weak likelihoods and test design. The complementary question here is the complete adjoint influence class, its sharp initial-state information loss and the explicit OU optimal representative. We do not claim weak GLS or optimized test selection themselves as new methods.

Finite-noise bias and uncertainty are also studied directly for WENDy-IRLS [@chawla2026]. Our bias/dispersion reporting concerns the specified correlated-noise experiment and does not replace these broader coverage studies.

'''
    at=front.index('Errors-in-variables effects')
    front=front[:at]+recent+front[at:]
    lead='''Boundary-vanishing weak equations can lose information even when an initial state is known, as recognized by orthogonality-condition estimation [@brunel2014] and recent weak identifiability work [@heitzman2026]. We quantify this effect: the closed class of boundary-free weak influences has exactly the local information available when the initial state is treated as a nuisance parameter in the full path likelihood. For an unknown initial state, this restriction causes no additional first-order information loss; for a known initial state, the lost information is the initial-sensitivity projection, and the retained information is its Schur complement. This is a local information statement, not an identity between finite-noise estimators.

'''
    at=front.index('## 1.');at=front.index('\n\n',at)+2
    front=front[:at]+lead+front[at:]
    statistical='''Statistical ODE estimation has long used smoothing to reduce repeated forward solves: spline least squares [@varah1982], nonparametric gradient matching [@brunel2008], measurement-error regression [@liang2008] and root-$n$ smooth-and-match estimation [@gugushvili2012]. Direct integral estimators have root-$n$ consistency for suitable linear-in-function parameterizations [@dattner2015]. Orthogonality-condition estimation already supplies weak moments, GMM weighting and boundary-information adaptations [@brunel2014]. Thus neither the qualitative boundary-loss idea nor optimally weighted moments originates here.

Generalized profiling can be asymptotically efficient without explicitly solving the ODE, under its approximation and i.i.d. sampling assumptions [@ramsay2007; @qi2010]. One-step likelihood refinement of an ODE pilot is also established [@dattner2018]; Dattner and Klaassen [@dattner2015] discuss such refinements as a direction beyond their direct-estimator result. We do not propose replacing these efficient full-observation methods. The question is the optimum available to boundary-free weak equations, which also provide tractable moments for library-based structural discovery. The present theorem assumes supplied identifiable structure; the one-step computation below is a likelihood reference on the same correlated-noise records. These small supplied-structure examples remain tractable for likelihood refinement; computational superiority of weak estimation is not inferred from them.

Gaussian regression geometry and small-noise efficiency belong to classical statistical theory [@parzen1961; @ibragimov1981; @kutoyants1994]. Here those principles are connected to the dynamical adjoint range. The distinction from prior qualitative boundary observations is the exact covariance-dual class, its sharp Schur-complement reference and the explicit OU representative.

'''
    at=front.index('## 2.');at=front.index('\n\n',at)+2
    front=front[:at]+statistical+front[at:]
    front=front.replace('Two experiments test', 'The experiments test')
    front=front.replace('Continuous-time regression with correlated errors has an established optimal-design theory.',
        'Continuous-time regression with correlated errors has an established optimal-design theory, including optimal signed measures with endpoint masses.')
    front=front.replace('### 4.2. Fisher information and weak efficiency','### 4.2. Continuous-record information and scalar weak ratios')
    front=front.replace('If $\\theta\\mapsto x_\\theta$ is differentiable',
        '**Proposition 1 (continuous-record information).** If $\\theta\\mapsto x_\\theta$ is differentiable',1)
    front=front.replace('Its inverse is the Cram',
        'The Gaussian mean-score identity proves the assertion: the directional score has variance equal to the squared Cameron-Martin norm of its sensitivity. Its inverse is the Cram',1)
    pieces=[]
    headings=[('## 1. Spaces and adjoint range','## 5. Boundary-free weak class and efficiency\n\n### 5.1. Exact boundary-free influence class'),
        ('## 2. Sharp covariance reference','### 5.2. Sharp information and optimal influence'),
        ('## 3. Explicit OU representative','### 5.3. Explicit optimal OU test'),
        ('## 4. GLS and the order of limits','### 5.4. Correct-covariance GLS and dense spans')]
    body=theory[theory.index('## 1.'):theory.index('## References')]
    for source,target in headings:body=body.replace(source,target)
    covariance=old[old.index('## 5.'):old.index('## 6.')]
    covariance=covariance.replace('Order $q$ has correlation time',
        'The order-one covariance supplies the GLS weights required by Corollary 1. Higher orders refine residual weights at finite amplitude; they do not alone remove coefficient bias. Order $q$ has correlation time')
    uncertainty=covariance[covariance.index('### 5.3.'):].replace('### 5.3.','## S10.')
    uncertainty=uncertainty.replace(r'\varepsilon=\sigma Z', r'\eta=\epsilon\eta_C').replace(r'\sigma P_0', r'\epsilon P_0').replace(r'O_Z(\sigma^2)', r'O_{\eta_C}(\epsilon^2)')
    covariance=covariance[:covariance.index('### 5.3.')].replace(r'\varepsilon',r'\eta_C')
    covariance=covariance.replace('## 5.','## 6.').replace('### 5.','### 6.')
    design='''## 7. Sampling and the scope of efficiency

Theorem 1 applies at fixed physical noise covariance and fixed horizon. Corollary 1 concerns a different limit: decreasing noise for each fixed test family, followed by dense test spans. Increasing the grid density improves the quadrature approximation; it does not replace either limit. Increasing the number of tests improves approximation of an admissible influence but cannot recover the initial-state information excluded by Proposition 2.

Independent records add their sensitivity information. Longer records add useful information only while parameter sensitivities remain appreciable. For scalar decay, information approaches a finite limit as duration increases. These distinctions and numerical acquisition and aggregation checks are detailed in S9. Averaging corrected weak matrices and stacking noisy regressions are different estimating equations; the latter can retain errors-in-variables bias [@fuller1987].

## 8. Tests of the efficiency prediction

The principal comparisons use supplied identifiable structure and known covariance shape; an additional sensitivity check fits the noise shape from data. They test feasible weak-estimation components, separately from the complete author pipelines archived in S6 and S9.6. No true parameter is supplied to a feasible pilot. Deterministic geometry calculations use true parameters only as local references. Protocols, seeds, source fingerprints and output archives are preserved in the public reproducibility archive.

@@NUMERICAL@@

## 9. Discussion and limitations

The adjoint-range characterization turns the initial-sensitivity constraint into an exact description of the closed influence class. Its sharp optimum follows from classical Gaussian projection, and the OU inverse covariance makes an optimal generalized test explicit. Endpoint atoms explain why compact tests generally approach the optimum without attaining it. The constructed scalar cutoff uses compact tests throughout; this article does not test an estimator with additional boundary terms.

Correct covariance weighting is optimal within a fixed span. Density in the covariance norm is required to approach the class bound. A well-conditioned representation changes computation, not the underlying span. An explicit fundamental-matrix inverse may itself be ill-conditioned for a contracting dynamical direction. The multiscale calculation therefore uses a fixed, data-independent basis transformation and reports retained dimensions and numerical relaxation.

The efficiency conclusion applies to the local GLS map under the regularity and risk assumptions of Corollary 1. Global covariance-iteration behavior, sparse support selection and simultaneous test-growth/noise-decay limits require additional analysis. Lemma 1 supplies a dense spline family, while the finite bump experiment tests approximation at three declared dimensions. Its iid comparison concerns this specific misspecification.

Exact residual covariance and marginal Hermite correction address distinct parts of the weak errors-in-variables problem. At finite noise, random design and pilot errors can still produce bias and excess dispersion. The numerical comparison therefore reports RMSE, bias, dispersion, stopping diagnostics and uncertainty, rather than equating a median error with a covariance bound. The separate 16-path complete-pipeline benchmark is retained as context; the matched 128-path WENDy comparison and the 2000-path small-noise weighting study address distinct practical questions.

The analysis assumes fully observed deterministic ODEs and additive stationary Gaussian measurement error with known positive covariance common across parameters. Unknown covariance, partial observation and process noise are outside the theorem proved here. The feasible vector implementation uses one fitted-trajectory solve to stabilize covariance weights. It therefore tests a statistical efficiency mechanism and does not establish a forward-solve-free or runtime-optimal algorithm. The result provides a local information reference for boundary-free weak estimation. The matched finite-noise comparison and the full-pipeline benchmark illustrate its practical scope; they do not establish a general algorithm ranking.

'''
    disclosures=old[old.index('## Declaration of'):]
    front=front.replace('defined in Section 5:', 'defined in Section 6:')
    front=front.replace(r'\varepsilon_j', r'\eta_{C,j}').replace(r'\varepsilon', r'\eta_C')
    front=front.replace('Initial conditions are known unless stated otherwise; unknown initial conditions would introduce nuisance sensitivities.', 'The initial state $u=x(0)$ is known in the baseline acquisition experiment; Sections 5 and 8 compare known and unknown $u$ explicitly. The symbol $Z$ denotes its fundamental sensitivity matrix, $h$ the observation mesh, $\\epsilon$ a scalar noise amplitude, and $C$ the continuous covariance operator. Finite weak covariance matrices are denoted by $\\Sigma$. In Corollary 1, $C$ is the unit-amplitude covariance shape and the physical covariance is $\\epsilon^2 C$.')
    front=front.replace('## 4. An estimator-independent fixed-horizon limit','## 4. Fixed-horizon information limits')
    (BASE/'local_uncertainty_appendix.md').write_text(uncertainty,encoding='utf-8')
    (BASE/'efficiency_article_template.md').write_text(front+body+'\n'+covariance+design+disclosures,encoding='utf-8')
    # Supporting evidence remains in a stable source instead of being discarded.
    supporting=old[old.index('### 7.1.'):old.index('## 8.')]
    supporting=re.sub(r'<!--.*?-->','',supporting)
    supporting=re.sub(r'### 7\.(\d+)\.',r'### S9.\1.',supporting)
    supporting=re.sub(r'\bSection 7\.(\d+)',r'Section S9.\1',supporting)
    supporting=re.sub(r'\bSections 7\.(\d+)',r'Sections S9.\1',supporting)
    supporting=re.sub(r'\bTable (\d+)',r'Table S\1',supporting)
    supporting=re.sub(r'\bFigure (\d+)',r'Figure S\1',supporting)
    supporting=supporting.replace('Section 4.3','Section 5.4')
    supporting=supporting.replace('no global shape optimum is claimed',
        'the compact-test infimum is now characterized by Theorem 2; these values concern only the listed finite spans')
    (BASE/'supporting_validation.md').write_text(supporting,encoding='utf-8')
    supplement=(BASE/'archive/pre_efficiency/supplement.md').read_text(encoding='utf-8')
    supplement=supplement[:supplement.index('## References')]
    supplement=re.sub(r'\bTable (\d+)',r'Table S\1',supplement)
    supplement=supplement.replace('Sections 7.7–7.9','Sections S9.7–S9.9').replace('Section 4.3','Section 5.4')
    supplement=supplement.replace('retained in the main text','retained in S9.8')
    (BASE/'efficiency_supplement_template.md').write_text(supplement,encoding='utf-8')
    print('Efficiency templates and supporting evidence prepared')


if __name__=='__main__':prepare()
