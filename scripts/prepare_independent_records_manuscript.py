"""Extend the rebuilt manuscript with replicate inference and factual author data."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'
TITLE='Weak-form estimation with correlated noise: information limits and independent-record inference'
OLD='Efficient weak-form estimation with correlated noise: the cost of boundary-free tests'


def main():
    p=BASE/'article_revised.md';text=p.read_text(encoding='utf-8');supp=(BASE/'supplement.md').read_text(encoding='utf-8')
    text=text.replace(OLD,TITLE)
    a=text.index('## Abstract');b=text.index('**Keywords:**',a)
    text=text[:a]+'''## Abstract

Boundary-free weak equations avoid differentiating noisy measurements but constrain the information available for dynamical estimation. For deterministic ODEs observed in stationary Ornstein-Uhlenbeck error, we characterize their closed adjoint influence class and its sharp local covariance bound: the inverse initial-state Schur complement. We construct an explicit generalized OU test with endpoint jumps; compact tests approach the bound as an infimum. Correct-covariance weak GLS is efficient within this class in a sequential small-noise, dense-span limit under stated regularity conditions. A finite Hermite expansion supplies exact polynomial residual covariance. In a separate regime of independent records at fixed noise, corrected moments yield consistent fixed-family GMM and square-root-record-count inference. Estimated correction variances contribute a nuisance influence to its sandwich covariance. Exact integration or controlled vanishing quadrature defects are required for the continuous-parameter asymptotics. Scalar, Van der Pol and Lorenz experiments examine test-space approximation, finite-noise bias and likelihood references. On 2000 shared small-noise Van der Pol paths, OU weighting lowers RMSE by 4.2% relative to iid weights. Repeated-record experiments at 5% noise assess estimated-covariance GMM and 95% interval coverage. The analysis separates single-record infill limits, local efficiency and replication-based inference.

'''+text[b:]
    comparison='''### 2.1. Complementarity with orthogonality-condition estimation

Brunel et al. [@brunel2014] use a nonparametric trajectory proxy, orthogonality conditions and GMM in an independent-observation, increasing-point experiment, with a formulation extending to delay equations. They already recognize boundary-information loss and describe adaptations using known boundary values. Our fully observed deterministic polynomial setting is narrower in model structure. It instead describes physical OU noise on a fixed horizon, the exact closed adjoint influence class, its sharp Schur-complement reference and an explicit optimal OU test. The small-noise dense-span result and the fixed-noise independent-record result below concern different asymptotic experiments from their increasing-point analysis. These results complement each other; neither a first use of weak GMM nor a strict dominance claim is intended.

'''
    at=text.index('## 3.');text=text[:at]+comparison+text[at:]
    section=r'''### 6.3. Fixed-noise inference from independent records

Let $n$ independent records share the same deterministic ODE trajectory on $[0,T]$. For a fixed identifiable polynomial structure and fixed boundary-free test family, the Hermite-corrected moments $m_r(\theta,\nu)=b_r-G_r(\nu)\theta$ satisfy $E m_r(\theta_0,\nu_0)=0$ at every finite noise amplitude, using exact integrals. This population moment identity, rather than an unbiased single-record ratio, permits replication-based inference.

**Theorem 3 (independent-record GMM; full statement and proof in S14).** Suppose $D=E G_r(\nu_0)$ has full column rank, moment covariance $\Omega$ is positive definite and the independent-record moments have finite second moments. For known covariance parameters and consistent weights $W_n\to W$, the closed-form estimate

$$
\widehat\theta_n=(\bar G_n^\top W_n\bar G_n)^{-1}\bar G_n^\top W_n\bar b_n
$$

is consistent and $\sqrt n(\widehat\theta_n-\theta_0)\Rightarrow N(0,P_W\Omega P_W^\top)$, where $P_W=(D^\top W D)^{-1}D^\top W$. The optimal fixed-family weight is $\Omega^{-1}$, and Proposition 3 gives this covariance with all Hermite orders. A consistent pilot and consistent covariance estimate yield feasible two-step GMM. This is an application of established GMM theory [@hansen1982; @newey1994].

For a regular same-data covariance estimate with influence $\xi_r$, correction parameters contribute through $H=E\partial_\nu m_r$: the effective covariance is $\widetilde\Omega=\operatorname{Cov}(m_r+H\xi_r)$. Its inverse is optimal for this nuisance-estimation scheme, and its sandwich gives asymptotic Wald coverage. Variance estimation inside Hermite correction cannot generally be ignored. Pair differences of records remove their shared mean and give consistent OU variances and range at fixed horizon, with independent pairs as inference clusters (S14.1). This replication result is distinct from lack of separate-parameter consistency under fixed-horizon infill of one OU path.

On a numerical grid the moment identity has a deterministic quadrature defect. Continuous-parameter consistency requires vanishing defect, and centered $\sqrt n$ inference requires that defect be $o(n^{-1/2})$. S14 states the fixed-family condition and the experiment reports its clean-grid error. Stacking noisy designs instead involves $E[G_r^\top m_r]$, which need not vanish even when $E m_r=0$; this explains the aggregation bias in Table S11. None of these statements guarantees finite-noise likelihood efficiency.

Replication uses averaged estimating equations, not an average of separately fitted coefficients. It therefore resolves consistency in this acquisition regime without removing the finite-noise bias of an arbitrary single-record estimator.

### 6.4. Temporal misspecification

For a fixed moment family with correct marginal correction, the independent-record variance under weights $W$ is $P_W\Omega_{\mathrm{true}}P_W^\top/n$. The same sandwich describes the leading small-noise single-record variance using the order-one residual covariance. Incorrect temporal weights can reduce efficiency while preserving the corrected mean identity; incorrect correction variances can invalidate that identity. In S12, assuming $\tau/2$ or $2\tau$ changes only weighting, whereas fitting variances changes correction as well. Replication allows consistent estimation of both OU range and variances under the shared-trajectory construction.

'''
    at=text.index('## 7.');text=text[:at]+section+text[at:]
    text=text.replace('A finite Hermite expansion supplies exact coupled polynomial residual weights.',
        'A finite Hermite expansion supplies exact polynomial residual weights. For independent recordings, exactly corrected moments yield consistent fixed-family GMM and square-root-record-count inference at fixed noise; estimating correction variances adds a nuisance influence to the sandwich covariance.')
    text=text.replace('Three dynamical systems and paired likelihood references test the resulting predictions.',
        'Three dynamical systems, paired likelihood references and repeated-record coverage experiments test the resulting predictions.')
    text=text.replace('Unknown covariance, partial observation and process noise are outside the theorem proved here.',
        'Unknown covariance is treated for independent copies in Theorem 3; the single-record efficiency theorem assumes known shape. Partial observation and process noise remain outside the analysis.')
    text=text.replace('Global covariance-iteration behavior, sparse support selection and simultaneous test-growth/noise-decay limits require additional analysis.',
        'Global covariance-iteration behavior and sparse support selection require additional analysis. A simultaneous limit with growing test spans and decreasing noise, including quantitative rates controlling Gram conditioning and remainders, remains open.')
    start=text.index('OpenAI Codex assisted');end=text.index('## Data availability',start)
    text=text[:start]+'''OpenAI Codex was used for literature searches, mathematical derivations, code development, experiment orchestration, numerical interpretation and drafting. Simulations were executed by archived source code, and figures were generated from recorded numerical outputs. This disclosure covers assistance in the research and writing process; AI tools are not credited as authors.

'''+text[end:]
    text=text.replace('The repository includes protocol configurations,',
        'For the independent-record study, all coefficient estimates, covariance summaries and sufficient-statistic fixtures are archived; full time series are regenerated from the frozen seed protocol. The repository includes protocol configurations,')
    text=text.replace('![Figure 4. Van der Pol likelihood refinement and weak weighting on the same 2000 paths at 0.5% noise. Bars are 95% independent-path bootstrap intervals.',
        '![Figure 4. Van der Pol on 2000 shared paths at 0.5% noise: one scoring step with known initial state, one or three joint steps with unknown initial state, and weak GLS with OU or iid weights. Whiskers give 95% bootstrap intervals for RMSE over independent paths, not parameter confidence intervals or error standard deviations.')
    proof=(BASE/'independent_records_theory.md').read_text(encoding='utf-8')
    proof=proof[proof.index('## S14.'):proof.index('## References')]
    supp=supp[:supp.index('## References')]+proof+'\n## References\n'
    report=ROOT/'reports/independent_records_results.json'
    if report.exists():
        data=json.loads(report.read_text(encoding='utf-8'));groups=data['groups']
        exp=r'''### 7.6. Independent-record GMM and interval calibration at 5% noise

One thousand independent Monte Carlo repetitions each generate 800 independent Van der Pol records sharing $u=(2,0)$, $T=20$, $h=0.01$ and $\tau=0.3$. Counts 50, 100, 200, 400 and 800 are prefixes within each repetition. Six fixed width-four tests give 12 moments for four coefficients. One two-step update uses empirical pair-influence covariance, including the estimated-variance contribution, with a diagonal relaxation equal to $1/(n/2)$ times the empirical diagonal. Pair differences estimate two variances and the common OU correlation time. Empirical weights estimate temporal dependence directly from replicated moment vectors; fitted tau is reported as a covariance-identification diagnostic and does not enter these weights. All repetitions are retained; no numerical failure is observed.

'''
        geom=data['geometry'];exp+=f"The clean-grid relative coefficient defect is {geom['clean_gls_relative_error']:.2e}, compared with the reported sampling error. For known correction variance the analytical square-root-$n$ relative SD is {geom['known_variance_gmm_sqrt_n_relative_sd']:.6f}; higher Hermite orders contribute {100*geom['higher_order_trace_fraction']:.3f}% of the residual covariance trace. Empirical full residual covariance agrees with Proposition 3 to relative Frobenius error {data['empirical_full_residual_covariance_relative_error']:.4f} across {data['independent_records_simulated']:,} independent records.\n\n"
        exp+=r'''Table 6. Van der Pol averaged-moment GMM with estimated noise at 5%, 1000 repetitions per count. Coverage refers to nominal 95% normal marginal intervals and a chi-square four-parameter ellipsoid. These finite-sample coverage estimates test the asymptotic theorem; nominal coverage is not guaranteed at every n. C1-C4 denote marginal coverage for the four coefficients. Binomial Wilson intervals and known-variance results are in S15.

| n | RMSE | $\sqrt n$ RMSE | C1 | C2 | C3 | C4 | Joint | Stacked RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''
        for n in data['protocol']['record_counts']:
            r=next(r for r in groups if r['n']==n and r['method']=='gmm_estimated');stack=next(r for r in groups if r['n']==n and r['method']=='stacked_known')
            cov=' | '.join(f'{100*v:.1f}%' for v in r['marginal_coverage'])
            exp+=f"| {n} | {r['relative_rmse']:.6f} | {r['sqrt_n_rmse']:.5f} | {cov} | {100*r['joint_coverage']:.1f}% | {stack['relative_rmse']:.6f} |\n"
        first=next(r for r in groups if r['n']==50 and r['method']=='gmm_estimated');last=next(r for r in groups if r['n']==800 and r['method']=='gmm_estimated')
        exp+=f"\nThe square-root-$n$ scaled RMSE decreases from {first['sqrt_n_rmse']:.3f} at 50 records to {last['sqrt_n_rmse']:.3f} at 800, approaching the fixed-family reference {geom['known_variance_gmm_sqrt_n_relative_sd']:.3f}. This supports the large-record scaling while exposing finite-sample weight and random-design effects. At 50 records the joint ellipsoid covers only {100*first['joint_coverage']:.1f}%; known-variance GMM has similarly poor coverage, so nuisance substitution alone does not explain it. At 800 records joint coverage is {100*last['joint_coverage']:.1f}% with Wilson interval [{100*last['joint_coverage_wilson_95'][0]:.1f}%, {100*last['joint_coverage_wilson_95'][1]:.1f}%], including 95%. Estimated tau has median {last['tau_quantiles'][1]:.5f} and central 90% range [{last['tau_quantiles'][0]:.5f}, {last['tau_quantiles'][2]:.5f}].\n\n"
        exp+='This six-test family is fixed for inference and is much less informative than the 480-window single-record family. It is not a dense-span efficiency experiment. Its stacking bias floor is larger than in Table S11 because the test family and acquisition protocol differ; the common mechanism is correlation between noisy designs and residuals. Higher-order covariance has a small numerical share here, but the fixed-noise theorem uses the full moment covariance rather than discarding it.\n\n'
        exp+='\n![Figure 5. Independent-record Van der Pol inference at fixed 5% noise. Left: coefficient RMSE with 95% repetition-bootstrap intervals; known-variance stacked regression is compared with averaged-moment GMM. The two GMM curves nearly coincide. Right: marginal and joint coverage with estimated variances and range; the dotted line is the nominal 95% level. Counts share records within each repetition. Coverage uncertainty is tabulated in S15.](figures/publication/independent_records.png)\n\n'
        at=text.index('## 8.');text=text[:at]+exp+text[at:]
        additional=r'''## S15. Independent-record protocol and finite-sample calibration

The protocol is frozen in `independent_records.json` and in each execution notebook before these outputs. Four shards partition repetition indices modulo four; seeds are 2700000-2700999. Every repetition simulates 800 independent records with common trajectory; record counts are prefixes, not independent experiments. Gaussian innovations and exact AR(1) recursion generate stationary OU observations at h=0.01. Marginal variances equal (0.05 times each clean-coordinate SD) squared. Known-variance GMM and stacking receive that variance; estimated-variance GMM receives only observations and the declared test functions. No true coefficient or clean trajectory enters fitted GMM weights. The clean path is used only for generating data and the separate covariance/quadrature audit. Twelve moments are fixed throughout: centers [2.1,5.2,8.3,11.4,14.5,17.6], support width four, polynomial bump order four.

Weights are the inverse empirical covariance of independent-pair influences at the initial averaged-moment ordinary least-squares estimate. The diagonal relaxation is empirical diagonal divided by the number of pairs; it vanishes asymptotically. One closed-form weighted update gives the final coefficients. The final residuals and variance-estimation influence supply a sandwich under the same initial weights. Standard errors use n/2 pair clusters, normal marginal critical value 1.959963984540054 and joint chi-square(4) quantile 0.95. Marginal intervals are individually 95%, not a simultaneous rectangular region. The ellipsoid supplies the separate joint coverage. No Student-t, bootstrap parameter interval or outcome-tuned finite-sample correction is used. Coverage Monte Carlo uncertainty uses binomial Wilson intervals; n-values and methods remain paired within repetition. The RMSE bootstrap uses 2000 independent-repetition resamples with seed 2702600. Finite empirical RMSE does not establish finite population moments of an unbounded random-design inverse.

Table S14. Independent-record known/estimated-variance averaged-moment GMM. The last column is the minimum-to-maximum marginal coverage across four coefficients; all individual coverages and their Wilson intervals are archived. All counts have 1000 retained repetitions.

| n | Correction | RMSE | Bias | Dispersion | Joint coverage (95% Wilson interval) | Marginal range |
|---:|---|---:|---:|---:|---|---|
'''
        for r in groups:
            if r['method']=='stacked_known':continue
            low,high=r['joint_coverage_wilson_95'];cov=f"{100*min(r['marginal_coverage']):.1f}-{100*max(r['marginal_coverage']):.1f}%"
            additional+=f"| {r['n']} | {'estimated' if r['method']=='gmm_estimated' else 'known'} | {r['relative_rmse']:.6f} | {r['relative_bias']:.6f} | {r['relative_dispersion']:.6f} | {100*r['joint_coverage']:.1f}% [{100*low:.1f}, {100*high:.1f}] | {cov} |\n"
        additional+='\nTable S15. OU range estimation from independent record pairs, the same repetitions. True tau=0.3. Variance error, marginal coverage Wilson intervals and standard-error calibration are recorded in the numerical report.\n\n| n | tau 5%, 50%, 95% quantiles | Marginal variance RMSE, two coordinates |\n|---:|---|---|\n'
        for r in groups:
            if r['method']=='gmm_estimated':additional+=f"| {r['n']} | {', '.join(f'{v:.5f}' for v in r['tau_quantiles'])} | {', '.join(f'{v:.2e}' for v in r['variance_rmse'])} |\n"
        additional+='\nAll sufficient-statistic input hashes and a complete first-repetition sufficient-statistic fixture per shard are archived. Aggregate residual sums and second moments permit a separate exact-covariance comparison without retaining every time series. Re-running the frozen seed protocol reconstructs every synthetic record.\n\n'
        supp=supp.replace('## References',additional+'## References')
        text=text.replace('finite-noise bias and dispersion are reported separately.', 'finite-noise bias and dispersion are reported separately.')
    start=text.index('## 8. Discussion');end=text.index('## Declaration',start)
    text=text[:start]+r'''## 8. Discussion and limitations

The exact adjoint range connects boundary-free weak equations to the initial-state nuisance geometry. Classical Gaussian projection then gives its sharp covariance reference, and the OU inverse covariance supplies an explicit optimal generalized test. For unknown initial state there is no additional first-order weak-class loss; known initial state permits information excluded by this class. Compact tests approach the optimum through endpoint layers. Lemma 1 proves spline density; the finite bump families are numerical approximations whose retained dimensions and conditioning remain relevant.

Single-record efficiency and independent-record consistency address different limits. Corollary 1 is sequential and requires regular pilots, controlled quadrature and risk assumptions. Exact Hermite correction and residual covariance do not remove finite-noise random-design bias from one fit. The matched WENDy result accordingly gives no evidence of superior practical accuracy of the weak component at 5%. Averaged-moment GMM instead has a fixed-noise record CLT, including fitted-variance influence. Its interval experiment shows strong undercoverage at small record counts and near-nominal coverage at 800; asymptotic coverage is not a finite-sample guarantee.

A simultaneous growing-span and decreasing-noise limit, with rates controlling Gram conditioning and remainders, remains open. Partial observation, process noise, global covariance-iteration behavior and support selection also need separate analysis. The results concern fully observed deterministic ODEs with supplied identifiable structure. The feasible vector implementation uses a fitted-trajectory solve for covariance weights, and no runtime superiority is inferred. These scope restrictions accompany a statistical reference for designing weak tests and for deciding when independent acquisition changes the inference problem.

'''+text[end:]
    p.write_text(text,encoding='utf-8');(BASE/'supplement.md').write_text(supp,encoding='utf-8')
    cover=BASE/'submission/cover_letter_draft.md';s=cover.read_text(encoding='utf-8').replace(OLD,TITLE)
    s=s.replace('The paper does not claim universal algorithm rankings or convergence at fixed noise.',
        'A separate independent-record theorem applies classical GMM theory to exactly corrected moments at fixed noise, including fitted-covariance influence and interval calibration. The paper does not claim universal algorithm rankings or single-record consistency from denser sampling.')
    cover.write_text(s,encoding='utf-8')
    print(json.dumps(dict(author_added=True,replicate_theorem_added=True,replicate_results_added=report.exists())))


if __name__=='__main__':main()
