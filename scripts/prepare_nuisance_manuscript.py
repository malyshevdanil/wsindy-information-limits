"""Put the nuisance-initial-state thesis and EJS metadata into the rebuilt paper."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'
TITLE='Boundary-free weak estimating equations under correlated noise: efficiency with nuisance initial states'


def main():
    p=BASE/'article_revised.md';text=p.read_text(encoding='utf-8');supp=(BASE/'supplement.md').read_text(encoding='utf-8')
    if TITLE in text:raise RuntimeError('Rebuild preceding manuscript stages before applying this transformation')
    old=re.search(r'title: "([^"]+)"',text)[1];text=text.replace(old,TITLE)
    text=text.replace('Corresponding author: malyshevdanil04@gmail.com','Corresponding author: malyshevdanil04@gmail.com\n\nORCID: https://orcid.org/0009-0007-0947-3866\n\n**MSC 2020:** Primary 62F12; secondary 62M10, 34A55.')
    start=text.index('## Abstract');end=text.index('**Keywords:**',start)
    text=text[:start]+'''## Abstract

Boundary-free weak equations eliminate initial-state directions from dynamical estimating equations. We show that this restriction is locally efficient when initial states are nuisance parameters and quantify its cost when they are known. For deterministic ODEs observed in stationary Ornstein-Uhlenbeck error, the closed adjoint influence class has the sharp covariance reference given by the inverse initial-state Schur complement. An explicit generalized optimal test has endpoint jumps, which compact tests approximate through boundary layers. Correct-covariance GLS approaches this reference in a sequential small-noise, dense-span limit. Discarding known initial-state information raises the local RMSE reference by 29% for Van der Pol and almost doubles it for Lorenz. On 2000 shared Van der Pol records, known-initial likelihood refinement has 23% lower RMSE than weak GLS, while OU weighting improves on iid weights by 4.2%. Exact Hermite covariance supports model weights for records with different initial states. Independent instruments preserve corrected weak moment identities at fixed noise and permit record-level sandwich inference without fitting every initial state in the estimating equation. Numerical tests assess span approximation, finite-noise accuracy and interval calibration. The contribution is a statistical theory and reproducible benchmark for weak estimation.

'''+text[end:]
    start=text.index('Data acquisition for dynamical identification');end=text.index('Sparse identification',start)
    text=text[:start]+'''The nuisance interpretation also connects single-record geometry to repeated measurements with different initial conditions. Each new record introduces its own initial state, while a common coefficient vector remains the target. Corrected weak identities hold conditionally on every initial state. Their local product-experiment information is the sum of the record-specific Schur complements. At finite noise, valid instruments must preserve those mean identities; merely inserting a noisy design into an otherwise optimal weight recreates errors-in-variables bias.

'''+text[end:]
    start=text.index('Our main result identifies');end=text.index('## 2.',start)
    text=text[:start]+'''The contribution is a theory and benchmark for boundary-free estimation with nuisance initial states. Theorems 1 and 2 give the fixed-horizon obstruction and the exact adjoint class with its sharp information reference. The OU construction supplies an explicit optimal generalized test; Corollary 1 explains sequential dense-span GLS efficiency. Proposition 3 supplies full polynomial residual covariance for computing model weights. Theorem 3 extends the conditional weak identity to heterogeneous records, using independent instruments to obtain fixed-noise consistency and a record-level sandwich. Its local information reference is the sum of the individual nuisance-adjusted informations. Established Gaussian projection, BLUE, GMM and data-fission theory provide the statistical tools [@dette2016; @hansen1982; @newey1994; @leiner2025]; the precise adjoint restriction and its quantitative dynamical cost form the central result.

'''+text[end:]
    # Shorten Theorem 1; retain its complete proof in the likelihood appendix.
    start=text.index('**Theorem 1 ');proofstart=text.index('**Proof.**',start);end=text.index('### 4.2.',proofstart)
    oldproof=text[proofstart:text.index('A constant estimator',proofstart)]
    theorem=r'''**Theorem 1 (fixed-horizon two-point bound; proof in S1.1).** Suppose the mean trajectories belong to $\mathcal H$ and share the stated positive OU covariance. For $\delta=\|\theta-\eta\|>0$ and $d=\|x_\theta-x_\eta\|_{\mathcal H}$, every finite-grid estimator satisfies

$$
\max_{v\in\{\theta,\eta\}}P_v\{\|\widehat\theta_h-v\|\geq\delta/2\}
\geq\Phi(-d/2)>0.
$$

Its maximal two-point squared risk is at least $(\delta^2/4)\Phi(-d/2)$. Thus dense-grid consistency throughout a nontrivial fixed-horizon family is excluded by Gaussian path equivalence [@feldman1958; @zhang2004].

'''
    text=text[:start]+theorem+text[end:]
    at=supp.index('## S2.');supp=supp[:at]+'### S1.1. Proof of Theorem 1\n\n'+oldproof+'\n'+supp[at:]
    # The old common-trajectory GMM theorem becomes a supporting lemma.
    start=text.index('### 6.3. Fixed-noise');end=text.index('### 6.4.',start)
    text=text[:start]+text[end:];text=text.replace('### 6.4.','### 6.3.')
    start=text.index('### 7.6.');end=text.index('## 8.',start)
    common_experiment=text[start:end];text=text[:start]+text[end:]
    common_experiment=common_experiment[common_experiment.index('![Figure 5.'):]
    common_experiment=common_experiment.replace('Figure 5.','Figure S8.')
    supp=supp.replace('## S14. Independent-record theorem and proof','## S14. Common-trajectory GMM reference and proof')
    s14=supp.index('## S14.');s15=supp.index('## S15.',s14)
    supp=supp[:s14]+supp[s14:s15].replace('Theorem 3','Lemma 2')+supp[s15:]
    supp=supp.replace('No Student-t, bootstrap parameter interval or outcome-tuned finite-sample correction is used.',
        'These rows use asymptotic critical values; S17 evaluates a separately frozen recentered bootstrap correction on identical records.')
    supp=supp.replace('Table S14. Independent-record known/estimated-variance',
        'Table S14. Common-trajectory known/estimated-variance')
    supp=supp.replace('## References',common_experiment+'\n## References')
    # Renumber the former experimental section before adding the new theory section.
    text=re.sub(r'(?m)^### 7\.(\d+)\.',r'### 8.\1.',text)
    text=text.replace('## 7. Tests','## 8. Tests').replace('## 8. Discussion','## 9. Discussion')
    for k in range(1,7):
        text=text.replace(f'Section 7.{k}',f'Section 8.{k}');supp=supp.replace(f'Section 7.{k}',f'Section 8.{k}')
    theory=(BASE/'heterogeneous_records_theory.md').read_text(encoding='utf-8')
    proof=theory[theory.index('## S16.'):theory.index('## References')]
    supp=supp.replace('## References',proof+'\n## References')
    section=r'''## 7. Independent records with nuisance initial states

### 7.1. Conditional moment identities and independent instruments

Let independent records have means $x_{\theta,u_r}$, with a common coefficient vector and a separate unknown initial state for each record. Hermite-corrected weak moments obey $E\{m_r(\theta_0,\nu_0)\mid u_r\}=0$. Consequently the initial-state distribution can remain unspecified. At finite noise, a same-record pilot used in $D_r^\top\Omega_r^{-1}$ can correlate with the residual and spoil that identity. We use an explicit independent-instrument construction to state a feasible fixed-noise result.

For known OU covariance $C$, draw independent $\zeta_r\sim N(0,C)$ and set $Y_r^\pm=Y_r\pm\zeta_r$. These independent Gaussian views have covariance $2C$ [@leiner2025]. Instruments $L_r^\pm$ use their own view and independent training data; moments use the opposite view with its doubled correction variance. Define $B_r=(L_r^+G_r^-+L_r^-G_r^+)/2$ and $a_r=(L_r^+b_r^-+L_r^-b_r^+)/2$.

**Assumptions (A1)–(A5).** Trajectories and sensitivities are regular on a compact initial-state set; records and randomizations are independent with known positive covariance; the instruments obey the stated view separation and uniform moment bounds; the average derivative has a nonsingular limit and the centered record scores satisfy a Lindeberg CLT; the average instrumented quadrature defect is $o(n^{-1/2})$. The precise numbered conditions are in S16.1.

**Theorem 3 (heterogeneous-record weak estimating equations).** Under (A1)–(A5), conditional on the independent training data,

$$
\widehat\theta_n=\left(\sum_r B_r\right)^{-1}\sum_r a_r
\xrightarrow{P}\theta_0,\qquad
\sqrt n(\widehat\theta_n-\theta_0)
\Rightarrow N(0,B^{-1}V_\psi B^{-\top}),
$$

where $B=\lim n^{-1}\sum E B_r$ and $V_\psi=\lim n^{-1}\sum\operatorname{Var}(a_r-B_r\theta_0)$. The empirical covariance of the coefficient-dimensional record scores consistently estimates its sandwich. The proof uses $E(L_r^+m_r^-\mid u_r)=0$ and the symmetric identity, followed by the standard estimating-equation CLT [@hansen1982; @newey1994]; see S16.1.

### 7.2. Model weights and the product information reference

For deterministic clean-trajectory instruments, $D_r=E G_r$ and Proposition 3's full covariance $\Omega_r$ yield $L_r=D_r^\top\Omega_r^{-1}$. The oracle fixed-family covariance is $(\sum_rD_r^\top\Omega_r^{-1}D_r)^{-1}$. The feasible experiment estimates these instrument maps on independent training trajectories. Its sandwich estimates only a four-dimensional record-score covariance, even when hundreds of weak moments supply the within-record model weights.

The joint local information for $\theta,u_1,\ldots,u_n$ has separate nuisance blocks. Taking their Schur complement gives

$$
\mathcal I_{e,n}=\sum_r\mathcal I_e(u_r).
$$

Theorem 2 gives this sharp boundary-free reference record by record. For fixed regular records, dense test spans recover it in the sequential small-noise limit. Symmetric view separation preserves the leading original-record influence when pilot instruments converge regularly, because the auxiliary noises cancel at first order (S16.2). At finite amplitude the full Hermite covariance and random instruments determine a different sandwich. A simultaneous record-count, span-growth and noise-decay rate theorem remains open.

### 7.3. Covariance calibration and the common-trajectory case

The known-covariance construction isolates the initial-state nuisance question. Estimated covariance requires a consistent calibration and its influence contribution. Separate fixed-horizon residual fits estimate a new initial state for every record; their residual quadratic forms can retain incidental-parameter bias. Thus consistency of a residual covariance estimate needs a bias correction or an independently validated calibration scheme (S16.3). The pair-difference estimator in S14 applies to equal means only.

If all records share one trajectory and covariance is known, their Gaussian path average is sufficient for the mean parameter and has covariance $C/n$. Path averaging therefore provides the natural route to Corollary 1. For nonlinear polynomial moments it also suppresses higher Hermite orders faster than averaging separate moments. The common-trajectory fixed-family GMM result is retained as Lemma 2 and a calibration diagnostic in S14–S17.

'''
    at=text.index('## 8.');text=text[:at]+section+text[at:]
    text=text.replace('Let $J\\in C^1',r'Along $x_{\theta,u}$, put $J(t)=D_xf(x_{\theta,u}(t),\theta)$ and $Z(t)=\partial_u x_{\theta,u}(t)$. Assume $J\in C^1')
    text=text.replace(r'\langle C^{-1}g,h\rangle=\langle g,h\rangle_{\mathcal H}',r'\langle C^{-1}g,v\rangle=\langle g,v\rangle_{\mathcal H}')
    text=text.replace(r'a function $h\in\mathcal H$ annihilates',r'a function $v\in\mathcal H$ annihilates').replace(r'$\mathcal D h=0$',r'$\mathcal D v=0$').replace(r'hence $h=Zc$',r'hence $v=Zc$')
    text=text.replace(r'has Gaussian variance $V_\phi=\iint\ell(u)C(u,v)\ell(v)dudv$.',r'has Gaussian variance'+ '\n\n$$\n'+r'V_\phi=\iint\ell(u)C(u,v)\ell(v)dudv.'+'\n$$\n\n')
    # Explicit sensitivity numbers replace a rounded tautology.
    start=text.index('A separate 128-path Van der Pol check');end=text.index('\n\n',start)
    text=text[:start]+r'''A separate 128-path Van der Pol check changes temporal weights to $\tau/2$ and $2\tau$, keeping the correction variance known. At 0.5% noise the RMSEs are 0.001747 and 0.001794, versus 0.001737 with the correct shape; at 5% they are 0.026242 and 0.026158, versus 0.025583. A stationary OU residual-profile fit using one known-initial scoring step gives RMSE 0.0017366 and 0.025183 at the two amplitudes. This check assesses a feasible residual plug-in; all fits converge (S12). Fixed-horizon infill identifies the microergodic combination $\sigma_j^2/\tau$, while separate variance and range consistency requires additional acquisition or calibration [@ying1991; @zhang2004].'''+text[end:]
    # Replace repeated defensive prose by direct scope statements.
    substitutions={
      'This is a local information statement, not an identity between finite-noise estimators.':'The equivalence concerns local information.',
      'Numerical resolution, feature bias, residual covariance and model selection interact, so reducing quadrature error does not by itself establish statistical consistency.':'Statistical consistency depends jointly on feature correction, residual covariance and the acquisition limit.',
      'We do not propose replacing these efficient full-observation methods.':'These methods provide efficient full-observation references.',
      'computational superiority of weak estimation is not inferred from them.':'the comparison concerns information and accuracy.',
      'it is not identical to a frozen-weight GLS update.':'its finite-noise objective differs from frozen-weight GLS.',
      'We do not claim weak GLS or optimized test selection themselves as new methods.':'Our contribution quantifies the attainable information of the complete weak class.',
      'does not replace these broader coverage studies.':'complements these broader coverage studies.',
      'it does not remove random design–response coupling.':'random design–response coupling remains in coefficient estimates.',
      "This efficiency concerns the Gaussian linearization; it is not the unguarded ratio's variance.":'This efficiency concerns the Gaussian linearization.',
      'rather than a pathwise pairing for every distribution.':'through Gaussian observations.',
      'It does not assert density of the finite multiscale bump families used in the experiments.':'The experiments use finite multiscale approximations.',
      'an analytical identity, not a guarantee of stable numerical inversion of $Z(t)$;':'an analytical identity whose numerical inversion depends on conditioning;',
      'rather than recovery of the known-initial score.':'under the unknown-initial information reference.',
      'but the released iid covariance is not the correct physical-OU covariance.':'and physical-OU covariance changes its temporal weights.',
      'it does not supply global optimization or covariance-calibration consistency for a released implementation.':'its assumptions concern a regular local optimizer and consistent covariance calibration.',
      'that are not supplied here.':'requiring a separate rate analysis.',
      'rather than an algorithm name or finite-noise support-selection procedure.':'with supplied identifiable structure.',
      'they do not alone remove coefficient bias.':'coefficient bias also depends on random-design terms.',
      'This population moment identity, rather than an unbiased single-record ratio, permits replication-based inference.':'This population moment identity permits replication-based inference.',
      'A clean-moment denominator diagnostic does not filter paths and is not used by the estimator.':'The clean-moment denominator is reported solely as a diagnostic; every path is retained.',
      'not a finite-noise bound on the nonlinear bounded estimator.':'for local information comparison.',
      'agreement of medians alone does not establish covariance efficiency.':'covariance efficiency is assessed separately by linearized variance.',
      'an approximation experiment, not a proof of density for an infinite sequence.':'a finite-span approximation experiment.',
      'it does not imply that every iid-weighted dense family must have a positive limiting gap.':'its magnitude depends on the span and weighting rule.',
      'The true trajectory is not supplied to this estimator.':'The covariance uses the fitted trajectory.',
      'these finite-cohort summaries do not assert finite unconditional moments.':'population moment claims require the risk conditions stated in Corollary 1.',
      'Intervals are not adjusted for multiple comparisons.':'Intervals are pointwise across comparisons.',
      'a likelihood refinement, not a computed global MLE.':'a local likelihood refinement.',
      'rather than evidence of a finite-noise bound violation.':'within the local reference uncertainty.',
      'as a follow-up, not as an independent replication or an unchanged one-step protocol.':'as a likelihood-refinement follow-up on the same observations.',
      'The theorem concerns the limiting influence and does not ensure efficiency of every finite-step implementation.':'The theorem concerns the limiting influence; the scoring implementation also needs a sufficiently accurate local initialization.',
      'contains zero and does not resolve an ordering.':'contains zero, leaving the ordering unresolved.',
      'Exact residual covariance does not cancel these terms.':'Their effects remain under exact residual covariance.',
      'This checks three-state polynomial coupling, not arbitrary Lorenz library selection or long-horizon chaotic inference.':'This checks three-state polynomial coupling with supplied structure over the stated horizon.',
      'AI tools are not credited as authors.':'Authorship remains with the named human author.',
      'Original author sources are downloaded from their pinned repositories rather than redistributed.':'Pinned download specifications retrieve original author sources.',
    }
    for before,after in substitutions.items():text=text.replace(before,after)
    start=text.index('## 9.');end=text.index('## Declaration',start)
    text=text[:start]+'''## 9. Discussion and limitations

The central result is a statistical reference for boundary-free weak equations when initial states are nuisance parameters. Their exact adjoint range retains the nuisance-adjusted local information. Known initial states make the excluded directions useful: the local RMSE cost is 29% for Van der Pol and almost a factor two for Lorenz. Rich test spans recover most available weak information; correct temporal weights add a smaller improvement. The full-path likelihood reference quantifies the value of observations excluded by the boundary-free class.

Different initial states link the product information geometry to conditional weak identities. Model covariance computes within-record weights; independent instruments preserve those identities at fixed noise, and the low-dimensional record score supplies its sandwich. A common trajectory instead favors averaging the Gaussian paths. The six-test common-trajectory experiment is consequently retained as a calibration diagnostic. Its bootstrap at 50 records overcovers and greatly expands the confidence region, motivating the richer model-weighted heterogeneous-record experiment.

Local efficiency requires regular pilots, controlled quadrature and compatible span enrichment. Moderate-noise performance depends additionally on random-design and nonlinear terms; the matched WENDy comparison leaves the accuracy ordering unresolved. Joint record/span/noise rates, consistent covariance calibration with heterogeneous fixed-horizon nuisance states, partial observation, process noise and support selection need further analysis. The results concern fully observed deterministic ODEs with supplied identifiable structure, and provide theory and a reproducible benchmark for evaluating weak estimators.

'''+text[end:]
    if (ROOT/'reports/record_bootstrap_results.json').exists():
        data=json.loads((ROOT/'reports/record_bootstrap_results.json').read_text());s=r'''## S17. Bootstrap calibration of the common-trajectory diagnostic

The 1000 original repetitions and all five record counts are reused exactly; SHA256 equality checks the regenerated sufficient statistics. The bootstrap protocol was fixed before its outputs: 999 resamples of independent record pairs, recentered moments, both GMM steps and correction variances reestimated, with marginal bootstrap-t quantiles and a studentized joint quadratic statistic [@hall1996]. Pair resampling preserves the covariance estimator's independent clusters. All 4,995,000 fits are retained. Batched fits agree with the scalar estimator within 1.3e-12. This is a supporting small-family diagnostic, conditional on common trajectories.

Table S16. Common-trajectory estimated-variance GMM: nominal 95% coverage before and after recentered pair bootstrap. Width ratios are medians across repetitions; the last column is the median four-dimensional region-volume ratio. Monte Carlo coverage uncertainty and individual marginal results are archived.

| n | Wald joint | Bootstrap joint | Bootstrap marginal range | Median width ratio range | Median volume ratio |
|---:|---:|---:|---:|---:|---:|
'''
        for g in data['groups']:
            s+=f"| {g['n']} | {100*g['wald_joint_coverage']:.1f}% | {100*g['joint_coverage']:.1f}% | {100*min(g['marginal_coverage']):.1f}–{100*max(g['marginal_coverage']):.1f}% | {min(g['median_marginal_width_ratio']):.2f}–{max(g['median_marginal_width_ratio']):.2f} | {g['median_joint_volume_ratio']:.2f} |\n"
        s+='\nAt 50 records the joint bootstrap region covers 99.9% and has median volume about 248 times the Wald region: the correction overcompensates for instability of estimating a 12-dimensional weight covariance from 25 independent pairs. Marginal coverage is 97.4–99.0%. At 200 records joint coverage is 94.7%, and at 800 it is 95.2%. Hence the bootstrap improves calibration at larger counts but the smallest-count construction remains excessively conservative.\n\n![Figure S9. Common-trajectory bootstrap diagnostic: coverage with Wilson Monte Carlo intervals and median marginal interval-width inflation. All methods reuse the identical 1000 repetitions.](figures/publication/record_bootstrap.png)\n\n'
        supp=supp.replace('## References',s+'## References')
    archive=json.loads((ROOT/'configs/public_archive.json').read_text())
    if archive.get('doi'):
        text=text.replace('a DOI has not been assigned.',
            'the code and data are preserved in the [Zenodo archive](https://doi.org/'+archive['doi']+').')
    p.write_text(text,encoding='utf-8');(BASE/'supplement.md').write_text(supp,encoding='utf-8')
    cover=BASE/'submission/cover_letter_draft.md'
    cover.write_text(f'''Dear Editors of the Electronic Journal of Statistics,

Please consider the manuscript “{TITLE}” by Danil Romanovich Malyshev.

The manuscript gives an exact adjoint characterization of boundary-free weak ODE influences under correlated Gaussian measurement error. Its sharp local information is the initial-state Schur complement. An explicit OU optimal generalized test makes the information cost computable, while exact polynomial residual covariance supplies model weights. A heterogeneous-record extension preserves conditional weak identities through independent instruments and connects them to the product nuisance-information reference. The contribution is statistical theory and a reproducible benchmark; established Gaussian projection, BLUE, GMM and data-fission methods are explicitly credited.

Scalar, Van der Pol and Lorenz simulations examine span approximation, temporal weighting, finite-noise accuracy and record-level inference. Complete adaptive benchmarks, the original common-trajectory calibration experiment and its conservative small-sample bootstrap remain in the supplement. The main manuscript is prepared in the IMS EJS class, with MSC 62F12, 62M10 and 34A55, and a versioned code/data archive.

Corresponding author: Danil Romanovich Malyshev, Independent researcher, Almaty, Kazakhstan; malyshevdanil04@gmail.com. ORCID: https://orcid.org/0009-0007-0947-3866.

Author approval, funding, conflict-of-interest and exclusive-submission declarations require factual completion before submission; this draft makes no assertion about them.
''',encoding='utf-8')
    print(json.dumps(dict(nuisance_frame=True,independent_record_section=7,experimental_section=8,target_journal='Electronic Journal of Statistics',orcid_added=True)))


if __name__=='__main__':main()
