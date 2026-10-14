"""Reframe the current manuscript around information cost and endpoint recovery."""
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def main():
    path=BASE/'article_revised.md';s=path.read_text(encoding='utf-8')
    if '## 7. Independent records' not in s:raise RuntimeError('Reframing already applied; edit the final text directly')
    independent=s[s.index('## 7. Independent records'):s.index('## 8.')]
    extension=s[s.index('### 8.6.'):s.index('### 8.8.')]
    archive=independent+'\n'+extension
    archive=archive.replace('Theorem 3','Theorem S16.1').replace('Corollary 2','Corollary S16.2')
    (BASE/'independent_record_extension_archive.md').write_text('# Independent-record extension retained outside the main narrative\n\n'+archive,encoding='utf-8')
    supppath=BASE/'supplement.md';supp=supppath.read_text(encoding='utf-8').replace('Theorem 3','Theorem S16.1').replace('Corollary 2','Corollary S16.2')
    extra='## S21. Independent-record and incidental-state numerical extensions\n\n'
    extra+=extension.replace('### 8.6.','### S21.1.').replace('### 8.7.','### S21.2.')
    supp=supp.replace('## References',extra+'\n## References');supppath.write_text(supp,encoding='utf-8')
    s=s.replace(independent,'').replace(extension,'')
    s=s.replace('Boundary-free weak estimating equations under correlated noise: efficiency with nuisance initial states',
                'The information cost of boundary-free weak estimation under correlated noise')
    a=s.index('## Abstract\n\n')+len('## Abstract\n\n');b=s.index('**Keywords:**',a)
    abstract=('Boundary-free weak equations retain the local information of a dynamical Gaussian mean experiment with unknown initial state. '
              'For deterministic ODEs observed in stationary Ornstein–Uhlenbeck error, we characterize their closed adjoint influence class and its sharp initial-state Schur-complement covariance reference. '
              'With a known initial state, the restriction raises the local RMSE reference by 29% for Van der Pol and almost doubles it for Lorenz. '
              'The extended known-initial class recovers full local information; this result quantifies established boundary adaptations. '
              'Universal time-only endpoint windows recover 97.1% of the weak-to-scoring RMSE gap on 1968 paired Van der Pol records at 0.5% noise. '
              'A 500-record ablation attributes 98.9% recovery to initial-trace tests alone; an initial-noise control alone recovers 14.3%. '
              'Exact polynomial Hermite covariance supplies model-based weights. For an added cubic coefficient, initial-state projection changes the local composite detection reference by only 0.68%, while the supplied adaptive WSINDy selects the term between two and four reference units at small noise. '
              'These are information and algorithmic benchmarks: covariance weights require a pilot model, and local detection references differ from finite-amplitude selection guarantees.\n\n')
    s=s[:a]+abstract+s[b:]
    old_start=s.index('The nuisance interpretation also connects');old_end=s.index('Sparse identification',old_start)
    s=s[:old_start]+('Known initial states also permit a remedy. Initial-trace conditions are established in orthogonality estimation [@brunel2014, Section 2.3]. '
                     'We characterize their complete local influence class, quantify finite-family recovery under physical temporal correlation, and distinguish their contribution from an initial-noise control. '
                     'The control reveals a correlated error whose covariance with interior moments can reduce variance. The ablation determines its practical contribution rather than presuming it drives recovery.\n\n')+s[old_end:]
    start=s.index('The contribution is a theory');end=s.index('## 2.',start)
    s=s[:start]+('The paper has three linked aims. Theorems 1 and 2 quantify the fixed-horizon limit and the exact information cost of boundary-free influences. '
                 'Theorem 3 describes full-information recovery through known-initial trace conditions or initial controls, separating class closure from finite-basis performance. '
                 'Proposition 3 supplies exact Hermite residual covariance for model weights and statistical testing. Simulations then distinguish boundary loss from the behavior of sparse selection. '
                 'Gaussian projection and correlated-regression BLUE supply established statistical tools [@dette2016]. Boundary adaptations are explicitly credited [@brunel2014]; the contribution is their exact OU information geometry and reproducible quantitative benchmarks. '
                 'Independent-record GMM, data fission and cross-fitting are retained as supplementary extensions (S14–S18 and S21).\n\n')+s[end:]
    pos=s.index('## 6.');s=s[:pos]+(BASE/'endpoint_class_theory.md').read_text(encoding='utf-8')+'\n\n'+s[pos:]
    start=s.index('### 6.3.');end=s.index('## 8.',start)
    s=s[:start]+('### 6.3. Temporal misspecification\n\n'
                 'For a single record, the leading small-noise covariance under weights $W$ is $P_W\\Omega_{\\mathrm{true},1}P_W^\\top$. '
                 'Weight misspecification changes efficiency; correction-variance misspecification can also change the moment mean. '
                 'In S12, assuming $\\tau/2$ or $2\\tau$ changes only weighting, whereas fitting variances changes correction as well. '
                 'The main experiments supply covariance. Calibration from paired differences requires equal mean trajectories (S14.1); the heterogeneous-record calibration question is treated separately in S16.3.\n\n')+s[end:]
    # Remove the short duplicate closure derivation from the experiment; cite the theorem.
    start=s.index('Here $f_\\theta^H$ denotes');end=s.index('We compare four branches',start)
    s=s[:start]+('Here $f_\\theta^H$ denotes coordinatewise Hermite correction under the supplied sensor variance. '
                 'Theorem 3 gives the full closed influence class when $u$ is known.\n\n')+s[end:]
    mapping={1:1,2:2,3:3,4:4,5:5,8:6,9:8}
    s=re.sub(r'^## 8\.', '## 7.', s, flags=re.M)
    s=re.sub(r'^## 9\.', '## 8.', s, flags=re.M)
    for old,new in mapping.items():
        s=s.replace(f'### 8.{old}.',f'### 7.{new}.')
        s=s.replace(f'Section 8.{old}',f'Section 7.{new}')
    ablation='''### 7.7. Endpoint-component attribution

The frozen ablation uses the same 500 records (indices 33–532). It is component attribution on already evaluated data. All new branches use identical Simpson quadrature and the same covariance-pilot rule; a matched boundary-free branch isolates quadrature from endpoint changes. The full branch and known-initial scoring are archived paired references. For a positive matched baseline-to-scoring gap, the prespecified material-contribution criterion is at least 50% gap recovery, with uncertainty reported separately.

| Added component | RMSE | Gap closed | Paired 95% interval |
|---|---:|---:|---:|
| Initial half-windows | 0.00132445 | 98.9% | [95.1%; 102.7%] |
| Terminal half-windows | 0.00164493 | −0.3% | [−1.8%; 1.2%] |
| Initial controls only | 0.00159795 | 14.3% | [7.0%; 21.4%] |
| All components | 0.00132453 | 98.9% | [95.1%; 102.6%] |
| Matched boundary-free | 0.00164395 | 0% | — |
| Known-initial scoring | 0.00132100 | 100% | — |

Initial traces alone recover essentially the entire finite-family gap. This is the mechanism already proposed by Brunel et al. [@brunel2014, Section 2.3], with its OU efficiency quantified here. The initial-noise controls provide a measurable standalone improvement, but fail the 50% material-contribution rule and add no resolved improvement over initial traces. Terminal traces alone give no resolved improvement, consistent with the unchanged closed class in Theorem 3. Component gains overlap; their percentages are not additive. The archived trapezoidal boundary-free RMSE is 0.00164343, so quadrature accounts for only 0.16% of the matched gap. All intervals use 5000 shared bootstrap resamples. Recovery percentages can exceed 100% through Monte Carlo variation; they are not clipped.

'''
    s=s.replace('### 7.8.',ablation+'### 7.8.')
    start=s.index('## 8. Discussion');end=s.index('## Declaration',start)
    discussion='''## 8. Discussion and limitations

Boundary-free equations retain exactly the local information available when the initial state is a nuisance parameter. Their restriction has a quantified cost when the state is known: 29% in the Van der Pol aggregate RMSE reference and almost a factor two for Lorenz. Rich spans and correct temporal weights approximate the sharp class reference. Theorem 3 supplies the complementary full-information class for known-initial boundary adaptations.

The recovery experiment and ablation distinguish class theory from method novelty. Universal endpoint windows recover 97.1% of the weak-to-scoring gap on 1968 records. On the 500-record ablation, initial traces alone explain essentially all that gain, reproducing the established Brunel adaptation. The correlated-noise initial control is useful separately, recovering 14.3% of the gap, but contributes no resolved improvement after initial traces. Thus the remedy is an established boundary mechanism with a quantitative OU benchmark. Covariance weights still require an integrated pilot model; computational superiority has not been demonstrated.

Information loss is direction specific. For the added cubic coefficient, initial-state projection changes the local composite detection reference by only 0.68%, despite its larger aggregate coefficient-estimation cost. At small noise the supplied MSTLS algorithm switches from zero to complete inclusion between two and four local reference units on the coarse grid. This brackets an algorithmic transition; it neither proves that a hard threshold alone causes it nor identifies a variance loss of four to sixteen. Statistical rejection and sparse selection have different objectives and require separate null calibration and power comparisons. Large finite perturbations also leave the null local Gaussian approximation.

The scope is fully observed deterministic polynomial ODEs with known positive sensor covariance. Small-noise influence efficiency, finite-noise feasible accuracy and support selection are separate claims. Exact residual covariance does not by itself give exact finite-noise estimator variance when weights and designs depend on the record. On the matched moderate-noise benchmark, likelihood is more accurate than the original boundary-free weak estimator. Independent-record GMM, fission and cross-fitting remain supplementary extensions. Joint span/noise rates, unknown-covariance calibration with heterogeneous states, partial observation and process noise require further analysis.

'''
    s=s[:start]+discussion+s[end:]
    path.write_text(s,encoding='utf-8')
    builder=ROOT/'scripts/build_ejs.py';code=builder.read_text(encoding='utf-8').replace('Weak equations with nuisance initial states','Information cost of boundary-free estimation');builder.write_text(code,encoding='utf-8')
    print('Main narrative reframed; independent-record material retained in supplement and archive')


if __name__=='__main__':main()
