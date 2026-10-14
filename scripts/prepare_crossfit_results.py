"""Use completed cross-fitting evidence; preserve external-pilot results in S18."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def main():
    report=json.loads((ROOT/'reports/crossfit_records_results.json').read_text());assert report['completed']
    validation=report.get('validation') or json.loads((ROOT/'reports/crossfit_run_validation.json').read_text());assert validation['full_inferential_runs']
    path=BASE/'article_revised.md';text=path.read_text(encoding='utf-8');p=report['protocol'];groups=report['groups'];geometry=report['geometry']
    if '### 8.6. Cross-fitted' in text:raise RuntimeError('Rebuild preceding manuscript stages before applying results')
    section=r'''### 8.6. Cross-fitted heterogeneous-record inference

Each of 400 repetitions contains 1050 independent Van der Pol records with different initial states, uniform on $[1.75,2.25]\times[-0.25,0.25]$. The first 800 records match the earlier conditional-pilot benchmark by input-array SHA256; the remaining 250 have a separate seed. Prefixes $N=50,100,200,400,800,1050$ use five fixed balanced folds. Each complement supplies an averaged six-test corrected-moment coefficient pilot at the same noise amplitude as its held-out fold. All $N$ records enter the final 480-window independent-view estimating equation. External training uses zero records. Noise covariance is known, with $T=20$, $h=0.01$ and $\tau=0.3$.

Full Proposition 3 covariance supplies model instruments. A data-independent, bounded degree-five polynomial cache approximates normalized model instruments on a coefficient cube of radius 0.06 around the declared parameter-box midpoint. Its 512 Sobol design points use a separate fixed seed. Fixed anchor-information factors restore the relative precision of the initial-state grid points when pooling. Coefficient-map coordinates are clipped to the cache cube; this changes weights while preserving the independent-view mean identity. Initial-state interpolation uses each instrument view's noisy first observation. Generating initial states enter only the information benchmark. This is a finite instrument approximation, with separate exact-model checks before Monte Carlo inference.

'''
    section+=f"The maximum checked trace-variance loss against exact within-record GLS is {100*geometry['maximum_variance_loss']:.3f}%, below the prespecified 5% threshold. This check covers the declared coefficient offsets and initial-state probes; it supplies neither a uniform global approximation theorem nor exact product-information efficiency for the finite cache.\n\n"
    if p.get('cache_protocol_file')=='crossfit_cache_v3.json':
        section=section.replace('degree-five','degree-seven').replace('512 Sobol','1024 Sobol')
    section+=r'''For each repetition and count, 999 bootstrap draws resample complete records within their original folds and refit every complement pilot. Copies of a physical record retain their original fold. A frozen-instrument comparator uses those identical stratified bootstrap counts, isolating the effect of pilot refitting on the critical values. Both use studentized marginal pivots and a joint quadratic pivot. All planned repetitions and bootstrap fits are retained.

Table 6. Five-fold heterogeneous-record inference at 5% noise, 400 repetitions. RMSE and its local reference are multiplied by $\sqrt N$. Bias is unscaled relative coefficient bias. Coverage refers to a joint 95% region; the bootstrap refits complement pilots.

| N | Scaled RMSE | Local reference | Bias | Wald joint | Bootstrap joint |
|---:|---:|---:|---:|---:|---:|
'''
    for r in groups:
        if r['level']==.05 and r['method']=='crossfit480':
            section+=f"| {r['n']} | {r['sqrt_n_rmse']:.5f} | {r['sqrt_n_local_reference']:.5f} | {r['relative_bias']:.6f} | {100*r['wald_joint_coverage']:.1f}% | {100*r['refit_joint_coverage']:.1f}% |\n"
    section+='\n'
    for level in p['levels']:
        first=next(r for r in groups if r['level']==level and r['method']=='crossfit480' and r['n']==50)
        last=next(r for r in groups if r['level']==level and r['method']=='crossfit480' and r['n']==1050);lo,hi=last['refit_joint_wilson_95']
        section+=f"At {100*level:g}% noise, scaled RMSE is {first['sqrt_n_rmse']:.6f} for N=50 and {last['sqrt_n_rmse']:.6f} for N=1050, against the latter local reference {last['sqrt_n_local_reference']:.6f}. At N=1050, joint coverage is {100*last['refit_joint_coverage']:.1f}% with pilot refitting (Wilson interval [{100*lo:.1f}%, {100*hi:.1f}%]) and {100*last['frozen_joint_coverage']:.1f}% with frozen instruments. The paired refit-minus-frozen coverage difference is {100*last['refit_minus_frozen_joint_coverage']:.1f} percentage points, with 95% interval [{100*last['refit_minus_frozen_paired_95'][0]:.1f}, {100*last['refit_minus_frozen_paired_95'][1]:.1f}].\n\n"
        deficient=[r['n'] for r in groups if r['level']==level and r['method']=='crossfit480' and r['refit_joint_wilson_95'][1]<.95]
        if deficient:
            section+=f"Joint undercoverage persists at {100*level:g}% for N={', '.join(map(str,deficient))}: the Wilson intervals lie below the nominal 95% level. Pilot refitting alone leaves this calibration question open; the finite instrument approximation and quadrature warrant separate sensitivity checks.\n\n"
        else:
            section+=f"At {100*level:g}%, no checked count has a joint-coverage Wilson interval entirely below 95%. These finite Monte Carlo intervals support a calibration assessment rather than an exact-coverage claim.\n\n"
    for r in report['paired_budget']:
        section+=f"For the same 1050-record budget at {100*r['level']:g}% noise, the split comparator trains on the last 1000 records and estimates from the first 50. Split minus cross-fitted relative RMSE is {r['split_minus_crossfit_rmse']:.6f}, with paired 95% interval [{r['paired_95'][0]:.6f}, {r['paired_95'][1]:.6f}].\n\n"
    section+=r'''![Figure 5. Heterogeneous-record cross-fitting with all records used for inference. Left: scaled RMSE and the summed nuisance-information reference, with a matched 1050-record split-budget comparator. Right: joint 95% bootstrap coverage, comparing refitted and frozen pilots under identical within-fold resampling. Error bars are Monte Carlo RMSE bootstrap intervals and binomial Wilson intervals. Fixed covariance, finite instrument approximation and quadrature apply throughout.](figures/publication/crossfit_gmm.png)

Cross-fitting changes the data budget and pilot construction. The previous independent, frozen 1000-record pilot benchmark remains in S18. It also differs in pilot amplitude and numerical instrument representation, so comparisons with that experiment measure the combined protocol change. The paired bootstrap comparison above isolates pilot refitting on each new point estimator. Finite-noise dispersion, instrument approximation and numerical quadrature remain separate from the local information reference; Corollary 2 supplies an asymptotic justification under its stated conditions. Complete bias, dispersion, marginal coverage, interval widths and cache-clipping counts are archived in S20.

'''
    start=text.index('### 8.6.');end=text.index('## 9.',start);text=text[:start]+section+text[end:]
    text=text.replace('The tagged snapshot identifies the version supporting this article;',
        'The tagged snapshot preserves the preceding benchmark; cross-fitting artifacts require deposition in a new versioned archive before submission;')
    path.write_text(text,encoding='utf-8')
    path=BASE/'supplement.md';supp=path.read_text(encoding='utf-8');assert '## S20.' not in supp
    details='''## S20. Cross-fitting protocol and complete calibration

The frozen crossfit_records.json and crossfit_cache.json protocols specify folds, coefficient-map support, Monte Carlo seeds, bootstrap seeds and thresholds. The execution manifests preserve source fingerprints and the instrument-cache SHA256. Source and protocol snapshots, observation identity, fixture estimates, sandwich covariances and every coverage indicator are validated. Known OU covariance and fixed coefficient structure are supplied. The coefficient pilot uses only original records in the complement of its held-out fold. Resampled copies stay in their original folds. Refitted and frozen bootstrap comparisons use identical counts and random seeds.

Table S18. Complete cross-fitted calibration, 400 repetitions per condition. Joint bootstrap coverage uses refitted pilots. Parentheses contain Wilson intervals; individual marginal coverage, width ratios, bias vectors, covariance matrices and pilot/cache clipping diagnostics are available in the machine-readable archive.

| Noise | N | RMSE | Bias | Dispersion | Wald joint | Refitted bootstrap joint (Wilson) | Frozen bootstrap joint |
|---:|---:|---:|---:|---:|---:|---|---:|
'''
    for r in groups:
        if r['method']!='crossfit480':continue
        lo,hi=r['refit_joint_wilson_95']
        details+=f"| {100*r['level']:g}% | {r['n']} | {r['relative_rmse']:.6f} | {r['relative_bias']:.6f} | {r['relative_dispersion']:.6f} | {100*r['wald_joint_coverage']:.1f}% | {100*r['refit_joint_coverage']:.1f}% ({100*lo:.1f}–{100*hi:.1f}%) | {100*r['frozen_joint_coverage']:.1f}% |\n"
    details+='\nThe pre-inference affine-cache failures are retained as computational diagnostics and contribute no inferential observations. The validated polynomial map remains a numerical approximation on a bounded cube; out-of-cube pilot coordinates are clipped only for evaluating that map. Centered continuous-parameter inference additionally requires negligible instrumented quadrature defect (A5). No Monte Carlo calibration of critical values or filtering of paths is used.\n\n'
    if p.get('cache_protocol_file')=='crossfit_cache_v3.json':
        details=details.replace('crossfit_records.json and crossfit_cache.json','repair480_records.json and crossfit_cache_v3.json')
    supp=supp.replace('## References',details+'## References');path.write_text(supp,encoding='utf-8')
    print('Completed cross-fitting evidence inserted')


if __name__=='__main__':main()
