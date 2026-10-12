"""Insert completed heterogeneous-record evidence into the nuisance-framed paper."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def main():
    data=json.loads((ROOT/'reports/heterogeneous_records_results.json').read_text());groups=data['groups'];text=(BASE/'article_revised.md').read_text(encoding='utf-8');supp=(BASE/'supplement.md').read_text(encoding='utf-8')
    if '### 8.6.' in text:raise RuntimeError('Rebuild nuisance manuscript before adding results')
    g=data['geometry'];p=data['protocol'];low=[r for r in groups if r['level']==.005 and r['method']=='fission480'];high=[r for r in groups if r['level']==.05 and r['method']=='fission480']
    g={**g,'training_coefficients':[round(v,5) for v in g['training_coefficients']]}
    section=r'''### 8.6. Heterogeneous records with rich model-weighted tests

Four hundred independent repetitions each contain 800 Van der Pol records, with initial states uniform on $[1.75,2.25]\times[-0.25,0.25]$, matching the acquisition range in S9.4. Counts 50, 100, 200, 400 and 800 are prefixes; 0.5% and 5% amplitudes share standard OU realizations. Marginal sensor variances use the fixed baseline trajectory calibration and $\tau=0.3$, with $T=20$ and $h=0.01$. Known covariance isolates the initial-state nuisance question. The initial states are hidden from feasible instruments.

A separate 1000-record training sample supplies one frozen coefficient pilot. Model trajectories at 25 initial-state grid points compute $D^\top\Omega^{-1}$ from the full Proposition 3 covariance. Bilinear interpolation at each fission view's noisy first observation supplies the initial-state instrument pilot; the weak moments use the independent opposite view with doubled correction variance, and the two directions are averaged. Thus finite-noise validity follows from instrument separation even with imperfect trajectory pilots. Inference is conditional on the frozen training sample: n counts inference records only, and the plotted information reference concerns those records rather than the total training-plus-inference budget.

The 120- and 480-window families are those in Section 8.2. Within-record weights are model based, while a four-dimensional empirical record-score covariance supplies the sandwich. A prespecified studentized bootstrap resamples complete record equations 999 times, with marginal bootstrap-t intervals and a joint quadratic critical value. An oracle 120-window reference receives generating coefficients and clean initial states; it is distinguished from both feasible constructions. No stacking comparator enters this experiment.

'''
    section+=f"The training pilot is {g['training_coefficients']}, with relative error {g['training_relative_error']:.5f}. The 9-by-9 sensitivity-information grid has maximum relative interpolation discrepancy {max(g['information_interpolation_relative_errors']):.4f} on 12 separately seeded checks. The record-specific information matrices enter the summed Schur-complement reference; they are used only for benchmarking.\n\n"
    section+=r'''Table 6. Heterogeneous initial states at 5% noise, 400 repetitions per count. All RMSE columns are scaled by $\sqrt n$. The local reference uses the sum of initial-state-adjusted information; finite-noise estimates also have nonlinear and random-instrument terms. Joint coverage is for nominal 95% studentized record-bootstrap regions. Monte Carlo Wilson intervals and marginal coverage are in S18.

| n | Oracle 120 | Feasible 120 | Feasible 480 | Local reference | Joint 120 | Joint 480 |
|---:|---:|---:|---:|---:|---:|---:|
'''
    for n in p['record_counts']:
        o=next(r for r in groups if r['n']==n and r['level']==.05 and r['method']=='oracle120');a=next(r for r in groups if r['n']==n and r['level']==.05 and r['method']=='fission120');b=next(r for r in groups if r['n']==n and r['level']==.05 and r['method']=='fission480')
        section+=f"| {n} | {o['sqrt_n_rmse']:.5f} | {a['sqrt_n_rmse']:.5f} | {b['sqrt_n_rmse']:.5f} | {b['sqrt_n_local_information_relative_sd']:.5f} | {100*a['bootstrap_joint_coverage']:.1f}% | {100*b['bootstrap_joint_coverage']:.1f}% |\n"
    section+=f"\nFor 480 tests at 0.5% noise, the scaled RMSE is {low[0]['sqrt_n_rmse']:.6f} at 50 records and {low[-1]['sqrt_n_rmse']:.6f} at 800, compared with local references {low[0]['sqrt_n_local_information_relative_sd']:.6f} and {low[-1]['sqrt_n_local_information_relative_sd']:.6f}. At 5% noise the corresponding scaled errors are {high[0]['sqrt_n_rmse']:.5f} and {high[-1]['sqrt_n_rmse']:.5f}; the finite-amplitude comparison is reported separately from local efficiency. At 50 records the 480-test joint bootstrap coverage is {100*high[0]['bootstrap_joint_coverage']:.1f}% with Wilson interval [{100*high[0]['bootstrap_joint_wilson_95'][0]:.1f}%, {100*high[0]['bootstrap_joint_wilson_95'][1]:.1f}%]. The model-weighted construction avoids estimating a moment-dimensional empirical weight covariance from those 50 records.\n\n"
    section+='![Figure 5. Different initial states: square-root-record-count-scaled RMSE against the summed nuisance-information reference, and nominal 95% joint record-bootstrap coverage. The oracle receives clean initial states; feasible instruments use independent views and a frozen training pilot. Left-panel whiskers resample independent Monte Carlo repetitions; right-panel whiskers are binomial Wilson intervals. The local Gaussian reference is a small-noise covariance scale. Training data are additional and fixed across repetitions; the reference counts inference records only.](figures/publication/heterogeneous_gmm.png)\n\n'
    section+='The feasible 480-test construction remains about 7% above the local RMSE reference at 0.5% noise and 28% above it at 5% for 800 records. The near-reference oracle uses clean initial states and generating coefficients; its performance separates span approximation from feasible random-instrument accuracy. The coefficient training pilot is fixed at 5% noise, so this experiment also retains a pilot discrepancy as the inference amplitude decreases. It checks finite-noise consistency and accuracy of the declared construction, rather than the regular-pilot small-noise efficiency limit. At 0.5% noise and 800 records, joint bootstrap coverage is 92.25%, and its Wilson interval excludes 95%. Consequently small-sample calibration at 5% improves substantially, while uniform calibration over noise amplitudes and record counts remains unresolved.\n\n'
    at=text.index('## 9.');text=text[:at]+section+text[at:]
    details=r'''## S18. Heterogeneous-record model weights and complete interval summaries

`heterogeneous_records.json` fixes 400 repetitions, four modulo-partitioned shards, seeds 3000000–3000399, independent coefficient-training seed 3100000, and record-bootstrap seed construction before outputs. Every repetition has 800 initial states uniform on the same square as S9.4. DOP853 integrates all generating trajectories jointly at tolerances 1e-10 and 1e-12. Exact stationary AR(1) recursion generates the original OU errors and independent auxiliary fission errors. The two amplitude conditions share initial states and standard noise arrays; all methods at an amplitude share observations. SHA256 hashes identify those input arrays.

The training sample has 1000 heterogeneous records at 5% noise and known marginal variance. Its averaged six-test moment system supplies the coefficient pilot, clipped to the declared fixed parameter box before model integration. It is generated independently and then frozen across all 400 repetitions. This experiment conditions on that pilot; it measures the accuracy of the inference-record equations and does not measure total-data efficiency of a training procedure. The instrument grid has initial positions [1.5,1.75,2,2.25,2.5] and velocities [-0.5,-0.25,0,0.25,0.5]. Model trajectories receive the training coefficients, not the generating coefficients. Bilinear interpolation clips noisy initial guesses to that grid; this changes instrument accuracy while preserving the independent-view mean identity. Full covariance uses twice the original marginal variances for feasible fission views, and the original covariance for the explicitly labeled oracle.

The basis-relative cutoff is 1e-10 and covariance relaxation is 1e-10 times its mean diagonal. All model weights contain every Hermite order. Local Fisher information is computed on a 9-by-9 grid over the actual acquisition square, checked on 12 independent points, and interpolated at each clean generating initial state for the benchmark only. Every reference uses the sum of coefficient Schur complements, not the Schur complement of an averaged trajectory. The exact boundary-free reference assumes arbitrary rich influences; the numerical span and instrument grid remain finite approximations.

The four-dimensional record equation is linear in coefficients. For each count its mean derivative is inverted once. Residual score covariance supplies the sandwich, with n independent clusters and no empirical high-dimensional moment weighting. The bootstrap resamples whole record equations, conditional on the independent training pilot and the original instrument construction. It recomputes the estimate and the four-dimensional sandwich for all 999 resamples. Exact identification of this reduced score makes its empirical moment mean zero at the original estimate, so recentering changes nothing beyond rounding. Studentized marginal and joint quantiles use the frozen confidence level 0.95. No resample or Monte Carlo repetition is discarded.

Table S17. Complete feasible heterogeneous-record RMSE and coverage. Coverage uses 400 repetitions per row; the parenthesized joint interval is the binomial Wilson interval. The marginal range summarizes four individual 95% intervals; individual results and interval-width ratios are archived.

| Noise | n | Tests | RMSE | Bias | Wald joint | Bootstrap joint (Wilson interval) | Bootstrap marginal range |
|---:|---:|---:|---:|---:|---:|---|---|
'''
    for r in groups:
        if not r['method'].startswith('fission'):continue
        lo,hi=r['bootstrap_joint_wilson_95'];m=r['bootstrap_marginal_coverage']
        details+=f"| {100*r['level']:g}% | {r['n']} | {r['method'][7:]} | {r['relative_rmse']:.6f} | {r['relative_bias']:.6f} | {100*r['wald_joint_coverage']:.1f}% | {100*r['bootstrap_joint_coverage']:.1f}% ({100*lo:.1f}–{100*hi:.1f}%) | {100*min(m):.1f}–{100*max(m):.1f}% |\n"
    details+='\nKnown OU covariance is used throughout this benchmark. The conditional moment theorem permits different initial states; it does not establish consistency of a covariance estimate formed from incidental-state residual fits. Estimated-shape inference requires the additional calibration conditions in S16.3.\n\n'
    supp=supp.replace('## References',details+'## References');(BASE/'article_revised.md').write_text(text,encoding='utf-8');(BASE/'supplement.md').write_text(supp,encoding='utf-8')
    print('Heterogeneous main table, figure and supplementary calibration added')


if __name__=='__main__':main()
