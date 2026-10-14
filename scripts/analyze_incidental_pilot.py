"""Apply only the prespecified pilot escalation criterion; no main claims."""
import argparse,json
from pathlib import Path
import numpy as np
from scripts.analyze_independent_records import wilson
ROOT=Path(__file__).resolve().parents[1]


def main(refined=False):
    p=json.loads((ROOT/'configs'/('incidental_pilot_v2.json' if refined else 'incidental_pilot.json')).read_text());rows=[]
    prefix='incidental_pilot_v2' if refined else 'incidental_pilot'
    for i in range(p['shards']):rows+=json.loads((ROOT/f'results/kaggle_{prefix}_{i}/rows.json').read_text())
    groups=[];rng=np.random.default_rng(4300000)
    for method in ('profile_likelihood','oracle_likelihood'):
        g=sorted([r for r in rows if r['method']==method],key=lambda r:r['repetition']);assert [r['repetition'] for r in g]==list(range(p['repetitions']))
        assert all(r['diagnostic']['success'] for r in g)
        e=(np.array([r['coefficients'] for r in g])-p['truth'])/2;bias=e.mean(axis=0);disp=float(np.sqrt(np.mean(np.sum((e-bias)**2,axis=1))))
        sd=float(np.mean([r['relative_information_sd'] for r in g]));count=sum(r['wald_joint_covered'] for r in g)
        ids=rng.integers(0,len(g),(5000,len(g)));bias_samples=np.linalg.norm(e[ids].mean(axis=1),axis=1)
        groups.append(dict(method=method,repetitions=len(g),n=p['records'],relative_bias=float(np.linalg.norm(bias)),mean_coefficient_bias=(2*bias).tolist(),relative_dispersion=disp,
            relative_rmse=float(np.sqrt(np.mean(np.sum(e**2,axis=1)))),bias_to_mean_fisher_sd=float(np.linalg.norm(bias)/sd),bias_to_empirical_dispersion=float(np.linalg.norm(bias)/disp),
            relative_bias_mc_bootstrap_95=np.quantile(bias_samples,[.025,.975]).tolist(),wald_joint_coverage=count/len(g),wald_joint_wilson_95=wilson(count,len(g)),
            median_runtime_seconds=float(np.median([r['diagnostic']['elapsed_seconds'] for r in g])),maximum_remaining_newton_step=max(r['diagnostic']['remaining_newton_step'] for r in g),
            fits_with_active_bounds=sum(r['diagnostic']['bound_active_coefficients']+r['diagnostic']['bound_active_initial_coordinates']>0 for r in g)))
    profile=groups[0];escalate=profile['bias_to_mean_fisher_sd']<p['escalation_threshold'];levels=[.05,.1,.2] if escalate else [.05]
    report=dict(completed=True,pilot_only=True,protocol=p,groups=groups,add_10_and_20_percent=escalate,planned_main_levels=levels,criterion='coefficient bias norm / mean Fisher norm-scale standard error',
        caveat='A 50-repetition pilot selects amplitudes; it is not pooled with the independent main Monte Carlo sample and does not establish inconsistency or a crossing.')
    (ROOT/'reports/incidental_pilot_results.json').write_text(json.dumps(report,indent=2))
    main_protocol={**p,'protocol':'shared_coefficient_incidental_comparison_main_v1','levels':levels,'record_counts':[50,100,200,400,800],
        'monte_carlo_repetitions':200,'repetitions_per_job':20,'job_partition':['amplitude','record_count','repetition_block'],
        'seed_base':4000000,'parametric_bootstrap_draws':50,'parametric_bootstrap_seed_base':50000000,
        'weak_bootstrap_draws':999,'weak_bootstrap_seed_base':600000000,'pilot_decision_source':'reports/incidental_pilot_results.json','covariance_known':True,
        'methods':['profile_likelihood','bias_corrected_likelihood','crossfit_fission_gmm','oracle_likelihood'],
        'status':'frozen amplitude decision; main execution requires validated matching weak instruments and optimizer smoke'}
    target=ROOT/'configs/incidental_main.json'
    if target.exists():assert json.loads(target.read_text())==main_protocol,'Preserve frozen main design'
    else:target.write_text(json.dumps(main_protocol,indent=2))
    print(json.dumps(report))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--refined',action='store_true');main(a.parse_args().refined)
