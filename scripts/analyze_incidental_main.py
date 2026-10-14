"""Analyze experiment A only after all four methods have verified matched inputs."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.validate_nuisance_runs import snapshot
from scripts.validate_crossfit_runs import validate
from scripts.analyze_independent_records import wilson
ROOT=Path(__file__).resolve().parents[1]


def main():
    compact=(ROOT/'configs/incidental_compact.json').exists();amended=(ROOT/'configs/incidental_trimmed.json').exists()
    p=json.loads((ROOT/'configs'/('incidental_compact.json' if compact else 'incidental_trimmed.json' if amended else 'incidental_main.json')).read_text());R=p['monte_carlo_repetitions']
    state=json.loads((ROOT/'reports/incidental_main_controller.json').read_text());assert state['likelihood_components_complete']
    families=p.get('weak_families',[f'repairA{i}' for i in range(len(p['levels']))] if (ROOT/'configs/repairA0_records.json').exists() else ['crossfitA'])
    for family in families:validate([f'{family}_records_{i}' for i in range(4)],full=True)
    rows=[];hashes={}
    for stage in state['downloaded']:
        if stage=='inc_main_smoke':continue
        snapshot(stage,'incidental_compact.json' if stage.startswith('inc_short_') else 'incidental_trimmed.json' if stage.startswith('inc_trim_') else 'incidental_main.json');folder=ROOT/f'results/kaggle_{stage}'
        rows+=json.loads((folder/'rows.json').read_text())
        for h in json.loads((folder/'input_hashes.json').read_text()):
            assert hashes.setdefault(h['repetition'],h['input_sha256'])==h['input_sha256']
        assert all(r['diagnostic']['success'] for r in json.loads((folder/'bootstrap_fits.json').read_text()))
    for family in families:
        for i in range(4):
            folder=ROOT/f'results/kaggle_{family}_records_{i}'
            for h in json.loads((folder/'input_hashes.json').read_text()):assert h['first_800_input_sha256']==hashes[h['repetition']]
            for r in json.loads((folder/'rows.json').read_text()):
                if r['method']=='crossfit'+str(p.get('weak_windows',480)):rows.append({**r,'method':'crossfit_fission_gmm'})
    methods=p['methods'] if compact else ['profile_likelihood','bias_corrected_likelihood','crossfit_fission_gmm','oracle_likelihood']
    truth=np.array(p['truth']);rng=np.random.default_rng(4600000);groups=[];paired=[]
    curvature=json.loads((ROOT/'reports/incidental_curvature_validation.json').read_text())['results']['results'][-1]
    predicted_bias=np.array(curvature['bias_coefficient'])
    bias_checks=[]
    for level in p['levels']:
        for n in p['record_counts']:
            errors={}
            cell_methods=[m for m in methods if m!='bias_corrected_likelihood' or not amended or any(c['level']==level and c['n']==n for c in p['bootstrap_correction_cells'])]
            for method in cell_methods:
                g=sorted([r for r in rows if r['method']==method and r['level']==level and r['n']==n],key=lambda r:r['repetition'])
                assert [r['repetition'] for r in g]==list(range(R))
                e=(np.array([r['coefficients'] for r in g])-truth)/2;bias=e.mean(axis=0);mse=np.sum(e*e,axis=1);errors[method]=mse
                ids=rng.integers(0,R,(5000,R));prefix='refit' if method=='crossfit_fission_gmm' else 'wald'
                count=sum(r[prefix+'_joint_covered'] for r in g);marginal=np.array([r[prefix+'_marginal_covered'] for r in g])
                group=dict(level=level,n=n,method=method,repetitions=R,relative_rmse=float(np.sqrt(mse.mean())),sqrt_n_rmse=float(np.sqrt(n*mse.mean())),
                    relative_bias=float(np.linalg.norm(bias)),coefficient_bias=(2*bias).tolist(),
                    coefficient_bias_mc_95=np.quantile(2*e[ids].mean(axis=1),[.025,.975],axis=0).tolist(),
                    relative_dispersion=float(np.sqrt(np.mean(np.sum((e-bias)**2,axis=1)))),
                    sqrt_n_rmse_mc_95=np.quantile(np.sqrt(n*mse[ids].mean(axis=1)),[.025,.975]).tolist(),
                    sqrt_n_local_reference=float(np.sqrt(n)*np.mean([r['local_information_relative_sd'] for r in g])),
                    joint_coverage=count/R,joint_wilson_95=wilson(count,R),marginal_coverage=marginal.mean(axis=0).tolist(),
                    marginal_wilson_95=[wilson(int(k),R) for k in marginal.sum(axis=0)],interval_method=prefix)
                if method=='crossfit_fission_gmm':
                    group.update(median_estimate_seconds=float(np.median([r['estimate_seconds'] for r in g])),
                        median_bootstrap_seconds=float(np.median([r['bootstrap_refit_seconds'] for r in g])),
                        median_shared_record_preprocessing_seconds=float(np.median([r['shared_preprocessing_seconds'] for r in g])),
                        shared_preprocessing_records=g[0]['shared_preprocessing_records'])
                    group['median_online_seconds_including_all_shared_preprocessing']=float(np.median([r['estimate_seconds']+r['shared_preprocessing_seconds'] for r in g]))
                    group['timing_caveat']='Shared preprocessing handles all generated records and n prefixes, not just this prefix; no unmeasured per-record scaling is asserted. Offline instrument-cache construction is separate.'
                else:
                    group['median_runtime_seconds']=float(np.median([r['diagnostic']['elapsed_seconds'] for r in g]))
                    if all('shared_weak_start_seconds' in r for r in g):
                        group['median_online_seconds_including_weak_start']=float(np.median([r['diagnostic']['elapsed_seconds']+r['shared_weak_start_seconds'] for r in g]))
                    else:group['timing_caveat']='Legacy block records optimizer time; its weak-start time was not recorded and is not imputed.'
                if method=='profile_likelihood':
                    coefficient_errors=2*e;prediction=level**2*predicted_bias;observed=coefficient_errors.mean(axis=0)
                    Vmean=np.cov(coefficient_errors,rowvar=False,ddof=1)/R;Vi=np.linalg.inv(Vmean)
                    discrepancy=observed-prediction;stat=float(discrepancy@Vi@discrepancy)
                    centered=coefficient_errors[ids].mean(axis=1)-observed
                    null_stats=np.einsum('bi,ij,bj->b',centered,Vi,centered)
                    direction=predicted_bias/np.linalg.norm(predicted_bias)
                    projections=coefficient_errors[ids].mean(axis=1)@direction
                    bias_checks.append(dict(level=level,n=n,predicted_coefficient_bias=prediction.tolist(),observed_coefficient_bias=observed.tolist(),
                        observed_bias_norm=group['relative_bias'],predicted_relative_bias=float(np.linalg.norm(prediction)/2),
                        vector_discrepancy_statistic=stat,centered_mc_bootstrap_pvalue=float((1+np.sum(null_stats>=stat))/(len(null_stats)+1)),
                        familywise_reference_alpha=.05/(len(p['levels'])*len(p['record_counts'])),
                        predicted_signed_projection=float(prediction@direction),observed_signed_projection=float(observed@direction),
                        signed_projection_mc_95=np.quantile(projections,[.025,.975]).tolist(),
                        interpretation='Vector comparison avoids upward bias of the norm of a noisy Monte Carlo mean. Finite-n and higher-noise remainder terms must be assessed before diagnosing a missing expansion term.'))
                assert np.isclose(group['relative_rmse']**2,group['relative_bias']**2+group['relative_dispersion']**2)
                groups.append(group)
            ids=rng.integers(0,R,(5000,R))
            for other in ('profile_likelihood','bias_corrected_likelihood'):
                if other not in errors:continue
                delta=np.sqrt(errors['crossfit_fission_gmm'][ids].mean(axis=1))-np.sqrt(errors[other][ids].mean(axis=1))
                paired.append(dict(level=level,n=n,comparison='weak_minus_'+other,
                    rmse_difference=float(np.sqrt(errors['crossfit_fission_gmm'].mean())-np.sqrt(errors[other].mean())),
                    paired_mc_95=np.quantile(delta,[.025,.975]).tolist()))
    crossings=[dict(level=level,first_checked_n=min(r['n'] for r in paired if r['level']==level and r['comparison']=='weak_minus_profile_likelihood' and r['paired_mc_95'][1]<0))
               for level in p['levels'] if any(r['level']==level and r['comparison']=='weak_minus_profile_likelihood' and r['paired_mc_95'][1]<0 for r in paired)]
    offline=[]
    for family in families:
        completion=ROOT/f'results/kaggle_{family}_geometry/completion.json'
        if completion.exists():offline.append(dict(family=family,cache_construction_seconds=json.loads(completion.read_text())['elapsed_seconds']))
    report=dict(completed=True,protocol=p,groups=groups,paired_comparisons=paired,weak_better_checked_counts=crossings,bias_prediction_checks=bias_checks,offline_instrument_costs=offline,
        all_inputs_matched=True,all_planned_fits_retained=True,
        caveat='Reported counts are a discrete grid, not a continuous crossing estimate. Both MLE methods use local bounded optimization; bias-corrected Wald regions retain the original Fisher covariance. Finite-sample coefficient estimates can be biased even when weak scores are exactly centered.')
    (ROOT/'reports/incidental_main_results.json').write_text(json.dumps(report,indent=2))
    fig,axes=plt.subplots(1,len(p['levels']),figsize=(5*len(p['levels']),4),squeeze=False)
    for ax,level in zip(axes[0],p['levels']):
        for method in methods:
            g=[r for r in groups if r['level']==level and r['method']==method]
            if not g:continue
            ax.plot([r['n'] for r in g],[r['sqrt_n_rmse'] for r in g],marker='o',label=method.replace('_',' '))
        g=[r for r in groups if r['level']==level and r['method']=='profile_likelihood']
        ax.plot([r['n'] for r in g],[r['sqrt_n_local_reference'] for r in g],'k--',label='nuisance information reference')
        ax.set(xlabel='Independent records n',ylabel='sqrt(n) × relative RMSE',title=f'{100*level:g}% noise')
        ax.grid(alpha=.25)
    axes[0,0].legend(fontsize=7);fig.tight_layout();folder=ROOT/'manuscript/figures/publication';folder.mkdir(exist_ok=True)
    fig.savefig(folder/'incidental_comparison.png',dpi=220);fig.savefig(folder/'incidental_comparison.pdf');plt.close(fig)
    print('Four-method matched incidental-state comparison validated and analyzed')


if __name__=='__main__':main()
