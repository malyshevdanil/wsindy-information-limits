"""Postprocess downloaded Kaggle outputs; no synthetic experiments run locally."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream))


def bootstrap_improvement(a,b,rng):
    a,b=np.asarray(a),np.asarray(b)
    indices=rng.integers(0,len(a),size=(4000,len(a)))
    estimates=1-np.median(b[indices],axis=1)/np.maximum(np.median(a[indices],axis=1),1e-15)
    return {'relative_median_improvement':float(1-np.median(b)/max(np.median(a),1e-15)),
            'bootstrap_95':np.quantile(estimates,[.025,.975]).tolist()}


def analyze():
    pilot=ROOT/'results/kaggle_pilot'
    covariance=ROOT/'results/kaggle_covariance'
    extended=ROOT/'results/kaggle_extended'
    rng=np.random.default_rng(981611)
    result={'status':'completed_exploratory_study','replication_unit':'independent noise seed',
            'bootstrap_replicates':4000,'bootstrap_seed':981611,'artifacts':{}}
    for name,folder,expected in [('pilot',pilot,4320),('covariance',covariance,640)]:
        rows=[json.loads(line) for line in (folder/'runs.jsonl').read_text().splitlines()]
        assert len(rows)==expected and len({r['key'] for r in rows})==expected
        result[name]={'attempts':len(rows),'fit_failures':sum(r['status']!='ok' for r in rows),
                      'forecast_failures':sum(r.get('forecast_status')!='ok' for r in rows)}
        manifest=json.loads((folder/'manifest.json').read_text())
        result['artifacts'][name]={'manifest':manifest,
            'download_zip_sha256':hashlib.sha256((ROOT/f'results/kaggle_{name}_download/wsindy_{name}_outputs.zip').read_bytes()).hexdigest()}
    pilot_rows=[json.loads(line) for line in (pilot/'runs.jsonl').read_text().splitlines()]
    result['representative_cells']=[r for r in read_csv(pilot/'summary.csv')
        if r['method']=='known_debias' and r['noise_level']=='0.05' and r['tau_over_width']=='0.1']
    paired=read_csv(pilot/'paired_grid_comparisons.csv')
    # Require both last adjacent refinements to meet the original +/-5% criterion.
    cells={tuple(r[k] for k in ('system','regime','tau_over_width','noise_level','method')) for r in paired}
    result['pilot']['two_successive_equivalence_cells']=sum(all(any(
        tuple(r[k] for k in ('system','regime','tau_over_width','noise_level','method'))==cell
        and float(r['h_coarse'])==h and r['equivalence_within_5pct_on_this_pair']=='True'
        for r in paired) for h in (.02,.01)) for cell in cells)
    result['pilot']['equivalence_cells_evaluated']=len(cells)
    diagnostics=[]
    for path in sorted(pilot.glob('diagnostic*.json')):
        d=json.loads(path.read_text())
        if d['cell'][1]!='physical_ou':
            continue
        chosen=d['selection']['h']
        record={'cell':d['cell'],'method':d['method'],'selection':d['selection']}
        rows=[r for r in pilot_rows if [r['system'],r['regime'],r['tau_over_width'],r['noise_level']]==d['cell']
              and r['method']==d['method'] and r['status']=='ok']
        if chosen is not None:
            a={r['seed']:r for r in rows if r['h']==chosen}
            b={r['seed']:r for r in rows if r['h']==.005}
            seeds=sorted(set(a)&set(b))
            if seeds:
                record.update(bootstrap_improvement([a[s]['coefficient_error'] for s in seeds],
                                                     [b[s]['coefficient_error'] for s in seeds],rng))
                record['selected_support_success']=sum(a[s]['exact_support'] for s in seeds)/len(seeds)
                record['finest_support_success']=sum(b[s]['exact_support'] for s in seeds)/len(seeds)
                record['fit_time_ratio_finest_to_selected']=float(np.median([b[s]['fit_seconds'] for s in seeds])/
                                                                 np.median([a[s]['fit_seconds'] for s in seeds]))
                record['record_point_reduction_ratio']=.04/.005 if chosen==.04 else chosen/.005
        diagnostics.append(record)
    result['diagnostics']=diagnostics
    result['diagnostic_selection_counts']=dict(Counter(str(r['selection']['h']) for r in diagnostics))
    result['covariance']['summaries']=read_csv(covariance/'summary.csv')
    result['covariance']['paired_comparisons']=read_csv(covariance/'paired_grid_comparisons.csv')
    result['scalar']=json.loads((covariance/'scalar/scalar_summary.json').read_text())
    scalar_rows=read_csv(covariance/'scalar/scalar_ratio.csv')
    scalar_groups={float(r['h']) for r in scalar_rows}
    comparisons=[]
    for coarse,fine in zip(sorted(scalar_groups,reverse=True),sorted(scalar_groups,reverse=True)[1:]):
        a={int(r['seed']):float(r['absolute_error']) for r in scalar_rows if float(r['h'])==coarse}
        b={int(r['seed']):float(r['absolute_error']) for r in scalar_rows if float(r['h'])==fine}
        seeds=sorted(set(a)&set(b))
        comparisons.append({'coarse':coarse,'fine':fine,**bootstrap_improvement([a[s] for s in seeds],[b[s] for s in seeds],rng)})
    result['scalar']['paired_comparisons']=comparisons
    result['extended']={name:json.loads((extended/file).read_text()) for name,file in
        [('manifest','manifest.json'),('completion','completion.json'),('residual','residual_summary.json'),
         ('repetitions','repetition_summary.json'),('reference','reference_convergence.json')]}
    rep=json.loads((extended/'repetitions.json').read_text())
    a=[r['coefficient_error'] for r in rep if r['M']==1]
    b=[r['coefficient_error'] for r in rep if r['M']==8]
    result['extended']['repetition_improvement_1_to_8']=bootstrap_improvement(a,b,rng)
    result['artifacts']['extended']={'download_zip_sha256':hashlib.sha256(
        (ROOT/'results/kaggle_extended_download/wsindy_extended_outputs.zip').read_bytes()).hexdigest()}
    replication=ROOT/'results/kaggle_replication'
    if (replication/'completion.json').exists():
        result['replication']={name:json.loads((replication/file).read_text()) for name,file in
            [('manifest','manifest.json'),('completion','completion.json'),('summary','summary.json'),
             ('population','population_check.json')]}
        replication_rows=json.loads((replication/'rows.json').read_text())
        result['replication']['paired_1_to_32']=[]
        for method in result['replication']['manifest']['methods']:
            a=[r['coefficient_error'] for r in replication_rows if r['method']==method and r['M']==1]
            b=[r['coefficient_error'] for r in replication_rows if r['method']==method and r['M']==32]
            result['replication']['paired_1_to_32'].append({'method':method,**bootstrap_improvement(a,b,rng)})
        result['artifacts']['replication']={'download_zip_sha256':hashlib.sha256(
            (ROOT/'results/kaggle_replication_download/wsindy_replication_outputs.zip').read_bytes()).hexdigest()}
    (ROOT/'reports/empirical_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    figures=ROOT/'manuscript/figures'
    figures.mkdir(exist_ok=True)
    for system_name in ('linear','vdp'):
        fig,axes=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
        for regime,marker in [('physical_ou','o'),('iid','s'),('fixed_index','^')]:
            selected=sorted([r for r in result['representative_cells'] if r['system']==system_name and r['regime']==regime],key=lambda r:float(r['h']))
            h=[float(r['h']) for r in selected]
            axes[0].loglog(h,[float(r['coefficient_error_median_successful_fits']) for r in selected],marker+'-',label=regime)
            axes[1].semilogx(h,[float(r['support_success_fraction_failure_counted']) for r in selected],marker+'-',label=regime)
        axes[0].set_ylabel('Median relative coefficient error')
        axes[1].set_ylabel('Exact support fraction (10 seeds)')
        axes[1].set_ylim(-.02,1.02)
        for ax in axes:
            ax.set_xlabel('Sampling step h');ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle(f'{system_name}: known variance correction; noise 5%; physical OU tau/width = 0.1')
        for ext in ('png','pdf'): fig.savefig(figures/f'{system_name}_sampling.{ext}',dpi=220)
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3.8),constrained_layout=True)
    summaries=result['scalar']['summary']
    for index,label in [(0,'Median'),(1,'90th percentile'),(2,'99th percentile')]:
        ax.semilogx([s['h'] for s in summaries],[s['absolute_error_quantiles_50_90_99'][index] for s in summaries],'o-',label=label)
    ax.set(xlabel='Sampling step h',ylabel='Absolute coefficient error',title='One-window ratio: 1,000 independent OU paths')
    ax.legend();ax.grid(alpha=.2)
    for ext in ('png','pdf'): fig.savefig(figures/f'scalar_ratio.{ext}',dpi=220)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
    residual=[r for r in result['extended']['residual'] if r['method']=='known_debias']
    levels=[r['noise_level'] for r in residual]
    axes[0].plot(levels,[r['exact_relative_frobenius_error'] for r in residual],'o-',label='Exact polynomial covariance')
    axes[0].plot(levels,[r['first_relative_frobenius_error'] for r in residual],'s-',label='First-order covariance')
    axes[0].set_ylabel('Relative empirical covariance discrepancy')
    axes[0].legend(fontsize=8)
    axes[1].plot(levels,[r['higher_order_trace_fraction'] for r in residual],'o-')
    axes[1].set_ylabel('Higher-order fraction of predicted trace')
    for ax in axes: ax.set_xlabel('Noise / clean state standard deviation');ax.grid(alpha=.2)
    fig.suptitle('Van der Pol: 4,000 residual realizations per cell')
    for ext in ('png','pdf'): fig.savefig(figures/f'covariance_validation.{ext}',dpi=220)
    plt.close(fig)
    if 'replication' in result:
        fig,ax=plt.subplots(figsize=(6.5,3.8),constrained_layout=True)
        for method in result['replication']['manifest']['methods']:
            selected=[r for r in result['replication']['summary'] if r['method']==method]
            ax.loglog([r['M'] for r in selected],[r['error_quantiles_25_50_75_90'][1] for r in selected],
                      'o-',label=method.replace('_',' '))
        ax.set(xlabel='Independent aligned recordings M',ylabel='Median relative coefficient error',
               title='Oracle support: averaging and stacking are different estimators')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
        for ext in ('png','pdf'): fig.savefig(figures/f'replication_comparison.{ext}',dpi=220)
        plt.close(fig)
    print(json.dumps({'pilot':result['pilot'],'covariance_attempts':result['covariance']['attempts'],
                      'diagnostics':result['diagnostic_selection_counts'],
                      'extended':result['extended']['completion']},indent=2))


if __name__=='__main__':
    analyze()
