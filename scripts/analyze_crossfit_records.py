"""Summarize all prespecified cross-fitting runs, including coverage and budget."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.analyze_independent_records import wilson
from scripts.validate_crossfit_runs import validate
ROOT=Path(__file__).resolve().parents[1]


def main(prefix='crossfit'):
    validation=validate([f'{prefix}_records_{i}' for i in range(4)],full=True)
    rows=[];pilots=[];manifests=[];archives={}
    for i in range(4):
        stage=f'{prefix}_records_{i}';folder=ROOT/f'results/kaggle_{stage}'
        rows+=json.loads((folder/'rows.json').read_text());pilots+=json.loads((folder/'pilot_rows.json').read_text());manifests.append(json.loads((folder/'manifest.json').read_text()))
        archives[stage]=hashlib.sha256((ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip').read_bytes()).hexdigest()
    p=manifests[0]['protocol'];geometry=manifests[0]['cache_geometry'];assert all(m['protocol']==p and m['cache_geometry']==geometry for m in manifests)
    R=p['monte_carlo_repetitions'];truth=np.array(p['truth']);rng=np.random.default_rng(3900000);groups=[]
    for level in p['levels']:
        for method in ('crossfit480','split1000_50'):
            for n in (p['record_counts'] if method=='crossfit480' else [1050]):
                g=sorted([r for r in rows if r['level']==level and r['method']==method and r['n']==n],key=lambda r:r['repetition']);assert len(g)==R
                error=(np.array([r['coefficients'] for r in g])-truth)/2;bias=error.mean(axis=0);mse=np.sum(error**2,axis=1);ids=rng.integers(0,R,(3000,R));V=np.array([r['covariance'] for r in g])
                joint=sum(r['wald_joint_covered'] for r in g)
                group=dict(level=level,n=n,method=method,repetitions=R,inference_records=g[0]['inference_records'],extra_training_records=g[0]['training_records_extra'],
                    relative_rmse=float(np.sqrt(mse.mean())),sqrt_n_rmse=float(np.sqrt(n*mse.mean())),relative_bias=float(np.linalg.norm(bias)),
                    relative_dispersion=float(np.sqrt(np.mean(np.sum((error-bias)**2,axis=1)))),
                    rmse_bootstrap_95=np.quantile(np.sqrt(mse[ids].mean(axis=1)),[.025,.975]).tolist(),
                    sqrt_n_local_reference=float(np.sqrt(n)*np.mean([r['local_information_relative_sd'] for r in g])),
                    wald_joint_coverage=joint/R,wald_joint_wilson_95=wilson(joint,R),
                    wald_marginal_coverage=np.mean([r['wald_marginal_covered'] for r in g],axis=0).tolist())
                if method=='crossfit480':
                    for prefix in ('refit','frozen'):
                        marginal=np.array([r[prefix+'_marginal_covered'] for r in g]);count=sum(r[prefix+'_joint_covered'] for r in g)
                        width=np.array([np.array(r[prefix+'_marginal_upper'])-r[prefix+'_marginal_lower'] for r in g]);waldwidth=2*1.959963984540054*np.sqrt(np.diagonal(V,axis1=1,axis2=2))
                        group.update({prefix+'_joint_coverage':count/R,prefix+'_joint_wilson_95':wilson(count,R),prefix+'_marginal_coverage':marginal.mean(axis=0).tolist(),
                            prefix+'_marginal_wilson_95':[wilson(k,R) for k in marginal.sum(axis=0)],prefix+'_median_width_ratio':np.median(width/waldwidth,axis=0).tolist()})
                    delta=np.array([int(r['refit_joint_covered'])-int(r['frozen_joint_covered']) for r in g])
                    group.update(refit_minus_frozen_joint_coverage=float(delta.mean()),refit_minus_frozen_paired_95=np.quantile(delta[ids].mean(axis=1),[.025,.975]).tolist())
                    pg=[r for r in pilots if r['level']==level and r['n']==n]
                    group.update(maximum_pilot_relative_error=max(r['maximum_relative_error'] for r in pg),pilot_box_clipped_coordinates=sum(r['clipped_coordinates'] for r in pg),cache_clipped_coordinates=sum(r['cache_clipped_coordinates'] for r in pg))
                assert np.isclose(group['relative_rmse']**2,group['relative_bias']**2+group['relative_dispersion']**2);groups.append(group)
    paired_budget=[]
    for level in p['levels']:
        a=sorted([r for r in rows if r['n']==1050 and r['level']==level and r['method']=='crossfit480'],key=lambda r:r['repetition'])
        b=sorted([r for r in rows if r['level']==level and r['method']=='split1000_50'],key=lambda r:r['repetition'])
        ea=np.array([r['relative_error']**2 for r in a]);eb=np.array([r['relative_error']**2 for r in b]);ids=rng.integers(0,R,(5000,R))
        delta=np.sqrt(eb[ids].mean(axis=1))-np.sqrt(ea[ids].mean(axis=1))
        paired_budget.append(dict(level=level,total_records=1050,split_minus_crossfit_rmse=float(np.sqrt(eb.mean())-np.sqrt(ea.mean())),paired_95=np.quantile(delta,[.025,.975]).tolist()))
    report=dict(completed=True,validation=validation,protocol=p,geometry=geometry,groups=groups,paired_budget=paired_budget,output_sha256=archives,
        all_planned_repetitions_retained=True,external_training_records=0,bootstrap_pilots_refitted=True,covariance_known=True,
        comparison_caveat='Previous frozen-pilot protocol also differs in pilot amplitude and instrument approximation; matched records do not isolate cross-fitting alone.')
    (ROOT/'reports/crossfit_records_results.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),sharex=True)
    for row,level in enumerate(p['levels']):
        g=[r for r in groups if r['level']==level and r['method']=='crossfit480'];ns=np.array([r['n'] for r in g]);v=np.array([r['sqrt_n_rmse'] for r in g]);ci=np.array([r['rmse_bootstrap_95'] for r in g])*np.sqrt(ns)[:,None]
        axes[row,0].errorbar(ns,v,yerr=[v-ci[:,0],ci[:,1]-v],fmt='o-',color='#0072B2',capsize=2,label='Five-fold cross-fitting')
        axes[row,0].plot(ns,[r['sqrt_n_local_reference'] for r in g],'k:',label='Local nuisance-information reference')
        b=next(r for r in groups if r['level']==level and r['method']=='split1000_50')
        axes[row,0].plot(1050,b['sqrt_n_rmse'],'^',color='#D55E00',label='1000 training + 50 inference')
        axes[row,0].set(ylabel=r'$\sqrt{N}$ relative RMSE',title=f'{100*level:g}% measurement noise')
        for prefix,color,label in [('refit','#0072B2','Refitted-pilot bootstrap'),('frozen','#009E73','Frozen-instrument bootstrap')]:
            v=np.array([r[prefix+'_joint_coverage'] for r in g]);ci=np.array([r[prefix+'_joint_wilson_95'] for r in g])
            axes[row,1].errorbar(ns,v,yerr=[v-ci[:,0],ci[:,1]-v],fmt='o-',color=color,capsize=2,label=label)
        axes[row,1].axhline(.95,color='grey',ls=':');axes[row,1].set(ylabel='Joint 95% coverage',ylim=(.65,1.005))
    for ax in axes.flat:ax.set(xscale='log',xlabel='Total records, N');ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(fontsize=5.5);axes[0,1].legend(fontsize=5.5);fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(ROOT/f'manuscript/figures/publication/crossfit_gmm.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig);print(json.dumps(dict(completed=True,groups=len(groups),paired_budget=paired_budget)))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--prefix',default='crossfit');main(a.parse_args().prefix)
