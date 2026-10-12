"""Summarize rich-window heterogeneous-record scores against local information."""
import hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.analyze_independent_records import wilson
ROOT=Path(__file__).resolve().parents[1]


def load(p):return json.loads(p.read_text(encoding='utf-8'))


def main():
    p=load(ROOT/'configs/heterogeneous_records.json');rows=[];hashes=[];geometries=[];archives={}
    for i in range(p['shards']):
        stage=f'heterogeneous_records_{i}';folder=ROOT/f'results/kaggle_{stage}'
        assert load(folder/'completion.json')['status']=='complete';rows+=load(folder/'rows.json');hashes+=load(folder/'input_hashes.json');geometries.append(load(folder/'geometry.json'))
        archives[stage]=hashlib.sha256((ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip').read_bytes()).hexdigest()
    assert all(g==geometries[0] for g in geometries);R=p['monte_carlo_repetitions'];assert sorted(r['repetition'] for r in hashes)==list(range(R))
    truth=np.array(p['truth']);rng=np.random.default_rng(3400000);groups=[]
    for level in p['levels']:
        for n in p['record_counts']:
            for method in p['methods']:
                g=sorted([r for r in rows if r['level']==level and r['n']==n and r['method']==method],key=lambda r:r['repetition'])
                assert [r['repetition'] for r in g]==list(range(R));q=np.array([r['coefficients'] for r in g]);e=(q-truth)/np.linalg.norm(truth);mean=e.mean(axis=0)
                mse=np.sum(e*e,axis=1);ids=rng.integers(0,R,(3000,R));rmse=float(np.sqrt(mse.mean()));V=np.array([r['covariance'] for r in g])
                wald=np.array([r['wald_marginal_covered'] for r in g]);joint=sum(r['wald_joint_covered'] for r in g)
                group=dict(level=level,n=n,method=method,repetitions=R,relative_rmse=rmse,sqrt_n_rmse=float(np.sqrt(n)*rmse),
                    relative_bias=float(np.linalg.norm(mean)),relative_dispersion=float(np.sqrt(np.mean(np.sum((e-mean)**2,axis=1)))),
                    rmse_bootstrap_95=np.quantile(np.sqrt(mse[ids].mean(axis=1)),[.025,.975]).tolist(),
                    mean_local_information_relative_sd=float(np.mean([r['local_information_relative_sd'] for r in g])),
                    sqrt_n_local_information_relative_sd=float(np.sqrt(n)*np.mean([r['local_information_relative_sd'] for r in g])),
                    sandwich_relative_sd=float(np.sqrt(np.trace(V.mean(axis=0)))/2),
                    wald_marginal_coverage=wald.mean(axis=0).tolist(),wald_joint_coverage=joint/R,wald_joint_wilson_95=wilson(joint,R))
                if method.startswith('fission'):
                    covered=np.array([r['bootstrap_marginal_covered'] for r in g]);count=sum(r['bootstrap_joint_covered'] for r in g)
                    width=np.array([np.array(r['bootstrap_marginal_upper'])-r['bootstrap_marginal_lower'] for r in g]);waldwidth=2*1.959963984540054*np.sqrt(np.diagonal(V,axis1=1,axis2=2))
                    group.update(bootstrap_marginal_coverage=covered.mean(axis=0).tolist(),bootstrap_joint_coverage=count/R,
                        bootstrap_joint_wilson_95=wilson(count,R),bootstrap_marginal_wilson_95=[wilson(k,R) for k in covered.sum(axis=0)],
                        median_interval_width_ratio=np.median(width/waldwidth,axis=0).tolist())
                groups.append(group)
    report=dict(completed=True,protocol=p,geometry=geometries[0],groups=groups,output_sha256=archives,
        all_planned_repetitions_retained=True,records_per_level=R*max(p['record_counts']),
        instrument_training_is_independent_and_fixed=True,covariance_known=True)
    (ROOT/'reports/heterogeneous_records_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.5),sharex=True)
    colors={'oracle120':'#777777','fission120':'#009E73','fission480':'#0072B2'}
    labels={'oracle120':'Oracle instruments, 120 tests','fission120':'Independent instruments, 120 tests','fission480':'Independent instruments, 480 tests'}
    for row,level in enumerate(p['levels']):
        for method in p['methods']:
            g=[r for r in groups if r['level']==level and r['method']==method];ns=np.array([r['n'] for r in g]);value=np.array([r['sqrt_n_rmse'] for r in g]);ci=np.array([r['rmse_bootstrap_95'] for r in g])*np.sqrt(ns)[:,None]
            axes[row,0].errorbar(ns,value,yerr=[value-ci[:,0],ci[:,1]-value],fmt='o-',color=colors[method],capsize=2,ms=3,label=labels[method])
            if method.startswith('fission'):
                v=np.array([r['bootstrap_joint_coverage'] for r in g]);c=np.array([r['bootstrap_joint_wilson_95'] for r in g])
                axes[row,1].errorbar(ns,v,yerr=[v-c[:,0],c[:,1]-v],fmt='o-',color=colors[method],capsize=2,ms=3,label=labels[method])
        g=[r for r in groups if r['level']==level and r['method']=='fission480']
        axes[row,0].plot([r['n'] for r in g],[r['sqrt_n_local_information_relative_sd'] for r in g],'k:',label='Local nuisance-information reference')
        axes[row,0].set(ylabel=r'$\sqrt{n}$ relative RMSE',title=f'{100*level:g}% measurement noise')
        axes[row,1].axhline(.95,color='grey',ls=':');axes[row,1].set(ylabel='Joint bootstrap 95% coverage',ylim=(.7,1.005))
    for ax in axes.flat:ax.set(xscale='log',xlabel='Independent records, n');ax.spines[['top','right']].set_visible(False)
    axes[0,0].legend(fontsize=5.7);axes[0,1].legend(fontsize=5.7);fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(ROOT/f'manuscript/figures/publication/heterogeneous_gmm.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig);print(json.dumps(dict(completed=True,groups=groups,geometry=geometries[0])))


if __name__=='__main__':main()
