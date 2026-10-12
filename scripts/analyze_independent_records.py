"""Summarize fixed-noise repeated-record inference, including coverage uncertainty."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]


def load(p):return json.loads(p.read_text(encoding='utf-8'))


def wilson(k,n):
    z=norm.ppf(.975);p=k/n;d=1+z*z/n
    c=(p+z*z/(2*n))/d;r=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(c-r),float(c+r)]


def main():
    p=load(ROOT/'configs/independent_records.json');stages=[f'independent_records_{i}' for i in range(p['shards'])]
    rows=[];hashes=[];pooled=[]
    for stage in stages:
        folder=ROOT/f'results/kaggle_{stage}';assert load(folder/'completion.json')['status']=='complete'
        rows+=load(folder/'rows.json');hashes+=load(folder/'input_hashes.json');pooled.append(load(folder/'pooled_residual_moments.json'))
    R=p['monte_carlo_repetitions'];assert sorted(r['repetition'] for r in hashes)==list(range(R))
    q=np.array(p['truth']);rng=np.random.default_rng(2702600);groups=[]
    for n in p['record_counts']:
        for method in p['comparators']:
            method={'known-variance averaged-moment GMM':'gmm_known','estimated-variance averaged-moment GMM':'gmm_estimated','known-variance stacked least squares':'stacked_known'}[method]
            group=sorted([r for r in rows if r['n']==n and r['method']==method],key=lambda r:r['repetition'])
            assert [r['repetition'] for r in group]==list(range(R))
            estimates=np.array([r['coefficients'] for r in group]);e=(estimates-q)/np.linalg.norm(q);mean=e.mean(axis=0)
            squared=np.sum(e*e,axis=1);ids=rng.integers(0,R,(2000,R));rmse=float(np.sqrt(squared.mean()))
            record=dict(n=n,method=method,repetitions=R,relative_rmse=rmse,sqrt_n_rmse=float(np.sqrt(n)*rmse),
                relative_bias=float(np.linalg.norm(mean)),relative_dispersion=float(np.sqrt(np.mean(np.sum((e-mean)**2,axis=1)))),
                rmse_bootstrap_95=np.quantile(np.sqrt(squared[ids].mean(axis=1)),[.025,.975]).tolist())
            if method.startswith('gmm'):
                covered=np.array([r['marginal_covered'] for r in group]);joint=sum(r['joint_covered'] for r in group)
                sd=np.array([r['standard_errors'] for r in group]);V=np.array([r['covariance'] for r in group])
                record.update(marginal_coverage=covered.mean(axis=0).tolist(),
                    marginal_coverage_wilson_95=[wilson(k,R) for k in covered.sum(axis=0)],
                    joint_coverage=joint/R,joint_coverage_wilson_95=wilson(joint,R),
                    empirical_sd=estimates.std(axis=0,ddof=1).tolist(),mean_estimated_sd=sd.mean(axis=0).tolist(),
                    sandwich_relative_rmse=float(np.sqrt(np.trace(V.mean(axis=0)))/np.linalg.norm(q)),
                    tau_quantiles=np.quantile([r['estimated_tau'] for r in group],[.05,.5,.95]).tolist(),
                    variance_mean=np.mean([r['estimated_variances'] for r in group],axis=0).tolist(),
                    variance_rmse=np.sqrt(np.mean((np.array([r['estimated_variances'] for r in group])-np.array(load(ROOT/f'results/kaggle_{stages[0]}/geometry.json')['sigma'])**2)**2,axis=0)).tolist())
            groups.append(record)
    geometry=load(ROOT/f'results/kaggle_{stages[0]}/geometry.json')
    for stage in stages:assert load(ROOT/f'results/kaggle_{stage}/geometry.json')==geometry
    count=sum(r['records'] for r in pooled);s=sum(np.array(r['sum']) for r in pooled);ss=sum(np.array(r['sum_outer']) for r in pooled)
    covariance=(ss-np.outer(s,s)/count)/(count-1);C=np.array(geometry['full_covariance'])
    report=dict(completed=True,protocol=p,groups=groups,geometry=geometry,
        independent_records_simulated=count,monte_carlo_repetitions=R,all_planned_repetitions_retained=True,
        empirical_full_residual_covariance_relative_error=float(np.linalg.norm(covariance-C)/np.linalg.norm(C)),
        empirical_mean_weak_residual=(s/count).tolist(),
        output_sha256={stage:hashlib.sha256((ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip').read_bytes()).hexdigest() for stage in stages},
        intervals='RMSE: 2000 repetition bootstrap draws, seed 2702600; coverage: binomial Wilson intervals; record counts paired within repetition')
    (ROOT/'reports/independent_records_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.2));colors=['#0072B2','#009E73','#D55E00']
    for m,label,color in zip(['gmm_known','gmm_estimated','stacked_known'],['GMM, known variance','GMM, estimated variance','Stacked LS, known variance'],colors):
        g=[r for r in groups if r['method']==m];ns=np.array([r['n'] for r in g]);v=np.array([r['relative_rmse'] for r in g]);ci=np.array([r['rmse_bootstrap_95'] for r in g])
        axes[0].errorbar(ns,v,yerr=[v-ci[:,0],ci[:,1]-v],fmt='o-',color=color,capsize=2,label=label)
    axes[0].set(xscale='log',yscale='log',xlabel='Independent records, n',ylabel='Relative coefficient RMSE');axes[0].legend(fontsize=6.5)
    g=[r for r in groups if r['method']=='gmm_estimated'];ns=[r['n'] for r in g]
    for j in range(4):axes[1].plot(ns,[r['marginal_coverage'][j] for r in g],'o-',lw=.9,ms=3,label=f'Coefficient {j+1}')
    axes[1].plot(ns,[r['joint_coverage'] for r in g],'s--',color='black',lw=.9,ms=3,label='Joint ellipsoid')
    minimum=min(min(r['marginal_coverage']+[r['joint_coverage']]) for r in g)
    axes[1].axhline(.95,color='grey',ls=':');axes[1].set(xscale='log',xlabel='Independent records, n',ylabel='Empirical 95% coverage',ylim=(max(0,minimum-.05),1.005));axes[1].legend(fontsize=6.5)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    for suffix in ['png','pdf']:fig.savefig(ROOT/f'manuscript/figures/publication/independent_records.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig);print(json.dumps(dict(records=count,groups=groups,geometry=geometry)))


if __name__=='__main__':main()
