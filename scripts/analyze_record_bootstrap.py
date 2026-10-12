"""Paired Monte Carlo coverage and interval-size audit for record bootstrap-t."""
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
    c=load(ROOT/'configs/record_bootstrap.json');original=load(ROOT/'reports/independent_records_results.json')
    rows=[];hashes=[];original_hashes={};original_rows={};archives={};completion=[]
    for i in range(c['shards']):
        stage=f'record_bootstrap_{i}';folder=ROOT/f'results/kaggle_{stage}'
        completed=load(folder/'completion.json');assert completed['status']=='complete';completion.append(completed)
        rows+=load(folder/'rows.json');hashes+=load(folder/'input_hashes.json')
        archives[stage]=hashlib.sha256((ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip').read_bytes()).hexdigest()
        old=ROOT/f'results/kaggle_independent_records_{i}'
        original_hashes.update({r['repetition']:r['sufficient_statistics_sha256'] for r in load(old/'input_hashes.json')})
        original_rows.update({(r['repetition'],r['n']):r for r in load(old/'rows.json') if r['method']=='gmm_estimated'})
    R=c['monte_carlo_repetitions'];assert sorted(r['repetition'] for r in hashes)==list(range(R))
    assert all(r['sufficient_statistics_sha256']==original_hashes[r['repetition']] for r in hashes)
    groups=[];truth=np.array(original['protocol']['truth'])
    for n in c['record_counts']:
        g=sorted([r for r in rows if r['n']==n],key=lambda r:r['repetition']);assert [r['repetition'] for r in g]==list(range(R))
        for r in g:
            old=original_rows[(r['repetition'],n)]
            assert np.array_equal(r['coefficients'],old['coefficients']) and np.array_equal(r['covariance'],old['covariance'])
            lower=np.array(r['marginal_lower']);upper=np.array(r['marginal_upper'])
            assert np.array_equal((truth>=lower)&(truth<=upper),r['marginal_covered'])
            assert r['joint_covered']==(r['joint_statistic']<=r['joint_critical_value'])
        marginal=np.array([r['marginal_covered'] for r in g]);joint=np.array([r['joint_covered'] for r in g])
        oldjoint=np.array([original_rows[(r['repetition'],n)]['joint_covered'] for r in g])
        widths=np.array([np.array(r['marginal_upper'])-r['marginal_lower'] for r in g])
        oldwidths=np.array([2*1.959963984540054*np.sqrt(np.diag(r['covariance'])) for r in g])
        critical=np.array([r['joint_critical_value'] for r in g])
        group=dict(n=n,repetitions=R,marginal_coverage=marginal.mean(axis=0).tolist(),
            marginal_coverage_wilson_95=[wilson(k,R) for k in marginal.sum(axis=0)],
            joint_coverage=float(joint.mean()),joint_coverage_wilson_95=wilson(joint.sum(),R),
            wald_joint_coverage=float(oldjoint.mean()),paired_joint_coverage_gain=float((joint.astype(int)-oldjoint).mean()),
            newly_covered=int(np.sum(joint&~oldjoint)),newly_missed=int(np.sum(~joint&oldjoint)),
            mean_marginal_width=widths.mean(axis=0).tolist(),median_marginal_width_ratio=np.median(widths/oldwidths,axis=0).tolist(),
            joint_critical_value_quantiles=np.quantile(critical,[.05,.5,.95]).tolist(),
            median_joint_volume_ratio=float(np.median((critical/9.487729036781154)**2)))
        groups.append(group)
    report=dict(completed=True,protocol=c,groups=groups,output_sha256=archives,
        observations_identical_to_original_sha256=True,point_estimates_and_covariances_identical=True,
        bootstrap_fits=sum(r['bootstrap_fits'] for r in completion),failed_or_discarded_draws=0,
        maximum_batch_fit_discrepancy=max(r['maximum_batch_fit_discrepancy'] for r in completion))
    (ROOT/'reports/record_bootstrap_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.2));ns=[r['n'] for r in groups]
    wald=[r for r in original['groups'] if r['method']=='gmm_estimated']
    for j in range(4):
        axes[0].plot(ns,[r['marginal_coverage'][j] for r in groups],'o-',lw=1,ms=3,label=f'Bootstrap, coefficient {j+1}')
    for g,label,color in ((wald,'Wald joint','#D55E00'),(groups,'Bootstrap joint','#0072B2')):
        v=np.array([r['joint_coverage'] for r in g]);ci=np.array([r['joint_coverage_wilson_95'] for r in g])
        axes[0].errorbar(ns,v,yerr=[v-ci[:,0],ci[:,1]-v],fmt='s--',capsize=2,color=color,label=label)
    axes[0].axhline(.95,color='grey',ls=':');axes[0].set(xscale='log',xlabel='Independent records, n',ylabel='Empirical 95% coverage',ylim=(.45,1.005));axes[0].legend(fontsize=6)
    for j in range(4):axes[1].plot(ns,[r['median_marginal_width_ratio'][j] for r in groups],'o-',ms=3,label=f'Coefficient {j+1}')
    axes[1].axhline(1,color='grey',ls=':');axes[1].set(xscale='log',xlabel='Independent records, n',ylabel='Median bootstrap / Wald interval width');axes[1].legend(fontsize=6)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(ROOT/f'manuscript/figures/publication/record_bootstrap.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig);print(json.dumps(report))


if __name__=='__main__':main()
