"""Summarize recorded likelihood, covariance and Lorenz checks."""
import json, hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts import analyze_efficiency_followup as follow
ROOT=Path(__file__).resolve().parents[1]

def load(p):return json.loads(p.read_text(encoding='utf-8'))

def main():
    follow.RNG=np.random.default_rng(2502600)
    stages=['efficiency_fisher_v2','efficiency_noise_v2','efficiency_lorenz']
    for s in stages:assert load(ROOT/f'results/kaggle_{s}/completion.json')['status']=='complete'
    fisher=load(ROOT/'results/kaggle_efficiency_fisher_v2/rows.json')
    weak=[]
    for i in range(4):weak+=load(ROOT/f'results/kaggle_efficiency_power_{i}/rows.json')
    groups={m:sorted([r for r in fisher if r['method']==m],key=lambda r:r['seed']) for m in ['fisher_known_initial','fisher_unknown_initial']}
    groups.update({m:sorted([r for r in weak if r['weight_model']==m],key=lambda r:r['seed']) for m in ['ou','iid']})
    for rows in groups.values():
        assert [r['seed'] for r in rows]==list(range(2000))
        assert [r['input_sha256'] for r in rows]==[r['input_sha256'] for r in groups['ou']]
    report={'fisher':{m:follow.summary(rows) for m,rows in groups.items()},
        'fisher_paired':{m:follow.paired(groups['ou'],groups[m]) for m in groups if m!='ou'},
        'discrete_information':load(ROOT/'results/kaggle_efficiency_fisher_v2/discrete_information.json')}
    refined=[]
    refinement_stages=['efficiency_refine_0','efficiency_refine_1']
    if all((ROOT/f'results/kaggle_{s}/completion.json').exists() for s in refinement_stages):
        for s in refinement_stages:
            assert load(ROOT/f'results/kaggle_{s}/completion.json')['status']=='complete'
            refined+=load(ROOT/f'results/kaggle_{s}/rows.json')
        report['refinement']={}
        for step in [2,3]:
            rows=sorted([r for r in refined if r['step']==step],key=lambda r:r['seed'])
            assert [r['seed'] for r in rows]==list(range(2000))
            assert [r['input_sha256'] for r in rows]==[r['input_sha256'] for r in groups['ou']]
            report['refinement'][f'joint_step_{step}']=follow.summary(rows)
        stages+=refinement_stages
    noise=load(ROOT/'results/kaggle_efficiency_noise_v2/rows.json');assert len(noise)==768
    original=[]
    for i in range(2):original+=load(ROOT/f'results/kaggle_efficiency_vdp_final_{i}/rows.json')
    report['noise']=[]
    for level in [.005,.05]:
        baseline=sorted([r for r in original if r['level']==level and r['windows']==480 and r['weight_model']=='ou'],key=lambda r:r['seed'])
        for m in ['tau_half','tau_double','fitted_shape']:
            rows=sorted([r for r in noise if r['level']==level and r['method']==m],key=lambda r:r['seed'])
            assert [r['seed'] for r in rows]==list(range(128))
            tau=np.array([r['estimated_tau'] for r in rows]);sigma=np.array([r['estimated_sigma'] for r in rows])
            report['noise'].append(dict(level=level,method=m,**follow.summary(rows),
                paired_vs_true_shape=follow.paired(baseline,rows),converged=sum(r['converged'] for r in rows),
                tau_quantiles=np.quantile(tau,[.05,.5,.95]).tolist(),
                sigma_median=np.median(sigma,axis=0).tolist(),
                microergodic_median=np.median(sigma*sigma/tau[:,None],axis=0).tolist(),
                tau_bound_hits=int(np.sum((tau<.02001)|(tau>2.999)))))
    lorenz=load(ROOT/'results/kaggle_efficiency_lorenz/rows.json');assert [r['seed'] for r in lorenz]==list(range(256))
    truth=np.array([10.,28.,8/3]);e=(np.array([r['coefficients'] for r in lorenz])-truth)/np.linalg.norm(truth)
    lengths=np.linalg.norm(e,axis=1);assert np.allclose(lengths,[r['relative_error'] for r in lorenz])
    ids=follow.RNG.integers(0,len(e),(5000,len(e)));mean=e.mean(axis=0)
    report['lorenz']=dict(n=256,relative_rmse=float(np.sqrt(np.mean(lengths**2))),relative_bias=float(np.linalg.norm(mean)),
        relative_dispersion=float(np.sqrt(np.mean(np.sum((e-mean)**2,axis=1)))),
        rmse_bootstrap_95=np.quantile(np.sqrt(np.mean(lengths[ids]**2,axis=1)),[.025,.975]).tolist(),
        converged=sum(r['converged'] for r in lorenz),pilot_clips=sum(r['pilot_clipped'] for r in lorenz),
        information=load(ROOT/'results/kaggle_efficiency_lorenz/information.json'),
        geometry=load(ROOT/'results/kaggle_efficiency_lorenz/geometry.json'))
    report['intervals']='5000 seed-paired bootstrap draws; RNG 2502600; all recorded paths retained'
    report['execution_sources']={s:load(ROOT/f'results/kaggle_{s}/manifest.json')['source_sha256'] for s in stages}
    report['output_sha256']={s:hashlib.sha256((ROOT/f'results/kaggle_{s}_download/wsindy_{s}_outputs.zip').read_bytes()).hexdigest() for s in stages}
    (ROOT/'reports/information_followup_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    fig,ax=plt.subplots(figsize=(6.4,3.3))
    methods=['fisher_known_initial','fisher_unknown_initial','ou','iid']
    labels=['One step\nknown initial','One step\nunknown initial','Weak GLS\nOU, 480 windows','Weak GLS\niid, 480 windows']
    plotted=dict(report['fisher'])
    if 'refinement' in report:
        methods.insert(2,'joint_step_3');labels.insert(2,'Three steps\nunknown initial');plotted.update(report['refinement'])
    x=np.arange(len(methods))
    values=np.array([plotted[m]['relative_rmse'] for m in methods]);ci=np.array([plotted[m]['rmse_bootstrap_95'] for m in methods])
    ax.errorbar(x,values*1000,yerr=np.array([values-ci[:,0],ci[:,1]-values])*1000,fmt='o',color='#0072B2',capsize=4,label='Empirical RMSE, 95% path interval')
    ax.axhline(report['discrete_information']['known']['relative_rmse']*1000,color='#D55E00',ls=':',label='Known-initial information reference')
    ax.axhline(report['discrete_information']['unknown']['relative_rmse']*1000,color='#009E73',ls='--',label='Unknown-initial information reference')
    ax.set_xticks(x,labels);ax.set_ylabel('Relative RMSE (×10⁻³)')
    ax.set_ylim(1.1,2.03);ax.legend(fontsize=7);ax.spines[['top','right']].set_visible(False);fig.tight_layout()
    for suffix in ['png','pdf']:fig.savefig(ROOT/f'manuscript/figures/publication/likelihood_comparison.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig)
    print(json.dumps({k:report[k] for k in ['fisher','noise','lorenz']}))

if __name__=='__main__':main()
