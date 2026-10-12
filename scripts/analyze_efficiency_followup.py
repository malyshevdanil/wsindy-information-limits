"""Analyze recorded paired paths; do not run numerical experiments locally."""
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pypandoc
ROOT=Path(__file__).resolve().parents[1]
RNG=np.random.default_rng(2402600)


def load(path):return json.loads(path.read_text(encoding='utf-8'))


def summary(rows):
    errors=(np.array([r['coefficients'] for r in rows])-np.array([1.,-1.,1.,-1.]))/2
    lengths=np.linalg.norm(errors,axis=1);mean=errors.mean(axis=0)
    indices=RNG.integers(0,len(rows),(5000,len(rows)))
    rmse=float(np.sqrt(np.mean(lengths**2)))
    bias=float(np.linalg.norm(mean));disp=float(np.sqrt(np.mean(np.sum((errors-mean)**2,axis=1))))
    assert abs(rmse**2-bias**2-disp**2)<1e-14
    assert np.allclose(lengths,[r['relative_error'] for r in rows],rtol=1e-14,atol=1e-15)
    return dict(n=len(rows),relative_rmse=rmse,relative_bias=bias,relative_dispersion=disp,
        bias_vector=mean.tolist(),median=float(np.median(lengths)),
        rmse_bootstrap_95=np.quantile(np.sqrt(np.mean(lengths[indices]**2,axis=1)),[.025,.975]).tolist())


def paired(left,right):
    assert [r['seed'] for r in left]==[r['seed'] for r in right]
    a=np.array([r['relative_error'] for r in left]);b=np.array([r['relative_error'] for r in right])
    ids=RNG.integers(0,len(a),(5000,len(a)))
    delta=np.sqrt(np.mean(b[ids]**2,axis=1))-np.sqrt(np.mean(a[ids]**2,axis=1))
    point=float(np.sqrt(np.mean(b*b))-np.sqrt(np.mean(a*a)))
    return dict(n=len(a),right_minus_left_rmse=point,paired_bootstrap_95=np.quantile(delta,[.025,.975]).tolist())


def figure():
    data=np.load(ROOT/'results/kaggle_efficiency_power_0/optimal_test_curves.npz')
    t=data['t'];fig,axes=plt.subplots(1,3,figsize=(7.5,2.8),gridspec_kw={'width_ratios':[2,1,1]})
    styles=[('exact','#0072B2',r'Generalized $\phi_*$','-'),('smooth','#009E73','Compact cutoff, 0.1','--'),('bell','#D55E00','Ordinary width-3 bell',':')]
    for ax in axes:
        for key,color,label,ls in styles:ax.plot(t,data[key],color=color,label=label,ls=ls,lw=1.5)
        # Zero extension makes the endpoint traces explicit, rather than connecting across a jump.
        ax.plot([-.12,0],[0,0],color='#0072B2');ax.plot([4,4.12],[0,0],color='#0072B2')
        for end,index in [(0,0),(4,-1)]:
            ax.plot([end,end],[0,data['exact'][index]],color='#0072B2',lw=1.4)
            ax.plot(end,0,'o',mfc='white',mec='#0072B2',ms=4)
            ax.plot(end,data['exact'][index],'o',color='#0072B2',ms=3)
        ax.axhline(0,color='.7',lw=.6);ax.set_xlabel('Time');ax.spines[['top','right']].set_visible(False)
    axes[0].set_xlim(-.12,4.12);axes[0].set_title('Full observation interval');axes[0].set_ylabel('Normalized test function')
    axes[0].legend(fontsize=6.5,loc='upper left')
    axes[1].set_xlim(-.015,.16);axes[1].set_title('Initial endpoint')
    axes[2].set_xlim(3.84,4.015);axes[2].set_title('Final endpoint')
    for ax in axes[1:]:
        keep=(t>=ax.get_xlim()[0])&(t<=ax.get_xlim()[1]);values=np.concatenate([data[k][keep] for k,_,_,_ in styles]+[np.array([0.])])
        margin=max(np.ptp(values)*.12,.01);ax.set_ylim(values.min()-margin,values.max()+margin)
    fig.tight_layout();folder=ROOT/'manuscript/figures/publication';folder.mkdir(parents=True,exist_ok=True)
    for suffix in ('png','pdf'):fig.savefig(folder/f'optimal_test_function.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig)


def main():
    stages=[f'efficiency_power_{i}' for i in range(4)]
    rows=[]
    for stage in stages:
        folder=ROOT/f'results/kaggle_{stage}'
        assert load(folder/'completion.json')['status']=='complete'
        rows+=load(folder/'rows.json')
    groups={m:sorted([r for r in rows if r['weight_model']==m],key=lambda r:r['seed']) for m in ('ou','iid')}
    for records in groups.values():assert [r['seed'] for r in records]==list(range(2000))
    assert all(r['input_sha256']==s['input_sha256'] for r,s in zip(groups['ou'],groups['iid']))
    old=[]
    for i in range(2):old+=load(ROOT/f'results/kaggle_efficiency_vdp_final_{i}/rows.json')
    lookup={(r['seed'],r['weight_model']):r for r in old if r['level']==.005 and r['windows']==480}
    discrepancy=max(np.linalg.norm(np.array(r['coefficients'])-np.array(lookup[r['seed'],r['weight_model']]['coefficients'])) for r in rows if r['seed']<128)
    assert discrepancy<1e-12
    report=dict(protocol=load(ROOT/'configs/efficiency_followup.json'),power={m:summary(groups[m]) for m in groups},
        power_paired=paired(groups['ou'],groups['iid']),additional_1872_paired=paired(groups['ou'][128:],groups['iid'][128:]),
        original_128_coefficient_max_discrepancy=float(discrepancy),all_power_converged=all(r['converged'] for r in rows),
        power_pilot_clips=sum(r['pilot_clipped'] for r in rows),power_max_iterations=max(r['iterations'] for r in rows),
        intervals='5000 seed-paired path bootstrap draws, RNG 2402600; fixed 2000 paths, including original 128; follow-up planned after inspection of original study')
    match_stages=[f'efficiency_match_{i}' for i in range(4)]
    if not all((ROOT/f'results/kaggle_{s}/completion.json').exists() for s in match_stages):
        match_stages=['efficiency_match'] if (ROOT/'results/kaggle_efficiency_match/completion.json').exists() else []
    if match_stages:
        matched=[];input_folders={}
        for stage in match_stages:
            matched_folder=ROOT/f'results/kaggle_{stage}'
            assert load(matched_folder/'completion.json')['status']=='complete'
            stage_rows=load(matched_folder/'rows.json');matched+=stage_rows
            input_folders.update({r['seed']:matched_folder for r in stage_rows})
        matched=sorted(matched,key=lambda r:r['seed'])
        assert [r['seed'] for r in matched]==list(range(128))
        assert all(r['status']=='ok' for r in matched)
        for r in matched:
            with np.load(input_folders[r['seed']]/f'inputs/seed_{r["seed"]:04d}.npz') as data:
                assert hashlib.sha256(data['t'].tobytes()+data['y'].tobytes()).hexdigest()==r['input_sha256']
        feasible={m:sorted([r for r in old if r['level']==.05 and r['windows']==480 and r['weight_model']==m],key=lambda r:r['seed']) for m in ('ou','iid')}
        report['matched']={'author_wendy':summary(matched),**{m:summary(feasible[m]) for m in feasible},
            'ou_minus_author':paired(matched,feasible['ou']),'iid_minus_author':paired(matched,feasible['iid']),
            'author_median_fit_seconds':float(np.median([r['seconds'] for r in matched]))}
        report['matched']['execution_stages']=match_stages
        report['matched']['unique_paths']=128
        report['matched']['duplicate_runs_not_counted']=True
        stages+=match_stages
    report['stages']=stages
    report['execution_sources']={s:load(ROOT/f'results/kaggle_{s}/manifest.json')['source_sha256'] for s in stages}
    report['output_sha256']={s:hashlib.sha256((ROOT/f'results/kaggle_{s}_download/wsindy_{s}_outputs.zip').read_bytes()).hexdigest() for s in stages}
    target=ROOT/'reports/efficiency_followup_results.json';target.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    figure()
    print(json.dumps({k:report[k] for k in ['power','power_paired','additional_1872_paired','original_128_coefficient_max_discrepancy']},ensure_ascii=True))
    if 'matched' in report:print(json.dumps(report['matched']))


if __name__=='__main__':main()
