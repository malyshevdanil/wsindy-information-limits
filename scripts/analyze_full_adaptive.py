"""Analyze completed full-author outputs; all simulation ran on Kaggle."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
NAMES={'fixed_raw_mstls':'Fixed raw MSTLS','fixed_hermite_mstls':'Fixed Hermite MSTLS',
       'wsindy':'Adaptive author WSINDy','wmsindy':'Adaptive author WmSINDy','wendy':'Adaptive author WENDy (known support)'}
NAMES.update({'fixed_raw_known_support':'Fixed raw, known support',
              'fixed_hermite_known_support':'Fixed Hermite, known support'})
KNOWN={'wendy','fixed_raw_known_support','fixed_hermite_known_support'}


def interval(values,rng):
    a=np.array(values,float)
    samples=np.median(a[rng.integers(0,len(a),(5000,len(a)))],axis=1)
    return np.quantile(samples,[.025,.975]).tolist()


def wilson(k,n):
    z=1.959963984540054;p=k/n;c=1+z*z/n
    mid=(p+z*z/(2*n))/c;half=z/c*np.sqrt(p*(1-p)/n+z*z/(4*n*n))
    return [mid-half,mid+half]


def main():
    cfg=json.loads((ROOT/'configs/full_adaptive.json').read_text());n=cfg['repeats']
    rows=[]
    for stage in ('adaptive_classical',*[f'adaptive_wm_{i}' for i in range(4)]):
        folder=ROOT/f'results/kaggle_{stage}'
        manifest=json.loads((folder/'manifest.json').read_text())
        if not manifest.get('completed'):raise ValueError('Uncompleted stage')
        if manifest['protocol']!=cfg:raise ValueError('Protocol mismatch')
        rows+=json.loads((folder/'rows.json').read_text())
    if len(rows)!=len(NAMES)*4*n:raise ValueError('Missing planned fits')
    hashes={}
    for row in rows:
        key=row['case'];hashes.setdefault(key,set()).add(row['input_sha256'])
    if any(len(h)!=1 for h in hashes.values()):raise ValueError('Unpaired algorithm inputs')
    grouped=[];rng=np.random.default_rng(925000)
    for regime in cfg['regimes']:
        for h in cfg['h']:
            for method,name in NAMES.items():
                a=[r for r in rows if (r['regime'],r['h'],r['method'])==(regime,h,method)]
                if len(a)!=n or len({r['seed'] for r in a})!=n:raise ValueError('Invalid record count')
                ok=[r for r in a if r['status']=='ok'];e=[r['error'] for r in ok]
                k=sum(bool(r.get('support_correct',False)) for r in a)
                grouped.append({'regime':regime,'h':h,'method':method,'name':name,'n':n,'completed':len(ok),
                                'failed':n-len(ok),'median':float(np.median(e)) if e else None,
                                'median_ci':interval(e,rng) if e else None,
                                'support_count':k if method not in KNOWN else None,
                                'support_ci':wilson(k,n) if method not in KNOWN else None})
    paired=[]
    for regime in cfg['regimes']:
        for method in NAMES:
            coarse={r['seed']:r for r in rows if (r['regime'],r['h'],r['method'])==(regime,max(cfg['h']),method)}
            fine={r['seed']:r for r in rows if (r['regime'],r['h'],r['method'])==(regime,min(cfg['h']),method)}
            seeds=[s for s in coarse if coarse[s]['status']==fine[s]['status']=='ok']
            diff=[fine[s]['error']-coarse[s]['error'] for s in seeds]
            paired.append({'regime':regime,'method':method,'paired_n':len(diff),
                           'median_fine_minus_coarse':float(np.median(diff)) if diff else None,
                           'ci':interval(diff,rng) if diff else None})
    comparisons=[]
    for regime in cfg['regimes']:
        for h in cfg['h']:
            for left,right in [('wsindy','fixed_hermite_mstls'),('wmsindy','wsindy'),
                               ('wmsindy','fixed_hermite_mstls'),('wendy','fixed_hermite_known_support'),
                               ('wendy','fixed_raw_known_support')]:
                a={r['seed']:r for r in rows if (r['regime'],r['h'],r['method'])==(regime,h,left)}
                b={r['seed']:r for r in rows if (r['regime'],r['h'],r['method'])==(regime,h,right)}
                seeds=[s for s in a if a[s]['status']==b[s]['status']=='ok']
                diff=[a[s]['error']-b[s]['error'] for s in seeds]
                comparisons.append({'regime':regime,'h':h,'left':left,'right':right,'paired_n':len(diff),
                                    'median_left_minus_right':float(np.median(diff)) if diff else None,
                                    'ci':interval(diff,rng) if diff else None})
    report={'protocol':cfg,'fit_count':len(rows),'paired_input_cases':len(hashes),'groups':grouped,'mesh_pairs':paired,
            'algorithm_pairs':comparisons,
            'bootstrap_replicates':5000,'bootstrap_seed':925000,'claim_scope':'16 paths per cell; supplied support for WENDy; one oscillator; no universal ranking'}
    (ROOT/'reports/full_adaptive_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    for (row,col),ax in np.ndenumerate(axes):
        regime=cfg['regimes'][col]
        for method,name in NAMES.items():
            if (method in KNOWN)!=(row==1):continue
            g=sorted([g for g in grouped if g['regime']==regime and g['method']==method],key=lambda g:g['h'])
            ax.errorbar([g['h'] for g in g],[g['median'] for g in g],
                        yerr=np.array([[g['median']-g['median_ci'][0] for g in g],[g['median_ci'][1]-g['median'] for g in g]]),
                        marker='o',capsize=3,label=name)
        ax.set_xscale('log');ax.set_yscale('log');ax.set_title(('Physical OU' if regime=='physical_ou' else 'Independent errors')+(' / known support' if row else ' / discovery'))
        ax.set_xlabel('Sampling interval h');ax.grid(alpha=.2);ax.set_xticks(cfg['h'],[str(h) for h in cfg['h']])
        ax.legend(fontsize=7,loc='best')
        if col==0:ax.set_ylabel('Median relative coefficient error')
    fig.tight_layout()
    figure=ROOT/'manuscript/figures/publication/full_adaptive.png'
    fig.savefig(figure,dpi=220,bbox_inches='tight');plt.close(fig)
    lines=['### 7.6. Complete adaptive author pipelines','',
      'A separate frozen comparison executes the complete pinned author WSINDy, WENDy and WmSINDy algorithms [@messenger2021; @bortz2023; @lopez2026]. WSINDy selects its test-function widths spectrally, places centers by the gradient, applies automatic SMAF smoothing and library scaling, and selects the MSTLS threshold. WENDy estimates noise, selects four adaptive radii and an SVD basis, and runs its full covariance iteration and stopping rules. WmSINDy jointly learns coefficients and noise with forward/backward RK4 and weak-residual losses, updates windows from the denoised data, and prunes and refits coefficients after each of eight cycles of 5000 Adam steps. Versions and compatibility adapters are detailed in Supplement S6.','',
      'The benchmark uses Van der Pol, five-percent state-scaled noise, fixed duration 20, correlation time 0.3, and 16 new independent paths per condition. Meshes 0.04 and 0.01 are nested. All methods receive identical observations verified by SHA256; discovery methods share the nine nonconstant monomials of degrees one through three. Fixed-window MSTLS is rerun on this library. WENDy receives the prespecified correct four-term structure and is assessed as parameter estimation. Hermite MSTLS uses the known calibration variance, whereas author pipelines estimate or reconstruct noise from observations.','',
      'Table 5. Discovery pipelines at h=0.01; 16 independent records per row. Median intervals are path bootstraps. Support counts include all planned fits.','',
      '| Error model | Method | Median error [95% interval] | Exact support | Failed fits |',
      '|---|---|---:|---:|---:|']
    for g in grouped:
        if g['h']!=min(cfg['h']) or g['method'] in KNOWN:continue
        med=f"{g['median']:.4f} [{g['median_ci'][0]:.4f}, {g['median_ci'][1]:.4f}]" if g['median'] is not None else 'n/a'
        support=f"{g['support_count']}/{n}" if g['support_count'] is not None else 'supplied'
        lines.append(f"| {'OU' if g['regime']=='physical_ou' else 'Independent'} | {g['name']} | {med} | {support} | {g['failed']}/{n} |")
    lines+=['','Table 6. Parameter estimation with supplied support at h=0.01; 16 independent records per row.','',
            '| Error model | Method | Median error [95% interval] | Failed fits |','|---|---|---:|---:|']
    for g in grouped:
        if g['h']!=min(cfg['h']) or g['method'] not in KNOWN:continue
        med=f"{g['median']:.4f} [{g['median_ci'][0]:.4f}, {g['median_ci'][1]:.4f}]" if g['median'] is not None else 'n/a'
        lines.append(f"| {'OU' if g['regime']=='physical_ou' else 'Independent'} | {g['name']} | {med} | {g['failed']}/{n} |")
    def finest(regime,method):
        return next(g for g in grouped if (g['regime'],g['h'],g['method'])==(regime,min(cfg['h']),method))
    ou_wm=finest('physical_ou','wmsindy');ou_ws=finest('physical_ou','wsindy')
    iid_wm=finest('iid','wmsindy');iid_fixed=finest('iid','fixed_hermite_mstls')
    effect=next(g for g in comparisons if (g['regime'],g['h'],g['left'],g['right'])==('physical_ou',min(cfg['h']),'wmsindy','wsindy'))
    interpretation=(f"On the finer OU grid, WmSINDy recovered support on {ou_wm['support_count']}/{n} paths with median error {ou_wm['median']:.4f}, "
        f"whereas adaptive WSINDy recovered {ou_ws['support_count']}/{n} with median {ou_ws['median']:.4f}. "
        f"The paired median error difference, WmSINDy minus WSINDy, was {effect['median_left_minus_right']:.4f} "
        f"with a descriptive 95% bootstrap interval [{effect['ci'][0]:.4f}, {effect['ci'][1]:.4f}]. "
        f"Under independent errors, WmSINDy recovered support on {iid_wm['support_count']}/{n} paths, "
        f"compared with {iid_fixed['support_count']}/{n} for fixed Hermite MSTLS. Its wide median interval reflects selection failures despite the full training budget. "
        "Thus the favorable OU result does not establish uniform superiority across observation regimes. These unadjusted comparisons are descriptive and limited to the frozen example.")
    lines+=['',interpretation,'',
      'The supplied-support comparison likewise concerns finite-sample accuracy: adaptive WENDy improves the OU coefficient median relative to these fixed-window fits, although its released noise model assumes independent observations. This empirical behavior does not establish calibrated uncertainty under OU error. Better window selection or reconstruction can improve accuracy within the finite information budget; it does not contradict the full-path infill bound.','',
      f"All {len(rows)} planned fits are archived; {sum(g['failed'] for g in grouped)} computational failures occurred. Error medians condition on completed fits. Wilson intervals, coarse-grid values, and paired mesh differences are retained in the reproducibility report. With 16 paths, support fractions have substantial binomial uncertainty. These conditions compare complete algorithms on one oscillator and do not establish a general ranking or an asymptotic plateau.",'',
      '![Figure 7. Complete adaptive author pipelines and matched fixed-window comparisons. Error bars are 95% path bootstrap intervals with n=16 per cell. WENDy uses known support; the discovery methods do not.](figures/publication/full_adaptive.png)','']
    section='\n'.join(lines)
    target=ROOT/'manuscript/article_revised.md';text=target.read_text(encoding='utf-8')
    start='<!-- FULL_ADAPTIVE_START -->';end='<!-- FULL_ADAPTIVE_END -->'
    block=start+'\n\n'+section+'\n'+end+'\n\n'
    if start in text:
        left=text.index(start);right=text.index(end)+len(end);text=text[:left]+block.rstrip()+text[right:]
    else:text=text.replace('## 8. Discussion and limitations',block+'## 8. Discussion and limitations')
    if 'Efficient weak-form estimation under correlated measurement error' not in text:
        target.write_text(text,encoding='utf-8')
    print(json.dumps({'fits':len(rows),'cases':len(hashes),'failures':sum(g['failed'] for g in grouped)}))


if __name__=='__main__':main()
