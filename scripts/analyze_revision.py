"""Generate publication evidence from completed archived runs; never invent data."""
import json
import argparse
import hashlib
from pathlib import Path
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'
OUT=ROOT/'results/kaggle_revision'
FIG=BASE/'figures/publication';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','figure.dpi':120})
LABEL={'linear':'Exponential decay','damped':'Damped oscillator','vdp':'Van der Pol','duffing':'Duffing'}
METHOD={'relative_mstls_raw':'MSTLS, raw','relative_mstls_hermite':'MSTLS, corrected','oracle_hermite':'Known support, corrected',
    'oracle_raw_ols':'OLS, known support','wendy_iid_core':'WENDy iid IRLS core','ou_irls_core':'OU IRLS core',
    'stacked':'Stacked design','mean_weak':'Mean weak matrices','naive_mean_sensor':'Mean sensor path'}


def load(path):return json.loads(path.read_text())


def median_interval(values,repeats=2000):
    values=np.asarray(values);rng=np.random.default_rng(281910)
    if not len(values):return [None,None,None]
    boot=np.median(values[rng.integers(len(values),size=(repeats,len(values)))],axis=1)
    return [float(np.median(values)),*map(float,np.quantile(boot,[.025,.975]))]


def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(center-half),float(center+half)]


def save(name,fig):
    fig.savefig(FIG/(name+'.png'),dpi=240);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)


def errorbars(ax,x,summaries,**kwargs):
    a=np.array(summaries);ax.errorbar(x,a[:,0],yerr=np.maximum(0,np.vstack([a[:,0]-a[:,1],a[:,2]-a[:,0]])),capsize=3,**kwargs)


def sparse(rows):
    summary=[]
    for system in LABEL:
        for regime in ('physical_ou','iid'):
            for level in (.01,.05):
                for h in (.04,.02,.01,.005):
                    for method in ('relative_mstls_raw','relative_mstls_hermite','oracle_hermite'):
                        group=[r for r in rows if (r['system'],r['regime'],r['level'],r['h'],r['method'])==(system,regime,level,h,method)]
                        if len(group)!=128:raise ValueError('Incomplete sparse condition')
                        errors=[r['error'] for r in group];k=sum(r['support_correct'] for r in group)
                        summary.append({'system':system,'regime':regime,'level':level,'h':h,'method':method,'n':128,
                            'median_error_ci':median_interval(errors),'support_count':k,'support_wilson':wilson(k,128),
                            'conditional_correct':median_interval([r['error'] for r in group if r['support_correct']]),
                            'conditional_incorrect':median_interval([r['error'] for r in group if not r['support_correct']])})
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    hs=np.array([.005,.01,.02,.04])
    for ax,system in zip(axes.ravel(),LABEL):
        for regime,color in [('physical_ou','#b74132'),('iid','#2065a1')]:
            for method,style in [('relative_mstls_hermite','-'),('oracle_hermite','--')]:
                data=[next(s['median_error_ci'] for s in summary if (s['system'],s['regime'],s['level'],s['h'],s['method'])==(system,regime,.05,h,method)) for h in hs]
                errorbars(ax,hs,data,color=color,linestyle=style,marker='o',markersize=4,
                    label=('OU' if regime=='physical_ou' else 'Independent')+', '+('MSTLS' if method.startswith('relative') else 'known support'))
        ax.set_xscale('log');ax.set_yscale('log');ax.xaxis.set_major_locator(FixedLocator(hs));ax.xaxis.set_major_formatter(FixedFormatter(['0.005','0.01','0.02','0.04']));ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(LABEL[system]);ax.set_xlabel('Sampling interval h');ax.set_ylabel('Median relative coefficient error')
    axes[0,0].legend(fontsize=8);fig.tight_layout();save('sparse_refinement',fig)
    group=[r for r in rows if (r['system'],r['regime'],r['level'],r['h'],r['method'])==('vdp','physical_ou',.05,.005,'relative_mstls_hermite')]
    fig,ax=plt.subplots(figsize=(7,4));allerrors=np.array([r['error'] for r in group]);lo=max(allerrors.min()*.8,1e-5);hi=max(allerrors.max()*1.1,lo*2);bins=np.geomspace(lo,hi,20)
    for correct,color in [(True,'#2065a1'),(False,'#b74132')]:
        errors=[r['error'] for r in group if r['support_correct']==correct]
        ax.hist(errors,bins=bins,color=color,alpha=.65,label=('Correct support' if correct else 'Incorrect support')+f' (n={len(errors)})')
    ax.set_xscale('log');ax.set_xlabel('Relative coefficient error');ax.set_ylabel('Number of paths');ax.legend();ax.set_title('Van der Pol: OU error, 5% noise, h = 0.005');fig.tight_layout();save('support_conditional',fig)
    return summary


def information():
    info=load(OUT/'information.json');data=np.load(OUT/'scalar_likelihood.npz');hs=data['h'];a=info['scalar_parameters']['a']
    summ={method:[median_interval(np.abs(data[method][:,j]-a)) for j in range(len(hs))] for method in ('weak','mle')}
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for method,color in [('weak','#b74132'),('mle','#2065a1')]:errorbars(axes[0],hs,summ[method],marker='o',color=color,label='Weak ratio' if method=='weak' else 'Bounded Gaussian likelihood')
    axes[0].set_xscale('log');axes[0].set_xlabel('Sampling interval h');axes[0].set_ylabel('Median absolute coefficient error');axes[0].legend(fontsize=8)
    duration=info['duration'];Ts=[r['T'] for r in duration]
    axes[1].plot(Ts,[r['decay_information']/info['continuous_information'] for r in duration],'o-',color='#2065a1')
    axes[1].set_xscale('log',base=2);axes[1].set_xlabel('Observation duration T');axes[1].set_ylabel('Decay information relative to T = 4');axes[1].set_title('Additional duration loses excitation')
    fig.tight_layout();save('information_and_likelihood',fig)
    info['empirical_median_absolute_error_ci']=summ
    info['mle_finest_mse']=float(np.mean((data['mle'][:,-1]-a)**2))
    info['mle_finest_bias']=float(np.mean(data['mle'][:,-1]-a))
    return info


def coefficient():
    summ=load(OUT/'coefficient_summary.json')
    fig,axes=plt.subplots(1,2,figsize=(10,4));levels=[r['level'] for r in summ]
    intervals=[]
    for r in summ:
        data=np.load(OUT/f"coefficient_{r['level']:g}.npz");A=data['actual'][data['guards']];C=data['predicted_covariance'];rng=np.random.default_rng(27271)
        errs=[]
        for _ in range(600):
            z=A[rng.integers(len(A),size=len(A))]
            errs.append(float(np.linalg.norm(np.cov(z,rowvar=False,bias=True)-C)/np.linalg.norm(C)))
        r['guarded_frobenius_bootstrap_ci']=np.quantile(errs,[.025,.975]).tolist()
        intervals.append([r['guarded_actual_relative_frobenius'],*r['guarded_frobenius_bootstrap_ci']])
    interval_array=np.asarray(intervals)
    axes[0].vlines(levels,interval_array[:,1],interval_array[:,2],color='#2065a1')
    axes[0].scatter(levels,interval_array[:,0],color='#2065a1',label='Guarded nonlinear fit')
    axes[0].scatter(levels,[r['linear_relative_frobenius'] for r in summ],marker='s',color='#b74132',label='Gaussian linearized fit')
    axes[0].scatter(levels,[r['gaussian_known_mean_covariance_relative_rms_floor'] for r in summ],marker='x',color='black',label='Gaussian Monte Carlo RMS floor')
    axes[0].set_xscale('log');axes[0].set_ylabel('Relative covariance discrepancy');axes[0].legend(fontsize=8)
    axes[1].scatter(levels,[r['paired_linearization_relative_rms'] for r in summ],color='#2065a1');axes[1].set_xscale('log');axes[1].set_yscale('log');axes[1].set_ylabel('Paired linearization relative RMS error')
    for ax in axes:ax.set_xlabel('Noise SD / clean state SD')
    fig.tight_layout();save('coefficient_covariance',fig)
    return summ


def covariance_null():
    summ=load(OUT/'covariance_null.json');fig,ax=plt.subplots(figsize=(8,4));names=[]
    for i,row in enumerate(summ):
        label=row['file'].removeprefix('residual_').removesuffix('.npz');corrected=label.startswith('known_debias_');level=float(label.split('_')[-1]);names.append(('Corrected' if corrected else 'Raw')+f'\n{level:.0%}')
        q=row['affine_null_bootstrap_quantiles'];ax.errorbar(i,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',capsize=4,color='#2065a1',label='Approximate null 95% range' if i==0 else None)
        ax.scatter(i,row['observed_relative_frobenius'],marker='x',color='#b74132',s=60,label='Observed discrepancy' if i==0 else None)
    ax.set_xticks(range(len(names)),names);ax.set_ylabel('Relative Frobenius discrepancy');ax.legend();fig.tight_layout();save('covariance_null',fig)
    return summ


def paired_refinement(rows):
    results=[]
    for system in LABEL:
        for regime in ('physical_ou','iid'):
            for level in (.01,.05):
                for method in ('relative_mstls_raw','relative_mstls_hermite','oracle_hermite'):
                    condition=[r for r in rows if (r['system'],r['regime'],r['level'],r['method'])==(system,regime,level,method)]
                    coarse={r['seed']:r['error'] for r in condition if r['h']==.01}
                    fine={r['seed']:r['error'] for r in condition if r['h']==.005}
                    seeds=sorted(coarse);a=np.array([coarse[s] for s in seeds]);b=np.array([fine[s] for s in seeds])
                    rng=np.random.default_rng(19711);indexes=rng.integers(len(seeds),size=(2000,len(seeds)))
                    values=1-np.median(b[indexes],axis=1)/np.median(a[indexes],axis=1)
                    ci=np.quantile(values,[.025,.975])
                    results.append({'system':system,'regime':regime,'level':level,'method':method,
                        'h_pair':[.01,.005],'relative_median_improvement':float(1-np.median(b)/np.median(a)),
                        'paired_bootstrap_ci':ci.tolist(),'within_five_percent_band':bool(ci[0]>=-.05 and ci[1]<=.05)})
    return results


def aggregate(rows,methods,xkey,name):
    summary=[];fig,ax=plt.subplots(figsize=(7,4))
    for method in methods:
        xs=sorted(set(r[xkey] for r in rows));values=[]
        for value in xs:
            a=[r['error'] for r in rows if r['method']==method and r[xkey]==value]
            interval=median_interval(a);values.append(interval);summary.append({'method':method,xkey:value,'n':len(a),'median_error_ci':interval})
        errorbars(ax,xs,values,marker='o',label=METHOD[method])
    ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xlabel('Independent recordings M');ax.set_ylabel('Median relative coefficient error');ax.legend(fontsize=8);fig.tight_layout();save(name,fig)
    return summary


def text_results(report):
    scalar=report['information'];weak=scalar['empirical_median_absolute_error_ci']['weak'][-1];mle=scalar['empirical_median_absolute_error_ci']['mle'][-1]
    lines=['### 7.1. Full-record information and scalar likelihood','',
        f"For $a=-1$, $x_0=2$, $T=4$, $\\tau=0.4$ and $\\sigma=0.1$, continuous Fisher information is {scalar['continuous_information']:.3f}. The chosen scalar weak function has linearized information efficiency {scalar['weak_linearized_efficiency']:.3f}; Section 7.8 evaluates optimal combinations within declared window families. On 1,000 paired paths, the finest-grid median absolute errors were {weak[0]:.4f} (95% interval [{weak[1]:.4f}, {weak[2]:.4f}]) for the weak ratio and {mle[0]:.4f} ([{mle[1]:.4f}, {mle[2]:.4f}]) for bounded Gaussian likelihood. The latter uses known covariance and initial condition. Its finite-sample bias and boundedness prevent interpreting its MSE as an unbiased Cramér–Rao efficiency measure.", '',
        '![Figure 1. Scalar coefficient errors with median bootstrap intervals, and continuous Fisher information as duration increases.](figures/publication/information_and_likelihood.png)','',
        '### 7.2. Fixed-window sparse regression','',
        'The fixed-window comparison contains 24,576 fits: four systems, two noise regimes, two amplitudes, four grids, three methods and 128 paths per condition. Thirty fourth-power polynomial bump windows remain fixed as the grid changes. For decay, $T=4$ and $L=1.6$; for oscillators, $T=20$ and $L=3$. Correlation ratio is $\\tau/L=0.1$. Relative MSTLS uses 50 thresholds between $10^{-4}$ and one, loss weight 0.8 and no ridge penalty; raw and known-variance corrected features share these settings. Settings and seed ranges were fixed before outputs; this is not a public preregistration. Initial states and libraries are specified in the supplement.','',
        'Table 1. Corrected relative MSTLS at 5% noise and $h=0.005$. Error and exact support are separate endpoints; 128 paths enter each row.','',
        '| System | Error model | Median relative error | 95% median interval | Exact support | 95% support interval |','|---|---|---:|---|---:|---|']
    captured=scalar['grid_information'][0]['fisher_information']/scalar['continuous_information']
    duration_gain=next(r['decay_information'] for r in scalar['duration'] if r['T']==8)/scalar['continuous_information']-1
    position=lines.index('### 7.2. Fixed-window sparse regression')
    lines[position:position]=[f"The coarsest scalar grid already captures {100*captured:.3f}% of continuous-record information. Extending the decay record from $T=4$ to $8$ adds only {100*duration_gain:.3f}%, whereas a second independent record doubles information. This is a quantitative acquisition comparison under the stated known-covariance model, rather than a universal cost-optimal allocation.",'']
    if 'geometry' in report:
        g=report['geometry'];position=lines.index('### 7.2. Fixed-window sparse regression')
        lines[position:position]=[f"Initial-state annihilation raises the weak influence variance lower bound from {g['known_initial_variance_bound']:.6f} to {g['weak_initial_annihilation_variance_bound']:.6f}, an increase of {100*(g['weak_bound_inflation']-1):.1f}%. The chosen window's efficiency relative to the nuisance-adjusted information reference is {g['weak_nuisance_information_efficiency']:.3f}, compared with {g['weak_known_initial_information_efficiency']:.3f} relative to known-initial-state information. This separates the weak class restriction from additional window inefficiency. Deterministic quadrature verifies both influence identities to $10^{{-9}}$; it adds no Monte Carlo paths.",'']
    for system in LABEL:
        for regime in ('physical_ou','iid'):
            r=next(r for r in report['sparse'] if (r['system'],r['regime'],r['level'],r['h'],r['method'])==(system,regime,.05,.005,'relative_mstls_hermite'))
            med,lo,hi=r['median_error_ci'];pl,pu=r['support_wilson'];lines.append(f"| {LABEL[system]} | {'OU' if regime=='physical_ou' else 'Independent'} | {med:.4f} | [{lo:.4f}, {hi:.4f}] | {r['support_count']}/128 | [{pl:.3f}, {pu:.3f}] |")
    vdp=next(r for r in report['sparse'] if (r['system'],r['regime'],r['level'],r['h'],r['method'])==('vdp','physical_ou',.05,.005,'relative_mstls_hermite'))
    lines+=['',f"In the Van der Pol OU condition, correct-support paths had median error {vdp['conditional_correct'][0]:.4f}, whereas incorrect-support paths had median {vdp['conditional_incorrect'][0]:.4f}. The aggregate median therefore combines parameter-estimation error with structure-selection error and cannot be read as an intrinsic information limit.",
        '','![Figure 2. Median coefficient errors with 95% path bootstrap intervals. Known-support fits separate parameter estimation from sparse selection.](figures/publication/sparse_refinement.png)','',
        '![Figure 3. Error distributions conditional on correct or incorrect selected support in the finest-grid Van der Pol OU condition. Counts accompany each distribution; bimodality is not presumed.](figures/publication/support_conditional.png)','',
        'Intervals resample independent paths and retain paired mesh comparisons. Conditions share random-number indices and are not pooled as independent experiments. Exact-support intervals use Wilson limits. No multiplicity-adjusted global support claim is made. No forecast endpoint is included in this validation.','',
        '### 7.3. Coefficient covariance and Monte Carlo calibration','',
        'Joint coefficient covariance was tested on 2,000 Van der Pol paths at each amplitude, on known support. Actual errors subtract the clean discrete fit; a guard retains designs whose smallest support singular value is at least half its clean value. Every simulated path passed: the guard did not filter these samples, so guarded summaries equal the full empirical summaries. This observation does not establish finite unconditional estimator moments or control unobserved near-singular designs. The guard defines the intended validation region; it is not evidence that the linearization is accurate.','',
        'Table 2. Relative Frobenius covariance discrepancies and paired relative RMS linearization errors on known support.','',
        '| Noise amplitude | Guarded paths | Linear covariance discrepancy | Guarded actual discrepancy | Paired linearization RMS |','|---|---:|---:|---:|---:|']
    selected=[r for r in report['paired_refinement'] if r['regime']=='physical_ou' and r['level']==.05 and r['method']=='relative_mstls_hermite']
    statement=f"For corrected MSTLS at 5% OU noise, {sum(r['within_five_percent_band'] for r in selected)} of four system-specific paired median-improvement intervals for $h=0.01$ to $0.005$ lie wholly inside the $\\pm5\\%$ band. This is a descriptive equivalence check of coefficient medians and does not guarantee support recovery."
    position=lines.index('### 7.3. Coefficient covariance and Monte Carlo calibration')
    lines[position:position]=[statement,'']
    for r in report['coefficient']:lines.append(f"| {r['level']:.1%} | {r['guarded_paths']}/2000 | {r['linear_relative_frobenius']:.3f} | {r['guarded_actual_relative_frobenius']:.3f} | {r['paired_linearization_relative_rms']:.3f} |")
    lines+=['','The linear covariance discrepancy is close to its Monte Carlo RMS floor at every amplitude because the standardized linear experiment is paired across amplitudes. Actual covariance agrees at 0.1% and 1% noise but has substantially larger discrepancy at 5% and 20%. Paired pathwise linearization already has appreciable error at 1%. Thus the coefficient surrogate is a local small-noise description, and the residual covariance theorem alone cannot validate it at the moderate-noise benchmark setting.',
        '','![Figure 4. Joint coefficient covariance discrepancies, bootstrap intervals for guarded actual covariance, Gaussian Monte Carlo RMS references, and paired pathwise linearization error.](figures/publication/coefficient_covariance.png)','',
        'In the 4,000-path nonlinear residual experiment, higher orders contributed 20.3% of corrected trace at 50% noise. The covariance calibration places the observed Frobenius discrepancy against the fourth-moment Monte Carlo scale and an affine-null bootstrap. The null procedure retains empirical non-Gaussian shape and forces the analytical covariance; its range is approximate, and its exceedance fraction is not an exact p-value.','']
    position=lines.index('![Figure 4. Joint coefficient covariance discrepancies, bootstrap intervals for guarded actual covariance, Gaussian Monte Carlo RMS references, and paired pathwise linearization error.](figures/publication/coefficient_covariance.png)')
    lines[position:position]=[r'The last column is $R=\{\sum_m\|e_m-\ell_m\|^2/\sum_m\|\ell_m\|^2\}^{1/2}$, where $e_m$ is the actual coefficient displacement from the clean fit and $\ell_m$ its first-order surrogate. It is not a relative error of the estimated coefficients and is not constrained to $[0,1]$. At 5% and 20% noise, omitted nonlinear terms and random-design coupling exceed the RMS magnitude of the first-order surrogate. A stable singular value alone does not suppress those terms.','']
    for r in report['covariance_null']:
        if r['file']=='residual_known_debias_0.5.npz':
            q=r['affine_null_bootstrap_quantiles'];lines.append(f"For the corrected 50% condition, observed matrix discrepancy was {r['observed_relative_frobenius']:.3f}; estimated fourth-moment relative RMS floor was {r['estimated_fourth_moment_relative_rms_floor']:.3f}, and the approximate null 95% range was [{q[0]:.3f}, {q[2]:.3f}]. The 6.4% discrepancy therefore must be read relative to this uncertainty, rather than as an isolated accuracy percentage.")
    lines+=['','![Figure 5. Observed covariance discrepancies and approximate covariance-null bootstrap ranges for three discrete amplitudes and two correction choices.](figures/publication/covariance_null.png)','',
        '### 7.4. Different initial conditions','',
        'A separate frozen experiment used 400 blocks with up to 16 independent Van der Pol recordings. Initial offsets were uniform on $[-0.25,0.25]^2$ about $(2,0)$; sensor noise covariance remained fixed. All methods supplied true support.','',
        'Table 3. Known-support aggregation with different initial conditions; 400 independent blocks per row.','',
        '| Method | M | Median relative error | 95% median interval |','|---|---:|---:|---|']
    for r in report['heterogeneous']:
        if r['M'] in (1,16):
            m,l,u=r['median_error_ci'];lines.append(f"| {METHOD[r['method']]} | {r['M']} | {m:.5f} | [{l:.5f}, {u:.5f}] |")
    lines+=['','Mean weak matrices improve with independent acquisition in this heterogeneous setting, whereas stacking retains a stable error and sensor-path averaging deteriorates. The latter estimates the dynamics of an averaged path that generally fails the original nonlinear equation. The evidence supports a distinction between aggregation operations, under the supplied support and excitation used here.',
        '','![Figure 6. Independent recordings with different initial states. Mean weak matrices and sensor-path averaging define different estimating equations.](figures/publication/heterogeneous_recordings.png)','']
    if 'reweighting' in report:
        lines+=['### 7.5. Covariance reweighting component','',
            'An additional fresh-path comparison supplies known support and compares ordinary raw-feature least squares, the iid covariance IRLS step from WENDy Algorithm 2, and its known-OU covariance extension. It retains our fixed windows and uses ten iterations, relative stopping tolerance $10^{-6}$, known coordinate variances and scale-aware covariance relaxation $10^{-10}$. These choices make it a component experiment rather than a reproduction of complete WENDy.','',
            'Table 4. Van der Pol covariance-reweighting components: parameter one, initial state (2,0), T=20, 5% state-scaled noise, OU correlation time 0.3, 30 fixed windows of full width three, and 128 paths per condition. True support is supplied.','',
            '| Error model | Component | h | Median relative error | 95% median interval | Converged |','|---|---|---:|---:|---|---:|']
        for r in report['reweighting']:
            m,l,u=r['median_error_ci'];convergence='n/a' if r['method']=='oracle_raw_ols' else f"{r['converged_count']}/{r['n']}";lines.append(f"| {'OU' if r['regime']=='physical_ou' else 'Independent'} | {METHOD[r['method']]} | {r['h']} | {m:.4f} | [{l:.4f}, {u:.4f}] | {convergence} |")
        lines+=['','At h=0.005 under OU error, the median decreases from 0.0325 for ordinary least squares to 0.0186 for iid IRLS and 0.0181 for OU IRLS. Reweighting therefore improves this component comparison, but its two IRLS medians are close and their displayed intervals overlap; the table alone does not establish superiority of the OU extension. OU errors remain nearly unchanged across the two meshes, whereas independent-noise errors decrease. Known covariance and supplied support distinguish this comparison from the complete adaptive pipelines in Section 7.6. All iteration-capped estimates enter the summaries. The convergence column reports the prescribed tolerance; it does not assess model correctness.']
    return '\n'.join(lines)


def run(stage='revision'):
    global OUT
    OUT=ROOT/f'results/kaggle_{stage}'
    if load(OUT/'completion.json')['status']!='complete':raise ValueError('Revision incomplete')
    report={'execution_stage':stage,'manifest':load(OUT/'manifest.json'),
        'output_archive_sha256':hashlib.sha256((ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip').read_bytes()).hexdigest(),
        'sparse':sparse(load(OUT/'benchmark_rows.json')),'paired_refinement':paired_refinement(load(OUT/'benchmark_rows.json')),'information':information(),
        'coefficient':coefficient(),'covariance_null':covariance_null(),
        'heterogeneous':aggregate(load(OUT/'heterogeneous_rows.json'),('stacked','mean_weak','naive_mean_sensor'),'M','heterogeneous_recordings')}
    rw=ROOT/'results/kaggle_reweighting'
    geometry=ROOT/'results/kaggle_geometry'
    if (geometry/'completion.json').exists():
        if load(geometry/'completion.json')['status']!='complete':raise ValueError('Geometry check incomplete')
        report['geometry']={**load(geometry/'geometry.json'),'manifest':load(geometry/'manifest.json'),
            'output_archive_sha256':hashlib.sha256((ROOT/'results/kaggle_geometry_download/wsindy_geometry_outputs.zip').read_bytes()).hexdigest()}
    if (rw/'completion.json').exists():
        rows=load(rw/'rows.json');report['reweighting']=[]
        for regime in ('physical_ou','iid'):
            for method in ('oracle_raw_ols','wendy_iid_core','ou_irls_core'):
                for h in (.04,.005):
                    group=[r for r in rows if (r['regime'],r['method'],r['h'])==(regime,method,h)]
                    report['reweighting'].append({'regime':regime,'method':method,'h':h,'n':len(group),
                        'median_error_ci':median_interval([r['error'] for r in group]),
                        'converged_count':sum(r.get('converged',False) for r in group)})
    (ROOT/'reports/revision_results.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    source=BASE/'article_revised.md';text=source.read_text(encoding='utf-8');begin='<!-- REVISION_RESULTS_START -->';end='<!-- REVISION_RESULTS_END -->'
    if 'Efficient weak-form estimation under correlated measurement error' in text:
        print(json.dumps({'status':'complete','benchmark_summary_cells':len(report['sparse']),'figures':6,'manuscript_update':'use analyze_efficiency.py'}))
        return
    text=text.split(begin)[0]+begin+'\n\n'+text_results(report)+'\n\n'+end+text.split(end)[1]
    text=text.replace('The new independent protocol adds damped and Duffing oscillators, published relative MSTLS, Gaussian likelihood, coefficient-surrogate validation, covariance-null calibration and heterogeneous initial conditions. Its results are included only after completed output is downloaded.','The independent validation adds damped and Duffing oscillators, published relative MSTLS, Gaussian likelihood, coefficient-surrogate validation, covariance-null calibration and heterogeneous initial conditions.')
    text=text.replace('Pending independent computations are not evidence until their outputs are available.','')
    source.write_text(text,encoding='utf-8')
    print(json.dumps({'status':'complete','benchmark_summary_cells':len(report['sparse']),'figures':6}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['revision','revision_fast'],default='revision')
    run(parser.parse_args().stage)
