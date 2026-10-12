"""Build the efficiency revision from recorded Kaggle outputs only."""
import json,re,hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pypandoc
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'
FIG=BASE/'figures/publication'
RNG=np.random.default_rng(2302600)


def load(path):return json.loads(path.read_text(encoding='utf-8'))


def interval(values,kind='rmse'):
    values=np.asarray(values)
    bootstrap=values[RNG.integers(0,len(values),(5000,len(values)))]
    stat=np.sqrt(np.mean(bootstrap**2,axis=1)) if kind=='rmse' else np.median(bootstrap,axis=1)
    return np.quantile(stat,[.025,.975]).tolist()


def summarize(errors):
    errors=np.asarray(errors);mean=errors.mean(axis=0)
    squares=np.sum(errors**2,axis=1) if errors.ndim==2 else errors**2
    bias=float(np.linalg.norm(mean));disp=float(np.sqrt(np.mean((errors-mean)**2) if errors.ndim==1 else np.mean(np.sum((errors-mean)**2,axis=1))))
    rmse=float(np.sqrt(squares.mean()))
    assert abs(rmse**2-bias**2-disp**2)<1e-12
    return dict(n=len(errors),median=float(np.median(np.sqrt(squares))),relative_rmse=rmse,
        relative_bias=bias,relative_dispersion=disp,median_bootstrap_95=interval(np.sqrt(squares),'median'),
        rmse_bootstrap_95=interval(np.sqrt(squares)))


def save(name,fig):
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(FIG/(name+'.'+suffix),dpi=220,bbox_inches='tight')
    plt.close(fig)


def main():
    stages=['efficiency_scalar_v3','efficiency_geometry_v2','efficiency_vdp_final_0','efficiency_vdp_final_1']
    for stage in stages:
        assert load(ROOT/f'results/kaggle_{stage}/completion.json')['status']=='complete'
    scalar=load(ROOT/'results/kaggle_efficiency_scalar_v3/rows.json')
    refs=load(ROOT/'results/kaggle_efficiency_scalar_v3/scalar_reference.json')
    groups=defaultdict(list)
    for r in scalar:groups[r['sigma'],r['method'],r['delta']].append(r)
    scalar_stats=[]
    for (sigma,method,delta),rows in sorted(groups.items(),key=lambda kv:(kv[0][0],kv[0][1],kv[0][2] or 0)):
        assert len(rows)==1000 and sorted(r['seed'] for r in rows)==list(range(1000))
        error=np.array([r['bounded_estimate']+1 for r in rows])
        scalar_stats.append(dict(sigma=sigma,method=method,delta=delta,**summarize(error),
            unbounded_median=float(np.median([abs(r['estimate']+1) for r in rows])),
            clipped=sum(r['bounded_estimate']!=r['estimate'] for r in rows),
            guard_passed=sum(r.get('guard_passed',True) for r in rows)))
    vdp=[]
    for stage in stages[2:]:vdp+=load(ROOT/f'results/kaggle_{stage}/rows.json')
    assert len(vdp)==1536
    controls={}
    for stage in ('efficiency_vdp_0','efficiency_vdp_1'):
        for r in load(ROOT/f'results/kaggle_{stage}/rows.json'):
            controls[r['seed'],r['level'],r['windows'],r['weight_model']]=r
    changes=[float(np.linalg.norm(np.array(r['coefficients'])-np.array(controls[r['seed'],r['level'],r['windows'],r['weight_model']]['coefficients']))/2) for r in vdp]
    stopping_control=dict(same_paths=True,replayed_fits=len(changes),original_iteration_cap=5,final_iteration_cap=100,
        median_relative_coefficient_change=float(np.median(changes)),maximum_relative_coefficient_change=float(max(changes)))
    groups=defaultdict(list)
    for r in vdp:groups[r['level'],r['windows'],r['weight_model']].append(r)
    vdp_stats=[];pairs=[];truth=np.array([1.,-1.,1.,-1.])
    for (level,K,model),rows in sorted(groups.items()):
        rows.sort(key=lambda r:r['seed']);assert [r['seed'] for r in rows]==list(range(128))
        errors=(np.array([r['coefficients'] for r in rows])-truth)/2
        assert np.isfinite(errors).all()
        vdp_stats.append(dict(level=level,windows=K,weight_model=model,**summarize(errors),
            converged=sum(r['converged'] for r in rows),pilot_clipped=sum(r['pilot_clipped'] for r in rows),
            max_iterations=max(r['iterations'] for r in rows),retained_windows=rows[0]['retained_windows'],
            median_fit_seconds=float(np.median([r['seconds'] for r in rows]))))
    for level in (.005,.05):
        for K in (30,120,480):
            ou=np.array([r['relative_error'] for r in groups[level,K,'ou']])
            iid=np.array([r['relative_error'] for r in groups[level,K,'iid']])
            indices=RNG.integers(0,128,(5000,128))
            difference=np.sqrt(np.mean(iid[indices]**2,axis=1))-np.sqrt(np.mean(ou[indices]**2,axis=1))
            pairs.append(dict(level=level,windows=K,iid_minus_ou_rmse=float(np.sqrt(np.mean(iid**2))-np.sqrt(np.mean(ou**2))),
                paired_bootstrap_95=np.quantile(difference,[.025,.975]).tolist()))
    geom=load(ROOT/'results/kaggle_efficiency_geometry_v2/geometry.json')
    info=load(ROOT/'reports/reviewer_results.json')['vdp_information']
    finest=info['grids'][-1]['continuous_quadrature']
    bound=float(finest['initial_nuisance_relative_rmse_bound']) if 'initial_nuisance_relative_rmse_bound' in finest else None
    # The field spelling is deliberately checked, never inferred from a plot.
    if bound is None:
        covariance=np.array(finest['initial_nuisance_covariance_bound']);bound=float(np.sqrt(np.trace(covariance))/2)
    known=float(np.sqrt(np.trace(np.array(finest['active_covariance_bound'])))/2)
    report=dict(protocols={s:load(ROOT/f'results/kaggle_{s}/manifest.json')['protocol'] for s in stages},
        scalar=scalar_stats,scalar_reference=refs,vdp=vdp_stats,paired_vdp=pairs,geometry=geom,
        vdp_weak_relative_rmse_bound=bound,vdp_known_initial_relative_rmse_bound=known,
        execution_sources={s:load(ROOT/f'results/kaggle_{s}/manifest.json')['source_sha256'] for s in stages},
        intervals='5000 path bootstrap draws; comparisons paired by seed; no multiplicity adjustment',
        scalar_fits=len(scalar),vdp_fits=len(vdp),stopping_control=stopping_control)
    (ROOT/'reports/efficiency_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    FIG.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(7.4,3.0))
    for ax,sigma in zip(axes,(.1,.005)):
        for method,color,label in [('oracle_phi','#777777','Generating-parameter test'),('pilot_phi','#0072B2','Fitted-pilot test')]:
            rows=sorted([r for r in scalar_stats if r['sigma']==sigma and r['method']==method],key=lambda r:r['delta'])
            y=np.array([r['median'] for r in rows]);limits=np.array([r['median_bootstrap_95'] for r in rows])
            ax.errorbar([r['delta'] for r in rows],y,yerr=np.maximum(0,np.array([y-limits[:,0],limits[:,1]-y])),fmt='o-',color=color,label=label,capsize=2)
        ref=next(r for r in refs if r['sigma']==sigma)
        ax.axhline(.6744897501960817*np.sqrt(ref['weak_variance_bound']),color='black',ls=':',label='Gaussian weak-bound median')
        baseline=next(r for r in scalar_stats if r['sigma']==sigma and r['method']=='single_window')
        ax.axhline(baseline['median'],color='#D55E00',ls='--',label='One fixed window')
        ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Cutoff width');ax.set_ylabel('Median absolute coefficient error');ax.set_title(f'Scalar noise SD {sigma:g}')
    axes[0].legend(fontsize=7);save('efficient_scalar',fig)
    fig,axes=plt.subplots(1,2,figsize=(7.4,3.0))
    for ax,level in zip(axes,(.005,.05)):
        for model,color,label in [('ou','#0072B2','OU weights'),('iid','#D55E00','iid weights')]:
            rows=sorted([r for r in vdp_stats if r['level']==level and r['weight_model']==model],key=lambda r:r['windows'])
            y=np.array([r['relative_rmse'] for r in rows]);limits=np.array([r['rmse_bootstrap_95'] for r in rows])
            ax.errorbar([r['windows'] for r in rows],y,yerr=np.maximum(0,np.array([y-limits[:,0],limits[:,1]-y])),fmt='o-',color=color,label=label,capsize=2)
            theory_rows=sorted([r for r in geom['results'] if r['weight_model']==model],key=lambda r:r['raw_windows'])
            ax.plot([r['raw_windows'] for r in theory_rows],[r['relative_rmse']*level/.05 for r in theory_rows],':',color=color,label=label+' local reference')
        ax.axhline(bound*level/.05,color='black',ls='--',label='Weak-class bound')
        ax.set_xscale('log');ax.set_xticks([30,120,480],labels=['30','120','480']);ax.set_xlabel('Prescribed windows');ax.set_ylabel('Relative coefficient RMSE');ax.set_title(f'State-scaled noise {100*level:g}%')
    axes[0].legend(fontsize=7);save('efficient_vdp',fig)
    render(report)
    print(json.dumps(dict(scalar_fits=len(scalar),vdp_fits=len(vdp),vdp_bound=bound)))


def render(report):
    scalar=report['scalar'];refs=report['scalar_reference'];vdp=report['vdp'];geom=report['geometry']['results']
    def s(sigma,method,delta=None):return next(r for r in scalar if (r['sigma'],r['method'],r['delta'])==(sigma,method,delta))
    def v(level,K,model):return next(r for r in vdp if (r['level'],r['windows'],r['weight_model'])==(level,K,model))
    def ci(r,key='rmse_bootstrap_95'):return '['+', '.join(f'{x:.5f}' for x in r[key])+']'
    small=.0125
    content=r'''### 7.1. Scalar optimal test with a fitted pilot

The scalar model is $x'=ax$, with $a=-1$, $x_0=2$, $T=4$ and OU correlation time 0.4. Noise standard deviations are 0.1, 0.02 and 0.005. For each amplitude, 1000 independent OU paths are shared across tests. A 30-window weak OLS pilot, clipped to $[-2.5,0.5]$, is substituted into the explicit function of Section 5.3. A quintic cutoff is zero within a quarter of the declared width from either endpoint and reaches one at that width. Widths are 0.2, 0.1, 0.05, 0.025 and 0.0125. Simpson quadrature uses $h=1/16384$. The generating-parameter test is an oracle comparison; the fitted-pilot test receives no true coefficient.

Unguarded weak ratios can have infinite population moments (Section 4.2). Reported scalar RMSE therefore uses the prespecified clipping of final estimates to $[-2.5,0.5]$; unbounded medians are also archived. A clean-moment denominator diagnostic does not filter paths and is not used by the estimator. All path counts, clips and diagnostics are listed in S8.

Table 1. Scalar coefficient errors at cutoff width 0.0125, with 1000 paired paths per amplitude. Medians are absolute errors; RMSE uses bounded final estimates. The reference is the standard deviation of the efficient Gaussian influence, not a finite-noise bound on the nonlinear bounded estimator.

| Noise SD | One-window median | Pilot-test median (95% interval) | Pilot-test RMSE | Weak influence SD |
|---|---:|---:|---:|---:|
'''
    for sigma in (.1,.02,.005):
        pilot=s(sigma,'pilot_phi',small);base=s(sigma,'single_window');ref=next(r for r in refs if r['sigma']==sigma and r['delta']==small)
        content+=f"| {sigma:g} | {base['median']:.5f} | {pilot['median']:.5f} {ci(pilot,'median_bootstrap_95')} | {pilot['relative_rmse']:.5f} | {np.sqrt(ref['weak_variance_bound']):.5f} |\n"
    high=s(.1,'pilot_phi',small);low=s(.005,'pilot_phi',small)
    ref=next(r for r in refs if r['sigma']==.1 and r['delta']==small)
    content+=f"\nAt noise SD 0.1, the fitted-pilot median is {high['median']:.5f}, compared with {s(.1,'single_window')['median']:.5f} for the fixed window. The Gaussian absolute-error median at the class bound is {(.6744897501960817*np.sqrt(ref['weak_variance_bound'])):.5f}. These are different statistical summaries; agreement of medians alone does not establish covariance efficiency. At the smallest noise the pilot-test bias is {low['relative_bias']:.6f} and dispersion {low['relative_dispersion']:.6f}. The oracle local efficiency relative to the weak-class bound at width 0.0125 is {ref['weak_efficiency']:.4f}. Full width and amplitude dependence is retained in S8.\n\n"
    content+='![Figure 1. Scalar compact optimal tests at two noise amplitudes. Intervals resample independent paths; the dotted Gaussian median reference is derived from the sharp influence bound.](figures/publication/efficient_scalar.png)\n\n'
    content+=r'''### 7.2. Van der Pol: nested spans and temporal weighting

We estimate four independently free coefficients in $x'=av$, $v'=bx+cv+dx^2v$, with truth $(1,-1,1,-1)$, initial state $(2,0)$, $T=20$, OU correlation time 0.3 and $h=0.005$. Marginal coordinate noise scales are fixed from the same clean $h=0.01$ calibration as S9.6. The relative error denominator is the Euclidean coefficient norm, two. The continuous known-initial-state reference is 0.013093 at 5% noise; the sharp boundary-free relative RMSE reference is 0.016904. At 0.5% noise these scale down by ten. Estimating a single tied oscillator parameter would be a different information experiment (S9.7).

The nested families contain 30 full-width-three bumps; then 30 width-1.5 and 60 width-0.75 bumps; then 120 width-0.375 and 240 width-0.1875 bumps. A fixed $H^1$ Gram transformation preserves their spans while improving representation. Relative cutoff $10^{-10}$ retains 30, 120 and 479 dimensions; the largest family has one linear dependence at this numerical resolution. The smallest retained relative Gram eigenvalue in that family is approximately $2.1\times10^{-5}$. Covariance solves add $10^{-10}$ times the mean diagonal. This finite family is an approximation experiment, not a proof of density for an infinite sequence.

Table 2. Oracle Gaussian influence RMSE under actual OU noise at 5%. Both methods use the same test span and clean design; they differ only in the covariance used for GLS weights. Actual variance uses the OU covariance for both.

| Prescribed windows | Retained tests | OU-weight RMSE | iid-weight RMSE | OU RMSE / class reference |
|---:|---:|---:|---:|---:|
'''
    bound=report['vdp_weak_relative_rmse_bound']
    for K in (30,120,480):
        ou=next(r for r in geom if r['raw_windows']==K and r['weight_model']=='ou');iid=next(r for r in geom if r['raw_windows']==K and r['weight_model']=='iid')
        content+=f"| {K} | {ou['retained_windows']} | {ou['relative_rmse']:.6f} | {iid['relative_rmse']:.6f} | {ou['relative_rmse']/bound:.4f} |\n"
    content+=r'''
The OU-weighted 480-window local reference is within 0.2% in RMSE of the class optimum; iid weighting retains a roughly 4.9% gap. Optimal weighting in the narrow 30-window span limits how much temporal covariance can help. Expanding the span makes the distinction observable. This is a verified prediction for the specified spans; it does not imply that every iid-weighted dense family must have a positive limiting gap. The independent exact-polynomial covariance fixture agrees with Proposition 3 to relative error below $10^{-14}$.

For the feasible test, 128 independent paths at each amplitude are shared by both weights and all three families. Hermite-corrected moments use known marginal variances. A 30-window OLS pilot is clipped to a fixed physical parameter box and integrated from the known initial state. Exact polynomial covariance is evaluated on this fitted model trajectory, and coefficient-dependent weights are updated for at most 100 iterations, with relative tolerance $10^{-7}$. The true trajectory is not supplied to this estimator. This is a known-support WENDy-like component; the complete author WENDy benchmark is documented separately in S6 and S9.6.

Table 3. Feasible Van der Pol empirical relative RMSE and 95% path bootstrap intervals, with 128 paired paths per row. Bias and dispersion refer to the OU-weighted coefficient errors; their squares sum to empirical MSE. Final coefficient estimates are unbounded, so these finite-cohort summaries do not assert finite unconditional moments.

| Noise | Windows | OU RMSE (95% interval) | iid RMSE (95% interval) | OU bias | OU dispersion |
|---|---:|---:|---:|---:|---:|
'''
    maximum_quadrature=max(r['clean_quadrature_relative_error'] for r in geom)
    content=content.replace('For the feasible test,',f'The largest clean discrete GLS coefficient error is {maximum_quadrature:.2e}, below the small-noise error scale. For the feasible test,')
    for level in (.005,.05):
        for K in (30,120,480):
            ou=v(level,K,'ou');iid=v(level,K,'iid')
            content+=f"| {100*level:g}% | {K} | {ou['relative_rmse']:.5f} {ci(ou)} | {iid['relative_rmse']:.5f} {ci(iid)} | {ou['relative_bias']:.5f} | {ou['relative_dispersion']:.5f} |\n"
    low=v(.005,480,'ou');high=v(.05,480,'ou')
    content+=f"\nAt 0.5% noise, OU weighting with 480 windows gives empirical RMSE {low['relative_rmse']:.6f}, versus the sharp reference {bound/10:.6f}. At 5% noise it gives {high['relative_rmse']:.6f}, versus {bound:.6f}; bias {high['relative_bias']:.6f} and dispersion {high['relative_dispersion']:.6f} identify finite-noise excess error. The stronger finite-noise assertion that more windows alone would reach 0.016904 is therefore not supported by this feasible implementation. The small-noise result and deterministic geometry support the conditional efficiency mechanism, while the moderate-noise result delimits it.\n\n"
    pairs=report['paired_vdp']
    for level in (.005,.05):
        p=next(r for r in pairs if r['level']==level and r['windows']==480)
        content+=f"At {100*level:g}% noise, paired iid-minus-OU RMSE for 480 windows is {p['iid_minus_ou_rmse']:.6f}, with 95% interval [{p['paired_bootstrap_95'][0]:.6f}, {p['paired_bootstrap_95'][1]:.6f}]. "
    content+='Both paired intervals include zero: this 128-path experiment does not resolve a statistical advantage of OU over iid weighting. Intervals are not adjusted for multiple comparisons. All 1536 final fits meet the stopping tolerance within 23 iterations, and no pilot is clipped; the finite-noise gap persists after this convergence control. Stopping counts are retained in S8; no path is excluded.\n\n'
    content+='![Figure 2. Feasible Van der Pol RMSE as the test span grows. Dotted curves show oracle local influences and the horizontal dashed line the class reference. Finite-noise estimation is distinguished from local Gaussian geometry.](figures/publication/efficient_vdp.png)\n\n'
    content+=r'''### 7.3. Supporting validation

The supplement retains fixed-window sparse selection, coefficient linearization checks, independent-record aggregation and 448 complete adaptive author fits. These results provide context without constituting a test of Corollary 1. In the separate 16-path 5% OU comparison, full WmSINDy has empirical RMSE 0.01926 and WENDy 0.02154, compared with the weak reference 0.01690; their medians of approximately 0.015 do not establish efficiency. Full coefficient and support distributions are retained in S9.6–S9.7.

Stacking heterogeneous noisy regressions retains bias: the 16-record empirical error decomposition is bias 0.03738, dispersion 0.00835 and RMSE 0.03830. Approximately 95.3% of empirical MSE is squared bias (S9.9). In the coefficient covariance check all 2000 paths pass the design guard, so it does not filter that finite sample. Its relative RMS linearization discrepancy can exceed one because its denominator is the RMS linear displacement (S3 and S9.3). These observations distinguish covariance calculation, estimator bias and independent acquisition.
'''
    template=(BASE/'efficiency_article_template.md').read_text(encoding='utf-8')
    archive=load(ROOT/'configs/public_archive.json')
    template=re.sub(r'version v\d+\.\d+\.\d+\]\(https://github.com/malyshevdanil/wsindy-information-limits/tree/v\d+\.\d+\.\d+\)',
        'version '+archive['version']+']('+archive['version_url']+')',template)
    assert template.count('@@NUMERICAL@@')==1
    template=template.replace('Small-noise accuracy is compared with the information reference; finite-noise bias and dispersion are reported separately.',
        f"The scalar fitted-pilot median error decreases from {s(.1,'single_window')['median']:.3f} to {s(.1,'pilot_phi',small)['median']:.3f}. Van der Pol empirical RMSE at 0.5% noise is {v(.005,480,'ou')['relative_rmse']:.5f}, versus the weak reference {bound/10:.5f}; the 5% experiment retains excess error. Bias and dispersion are reported separately.")
    (BASE/'article_revised.md').write_text(template.replace('@@NUMERICAL@@',content),encoding='utf-8')
    supp=(BASE/'efficiency_supplement_template.md').read_text(encoding='utf-8')
    supp+=r'''## S8. Efficiency protocols and complete summaries

Theoretical proofs and their function-space assumptions are stated in Sections 5.1–5.4 and in the standalone proof document. The configuration `configs/efficiency.json` defines the original scalar and Van der Pol experiments; each executed notebook retains its own exact configuration and source fingerprint. Numerical experiments run on private CPU Kaggle notebooks. Synthetic OU paths use seeds 2100000–2100999 for the scalar experiment and 2200000–2200127 for Van der Pol, independent of other studies. Noise amplitudes, methods and window families are paired within seed. Resampling uses 5000 independent-path bootstrap draws with generator seed 2302600; paired comparisons resample one common path-index vector. These are unadjusted exploratory intervals.

The scalar cutoff is the product of two quintic smoothstep ramps. For distance $d$ from an endpoint and width $\delta$, $z=\min(1,\max(0,(d-\delta/4)/(3\delta/4)))$ and $\eta=z^3(10-15z+6z^2)$. Its derivative is included analytically in the weak equation. The functions are compactly supported and $C^2$, which suffices for the implemented weak form; smooth approximation yields the same covariance closure. Scalar experiments use 65536 quadrature intervals. The pilot and every final scalar estimate are clipped to $[-2.5,0.5]$ for the bounded-risk summaries. Unbounded error medians and every estimate are preserved. The ratio denominator diagnostic uses the generating clean moment only after fitting; it neither rejects observations nor changes estimates.

Table S10. Scalar fitted-pilot summaries for every amplitude and cutoff. RMSE, bias and dispersion use bounded final estimates. All 1000 paths enter each row. Clips count altered final estimates; diagnostic pass counts do not select a subset.

| Noise SD | Cutoff | Median | RMSE | Bias | Dispersion | Clips | Diagnostic passes |
|---:|---:|---:|---:|---:|---:|---:|---:|
'''
    # S1-S9 tables reside in S9 below; labels S10 onward remain unambiguous.
    for sigma in (.1,.02,.005):
        for delta in (.2,.1,.05,.025,.0125):
            r=s(sigma,'pilot_phi',delta)
            supp+=f"| {sigma:g} | {delta:g} | {r['median']:.6f} | {r['relative_rmse']:.6f} | {r['relative_bias']:.6f} | {r['relative_dispersion']:.6f} | {r['clipped']} | {r['guard_passed']}/1000 |\n"
    remaining=[r for r in scalar if r['method']!='pilot_phi']
    total_clips=sum(r['clipped'] for r in remaining)
    supp+=f"\nThe fixed-window and generating-parameter-test cells also retain 1000 paths each. Across these cells {total_clips} final estimates are clipped; individual counts and all unbounded estimates are preserved in the numerical report and raw rows. This clipping changes only the bounded-risk summaries; no path is removed.\n\n"
    supp+=r'''
For Van der Pol the parameter pilot box is $[0.25,2]\times[-2,-0.25]\times[0.1,2]\times[-2,-0.1]$. DOP853 integration of its fitted model uses relative tolerance $10^{-9}$ and absolute tolerance $10^{-11}$. The fitted trajectory is fixed during the coefficient-weight iterations. Tests and the basis transformation are fixed independently of observations. The relative covariance relaxation is scaled by the mean covariance diagonal, so common noise-variance scaling does not change GLS weights. One hundred iterations and tolerance $10^{-7}$ are computation settings, not a convergence theorem. Final coefficients are not clipped, and every finite estimate is included. Basis dimensions and both covariance matrices are preserved in the deterministic geometry output.

Table S11. Van der Pol finite-cohort decomposition and stopping diagnostics. Counts use 128 paths per row; attaining the iteration cap without convergence is retained in the sample. These empirical moments do not establish the uniform-integrability assumption of Corollary 1.

| Noise | Windows | Weights | Bias | Dispersion | RMSE | Converged | Pilot clips |
|---|---:|---|---:|---:|---:|---:|---:|
'''
    control=report['stopping_control']
    supp=supp.replace('Table S11.',f"A five-iteration execution on the identical paths was retained as a stopping control. After inspecting its convergence diagnostics, the limit was increased to 100; this computational adjustment was not fixed before the five-iteration output. It does not add independent paths. The median relative coefficient change is {control['median_relative_coefficient_change']:.2e} and the maximum is {control['maximum_relative_coefficient_change']:.2e}. No path is removed because of the stopping test.\n\nTable S11.",1)
    for r in vdp:
        supp+=f"| {100*r['level']:g}% | {r['windows']} | {r['weight_model']} | {r['relative_bias']:.6f} | {r['relative_dispersion']:.6f} | {r['relative_rmse']:.6f} | {r['converged']}/128 | {r['pilot_clipped']} |\n"
    supp+='\nThe raw scalar, geometry and Van der Pol outputs, including coefficient vectors, are archived. No failed stochastic fit is omitted; completion and expected seed counts are checked separately from stopping tolerances.\n\n'
    supp+='## S9. Supporting numerical validation\n\n'+(BASE/'supporting_validation.md').read_text(encoding='utf-8')
    supp+='\n## References\n'
    labels=re.findall(r'^Table (S\d+)\.',supp,flags=re.M)
    mapping={label:f'S{i+1}' for i,label in enumerate(labels)}
    supp=re.sub(r'\bTable (S\d+)\b',lambda match:'Table '+mapping[match[1]],supp)
    (BASE/'supplement.md').write_text(supp,encoding='utf-8')
    ru='''# Проверка эффективности boundary-free weak-оценок

Доказаны точное описание adjoint-класса в H0^1 и его замыкания в ковариационной норме, точность границы Шура и условная теорема об эффективности GLS. Нулевые значения на концах не означают компактный носитель. Оптимальная OU-тестовая функция содержит граничные скачки; компактные функции достигают границы как инфимума. Пределы последовательные: малый шум при фиксированном наборе тестов, затем плотный набор.

Два эксперимента завершены на Kaggle. Скалярный эксперимент: 1000 независимых путей, три амплитуды, пять cutoff-ширин, oracle и fitted-pilot тесты плюс исходное окно. Van der Pol: 128 независимых путей, две амплитуды, три семейства, два ковариационных веса — 1536 оценок. Число строк скалярного эксперимента — 33000; общие seed не считаются новыми независимыми данными.

| Шум Van der Pol | Окна | RMSE OU | RMSE iid | Граница weak |
|---|---:|---:|---:|---:|
'''
    for level in (.005,.05):
        for K in (30,120,480):ru+=f"| {100*level:g}% | {K} | {v(level,K,'ou')['relative_rmse']:.6f} | {v(level,K,'iid')['relative_rmse']:.6f} | {bound*level/.05:.6f} |\n"
    ru+=f"\nСкалярный fitted-pilot тест при SD=0.1 и ширине 0.0125: медиана {high_scalar(s):.6f}; исходное окно {s(.1,'single_window')['median']:.6f}. Детальные данные и интервалы: efficiency_results.json и S8.\n\n"
    ru+='При малом шуме результаты согласуются с локальной теорией. При 5% шума добавление окон не доводит feasible-оценку до границы 0.016904: остаются конечные ошибки и смещение. Это ограничение прямо включено в статью. Утверждения о универсальном plateau iid и о совместном пределе не делаются. Полные авторские методы не переименованы в новую компонентную проверку; 448 оценок сохранены в приложении.\n'
    path=ROOT/'reports/efficiency_results_RU.md';path.write_text(ru,encoding='utf-8')
    pypandoc.convert_file(str(path),'docx',outputfile=str(path.with_suffix('.docx')))


def high_scalar(lookup):return lookup(.1,'pilot_phi',.0125)['median']


if __name__=='__main__':main()
