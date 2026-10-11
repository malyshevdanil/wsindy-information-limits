"""Postprocess completed reviewer calculations and write publication sections."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pypandoc
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/kaggle_reviewer'


def load(p):return json.loads(p.read_text(encoding='utf-8'))


def decomposition(estimates,truth):
    errors=(np.asarray(estimates)-truth).reshape(len(estimates),-1)/np.linalg.norm(truth)
    mean=errors.mean(axis=0);bias2=float(mean@mean)
    variance=float(np.mean(np.sum((errors-mean)**2,axis=1)))
    mse=float(np.mean(np.sum(errors**2,axis=1)))
    assert abs(mse-bias2-variance)<1e-12
    adjusted=max(0.,bias2-variance/(len(errors)-1))
    return {'n':len(errors),'relative_bias':float(np.sqrt(bias2)),
            'relative_dispersion':float(np.sqrt(variance)),'relative_rmse':float(np.sqrt(mse)),
            'empirical_squared_bias_fraction':bias2/mse,
            'mean_sampling_adjusted_squared_bias_fraction':adjusted/mse,
            'normalized_mean_error':mean.tolist()}


def replace_section(text,begin,end,section,anchor):
    if begin in text:
        left=text.index(begin);right=text.index(end,left)+len(end)
        return text[:left]+begin+'\n\n'+section+'\n\n'+end+text[right:]
    return text.replace(anchor,begin+'\n\n'+section+'\n\n'+end+'\n\n'+anchor)


def main():
    assert load(OUT/'completion.json')['status']=='complete'
    info=load(OUT/'vdp_information.json');windows=load(OUT/'window_families.json')
    replay=load(OUT/'heterogeneous_coefficients.json');truth=np.array(replay['truth'])
    archived=load(ROOT/'results/kaggle_revision_fast/heterogeneous_rows.json')
    old={(r['block'],r['M'],r['method']):r['error'] for r in archived}
    assert len(replay['rows'])==len(old)==3600
    difference=max(abs(r['error']-old[(r['block'],r['M'],r['method'])]) for r in replay['rows'])
    assert difference<1e-8
    groups=[]
    for method in ('stacked','mean_weak','naive_mean_sensor'):
        for M in (1,4,16):
            rows=[r for r in replay['rows'] if (r['method'],r['M'])==(method,M)]
            groups.append({'method':method,'M':M,**decomposition([r['coefficients'] for r in rows],truth)})
    author_rows=[]
    for stage in ('adaptive_classical',*[f'adaptive_wm_{i}' for i in range(4)]):
        author_rows+=load(ROOT/f'results/kaggle_{stage}/rows.json')
    fine_input=next((ROOT/'results/kaggle_adaptive_classical/inputs').glob('*h0.01000.npz'))
    with np.load(fine_input) as z:
        full_truth=z['truth'];sigma=.05*np.std(z['clean'],axis=0)
    assert np.allclose(sigma,info['sigma_fixed'],rtol=1e-10,atol=1e-12)
    methods=[];rng=np.random.default_rng(1926000)
    for method in ('wsindy','wmsindy','wendy','fixed_hermite_known_support'):
        rows=[r for r in author_rows if (r['method'],r['regime'],r['h'])==(method,'physical_ou',.01)]
        assert len(rows)==16 and all(r['status']=='ok' for r in rows)
        result=decomposition([r['coefficients'] for r in rows],full_truth)
        squares=np.array([r['error']**2 for r in rows])
        boot=np.sqrt(np.mean(squares[rng.integers(0,len(rows),(5000,len(rows)))],axis=1))
        methods.append({'method':method,**result,'relative_rmse_bootstrap_95':np.quantile(boot,[.025,.975]).tolist()})
    geometry=load(ROOT/'results/kaggle_geometry/geometry.json')
    I=geometry['coefficient_information_known_initial'];Ie=geometry['coefficient_information_initial_nuisance']
    opt=[]
    for row in windows['results']:
        opt.append({**row,'known_initial_efficiency':row['optimal_linear_information']/I,
                    'weak_class_efficiency':row['optimal_linear_information']/Ie})
        assert row['optimal_linear_information']<=Ie*(1+1e-6)
    report={'completion':load(OUT/'completion.json'),'manifest':load(OUT/'manifest.json'),
        'output_archive_sha256':hashlib.sha256((ROOT/'results/kaggle_reviewer_download/wsindy_reviewer_outputs.zip').read_bytes()).hexdigest(),
        'vdp_information':info,'optimal_window_families':opt,'heterogeneous_bias_variance':groups,
        'replay_maximum_error_difference':difference,'adaptive_ou_finest_bias_variance':methods,
        'scope':'All numerical experiments execute on Kaggle. Replay adds coefficient preservation, not independent evidence. Information and empirical RMSE are comparable scales, not proof of efficiency of biased estimators.'}
    (ROOT/'reports/reviewer_results.json').write_text(json.dumps(report,indent=2))
    known=info['grids'][-1]['continuous_quadrature']['relative_rmse_bound']
    weak=info['grids'][-1]['continuous_quadrature']['weak_initial_annihilation_relative_rmse_bound']
    mu=info['grids'][-1]['continuous_quadrature']['mu_only_standard_deviation_bound']
    section=r'''### 7.7. Van der Pol information reference

The OU experiment in Section 7.6 has $x'=a v$ and $v'=b x+c v+d x^2v$, with four independently free active coefficients $(a,b,c,d)=(1,-1,1,-1)$. The initial state (2,0) is known. We hold the calibrated coordinate noise variances fixed when differentiating the mean, rather than re-estimating state-scaled variance at perturbed parameters. The calibration is exactly that of the benchmark fine grid. Coefficient and initial-state sensitivities solve $S'=JS+\partial_\theta f$ and $Z'=JZ$, with $S(0)=0$ and $Z(0)=I$. OU innovations give discrete information; the continuous norm in Section 4.1 gives its limit.

For four active coefficients, the relative RMS reference is $\{\operatorname{tr}(\mathcal I_{\theta\theta}^{-1})\}^{1/2}/\|\theta\|_2$. The weak initial-state-annihilation reference replaces information by its Schur complement. Both are local covariance bounds, rather than bounds on error medians or on every biased estimator.

Table 7. Van der Pol OU relative RMS information references for the exact Section 7.6 calibration. Four active coefficients are independently free. The second column assumes known initial state; the third also imposes weak initial-sensitivity annihilation.

| Observation | Known-initial reference | Weak influence reference |
|---|---:|---:|'''
    for h in (.04,.01):
        row=next(r for r in info['grids'] if r['h']==h)['discrete']
        section+=f"\n| Grid h={h} | {row['relative_rmse_bound']:.6f} | {row['weak_initial_annihilation_relative_rmse_bound']:.6f} |"
    section+=f'\n| Continuous record | {known:.6f} | {weak:.6f} |\n\n'
    wms=next(r for r in methods if r['method']=='wmsindy');wendy=next(r for r in methods if r['method']=='wendy')
    section+=f"For the same 16 OU records, WmSINDy has relative RMSE {wms['relative_rmse']:.5f} (descriptive 95% bootstrap interval [{wms['relative_rmse_bootstrap_95'][0]:.5f}, {wms['relative_rmse_bootstrap_95'][1]:.5f}]) and WENDy {wendy['relative_rmse']:.5f} ([{wendy['relative_rmse_bootstrap_95'][0]:.5f}, {wendy['relative_rmse_bootstrap_95'][1]:.5f}]). Their RMSEs are respectively {wms['relative_rmse']/known:.2f} and {wendy['relative_rmse']/known:.2f} times the known-initial local reference. This scale comparison is more informative than comparing their medians near 0.015 with a variance bound. It is not an efficiency certificate: sparse selection and nonlinear bias are present, the sample is small, and the regular locally unbiased assumption is unverified. WmSINDy also includes trajectory-prediction terms, so it is not restricted to boundary-free weak influences.\n\n"
    section+=f"A one-parameter Van der Pol model fixes a=1 and b=-1 and ties c=mu, d=-mu. Its continuous standard-deviation reference for mu is {mu:.6f}. That substantially stronger constraint changes the estimation problem and must not be used to judge four-coefficient discovery. Central finite differences check the four sensitivities with relative discrepancies below {max(info['finite_difference_relative_errors']):.1e}; refining information quadrature from h=0.0025 to 0.00125 checks the reported reference.\n\n"
    section+='### 7.8. Optimal efficiency within declared window spans\n\n'
    section+=r'The scalar parameters are those of Section 7.1. We optimize the covariance-weighted linear combination in Section 4.3, including every cross-window covariance. Family F30 contains 30 fourth-power bump windows of full width 1.6; centers range uniformly from 0.805 to 3.195. Family F120 contains 60 centers for each of widths 1.6 and 3.0, with centers 0.005 inside the admissible support range. These are declared finite families, rather than an exhaustive search over all shapes and widths. The reference single window has width three and center two.'+'\n\n'
    section+='Table 8. Oracle optimal local influence variance and efficiency within each scalar window span. Information references are 144.187 for known initial state and 121.752 after initial-state projection. The table uses 8192 quadrature intervals and relative covariance-eigenvalue cutoff 1e-12.\n\n| Family | Retained rank | Variance | Known-initial efficiency | Weak-class efficiency |\n|---|---:|---:|---:|---:|'
    names={'reference_single':'Single window','fixed_30':'F30 (30 windows)','enriched_60':'F120 (120 windows)'}
    for r in opt:
        if r['grid_intervals']==8192 and r['eigenvalue_cutoff']==1e-12:
            section+=f"\n| {names[r['family']]} | {r['rank']} | {r['variance']:.6f} | {r['known_initial_efficiency']:.4f} | {r['weak_class_efficiency']:.4f} |"
    change=max(abs(next(q['variance'] for q in opt if (q['family'],q['grid_intervals'],q['eigenvalue_cutoff'])==(r['family'],4096,r['eigenvalue_cutoff']))/r['variance']-1) for r in opt if r['grid_intervals']==8192)
    conservative=next(r for r in opt if (r['family'],r['grid_intervals'],r['eigenvalue_cutoff'])==('enriched_60',8192,1e-10))
    section+=f"\n\nThus 0.063 describes one window, not the optimal efficiency of the weak family. F120 attains {next(r['known_initial_efficiency'] for r in opt if (r['family'],r['grid_intervals'],r['eigenvalue_cutoff'])==('enriched_60',8192,1e-12)):.3f} of known-initial information in the retained span. Tightening the cutoff from 1e-10 to 1e-12 raises retained rank from {conservative['rank']} to 68; at 1e-10 its known-initial efficiency is {conservative['known_initial_efficiency']:.3f}. This sensitivity is disclosed rather than treated as an exact infinite-dimensional optimum. Doubling quadrature intervals changes every reported variance by less than {100*change:.4f}%. All optimal influences satisfy coefficient normalization and initial-state annihilation numerically. Their efficiencies remain below the initial-state-projected ceiling.\n\n"
    section+='### 7.9. Bias and dispersion of aggregation\n\n'
    section+=r'The coefficient-preserving execution replays exactly the 400 independent blocks in Table 3; it adds no independent samples. For normalized coefficient errors, the empirical identity states that squared RMSE equals squared bias plus squared dispersion. Table 9 reports the norm of the sample mean error as bias and the square root of the second term as dispersion. Bias here is an empirical estimate, not a proved population limit. We also subtract sample-mean sampling variance when estimating the squared-bias fraction; the adjustment is small for 400 blocks.'+'\n\n'
    section+='Table 9. Relative coefficient bias, dispersion and RMS error for heterogeneous Van der Pol aggregation; the same 400 blocks as Table 3. Values are dimensionless. Bias fraction uses the unadjusted empirical squared mean error.\n\n| Method | M | Bias | Dispersion | RMSE | Squared-bias fraction |\n|---|---:|---:|---:|---:|---:|'
    labels={'stacked':'Stacked design','mean_weak':'Mean weak matrices','naive_mean_sensor':'Mean sensor path'}
    for r in groups:
        if r['M'] in (1,16):
            section+=f"\n| {labels[r['method']]} | {r['M']} | {r['relative_bias']:.5f} | {r['relative_dispersion']:.5f} | {r['relative_rmse']:.5f} | {r['empirical_squared_bias_fraction']:.1%} |"
    single=next(r for r in groups if (r['method'],r['M'])==('stacked',1));stack=next(r for r in groups if (r['method'],r['M'])==('stacked',16))
    section+=f"\n\nFor stacking, dispersion decreases from {single['relative_dispersion']:.5f} at M=1 to {stack['relative_dispersion']:.5f} at M=16, while bias changes from {single['relative_bias']:.5f} to {stack['relative_bias']:.5f}. The estimated bias contribution is {100*single['empirical_squared_bias_fraction']:.1f}% of single-record MSE and {100*stack['empirical_squared_bias_fraction']:.1f}% at M=16; the sample-mean-adjusted fractions are {100*single['mean_sampling_adjusted_squared_bias_fraction']:.1f}% and {100*stack['mean_sampling_adjusted_squared_bias_fraction']:.1f}%. The stacking floor is therefore dominated by systematic random-design error as acquisition increases. Whether the single-record error is mainly bias must be judged by this decomposition, not by the two median norms. Mean weak matrices reduce both errors under the stated excitation, whereas sensor-path averaging produces nonlinear averaging bias.\n"
    path=ROOT/'manuscript/article_revised.md';text=path.read_text(encoding='utf-8')
    text=replace_section(text,'<!-- REVIEWER_RESULTS_START -->','<!-- REVIEWER_RESULTS_END -->',section,'## 8. Discussion and limitations')
    phrase='The evidence supports a distinction between aggregation operations, under the supplied support and excitation used here.'
    cross_reference=' Coefficient-level bias and dispersion are reported in Section 7.9; the near-equal single-record and stacked median errors alone do not identify the bias fraction.'
    text=text.replace(cross_reference,'')
    text=text.replace(phrase,phrase+cross_reference)
    if 'Efficient weak-form estimation under correlated measurement error' not in text:
        path.write_text(text,encoding='utf-8')
    ru=['# Дополнительные расчёты по замечаниям к статье','',f'Все расчёты выполнены на Kaggle. Повтор 400 блоков совпадает с прежними нормами ошибок: максимальное отклонение {difference:.3e}; это не новая независимая выборка.','',
        f'Информационная граница Van der Pol: относительная RMS-шкала {known:.6f} при известном начальном состоянии и {weak:.6f} для слабых влияний с аннулированием чувствительности к начальному состоянию.',
        f'WmSINDy: RMSE {wms["relative_rmse"]:.6f}; WENDy: {wendy["relative_rmse"]:.6f}. Это сравнение шкал, а не доказательство эффективности смещённых оценок.','',
        '## Оптимальные комбинации окон','', '| Семейство | Эффективность относительно полного пути | Эффективность относительно слабого класса |','|---|---:|---:|']
    for r in opt:
        if r['grid_intervals']==8192 and r['eigenvalue_cutoff']==1e-12:
            ru.append(f"| {names[r['family']]} | {r['known_initial_efficiency']:.5f} | {r['weak_class_efficiency']:.5f} |")
    ru+=['','## Смещение и разброс','', '| Метод | M | Смещение | Разброс | RMSE | Доля смещения в MSE |','|---|---:|---:|---:|---:|---:|']
    for r in groups:ru.append(f"| {labels[r['method']]} | {r['M']} | {r['relative_bias']:.6f} | {r['relative_dispersion']:.6f} | {r['relative_rmse']:.6f} | {r['empirical_squared_bias_fraction']:.1%} |")
    path=ROOT/'reports/reviewer_results_RU.md';path.write_text('\n'.join(ru)+'\n',encoding='utf-8')
    pypandoc.convert_file(str(path),'docx',format='markdown',outputfile=str(path.with_suffix('.docx')))
    print(json.dumps({'known_initial_cr_rmse':known,'weak_class_cr_rmse':weak,
        'wm_rmse':wms['relative_rmse'],'wendy_rmse':wendy['relative_rmse'],
        'single_bias_fraction':single['empirical_squared_bias_fraction'],'stack16_bias_fraction':stack['empirical_squared_bias_fraction'],
        'replay_max_error_difference':difference}))


if __name__=='__main__':main()
