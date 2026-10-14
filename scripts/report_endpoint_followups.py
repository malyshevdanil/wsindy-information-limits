"""Render only complete, audited endpoint/selection summaries; no simulations."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'manuscript/figures/publication'


def general_report(system):
    path=ROOT/f'reports/boundary_generalization_{system}_results.json'
    if not path.exists():return
    r=json.loads(path.read_text());assert r['completed']
    lines=[f'# Универсальная граничная форма: {system}',
           '\nВсе запланированные пути сохранены; исходники, конфигурация и парные хеши проверены.',
           '\n| Метод | RMSE | Норма смещения | Разброс |', '|---|---:|---:|---:|']
    for g in r['groups']:
        lines.append(f"| {g['method']} | {g['relative_rmse']:.8f} | {g['relative_bias']:.8f} | {g['relative_dispersion']:.8f} |")
    lo,hi=r['fraction_ci95']
    lines.extend([f"\nЗакрыто {100*r['fraction_gap_closed']:.1f}% парного разрыва; 95% интервал [{100*lo:.1f}%; {100*hi:.1f}%].",
                  f"Замороженный критерий ≥50%: {'пройден' if r['pass_frozen_point_rule'] else 'не пройден'}.",
                  '\nОриентиры первого порядка и полной эрмитовой ковариации относятся к фиксированному чистому дизайну. Их разность не является точным разложением дисперсии реализуемой оценки.'])
    (ROOT/f'reports/boundary_generalization_{system}_results_RU.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def ablation_figure():
    r=json.loads((ROOT/'reports/boundary_ablation_results.json').read_text())
    # Use stored summaries only, without reevaluating the frozen bootstrap.
    groups=r['groups']
    lookup={g['method']:g for g in groups}
    names=['initial_only','terminal_only','control_only','universal_boundary']
    labels=['Initial traces','Terminal traces','Initial-noise controls','All components']
    vals=[];intervals=[]
    for name in names:
        g=lookup[name];vals.append(100*g['fraction_matched_gap_closed']);intervals.append(100*np.array(g['fraction_closed_ci95']))
    fig,ax=plt.subplots(figsize=(6.8,3.2))
    vals=np.array(vals);intervals=np.array(intervals)
    ax.errorbar(vals,np.arange(4),xerr=np.array([vals-intervals[:,0],intervals[:,1]-vals]),fmt='o',capsize=4,color='#245b8a')
    ax.set_yticks(np.arange(4),labels);ax.invert_yaxis();ax.axvline(0,color='.7',lw=.8)
    ax.axvline(50,color='.6',ls='--',lw=.8);ax.set_xlabel('Matched weak-to-scoring RMSE gap recovered (%)')
    ax.grid(axis='x',alpha=.2);fig.tight_layout();FIG.mkdir(parents=True,exist_ok=True)
    for suffix in ['png','pdf']:fig.savefig(FIG/f'endpoint_ablation.{suffix}',dpi=200)
    plt.close(fig)


def wald_report():
    path=ROOT/'reports/support_wald_results.json'
    if not path.exists():return
    r=json.loads(path.read_text());assert r['completed'];groups=r['groups']
    lines=['# Wald-отбор и авторский MSTLS: парное сравнение',
           '\nЗаранее выбран один коэффициент x³ во втором уравнении. Это не процедура множественного отбора библиотеки и не вывод после отбора.',
           '\n| Шум | Метод | Нулевое включение | Wilson 95% | Калибровка α=5% |', '|---|---|---:|---|---|']
    for g in r['calibration']:
        lo,hi=g['wilson_ci95'];gate='не применимо' if g['nominal_alpha'] is None else ('да' if g['calibration_gate_pass'] else 'нет')
        lines.append(f"| {100*g['level']:g}% | {g['method']} | {100*g['null_rate']:.0f}% | [{100*lo:.1f}%; {100*hi:.1f}%] | {gate} |")
    lines.extend(['\n| Шум | Множитель локального двухточечного ориентира | Wald | MSTLS | Парная разность | 95% интервал разности |','|---|---:|---:|---:|---:|---|'])
    fig,axes=plt.subplots(1,2,figsize=(8.2,3.5),sharey=True)
    for ax,level in zip(axes,[.005,.05]):
        for method,label,color in [('weak_gls_wald','Wald, nominal 5%','#245b8a'),('author_mstls','Author MSTLS','#bb5638')]:
            sample=[g for g in groups if g['level']==level and g['method']==method]
            xx=np.abs([g['multiplier'] for g in sample]);yy=np.array([g['rate'] for g in sample]);ci=np.array([g['wilson_ci95'] for g in sample])
            ax.errorbar(xx,yy,yerr=np.maximum(0,np.array([yy-ci[:,0],ci[:,1]-yy])),marker='o',capsize=2,label=label,color=color,lw=1)
        diff=[g for g in groups if g['level']==level and g['method']=='paired_difference']
        ax.plot(np.abs([r['protocol']['multipliers'][g['delta_index']] for g in diff]),[g['null_local_gaussian_wald_power_reference'] for g in diff],ls='--',color='.35',label='Null-local Gaussian Wald reference')
        for g in diff:
            w=next(v for v in groups if v['level']==level and v['delta_index']==g['delta_index'] and v['method']=='weak_gls_wald')
            m=next(v for v in groups if v['level']==level and v['delta_index']==g['delta_index'] and v['method']=='author_mstls')
            lo,hi=g['difference_ci95'];lines.append(f"| {100*level:g}% | {abs(w['multiplier']):g} | {100*w['rate']:.0f}% | {100*m['rate']:.0f}% | {100*g['wald_minus_mstls']:.0f} п.п. | [{100*lo:.0f}; {100*hi:.0f}] п.п. |")
        cal=next(g for g in r['calibration'] if g['level']==level and g['method']=='weak_gls_wald')
        ax.set_title(f"{100*level:g}% noise (null gate {'passed' if cal['calibration_gate_pass'] else 'failed'})")
        ax.set_xlabel('Magnitude / two-point information reference');ax.set_ylim(-.03,1.04);ax.grid(alpha=.2)
    axes[0].set_ylabel('Rejection / inclusion frequency');axes[1].legend(fontsize=7,loc='lower right');fig.tight_layout()
    for suffix in ['png','pdf']:fig.savefig(FIG/f'support_wald_mstls.{suffix}',dpi=200)
    plt.close(fig)
    lines.extend(['\nЕсли нулевой Wilson-интервал Wald не содержит 5%, сравнение мощности не подтверждает калиброванную эффективность.',
                  '\nТочная эрмитова ковариация невязок в пилотной точке не является точной конечномерной дисперсией оценки со случайным дизайном. Двухточечный информационный ориентир и двусторонний Wald-тест используют разные критерии.',
                  '\nВерсия v2 исправляет детерминированный дефект квадратуры по пилотной модели. Порог чистого контроля, уровни значимости, данные и настройки MSTLS сохранены. Старые версии архивированы.'])
    (ROOT/'reports/support_wald_results_RU.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    for system in ['vdp5','lorenz']:general_report(system)
    ablation_figure();wald_report()
