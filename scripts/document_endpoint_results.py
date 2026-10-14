"""Synchronize supplementary tables and protocol audit from complete reports."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    w=json.loads((ROOT/'reports/support_wald_results.json').read_text());assert w['completed']
    lines=['### S23.1. Complete rejection intervals and paired differences','',
        'Table S19. All targeted-test and sparse-inclusion frequencies. Each count uses one hundred common noise paths. Intervals for individual frequencies are Wilson 95%; difference intervals use 5000 paired bootstrap resamples. The nominal Wald size is 5%, with a failed null-calibration gate at 5% noise. MSTLS is a selection procedure without a nominal size parameter.','',
        '| Noise | Magnitude | Wald (Wilson 95%) | MSTLS (Wilson 95%) | Difference, percentage points (paired 95%) |',
        '|---|---:|---|---|---|']
    for level in [.005,.05]:
        for di,mult in enumerate(w['protocol']['multipliers']):
            group={g['method']:g for g in w['groups'] if g['level']==level and g['delta_index']==di}
            def rate(method):
                g=group[method];lo,hi=g['wilson_ci95'];return f"{g['included']}/100 ({100*lo:.1f}–{100*hi:.1f}%)"
            g=group['paired_difference'];lo,hi=g['difference_ci95']
            lines.append(f"| {100*level:g}% | {abs(mult):g} | {rate('weak_gls_wald')} | {rate('author_mstls')} | {100*g['wald_minus_mstls']:.0f} ({100*lo:.1f}–{100*hi:.1f}) |")
    p=ROOT/'manuscript/supplement.md';t=p.read_text(encoding='utf-8')
    if '### S23.1.' in t:t=t.split('### S23.1.',1)[0].rstrip()+'\n\n## References\n'
    assert t.count('## References')==1;t=t.replace('## References','\n'.join(lines)+'\n\n## References')
    p.write_text(t,encoding='utf-8')
    p=ROOT/'reports/boundary_protocol_decision.json';d=json.loads(p.read_text());d['endpoint_followup_completion']={
        'all_declared_primary_studies_complete':True,'ablation_paths':500,'vdp5_paths':128,'lorenz_paths':256,
        'wald_mstls_paired_paths_by_grid':1600,'full_paired_input_fingerprints_verified':True,
        'source_and_embedded_configuration_audit':'reports/boundary_ablation_results.json; reports/boundary_followup_status.json',
        'scientific_configurations':['configs/boundary_ablation.json','configs/boundary_generalization.json','configs/support_wald_v2.json'],
        'packaging_changes':'ablation v2 reduces embedded configs after source-size rejection; generalization uses shorter endpoint slugs after API errors; scientific protocols unchanged',
        'numerical_amendment':w['protocol']['numerical_amendment'],
        'outcome_scope':'null gate passes at 0.5%, fails at 5%; universal endpoint recovery passes Lorenz, fails VdP5; no outcome-based repairs',
        'legacy_queues_remain_held':True}
    p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
    p=ROOT/'reports/unstarted_followup_queue.json';q=json.loads(p.read_text())
    q['final_verified_primary_states']={r['stage']:r['status'] for r in json.loads((ROOT/'reports/boundary_followup_status.json').read_text())}
    assert all(s=='COMPLETE' for s in q['final_verified_primary_states'].values())
    q['all_declared_primary_studies_complete']=True;p.write_text(json.dumps(q,indent=2)+'\n',encoding='utf-8')
    lines=['# Граничная форма и отбор структуры: завершённый раунд','',
        'Все десять основных заданий Kaggle завершены. Источник, замороженные конфигурации и парные хеши проверены; пути не исключались.',
        '', '| Проверка | Результат | Вывод |','|---|---|---|',
        '| Абляция, 500 путей | Начальные полуокна: 98,9%; контроль Y(0)-u: 14,3%; терминальные: -0,3% | Главный механизм уже описан Brunel; контроль полезен отдельно |',
        '| Лоренц, 256 путей | RMSE 0,00040994 → 0,00021686; MLE 0,00020784 | Закрыто 95,5% [92,8%; 98,0%]; критерий пройден |',
        '| VdP, 5%, 128 путей | RMSE 0,0255965 → 0,0218790; MLE 0,0142088 | Закрыто 32,6% [18,3%; 45,8%]; критерий не пройден |',
        '| Wald, 0,5% | На одном информационном ориентире: 91/100 против 0/100 MSTLS | Согласуется с локальным прогнозом 90,8%; нулевой gate пройден |',
        '| Wald, 5% | Ложные срабатывания 23/100 [15,8%; 32,2%] | Калибровка не пройдена; преимущества по мощности заявлять нельзя |','',
        'Критерий калибровки проверяет совместимость Wilson-интервала с 5%, а не точность размера теста. 100 нулевых путей дают широкий интервал. Wald — один заранее выбранный коэффициент, не множественный отбор.',
        '', 'Теорема 3 описывает замыкание расширенного класса и полную информационную границу при известном начальном состоянии. Статья перестроена под цену, восстановление и отбор; материалы о независимых записях перенесены в приложение. Нативный основной PDF использует imsart[ejs].','',
        'Основные артефакты: manuscript/article_ejs.pdf, manuscript/article_revised.docx, manuscript/supplement.pdf. Полные парные оценки: data/endpoint_followups/; SHA256 manifest сохраняет исходные байты. Текущий публичный DOI ещё относится к предыдущей версии; новое архивирование и авторские декларации нужны перед подачей.']
    (ROOT/'reports/endpoint_followups_summary_RU.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':main()
