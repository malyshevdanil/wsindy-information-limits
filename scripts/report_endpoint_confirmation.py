"""Render the complete paired confirmation; never consume partial run rows."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    result = json.loads((ROOT / 'reports/endpoint_confirmation_results.json').read_text())
    assert result['completed'] and result['paired_input_fingerprints_verified']
    groups = {g['method']: g for g in result['groups']}
    assert all(g['paths'] == 1968 for g in groups.values())
    order = ['boundary_free_480', 'universal_boundary', 'model_adjoint', 'known_initial_scoring']
    labels = ['Boundary-free, 480', 'Universal + endpoints', 'Model-adjoint', 'Known-initial scoring']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    values = np.array([groups[m]['relative_rmse'] for m in order])
    ci = np.array([groups[m]['rmse_bootstrap_ci95'] for m in order])
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.errorbar(values * 1000, np.arange(4),
                xerr=np.maximum(0., np.vstack((values-ci[:, 0], ci[:, 1]-values))) * 1000,
                fmt='o', capsize=3, color='#24476b')
    ax.set_yticks(np.arange(4), labels)
    ax.invert_yaxis()
    ax.set_xlabel('Relative RMSE (10$^{-3}$); paired-bootstrap 95% intervals')
    ax.set_title('Van der Pol, known initial state, 0.5% OU noise; 1968 paired records')
    ax.grid(axis='x', alpha=.2)
    fig.tight_layout()
    destination = ROOT / 'reports/figures'
    destination.mkdir(exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(destination / f'endpoint_confirmation.{ext}', dpi=200)
    criteria = result['criteria']
    lines = ['# Подтверждающий граничный эксперимент', '',
             'Полная выборка: пути 33–2000 (1968 записей), Van der Pol, '
             'T=20, h=0.005, τ=0.3, шум 0.5%, известное начальное состояние.', '',
             'Все четыре ветки используют побитово совпадающие наблюдения. '
             'Проверены замороженные исходники, конфигурация и контроль квадратуры.', '',
             '| Метод | RMSE | 95% MC-интервал | Норма смещения | Разброс |',
             '|---|---:|---:|---:|---:|']
    for method, label in zip(order, labels):
        g = groups[method]
        low, high = g['rmse_bootstrap_ci95']
        lines.append(f"| {label} | {g['relative_rmse']:.8f} | [{low:.8f}; {high:.8f}] | "
                     f"{g['relative_bias']:.8f} | {g['relative_dispersion']:.8f} |")
    ratio = criteria['model_adjoint_ratio']
    ci_ratio = criteria['model_adjoint_ratio_ci95']
    closed = criteria['universal_fraction_gap_closed']
    lines += ['', f"Модельная сопряжённая ветка / scoring: {ratio:.5f}; "
              f"95% парный интервал [{ci_ratio[0]:.5f}; {ci_ratio[1]:.5f}]. "
              f"Критерий ≤1.05: {'выполнен' if criteria['model_adjoint_pass'] else 'не выполнен'}."]
    if closed is not None:
        low, high = criteria['gap_closed_ci95_conditional_positive_denominator']
        lines += ['', f"Универсальная ветка закрывает {closed*100:.2f}% RMSE-разрыва; "
                  f"95% парный интервал [{low*100:.2f}%; {high*100:.2f}%]. "
                  f"Критерий ≥50%: {'выполнен' if criteria['universal_pass'] else 'не выполнен'}."]
    else:
        lines += ['', 'Исходный RMSE-разрыв неположителен; критерий доли не оценивается.']
    lines += ['', 'Интервалы получены 5000 парными bootstrap-повторениями по записям. '
              'Первичное решение принимается по заранее замороженным точечным критериям.', '',
              'Модельная сопряжённая ветка реализует направления эффективного score в первом порядке. '
              'Её совпадение со scoring является проверкой реализации, а не новым преимуществом weak-метода.', '',
              'Универсальные окна строятся только по времени: 480 внутренних bump-функций и '
              '12 граничных полуокон, после H¹-ортогонализации — 491 независимый скалярный тест; '
              'добавлены два момента Y(0)−u. Ковариационные веса по-прежнему используют '
              'интегрирование пилотной ОДУ. Результат не доказывает полностью безмодельную или '
              'не требующую решения ОДУ оценку.', '',
              'Архивный likelihood-ориентир — один шаг Fisher scoring при известном u, '
              'а не полностью оптимизированная MLE. Архивные weak/scoring результаты на этих '
              'наблюдениях существовали до подтверждения; граничные ветки здесь оценены впервые.', '',
              'Абсолютный порог 0.00140 сохранён как историческое решение скрининга; '
              'его замена парными критериями зафиксирована до новых вычислений в '
              'boundary_protocol_decision.json. Неудача конечного универсального базиса '
              'ограничивает этот базис и эти веса, а не весь класс универсальных тестов.', '',
              'Время включает общий 30-оконный пилот; архивное время scoring отсутствует.']
    for method in ['model_adjoint', 'universal_boundary']:
        seconds = groups[method].get('median_online_seconds_with_pilot')
        if seconds is not None:
            lines.append(f"- {method}: медиана {seconds:.3f} с на запись.")
    (ROOT / 'reports/endpoint_confirmation_results_RU.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('Complete confirmation report and figure written')


if __name__ == '__main__':
    main()
