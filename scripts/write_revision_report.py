"""Write an honest Russian revision report from validated completed results."""
import json
from pathlib import Path
import pypandoc
ROOT=Path(__file__).resolve().parents[1]
report=json.loads((ROOT/'reports/revision_results.json').read_text())
info=report['information'];weak=info['empirical_median_absolute_error_ci']['weak'][-1];mle=info['empirical_median_absolute_error_ci']['mle'][-1]
lines=['# Результаты доработки по замечаниям','',
    'Целевой журнал — Applied Mathematics and Computation из исходного Word-плана. Подготовлена новая английская статья, supplement, native OMML Word, LaTeX, PDF, highlights, cover letter и поля авторских деклараций. Исходный план сохранён, старая статья архивирована.','',
    '## Математическое усиление','',
    'Выведена явная Cameron–Martin норма стационарного OU-процесса, включая граничные члены. Доказаны эквивалентность законов гладких ODE-траекторий, количественная двухточечная граница для любого оценивания и невозможность состоятельности одновременно во всём семействе. Это не запрещает состоятельность в одной изолированной точке. Классический принцип отнесён к литературе, а не заявлен новым.','',
    f"Для скалярного затухания непрерывная информация Фишера равна {info['continuous_information']:.3f}, эффективность локального Gaussian weak-суррогата — {info['weak_linearized_efficiency']:.3f}. Выведены точная информация на сетке, граница Cramér–Rao и различия между h, T и M. Для затухания увеличение T не даёт неограниченного роста информации; независимые записи дают сумму информации.",'',
    'Отдельно доказана потеря информации об исходном состоянии в boundary-free weak-влияниях. Интегрирование по частям даёт нулевой отклик на initial-state sensitivities и нижнюю границу covariance через Schur complement совместной информации. В скалярном примере она повышает нижнюю границу локальной дисперсии на 18.4% даже при известном исходном состоянии. Детерминированная проверка geometry на Kaggle завершена; она не добавляет Monte Carlo paths и не утверждает конечной дисперсии unguarded ratio.','',
    '## Новая численная проверка','',
    'На Kaggle завершено 24 576 свежих fits relative MSTLS/oracle-support: четыре системы, два режима шума, две амплитуды, четыре сетки, 128 paths на условие. Дополнительно: 1 000 scalar paths с exact Gaussian likelihood, 2 000 paths на каждый из четырёх уровней для joint coefficient covariance, шесть ковариационных ячеек по 250 affine-null bootstrap повторов, 400 блоков heterogeneous initial conditions и 1 536 fits ковариационного IRLS-компонента. Не все endpoints независимы: сетки/методы и некоторые условия специально используют common random numbers.','',
    f"На самой частой scalar-сетке медиана абсолютной ошибки weak ratio — {weak[0]:.4f}, bounded Gaussian likelihood — {mle[0]:.4f}. Это сравнение не превращает likelihood MSE в несмещённую CR-efficiency и не приписывает unguarded ratio конечную дисперсию.",'',
    '| Уровень шума | Guarded paths | Расхождение actual covariance | Paired linearization RMS |','|---|---:|---:|---:|']
for r in report['coefficient']:lines.append(f"| {r['level']:.1%} | {r['guarded_paths']}/2000 | {r['guarded_actual_relative_frobenius']:.3f} | {r['paired_linearization_relative_rms']:.3f} |")
lines+=['','Actual covariance условна на singular-value guard. Это не доказательство существования безусловных моментов оценивания.','',
    'В Van der Pol при 5% OU-шуме исправленный MSTLS восстановил точную структуру в 81 из 128 случаев. При правильной структуре медианная ошибка — 0.0340, при неправильной — 0.7189; гистограмма показывает, почему один общий показатель скрывает разные исходы. Для corrected residual при 50% шуме расхождение covariance 0.064 сопоставимо с fourth-moment RMS floor 0.063 и approximate-null 95% range [0.053, 0.072].','',
    '## Что исправлено в представлении','',
    'Аннотация переписана вокруг информационного предела, а не числа runs. Неподтверждённая grid-диагностика убрана из abstract/contributions и сохранена как отрицательный результат в supplement. Добавлены 19 проверенных DOI, fixed-domain/EiV контекст, нумерованные ссылки и литературные положения внутри аргументов. Новые графики имеют человеческие подписи, интервалы и support-conditional распределения; таблицы разделяют coefficient error и support recovery. Формулы Word являются native OMML, raw TeX проверяется автоматически.','',
    '## Оставшиеся границы готовности','',
    'Опубликованный relative MSTLS регрессионный блок сверён с авторским кодом, но полный adaptive WSINDy не воспроизведён. IRLS-проверка относится к шагу Algorithm 2 WENDy и его OU-расширению; окна, остановка и регуляризация отличаются. Полный WmSINDy с совместным noise reconstruction не запущен. Поэтому статья не ранжирует эти pipelines и не утверждает превосходство над ними.','',
    'Перед подачей нужны реальные авторы/аффилиации, декларации и человеческая проверка математики, кода и выводов. Квартиль проверяется по конкретной базе, году и категории; подготовка под журнал не гарантирует Q1-принятие. В журнал ничего не отправлено, открытый code/data DOI не опубликован.']
target=ROOT/'reports/q1_revision_results_RU.md';target.write_text('\n'.join(lines),encoding='utf-8')
pypandoc.convert_file(str(target),'docx',format='markdown+tex_math_dollars',outputfile=str(target.with_suffix('.docx')))
print('Russian revision report built')
