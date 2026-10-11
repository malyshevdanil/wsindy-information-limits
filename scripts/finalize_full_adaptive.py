"""Update method descriptions only after full benchmark validation and analysis."""
import json
from pathlib import Path
import pypandoc
from docx import Document
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'


def main():
    report=json.loads((ROOT/'reports/full_adaptive_results.json').read_text())
    if report['fit_count']!=448:raise ValueError('Incomplete full-pipeline benchmark')
    path=BASE/'article_revised.md';text=path.read_text(encoding='utf-8')
    old='WmSINDy combines weak regression with simultaneous dynamics and noise reconstruction [@lopez2026]. Its adaptive reconstruction is more extensive than the relative-magnitude MSTLS regression component implemented for the independent validation. The latter compares regression on fixed weak windows rather than complete discovery pipelines.'
    new='WmSINDy combines weak regression with simultaneous dynamics and noise reconstruction [@lopez2026]. We distinguish the large fixed-window regression-component study from a separate paired execution of the complete adaptive author WSINDy, WENDy and WmSINDy pipelines. The latter retains noise estimation or reconstruction and adaptive test-function selection; WENDy supplies structure for parameter estimation.'
    text=text.replace(old,new)
    old='The sparse comparison concerns a published regression component with fixed windows. Complete adaptive WSINDy, WENDy and WmSINDy pipelines are outside the executed benchmark, so the study cannot rank them. Oracle states and supplied support are verification devices.'
    new='The large sparse comparison concerns a published regression component with fixed windows. The additional full adaptive author-pipeline benchmark uses one oscillator, two meshes and 16 paths per condition. It establishes executable comparisons under those conditions, with substantial support uncertainty, but cannot establish universal method rankings or an asymptotic plateau. WENDy and its fixed-window parameter comparators receive supplied structure; discovery results are reported separately. Oracle states and supplied support are verification devices.'
    text=text.replace(old,new)
    text=text.replace('Numerical comparisons examine full Gaussian likelihood, sparse weak regression, coefficient covariance and independent recordings.',
                      'Numerical comparisons examine full Gaussian likelihood, sparse weak regression, complete adaptive author pipelines, coefficient covariance and independent recordings.')
    path.write_text(text,encoding='utf-8')
    supplement=BASE/'supplement.md';text=supplement.read_text(encoding='utf-8')
    section='''## S6. Complete adaptive author-pipeline reproduction

A separate protocol, `configs/full_adaptive.json`, was frozen before comparative outputs. It uses Van der Pol with parameter one, initial state (2,0), duration 20, noise amplitude five percent, OU correlation time 0.3 or independent errors, and meshes 0.04 and 0.01. Seeds 910000--910015 give 16 independent records per cell. All algorithms receive the same fine-grid error realization and its coarse restriction; SHA256 checks verify paired inputs. This seed range is separate from the earlier studies. Discovery libraries have nine nonconstant monomials of total degree one through three, matching the released WmSINDy example. Fixed-window discovery and supplied-support comparators are rerun on that library, rather than reusing ten-feature results.

Pinned author repositories are WSINDy_ODE commit `72d8d6e5be2b80bcfe5015dced399ace1ceca43f`, WENDy commit `65c52a893b1f693cae3cd17e21854665ea814b9d`, and the WmSINDy Van der Pol example at commit `f1a781009a5c6adfd5650bec63ef30cc0673ae97` [@messenger2021; @bortz2023; @lopez2026]. Files or repository ZIPs are fetched from those revisions and verified with the SHA256 values in `configs/author_pins.json`. Core author functions are executed without rewriting the statistical algorithms. Author sources are excluded from distributed archives; download instructions and our adapters are provided.

WSINDy uses FFT-corner widths and polynomial degrees, gradient-adapted placement with relax_AG=1, automatic SMAF smoothing initialized as in the released driver, scaling parameter two, 100 thresholds from 1e-4 to one, zero ridge penalty and loss weight 0.8. The adaptive placement is explicitly enabled relative to the released driver's uniform-placement default. GLS remains disabled, as in that driver. Selected widths, degrees, centers, smoothing-filter lengths and losses are saved.

WENDy uses the released set-and-forget settings: estimated noise, the exp(-9/(1-x²)) bump, automatic minimum radius and four dyadic multiples, SVD corner truncation, maximum 100 iterations, relative tolerance 1e-10, diagonal relaxation 1e-10, and normality checks after iteration ten with cutoff 1e-4. Its four-term structure is prespecified. Actual coefficients are withheld: the diagnostic true_vec is replaced by a unit vector and use_true=0. The returned normality value and iterates are archived. These are the released iid-noise equations; under physical OU error they are evaluated as a misspecified observation model, without claiming a full adaptive OU-WENDy extension.

The MATLAB functions execute under GNU Octave 8.4.0 with statistics and signal packages. Separate compatibility functions provide exact analytic polynomial Jacobians in the declared feature order, exact values and first derivatives of the specified bump with author normalization, and the two-input pagemtimes operation used in covariance construction. These replace Symbolic Toolbox and missing host operations. They do not change radius selection, SVD, covariance updates, regression or stopping logic. Bitwise agreement with a proprietary MATLAB runtime is not asserted.

WmSINDy executes the author Python initializer, RK4 predictions, weak-residual loss, gradient updates, adaptive window function, pruning and least-squares refits. Every completed fit performs eight outer cycles of 5000 Adam steps (40,000 total), with threshold 0.15, learning rate 0.001, epsilon 1e-7, q=1 and decay 0.9. Noise is initialized to zero. Legacy Adam is used with TensorFlow 2.20.0 / tf-keras 2.20.1. Slots for new coefficient variables are created outside traced functions, and the author loss functions are retraced after each window update so changed variable-container shapes are actually used. Adam is initialized independently for each recording; noise-variable slots and iteration count continue through its eight cycles. The released all-valid-convolution-row calculation is retained. Loss checkpoints, eight coefficient iterates, window widths, spectral corners, optimizer counts and the reconstructed noise signal are archived.

The compatibility control uses the author's Van der Pol example: parameter 0.5, initial state (-2,1), duration ten, 1000 linspace points, iid Gaussian amplitude 25 percent and seed one. Its nominal WmSINDy prediction step 0.01 is retained although the linspace spacing is 10/999. The benchmark instead uses its exact grid spacing. Clean trajectories are integrated with DOP853 at relative and absolute tolerances 1e-12 for this control. The single control checks execution and recovery and does not reproduce the published 50-run aggregate figure. Reported comparison coefficient errors use the Frobenius norm.

The comparison contains 448 planned fits: four discovery algorithms and three supplied-support parameter estimators, four observation cells and 16 records. Discovery success is exact support equality; supplied-support estimators receive no structure-discovery score. Failure counts remain in success denominators. Error medians and bootstrap intervals condition on completed fits; all failures and their logs are retained. Mesh differences are paired by recording. Wilson limits accompany support counts in `reports/full_adaptive_results.json`. Computational runtimes are execution diagnostics, not cross-language speed claims. With one oscillator and 16 paths per cell, these results do not justify universal rankings or fine precision of support probabilities.

'''
    if '## S6. Complete adaptive author-pipeline reproduction' not in text:text=text.replace('## References',section+'## References')
    for pin in json.loads((ROOT/'configs/author_pins.json').read_text()).values():
        text=text.replace('`'+pin['commit']+'`','`'+pin['commit'][:12]+'`')
    text=text.replace('verified with the SHA256 values in `configs/author_pins.json`.',
                      'verified with the full revision identifiers and SHA256 values in `configs/author_pins.json`.')
    supplement.write_text(text,encoding='utf-8')
    lines=['# Результаты полного адаптивного воспроизведения','',
           'На Kaggle завершены полные авторские WSINDy, WENDy и WmSINDy с закреплёнными версиями и проверенными SHA256. WmSINDy прошёл 8 × 5000 шагов в каждом успешном fit, включая обновление окон, восстановление шума и повторную регрессию.','',
           'Сравнение: Ван-дер-Поль, T=20, шум 5%, OU с τ=0.3 и независимый шум, h=0.04/0.01, 16 независимых записей в каждой ячейке. Всего 448 fits; наблюдения между методами совпадают побитово. WENDy и два параметрических сравнения получают заранее заданную структуру; их результаты не оценивают обнаружение структуры.','',
           '| Шум | h | Метод | Медианная ошибка | 95% интервал медианы | Структура | Сбои |',
           '|---|---:|---|---:|---|---:|---:|']
    for g in report['groups']:
        support='задана' if g['support_count'] is None else f"{g['support_count']}/{g['n']}"
        median='n/a' if g['median'] is None else f"{g['median']:.5f}"
        ci='n/a' if g['median_ci'] is None else f"[{g['median_ci'][0]:.5f}, {g['median_ci'][1]:.5f}]"
        lines.append(f"| {g['regime']} | {g['h']} | {g['name']} | {median} | {ci} | {support} | {g['failed']}/{g['n']} |")
    lines+=['','Автоматические проверки подтверждают версии, полный бюджет обучения и парность входов. Интервалы рассчитаны bootstrap по записям (5000 повторов), интервалы долей — Wilson. Результаты ограничены одной системой и n=16. Не заявляются всеобщая иерархия методов, точная MATLAB/Octave идентичность или гарантированное принятие в Q1.','',
            'Обновлены раздел 7.6 статьи, две таблицы, Figure 7, supplement S6 и воспроизводимый пакет. Исходники авторов скачиваются по ссылкам и не распространяются в архиве.','',
            'Перед подачей остаются реальные авторы, аффилиации, декларации, человеческая проверка и открытый архив кода/данных.']
    ru=ROOT/'reports/full_adaptive_results_RU.md';ru.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    pypandoc.convert_file(str(ru),'docx',outputfile=str(ru.with_suffix('.docx')))
    path=ROOT/'reports/q1_revision_results_RU.md';text=path.read_text(encoding='utf-8')
    old='Опубликованный relative MSTLS регрессионный блок сверён с авторским кодом, но полный adaptive WSINDy не воспроизведён. IRLS-проверка относится к шагу Algorithm 2 WENDy и его OU-расширению; окна, остановка и регуляризация отличаются. Полный WmSINDy с совместным noise reconstruction не запущен. Поэтому статья не ранжирует эти pipelines и не утверждает превосходство над ними.'
    new='Полные adaptive WSINDy, WENDy и WmSINDy воспроизведены на закреплённых авторских версиях с описанными адаптерами Octave/TensorFlow. Завершены 448 парных fits, включая полный бюджет 40 000 шагов WmSINDy в каждом успешном запуске. Сравнение ограничено Ван-дер-Полем и 16 записями на условие; WENDy получает заданную структуру. Подробные результаты — reports/full_adaptive_results_RU.md и раздел 7.6 статьи. Компонентные проверки остаются отдельными экспериментами.'
    interim='Полные адаптивные авторские WSINDy, WENDy и WmSINDy успешно воспроизведены на контрольном примере Ван-дер-Поля с 25% независимого шума. WmSINDy прошёл все 40 000 шагов и обновления окон. Результаты и версии проверены; см. reports/full_adaptive_control_RU.md. Отдельная парная серия из 448 fits ещё выполняется на Kaggle. До её завершения сравнительные выводы статьи остаются основанными на ранее завершённых компонентных экспериментах; один контрольный пример не устанавливает статистического превосходства.'
    text=text.replace(old,new).replace(interim,new);path.write_text(text,encoding='utf-8')
    pypandoc.convert_file(str(path),'docx',outputfile=str(path.with_suffix('.docx')))
    path=ROOT/'reports/full_adaptive_protocol_RU.md';text=path.read_text(encoding='utf-8')
    text=text.replace('Проверка окружения пройдена. Полные вычислительные результаты добавляются только после успешного завершения и проверки выходных файлов Kaggle.',
                      'Полные запуски завершены; проверены SHA256 исходников и входов, 448 сравнительных fits и полный бюджет WmSINDy. Результаты находятся в reports/full_adaptive_results_RU.md; артефакты — в results/kaggle_adaptive_*.')
    path.write_text(text,encoding='utf-8')
    print('Full author-pipeline descriptions and Russian reports updated.')


if __name__=='__main__':main()
