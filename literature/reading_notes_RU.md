# Обзор литературы и карта вклада

Даты проверки: 10–11 октября 2026 года. Это целевой обзор первичных источников, а не завершённый систематический обзор. Исходное задание прочитано в `Sampling_Saturation_WSINDy_Research_Plan_RU.docx` и его Markdown-источнике. Имена методов ниже не означают, что их авторские реализации уже воспроизведены.

## Проверенные источники

| ID | Работа и версия | Прочитано / проверено | Положение, которое используется в тексте |
|---|---|---|---|
| R1 | Brunton, Proctor, Kutz. Discovering governing equations from data by sparse identification of nonlinear dynamics. PNAS 113, 3932–3937 (2016). [DOI](https://doi.org/10.1073/pnas.1517384113), [открытый текст](https://pmc.ncbi.nlm.nih.gov/articles/PMC4839439/) | Постановка и принцип sparse library | Разреженная библиотека задаёт одновременно структуру и коэффициенты модели. |
| R2 | Messenger, Bortz. Weak SINDy: Galerkin-Based Data-Driven Model Selection. MMS 19(3), 1474–1497 (2021). [DOI](https://doi.org/10.1137/20M1343166), [arXiv v3](https://arxiv.org/abs/2005.04339v3), [текст](https://pmc.ncbi.nlm.nih.gov/articles/PMC10795802/) | Слабая постановка, действие окон на шум; авторский код найден | Перенос производной на тестовую функцию и учёт ковариации уже присутствуют в слабой идентификации. |
| R3 | Messenger, Bortz. Asymptotic consistency of the WSINDy algorithm in the limit of continuum data. IMA JNA 45(6), 3264–3312 (2025), online 12.12.2024. [DOI](https://doi.org/10.1093/imanum/drae086), [arXiv v1](https://arxiv.org/abs/2211.16000v1), [полный препринт](https://arxiv.org/html/2211.16000) | Assumptions, библиотечное смещение, consistency, Discussion III; журнал и препринт проверены отдельно | Результат зависит от модели шума; Discussion III оставляет коррелированный измерительный шум будущей работе. Корреляция после фильтрации iid-данных уже анализируется и отличается от исходного OU-шума с фиксированным τ. |
| R4 | Bortz, Messenger, Dukic. Direct Estimation of Parameters in ODE Models Using WENDy: Weak-form Estimation of Nonlinear Dynamics. BMB 85, 110 (2023). [DOI](https://doi.org/10.1007/s11538-023-01208-6), [текст](https://pmc.ncbi.nlm.nih.gov/articles/PMC10556152/) | §2.1–2.3, EiV, Taylor residual, IRLS | Ошибки в G и b связаны. Оценка параметров заданной структуры отличается от выбора sparse support. |
| R5 | López, Naranjo, Salazar, Moore. Weak-form modified sparse identification of nonlinear dynamics. [arXiv v1](https://arxiv.org/abs/2410.17838v1), [PDF](https://arxiv.org/pdf/2410.17838v1), [журнальный DOI](https://doi.org/10.1016/j.jcp.2025.114410) | Методы и эксперименты препринта; DOI подтверждён README авторского кода, страница издателя не прочитана | Совместное восстановление динамики и шума — существующий подход. В доступном препринте нельзя автоматически искать все эксперименты журнальной версии. |
| R5-code | [Авторский WmSINDy](https://github.com/CristianLopez85/WmSINDy) | Список файлов, README | Есть `Fig31_Colored_noise` и `Fig32_Correlated Noise`. Это подтверждает близость baseline, но не эквивалентность его шума нашей модели фиксированного τ. |
| R6 | Chen, Fasel, Bizyaeva. Fourier Weak SINDy: Spectral Test Function Selection for Robust Model Identification. [arXiv v2, 27.07.2026](https://arxiv.org/abs/2604.20141v2), [полный текст](https://arxiv.org/html/2604.20141v2) | Метод спектральных окон, обсуждение baselines, metadata | Ортогональные синусоидальные окна выбираются через multitaper spectrum. Карточка v2 сообщает публикацию в proceedings L4DC 2026; отдельная страница proceedings в этом проходе не найдена. В исходном плане был указан только v1. |
| R7 | Isserlis. On a formula for the product-moment coefficient of any order of a normal frequency distribution in any number of variables. Biometrika 12(1–2), 134–139 (1918). [DOI](https://doi.org/10.1093/biomet/12.1-2.134) | Библиографическая карточка; наш вывод через generating function дан полностью | Гауссовские моменты и Hermite orthogonality — известный инструмент, а не заявка на новизну. |
| R8 | Uhlenbeck, Ornstein. On the Theory of the Brownian Motion. Physical Review 36, 823–841 (1930). [DOI](https://doi.org/10.1103/PhysRev.36.823) | Карточка первичного издателя | Исторический источник OU-модели; её ковариация и используемая рекурсия выведены/определены в нашей постановке. |

## Проверка новизны

Поисковые запросы выполнены через web search 10.10.2026:

- `WSINDy correlated measurement noise infill`;
- `WSINDy sampling correlated noise saturation`;
- `weak form equation discovery fixed horizon`;
- `WSINDy sampling saturation`;
- `weak-form Hermite correlated identification`;
- точные названия R2–R6 и просмотр авторских репозиториев.

Готового решения всей постановки в просмотренном наборе не обнаружено. Это ограниченный результат поиска. Citation chaining, проверка журнальной версии R5 и поиск общего polynomial Gaussian errors-in-variables literature ещё нужны перед заявлением новизны и выбором журнала.

| Возможное заявление | Решение |
|---|---|
| Впервые применили слабую форму к коррелированному шуму | Исключено: близкая литература и авторский код уже существуют. |
| Впервые использовали GLS или подобрали окна | Исключено: R2, R4, R6. |
| Сгущение сетки не снижает дисперсию OU-среднего до нуля | Стандартный подготовительный расчёт. |
| Дана несостоятельность конкретной слабой ratio-оценки на фиксированном горизонте | Наше доказательство оформлено; оригинальность общего статистического механизма не заявляется. |
| Точная covariance слабого polynomial residual как конечная сумма OU-ядер с временами τ/q | Получен самостоятельный вывод в `theory.md`; перспектива вклада требует более широкого поиска аналогов. |
| Критерий grid coarsening сохраняет качество sparse discovery | Гипотеза и прототип. Выполненный Kaggle-пилот не подтверждает критерий во всех случаях; нужны доработка и сильные baselines. |

## Как положения литературы встроены в статью

В `manuscript/article_en.md` ссылки стоят непосредственно в Introduction, Related work и обсуждении. R3 обосновывает различение sampling models; R4 — анализ совместных ошибок; R5 — обязательность сильного baseline; R6 — фиксацию окон в физических единицах; R7 — происхождение гауссовской алгебры. Собственные теоремы сопровождаются доказательствами и явно отделены от эмпирических гипотез.

## Следующий проход чтения

1. Сверить актуальную журнальную R5 с v1 и расширить прочитанный ниже набор noise-generation scripts.
2. Пройти цитирующие R3/R4 работы по fixed-domain covariance, estimating equations и polynomial bias correction.
3. Проверить оригинальность точной Hermite covariance до использования слова novel.
4. Воспроизвести авторский WSINDy/MSTLS и выбранные WmSINDy-ячейки; наш STLS-прототип пока не является их репликацией.

## Уточнение по авторскому коду, 11.10.2026

Проверен commit `f1a781009a5c6adfd5650bec63ef30cc0673ae97`; пути и SHA-256 просмотренных файлов сохранены в `author_code_audit.json`.

В [Fig32 / Fig16_Wm_corr.py](https://github.com/CristianLopez85/WmSINDy/blob/f1a781009a5c6adfd5650bec63ef30cc0673ae97/Fig32_Correlated%20Noise/Fig16_Wm_corr.py) каждый момент получает новый многомерный Gaussian vector. Матрица Q создаёт корреляцию между координатами состояния; независимые вызовы генератора не задают корреляцию во времени. Поэтому имя папки само по себе было недостаточно для классификации шума.

В [Fig31 / Fig15_Wm_pink.py](https://github.com/CristianLopez85/WmSINDy/blob/f1a781009a5c6adfd5650bec63ef30cc0673ae97/Fig31_Colored_noise/Fig15_Wm_pink.py) временной цветной шум получают спектральным фильтром с амплитудой ∝1/√f. Частоты в генераторе заданы в отсчётах; скрипт использует фиксированный dt=0.01. По этому файлу не установлена серия OU-инфилл с фиксированным физическим τ. Он остаётся релевантным цветным baseline. Данный вывод относится к двум прочитанным скриптам, а не ко всему репозиторию или журнальной статье.
