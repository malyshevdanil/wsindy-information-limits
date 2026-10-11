# Efficient weak-form estimation under correlated measurement error

Исследование по Word-плану; журнал из §21 плана — **Applied Mathematics and Computation**. Новая редакция описывает полный boundary-free adjoint-класс, доказывает точность границы Шура и условную эффективность GLS в последовательном малошумовом пределе. Прежние вычисления сохранены в supplement и неизменных архивах исполнения.

- [Статья: Word](manuscript/article_revised.docx), [PDF](manuscript/article_revised.pdf), [исходный текст](manuscript/article_revised.md), [LaTeX](manuscript/article_revised.tex).
- [Supplement: PDF](manuscript/supplement.pdf), [Word](manuscript/supplement.docx), [исходный текст](manuscript/supplement.md).
- [Строгие доказательства weak-эффективности](manuscript/weak_efficiency_theory.md), [PDF](manuscript/weak_efficiency_theory.pdf).
- [Результаты двух экспериментов](reports/efficiency_results_RU.md), [Word-отчёт](reports/efficiency_results_RU.docx), [численные результаты и интервалы](reports/efficiency_results.json), [проверка](reports/efficiency_validation.json).
- [Литературные основания и границы новизны](literature/efficiency_notes_RU.md), [проверка 19 DOI](literature/publication_metadata_audit.json).
- [Полные авторские методы: 448 оценок](reports/full_adaptive_results_RU.md), [проверка](reports/full_adaptive_validation.json).
- [Информация Van der Pol, оптимальные семейства и bias/dispersion](reports/reviewer_results_RU.md).
- [Пакет журнала и авторские поля](manuscript/submission/README_RU.md).

Скалярный тест с fitted pilot при SD=0.1 снизил медиану ошибки с 0.2244 до 0.0598 на 1000 путях. Для Van der Pol на 128 путях OU-RMSE с 480 окнами при шуме 0.5% равна 0.001737 при weak-границе 0.001690; при 5% — 0.025583 при границе 0.016904. В локальном Gaussian расчёте OU-RMSE снижается 0.019605 → 0.017043 → 0.016938. Парные интервалы Монте-Карло для преимущества OU над iid включают ноль. Финальная серия состоит из 1536 оценок; все достигли критерия остановки. Пять-итерационный контроль использует те же пути и не добавляет независимых данных.

Эксперименты выполняются на Kaggle. Локально проводится постобработка записанных результатов и сборка документов. Снимки notebook сохраняют исполнявшийся код и конфигурацию; изменения текущих исходников не переписывают их. Все методы и амплитуды сравниваются с парными seed, которые нельзя считать независимыми репликациями. Протокол `configs/efficiency.json` и каждый manifest явно задают параметры своего исполнения. Новые IRLS-проверки относятся к WENDy-подобной компоненте с заданной структурой; полный авторский WENDy остаётся отдельным экспериментом.

Публичный архив: [GitHub v1.1.0](https://github.com/malyshevdanil/wsindy-information-limits/tree/v1.1.0). Точный commit и проверка всех публичных blob записаны в `reports/public_archive_validation.json`. `release_manifest.json` в публичном снимке фиксирует SHA256 каждого файла; DOI не присвоен. Исходные авторские библиотеки скачиваются по закреплённым версиям и не перераспределяются; наши адаптеры включены.

После распаковки архивов исходников и данных в одну папку:

```powershell
python -m pytest -q tests
python scripts/analyze_revision.py --stage revision_fast
python scripts/analyze_full_adaptive.py
python scripts/analyze_reviewer_validation.py
python scripts/prepare_efficiency_manuscript.py
python scripts/analyze_efficiency.py
python scripts/build_publication.py
python scripts/validate_artifacts.py
python scripts/validate_efficiency.py --publication
python scripts/validate_reviewer.py --publication
python scripts/package_publication.py
python scripts/validate_full_adaptive.py --publication
```

Постобработка прежних исследований создаёт их численные отчёты и фигуры, не перезаписывая новую статью. Стабильные шаблоны и supporting validation включены в архив исходников. Численные исследования повторно запускаются только на Kaggle; архивы outputs проверяются по CRC, SHA256 и fingerprint исполненного notebook. После установки `requirements-publication.txt` нужен portable Tectonic (`scripts/install_publication_tools.py` на Windows). На других ОС следует указать соответствующий Tectonic binary в сборщике. Формулы Word — редактируемые OMML. `article_en.*` совпадают с новой статьёй.

В `dist/` находятся отдельные ZIP рукописей, исходников и данных с `manifest.json`. Обе версии статьи до переработки сохранены в `manuscript/archive/`. Авторские сведения остаются для заполнения; человеческая проверка содержания и деклараций перед подачей ещё требуется. Статья в журнал не отправлена.
