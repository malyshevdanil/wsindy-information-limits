# Выполненные исследования на Kaggle

Пользователь подтвердил аккаунт malyshevdanil. Все notebooks приватные, CPU, без доступа к интернету во время исполнения. Credentials используются только локальным API-клиентом и не включены в исходники, документы или ZIP-пакеты. Перед экспериментами запускаются correctness tests; исходные серии содержат семь тестов, revision-fast — девять.

| Notebook | Статус | Выходы |
|---|---|---|
| [Smoke](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-smoke) | COMPLETE | 4 установочных fit; не используются как научные результаты |
| [Pilot](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-pilot) | COMPLETE | 4 320 test-fit, независимые train/validation/calibration |
| [Covariance](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-covariance) | COMPLETE | 640 fit на 80 seed; 1 000 scalar paths × 6 сеток |
| [Extended](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-extended) | COMPLETE | 24 000 residual vectors; 1 600 oracle-support fit; reference refinement |
| [Replication](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-replication) | COMPLETE | 7 200 oracle-support fit; 400 блоков по 32 независимые записи |
| [Revision-fast](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-revision-fast) | COMPLETE | 24 576 fit; 1 000 scalar paths; 2 000 coefficient paths на каждый из четырёх уровней; covariance-null и heterogeneous-IC исследования |
| [Reweighting](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-reweighting) | COMPLETE | 1 536 known-support fits OLS / WENDy iid IRLS core / OU IRLS core |
| [Geometry](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-geometry) | COMPLETE | Детерминированная проверка потери initial-state information; обе weak-influence identities проверены до 1e-9; новых случайных paths нет |
| [Revision, исходная медленная реализация](https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-revision) | RUNNING при последней проверке 11.10.2026 | Дублирует revision-fast с теми же seeds; не используется как независимое подтверждение |

Архивы скачаны в results/kaggle_STAGE_download/wsindy_STAGE_outputs.zip; CRC проверены, распакованные данные находятся в results/kaggle_STAGE. ZIP-хеши и manifests включены в reports/empirical_results.json. Smoke распакован CLI в results/kaggle_smoke_download/results_smoke. Логи сохраняют фактическое выполнение и версии; некоторые stdout строки дублируются notebook-логгером, число запусков считается по уникальным JSONL keys.

В Kaggle использовались Python 3.13.15, NumPy 2.1.3, SciPy 1.16.3 и Matplotlib 3.10.0. Время специализированного extended runner — 102.7 секунды, replication runner — 45.3 секунды, revision-fast — 256.8 секунды, reweighting — 57.8 секунды. Это время runners на Kaggle CPU, без очереди и упаковки; оно не является end-to-end сравнением алгоритмов. Notebook version устанавливается по metadata выполненного manifest, а не предполагается равным одному после повторных push.

## Воспроизведение

Сохранённые notebooks — неизменяемые execution snapshots. Их встроенные исходники могут отличаться от новых скриптов оформления и постобработки. scripts/validate_artifacts.py проверяет соответствие снимка source hash в скачанном manifest и отдельно сообщает изменения рабочего дерева. Для повторения точного эксперимента используют сохранённый notebook. Новый эксперимент следует выпускать отдельной версией с новой конфигурацией и выходной папкой.

API-клиент работает через python -m kaggle. При уже настроенном доступе:

```powershell
python -m kaggle kernels status malyshevdanil/wsindy-sampling-saturation-pilot
python scripts/download_kaggle_archive.py pilot
python scripts/analyze_completed_study.py
python scripts/analyze_revision.py --stage revision_fast
python scripts/build_publication.py
python scripts/validate_artifacts.py
```

Downloader повторяет временные сетевые ошибки по отдельным страницам, скачивает только исследовательский ZIP и проверяет CRC. Он не выводит credentials или подписанные URL. Все локальные команды выше выполняют чтение/анализ сохранённых выходов; Monte Carlo работает внутри Kaggle.

Каждый notebook содержит ZIP явно выбранных src/scripts/tests/configs/requirements в base64. Исходный Word и локальные results туда не включаются. Автоматический межсессионный resume не реализован: продолжение требует отдельно восстановить checkpoint, совпадающий по config/source hash.

Официальные сведения о CLI: [authentication](https://github.com/Kaggle/kaggle-cli/blob/main/docs/README.md#authentication), [kernels](https://github.com/Kaggle/kaggle-cli/blob/main/docs/kernels.md), [metadata](https://github.com/Kaggle/kaggle-cli/blob/main/docs/kernels_metadata.md).
