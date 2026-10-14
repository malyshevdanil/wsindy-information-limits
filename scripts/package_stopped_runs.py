"""Package recovered outputs and redacted logs for review."""
import hashlib,json,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    report=json.loads((ROOT/'reports/stopped_runs_collection.json').read_text());assert report['collection_complete']
    lines=['# Результаты остановленных запусков','',
           'Все пять отменённых запусков собраны. Отпечатки исходников и замороженных конфигураций совпадают; файлы снабжены SHA256. Автоматическая отправка заданий отключена.',
           '', '| Запуск | Сохранено |', '|---|---|']
    for run in report['runs']:
        if run.get('method_counts'):
            text=', '.join(f'{m}: {n}' for m,n in run['method_counts'].items())+f"; бутстреп-подгонок: {run['saved_bootstrap_refits']}"
        elif run.get('saved_cache_checkpoints'):text=f"{run['saved_cache_checkpoints']} блоков кеша"
        else:text=f"{len(run.get('checks',[]))} проверок геометрии; полного кеша нет"
        lines.append(f"| `{run['stage']}` | {text} |")
    lines+=['','Сохранено 23 профильные и 23 оракульные оценки, 21 завершённая поправка смещения и 1125 бутстреп-подгонок. В двух последних повторениях закончились только 35/50 и 40/50 бутстреп-подгонок. Выборка не образует завершённый Monte Carlo эксперимент; итоговые RMSE и накрытие по ней не публикуются.',
            '', 'Прежняя `crossfitA_geometry` сохранила 39 проверок при 5% шума; максимальная потеря дисперсии 8,8832% превышает порог 5%. Геометрия не допускается к основным расчётам.',
            '', 'Отдельные детерминированные расчёты успели завершиться. Для оракульных 480 тестов с fission прогноз n* ≈ 2668 / 2185 / 2029 при 5% / 10% / 20% шума. Это ориентир с истинными начальными состояниями и малошумовым приближением MLE, а не прогноз реализуемой cross-fitted оценки.',
            '', 'Поправка 1/n меняет прогноз нормы смещения при n=800 примерно на 0,34%; расхождение пилота она не объясняет. Проверка полного вектора остаётся предметом основного сравнения.',
            '', 'После остановки первоначальных пяти запусков три задания успели занять освободившиеся слоты до отключения локальных контроллеров: `repairA0_geometry`, `repairA1_geometry`, `cache_fix_v4`. На момент проверки они ещё работали; автоматическая очередь их больше не перезапускает.']
    summary=ROOT/'reports/stopped_runs_RU.md';summary.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    files=[summary,ROOT/'reports/stopped_runs_collection.json']
    for run in report['runs']:
        files.append(ROOT/run['log_file']);files.extend(ROOT/f['path'] for f in run['outputs'])
    for name in ['incidental_crossing_forecast.json','incidental_forecast_exact.json','incidental_finite_n_bias.json','pilot_curvature_bias_check.json']:
        files.append(ROOT/'reports'/name)
    archive=ROOT/'dist/stopped_kaggle_runs.zip';archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for file in sorted(set(files)):z.write(file,file.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
    (ROOT/'reports/stopped_runs_archive.json').write_text(json.dumps(dict(path=str(archive.relative_to(ROOT)),bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest()),indent=2))
    print(f'Archived {len(set(files))} files, {archive.stat().st_size} bytes')


if __name__=='__main__':main()
