"""Read-only collection of the four authorized ablation notebooks."""
import json, sys, time, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.download_kaggle_archive import download
    from scripts.collect_stopped_runs import snapshot
    from scripts.analyze_boundary_ablation import main as analyze
    api = KaggleApi(); api.authenticate()
    while True:
        records = []
        for i in range(4):
            stage = f'boundary_ablate_v2_{i}'
            try:
                response = api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                state = getattr(response.status, 'name', str(response.status)); r = {'stage': stage, 'status': state}
                print(stage+': '+state, flush=True)
                if state == 'COMPLETE':
                    download(stage); folder = ROOT/f'results/kaggle_{stage}'; folder.mkdir(exist_ok=True)
                    with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                        assert z.testzip() is None and all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                        z.extractall(folder)
                    r['audit'] = snapshot(stage, json.loads((folder/'manifest.json').read_text()))
                    r['completion'] = json.loads((folder/'completion.json').read_text())
                elif state == 'ERROR':
                    target = ROOT/f'results/failed_{stage}'; target.mkdir(exist_ok=True)
                    api.kernels_output('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'),str(target),quiet=True)
            except Exception as error:
                r = {'stage': stage, 'read_error': type(error).__name__}
            records.append(r)
        (ROOT/'reports/boundary_ablation_status.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
        if any(r.get('status') in ('ERROR','CANCELED','CANCEL_ACKNOWLEDGED') for r in records):
            print('Failed/cancelled run preserved; no automatic resubmission'); return
        if all(r.get('completion',{}).get('status') == 'complete' for r in records):
            analyze(); return
        time.sleep(45)


if __name__ == '__main__': main()
