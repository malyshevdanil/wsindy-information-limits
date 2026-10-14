"""Read-only collector; never submits, cancels or restarts Kaggle notebooks."""
import argparse, json, sys, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main(status_only=False):
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.download_kaggle_archive import download
    from scripts.collect_stopped_runs import snapshot
    api = KaggleApi(); api.authenticate(); statuses = []
    stages = [f'boundary_confirm_{i}' for i in range(4)]+['support_extra']
    for stage in stages:
        try:
            response = api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_', '-'))
            status = getattr(response.status, 'name', str(response.status)); record = {'stage': stage, 'status': status}
            print(stage+': '+status, flush=True)
            if status == 'COMPLETE' and not status_only:
                download(stage); folder = ROOT/f'results/kaggle_{stage}'; folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None
                    assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                    z.extractall(folder)
                manifest = json.loads((folder/'manifest.json').read_text())
                record['snapshot'] = snapshot(stage, manifest)
                record['completion'] = json.loads((folder/'completion.json').read_text())
                assert record['completion']['status'] == 'complete'
            elif status == 'ERROR' and not status_only:
                dest = ROOT/f'results/failed_{stage}'; dest.mkdir(exist_ok=True)
                api.kernels_output('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_', '-'),
                                   str(dest), file_pattern='(^|/)results_'+stage+'/', quiet=True)
                for path in dest.rglob('completion.json'): record['failure'] = json.loads(path.read_text())
        except Exception as exc:
            record = {'stage': stage, 'collection_error_type': type(exc).__name__}
            print(stage+': read-only audit error '+type(exc).__name__, flush=True)
        statuses.append(record)
        (ROOT/'reports/endpoint_confirmation_status.json').write_text(json.dumps(statuses, indent=2), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--status-only', action='store_true')
    main(parser.parse_args().status_only)
