"""Submit only explicitly named frozen endpoint follow-up notebooks."""
import argparse, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main(stage, api=None):
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.collect_stopped_runs import redact
    folder = ROOT/'kaggle'/stage
    meta = json.loads((folder/'kernel-metadata.json').read_text())
    if api is None:
        api = KaggleApi()
        try:api.authenticate()
        except Exception as error:
            print(stage+': authentication transport failed '+type(error).__name__,flush=True)
            return False
    try:
        status = api.kernels_status(meta['id'])
        state = getattr(status.status, 'name', str(status.status))
        print(f'{stage}: existing {state}; no duplicate submission', flush=True)
        return
    except Exception as error:
        response = getattr(error, 'response', None)
        absent_new_slug = isinstance(error, ValueError) and str(error).startswith('Cannot access kernel')
        if not absent_new_slug and (response is None or response.status_code != 404):
            print(f'{stage}: status unavailable ({type(error).__name__}); not submitting')
            print(redact(str(error))[:300])
            return
    try:
        response=api.kernels_push(str(folder))
        receipt={'stage':stage,'ref':response.ref,'version_number':response.version_number,'error':redact(response.error)}
        path=ROOT/'reports/endpoint_submission_receipts.json'
        history=json.loads(path.read_text()) if path.exists() else []
        history.append(receipt);path.write_text(json.dumps(history,indent=2),encoding='utf-8')
        if response.error:
            print(stage+': not started; '+redact(response.error),flush=True)
            return False
        print(stage+': submitted version '+str(response.version_number),flush=True)
        return True
    except Exception as error:
        response = getattr(error, 'response', None)
        print(stage+': submission failed '+type(error).__name__, flush=True)
        if response is not None: print(redact(response.text)[:1800], flush=True)
        elif hasattr(error,'doc'):print('Non-JSON server response: '+redact(error.doc)[:500],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage'); args = parser.parse_args()
    main(args.stage)
