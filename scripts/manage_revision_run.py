"""Small retrying Kaggle controller; never print credentials or signed URLs."""
import argparse
import json
from pathlib import Path
import time
from kaggle.api.kaggle_api_extended import KaggleApi
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('action',choices=['push','status'])
parser.add_argument('--stage',default='revision',choices=['revision','revision_fast','reweighting','geometry','reviewer','adaptive_probe','adaptive_smoke','adaptive_octave_smoke','adaptive_classical','adaptive_wm',*[f'adaptive_wm_{i}' for i in range(4)]])
args=parser.parse_args();action=args.action;stage=args.stage
for attempt in range(6):
    try:
        api=KaggleApi();api.authenticate()
        if action=='push':
            response=api.kernels_push(str(ROOT/'kaggle'/stage))
            error = response.get('error') if isinstance(response,dict) else getattr(response,'error',None)
            if error:
                print(json.dumps({'action':'push','rejected':True,'reason':error}),flush=True)
                raise SystemExit(1)
            print(json.dumps({'action':'push','status':'response_received','url':'https://www.kaggle.com/code/malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-')}),flush=True)
        else:
            print(api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-')),flush=True)
        break
    except Exception as error:
        print(json.dumps({'attempt':attempt+1,'error_type':type(error).__name__}),flush=True)
        if attempt==5:raise SystemExit(1)
        time.sleep(2)
