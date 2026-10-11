"""Durable local orchestration of Kaggle runs and verified publication updates.

No simulations execute locally. Credentials are read by KaggleApi only. This
worker can survive the editor/chat turn and resumes from completed archives.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess
import sys
import time
from kaggle.api.kaggle_api_extended import KaggleApi
ROOT=Path(__file__).resolve().parents[1]
STAGES=['adaptive_smoke','adaptive_classical',*[f'adaptive_wm_{i}' for i in range(4)]]


def command(*args):
    r=subprocess.run([sys.executable,*args],cwd=ROOT,timeout=600)
    if r.returncode:raise RuntimeError('Local preparation/postprocessing command failed')


def record(statuses,phase,**extra):
    obj={'updated_utc':datetime.now(timezone.utc).isoformat(),'phase':phase,'stages':statuses,**extra}
    p=ROOT/'reports/full_adaptive_job_status.json';temp=p.with_suffix('.tmp')
    temp.write_text(json.dumps(obj,indent=2),encoding='utf-8');temp.replace(p)


def main():
    statuses={s:'UNKNOWN' for s in STAGES};api=KaggleApi();api.authenticate()
    record(statuses,'monitoring')
    for stage in STAGES:
        if not (ROOT/'kaggle'/stage/'kernel-metadata.json').exists():command('scripts/build_adaptive_kaggle.py',stage)
    while True:
        for stage in STAGES:
            folder=ROOT/f'results/kaggle_{stage}'
            if (folder/'manifest.json').exists():
                manifest=json.loads((folder/'manifest.json').read_text())
                if not manifest.get('completed'):raise ValueError('Invalid completed output')
                statuses[stage]='DOWNLOADED';continue
            slug='malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-')
            try:
                response=api.kernels_status(slug)
                status=response.get('status') if isinstance(response,dict) else response.status
                status=str(getattr(status,'name',status)).upper().rsplit('.',1)[-1];statuses[stage]=status
                if status=='COMPLETE':
                    try:
                        command('scripts/download_kaggle_archive.py',stage)
                        command('scripts/unpack_kaggle_outputs.py',stage)
                        statuses[stage]='DOWNLOADED'
                    except (RuntimeError,subprocess.TimeoutExpired) as download_error:
                        statuses[stage]='DOWNLOAD_RETRY'
                        print(json.dumps({'stage':stage,'download_retry':type(download_error).__name__}),flush=True)
                elif status in ('ERROR','CANCELLED'):
                    record(statuses,'execution_failed',failed_stage=stage)
                    raise RuntimeError('Kaggle execution failed; inspect preserved logs')
            except RuntimeError:raise
            except Exception as error:
                code=getattr(getattr(error,'response',None),'status_code',None)
                missing=(code==404 or isinstance(error,ValueError) and 'kernels.get' in str(error))
                if missing and stage.startswith('adaptive_wm_') and statuses[stage] in ('UNKNOWN','WAITING_FOR_KAGGLE_SLOT'):
                    try:
                        response=api.kernels_push(str(ROOT/'kaggle'/stage))
                        rejection=response.get('error') if isinstance(response,dict) else getattr(response,'error',None)
                        statuses[stage]='WAITING_FOR_KAGGLE_SLOT' if rejection else 'QUEUED'
                        print(json.dumps({'stage':stage,'push_rejection':rejection,'state':statuses[stage]}),flush=True)
                    except Exception as push_error:
                        print(json.dumps({'stage':stage,'transient_error':type(push_error).__name__}),flush=True)
                else:
                    print(json.dumps({'stage':stage,'transient_error':type(error).__name__}),flush=True)
            record(statuses,'monitoring')
        print(json.dumps(statuses),flush=True)
        if all(v=='DOWNLOADED' for v in statuses.values()):break
        time.sleep(45)
    record(statuses,'validating')
    command('scripts/validate_full_adaptive.py')
    command('scripts/analyze_full_adaptive.py')
    command('scripts/finalize_full_adaptive.py')
    command('scripts/build_publication.py')
    command('scripts/validate_artifacts.py')
    command('scripts/package_publication.py')
    command('scripts/validate_full_adaptive.py','--publication')
    record(statuses,'complete',results='reports/full_adaptive_results.json',article='manuscript/article_revised.docx')
    print('Complete adaptive reproduction and publication update finished.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.parse_args();main()
