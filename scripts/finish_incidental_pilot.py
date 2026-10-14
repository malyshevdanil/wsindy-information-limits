"""Recoverable Kaggle smoke-to-pilot queue; freeze escalation from pilot only."""
import argparse,hashlib,json,subprocess,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main(refined=False):
    config='incidental_pilot_v2.json' if refined else 'incidental_pilot.json'
    ledger=ROOT/('reports/incidental_pilot_v2_controller.json' if refined else 'reports/incidental_pilot_controller.json')
    state=json.loads(ledger.read_text()) if ledger.exists() else dict(launched=[],downloaded=(['incidental_smoke'] if refined else []));last={}
    def save():ledger.write_text(json.dumps(state,indent=2))
    while len(state['downloaded'])<3:
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            planned=['incidental_smoke'] if 'incidental_smoke' not in state['downloaded'] else ([f'incidental_pilot_v2_{i}' for i in range(2)] if refined else ['incidental_pilot_0','incidental_pilot_1'])
            for stage in planned:
                if stage in state['launched']:continue
                if not (ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').exists():subprocess.run([sys.executable,'scripts/build_efficiency_kaggle.py',stage],cwd=ROOT,check=True)
                try:
                    api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'));exists=True
                except ValueError:
                    api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match');exists=False
                if not exists:
                    response=api.kernels_push(str(ROOT/'kaggle'/stage));error=getattr(response,'error',None)
                    if error:raise RuntimeError('Kaggle rejected pilot submission')
                state['launched'].append(stage);save();print(json.dumps(dict(stage=stage,execution_adopted=exists)),flush=True)
            statuses={}
            for stage in state['launched']:
                if stage in state['downloaded']:continue
                response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                value=response.get('status') if isinstance(response,dict) else response.status;statuses[stage]=getattr(value,'name',str(value)).split('.')[-1].upper()
            if statuses!=last:print(json.dumps(statuses),flush=True);last=statuses
            for stage,status in statuses.items():
                if status in ('ERROR','FAILED','CANCELED','CANCELLED'):
                    state['failure']=dict(stage=stage,status=status);save();print('Incidental pilot failed; promotion stopped',flush=True);return
                if status!='COMPLETE':continue
                from scripts.download_kaggle_archive import download
                download(stage);folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as archive:
                    assert archive.testzip() is None
                    for name in archive.namelist():assert (folder/name).resolve().is_relative_to(folder.resolve())
                    archive.extractall(folder)
                from scripts.validate_nuisance_runs import snapshot
                snapshot(stage,config);rows=json.loads((folder/'rows.json').read_text())
                assert all(r['diagnostic']['success'] and r['diagnostic']['remaining_newton_step']<1e-5 for r in rows)
                if stage=='incidental_smoke':assert all(c['success'] and c['maximum_coefficient_difference']<=1e-5 for c in json.loads((folder/'technical_checks.json').read_text()))
                state['downloaded'].append(stage);save()
        except AssertionError:
            state['validation_failure']='frozen optimizer diagnostic failed; refinement required';save();return
        except Exception as error:
            print(json.dumps(dict(pilot_controller_retry=True,error_type=type(error).__name__)),flush=True)
        if len(state['downloaded'])<3:time.sleep(45)
    subprocess.run([sys.executable,'-m','scripts.analyze_incidental_pilot']+(['--refined'] if refined else []),cwd=ROOT,check=True)
    state['pilot_analysis_complete']=True;save()


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--refined',action='store_true');main(a.parse_args().refined)
