"""Queue bounded-duration likelihood/bootstrap jobs after validated amplitude pilot.

This ledger tracks likelihood components, not completion of the matched weak
comparison. It preserves each execution snapshot and stops on failed refits.
"""
import json,subprocess,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    hold=ROOT/'reports/kaggle_submission_hold.json'
    if hold.exists() and json.loads(hold.read_text()).get('enabled'):
        print('Automatic Kaggle submissions held; read-only collection remains available',flush=True);return
    compact=(ROOT/'configs/incidental_compact.json').exists();trimmed=(ROOT/'configs/incidental_trimmed.json').exists()
    config=ROOT/('configs/incidental_compact.json' if compact else 'configs/incidental_trimmed.json' if trimmed else 'configs/incidental_main.json');ledger=ROOT/'reports/incidental_main_controller.json'
    state=json.loads(ledger.read_text()) if ledger.exists() else dict(launched=[],downloaded=[])
    def save():ledger.write_text(json.dumps(state,indent=2))
    while not config.exists():time.sleep(45)
    p=json.loads(config.read_text());blocks=(p['monte_carlo_repetitions']+p['repetitions_per_job']-1)//p['repetitions_per_job']
    stem='inc_short' if compact else 'inc_trim' if trimmed else 'inc_main'
    jobs=[f'{stem}_l{li}_n{n}_b{b:02d}' for li in range(len(p['levels'])) for n in reversed(p['record_counts']) for b in range(blocks)]
    if trimmed:
        old=json.loads((ROOT/'reports/incidental_main_controller_pre_budget.json').read_text())
        state['legacy_jobs_outside_amended_queue']=[s for s in old['launched'] if s!='inc_main_smoke']
        # Recompute all amended cells on unchanged seeds. Legacy expensive jobs
        # are never adopted, relaunched, or mistaken for amended completion.
        state['launched']=[s for s in state['launched'] if s=='inc_main_smoke' or s.startswith(stem+'_')]
        state['downloaded']=[s for s in state['downloaded'] if s=='inc_main_smoke' or s.startswith(stem+'_')]
        state['protocol_file']=config.name;save()
    while len(state['downloaded'])<len(jobs)+1:
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            active=[s for s in state['launched'] if s not in state['downloaded']]
            for stage in active:
                response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                value=response.get('status') if isinstance(response,dict) else response.status
                status=getattr(value,'name',str(value)).split('.')[-1].upper()
                if status=='CANCEL_ACKNOWLEDGED':
                    submission=api.kernels_push(str(ROOT/'kaggle'/stage))
                    print(stage+' waiting for CPU slot' if getattr(submission,'error',None) else stage+' resubmitted frozen snapshot',flush=True)
                    continue
                if status in ('ERROR','FAILED','CANCELED','CANCELLED'):
                    state['failure']=dict(stage=stage,status=status);save();return
                if status!='COMPLETE':continue
                from scripts.download_kaggle_archive import download
                download(stage);folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None
                    assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                    z.extractall(folder)
                from scripts.validate_nuisance_runs import snapshot
                snapshot(stage,'incidental_compact.json' if stage.startswith('inc_short_') else 'incidental_trimmed.json' if stage.startswith('inc_trim_') else 'incidental_main.json')
                c=json.loads((folder/'completion.json').read_text());assert c['status']=='complete' and c['discarded_records_or_fits']==0
                rows=json.loads((folder/'rows.json').read_text());boot=json.loads((folder/'bootstrap_fits.json').read_text())
                assert len(boot)==c['bootstrap_refits']
                assert len(rows)==(3 if boot else 2)*c['repetitions']
                assert all(r['diagnostic']['success'] for r in rows+boot)
                state['downloaded'].append(stage);save();print(stage+' validated',flush=True)
            candidates=['inc_main_smoke'] if 'inc_main_smoke' not in state['downloaded'] else jobs
            if trimmed and not (ROOT/'reports/incidental_forecast_exact.json').exists():
                candidates=[]
                state['waiting_for_analytic_forecast']=True;save()
            else:state.pop('waiting_for_analytic_forecast',None)
            for stage in candidates:
                if stage in state['launched']:continue
                if len(set(state['launched'])-set(state['downloaded']))>=2:break
                if not (ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').exists():
                    subprocess.run([sys.executable,'scripts/build_efficiency_kaggle.py',stage],cwd=ROOT,check=True)
                try:api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                except ValueError:
                    api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match')
                    response=api.kernels_push(str(ROOT/'kaggle'/stage))
                    if getattr(response,'error',None):raise RuntimeError('Submission rejected')
                state['launched'].append(stage);save();print(stage+' submitted or adopted',flush=True)
        except (AssertionError,FileNotFoundError) as error:
            state['validation_failure']=type(error).__name__;save();return
        except Exception as error:print(json.dumps(dict(retry=True,error_type=type(error).__name__)),flush=True)
        if len(state['downloaded'])<len(jobs)+1:time.sleep(45)
    state['likelihood_components_complete']=True
    state['matched_weak_component_pending']=True;save()


if __name__=='__main__':main()
