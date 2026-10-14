"""Download and validate the repair probe before any full cache promotion."""
import argparse,json,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main(large_noise=False,v4=False,compact=False):
    hold=ROOT/'reports/kaggle_submission_hold.json'
    if hold.exists() and json.loads(hold.read_text()).get('enabled'):
        print('Automatic Kaggle submissions held',flush=True);return
    stage='cache120a' if compact else 'cache_fix_v4' if v4 else 'cache_fix_a' if large_noise else 'cache_fix_probe'
    config='cache120_main_probe.json' if compact else 'cache_repair_probe_v4.json' if v4 else 'cache_repair_probe_A.json' if large_noise else 'cache_repair_probe.json'
    report=ROOT/('reports/cache120_main_probe_validation.json' if compact else 'reports/cache_repair_probe_v4_validation.json' if v4 else 'reports/cache_repair_probe_A_validation.json' if large_noise else 'reports/cache_repair_probe_validation.json')
    while not report.exists():
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            try:
                response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
            except Exception as error:
                # Verify the account against a known kernel before treating a
                # missing/rejected new entry as an unsubmitted queued probe.
                if not (v4 or compact):raise
                api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match')
                submitted=api.kernels_push(str(ROOT/'kaggle'/stage))
                print('Probe waiting for CPU slot' if getattr(submitted,'error',None) else 'Probe submitted',flush=True)
                time.sleep(45);continue
            value=response.get('status') if isinstance(response,dict) else response.status
            status=getattr(value,'name',str(value)).split('.')[-1].upper()
            if (v4 or compact) and status=='CANCEL_ACKNOWLEDGED':
                submitted=api.kernels_push(str(ROOT/'kaggle'/stage))
                print('Probe waiting for CPU slot' if getattr(submitted,'error',None) else 'Probe resubmitted unchanged snapshot',flush=True)
                time.sleep(45);continue
            if status in ('ERROR','FAILED','CANCELED','CANCELLED'):
                report.write_text(json.dumps(dict(passed=False,stage=stage,status=status)))
                print('Repair probe failed; full cache rebuild blocked',flush=True);return
            if status=='COMPLETE':
                from scripts.download_kaggle_archive import download
                download(stage);folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None
                    assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                    z.extractall(folder)
                from scripts.validate_nuisance_runs import snapshot
                verified=snapshot(stage,config);data=json.loads((folder/'probe.json').read_text())
                assert data['maximum_variance_loss']<=data['threshold']==.05
                report.write_text(json.dumps(dict(passed=True,stage=stage,execution_snapshot=verified,probe=data),indent=2))
                print('Repair probe validated; full geometry promotion permitted',flush=True);return
        except AssertionError:
            report.write_text(json.dumps(dict(passed=False,stage=stage,validation_failure=True)));return
        except Exception as error:print(json.dumps(dict(retry=True,error_type=type(error).__name__)),flush=True)
        time.sleep(45)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--large-noise',action='store_true');a.add_argument('--v4',action='store_true');a.add_argument('--compact',action='store_true')
    args=a.parse_args();main(args.large_noise,args.v4,args.compact)
