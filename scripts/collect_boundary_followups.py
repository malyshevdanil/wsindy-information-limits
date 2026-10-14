"""Read-only collection for the six declared generalization/selection studies."""
import argparse,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
STAGES=['endpoint_vdp5','endpoint_lorenz']+[f'support_wald_v2_{i}' for i in range(4)]


def main(status_only=False):
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.download_kaggle_archive import download
    from scripts.collect_stopped_runs import snapshot
    api=KaggleApi();api.authenticate(); records=[]
    for stage in STAGES:
        try:
            response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
            state=getattr(response.status,'name',str(response.status));r={'stage':stage,'status':state}
            print(stage+': '+state,flush=True)
            if state=='COMPLETE' and not status_only:
                download(stage,api=api);folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None and all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist());z.extractall(folder)
                r['audit']=snapshot(stage,json.loads((folder/'manifest.json').read_text()))
                r['completion']=json.loads((folder/'completion.json').read_text())
            elif state=='ERROR' and not status_only:
                target=ROOT/f'results/failed_{stage}';target.mkdir(exist_ok=True)
                api.kernels_output('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'),str(target),quiet=True)
        except Exception as error:r={'stage':stage,'read_error':type(error).__name__};print(stage+': read error '+type(error).__name__,flush=True)
        records.append(r)
    target=ROOT/'reports/boundary_followup_status.json'
    temporary=target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(records,indent=2),encoding='utf-8')
    temporary.replace(target)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--status-only',action='store_true');args=parser.parse_args();main(args.status_only)
