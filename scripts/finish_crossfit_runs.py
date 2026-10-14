"""Run geometry, a technical smoke check, then four frozen crossfit shards."""
import argparse,hashlib,json,subprocess,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class StageFailure(RuntimeError):pass


def main():
    ledger=ROOT/'reports/crossfit_run_controller_v2.json'
    state=json.loads(ledger.read_text()) if ledger.exists() else dict(launched=[],downloaded=[])
    last={}
    def save():ledger.write_text(json.dumps(state,indent=2))
    def build(stage):
        if not (ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').exists():
            subprocess.run([sys.executable,'scripts/build_efficiency_kaggle.py',stage],cwd=ROOT,check=True)
    def push(api,stage):
        build(stage)
        # Recover a successful submission whose HTTP response was lost.
        try:
            existing=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
            value=existing.get('status') if isinstance(existing,dict) else existing.status
            status=getattr(value,'name',str(value)).split('.')[-1].upper()
            if status in ('RUNNING','QUEUED','COMPLETE','ERROR'):
                state['launched'].append(stage);save();print(json.dumps(dict(stage=stage,existing_execution_adopted=True)),flush=True);return
        except Exception as error:
            status=getattr(getattr(error,'response',None),'status_code',None)
            if isinstance(error,ValueError) and str(error).startswith('Cannot access kernel'):
                # Kaggle conceals absent private slugs with 403. Check access
                # to a known owner notebook before creating the new slug.
                api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match')
            elif status!=404 and not any(s in str(error).lower() for s in ('404','not found','notfound')):raise
        response=api.kernels_push(str(ROOT/'kaggle'/stage))
        error=response.get('error') if isinstance(response,dict) else getattr(response,'error',None)
        if error:raise RuntimeError('Kaggle push rejected: '+str(error))
        state['launched'].append(stage);save();print(json.dumps(dict(stage=stage,launched=True)),flush=True)
    def download(stage):
        from scripts.download_kaggle_archive import download
        download(stage)
        target=(ROOT/f'results/kaggle_{stage}').resolve();target.mkdir(exist_ok=True)
        with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
            assert z.testzip() is None
            for info in z.infolist():assert (target/info.filename).resolve().is_relative_to(target)
            z.extractall(target)
        assert json.loads((target/'completion.json').read_text())['status']=='complete'
        state['downloaded'].append(stage);save();print(stage+' downloaded and extracted',flush=True)
    while len(state['downloaded'])<7:
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            if 'crossfit_poly_probe_v2' not in state['launched']:push(api,'crossfit_poly_probe_v2')
            statuses={}
            for stage in [s for s in state['launched'] if s not in state['downloaded']]:
                response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                value=response.get('status') if isinstance(response,dict) else response.status
                statuses[stage]=getattr(value,'name',str(value)).split('.')[-1].upper()
            if statuses!=last:print(json.dumps(statuses),flush=True);last=statuses
            for stage,value in statuses.items():
                if value=='COMPLETE':download(stage)
                elif value in ('ERROR','FAILED','CANCELED','CANCELLED'):
                    raise StageFailure('Required crossfit stage failed: '+stage)
            if 'crossfit_poly_probe_v2' in state['downloaded'] and 'crossfit_geometry_v2' not in state['launched']:
                probe=json.loads((ROOT/'results/kaggle_crossfit_poly_probe_v2/probe.json').read_text())
                if probe['maximum_variance_loss']>json.loads((ROOT/'configs/crossfit_records.json').read_text())['maximum_cache_variance_loss']:
                    raise StageFailure('Polynomial probe failed prespecified accuracy threshold')
                push(api,'crossfit_geometry_v2')
            if 'crossfit_geometry_v2' in state['downloaded'] and 'crossfit_records_smoke' not in state['launched']:
                folder=ROOT/'results/kaggle_crossfit_geometry_v2'
                geometry=json.loads((folder/'geometry.json').read_text());assert geometry['cache_validated']
                (ROOT/'configs/crossfit_instruments.npz').write_bytes((folder/'crossfit_instruments.npz').read_bytes())
                geometry['cache_sha256']=hashlib.sha256((folder/'crossfit_instruments.npz').read_bytes()).hexdigest()
                (ROOT/'configs/crossfit_cache_geometry.json').write_text(json.dumps(geometry,indent=2))
                old=json.loads((ROOT/'reports/heterogeneous_records_results.json').read_text())['geometry']
                info={k:old[k] for k in ('information_grid','information_x','information_v','information_interpolation_relative_errors')}
                (ROOT/'configs/crossfit_information.json').write_text(json.dumps(info,indent=2))
                push(api,'crossfit_records_smoke')
            if 'crossfit_records_smoke' in state['downloaded']:
                if not state.get('smoke_validated'):
                    from scripts.validate_crossfit_runs import validate
                    try:validate(['crossfit_records_smoke'])
                    except Exception:raise StageFailure('Cross-fitting smoke validation failed; full inferential runs were not launched') from None
                    state['smoke_validated']=True;save()
                for i in range(4):
                    stage=f'crossfit_records_{i}'
                    if stage not in state['launched']:push(api,stage)
        except StageFailure:raise
        except Exception as error:
            print(json.dumps(dict(controller_retry=True,error_type=type(error).__name__)),flush=True)
        if len(state['downloaded'])<7:time.sleep(45)
    print('All cross-fitting runs downloaded',flush=True)
    if not state.get('analysis_and_manuscript_completed'):
        for script in ('analyze_crossfit_records.py','prepare_crossfit_theory.py','prepare_crossfit_results.py','clean_v13_manuscript.py','build_publication.py','build_ejs.py'):
            subprocess.run([sys.executable,'-m','scripts.'+Path(script).stem],cwd=ROOT,check=True)
        subprocess.run([sys.executable,str(ROOT/'scripts/validate_efficiency.py'),'--publication'],cwd=ROOT,check=True)
        state['analysis_and_manuscript_completed']=True;save()


if __name__=='__main__':main()
