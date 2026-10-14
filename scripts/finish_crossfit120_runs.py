"""Resume the matched 120-test cross-fitting and oracle fission-price protocol."""
import argparse,hashlib,json,subprocess,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main(main_a=False,repair_family=None):
    hold=ROOT/'reports/kaggle_submission_hold.json'
    if hold.exists() and json.loads(hold.read_text()).get('enabled'):
        print('Automatic Kaggle submissions held',flush=True);return
    prefix=repair_family if repair_family else 'crossfitA' if main_a else 'crossfit120'
    config=repair_family+'_records.json' if repair_family else 'incidental_weak.json' if main_a else 'crossfit120_records.json'
    if repair_family:
        gate=ROOT/('reports/cache120_main_probe_validation.json' if repair_family=='repairAb' else 'reports/cache_repair_probe_A_validation.json' if repair_family in ('repairA1','repairA2') else 'reports/cache_repair_probe_groups.json')
        while not gate.exists():time.sleep(45)
        audit=json.loads(gate.read_text())
        passed=(audit['passed'] if repair_family in ('repairA1','repairA2','repairAb') else audit['families']['repair120' if repair_family in ('repairBa','repairAa') else repair_family]['passed'])
        if not passed:
            print('Repair accuracy probe failed; geometry promotion stopped',flush=True);return
    if main_a:
        while not (ROOT/'configs/incidental_main.json').exists():time.sleep(45)
        from scripts.prepare_incidental_weak import prepare
        prepare()
    ledger=ROOT/('reports/'+prefix+'_controller.json')
    protocol=json.loads((ROOT/'configs'/config).read_text())
    geometry_stage=protocol.get('geometry_source_stage',prefix+'_geometry')
    state=json.loads(ledger.read_text()) if ledger.exists() else dict(launched=([] if main_a or repair_family else [prefix+'_geometry']),downloaded=[])
    def save():ledger.write_text(json.dumps(state,indent=2))
    def push(api,stage):
        if not (ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').exists():
            subprocess.run([sys.executable,'scripts/build_efficiency_kaggle.py',stage],cwd=ROOT,check=True)
        try:api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-').lower())
        except ValueError:
            api.kernels_status('malyshevdanil/wsindy-sampling-saturation-efficiency-match')
            response=api.kernels_push(str(ROOT/'kaggle'/stage))
            if getattr(response,'error',None):raise RuntimeError('Submission rejected')
        state['launched'].append(stage);save();print(stage+' submitted or adopted',flush=True)
    save()
    while len(state['downloaded'])<6:
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate()
            if geometry_stage not in state['launched']:
                if geometry_stage!=prefix+'_geometry':
                    api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+geometry_stage.replace('_','-').lower())
                    state['launched'].append(geometry_stage);save()
                else:push(api,geometry_stage)
            for stage in list(state['launched']):
                if stage in state['downloaded']:continue
                foreign_geometry=stage==geometry_stage and geometry_stage!=prefix+'_geometry'
                if foreign_geometry:
                    owner_ledger=ROOT/('reports/'+geometry_stage.split('_')[0]+'_controller.json')
                    if not owner_ledger.exists() or stage not in json.loads(owner_ledger.read_text()).get('downloaded',[]):continue
                response=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-').lower())
                value=response.get('status') if isinstance(response,dict) else response.status
                status=getattr(value,'name',str(value)).split('.')[-1].upper()
                if status=='CANCEL_ACKNOWLEDGED':
                    # Kaggle creates a kernel entry even when CPU quota rejects
                    # execution. Retry the identical frozen notebook later.
                    submitted=api.kernels_push(str(ROOT/'kaggle'/stage))
                    error=getattr(submitted,'error',None)
                    if error:
                        state['waiting_for_cpu_slot']=stage;save()
                        print(stage+' waiting for CPU slot',flush=True)
                    else:
                        state.pop('waiting_for_cpu_slot',None);save()
                        print(stage+' resubmitted unchanged snapshot',flush=True)
                    continue
                if status in ('ERROR','FAILED','CANCELED','CANCELLED'):
                    state['failure']=dict(stage=stage,status=status);save();return
                if status!='COMPLETE':continue
                from scripts.download_kaggle_archive import download
                if not foreign_geometry:download(stage)
                folder=ROOT/f'results/kaggle_{stage}';folder.mkdir(exist_ok=True)
                with zipfile.ZipFile(ROOT/f'results/kaggle_{stage}_download/wsindy_{stage}_outputs.zip') as z:
                    assert z.testzip() is None
                    assert all((folder/n).resolve().is_relative_to(folder.resolve()) for n in z.namelist())
                    z.extractall(folder)
                from scripts.validate_nuisance_runs import snapshot
                snapshot(stage,protocol.get('geometry_protocol_file',config) if stage==geometry_stage else config)
                assert json.loads((folder/'completion.json').read_text())['status']=='complete'
                state['downloaded'].append(stage);save();print(stage+' validated snapshot',flush=True)
            if geometry_stage in state['downloaded'] and prefix+'_records_smoke' not in state['launched']:
                folder=ROOT/('results/kaggle_'+geometry_stage)
                geometry=json.loads((folder/'geometry.json').read_text());assert geometry['cache_validated']
                cache=(folder/'crossfit_instruments.npz').read_bytes();geometry['cache_sha256']=hashlib.sha256(cache).hexdigest()
                (ROOT/('configs/'+prefix+'_instruments.npz')).write_bytes(cache)
                (ROOT/('configs/'+prefix+'_cache_geometry.json')).write_text(json.dumps(geometry,indent=2))
                old=json.loads((ROOT/'reports/heterogeneous_records_results.json').read_text())['geometry']
                info={k:old[k] for k in ('information_grid','information_x','information_v','information_interpolation_relative_errors')}
                (ROOT/'configs/crossfit_information.json').write_text(json.dumps(info,indent=2))
                push(api,prefix+'_records_smoke')
            if prefix+'_records_smoke' in state['downloaded']:
                from scripts.validate_crossfit_runs import validate
                if not state.get('smoke_validated'):
                    validate([prefix+'_records_smoke']);state['smoke_validated']=True;save()
                for i in range(4):
                    stage=f'{prefix}_records_{i}'
                    if stage not in state['launched']:push(api,stage)
        except (AssertionError,FileNotFoundError) as error:
            state['validation_failure']=type(error).__name__;save();return
        except Exception as error:print(json.dumps(dict(retry=True,error_type=type(error).__name__)),flush=True)
        if len(state['downloaded'])<6:time.sleep(45)
    from scripts.validate_crossfit_runs import validate
    validate([f'{prefix}_records_{i}' for i in range(4)],full=True)
    state['full_validation_complete']=True;save()
    if not main_a and not (repair_family and repair_family.startswith('repairA')):
        analyzer='scripts.analyze_crossfit_records' if repair_family=='repair480' else 'scripts.analyze_crossfit120'
        subprocess.run([sys.executable,'-m',analyzer]+(['--prefix',repair_family] if repair_family else []),cwd=ROOT,check=True)
        state['analysis_complete']=True;save()
        if repair_family=='repair480':
            for module in ('prepare_crossfit_theory','prepare_crossfit_results','clean_v13_manuscript','build_publication','build_ejs'):
                subprocess.run([sys.executable,'-m','scripts.'+module],cwd=ROOT,check=True)
            subprocess.run([sys.executable,'-m','scripts.validate_efficiency','--publication'],cwd=ROOT,check=True)
            state['analysis_and_manuscript_completed']=True;save()


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--main-a',action='store_true');a.add_argument('--repair-family',choices=['repair120','repair480','repairA0','repairA1','repairA2','repairAa','repairAb','repairBa'])
    args=a.parse_args();main(args.main_a,args.repair_family)
