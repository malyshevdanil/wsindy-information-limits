"""Run only the frozen independent-validation queue; never restart failed jobs."""
import json,sys,time,zipfile,hashlib
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.submit_endpoint_followup import main as submit
from scripts.download_kaggle_archive import download
from scripts.collect_stopped_runs import snapshot

def collect(stage,api):
    folder=ROOT/'results'/stage
    if not (folder/'completion.json').exists():
        download(stage,api=api);folder.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(ROOT/'results'/f'kaggle_{stage}_download'/f'wsindy_{stage}_outputs.zip') as archive:
            assert archive.testzip() is None
            assert all((folder/name).resolve().is_relative_to(folder.resolve()) for name in archive.namelist())
            archive.extractall(folder)
    manifest=json.loads((folder/'manifest.json').read_text());proof=snapshot(stage,manifest)
    assert manifest['protocol']==json.loads((ROOT/'configs/independent_validation_v2.json').read_text())
    completion=json.loads((folder/'completion.json').read_text());assert completion['status']=='complete'
    if stage.startswith('iq_smoke'):
        assert completion['rows']==2
        assert len(json.loads((folder/'geometry_rows.json').read_text()))==2
        assert len(json.loads(next(folder.glob('bootstrap_c*.json')).read_text())['targets'])==4
    elif stage.startswith('iq_geo_'):assert completion['rows']==200
    elif stage.startswith('iq_nom_'):assert completion['rows']==1000
    else:
        assert completion['rows']==40
        boot=list(folder.glob('bootstrap_c*.json'));assert len(boot)==20
        assert all(len(json.loads(path.read_text())['targets'])==399 for path in boot)
    return {'source_audit':proof,'completion':completion}

def main():
    # The SDK otherwise sends API requests without a timeout. Scope this
    # timeout to this standalone controller process, including archive reads.
    import requests
    original_send=requests.Session.send
    def bounded_send(session,request,**kwargs):
        kwargs.setdefault('timeout',(15,60))
        return original_send(session,request,**kwargs)
    requests.Session.send=bounded_send
    for attempt in range(10):
        try:
            from kaggle.api.kaggle_api_extended import KaggleApi
            api=KaggleApi();api.authenticate();break
        except Exception as error:
            print('startup authentication retry '+str(attempt+1)+': '+type(error).__name__,flush=True)
            if attempt==9:raise RuntimeError('Kaggle authentication unavailable after retries') from None
            time.sleep(30)
    protocol=ROOT/'configs/independent_validation_v2.json'
    decision=json.loads((ROOT/'reports/independent_validation_decision_v2.json').read_text())
    assert hashlib.sha256(protocol.read_bytes()).hexdigest()==decision['protocol_sha256']
    smoke=collect('iq_smoke_v2',api)
    (ROOT/'reports/independent_smoke_v2_validation.json').write_text(json.dumps(smoke,indent=2))
    stages=[f'iq_nom_{i}' for i in range(4)]+[f'iq_geo_{i}' for i in range(5)]+[f'iq_bootn_{i:02d}' for i in range(50)]+[f'iq_bootp_{i:02d}' for i in range(10)]
    path=ROOT/'reports/independent_validation_controller.json'
    state=json.loads(path.read_text()) if path.exists() else {'stages':{s:{'status':'NOT_SUBMITTED'} for s in stages},'queue_complete':False}
    def save():
        state['updated_at_utc']=datetime.now(timezone.utc).isoformat();path.write_text(json.dumps(state,indent=2))
    save()
    while True:
        active=0
        for stage,record in state['stages'].items():
            if record['status'] in ['VERIFIED_COMPLETE','HALTED_ERROR']:continue
            if record['status']=='NOT_SUBMITTED':continue
            slug='malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-')
            try:
                print(stage+': checking remote state',flush=True)
                response=api.kernels_status(slug);status=getattr(response.status,'name',str(response.status))
                if status=='COMPLETE':record.update(status='VERIFIED_COMPLETE',**collect(stage,api))
                elif status in ['ERROR','CANCELLED','CANCELED']:
                    record.update(status='HALTED_ERROR',remote_status=status);state['halted']=True;save();return
                else:record['status']=status;active+=1
                save()
            except Exception as error:
                response=getattr(error,'response',None)
                absent=(isinstance(error,ValueError) and str(error).startswith('Cannot access kernel')) or (response is not None and response.status_code==404)
                if absent and record['status']=='RECONCILE':
                    record['status']='NOT_SUBMITTED';continue
                record['last_transport_error']=type(error).__name__
                active+=record['status']!='NOT_SUBMITTED'
        save()
        if all(r['status']=='VERIFIED_COMPLETE' for r in state['stages'].values()):
            from scripts.analyze_independent_validation import analyze
            analyze();state['queue_complete']=True;save();return
        if not state.get('auto_submission_enabled',True):
            if active==0:
                from scripts.analyze_independent_validation import analyze
                analyze();state['collection_complete_for_submitted']=True;save();return
            print(json.dumps({'mode':'collect_only','active':active,'new_submissions':False}),flush=True)
            time.sleep(45);continue
        for stage,record in state['stages'].items():
            if active>=5:break
            if record['status']=='NOT_SUBMITTED':
                result=submit(stage,api=api)
                if result is True:record['status']='SUBMITTED';active+=1;save()
                elif result is None:
                    # An accepted push can lose its response, or this slug may
                    # already exist. Reconcile before allocating another slot.
                    record['status']='RECONCILE';save();break
                else:
                    # Includes CPU-capacity rejection: do not try every other
                    # stage against the same full pool in this polling cycle.
                    save();break
        print(json.dumps({'active':active,'verified':sum(r['status']=='VERIFIED_COMPLETE' for r in state['stages'].values()),'total':len(stages)}),flush=True)
        time.sleep(45)

if __name__=='__main__':main()
