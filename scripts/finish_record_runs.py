"""Download matched fits and schedule replicated-record jobs as slots free."""
import concurrent.futures,json,subprocess,sys,time
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
ROOT=Path(__file__).resolve().parents[1]
OLD=[f'efficiency_match_{i}' for i in range(4)]
NEW=[f'independent_records_{i}' for i in range(4)]
STATE=ROOT/'reports/record_run_controller.json'


def status(stage):
    api=KaggleApi();api.authenticate()
    try:
        r=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
        v=r.get('status') if isinstance(r,dict) else r.status
        return stage,getattr(v,'name',str(v)).split('.')[-1].upper()
    except Exception as e:return stage,'READ_ERROR:'+type(e).__name__


def call(script,stage):subprocess.run([sys.executable,'scripts/'+script+'.py',stage],cwd=ROOT,check=True)


def main():
    state=json.loads(STATE.read_text()) if STATE.exists() else dict(launched=[],downloaded=[])
    last={}
    while len(state['downloaded'])<len(OLD+NEW):
        check=[s for s in OLD+state['launched'] if s not in state['downloaded']]
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:states=dict(pool.map(status,check))
        if states!=last:print(json.dumps(dict(statuses=states)),flush=True);last=states
        for stage,value in states.items():
            if value=='COMPLETE':
                call('download_kaggle_archive',stage);call('unpack_kaggle_outputs',stage);state['downloaded'].append(stage)
            elif value in ('ERROR','FAILED','CANCELED','CANCELLED'):raise RuntimeError('Required stage failed: '+stage)
        # The fifth slot remains occupied by the original unsplit author job.
        available=sum(s in state['downloaded'] for s in OLD)-len(state['launched'])+sum(s in state['downloaded'] for s in NEW)
        for stage in [s for s in NEW if s not in state['launched']][:max(available,0)]:
            api=KaggleApi();api.authenticate();r=api.kernels_push(str(ROOT/'kaggle'/stage))
            error=r.get('error') if isinstance(r,dict) else getattr(r,'error',None)
            if error:
                print(json.dumps(dict(stage=stage,push_rejected=str(error))),flush=True);break
            state['launched'].append(stage);print(json.dumps(dict(stage=stage,launched=True)),flush=True)
        STATE.write_text(json.dumps(state,indent=2),encoding='utf-8')
        if len(state['downloaded'])<len(OLD+NEW):time.sleep(45)
    print(json.dumps(dict(all_matched_and_record_runs_downloaded=True)),flush=True)


if __name__=='__main__':main()
