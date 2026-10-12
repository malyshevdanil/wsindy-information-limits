"""Finish already authorized computations and build reviewable publication files."""
import concurrent.futures,json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STAGES=[*[f'efficiency_match_{i}' for i in range(4)],*[f'independent_records_{i}' for i in range(4)]]
def status(stage):
    from kaggle.api.kaggle_api_extended import KaggleApi
    try:
        api=KaggleApi();api.authenticate();r=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
        value=r.get('status') if isinstance(r,dict) else r.status
        return stage,getattr(value,'name',str(value)).split('.')[-1]
    except Exception as e:return stage,'read_error:'+type(e).__name__
def call(*args):subprocess.run([sys.executable,*args],cwd=ROOT,check=True)
def main():
    downloaded=set();last={}
    while len(downloaded)<len(STAGES):
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:states=dict(pool.map(status,[s for s in STAGES if s not in downloaded]))
        if states!=last:print(json.dumps({'statuses':states}),flush=True);last=states
        for stage,state in states.items():
            if state.upper()=='COMPLETE':
                call('scripts/download_kaggle_archive.py',stage);call('scripts/unpack_kaggle_outputs.py',stage);downloaded.add(stage)
            elif state.upper() in ['ERROR','FAILED','CANCELLED','CANCELED']:raise RuntimeError('Required author stage did not complete: '+stage)
        if len(downloaded)<len(STAGES):time.sleep(30)
    for module in ['prepare_efficiency_manuscript','analyze_efficiency','analyze_efficiency_followup','analyze_information_followup',
                   'prepare_efficiency_followup_manuscript','analyze_independent_records','prepare_independent_records_manuscript',
                   'validate_efficiency_followup','validate_independent_records','build_publication']:
        call('-m','scripts.'+module)
    call('scripts/validate_artifacts.py');call('scripts/validate_efficiency.py','--publication');call('scripts/validate_reviewer.py','--publication')
    call('scripts/package_publication.py');call('scripts/validate_full_adaptive.py','--publication')
    call('scripts/reproduce_efficiency_release.py');call('scripts/package_publication.py')
    print(json.dumps({'completed':True,'built_and_reproduced':True,'public_push_pending':True}),flush=True)
if __name__=='__main__':main()
