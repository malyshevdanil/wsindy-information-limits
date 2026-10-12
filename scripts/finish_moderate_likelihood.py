"""Launch the frozen likelihood check when a required heterogeneous job frees a slot."""
import json,subprocess,sys,time
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
from scripts.finish_record_runs import status
ROOT=Path(__file__).resolve().parents[1]


def main():
    marker=ROOT/'reports/moderate_likelihood_controller.json';state=json.loads(marker.read_text()) if marker.exists() else dict(launched=False)
    while not state['launched']:
        available=any(status(f'heterogeneous_records_{i}')[1]=='COMPLETE' for i in range(4))
        if available:
            api=KaggleApi();api.authenticate();r=api.kernels_push(str(ROOT/'kaggle/moderate_likelihood'));error=r.get('error') if isinstance(r,dict) else getattr(r,'error',None)
            if not error:state['launched']=True;marker.write_text(json.dumps(state));print('Moderate likelihood launched',flush=True)
        if not state['launched']:time.sleep(45)
    last=''
    while True:
        _,value=status('moderate_likelihood')
        if value!=last:print('Moderate likelihood '+value,flush=True);last=value
        if value=='COMPLETE':
            for script in ('download_kaggle_archive','unpack_kaggle_outputs'):subprocess.run([sys.executable,'scripts/'+script+'.py','moderate_likelihood'],cwd=ROOT,check=True)
            state['downloaded']=True;marker.write_text(json.dumps(state));break
        if value in ('FAILED','ERROR','CANCELED','CANCELLED'):raise RuntimeError('Moderate likelihood failed')
        time.sleep(45)


if __name__=='__main__':main()
