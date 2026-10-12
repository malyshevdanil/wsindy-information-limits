"""Resume downloads of the four immutable bootstrap jobs."""
import concurrent.futures,json,subprocess,sys,time
from pathlib import Path
from scripts.finish_record_runs import status
ROOT=Path(__file__).resolve().parents[1]


def main():
    stages=[f'record_bootstrap_{i}' for i in range(4)];state=ROOT/'reports/bootstrap_run_controller.json'
    done=json.loads(state.read_text()) if state.exists() else [];last={}
    while len(done)<len(stages):
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            states=dict(pool.map(status,[s for s in stages if s not in done]))
        if states!=last:print(json.dumps(states),flush=True);last=states
        for stage,value in states.items():
            if value=='COMPLETE':
                for script in ('download_kaggle_archive','unpack_kaggle_outputs'):
                    subprocess.run([sys.executable,'scripts/'+script+'.py',stage],cwd=ROOT,check=True)
                done.append(stage);state.write_text(json.dumps(done,indent=2))
            elif value in ('ERROR','FAILED','CANCELED','CANCELLED'):raise RuntimeError('Bootstrap stage failed: '+stage)
        if len(done)<len(stages):time.sleep(45)
    print('All bootstrap results downloaded',flush=True)


if __name__=='__main__':main()
