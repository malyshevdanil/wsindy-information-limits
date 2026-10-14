"""Start only three declared, never-run studies as authorized CPU slots free.

Never restarts ERROR/CANCELED/COMPLETE studies or any legacy queue.
"""
import json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi
    from scripts.submit_endpoint_followup import main as submit
    from scripts.collect_boundary_followups import STAGES
    api=KaggleApi();api.authenticate();pending=['support_wald_v2_1','support_wald_v2_2','support_wald_v2_3']
    busy_candidates=STAGES+['support_wald_0','support_wald_1']
    while pending:
        states={}
        (ROOT/'reports/unstarted_followup_queue.json').write_text(json.dumps({'pending':pending,'legacy_queues_untouched':True},indent=2),encoding='utf-8')
        for stage in busy_candidates:
            try:
                r=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'))
                states[stage]=getattr(r.status,'name',str(r.status))
            except Exception as error:
                response=getattr(error,'response',None)
                absent=(response is not None and response.status_code==404) or (isinstance(error,ValueError) and str(error).startswith('Cannot access kernel'))
                states[stage]='NO_RUN' if absent else 'TRANSPORT_UNCERTAIN'
        pending=[s for s in pending if states.get(s) not in ('RUNNING','QUEUED','COMPLETE','ERROR','CANCELED','CANCEL_ACKNOWLEDGED')]
        busy=sum(v in ('RUNNING','QUEUED','TRANSPORT_UNCERTAIN') for v in states.values())
        if busy<5:
            for stage in list(pending):
                if states.get(stage)=='NO_RUN' and submit(stage,api):
                    pending.remove(stage);busy+=1
                    if busy>=5:break
        (ROOT/'reports/unstarted_followup_queue.json').write_text(json.dumps({'pending':pending,'states':states,'legacy_queues_untouched':True},indent=2),encoding='utf-8')
        print('Unstarted authorized studies remaining: '+', '.join(pending),flush=True)
        if pending:time.sleep(45)
    print('All declared studies have runs; no restart performed',flush=True)


if __name__=='__main__':main()
