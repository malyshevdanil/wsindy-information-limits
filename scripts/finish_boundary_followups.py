"""Read-only completion; analyze each complete study without inspecting partial samples."""
import json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from scripts.collect_boundary_followups import main as collect
    from scripts.analyze_boundary_followups import general,wald
    done=set()
    while True:
        try:
            collect()
        except Exception as error:
            print('Collection transport failure: '+type(error).__name__+'; retrying read-only',flush=True)
            time.sleep(45)
            continue
        records=json.loads((ROOT/'reports/boundary_followup_status.json').read_text())
        complete={r['stage'] for r in records if r.get('completion',{}).get('status')=='complete'}
        for system in ['vdp5','lorenz']:
            if system not in done and 'endpoint_'+system in complete:general(system);done.add(system)
        if 'wald' not in done and all(f'support_wald_v2_{i}' in complete for i in range(4)):wald();done.add('wald')
        if len(done)==3:print('All full-sample studies analyzed',flush=True);return
        if any(r.get('status') in ('ERROR','CANCELED','CANCEL_ACKNOWLEDGED') for r in records):
            print('Failed/cancelled study preserved; no automatic resubmission',flush=True);return
        time.sleep(45)


if __name__=='__main__':main()
