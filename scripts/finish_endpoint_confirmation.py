"""Finish the five authorized studies by read-only polling and local analysis.

No push, cancellation, retry submission or legacy queue activation is performed.
"""
import json, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))


def main():
    from scripts.collect_endpoint_confirmation import main as collect
    from scripts.analyze_endpoint_confirmation import main as boundary
    from scripts.report_endpoint_confirmation import main as boundary_report
    from scripts.analyze_support_extra_term import main as support
    print('Read-only finisher; no notebook submission or resubmission', flush=True)
    boundary_done = False; support_done = False
    while True:
        collect()
        statuses = json.loads((ROOT/'reports/endpoint_confirmation_status.json').read_text())
        complete = {r['stage'] for r in statuses if r.get('status') == 'COMPLETE' and r.get('completion', {}).get('status') == 'complete'}
        if not boundary_done and all(f'boundary_confirm_{i}' in complete for i in range(4)):
            boundary(); boundary_report(); boundary_done = True
        if not support_done and 'support_extra' in complete:
            support(); support_done = True
        if boundary_done and support_done:
            print('All studies complete; full paired and inclusion reports written', flush=True)
            return
        if any(r.get('status') in ('ERROR', 'CANCEL_ACKNOWLEDGED', 'CANCELED', 'CANCELLED') for r in statuses):
            print('Stopped on a failed or cancelled study; preserved available outputs', flush=True)
            return
        time.sleep(45)


if __name__ == '__main__': main()
