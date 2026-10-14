"""Wait for all matching components before running the four-method analysis."""
import json,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    while True:
        a=ROOT/'reports/incidental_main_controller.json'
        compact=ROOT/'configs/incidental_compact.json'
        families=json.loads(compact.read_text())['weak_families'] if compact.exists() else ['repairA0','repairA1','repairA2'] if (ROOT/'configs/repairA0_records.json').exists() else ['crossfitA']
        try:
            sa=json.loads(a.read_text()) if a.exists() else {}
            states=[json.loads(path.read_text()) if path.exists() else {} for path in [ROOT/('reports/'+family+'_controller.json') for family in families]]
        except (OSError,json.JSONDecodeError):
            time.sleep(45);continue
        if sa.get('failure') or sa.get('validation_failure') or any(s.get('failure') or s.get('validation_failure') for s in states):
            print('A component failed; combined analysis stopped',flush=True);return
        if sa.get('likelihood_components_complete') and all(s.get('full_validation_complete') for s in states):break
        time.sleep(45)
    subprocess.run([sys.executable,'-m','scripts.analyze_incidental_main'],cwd=ROOT,check=True)


if __name__=='__main__':main()
