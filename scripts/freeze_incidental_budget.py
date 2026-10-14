"""Freeze a user-directed budget amendment without reading main outcomes."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    original=ROOT/'configs/incidental_main.json';target=ROOT/'configs/incidental_trimmed.json'
    assert not target.exists(), 'Never overwrite a frozen protocol'
    p=json.loads(original.read_text());p.update(protocol='shared_coefficient_incidental_comparison_budget_v2',
        previous_protocol='incidental_main.json',bootstrap_correction_cells=[{'level':.2,'n':400},{'level':.2,'n':800}],
        revision_reason='User-directed allocation before inspection of main outcomes: retain B=50 correction only at 20% noise, n=400/800; other cells retain profile and oracle likelihood',
        status='frozen budget amendment; physical seeds and optimizer settings unchanged')
    target.write_text(json.dumps(p,indent=2))
    state=json.loads((ROOT/'reports/incidental_main_controller.json').read_text())
    report=dict(decision_time_utc=datetime.now(timezone.utc).isoformat(),
        original_protocol_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
        amended_protocol_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        evidence_used=['validated 50-repetition pilot','validated deterministic curvature','user budget instruction'],
        main_outcomes_inspected_for_decision=False,already_launched=state['launched'],
        already_downloaded=state['downloaded'],new_correction_cells=p['bootstrap_correction_cells'],
        new_planned_bootstrap_refits=2*200*50,original_planned_bootstrap_refits=3*5*200*50,
        existing_jobs='Do not retroactively relabel existing executions. Two already-running old 5%-noise jobs may finish; their corrected estimates are outside the amended comparison.',
        future_queue='No further full-bootstrap execution outside the two amended cells')
    (ROOT/'reports/incidental_budget_decision.json').write_text(json.dumps(report,indent=2))
    (ROOT/'reports/incidental_main_controller_pre_budget.json').write_text(json.dumps(state,indent=2))
    print('Budget amendment frozen without inspecting main outcomes')


if __name__=='__main__':main()
