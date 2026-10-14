"""Check signed pilot bias against the expansion, using no main outcomes."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    from scripts.validate_nuisance_runs import snapshot
    records=[];audits=[]
    for i in range(2):
        stage=f'incidental_pilot_v2_{i}'
        audits.append(snapshot(stage,'incidental_pilot_v2.json'))
        records+=json.loads((ROOT/f'results/kaggle_{stage}/rows.json').read_text())
    records=sorted([r for r in records if r['method']=='profile_likelihood'],key=lambda r:r['repetition'])
    assert [r['repetition'] for r in records]==list(range(50))
    truth=np.array(json.loads((ROOT/'configs/incidental_pilot_v2.json').read_text())['truth'])
    a=np.array(json.loads((ROOT/'reports/incidental_curvature_validation.json').read_text())['results']['results'][-1]['bias_coefficient'])
    e=np.array([r['coefficients'] for r in records])-truth;observed=e.mean(axis=0);prediction=.05**2*a
    V=np.cov(e,rowvar=False,ddof=1)/len(e);Vi=np.linalg.inv(V);delta=observed-prediction
    rng=np.random.default_rng(4700000);ids=rng.integers(0,len(e),(20000,len(e)))
    boot=e[ids].mean(axis=1);centered=boot-observed
    stat=float(delta@Vi@delta);null=np.einsum('bi,ij,bj->b',centered,Vi,centered)
    direction=a/np.linalg.norm(a)
    report=dict(completed=True,pilot_only=True,main_outcomes_used=False,execution_audits=audits,
        predicted_coefficient_bias=prediction.tolist(),observed_coefficient_bias=observed.tolist(),
        predicted_relative_bias=float(np.linalg.norm(prediction)/2),observed_relative_bias=float(np.linalg.norm(observed)/2),
        bias_norm_ratio_descriptive=float(np.linalg.norm(observed)/np.linalg.norm(prediction)),
        vector_discrepancy_statistic=stat,centered_bootstrap_pvalue=float((1+np.sum(null>=stat))/(len(null)+1)),
        coefficient_bias_mc_95=np.quantile(boot,[.025,.975],axis=0).tolist(),
        predicted_signed_projection=float(prediction@direction),observed_signed_projection=float(observed@direction),
        signed_projection_mc_95=np.quantile(boot@direction,[.025,.975]).tolist(),
        caveat='Pilot has 50 repetitions; a norm of an estimated mean is upward biased. Agreement here does not validate the fixed-noise or finite-n remainder. Repeat the vector check on the amended main comparison.')
    finite=ROOT/'reports/incidental_finite_n_bias.json'
    if finite.exists():
        expansion=json.loads(finite.read_text());assert expansion['completed']
        corrected=next(r for r in expansion['rows'][-1]['predictions'] if r['n']==800 and r['level']==.05)
        prediction_n=np.array(corrected['coefficient_bias']);delta_n=observed-prediction_n;stat_n=float(delta_n@Vi@delta_n)
        report['finite_n_comparison']=dict(predicted_coefficient_bias=prediction_n.tolist(),predicted_relative_bias=float(np.linalg.norm(prediction_n)/2),
            relative_change_of_bias_norm=float(np.linalg.norm(prediction_n)/np.linalg.norm(prediction)-1),
            vector_discrepancy_statistic=stat_n,centered_bootstrap_pvalue=float((1+np.sum(null>=stat_n))/(len(null)+1)),
            conclusion='The computed 1/n correction is too small to explain the pilot discrepancy; finite-noise remainders and Monte Carlo fluctuation remain unresolved.')
    (ROOT/'reports/pilot_curvature_bias_check.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:report[k] for k in ('bias_norm_ratio_descriptive','centered_bootstrap_pvalue','predicted_signed_projection','observed_signed_projection','signed_projection_mc_95')}))


if __name__=='__main__':main()
