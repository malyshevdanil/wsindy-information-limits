"""Matched separation of finite-span, fission and feasible-pilot errors."""
import argparse,json
from pathlib import Path
import numpy as np
from scripts.validate_crossfit_runs import validate
from scripts.analyze_independent_records import wilson
ROOT=Path(__file__).resolve().parents[1]


def main(prefix='crossfit120'):
    validate([f'{prefix}_records_{i}' for i in range(4)],full=True)
    p=json.loads((ROOT/f'configs/{prefix}_records.json').read_text());R=p['monte_carlo_repetitions']
    rows=[r for i in range(4) for r in json.loads((ROOT/f'results/kaggle_{prefix}_records_{i}/rows.json').read_text())]
    truth=np.asarray(p['truth']);rng=np.random.default_rng(4500000);groups=[];paired=[]
    for level in p['levels']:
        for n in p['record_counts']:
            errors={}
            for method in ('oracle120_original','oracle120_fission','crossfit120'):
                g=sorted([r for r in rows if r['method']==method and r['n']==n and r['level']==level],key=lambda r:r['repetition'])
                assert [r['repetition'] for r in g]==list(range(R))
                e=(np.array([r['coefficients'] for r in g])-truth)/2;mse=np.sum(e*e,axis=1);errors[method]=mse
                bias=e.mean(axis=0);ids=rng.integers(0,R,(5000,R));count=sum(r['wald_joint_covered'] for r in g)
                group=dict(level=level,n=n,method=method,repetitions=R,relative_rmse=float(np.sqrt(mse.mean())),
                    sqrt_n_rmse=float(np.sqrt(n*mse.mean())),relative_bias=float(np.linalg.norm(bias)),
                    coefficient_bias=(2*bias).tolist(),relative_dispersion=float(np.sqrt(np.mean(np.sum((e-bias)**2,axis=1)))),
                    sqrt_n_reference=float(np.sqrt(n)*np.mean([r['local_information_relative_sd'] for r in g])),
                    rmse_mc_95=np.quantile(np.sqrt(mse[ids].mean(axis=1)),[.025,.975]).tolist(),
                    wald_joint_coverage=count/R,wald_joint_wilson_95=wilson(count,R))
                if method=='crossfit120':
                    count=sum(r['refit_joint_covered'] for r in g);marginal=np.array([r['refit_marginal_covered'] for r in g])
                    group.update(bootstrap_joint_coverage=count/R,bootstrap_joint_wilson_95=wilson(count,R),
                        bootstrap_marginal_coverage=marginal.mean(axis=0).tolist(),
                        bootstrap_marginal_wilson_95=[wilson(int(k),R) for k in marginal.sum(axis=0)])
                assert np.isclose(group['relative_rmse']**2,group['relative_bias']**2+group['relative_dispersion']**2)
                groups.append(group)
            ids=rng.integers(0,R,(5000,R))
            for name,a,b in (('fission','oracle120_fission','oracle120_original'),('feasible_instruments','crossfit120','oracle120_fission')):
                delta=np.sqrt(errors[a][ids].mean(axis=1))-np.sqrt(errors[b][ids].mean(axis=1))
                paired.append(dict(level=level,n=n,component=name,rmse_difference=float(np.sqrt(errors[a].mean())-np.sqrt(errors[b].mean())),
                    paired_mc_95=np.quantile(delta,[.025,.975]).tolist()))
    report=dict(completed=True,protocol=p,groups=groups,paired_decomposition=paired,
        caveat='Deterministic oracle maps use the same 25-point initial-state interpolation; fission comparison isolates the two-view moments within that shared numerical family. Feasible instruments also use noisy initial pilots and cross-fitted coefficient pilots.')
    (ROOT/'reports/crossfit120_results.json').write_text(json.dumps(report,indent=2))
    print('Matched 120-test and fission decomposition complete')


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--prefix',default='crossfit120');main(a.parse_args().prefix)
