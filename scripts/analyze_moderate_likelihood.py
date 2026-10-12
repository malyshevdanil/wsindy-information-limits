"""Pair moderate-noise likelihood coefficients with the completed 128-path cohort."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    folder=ROOT/'results/kaggle_moderate_likelihood';load=lambda p:json.loads(p.read_text(encoding='utf-8'))
    assert load(folder/'completion.json')['status']=='complete';rows=load(folder/'rows.json');p=load(folder/'manifest.json')['protocol'];truth=np.array([1,-1,1,-1])
    matched=[];ou=[]
    for i in range(4):matched+=load(ROOT/f'results/kaggle_efficiency_match_{i}/rows.json');ou+=load(ROOT/f'results/kaggle_efficiency_power_{i}/rows.json')
    # The 5% weak estimates are retained in the original 128-path final stages.
    ou=[]
    for i in range(2):ou+=load(ROOT/f'results/kaggle_efficiency_vdp_final_{i}/rows.json')
    cohorts={}
    for method in ('known_initial','unknown_initial'):
        g=sorted([r for r in rows if r['method']==method],key=lambda r:r['seed']);assert [r['seed'] for r in g]==list(range(128));cohorts[method]=g
    author={r['seed']:r for r in matched}
    assert set(author)==set(range(128))
    assert all(r['input_sha256']==author[r['seed']]['input_sha256'] for r in rows)
    rng=np.random.default_rng(3500000);ids=rng.integers(0,128,(5000,128));groups=[]
    for method,g in cohorts.items():
        q=np.array([r['coefficients'] for r in g]);e=(q-truth)/2;mean=e.mean(axis=0);squared=np.sum(e*e,axis=1)
        groups.append(dict(method=method,paths=128,relative_rmse=float(np.sqrt(squared.mean())),relative_bias=float(np.linalg.norm(mean)),
            relative_dispersion=float(np.sqrt(np.mean(np.sum((e-mean)**2,axis=1)))),rmse_bootstrap_95=np.quantile(np.sqrt(squared[ids].mean(axis=1)),[.025,.975]).tolist(),
            bound_active_paths=sum(any(r['bound_active']) for r in g),maximum_nfev=max(r['nfev'] for r in g),maximum_optimality=max(r['optimality'] for r in g)))
    aq=np.array([author[i]['coefficients'] for i in range(128)]);ae=np.sum(((aq-truth)/2)**2,axis=1);paired=[]
    for method,g in cohorts.items():
        e=np.sum(((np.array([r['coefficients'] for r in g])-truth)/2)**2,axis=1)
        paired.append(dict(method=method,author_minus_likelihood_rmse=float(np.sqrt(ae.mean())-np.sqrt(e.mean())),
            bootstrap_95=np.quantile(np.sqrt(ae[ids].mean(axis=1))-np.sqrt(e[ids].mean(axis=1)),[.025,.975]).tolist()))
    result=dict(completed=True,protocol=p,groups=groups,paired_author_comparisons=paired,observation_hashes_match_author_128=True,
        output_sha256=hashlib.sha256((ROOT/'results/kaggle_moderate_likelihood_download/wsindy_moderate_likelihood_outputs.zip').read_bytes()).hexdigest(),
        scope='bounded local path likelihood optimization; no global optimum guarantee')
    (ROOT/'reports/moderate_likelihood_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))


if __name__=='__main__':main()
