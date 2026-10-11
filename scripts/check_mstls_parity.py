"""Read-only parity check against a separately downloaded, pinned author file.

The author source is excluded from the publication archive. It is inspected and
used privately as a reference; no full WSINDy or WmSINDy replication is claimed.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.mstls import relative_mstls
from src.systems import system,trajectory
from src.weak_features import weak_matrices


def run():
    source=ROOT/'results/baseline_sources/utils_WSINDy.py'
    if not source.exists():
        import requests
        source.parent.mkdir(parents=True,exist_ok=True)
        url='https://raw.githubusercontent.com/CristianLopez85/WmSINDy/f1a781009a5c6adfd5650bec63ef30cc0673ae97/Fig32_Correlated%20Noise/utils_WSINDy.py'
        response=requests.get(url,timeout=(20,60));response.raise_for_status()
        if hashlib.sha256(response.content).hexdigest()!='8e5692ed15aea69cd82f6a78fda4bf0e84d5293fb227b4375485e4277aff1997':
            raise ValueError('Author reference fingerprint mismatch')
        source.write_bytes(response.content)
    spec=importlib.util.spec_from_file_location('author_wsindy_reference',source)
    reference=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    records=[]
    thresholds=np.logspace(-4,0,50)
    for name in ('linear','damped','vdp','duffing'):
        _,powers,truth=system(name)
        t=np.linspace(0,20,1001)
        x=trajectory(name,t)
        # Deterministic analytic perturbation, not an empirical noise experiment.
        y=x+.02*np.column_stack([np.sin((j+1)*3.17*t+.13*j) for j in range(x.shape[1])])
        G,b=weak_matrices(t,y,powers,np.linspace(1.505,18.495,30),3.)
        ours,choices=relative_mstls(G,b,thresholds)
        for equation in range(b.shape[1]):
            ref,loss=reference.wsindy_pde_RGLS_seq(thresholds,0.,
                np.column_stack([b[:,equation],G]),1,np.ones(G.shape[1]+1),.8)
            delta=float(np.max(np.abs(ours[:,equation]-ref[:,0])))
            assert np.allclose(ours[:,equation],ref[:,0],rtol=1e-11,atol=1e-12)
            assert abs(choices[equation]['threshold']-thresholds[np.argmin(loss[0])])<1e-12
            records.append({'system':name,'equation':equation,'max_coefficient_difference':delta})
    report={'scope':'relative MSTLS regression only, gamma=0, unit rescaling',
            'source_url':'https://github.com/CristianLopez85/WmSINDy/blob/f1a781009a5c6adfd5650bec63ef30cc0673ae97/Fig32_Correlated%20Noise/utils_WSINDy.py',
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'reference_redistributed':False,'fixtures':'analytic trajectories and deterministic perturbation',
            'records':records,'status':'passed'}
    (ROOT/'reports/mstls_parity.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))


if __name__=='__main__':run()
