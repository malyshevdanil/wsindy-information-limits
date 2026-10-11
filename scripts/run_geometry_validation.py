"""Deterministic information check, executed on Kaggle rather than locally."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.run_experiments import source_hash
from src.weak_features import test_functions,quadrature_weights
from src.covariance import exponential_apply


def run(out):
    out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    cfg=json.loads((ROOT/'configs/geometry.json').read_text())
    a,x0,sigma,tau,T=[cfg[k] for k in ('a','x0','sigma','tau','T')]
    def integral(function):
        return quad(function,0,T,epsabs=cfg['quad_epsabs'],epsrel=cfg['quad_epsrel'])[0]
    # Columns correspond to coefficient a and initial amplitude x0.
    coefficient=lambda t:x0*np.exp(a*t)*(1+(a+1/tau)*t)
    initial=lambda t:(a+1/tau)*np.exp(a*t)
    factor=tau/(2*sigma*sigma)
    Iaa=factor*integral(lambda t:coefficient(t)**2)
    Iau=factor*integral(lambda t:coefficient(t)*initial(t))
    Iuu=1/(sigma*sigma)+factor*integral(lambda t:initial(t)**2)
    effective=Iaa-Iau*Iau/Iuu
    n=round(T/cfg['reference_h']);t=np.linspace(0,T,n+1)
    phi,der=test_functions(t,[cfg['window_center']],cfg['window_width'])
    weights=quadrature_weights(t);x=x0*np.exp(a*t);s=t*x;q=x/x0
    D=float((phi[0]*weights)@x)
    ell=(-der[0]-a*phi[0])/D;weighted=ell*weights
    variance=float(sigma*sigma*weighted@exponential_apply(weighted,np.exp(-(t[1]-t[0])/tau)))
    annihilation=float(weighted@q);normalization=float(weighted@s)
    assert abs(annihilation)<1e-9 and abs(normalization-1)<1e-9
    assert variance*effective>=1-1e-8
    result={'parameters':cfg,'joint_information':[[Iaa,Iau],[Iau,Iuu]],
            'coefficient_information_known_initial':Iaa,
            'coefficient_information_initial_nuisance':effective,
            'known_initial_variance_bound':1/Iaa,
            'weak_initial_annihilation_variance_bound':1/effective,
            'weak_bound_inflation':Iaa/effective,
            'weak_linearized_variance_reference_grid':variance,
            'weak_known_initial_information_efficiency':1/(variance*Iaa),
            'weak_nuisance_information_efficiency':1/(variance*effective),
            'initial_sensitivity_annihilation':annihilation,
            'coefficient_sensitivity_normalization':normalization,
            'scope':'Local Gaussian influence covariance; not variance of unguarded nonlinear ratio'}
    (out/'geometry.json').write_text(json.dumps(result,indent=2))
    (out/'manifest.json').write_text(json.dumps({'stage':'geometry','protocol':cfg,'source_sha256':source_hash()},indent=2))
    (out/'completion.json').write_text(json.dumps({'status':'complete','elapsed_seconds':time.perf_counter()-start},indent=2))
    print(json.dumps({'status':'complete','known_initial_information':Iaa,'weak_class_information_bound':effective}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
