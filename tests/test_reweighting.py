"""Deterministic invariance checks for the covariance-weighted estimating step."""
import numpy as np
from src.systems import system,trajectory
from src.weak_irls import fit_irls


def test_irls_common_variance_scale_invariance():
    t=np.linspace(0,4,201);_,powers,truth=system('linear')
    y=trajectory('linear',t)+.01*np.sin(3.17*t)[:,None]
    centers=np.linspace(.805,3.195,12)
    for tau in (None,.4):
        first,meta=fit_irls(t,y,powers,truth!=0,centers,1.6,[.01],tau)
        second,_=fit_irls(t,y,powers,truth!=0,centers,1.6,[.05],tau)
        np.testing.assert_allclose(first,second,rtol=1e-10,atol=1e-12)
        assert np.array_equal(first!=0,truth!=0)
        assert meta['iterations']>=1
