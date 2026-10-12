"""Averaged weak moments with independent-pair nuisance influence correction."""
import numpy as np
from scipy.linalg import cho_factor,cho_solve


def corrected_design(raw,variance_x):
    A=raw.copy()
    A[...,1::2,3]-=variance_x*A[...,1::2,2]
    return A


def covariance(values):
    centered=values-values.mean(axis=0)
    return centered.T@centered/(len(values)-1)


def pair_influences(A,b,q,variance_statistics,variance_x,estimate_variance):
    residual=b-np.einsum('bmp,p->bm',A,q)
    if estimate_variance:
        # d E[b-Gq]/d variance_x = q_d E[integral phi*y_2].
        H=np.zeros(A.shape[1]);H[1::2]=q[3]*A.mean(axis=0)[1::2,2]
        residual=residual+(variance_statistics[:,0]-variance_x)[:,None]*H
    return residual


def fit_mean(raw_pair,b_pair,variance_statistics,variance,estimate_variance):
    """One two-step GMM update; all randomness is clustered by record pair."""
    B=len(raw_pair);A=corrected_design(raw_pair,variance[0]);meanA=A.mean(axis=0);meanb=b_pair.mean(axis=0)
    q0=np.linalg.lstsq(meanA,meanb,rcond=None)[0]
    psi=pair_influences(A,b_pair,q0,variance_statistics,variance[0],estimate_variance)
    Omega=covariance(psi)
    regularized=Omega+np.diag(np.diag(Omega))/B
    chol=cho_factor(regularized,lower=True,check_finite=False)
    WA=cho_solve(chol,meanA,check_finite=False)
    P=np.linalg.solve(meanA.T@WA,WA.T);q=P@meanb
    final=pair_influences(A,b_pair,q,variance_statistics,variance[0],estimate_variance)
    V=P@covariance(final)@P.T/B
    return q,(V+V.T)/2,dict(weight_condition=float(np.linalg.cond(regularized)),design_condition=float(np.linalg.cond(meanA)))
