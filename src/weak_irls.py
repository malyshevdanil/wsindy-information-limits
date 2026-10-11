"""WENDy Algorithm-2 reweighting core on a supplied support and fixed windows.

The complete multiscale test-function selection, variance filter and stopping
diagnostics of author WENDy are not reproduced. Two noise models are compared:
the iid residual covariance from the published core, and its known-OU extension.
"""
import numpy as np
from scipy.linalg import solve_triangular
from .weak_features import test_functions,quadrature_weights,derivative_polynomial
from .covariance import exponential_apply


def fit_irls(t,y,powers,support,centers,width,sigma,tau=None,iterations=10,tolerance=1e-6):
    from .weak_features import weak_matrices
    G,b=weak_matrices(t,y,powers,centers,width)
    K,d=b.shape;pairs=[(j,i) for i in range(d) for j in np.flatnonzero(support[:,i])]
    A=np.zeros((K*d,len(pairs)))
    for col,(j,i) in enumerate(pairs):A[np.arange(i,K*d,d),col]=G[:,j]
    response=b.ravel();q=np.linalg.lstsq(A,response,rcond=None)[0]
    phi,der=test_functions(t,centers,width);w=quadrature_weights(t)
    sigma=np.broadcast_to(sigma,(d,));h=t[1]-t[0]
    history=[]
    for index in range(iterations):
        coefficients=np.zeros_like(support,dtype=float)
        for value,(j,i) in zip(q,pairs):coefficients[j,i]=value
        C=np.zeros((K*d,K*d))
        for j in range(d):
            alpha=tuple(int(k==j) for k in range(d))
            jac=derivative_polynomial(y,powers,coefficients,alpha)
            c=-phi[:,:,None]*jac[None,:,:];c[:,:,j]-=der
            U=(c*w[None,:,None]).transpose(1,0,2).reshape(len(t),K*d)
            applied=U if tau is None else exponential_apply(U,np.exp(-h/tau))
            C+=sigma[j]**2*(U.T@applied)
        # Scale-aware covariance relaxation, explicitly different from author default.
        relaxation=1e-10*max(np.trace(C)/(K*d),np.finfo(float).tiny)
        chol=np.linalg.cholesky((C+C.T)/2+relaxation*np.eye(K*d))
        AW=solve_triangular(chol,A,lower=True)
        bw=solve_triangular(chol,response,lower=True)
        qnew=np.linalg.lstsq(AW,bw,rcond=None)[0]
        change=float(np.linalg.norm(qnew-q)/max(np.linalg.norm(q),1e-30))
        history.append(change);q=qnew
        if change<tolerance:break
    estimate=np.zeros_like(support,dtype=float)
    for value,(j,i) in zip(q,pairs):estimate[j,i]=value
    return estimate,{'iterations':len(history),'converged':history[-1]<tolerance,'last_change':history[-1]}
