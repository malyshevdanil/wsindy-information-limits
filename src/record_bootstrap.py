"""Recentered pair-cluster bootstrap for the two-step averaged-moment estimator."""
import numpy as np


def batch_fit(raw,b,vs):
    """The scalar fit_mean algorithm, batched without altering its regularization."""
    B=raw.shape[1];v=vs.mean(axis=1)[:,0]
    A=raw.copy();A[:,:,1::2,3]-=v[:,None,None]*A[:,:,1::2,2]
    meanA=A.mean(axis=1);meanb=b.mean(axis=1)
    U,s,Vh=np.linalg.svd(meanA,full_matrices=False)
    q0=np.einsum('bij,bj->bi',Vh.transpose(0,2,1),np.einsum('bmi,bm->bi',U,meanb)/s)

    def influences(q):
        r=b-np.einsum('brmp,bp->brm',A,q)
        H=np.zeros(meanb.shape);H[:,1::2]=q[:,3,None]*meanA[:,1::2,2]
        return r+(vs[:,:,0]-v[:,None])[:,:,None]*H[:,None,:]

    def cov(r):
        r=r-r.mean(axis=1,keepdims=True)
        return np.einsum('brm,brk->bmk',r,r)/(B-1)

    Omega=cov(influences(q0));regularized=Omega.copy()
    ix=np.arange(Omega.shape[1]);regularized[:,ix,ix]+=Omega[:,ix,ix]/B
    WA=np.linalg.solve(regularized,meanA)
    P=np.linalg.solve(np.einsum('bmp,bmq->bpq',meanA,WA),WA.transpose(0,2,1))
    q=np.einsum('bpm,bm->bp',P,meanb)
    V=np.einsum('bpm,bmk,bqk->bpq',P,cov(influences(q)),P)/B
    return q,(V+V.transpose(0,2,1))/2


def bootstrap_t(raw,b,vs,q,V,rng,draws=999,batch_size=64):
    """Resample whole independent pairs, impose the fitted moment null, studentize."""
    B=len(raw);v=vs[:,0].mean();A=raw.copy();A[:,1::2,3]-=v*A[:,1::2,2]
    offset=b.mean(axis=0)-A.mean(axis=0)@q
    marginal=[];joint=[];minimum_eigenvalue=np.inf
    for start in range(0,draws,batch_size):
        indices=rng.integers(0,B,size=(min(batch_size,draws-start),B))
        qb,Vb=batch_fit(raw[indices],b[indices]-offset,vs[indices])
        error=qb-q;sd=np.sqrt(np.diagonal(Vb,axis1=1,axis2=2))
        eigenvalues=np.linalg.eigvalsh(Vb);minimum_eigenvalue=min(minimum_eigenvalue,float(eigenvalues.min()))
        if eigenvalues.min()<=0 or not np.isfinite(qb).all():raise RuntimeError('Invalid bootstrap fit; no draws discarded')
        marginal.append(error/sd)
        joint.append(np.einsum('bi,bi->b',error,np.linalg.solve(Vb,error[...,None])[...,0]))
    marginal=np.concatenate(marginal);joint=np.concatenate(joint)
    return dict(marginal_t_quantiles=np.quantile(marginal,[.025,.975],axis=0).tolist(),
        joint_critical_value=float(np.quantile(joint,.95)),bootstrap_draws=draws,
        minimum_bootstrap_covariance_eigenvalue=minimum_eigenvalue)
