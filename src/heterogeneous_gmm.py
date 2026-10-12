"""Four-dimensional independent-record estimating equations and studentization."""
import numpy as np


def interpolate_grid(grid,xs,ys,u):
    u=np.asarray(u);i=np.clip(np.searchsorted(xs,u[:,0],side='right')-1,0,len(xs)-2)
    j=np.clip(np.searchsorted(ys,u[:,1],side='right')-1,0,len(ys)-2)
    ax=np.clip((u[:,0]-np.asarray(xs)[i])/(np.asarray(xs)[i+1]-np.asarray(xs)[i]),0,1)
    ay=np.clip((u[:,1]-np.asarray(ys)[j])/(np.asarray(ys)[j+1]-np.asarray(ys)[j]),0,1)
    shape=(len(u),)+(1,)*(grid.ndim-2);ax=ax.reshape(shape);ay=ay.reshape(shape)
    return (1-ax)*(1-ay)*grid[i,j]+ax*(1-ay)*grid[i+1,j]+(1-ax)*ay*grid[i,j+1]+ax*ay*grid[i+1,j+1]


def record_fit(B,a):
    meanB=B.mean(axis=-3);meana=a.mean(axis=-2)
    q=np.linalg.solve(meanB,meana[...,None])[...,0]
    residual=a-np.einsum('...rij,...j->...ri',B,q)
    residual=residual-residual.mean(axis=-2,keepdims=True);n=len(B) if B.ndim==3 else B.shape[1]
    C=np.einsum('...ri,...rj->...ij',residual,residual)/(n-1)
    inverse=np.linalg.inv(meanB);V=np.einsum('...ik,...kl,...jl->...ij',inverse,C,inverse)/n
    return q,(V+V.swapaxes(-1,-2))/2


def record_bootstrap(B,a,q,rng,draws=999):
    n=len(B);pivots=[];joint=[]
    for start in range(0,draws,64):
        ids=rng.integers(0,n,(min(64,draws-start),n));qb,Vb=record_fit(B[ids],a[ids]);error=qb-q
        if not np.isfinite(qb).all() or np.linalg.eigvalsh(Vb).min()<=0:raise RuntimeError('Invalid resample; no discarded fits')
        pivots.append(error/np.sqrt(np.diagonal(Vb,axis1=-2,axis2=-1)))
        joint.append(np.einsum('bi,bi->b',error,np.linalg.solve(Vb,error[...,None])[...,0]))
    return np.quantile(np.concatenate(pivots),[.025,.975],axis=0),float(np.quantile(np.concatenate(joint),.95))
