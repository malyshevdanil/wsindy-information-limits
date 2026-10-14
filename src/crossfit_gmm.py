"""Record-wise cross-fitted pilots and bootstrap with pilot refitting.

Instrument caches use either an affine map or a supplied bounded polynomial
basis in the coefficient pilot. Records are assigned to fixed folds before
any observed data are used. Bootstrap copies keep their original fold.
"""
import numpy as np
from src.crossfit_cache import features


def pilot_from_means(A,b,box):
    q=np.linalg.solve(np.einsum('...mi,...mj->...ij',A,A),
                      np.einsum('...mi,...m->...i',A,b)[...,None])[...,0]
    return np.clip(q,np.asarray(box)[:,0],np.asarray(box)[:,1])


def crossfit_equations(A,b,Bcache,acache,anchor,box,folds,counts=None,cache_exponents=None,cache_radius=None):
    """Return each bootstrap estimate, sandwich, and complement-fold pilots.

    counts has shape (draws, records); absent counts means one copy per record.
    This routine refits pilots for each draw and excludes all copies of a
    held-out record by keeping bootstrap resampling within original folds.
    """
    n=len(A);folds=np.asarray(folds);anchor=np.asarray(anchor)
    if counts is None:counts=np.ones((1,n))
    counts=np.asarray(counts,dtype=float);total=counts.sum(axis=1)
    sumA=np.einsum('dr,rmi->dmi',counts,A);sumb=counts@b
    aggregateB=np.zeros((len(counts),4,4));aggregatea=np.zeros((len(counts),4))
    raw_second=np.zeros_like(aggregateB);pilots=[];fold_terms=[]
    for k in np.unique(folds):
        take=folds==k;c=counts[:,take];mass=c.sum(axis=1)
        train_mass=total-mass
        if np.any(train_mass<=0):raise ValueError('Empty training fold complement')
        trainA=(sumA-np.einsum('dr,rmi->dmi',c,A[take]))/train_mass[:,None,None]
        trainb=(sumb-c@b[take])/train_mass[:,None]
        pilot=pilot_from_means(trainA,trainb,box);pilots.append(pilot)
        affine=(np.column_stack((np.ones(len(counts)),pilot-anchor)) if cache_exponents is None
                else features(pilot,anchor,cache_radius,cache_exponents))
        # Record equations, evaluated at each complement pilot.
        B=np.einsum('dp,prij->drij',affine,Bcache[:,take]);a=np.einsum('dp,pri->dri',affine,acache[:,take])
        aggregateB+=np.einsum('dr,drij->dij',c,B);aggregatea+=np.einsum('dr,dri->di',c,a)
        fold_terms.append((c,B,a))
    q=np.linalg.solve(aggregateB,aggregatea[...,None])[...,0]
    for c,B,a in fold_terms:
        residual=a-np.einsum('drij,dj->dri',B,q)
        raw_second+=np.einsum('dr,dri,drj->dij',c,residual,residual)
    inv=np.linalg.inv(aggregateB)
    V=np.einsum('dik,dkl,djl->dij',inv,raw_second,inv)*total[:,None,None]/(total-1)[:,None,None]
    V=(V+V.swapaxes(-1,-2))/2
    return q,V,np.stack(pilots,axis=1)


def stratified_counts(rng,folds,draws):
    folds=np.asarray(folds);result=np.zeros((draws,len(folds)),dtype=np.int32)
    for k in np.unique(folds):
        ids=np.flatnonzero(folds==k)
        result[:,ids]=rng.multinomial(len(ids),np.full(len(ids),1/len(ids)),size=draws)
    return result


def crossfit_bootstrap(A,b,Bcache,acache,anchor,box,folds,q,rng,draws=999,cache_exponents=None,cache_radius=None,refit_pilots=True):
    if not refit_pilots:
        _,_,pilots=crossfit_equations(A,b,Bcache,acache,anchor,box,folds,cache_exponents=cache_exponents,cache_radius=cache_radius)
        weights=(np.column_stack((np.ones(len(folds)),pilots[0,folds]-anchor)) if cache_exponents is None
                 else features(pilots[0,folds],anchor,cache_radius,cache_exponents))
        Bcache=np.einsum('rp,prij->rij',weights,Bcache)[None]
        acache=np.einsum('rp,pri->ri',weights,acache)[None]
        cache_exponents=np.zeros((1,4),dtype=int);cache_radius=1.
    pivots=[];joint=[]
    for start in range(0,draws,32):
        counts=stratified_counts(rng,folds,min(32,draws-start))
        qb,Vb,_=crossfit_equations(A,b,Bcache,acache,anchor,box,folds,counts,cache_exponents,cache_radius)
        if not np.isfinite(qb).all() or np.linalg.eigvalsh(Vb).min()<=0:
            raise RuntimeError('Invalid cross-fit bootstrap fit; no draws discarded')
        error=qb-q
        pivots.append(error/np.sqrt(np.diagonal(Vb,axis1=-2,axis2=-1)))
        joint.append(np.einsum('di,di->d',error,np.linalg.solve(Vb,error[...,None])[...,0]))
    return np.quantile(np.concatenate(pivots),[.025,.975],axis=0),float(np.quantile(np.concatenate(joint),.95))
