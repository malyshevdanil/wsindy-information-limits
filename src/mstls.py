"""Independent implementation of relative-magnitude MSTLS regression.

Definitions follow Messenger--Bortz's regression rule and the published author
Python translation in WmSINDy, with gamma=0, unit rescaling, alpha_loss=.8.
This is the regression component on supplied weak matrices, not the adaptive
test-function/data-denoising pipeline of either complete author algorithm.
"""
import numpy as np


def relative_mstls(G, b, thresholds=None, alpha_loss=.8):
    G=np.asarray(G,float)
    b=np.asarray(b,float)
    if b.ndim==1:
        b=b[:,None]
    thresholds=np.logspace(-4,0,50) if thresholds is None else np.asarray(thresholds)
    if np.any(thresholds<=0) or not 0<alpha_loss<1:
        raise ValueError('Need positive thresholds and alpha_loss in (0,1)')
    column_norms=np.linalg.norm(G,axis=0)
    valid=column_norms>0
    answer=np.zeros((G.shape[1],b.shape[1]))
    choices=[]
    for equation in range(b.shape[1]):
        rhs=b[:,equation]
        base=np.linalg.lstsq(G,rhs,rcond=-1)[0]
        reference_norm=np.linalg.norm(G@base)
        if reference_norm==0:
            choices.append({'threshold':float(thresholds[0]),'loss':0.})
            continue
        scales=np.divide(np.linalg.norm(rhs),column_norms,
                         out=np.full_like(column_norms,np.inf),where=valid)
        solved={tuple(valid):np.linalg.lstsq(G[:,valid],rhs,rcond=-1)[0]}
        candidates=[]
        for threshold in thresholds:
            coefficients=base.copy()
            active=valid.copy()
            lower=threshold*np.maximum(1.,scales)
            upper=np.minimum(1.,scales)/threshold
            for _ in range(G.shape[1]):
                keep=valid & (np.abs(coefficients)>=lower) & (np.abs(coefficients)<=upper)
                if np.array_equal(keep,active):
                    break
                active=keep
                coefficients[:]=0.
                if not active.any():
                    break
                key=tuple(active)
                if key not in solved:
                    solved[key]=np.linalg.lstsq(G[:,active],rhs,rcond=-1)[0]
                coefficients[active]=solved[key]
            loss=(2*alpha_loss*np.linalg.norm(G@(coefficients-base))/reference_norm
                  +2*(1-alpha_loss)*np.count_nonzero(coefficients)/len(coefficients))
            candidates.append((float(loss),float(threshold),coefficients))
        loss,threshold,coefficients=min(candidates,key=lambda item:item[0])
        answer[:,equation]=coefficients
        choices.append({'threshold':threshold,'loss':loss})
    return answer,choices
