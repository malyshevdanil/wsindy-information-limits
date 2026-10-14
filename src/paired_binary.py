"""Conservative finite-sample intervals for a paired binary risk difference."""
import numpy as np
from scipy.stats import beta


def clopper_pearson(k,n,alpha):
    if n<1 or not 0<=k<=n or not 0<alpha<1:raise ValueError('Invalid binomial inputs')
    lower=0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1))
    upper=1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))
    return lower,upper


def paired_difference_interval(first,second,alpha=.05):
    """Each discordant probability gets a 1-alpha/2 CP interval.

    Union bound, not independence of discordant counts, gives joint coverage
    >=1-alpha and hence coverage for p10-p01. This is pointwise per cell.
    """
    first=np.asarray(first);second=np.asarray(second)
    if first.ndim!=1 or first.shape!=second.shape or len(first)==0:raise ValueError('Unequal or empty paired vectors')
    if not np.isin(first,[0,1]).all() or not np.isin(second,[0,1]).all():raise ValueError('Outcomes must be binary')
    n=len(first);n10=int(np.sum((first==1)&(second==0)));n01=int(np.sum((first==0)&(second==1)))
    n11=int(np.sum((first==1)&(second==1)));n00=n-n10-n01-n11
    l10,u10=clopper_pearson(n10,n,alpha/2);l01,u01=clopper_pearson(n01,n,alpha/2)
    return {'n':n,'n00':n00,'n01':n01,'n10':n10,'n11':n11,
        'difference':(n10-n01)/n,'ci95':[max(-1.,l10-u01),min(1.,u10-l01)],
        'method':'Bonferroni-projected Clopper-Pearson discordant-cell intervals',
        'coverage_at_least':1-alpha,'simultaneous_across_grid':False}
