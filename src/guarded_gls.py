"""Bounded vanishing-ridge variant; archived experiment solvers are unchanged."""
import numpy as np
from scipy.linalg import cho_factor,cho_solve


def guarded_gls(A,b,covariance,epsilon,bounds):
    """Full-rank polynomial weak GLS with rho_epsilon=epsilon**2.

    bounds is a fixed compact coefficient box. Covariance is physical, so
    divide by epsilon**2 before ridge and numerical guards. This variant
    is analyzed in S27; it has not been used in archived Monte Carlo runs.
    Returns a bounded estimate and execution diagnostics, not an exact SE.
    """
    A=np.asarray(A,dtype=float);b=np.asarray(b,dtype=float);covariance=np.asarray(covariance,dtype=float)
    bounds=np.asarray(bounds,dtype=float)
    if A.ndim!=2 or bounds.shape!=(A.shape[1],2) or b.shape!=(A.shape[0],) or covariance.shape!=(len(b),len(b)):
        raise ValueError('Incompatible GLS inputs')
    if not np.isfinite(epsilon) or not 0<epsilon<=1 or not np.isfinite(bounds).all() or np.any(bounds[:,0]>=bounds[:,1]):
        raise ValueError('Amplitude in (0,1] and finite ordered coefficient bounds required')
    fallback=bounds.mean(axis=1)
    def failed(reason):return fallback.copy(),{'fallback':True,'reason':reason,'relative_ridge':float(epsilon**2)}
    if not all(np.isfinite(v).all() for v in [A,b,covariance]):return failed('nonfinite_input')
    scales=np.linalg.norm(A,axis=0)
    if np.any(scales==0):return failed('zero_design_column')
    shape=(covariance+covariance.T)/(2*epsilon**2);average=float(np.trace(shape)/len(shape))
    if average<=0:return failed('invalid_covariance_scale')
    threshold=epsilon**2*1e-8
    if np.linalg.eigvalsh(shape)[0]<=threshold*average:return failed('covariance_guard')
    relative_ridge=epsilon**2;ridge=relative_ridge*average
    try:
        chol=cho_factor(shape+ridge*np.eye(len(shape)),lower=True,check_finite=False)
        design=A/scales;weighted=cho_solve(chol,design,check_finite=False);normal=design.T@weighted
        values=np.linalg.eigvalsh(normal)
        if values[0]<=threshold*values[-1]:return failed('design_guard')
        estimate=np.linalg.solve(normal,weighted.T@b)/scales
    except np.linalg.LinAlgError:return failed('linear_algebra_guard')
    if not np.isfinite(estimate).all():return failed('nonfinite_estimate')
    bounded=np.clip(estimate,bounds[:,0],bounds[:,1])
    return bounded,{'fallback':False,'clipped':bool(np.any(bounded!=estimate)),
        'relative_ridge':float(relative_ridge),'normalized_covariance_ridge':float(ridge),
        'physical_covariance_ridge':float(epsilon**2*ridge)}


def fit_guarded_polynomial(t,y,phi,der,pilot_phi,pilot_der,powers,epsilon,
        sigma_shape,tau,core,pilot_bounds,final_bounds,integrate_pilot):
    """Full deterministic-window pipeline of Proposition S27.1.

    integrate_pilot(t, core_coefficients) is a supplied ODE solver with its
    known initial state fixed. Research uses of this pipeline run on Kaggle.
    Actual archived runners still use their original solver, not this one.
    """
    from .polynomial_gls import moments,covariance
    sigma=epsilon*np.asarray(sigma_shape);dimension=np.asarray(y).shape[1]
    if not np.isfinite(tau) or tau<=0 or np.asarray(sigma_shape).shape!=(dimension,) or np.any(np.asarray(sigma_shape)<=0):
        raise ValueError('Positive common OU range and coordinate noise scales required')
    if np.any(np.asarray(phi)[:,[0,-1]]!=0.) or np.any(np.asarray(pilot_phi)[:,[0,-1]]!=0.):
        raise ValueError('This variant requires boundary-free test traces')
    pilot_bounds=np.asarray(pilot_bounds,dtype=float);final_bounds=np.asarray(final_bounds,dtype=float)
    if not np.isfinite(epsilon) or not 0<epsilon<=1 or final_bounds.shape!=(len(powers)*dimension,2) or not np.isfinite(final_bounds).all() or np.any(final_bounds[:,0]>=final_bounds[:,1]):
        raise ValueError('Valid amplitude and final compact box required')
    if pilot_bounds.shape!=(len(core),2) or not np.isfinite(pilot_bounds).all() or np.any(pilot_bounds[:,0]>=pilot_bounds[:,1]):
        raise ValueError('Finite ordered pilot box required')
    def fallback(reason):return final_bounds.mean(axis=1),{'fallback':True,'reason':reason}
    A0,b0=moments(t,y,pilot_phi,pilot_der,powers,sigma);B=A0[:,core]
    try:
        scales=np.linalg.norm(B,axis=0)
        if np.any(scales==0):return fallback('pilot_design_guard')
        Bs=B/scales;normal=Bs.T@Bs;values=np.linalg.eigvalsh(normal)
        if values[0]<=1e-8*epsilon**2*values[-1]:return fallback('pilot_design_guard')
        pilot=np.linalg.solve(normal,Bs.T@b0)/scales
        pilot=np.clip(pilot,pilot_bounds[:,0],pilot_bounds[:,1]);model=integrate_pilot(t,pilot)
        coef=np.zeros(len(powers)*dimension);coef[core]=pilot
        A,b=moments(t,y,phi,der,powers,sigma)
        Am,bm=moments(t,model,phi,der,powers,np.zeros(dimension));b=b-(bm-Am@coef)
        C=covariance(t,model,coef.reshape(dimension,len(powers)).T,phi,der,powers,sigma,tau)
        estimate,diag=guarded_gls(A,b,C,epsilon,final_bounds)
        return estimate,{**diag,'pilot':pilot.tolist()}
    except (np.linalg.LinAlgError,RuntimeError,FloatingPointError,ValueError,OverflowError):return fallback('pilot_or_covariance_failure')
