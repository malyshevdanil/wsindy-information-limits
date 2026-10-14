"""Full-library guarded weak GLS using smoothed observations for instruments."""
import numpy as np
from scipy.linalg import eigh
from .weak_features import library,derivative_polynomial
from .polynomial_gls import moments,covariance
from .boundary_weak import quadrature_weights
from .guarded_gls import guarded_gls


def dyadic_cubic_tests(t,level):
    """Independent C2 cardinal cubic splines, with zero endpoint traces."""
    n=2**level;step=float(t[-1])/n;centers=np.arange(2,n-1)*step
    z=(t[None,:]-centers[:,None])/step;r=np.abs(z);sign=np.sign(z)
    phi=np.where(r<1,(4-6*r*r+3*r**3)/6,np.where(r<2,(2-r)**3/6,0.))
    der=np.where(r<1,sign*(-12*r+9*r*r)/6,np.where(r<2,-sign*(2-r)**2/2,0.))/step
    return phi,der


def interpolation_kernel(h,bandwidth):
    """Convolution of compact polynomial kernel with linear interpolation hats.

    Split at hat/kernel knots; five Gauss nodes integrate each degree-nine
    polynomial piece exactly in real arithmetic. Derivative weights refer
    to the derivative of the same continuous convolution.
    """
    nodes,weights=np.polynomial.legendre.leggauss(5);m=int(np.ceil(bandwidth/h))+1
    kernel=np.zeros(2*m+1);derivative=np.zeros_like(kernel)
    for index,lag in enumerate(range(-m,m+1)):
        cuts=sorted(set([-h,0.,h]+[q for q in [lag*h-bandwidth,lag*h+bandwidth] if -h<q<h]))
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            s=(lo+hi)/2+(hi-lo)/2*nodes;z=(lag*h-s)/bandwidth
            base=np.maximum(1-z*z,0);hat=1-np.abs(s)/h;factor=315/256
            kernel[index]+=np.dot(weights,hat*factor*base**4/bandwidth)*(hi-lo)/2
            derivative[index]+=np.dot(weights,hat*factor*(-8*z)*base**3/bandwidth**2)*(hi-lo)/2
    mass=kernel.sum();kernel/=mass;derivative/=mass
    kernel=(kernel+kernel[::-1])/2;derivative=(derivative-derivative[::-1])/2
    return kernel,derivative,m


def proxy_and_tests(t,y,epsilon,powers,universal_phi,universal_der):
    h=float(t[1]-t[0]);bandwidth=.25*epsilon**.1
    kernel,derivative,m=interpolation_kernel(h,bandwidth)
    padded=np.pad(y,((m,m),(0,0)),mode='reflect')
    proxy=np.column_stack([np.convolve(padded[:,j],kernel,mode='valid') for j in range(y.shape[1])])
    proxy_der=np.column_stack([np.convolve(padded[:,j],derivative,mode='valid') for j in range(y.shape[1])])
    center=float(t[-1])/2;radius=center-.5;z=(t-center)/radius
    base=np.maximum(1-z*z,0);w=base**4;wp=np.where(np.abs(z)<1,-8*z*base**3/radius,0.)
    features=library(proxy,powers,np.zeros(y.shape[1]));gradient=np.zeros_like(features)
    for j in range(y.shape[1]):
        alpha=[0]*y.shape[1];alpha[j]=1
        gradient+=derivative_polynomial(proxy,powers,np.eye(len(powers)),tuple(alpha))*proxy_der[:,j,None]
    pilot_phi=(w[:,None]*features).T;pilot_der=(wp[:,None]*features+w[:,None]*gradient).T
    raw_phi=np.concatenate([pilot_phi,universal_phi]);raw_der=np.concatenate([pilot_der,universal_der])
    weights=quadrature_weights(t);gram=(raw_phi*weights)@raw_phi.T+(raw_der*weights)@raw_der.T
    values,V=eigh((gram+gram.T)/2);cutoff=1e-10*epsilon**.1;keep=values>cutoff*values[-1]
    transform=(V[:,keep]/np.sqrt(values[keep])).T
    return proxy,pilot_phi,pilot_der,transform@raw_phi,transform@raw_der


def fit(t,y,epsilon,sigma_shape,tau,powers,universal_phi,universal_der,bounds):
    """Requires no ODE integration or supplied generating support.

    Finite-noise unbiasedness is not asserted for data-dependent instruments.
    All coefficients in the supplied polynomial library are fitted.
    """
    bounds=np.asarray(bounds,dtype=float);fallback=bounds.mean(axis=1)
    if not 0<epsilon<=1 or t[-1]<=1 or bounds.shape!=(len(powers)*y.shape[1],2) or not np.isfinite(bounds).all() or np.any(bounds[:,0]>=bounds[:,1]):
        raise ValueError('Valid amplitude, horizon and compact coefficient box required')
    try:
        proxy,pilot_phi,pilot_der,phi,der=proxy_and_tests(t,y,epsilon,powers,universal_phi,universal_der)
        sigma=epsilon*np.asarray(sigma_shape)
        A0,b0=moments(t,y,pilot_phi,pilot_der,powers,sigma)
        scales=np.linalg.norm(A0,axis=0)
        if np.any(scales==0):return fallback,{'fallback':True,'reason':'pilot_zero_column'}
        design=A0/scales;normal=design.T@design;values=np.linalg.eigvalsh(normal)
        if values[0]<=1e-8*epsilon**2*values[-1]:return fallback,{'fallback':True,'reason':'pilot_rank_guard'}
        raw=np.linalg.solve(normal,design.T@b0)/scales;pilot=np.clip(raw,bounds[:,0],bounds[:,1])
        A,b=moments(t,y,phi,der,powers,sigma)
        C=covariance(t,proxy,pilot.reshape(y.shape[1],len(powers)).T,phi,der,powers,sigma,tau)
        estimate,diag=guarded_gls(A,b,C,epsilon,bounds)
        return estimate,{**diag,'pilot_clipped':bool(np.any(pilot!=raw)),
            'retained_tests':len(phi),'raw_tests':len(pilot_phi)+len(universal_phi),
            'bandwidth':float(.25*epsilon**.1)}
    except (np.linalg.LinAlgError,ValueError,FloatingPointError,OverflowError):
        return fallback,{'fallback':True,'reason':'instrument_or_covariance_failure'}
