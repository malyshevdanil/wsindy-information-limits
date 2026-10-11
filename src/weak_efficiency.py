"""Optimal boundary-free influences and covariance-weighted weak components."""
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.linalg import cho_factor,cho_solve,eigh
from .weak_features import test_functions,quadrature_weights
from .covariance import exponential_apply


def scalar_information(a,x0,T,tau,sigma):
    c=a+1/tau;f=tau/(2*sigma*sigma)
    F=[quad(lambda t,q=q:t**q*np.exp(2*a*t),0,T,epsabs=1e-11,epsrel=1e-11)[0] for q in range(3)]
    Iaa=f*x0*x0*(F[0]+2*c*F[1]+c*c*F[2])
    Iau=f*x0*c*(F[0]+c*F[1]);Iuu=1/sigma**2+f*c*c*F[0]
    return Iaa,Iau,Iuu,Iaa-Iau*Iau/Iuu


def scalar_optimal_test(t,a,x0,T,tau,sigma,delta):
    """Smooth the optimal adjoint test; no endpoint observations enter the fit."""
    Iaa,Iau,Iuu,Ieff=scalar_information(a,x0,T,tau,sigma)
    k=Iau/Iuu
    A=(1-tau*tau*a*a)*x0/(2*sigma*sigma*tau)
    B=(-(1-tau*tau*a*a)*k-2*tau*tau*a*x0)/(2*sigma*sigma*tau)
    atom0=(-k-tau*(x0-a*k))/(2*sigma*sigma)
    if abs(a*T)<1e-4:
        F0=sum((2*a)**q*t**(q+1)/(__import__('math').factorial(q)*(q+1)) for q in range(10))
        F1=sum((2*a)**q*t**(q+2)/(__import__('math').factorial(q)*(q+2)) for q in range(10))
    else:
        F0=np.expm1(2*a*t)/(2*a)
        F1=(np.exp(2*a*t)*(2*a*t-1)+1)/(4*a*a)
    phi=-np.exp(-a*t)*(atom0+A*F1+B*F0)/Ieff
    density=np.exp(a*t)*(A*t+B)/Ieff
    derivative=-a*phi-density
    def ramp(distance):
        z=np.clip((distance-delta/4)/(3*delta/4),0,1)
        eta=z**3*(10-15*z+6*z*z)
        prime=30*z*z*(1-z)**2/(3*delta/4)
        return eta,prime
    left,ld=ramp(t);right,rd=ramp(T-t)
    cutoff=left*right;cutoff_derivative=ld*right-left*rd
    return cutoff*phi,cutoff*derivative+cutoff_derivative*phi


def nested_bumps(t,count,cutoff=1e-10):
    """A fixed, data-independent H1 basis for nested multiscale bump spans."""
    schedule=[(30,3.),(30,1.5),(60,.75),(120,.375),(240,.1875)]
    phis=[];ders=[];used=0;T=t[-1]
    for n,L in schedule:
        if used>=count:break
        phi,der=test_functions(t,np.linspace(L/2+.0025,T-L/2-.0025,n),L)
        phis.extend(phi);ders.extend(der);used+=n
    assert used==count
    phi=np.array(phis);der=np.array(ders);w=quadrature_weights(t)
    gram=(phi*w)@phi.T+(der*w)@der.T
    vals,V=eigh((gram+gram.T)/2);keep=vals>cutoff*vals.max()
    transform=(V[:,keep]/np.sqrt(vals[keep])).T
    return transform@phi,transform@der,{'raw_windows':count,'retained_windows':int(keep.sum()),
        'relative_cutoff':cutoff,'minimum_retained_relative_eigenvalue':float(vals[keep].min()/vals.max())}


def vdp_pilot_trajectory(t,q,initial,box):
    constrained=np.clip(q,np.array(box)[:,0],np.array(box)[:,1])
    a,b,c,d=constrained
    def rhs(_,state):
        x,v=state;return [a*v,b*x+c*v+d*x*x*v]
    sol=solve_ivp(rhs,(0.,float(t[-1])),initial,t_eval=t,method='DOP853',rtol=1e-9,atol=1e-11)
    if not sol.success:raise RuntimeError('Pilot trajectory integration failed')
    return sol.y.T,constrained


def vdp_moments(t,y,phi,der,sigma):
    w=quadrature_weights(t);x,v=y.T
    K=len(phi);A=np.zeros((2*K,4));b=-(der*w)@y
    A[0::2,0]=(phi*w)@v
    A[1::2,1]=(phi*w)@x;A[1::2,2]=(phi*w)@v
    A[1::2,3]=(phi*w)@((x*x-sigma[0]**2)*v)
    return A,b.ravel()


def covariance_blocks(t,model,phi,der,sigma,tau):
    """Precompute the exact Prop3 covariance polynomial in four coefficients."""
    w=quadrature_weights(t);h=t[1]-t[0];x,v=model.T
    arrays={'A':(phi*w).T,'B':(der*w).T,'V':(phi*(x*v*w)).T,'X':(phi*(x*x*w)).T}
    def apply(U,order=1):
        return U if tau is None else exponential_apply(U,np.exp(-order*h/tau))
    pairs=['AA','AB','BB','AV','BV','VV','AX','BX','XX']
    applied={k:apply(U) for k,U in arrays.items()}
    blocks={pair:arrays[pair[0]].T@applied[pair[1]] for pair in pairs}
    # The only nonzero corrected higher-order coefficients are d*v, 2*d*x, d.
    Fv=(phi*(v*w)).T;Fx=(phi*(x*w)).T;Fa=arrays['A']
    higher=(2*sigma[0]**4*(Fv.T@apply(Fv,2))+
            4*sigma[0]**2*sigma[1]**2*(Fx.T@apply(Fx,2))+
            2*sigma[0]**4*sigma[1]**2*(Fa.T@apply(Fa,3)))
    return blocks,(higher+higher.T)/2


def assemble_covariance(q,blocks,higher,sigma):
    a,b,c,d=q;B=blocks;sx,sv=sigma**2
    C11=sx*B['BB']+sv*a*a*B['AA']
    C12=sx*(b*B['AB'].T+2*d*B['BV'])+sv*a*(B['AB']+c*B['AA']+d*B['AX'])
    C22=sx*(b*b*B['AA']+2*b*d*(B['AV']+B['AV'].T)+4*d*d*B['VV'])
    C22+=sv*(B['BB']+c*(B['AB']+B['AB'].T)+d*(B['BX']+B['BX'].T)+
             c*c*B['AA']+c*d*(B['AX']+B['AX'].T)+d*d*B['XX'])
    C22+=d*d*higher
    K=len(C11);C=np.empty((2*K,2*K))
    C[0::2,0::2]=C11;C[0::2,1::2]=C12;C[1::2,0::2]=C12.T;C[1::2,1::2]=C22
    return (C+C.T)/2


def gls_map(A,C,relaxation=1e-10):
    level=relaxation*max(float(np.trace(C)/len(C)),np.finfo(float).tiny)
    chol=cho_factor(C+level*np.eye(len(C)),lower=True,check_finite=False)
    WA=cho_solve(chol,A,check_finite=False)
    information=A.T@WA
    return np.linalg.solve(information,WA.T),level


def fit_covariance_iteration(A,b,q0,blocks,higher,sigma,iterations,tolerance,relaxation):
    q=q0.copy();history=[]
    for _ in range(iterations):
        C=assemble_covariance(q,blocks,higher,sigma)
        P,level=gls_map(A,C,relaxation);new=P@b
        change=float(np.linalg.norm(new-q)/max(np.linalg.norm(q),1e-30))
        history.append(change);q=new
        if change<tolerance:break
    return q,{'iterations':len(history),'converged':history[-1]<tolerance,'last_change':history[-1],
              'relaxation':level}
