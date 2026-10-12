"""Exact discrete stationary OU score and ODE sensitivities."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


def vdp_sensitivities(t,q,initial):
    a,b,c,d=q
    def rhs(_,state):
        x,v=state[:2];U=state[2:].reshape(2,6)
        J=np.array([[0,a],[b+2*d*x*v,c+d*x*x]])
        forcing=np.array([[v,0,0,0],[0,x,v,x*x*v]])
        derivative=J@U;derivative[:,:4]+=forcing
        return np.r_[a*v,b*x+c*v+d*x*x*v,derivative.ravel()]
    U=np.zeros((2,6));U[:,4:]=np.eye(2)
    sol=solve_ivp(rhs,(0,float(t[-1])),np.r_[initial,U.ravel()],t_eval=t,method='DOP853',rtol=1e-10,atol=1e-12)
    if not sol.success:raise RuntimeError('Sensitivity integration failed')
    return sol.y[:2].T,sol.y[2:].T.reshape(len(t),2,6)


def ou_whiten(values,h,tau,sigma):
    rho=np.exp(-h/tau);scale=np.asarray(sigma)
    if values.ndim==3:scale=scale[:,None]
    result=np.empty_like(values);result[0]=values[0]/scale
    result[1:]=(values[1:]-rho*values[:-1])/(scale*np.sqrt(-np.expm1(-2*h/tau)))
    return result


def one_step_vdp(t,y,q,initial,tau,sigma,unknown_initial=False):
    model,U=vdp_sensitivities(t,q,initial)
    if not unknown_initial:U=U[:,:,:4]
    W=ou_whiten(U,t[1]-t[0],tau,sigma).reshape(-1,U.shape[2])
    residual=ou_whiten(y-model,t[1]-t[0],tau,sigma).ravel()
    scale=np.linalg.norm(W,axis=0)
    delta=np.linalg.lstsq(W/scale,residual,rcond=None)[0]/scale
    result=q+delta[:4]
    return result,dict(initial_estimate=(np.asarray(initial)+delta[4:]).tolist() if unknown_initial else list(initial),
        scaled_information_condition=float(np.linalg.cond((W/scale).T@(W/scale))))


def profile_ou_shape(residual,h,bounds):
    """Stationary OU profile likelihood of fitted-mean residuals, not joint MLE."""
    N=len(residual)
    def evaluate(logtau,return_variance=False):
        tau=np.exp(logtau);rho=np.exp(-h/tau);den=-np.expm1(-2*h/tau)
        variance=(residual[0]**2+np.sum((residual[1:]-rho*residual[:-1])**2,axis=0)/den)/N
        if return_variance:return variance
        return float(N*np.log(variance).sum()+residual.shape[1]*(N-1)*np.log(den))
    fitted=minimize_scalar(evaluate,bounds=np.log(bounds),method='bounded',options={'xatol':1e-7})
    if not fitted.success:raise RuntimeError('OU residual-profile fit failed')
    return float(np.exp(fitted.x)),np.sqrt(evaluate(fitted.x,True))
