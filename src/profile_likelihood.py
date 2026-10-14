"""Shared-coefficient OU likelihood with record-specific initial states.

Block Gauss--Newton avoids forming a (4+2R)-square dense matrix. Fits are
bounded local likelihood optimizations; convergence is recorded explicitly.
"""
import time
import numpy as np
from scipy.integrate import solve_ivp


def ensemble_sensitivities(t,q,u,rtol=1e-9,atol=1e-11):
    R=len(u);a,b,c,d=q;initial=np.zeros((R,14));initial[:,:2]=u
    initial[:,2:].reshape(R,2,6)[:,:,4:]=np.eye(2)
    def rhs(_,state):
        state=state.reshape(R,14);x,v=state[:,:2].T;U=state[:,2:].reshape(R,2,6)
        derivative=np.empty_like(U);derivative[:,0]=a*U[:,1]
        derivative[:,1]=(b+2*d*x*v)[:,None]*U[:,0]+(c+d*x*x)[:,None]*U[:,1]
        derivative[:,0,0]+=v;derivative[:,1,1]+=x;derivative[:,1,2]+=v;derivative[:,1,3]+=x*x*v
        return np.column_stack((a*v,b*x+c*v+d*x*x*v,derivative.reshape(R,12))).ravel()
    sol=solve_ivp(rhs,(0.,float(t[-1])),initial.ravel(),t_eval=t,method='DOP853',rtol=rtol,atol=atol)
    if not sol.success:raise RuntimeError('Ensemble sensitivity integration failed')
    state=sol.y.T.reshape(len(t),R,14).transpose(1,0,2)
    return state[:,:,:2],state[:,:,2:].reshape(R,len(t),2,6)


def whiten_records(values,h,tau,sigma):
    rho=np.exp(-h/tau);scale=np.asarray(sigma).reshape((1,1,2)+(1,)*(values.ndim-3))
    result=np.empty_like(values);result[:,0:1]=values[:,0:1]/scale
    result[:,1:]=(values[:,1:]-rho*values[:,:-1])/(scale*np.sqrt(-np.expm1(-2*h/tau)))
    return result


def ensemble_hessian(t,q,u,rtol=1e-10,atol=1e-12):
    """Full second sensitivities in (four coefficients, two initial states)."""
    R=len(u);a,b,c,d=q;initial=np.zeros((R,86));initial[:,:2]=u
    initial[:,2:14].reshape(R,2,6)[:,:,4:]=np.eye(2)
    basis=np.eye(6)
    def rhs(_,state):
        state=state.reshape(R,86);x,v=state[:,:2].T;U=state[:,2:14].reshape(R,2,6);H=state[:,14:].reshape(R,2,6,6)
        jx=b+2*d*x*v;jv=c+d*x*x
        dU=np.empty_like(U);dU[:,0]=a*U[:,1];dU[:,1]=jx[:,None]*U[:,0]+jv[:,None]*U[:,1]
        dU[:,0,0]+=v;dU[:,1,1]+=x;dU[:,1,2]+=v;dU[:,1,3]+=x*x*v
        dH=np.empty_like(H)
        dH[:,0]=a*H[:,1]+np.einsum('i,rj->rij',basis[0],U[:,1])+np.einsum('ri,j->rij',U[:,1],basis[0])
        source=2*d*v[:,None,None]*np.einsum('ri,rj->rij',U[:,0],U[:,0])
        source+=2*d*x[:,None,None]*(np.einsum('ri,rj->rij',U[:,0],U[:,1])+np.einsum('ri,rj->rij',U[:,1],U[:,0]))
        mixed=np.einsum('ri,j->rij',U[:,0],basis[1])+np.einsum('ri,j->rij',U[:,1],basis[2])
        mixed+=np.einsum('ri,j->rij',2*(x*v)[:,None]*U[:,0]+(x*x)[:,None]*U[:,1],basis[3])
        dH[:,1]=jx[:,None,None]*H[:,0]+jv[:,None,None]*H[:,1]+source+mixed+mixed.swapaxes(1,2)
        return np.column_stack((a*v,b*x+c*v+d*x*x*v,dU.reshape(R,12),dH.reshape(R,72))).ravel()
    sol=solve_ivp(rhs,(0.,float(t[-1])),initial.ravel(),t_eval=t,method='DOP853',rtol=rtol,atol=atol)
    if not sol.success:raise RuntimeError('Full sensitivity Hessian integration failed')
    state=sol.y.T.reshape(len(t),R,86).transpose(1,0,2)
    return state[:,:,:2],state[:,:,2:14].reshape(R,len(t),2,6),state[:,:,14:].reshape(R,len(t),2,6,6)


def initial_curvature(t,q,u,rtol=1e-9,atol=1e-11):
    """Second derivatives of the trajectory with respect to its initial state."""
    R=len(u);a,b,c,d=q;initial=np.zeros((R,14));initial[:,:2]=u
    initial[:,2:6].reshape(R,2,2)[:]=np.eye(2)
    def rhs(_,state):
        state=state.reshape(R,14);x,v=state[:,:2].T;Z=state[:,2:6].reshape(R,2,2);H=state[:,6:].reshape(R,2,2,2)
        jx=b+2*d*x*v;jv=c+d*x*x;dZ=np.empty_like(Z);dZ[:,0]=a*Z[:,1];dZ[:,1]=jx[:,None]*Z[:,0]+jv[:,None]*Z[:,1]
        dH=np.empty_like(H);dH[:,0]=a*H[:,1]
        source=2*d*v[:,None,None]*np.einsum('ri,rj->rij',Z[:,0],Z[:,0])
        source+=2*d*x[:,None,None]*(np.einsum('ri,rj->rij',Z[:,0],Z[:,1])+np.einsum('ri,rj->rij',Z[:,1],Z[:,0]))
        dH[:,1]=jx[:,None,None]*H[:,0]+jv[:,None,None]*H[:,1]+source
        return np.column_stack((a*v,b*x+c*v+d*x*x*v,dZ.reshape(R,4),dH.reshape(R,8))).ravel()
    sol=solve_ivp(rhs,(0.,float(t[-1])),initial.ravel(),t_eval=t,method='DOP853',rtol=rtol,atol=atol)
    if not sol.success:raise RuntimeError('Initial-curvature integration failed')
    return sol.y.T.reshape(len(t),R,14).transpose(1,0,2)[:,:,6:].reshape(R,len(t),2,2,2)


def normal_blocks(W,residual):
    W=W.reshape(len(W),-1,W.shape[-1]);r=residual.reshape(len(W),-1)
    return np.einsum('rmi,rmj->rij',W,W,optimize=True),np.einsum('rmi,rm->ri',W,r,optimize=True)


def block_step(info,score,unknown_initial=True,damping=0.):
    H=info.copy();ids=np.arange(H.shape[-1]);H[:,ids,ids]+=damping*np.maximum(np.diagonal(info,axis1=1,axis2=2),1e-12)
    if not unknown_initial:
        S=H[:,:4,:4].sum(axis=0);g=score[:,:4].sum(axis=0);return -np.linalg.solve(S,g),None,S,g
    C=H[:,:4,4:];N=H[:,4:,4:];NiC=np.linalg.solve(N,C.swapaxes(1,2));Nig=np.linalg.solve(N,score[:,4:,None])[...,0]
    S=(H[:,:4,:4]-C@NiC).sum(axis=0);g=(score[:,:4]-np.einsum('rij,rj->ri',C,Nig)).sum(axis=0)
    dq=-np.linalg.solve(S,g);du=-Nig-np.einsum('rij,j->ri',NiC,dq)
    return dq,du,(S+S.T)/2,g


def fit_profile(t,y,q0,u0,tau,sigma,config,known_initial=False):
    start=time.perf_counter();q=np.clip(q0,config['coefficient_lower'],config['coefficient_upper']);u=np.asarray(u0).copy()
    if not known_initial:u=np.clip(u,config['initial_lower'],config['initial_upper'])
    damping=config['initial_damping'];evaluations=0;success=False;reason='maximum_iterations';history=[]
    def evaluate(q,u):
        nonlocal evaluations
        model,U=ensemble_sensitivities(t,q,u,config['ode_rtol'],config['ode_atol']);evaluations+=1
        r=whiten_records(model-y,t[1]-t[0],tau,sigma);W=whiten_records(U,t[1]-t[0],tau,sigma)
        I,g=normal_blocks(W,r);return float(np.sum(r*r)/2),I,g
    cost,I,g=evaluate(q,u)
    for iteration in range(config['maximum_iterations']):
        dq,du,S,profile_score=block_step(I,g,not known_initial,damping)
        if max(np.max(np.abs(dq)),0. if du is None else np.max(np.abs(du)))<config['step_tolerance']:
            success=True;reason='step_tolerance';break
        accepted=False
        for j in range(config['maximum_line_search']):
            alpha=2.**(-j);qn=np.clip(q+alpha*dq,config['coefficient_lower'],config['coefficient_upper'])
            un=u if known_initial else np.clip(u+alpha*du,config['initial_lower'],config['initial_upper'])
            cn,In,gn=evaluate(qn,un)
            if cn<cost:
                change=cost-cn;previous=cost;q,u,cost,I,g=qn,un,cn,In,gn;accepted=True;damping=max(damping*.2,1e-12)
                history.append(cost)
                if change<=config['objective_tolerance']*(1+previous):
                    # An objective stop is reported, with score diagnostics.
                    success=True;reason='objective_tolerance'
                break
        if success:break
        if not accepted:
            damping*=10.
            if damping>config['maximum_damping']:reason='line_search_failed';break
    dq,du,S,profile_score=block_step(I,g,not known_initial,0.)
    V=np.linalg.inv(S);active_q=((q<=np.asarray(config['coefficient_lower'])+1e-8)|(q>=np.asarray(config['coefficient_upper'])-1e-8))
    active_u=np.zeros_like(u,dtype=bool) if known_initial else ((u<=np.asarray(config['initial_lower'])+1e-8)|(u>=np.asarray(config['initial_upper'])-1e-8))
    return q,u,V,dict(success=success,reason=reason,iterations=iteration+1,evaluations=evaluations,objective=cost,
        remaining_newton_step=float(max(np.max(np.abs(dq)),0. if du is None else np.max(np.abs(du)))),
        standardized_profile_score=float(np.max(np.abs(profile_score)/np.sqrt(np.diag(S)))),
        bound_active_coefficients=int(active_q.sum()),bound_active_initial_coordinates=int(active_u.sum()),elapsed_seconds=time.perf_counter()-start,objective_history=history)
