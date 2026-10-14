import numpy as np
from src.profile_likelihood import block_step,ensemble_sensitivities,whiten_records
from src.path_scoring import vdp_sensitivities,ou_whiten
from src.profile_likelihood import ensemble_hessian,initial_curvature


def test_full_hessian_matches_mixed_sensitivity_differences_and_initial_block():
    t=np.linspace(0,.3,11);q=np.array([1.,-1.,1.,-1.]);u=np.array([[2.,0.],[1.8,.2]])
    x,U,H=ensemble_hessian(t,q,u)
    np.testing.assert_allclose(H,H.swapaxes(-1,-2),atol=1e-12)
    np.testing.assert_allclose(H[...,4:,4:],initial_curvature(t,q,u),atol=1e-9)
    step=1e-5
    for coordinate in (0,3,4):
        qp=q.copy();qm=q.copy();up=u.copy();um=u.copy()
        if coordinate<4:qp[coordinate]+=step;qm[coordinate]-=step
        else:up[:,coordinate-4]+=step;um[:,coordinate-4]-=step
        _,plus=ensemble_sensitivities(t,qp,up,1e-11,1e-13)
        _,minus=ensemble_sensitivities(t,qm,um,1e-11,1e-13)
        np.testing.assert_allclose((plus-minus)/(2*step),H[...,coordinate],rtol=2e-6,atol=2e-8)


def test_block_step_matches_dense_joint_normal_equation():
    rng=np.random.default_rng(193);W=rng.normal(size=(3,20,6));r=rng.normal(size=(3,20))
    I=np.einsum('rmi,rmj->rij',W,W);g=np.einsum('rmi,rm->ri',W,r)
    dq,du,S,_=block_step(I,g)
    dense=np.zeros((60,10))
    for k in range(3):dense[k*20:(k+1)*20,:4]=W[k,:,:4];dense[k*20:(k+1)*20,4+2*k:6+2*k]=W[k,:,4:]
    delta=np.linalg.lstsq(dense,-r.ravel(),rcond=None)[0]
    np.testing.assert_allclose(dq,delta[:4],atol=1e-13);np.testing.assert_allclose(du.ravel(),delta[4:],atol=1e-13)
    np.testing.assert_allclose(np.linalg.inv(S),np.linalg.inv(dense.T@dense)[:4,:4],atol=1e-13)


def test_ensemble_sensitivities_and_whitening_match_single_records():
    t=np.linspace(0,1,31);q=[1,-1,1,-1];u=np.array([[2.,0.],[1.8,.2]])
    x,U=ensemble_sensitivities(t,q,u,1e-10,1e-12)
    for k in range(2):
        xe,Ue=vdp_sensitivities(t,q,u[k]);np.testing.assert_allclose(x[k],xe,atol=1e-8);np.testing.assert_allclose(U[k],Ue,atol=1e-8)
        np.testing.assert_allclose(whiten_records(U,t[1]-t[0],.3,[.1,.2])[k],ou_whiten(U[k],t[1]-t[0],.3,np.array([.1,.2])),atol=1e-13)
