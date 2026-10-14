import numpy as np
from src.crossfit_gmm import crossfit_equations,stratified_counts
from src.heterogeneous_gmm import record_fit
from src.crossfit_cache import features,powers


def test_polynomial_cache_matches_explicit_record_equations():
    A,b,B,a,folds=fixture();rng=np.random.default_rng(121)
    exponents=powers(2);anchor=np.zeros(4);radius=2.
    B=np.zeros((len(exponents),20,4,4));B[0]=np.eye(4)
    a=rng.normal(size=(len(exponents),20,4))*.01
    q,V,p=crossfit_equations(A,b,B,a,anchor,[[-100,100]]*4,folds,cache_exponents=exponents,cache_radius=radius)
    weights=features(p[0,folds],anchor,radius,exponents)
    Be=np.einsum('rp,prij->rij',weights,B);ae=np.einsum('rp,pri->ri',weights,a)
    qe,Ve=record_fit(Be,ae)
    np.testing.assert_allclose(q[0],qe);np.testing.assert_allclose(V[0],Ve)
    np.testing.assert_array_equal(features([100,-100,0,0],anchor,radius,exponents),features([2,-2,0,0],anchor,radius,exponents))
    counts=stratified_counts(rng,folds,3)
    qb,Vb,pb=crossfit_equations(A,b,B,a,anchor,[[-100,100]]*4,folds,counts,exponents,radius)
    for draw,c in enumerate(counts):
        ids=np.repeat(np.arange(20),c)
        qe,Ve,pe=crossfit_equations(A[ids],b[ids],B[:,ids],a[:,ids],anchor,[[-100,100]]*4,folds[ids],cache_exponents=exponents,cache_radius=radius)
        np.testing.assert_allclose(qb[draw],qe[0],atol=1e-13)
        np.testing.assert_allclose(Vb[draw],Ve[0],atol=1e-13)
        np.testing.assert_allclose(pb[draw],pe[0],atol=1e-13)


def fixture():
    rng=np.random.default_rng(871)
    A=np.tile(np.eye(4),(20,1,1));b=rng.normal(size=(20,4))
    B=np.zeros((5,20,4,4));B[0]=np.eye(4)
    a=np.zeros((5,20,4));a[0]=rng.normal(size=(20,4));a[1:]=rng.normal(size=(4,20,4))*.1
    return A,b,B,a,np.arange(20)%5


def test_holdout_copies_excluded_and_explicit_estimate_matches():
    A,b,B,a,folds=fixture();anchor=np.zeros(4);box=[[-100,100]]*4
    q,V,pilots=crossfit_equations(A,b,B,a,anchor,box,folds)
    explicit_a=np.empty((20,4))
    for k in range(5):
        expected=b[folds!=k].mean(axis=0)
        np.testing.assert_allclose(pilots[0,k],expected)
        ids=folds==k;explicit_a[ids]=a[0,ids]+np.einsum('p,pri->ri',expected,a[1:,ids])
    qe,Ve=record_fit(B[0],explicit_a)
    np.testing.assert_allclose(q[0],qe);np.testing.assert_allclose(V[0],Ve)
    changed=b.copy();changed[folds==0]+=1000
    _,_,p2=crossfit_equations(A,changed,B,a,anchor,box,folds)
    np.testing.assert_allclose(pilots[0,0],p2[0,0],rtol=0,atol=1e-12)


def test_bootstrap_batch_equals_expanded_record_calculation():
    A,b,B,a,folds=fixture();counts=stratified_counts(np.random.default_rng(17),folds,3)
    q,V,p=crossfit_equations(A,b,B,a,np.zeros(4),[[-100,100]]*4,folds,counts)
    for draw,c in enumerate(counts):
        ids=np.repeat(np.arange(20),c)
        qe,Ve,pe=crossfit_equations(A[ids],b[ids],B[:,ids],a[:,ids],np.zeros(4),[[-100,100]]*4,folds[ids])
        np.testing.assert_allclose(q[draw],qe[0],atol=1e-13)
        np.testing.assert_allclose(V[draw],Ve[0],atol=1e-13)
        np.testing.assert_allclose(p[draw],pe[0],atol=1e-13)
        for k in range(5):assert c[folds==k].sum()==4
