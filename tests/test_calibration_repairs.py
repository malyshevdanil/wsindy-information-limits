import unittest
import numpy as np
from scipy.stats import multinomial
from src.paired_binary import paired_difference_interval
from src.guarded_gls import guarded_gls


class CalibrationRepairs(unittest.TestCase):
    def test_paired_extremes_have_width_and_swap_symmetry(self):
        a=np.ones(100);b=np.zeros(100)
        forward=paired_difference_interval(a,b);reverse=paired_difference_interval(b,a)
        self.assertGreater(forward['ci95'][0],.9);self.assertLess(forward['ci95'][0],1.)
        self.assertEqual(forward['ci95'][1],1.)
        np.testing.assert_allclose(reverse['ci95'],-np.array(forward['ci95'])[::-1])
        equal=paired_difference_interval(a,a)
        self.assertLess(equal['ci95'][0],0);self.assertGreater(equal['ci95'][1],0)

    def test_exact_small_sample_coverage_by_enumeration(self):
        # Deterministic enumeration, not a simulation or asymptotic check.
        n=6
        for p10,p01 in [(0.,0.),(1.,0.),(.1,.2),(.7,.1),(.49,.49)]:
            coverage=0.;delta=p10-p01
            for n10 in range(n+1):
                for n01 in range(n-n10+1):
                    first=np.array([1]*n10+[0]*n01+[0]*(n-n10-n01))
                    second=np.array([0]*n10+[1]*n01+[0]*(n-n10-n01))
                    lo,hi=paired_difference_interval(first,second)['ci95']
                    if lo-1e-12<=delta<=hi+1e-12:
                        coverage+=multinomial.pmf([n10,n01,n-n10-n01],n,[p10,p01,1-p10-p01])
            self.assertGreaterEqual(coverage,.95-1e-12)

    def test_guarded_exact_fit_and_ridge_scaling(self):
        A=np.array([[1.,0.],[0.,1.],[1.,1.]]);truth=np.array([.2,-.4]);bounds=[[-1.,1.],[-1.,1.]]
        for epsilon in [.1,.01]:
            estimate,diag=guarded_gls(A,A@truth,epsilon**2*np.eye(3),epsilon,bounds)
            np.testing.assert_allclose(estimate,truth,atol=1e-12)
            self.assertAlmostEqual(diag['physical_covariance_ridge'],epsilon**4)
            self.assertFalse(diag['fallback'])
        estimate,diag=guarded_gls(np.ones((3,2)),np.ones(3),np.eye(3)*.01,.1,bounds)
        self.assertTrue(diag['fallback']);np.testing.assert_equal(estimate,[0.,0.])
        estimate,_=guarded_gls(A,A@np.array([20.,-30.]),np.eye(3)*.01,.1,bounds)
        np.testing.assert_equal(estimate,[1.,-1.])

    def test_fixed_relative_ridge_limit_and_covariance_identity(self):
        from src.weak_efficiency import gls_map
        A=np.array([[1.],[2.],[3.]])
        shape=np.diag([1.,4.,9.])
        ideal,_=gls_map(A,shape,0.)
        regularized,_=gls_map(A,shape,.2)
        scaled,_=gls_map(A,.01*shape,.2)
        np.testing.assert_allclose(regularized,scaled,atol=1e-12)
        self.assertGreater(np.linalg.norm(regularized-ideal),.001)
        difference=regularized-ideal
        excess=regularized@shape@regularized.T-ideal@shape@ideal.T
        np.testing.assert_allclose(excess,difference@shape@difference.T,atol=1e-12)
        self.assertGreater(float(excess[0,0]),0)

    def test_full_pipeline_failed_pilot_returns_bounded_fallback(self):
        from src.guarded_gls import fit_guarded_polynomial
        t=np.array([0.,.5,1.]);y=np.array([[1.],[.8],[.6]])
        phi=np.array([[0.,1.,0.]]);der=np.array([[1.,0.,-1.]])
        called=[]
        def failed_model(t,pilot):
            called.append(True)
            raise RuntimeError('Deliberate test fixture, no ODE solve')
        estimate,diag=fit_guarded_polynomial(t,y,phi,der,phi,der,[(1,)],
            .05,[1.],.3,[0],[[-2.,2.]],[[-3.,3.]],failed_model)
        self.assertEqual(len(called),1)
        self.assertTrue(diag['fallback']);np.testing.assert_equal(estimate,[0.])


if __name__=='__main__':unittest.main()
