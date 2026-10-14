"""Algebra and basis checks; no trajectory simulation."""
import unittest
import numpy as np
from src.self_instrumented_gls import interpolation_kernel,dyadic_cubic_tests

class ConstructionChecks(unittest.TestCase):
    def test_exact_kernel_second_moment_in_smoothing_lemma(self):
        import sympy as sp
        s=sp.symbols('s')
        K=sp.Rational(315,256)*(1-s**2)**4
        self.assertEqual(sp.integrate(K,(s,-1,1)),1)
        self.assertEqual(sp.integrate(s*K,(s,-1,1)),0)
        self.assertEqual(sp.integrate(s**2*K,(s,-1,1)),sp.Rational(1,11))

    def test_interpolation_reproduces_constant_and_linear_derivative(self):
        for h,b in [(.005,.147),(.001,.1)]:
            k,d,m=interpolation_kernel(h,b);lag=np.arange(-m,m+1)*h
            self.assertAlmostEqual(k.sum(),1.,places=13)
            self.assertAlmostEqual(d.sum(),0.,places=12)
            self.assertAlmostEqual(d@(-lag),1.,places=12)
            np.testing.assert_allclose(k,k[::-1],atol=1e-14)
            np.testing.assert_allclose(d,-d[::-1],atol=1e-14)

    def test_nested_cubic_refinement_and_boundary_traces(self):
        t=np.linspace(0,20,2001)
        coarse,_=dyadic_cubic_tests(t,4);fine,der=dyadic_cubic_tests(t,5)
        mask=np.array([1,4,6,4,1])/8
        for index,j in enumerate(range(2,15)):
            inds=np.arange(2*j-2,2*j+3)-2
            np.testing.assert_allclose(coarse[index],mask@fine[inds],atol=3e-15)
        self.assertTrue(np.all(fine[:,[0,-1]]==0))
        self.assertTrue(np.all(der[:,[0,-1]]==0))

if __name__=='__main__':unittest.main()
