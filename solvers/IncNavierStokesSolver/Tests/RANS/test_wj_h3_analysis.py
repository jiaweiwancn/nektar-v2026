# SPDX-License-Identifier: MIT
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj_h2_profiles import constitutive
from wj_h3_analysis import algebra,principal


class PassiveAlgebra(unittest.TestCase):
    def test_full_three_dimensional_batch_matches_independent_scalar_algebra(self):
        rng=np.random.default_rng(3147);G=rng.normal(size=(11,3,3));k=rng.uniform(.1,2,11);e=rng.uniform(.1,2,11)
        R,nu,N=algebra(k,e,G)
        for i in range(len(k)):
            r,n,m=constitutive(k[i],e[i],G[i]);assert_allclose(R[i],r,rtol=2e-12,atol=1e-13)
            assert_allclose([nu[i],N[i]],[n,m],rtol=2e-12,atol=1e-13)

    def test_rank_one_quartic_and_trace_invariance(self):
        G=np.zeros((3,3));G[0,2]=3.5
        p=principal(1.,1.,G)
        self.assertLess(p['in_plane'],0)
        xi=np.array([1.,0,-7])/np.sqrt(50);a=np.array([7.,0,1])/np.sqrt(50);h=np.outer(a,xi)
        R,nu,N=constitutive(1.,1.,G)
        witness=-np.sum(h*constitutive(1.,1.,G.astype(complex)+1e-27j*h,N)[0].imag)/1e-27/nu
        assert_allclose(witness,-.00633378103539413,atol=1e-12,rtol=0)
        self.assertLessEqual(p['in_plane'],witness+1e-14)
        q=principal(1.,1.,G+np.eye(3)*.7)
        assert_allclose([q['in_plane'],q['transverse']],[p['in_plane'],p['transverse']],atol=2e-14,rtol=0)
        assert_allclose(principal(1.,1.,np.zeros((3,3)))['in_plane'],1.,atol=2e-14,rtol=0)

if __name__=='__main__':unittest.main()
