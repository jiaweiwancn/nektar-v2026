# SPDX-License-Identifier: MIT
import itertools,unittest
import numpy as np
from wj_h12_phase import identity
from wj_h22_linear_roundoff import dev_bound,propagated_identity

class LinearRoundoff(unittest.TestCase):
    def test_deviatoric_bound_contains_all_diagonal_sign_corners(self):
        b=np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])*1e-14;limit=dev_bound(b)
        for signs in itertools.product((-1,1),repeat=3):
            e=b.copy()
            for j,s in enumerate(signs):e[j,j]*=s
            d=e-np.eye(3)*np.trace(e)/3
            self.assertTrue(np.all(abs(d)<=limit+1e-28))

    def test_input_certificate_must_propagate_without_changing_128_eps(self):
        e=np.diag([0.,0.,2e-14]);g=identity(e,[e],1.)
        tau=5.;response=tau*(e-np.eye(3)*np.trace(e)/3)
        with self.assertRaises(ValueError):identity(response,[response],1.)
        inherited=abs(tau)*dev_bound(g['bound'])
        proof=propagated_identity(response,[response],inherited,1.)
        self.assertGreater(proof['direct_bound_exceeded_count'],0)
        direct=128*np.finfo(float).eps*(1+abs(response))
        np.testing.assert_array_equal(proof['bound'],direct+inherited)
        with self.assertRaises(ValueError):propagated_identity(response+1e-8,[response],inherited,1.)

if __name__=='__main__':unittest.main()
