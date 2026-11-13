#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused invariant-branch checks for the phase-I trial tensor evaluator."""
import os,unittest
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from search_cao_tamura_budget_h2d import tensors
from validate_cls_epawf_m1 import indexed

class BranchAwareTrial(unittest.TestCase):
    def test_fresh_scalar_and_max_branch(self):
        for k,e in ((.0176,.008),(.0352,.004),(.0088,.016)):
            for rotation in (-1e-5,0,1e-5):
                g=np.zeros((3,3));g[0,2]=1.3;g[2,0]=rotation
                got=tensors(g[None],np.array([k]),np.array([e]),np.array([False]),np.array([0.]),np.array([[.03,0,.17]]),.0743,.440553905186659)[0]
                c=indexed(g,k,e,'physical');S=g+g.T;sr=np.zeros((3,3));sr[0,2]=sr[2,0]=.0743/(.440553905186659*.17)
                want=2*k*np.eye(3)/3-float(c['nut'].real)*S+c['quadratic'].real+.4*c['cubic'].real-.440553905186659*.0743*.17*(S-sr)
                self.assertLess(np.max(abs(got-want)),1e-13)
    def test_one_sided_derivatives_at_switch(self):
        # Pure shear has q_S=q_Omega. Use the actual max on each side;
        # no unsupported common central derivative is asserted at the switch.
        g=np.zeros((3,3));g[0,2]=1.3;d=np.zeros((3,3));d[2,0]=1.
        def R(h):return tensors((g+h*d)[None],np.array([.0176]),np.array([.008]),np.array([False]),np.array([0.]),np.array([[.03,0,.17]]),.0743,.440553905186659)[0]
        for sign in (-1,1):
            h=sign*1e-6;coarse=(R(h)-R(0))/h;fine=(R(h/2)-R(0))/(h/2)
            self.assertLess(np.max(abs(coarse-fine)),1e-7)
if __name__=='__main__':unittest.main()
