#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent integral fixtures for external versus shared-face accounting."""
import copy
from types import SimpleNamespace
import unittest
from analyze_rans_physical_budgets import whole_patch


class PhysicalBudgetAccounting(unittest.TestCase):
    def setUp(self):
        self.ref=SimpleNamespace(K=2.,us=1.,H=2.)
        def part(storage, source, flux):
            return dict(volume=4.,storage_backward_difference=storage,stage_source=source,
                instantaneous_source=source,stage_physical_boundary_flux=flux,
                instantaneous_physical_boundary_flux=flux)
        self.receipt=dict(physical_budget=dict(equations=dict(
            inner_k=part(10.,5.,3.),outer_k=part(20.,11.,4.)),physical_faces=[
            dict(axis=2,coordinate=[0.,0.,2.],weight=1.,k_frozen_flux=1.,k_current_flux=1.),
            dict(axis=2,coordinate=[0.,0.,2.],weight=1.,k_frozen_flux=-2.,k_current_flux=-2.)]))

    def test_external_integral_defect(self):
        # Storage 30, source 16, regional flux 7 with internal sum -1:
        # physical external flux is 8, hence defect 6 and rate scale 8.
        row=whole_patch(self.receipt,self.ref)
        self.assertEqual(row['stage_physical_boundary_flux'],1.)
        self.assertEqual(row['frozen_stage_physical_defect'],.75)
        self.assertEqual(row['interface_frozen_outward_flux_sum'],-.125)

    def test_internal_flux_is_counted_once(self):
        changed=copy.deepcopy(self.receipt)
        q=changed['physical_budget'];q['physical_faces'][0]['k_frozen_flux']+=1
        q['physical_faces'][0]['k_current_flux']+=1
        q['equations']['inner_k']['stage_physical_boundary_flux']+=1
        q['equations']['inner_k']['instantaneous_physical_boundary_flux']+=1
        a,b=whole_patch(self.receipt,self.ref),whole_patch(changed,self.ref)
        self.assertEqual(a['frozen_stage_physical_defect'],b['frozen_stage_physical_defect'])
        self.assertNotEqual(a['interface_frozen_outward_flux_sum'],b['interface_frozen_outward_flux_sum'])

    def test_physical_flux_changes_defect_but_reaction_does_not(self):
        changed=copy.deepcopy(self.receipt)
        changed['prescribed_row_reaction']=1e99
        self.assertEqual(whole_patch(self.receipt,self.ref),whole_patch(changed,self.ref))
        changed['physical_budget']['equations']['outer_k']['stage_physical_boundary_flux']+=1
        self.assertEqual(whole_patch(changed,self.ref)['frozen_stage_physical_defect'],.625)


if __name__=='__main__':unittest.main()
