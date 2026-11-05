#!/usr/bin/env python3
"""Regression tests for binding the actual selected R11 native candidate."""
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('r11',Path(__file__).with_name('experiment.py'))
r11=importlib.util.module_from_spec(spec)
spec.loader.exec_module(r11)


def trial(f, state, residual, accepted=True):
    return dict(comparison='Nf',fraction=f,state={'sha256':state},residual={'sha256':residual},
                classification='accepted' if accepted else 'residual_non_descent')


class Selection(unittest.TestCase):
    def test_zero(self):
        self.assertIsNone(r11.selected_trial(dict(name='Nf',accepted=False),[]))

    def test_one(self):
        t=trial(.125,'a','r')
        p=dict(name='Nf',accepted=True,fraction=.125,selected_state=t['state'],selected_residual=t['residual'])
        self.assertIs(r11.selected_trial(p,[t]),t)

    def test_multiple_bind_selected_not_first(self):
        first=trial(.125,'a','r')
        chosen=trial(.03125,'b','s')
        last=trial(.0625,'c','t')
        p=dict(name='Nf',accepted=True,fraction=.03125,selected_state=chosen['state'],selected_residual=chosen['residual'])
        self.assertIs(r11.selected_trial(p,[first,chosen,last]),chosen)

    def test_changed_residual_rejected(self):
        t=trial(.125,'a','r')
        p=dict(name='Nf',accepted=True,fraction=.125,selected_state=t['state'],selected_residual={'sha256':'different'})
        with self.assertRaises(ValueError):
            r11.selected_trial(p,[t])

    def test_rejected_selected_state_rejected(self):
        t=trial(.125,'a','r',False)
        p=dict(name='Nf',accepted=True,fraction=.125,selected_state=t['state'],selected_residual=t['residual'])
        with self.assertRaises(ValueError):
            r11.selected_trial(p,[t])


if __name__=='__main__':
    unittest.main()
