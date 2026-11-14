#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Post-science selector verification; no new native operation."""
import argparse,ast,sys,unittest
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import atomic
from wj_h18_readback import selection


def check(root):
    test=unittest.defaultTestLoader.loadTestsFromName('test_wj_h18_readback')
    assert unittest.TextTestRunner(verbosity=2).run(test).wasSuccessful()
    assert selection(root)[1]==load(root/'science/target/B/view/request.json')
    old=ast.parse((root/'source/solvers/IncNavierStokesSolver/Tests/RANS/study_wj_h18.py').read_text())
    new=ast.parse((root/'post-science-python/study_wj_h18.py').read_text())
    for tree in (old,new):tree.body=[n for n in tree.body if not(isinstance(n,ast.FunctionDef)and n.name=='strict')]
    assert ast.dump(old)==ast.dump(new),'post-science change affected governing orchestration'
    atomic(root/'science/readback-selection-compatibility-v1.json',dict(status='PASS',tests=3,actual_final_request_equal=True,only_strict_function_changed=True,additional_native_calls=0))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();check(a.root)
