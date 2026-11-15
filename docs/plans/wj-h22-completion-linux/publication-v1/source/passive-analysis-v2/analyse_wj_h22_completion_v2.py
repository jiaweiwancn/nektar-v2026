#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned passive rerun on the unchanged seven-attempt composite index."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS
os.environ.update(THREADS)
from wj_h15_observer import identity
from wj_h22_phase import run

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    run(a.root,out_name='phase-analysis-v2',source_index=identity(a.root/'science/composite-source-index.json'))
