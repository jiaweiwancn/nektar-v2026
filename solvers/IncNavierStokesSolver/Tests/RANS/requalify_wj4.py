#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Requalify immutable native observations; no probe or governing launch."""
import argparse
import os
from pathlib import Path
import time
from wj_immutable import Ledger,read,write_once,artifacts,verify_artifacts
from wj4_qualification import wj4_gates

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    start=time.time();root=a.study.resolve()/a.policy;scripts=Path(__file__).resolve().parent
    inputs=artifacts([root/'summary.json',root/'design.json',Path(__file__),scripts/'wj4_qualification.py',scripts/'wj_qualification.py'])
    ledger=Ledger(Path(os.environ['WJ4_PACKAGE_ROOT']),a.policy);label=str(a.output.resolve());ledger.reserve(label,dict(inputs=inputs),0,30)
    try:
        result=wj4_gates(read(root/'summary.json'),read(root/'design.json'),ledger.totals());verify_artifacts(inputs)
        write_once(a.output,dict(schema='WJ4CompleteGateReceiptV1',policy=a.policy,inputs=inputs,governing_attempts=0,qualification=result))
    finally:ledger.finish(label,0,time.time()-start,str(a.output))
    print(a.output)

if __name__=='__main__':main()
