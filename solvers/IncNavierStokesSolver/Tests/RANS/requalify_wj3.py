#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned, read-only gate evaluation of retained WJ-3 publication evidence."""
import argparse
from pathlib import Path
from wj_qualification import wj3_gates
from wj_immutable import read,artifacts,file_hash,write_once,digest,verify_artifacts

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    inputs=artifacts([a.input/(policy+'.json') for policy in ('EP','Budget')]+[a.input/name for name in ('publication.json','receipt.json','manifest.json')])
    result=dict(schema='WJ3RequalificationV1',inputs=inputs,qualification_code=artifacts([Path(__file__),Path(__file__).with_name('wj_qualification.py')]),
        policies={policy:wj3_gates(read(a.input/(policy+'.json'))) for policy in ('EP','Budget')},governing_solves=0,original_receipts_unchanged=True)
    write_once(a.output/('wj3-requalified-'+digest(result['qualification_code'])[:12]+'.json'),result);verify_artifacts(inputs)
    print({p:r['status'] for p,r in result['policies'].items()})

if __name__=='__main__':main()
