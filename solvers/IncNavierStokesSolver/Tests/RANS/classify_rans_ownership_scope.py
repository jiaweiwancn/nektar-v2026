#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Emit an independent required-scope receipt and fail the process on unmet gates."""
import argparse
from pathlib import Path
from qualify_rans_preservation import load
from rans_ownership_contract import classify_scope
from rans_qualification import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:result=classify_scope(load(a.input))
    except (ValueError,KeyError,TypeError) as e:result=dict(status='fail',failure=str(e))
    require_new=not a.output.exists()
    if not require_new:raise FileExistsError('immutable required-scope output already exists')
    write_json(a.output,result);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
