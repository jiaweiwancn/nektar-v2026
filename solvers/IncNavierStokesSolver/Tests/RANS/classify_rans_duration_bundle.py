#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CLI aggregate classifier: required numerical/coverage failures return nonzero."""
import argparse,json
from pathlib import Path
from rans_duration_contract import aggregate


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--request',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    try:result=aggregate(json.loads(a.request.read_text()),json.loads(a.receipt.read_text()))
    except Exception as e:result=dict(status='fail',failure=repr(e))
    a.output.write_text(json.dumps(result,indent=2)+'\n');return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
