#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded summaries of retained CT8 preparations, without printing packets."""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_request_manifest import sha

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r={}
    for label in ('preflight-v1','preflight-v2','preflight-v3'):
        root=a.root/label;row={}
        f=root/'failure.json'
        if f.exists():row['failure']=json.loads(f.read_text())
        f=root/'operator-p3/operator.json'
        if f.exists():
            q=json.loads(f.read_text());row.update(packet_sha256=sha(f),status=q['status'],physical_admissibility=q.get('physical_admissibility'),failure=q.get('failure'))
            audit=q.get('observation',{}).get('audit_alpha_zero',{})
            row['domain_and_residual_primitives']={k:v for k,v in audit.items() if any(s in k for s in ('covariance','radicand','domain','admissible','minimum','divergence'))}
        r[label]=row
    f=a.root/'reference-v2/results.json';q=json.loads(f.read_text());r['central']={k:q[k] for k in ('central_fit','central_root','central_column','extraction_summary','native')}
    write_json(a.output,r);return 0
if __name__=='__main__':raise SystemExit(main())
