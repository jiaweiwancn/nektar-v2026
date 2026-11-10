#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Conservative live scientific-process thread inventory before allocation."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    p=argparse.ArgumentParser();p.add_argument('--reserve',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rows=[]
    for line in subprocess.check_output(['ps','-eo','pid=,ppid=,nlwp=,pcpu=,args='],text=True).splitlines():
        pid,parent,n,cpu,args=line.strip().split(None,4)
        command=args.split()[0].split('/')[-1]
        if any(command.startswith(c) for c in ('IncNavierStokes','RANS','python','cc1','orted','mpirun','mpiexec','lto1','f951')):
            rows.append(dict(pid=int(pid),parent=int(parent),threads=int(n),cpu=float(cpu),args=args))
    counted=sum(r['threads'] for r in rows);value=dict(epoch=time.time(),existing_threads=counted,additional_reserved=a.reserve,total=counted+a.reserve,limit=72,allowed=counted+a.reserve<=72,processes=rows)
    if a.output.exists():raise ValueError('resource receipt is write-once')
    a.output.write_text(json.dumps(value,indent=2)+'\n');print(json.dumps({k:v for k,v in value.items() if k!='processes'}))
    if not value['allowed']:raise SystemExit(2)

if __name__=='__main__':main()
