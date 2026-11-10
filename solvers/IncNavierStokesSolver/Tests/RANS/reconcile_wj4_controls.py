#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Append missing conservative legacy-control units; never rewrite receipts."""
import argparse
from pathlib import Path
import re
import time
from wj_immutable import Ledger,read,write_once,artifacts,verify_artifacts,file_hash
from wj4_validate import legacy_control_budget


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
    output=root/'legacy-control-accounting-supplement.json'
    if output.exists():verify_artifacts(read(output)['inputs']);print(output);return
    rows=[];inputs=[];commands=[]
    for directory,label in ((root,'qualification'),(root/'readback-runtime','observer2')):
        inventory=directory/(label+'-controls-inventory.json');accounting=directory/(label+'-controls-accounting.json')
        if read(accounting)['attempt_index']:raise ValueError('controls already have budget units; do not double charge')
        command=next(r for r in read(directory/'command-receipts.json') if r['label']=='ctest-'+label+'-controls')
        log=Path(command['log'])
        if command['exit_code']!=0 or file_hash(log)!=command['log_sha256'] or not re.search(r'100% tests passed, 0 tests failed out of 24',log.read_text()):
            raise ValueError('successful original complete controls are required')
        value=legacy_control_budget(read(inventory)['tests'])
        rows.extend(dict(runtime_root=str(directory),**r) for r in value);commands.append(command)
        inputs.extend([inventory,accounting,log,*[Path(r['xml']) for r in value]])
    units=sum(r['attempts'] for r in rows);before=artifacts(inputs);ledgers=[Ledger(root,p) for p in ('EP','Budget')]
    label=str(output);start=time.time();reserved=[]
    try:
        for l in ledgers:l.reserve(label,dict(inputs=before,code=artifacts([Path(__file__)])),units,30);reserved.append(l)
        write_once(output,dict(schema='WJ4LegacyControlSupplementV1',inputs=before,original_command_entries=commands,controls=rows,
            extra_conservative_units_per_policy=units,native_WJ_sweep_attempts_added=0,
            scope='ordinary legacy transport/evolution timesteps and one standalone stage invocation; prior wall time already charged; original receipts/reservations unchanged; future wrapper reserves these units before launch'))
        verify_artifacts(before)
    finally:
        for l in reserved:l.finish(label,units if output.exists() else 0,time.time()-start,str(output))
    print(output)


if __name__=='__main__':main()
