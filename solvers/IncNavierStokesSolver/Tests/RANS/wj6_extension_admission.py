#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Write-once joint T3 admission; unlaunched phases never become PASS."""
import argparse
import json
from pathlib import Path
import time
from wj_immutable import Ledger,artifacts,atomic,file_identity,file_hash,write_once
from wj6_qualification import common_extension_authorized
from wj_qualification import component


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--study',type=Path,required=True);a=p.parse_args()
    started=time.time();root=a.root.resolve();summaries={n:json.loads((a.study/n/'summary.json').read_text()) for n in ('EP','Budget')}
    inputs={};states={}
    for policy,s in summaries.items():
        for phase in ('primary','transactions'):
            if s['components'][phase]['status']!='PASS':raise ValueError('joint extension admission requires completed Tx/transaction qualification')
        path=a.study/policy/'T3-resource-projection.json';projection=json.loads(path.read_text())
        if s['components']['extension_resource']['projection']['sha256']!=file_hash(path):raise ValueError('extension projection pointer differs')
        inputs[policy]=file_identity(path);states[policy]=projection['status']
    allowed=common_extension_authorized({p:s['components'] for p,s in summaries.items()})
    status='PASS' if allowed else 'RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in states.values() else 'UNRESOLVED'
    output=root/'extension-admission.json'
    write_once(output,dict(schema='WJ6JointExtensionAdmissionV1',status=status,authorized=allowed,
        resource_states=states,projections=inputs,frozen_target=dict(T=.004096,coarse=128,fine=256,baseline=128),
        governing_attempts=0,code=artifacts([Path(__file__)])))
    for policy,s in summaries.items():
        ledger=Ledger(root,policy,cap=800);label='joint-T3-admission'
        ledger.reserve(label,dict(admission=file_identity(output)),0,30)
        if not allowed:s['components']['extension']=component(status,'frozen T3 target retained; no extension launched because joint resource admission did not pass',admission=file_identity(output))
        s['extension_admission']=file_identity(output);atomic(a.study/policy/'summary.json',s)
        ledger.finish(label,0,time.time()-started,str(output))
    print(output)


if __name__=='__main__':main()
