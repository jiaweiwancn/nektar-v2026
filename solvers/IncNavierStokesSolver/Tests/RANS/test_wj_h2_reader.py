#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Passive producer/earlier-checkpoint fixtures on retained H2 native controls."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from measure_wj_h2 import checkpoint_for
from validate_wj_h1 import run, save, sha, require
from run_wj_h2 import load


def vector_hash(values):
    return hashlib.sha256(np.asarray(values, dtype='<f8').tobytes()).hexdigest()


def check(probe, summary, output):
    output=Path(output);output.mkdir(exist_ok=False)
    original=load(summary);runs=original['runs']
    rp=Path(runs['two']['report']);report=load(rp);xml=rp.with_name('case.xml')
    record=dict(report=str(rp),xml=str(xml))
    cp,payload=checkpoint_for(record,1)
    packet=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),
        native_report_sha256=sha(rp),native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),
        physical_observations=[[4.,0.,1.]],assembly=False,weak=False,reuse=False)
    immutable={str(p):sha(p) for p in (rp,Path(cp)/'State.bin')}
    receipts={}
    def observe(label,request,reason=None):
        folder=output/label;folder.mkdir();inp=folder/'packet.json';out=folder/'result.json';save(inp,request)
        receipts[label]=run([Path(probe).resolve(),xml,inp,out],folder,timeout=180,expected_failure=reason)
        return json.loads(out.read_text()) if reason is None else None
    earlier=observe('earlier',packet)
    audit=report['wj_audit_packets'][1]
    require(earlier['selected_absolute_step']==1 and earlier['producer_absolute_step']==2,'wrong earlier producer index')
    require(earlier['state_sha256']==vector_hash(payload['state']) and
        earlier['callback_history_sha256']==audit['stage_history_sha256'],'earlier callback/state differs')
    require(earlier['initialization_work']==[0,0] and earlier['governing_solves']==0,'reader performed governing work')
    require(earlier['checkpoint_native_history_sha256']==[vector_hash(v) for v in payload['native_history']],
        'reader did not retain native checkpoint histories')
    bad=dict(packet,native_report_sha256='0'*64)
    observe('wrong-report-hash',bad,'open producing report SHA256 differs')
    other=Path(runs['one']['report']);bad=dict(packet,native_report=str(other),native_report_sha256=sha(other),
        native_checkpoint=report['checkpoint'],checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin'))
    observe('wrong-producer',bad,'open checkpoint outside producer accepted lineage')
    tampered=copy.deepcopy(report);tampered['wj_audit_packets'][1]['callback_history'][0]+=1.
    fake=output/'altered-callback.json';save(fake,tampered)
    observe('wrong-callback',dict(packet,native_report=str(fake),native_report_sha256=sha(fake)),
        'open accepted audit/callback binding differs')
    tampered=copy.deepcopy(report);tampered['native_history'][0][0]+=1.
    fake=output/'altered-native-history.json';save(fake,tampered)
    terminal=dict(packet,native_report=str(fake),native_report_sha256=sha(fake),native_checkpoint=report['checkpoint'],
        checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin'))
    observe('wrong-native-history',terminal,'open native GLM history differs')
    failed=Path(runs['rejected']['report']);d=load(failed)
    predecessor=dict(packet,native_report=str(failed),native_report_sha256=sha(failed),
        native_checkpoint=d['checkpoint'],checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'))
    observe('failed-role-required',predecessor,'requires an explicit diagnostic role')
    pred=observe('accepted-predecessor',dict(predecessor,allow_failed_producer=True))
    require(pred['state_sha256']==earlier['state_sha256'] and pred['producer_failure'],
        'failed producer promoted rejected state or lost diagnostic role')
    require(all(sha(p)==h for p,h in immutable.items()),'passive fixtures altered original producer bytes')
    result=dict(status='PASS',governing_attempts=0,original_summary=str(summary),original_summary_sha256=sha(summary),
        selected_earlier_step=1,producer_step=2,unchanged=immutable,receipts=receipts)
    save(output/'summary.json',result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--control-summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();check(a.probe,a.control_summary,a.output);print(json.dumps(dict(status='PASS',governing_attempts=0)))
if __name__=='__main__':main()
