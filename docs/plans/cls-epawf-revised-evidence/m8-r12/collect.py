#!/usr/bin/env python3
"""Bind the R12 selector's exact chosen tuple and unchanged native history."""
import argparse
import hashlib
import importlib.util
import struct
from pathlib import Path
import platform
import sys
sys.dont_write_bytecode=True

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert platform.system()=='Linux'
    e=a.epoch.resolve();evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r12_collect')
    receipt=h.load(e/'comparison-run/receipt.json')
    assert receipt['status']=='completed' and receipt['runtime_unchanged']
    C=h.load(h.bound(receipt['durable_collection']['file']))
    assert C['terminal']['terminal']=='r12_fixed_comparison_completed'
    base=C['completed_iterates'][0]
    assert base['state']['sha256']=='eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c'
    assert base['history']['sha256']=='b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71'
    matrix,=[x for x in C['jacobians'] if x['status']=='complete']
    assert matrix['event']['matrix']['sha256']=='2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6'
    events=module(evidence/'m8-r7/experiment.py','r7_events').events(C)
    selector=module(evidence/'m8-r11/experiment.py','r11_selector')
    trials=[x for x in events if x['kind']=='r8_trial'];policies=[x for x in events if x['kind']=='r11_policy']
    chosen={x['name']:selector.selected_trial(x,trials) for x in policies}
    decision,=[x for x in events if x['kind']=='r11_decision'];last=C['last_snapshot']
    retained=chosen[decision['retained_policy']]
    assert all(retained[k]['sha256']==last[k]['sha256'] for k in ('state','residual'))
    report=h.load(h.bound(receipt['native_report']))
    def digest(v):return hashlib.sha256(struct.pack('<'+str(len(v))+'d',*v)).hexdigest()
    start=report['stage_test_initial_native_history'];h0=base['history']['sha256']
    assert report['time']==report['steps']==0 and digest(report['state'])==digest(start['state'])==h0
    assert [digest(v) for v in start['history']]==[digest(v) for v in report['native_history']]
    assert start['times']==report['native_history_times']
    level,=[x for x in report['stage_test_levels'] if x.get('attempted') and digest(x['state'])==last['state']['sha256']]
    assert digest(level['history'])==h0 and digest(level['residual'])==last['residual']['sha256'] and not level['native_history_advanced']
    assert all(not s['history_advanced'] for s in C['completed_iterates'])
    h.save(a.output,dict(receipt=h.pin(e/'comparison-run/receipt.json'),collection=receipt['durable_collection']['file'],
        initial=base,selected=last,matrix=matrix['event']['matrix'],events=events,policies=policies,chosen=chosen,
        decision=decision,native_report=receipt['native_report'],native_history_verified=True))

if __name__=='__main__':main()
