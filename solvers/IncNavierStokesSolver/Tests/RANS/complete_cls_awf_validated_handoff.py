#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete a prepared handoff after infrastructure failure, without a PDE replay.

Requires the exact successful checkpoint audit and isolated-restart proof,
unchanged parent bytes, stopped old owners, and a positively acquired lock.
It never accepts a numerical rejection or invents a prefix/history.
"""
import argparse
from pathlib import Path
import shutil
from handoff_cls_awf_runtime import launch_validated_successor,observe_successor,matches
from continue_cls_epawf_pimple import atomic,header,process_identity,assert_native_stopped
from run_cls_epawf_picard import load,sha,require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--block-steps',type=int,choices=(8,32),default=32)
    p.add_argument('--observe-existing',action='store_true',help='Only finish the receipt for an already launched, sealed successor; never launch or signal a process')
    p.add_argument('--timeout',type=float,default=7200);p.add_argument('--targets',type=int,nargs='+',default=[32,160,800,1600,3200,8000]);a=p.parse_args()
    source=a.output/'handoff.json';out=load(source);a.campaign=Path(out['parent_campaign']);old=load(a.campaign/'progress.json')
    require(sha(a.solver)==out['solver_sha256'],'validated successor executable differs')
    require(out['complete_corrected_audit']['accepted_awf_hard_admissible'][0]==1 and
        out['isolated_restart_equivalence']['bitwise_checkpoint_payload'] and out['verification_governing_sweeps']==2,'missing successful audit/restart proof')
    for key in ('old_driver_identity','old_supervisor_identity'):
        expected=out[key];require(not matches(process_identity(expected['pid']),expected),'old owner still alive')
    assert_native_stopped(old);inputs=a.output/'parent';h=header(inputs/'checkpoint',True)
    require(h['sha256']==out['parent_checkpoint_sha256']==old['checkpoint_sha256'] and h['steps']==old['accepted_steps'],'validated parent checkpoint differs')
    def save(**changes):out.update(changes);atomic(source,out)
    if a.observe_existing:
        require(out['status']=='successor_launched','only an already launched successor is eligible for observation')
        seal=load(a.campaign/'runtime-handoff.json')
        require(seal['successor']==str(a.output) and seal['parent_checkpoint_sha256']==h['sha256'] and
            seal['successor_supervisor_identity']==out['successor_supervisor_identity'],'successor seal differs')
        require(matches(process_identity(out['successor_supervisor_identity']['pid']),out['successor_supervisor_identity']),'recorded successor supervisor is not live')
        current=load(a.output/'trajectory/progress.json')
        require(current['executable_sha256']==out['solver_sha256'] and current['prefix_steps']==h['steps'] and
            current['audit_revision']=='NativeElementSideV2' and current['evaluation_reuse']=='ExactStepV1','live successor configuration differs')
        observe_successor(a,h,out['expected_first_successor_checkpoint_sha256'],save);return
    require(out['status']=='handoff_failed' and out.get('old_owners_terminated') is True,'only a validated retired-owner infrastructure interruption is eligible')
    require(not (a.output/'trajectory').exists() and not (a.campaign/'runtime-handoff.json').exists(),'successor/seal already exists; inspect actual state')
    copy=a.output/'handoff-before-infrastructure-completion.json';require(not copy.exists(),'completion was already attempted');shutil.copy2(source,copy)
    save(status='resuming_validated_handoff',previous_infrastructure_failure=out.get('failure'),failure=None)
    launch_validated_successor(a,old,h,inputs,out['expected_first_successor_checkpoint_sha256'],save)

if __name__=='__main__':main()
