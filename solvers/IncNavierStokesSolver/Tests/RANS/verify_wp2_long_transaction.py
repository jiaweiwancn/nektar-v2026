#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bitwise 4+4 replay of eight steps in the executed WP2 transit trajectory."""
import argparse
from pathlib import Path
from cao_tamura_h2d_state import load
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_epawf_picard import require,sha,state_hash
from run_cao_tamura_budget_h2f_startup import verify_retained_rollback
from wp2_native import Run
from wp2_qualification import verify_report


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','h1','preflight','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    a=p.parse_args();origin=a.root/f'preservation-coarse-{a.model}/step60';ordinary=load(origin/'case.evolution.json')
    parents=[x for x in origin.glob('checkpoint*') if x.is_dir() and checkpoint(x)['steps']==52]
    require(len(parents)==1,'unique actual accepted step52 checkpoint required')
    parent=parents[0];before=sha(parent/'State.bin');final=Path(ordinary['checkpoint']);final_before=sha(final/'State.bin')
    run=Run(a.root,a.solver,a.h1,a.preflight,a.output,a.model,4,16,'RadialEnergyV1')
    K=.017598891430879275 if a.model=='CLS' else .01838575534504143
    try:
        first=run.native('split4',4,.02,parent);verify_report(first,.07426793792419734,K)
        second=run.native('restart4',4,.02,Path(first['checkpoint']));verify_report(second,.07426793792419734,K)
        require(same_history(ordinary,second),'transit-time 4+4 state/history differs from uninterrupted path')
        for key in ['fixed_reference','fixed_defect']:
            require(state_hash(ordinary[key])==state_hash(second[key]),'split changed '+key)
        stages=['history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256']
        for expected,actual in zip(ordinary['pimple_attempts'][-8:],first['pimple_attempts']+second['pimple_attempts']):
            require(all(expected[k]==actual[k] for k in stages),'split governing-stage hashes differ')
        # Compare rollback to the exact restart-initialized accepted context,
        # including its ephemeral stage cache, not an earlier process cache.
        readback=run.native('parent-readback',0,.02,final)
        bad=run.native('post-assignment-reject',1,.02,final,reject=True)
        verify_retained_rollback(readback,bad,final,final_before)
        require(sha(parent/'State.bin')==before and sha(final/'State.bin')==final_before,'transaction rewrote parent bytes')
        run.result.update(status='passed_transit_transaction',actual_final_time=ordinary['time'],
            split_bitwise=True,governing_stage_hashes_bitwise=True,post_assignment_rollback=True,
            parent_state_history_reference_defect_boundary_preserved=True,accepted_step_executions=8,deliberate_rejections=1,
            parent_checkpoint=str(parent),parent_checkpoint_sha256=before,rollback_parent=str(final),rollback_parent_sha256=final_before)
    except Exception as e:run.result.update(status='terminal_stop',failure=str(e));raise
    finally:run.save()


if __name__=='__main__':main()
