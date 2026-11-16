#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Explicit same-checkpoint recovery milestone, without resetting WP2 budgets.

The optional transaction wait reserves the spare native slot for sequential
terminal probes and one continuation. It never retries a rejected physical leg.
"""
import argparse,os,subprocess,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,THREADS,sha


def main():
    p=argparse.ArgumentParser()
    for key in ['root','solver','probe','h1','preflight','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    p.add_argument('--endpoint',type=float,required=True)
    p.add_argument('--after-transactions',action='store_true')
    a=p.parse_args();a.output.mkdir(exist_ok=False);here=Path(__file__).parent
    record=dict(version='WP2BoundedFollowupV1',status='starting',model=a.model,endpoint=a.endpoint,
                solver_sha256=sha(a.solver),commands=[],origin_sha256=sha(a.root/'origin.json'))
    def save():atomic(a.output/'progress.json',record)
    def remaining():
        seconds=43200-(time.time()-load(a.root/'origin.json')['started_epoch'])
        require(seconds>0,'original WP2 wall envelope exhausted');return seconds
    def execute(name,command):
        record['status']=name;record['commands'].append(command);save()
        with (a.output/(name+'.log')).open('x') as log:
            result=subprocess.run(command,env=dict(os.environ,**THREADS),stdout=log,stderr=subprocess.STDOUT,timeout=remaining())
        require(result.returncode==0,'followup stopped after '+name+'; no later invocation permitted')
    common=['--root',str(a.root),'--solver',str(a.solver),'--h1',str(a.h1),'--preflight',str(a.preflight)]
    save()
    try:
        if a.after_transactions:
            record['status']='waiting_for_existing_transit_transactions';save()
            while True:
                remaining();statuses=[load(a.root/f'transit-transaction-{m}/progress.json')['status'] for m in ['CLS','standard']]
                require('terminal_stop' not in statuses,'existing transit transaction failed')
                if all(s=='passed_transit_transaction' for s in statuses):break
                time.sleep(15)
            execute('terminal_seen_grid_checks',['python3',str(here/'run_wp2_final_checks.py'),*common,
                    '--probe',str(a.probe),'--terminal-only','--record','terminal-checks-after-transactions.json'])
            require(load(a.root/'terminal-checks-after-transactions.json')['status']=='completed_bounded_final_checks','terminal checks incomplete')
        before=a.root/f'preservation-recovery-{a.model}/progress.json'
        prior=load(before)
        require(prior['status']=='completed_bounded_trajectory' and prior['actual_time']<a.endpoint,'recovery milestone is not a new compatible continuation')
        require(a.endpoint<=4*1.1942387693600194,'declared recovery cap')
        record['parent_progress_sha256']=sha(before);save()
        execute('same_checkpoint_continuation',['python3',str(here/'run_wp2_preservation.py'),*common,
                '--output',str(before.parent),'--model',a.model,'--mode','recover','--level','4','--dt','.01',
                '--endpoint',str(a.endpoint),'--resume'])
        receipt=load(before);require(receipt['status']=='completed_bounded_trajectory','recovery did not reach the bounded endpoint')
        execute('read_only_return_assessment',['python3',str(here/'audit_wp2_recovery_leg.py'),'--root',str(a.root),
                '--output',str(a.output/'assessment'),'--model',a.model])
        result=a.output/'assessment/result.json'
        record.update(status='completed_bounded_followup',actual_time=receipt['actual_time'],
                      assessment_sha256=sha(result),recovery_passed=load(result)['passed'])
    except Exception as error:record.update(status='terminal_stop',failure=str(error));raise
    finally:save()


if __name__=='__main__':main()
