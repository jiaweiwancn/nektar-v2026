#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish read-only WP2 ledgers after the two explicitly bounded native legs."""
import argparse,os,subprocess,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,THREADS,sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root;here=Path(__file__).parent
    record=r/'final-readonly-completion.json';require(not record.exists(),'final assessment already has an owner')
    result=dict(version='WP2FinalReadonlyCompletionV1',status='waiting_for_bounded_native_legs',models={},commands=[])
    def save():atomic(record,result)
    def remaining():
        seconds=43200-(time.time()-load(r/'origin.json')['started_epoch']);require(seconds>0,'original WP2 wall envelope exhausted');return seconds
    save()
    try:
        for model in ['CLS','standard']:
            owner=r/f'recovery-{model}-followup-T4p2/progress.json'
            while True:
                remaining();d=load(owner)
                if d['status'] in ['completed_bounded_followup','terminal_stop']:break
                time.sleep(15)
            result['models'][model]=dict(owner=str(owner),owner_sha256=sha(owner),status=d['status'],failure=d.get('failure'))
            if d['status']=='completed_bounded_followup':
                snapshot=Path(d['completed_trajectory_snapshot']['path'])
                require(sha(snapshot)==d['completed_trajectory_snapshot']['sha256'],'completed trajectory snapshot changed')
                cmd=['python3',str(here/'audit_wp2_window_envelope.py'),'--root',str(r),'--output',str(r/f'recovery-{model}-final-all-pairs'),
                     '--models',model,'--phases','recovery','--trajectory',str(snapshot)]
                result['commands'].append(cmd);result['status']='read_only_window_'+model;save()
                with (r/f'recovery-{model}-final-all-pairs.log').open('x') as log:
                    q=subprocess.run(cmd,env=dict(os.environ,**THREADS),stdout=log,stderr=subprocess.STDOUT,timeout=remaining())
                require(q.returncode==0,'read-only final window failed; no later automatic invocation')
            save()
        cmd=['python3',str(here/'collect_wp2_closure.py'),'--root',str(r),'--checkpoints',str(r/'checkpoint-evidence'),'--output',str(r/'final-status.json')]
        result['commands'].append(cmd);result['status']='collecting_closure_ledger';save()
        with (r/'final-closure-collector.log').open('x') as log:q=subprocess.run(cmd,env=dict(os.environ,**THREADS),stdout=log,stderr=subprocess.STDOUT,timeout=remaining())
        require(q.returncode==0,'closure ledger failed')
        result.update(status='completed_read_only_assessment',closed=load(r/'final-status.json')['closed'],final_status_sha256=sha(r/'final-status.json'))
    except Exception as error:result.update(status='terminal_stop',failure=str(error));raise
    finally:save()


if __name__=='__main__':main()
