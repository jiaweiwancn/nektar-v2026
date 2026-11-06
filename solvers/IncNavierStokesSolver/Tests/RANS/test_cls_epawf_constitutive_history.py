#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""V2 active restart, transaction rollback, and refusal of legacy fingerprints."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,write_json,require,state_hash,sha,set_parameter,write
from run_cls_epawf_pimple import prepare,check_report
from continue_cls_epawf_pimple import header
from compare_cls_epawf_pimple import same_history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--legacy-checkpoint',type=Path,required=True);a=p.parse_args()
    a.solver=a.solver.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    summary=dict(schema=1,status='running',cases={})
    def case(name,steps,restart=None,hard=0):
        folder=a.output/name
        prepare(folder,.000625,steps,True,restart=restart,checkpoint=True,k_bound='MassProjection',floor_relative=1e-8,
            snapshots=True,hard_reject_step=hard,inventory='ConstitutiveV2')
        tree=ET.parse(folder/'case.xml').getroot();set_parameter(tree,'IO_CheckSteps',1);write(tree,folder/'case.xml')
        receipt=run([str(a.solver),'case.xml'],folder,'solver.log',14400)
        write_json(folder/'receipt.json',receipt)
        r=load(folder/'case.evolution.json')
        if r.get('initialization_complete') is not False:check_report(r,bool(restart))
        summary['cases'][name]=dict(exit_code=receipt['exit_code'],failure=r.get('failure'),steps=r.get('steps'),time=r.get('time'))
        write_json(a.output/'summary.json',summary)
        return r
    full=case('uninterrupted',6)
    require(full['steps']==6 and not full['failure'],'short V2 transaction trajectory did not complete')
    cps={header(p.parent)['steps']:p.parent for p in (a.output/'uninterrupted').glob('checkpoint*/State.bin')}
    require(4 in cps and 5 in cps and 6 in cps,'missing actual native publications')
    resumed=case('restart4',2,cps[4]);same_history(full,resumed)
    require(any(v['bound_activations'] for v in resumed['pimple_attempts']),'restart interval did not activate correction')
    rejected=case('reject6',2,cps[4],6)
    require(rejected['steps']==5 and 'hard-audit rejection' in rejected['failure'],'deliberate rejection did not occur after active correction')
    cp=header(cps[5],True)
    require(state_hash(rejected['state'])==state_hash(cp['coupled_state']),'rollback state differs')
    require([state_hash(v) for v in rejected['native_history']]==[state_hash(v) for v in cp['history']],'rollback native history differs')
    require(rejected['native_history_times']==cp['history_times'] and rejected['time']==cp['time'],'rollback times differ')
    candidate=load(Path(full['pimple_attempts'][-1]['k_snapshot']))
    require(state_hash(rejected['rejected_state'])==state_hash(candidate['corrected']),'rejection saved a previous state')
    require(rejected['pimple_attempts'][-1]['bound_activations']==1,'rejected step did not contain active bounding')
    legacy_hash=sha(a.legacy_checkpoint/'State.bin')
    mismatch=case('legacy-refused',0,a.legacy_checkpoint)
    require('checkpoint signature mismatch' in mismatch['failure'],'legacy fingerprint was not refused')
    require(legacy_hash==sha(a.legacy_checkpoint/'State.bin'),'legacy checkpoint changed')
    summary.update(status='passed',bitwise_restart_state_time_all_histories=True,
        activated_rejection_rollback=True,actual_rejected_candidate_saved=True,legacy_checkpoint_refused_unchanged=True)
    write_json(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
