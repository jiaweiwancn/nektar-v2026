#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One bounded local time-sensitivity control; no global recovery promotion."""
import argparse,time
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha,state_hash
from diagnose_cls_accepted_covariance import checkpoint
from continue_cls_epawf_pimple import atomic
from wp2_native import Run
from run_cao_tamura_budget_wp2 import charge
from wp2_budget import wall_limit


def main():
    parser=argparse.ArgumentParser()
    for name in ['root','solver','h1','preflight','output']:parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args();parent=a.root/'preservation-recovery-standard/step336/checkpoint.32'
    seed=checkpoint(parent);before=sha(parent/'State.bin');out=Run(a.root,a.solver,a.h1,a.preflight,a.output,'standard',4,16,'RadialEnergyV1')
    started=time.monotonic();out.result.update(status='validating_local_import',parent_checkpoint=str(parent),parent_sha256=before,
        physical_t0=seed['time'],maximum_fine_steps=68,dt=.005,cumulative_wall_limit=wall_limit(a.root));out.save()
    try:
        charge(a.root,'late standard local dt control',timestep_trials=1,boundary_calls=4)
        ordinary=out.native('ordinary-one',1,.01,parent=parent)
        imported=out.native('imported-one',1,.01,ivp_parent=parent)
        archived=load(a.root/'preservation-recovery-standard/step368/case.evolution.json')
        checks={}
        for name in ['history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256']:
            checks[name]=ordinary['pimple_attempts'][0][name]==imported['pimple_attempts'][0][name]==archived['pimple_attempts'][0][name]
        for name in ['state','fixed_reference','fixed_defect','native_history','native_history_times']:
            if name=='native_history':same=[state_hash(x) for x in ordinary[name]]==[state_hash(x) for x in imported[name]]
            else:same=state_hash(ordinary[name])==state_hash(imported[name])
            checks[name]=same
        checks.update(archived_first_state=state_hash(ordinary['state'])==state_hash(archived['accepted_states'][1]),
            exact_seed=state_hash(imported['accepted_states'][0])==state_hash(seed['coupled']),
            no_initialization_Lift=imported['initialization_lift_calls_iterations']==[0,0],
            actual_time=ordinary['time']==imported['time'],distinct_counters=ordinary['steps']==337 and imported['steps']==1,
            candidate_differs_from_reference=state_hash(imported['accepted_states'][0])!=state_hash(imported['fixed_reference']))
        require(all(checks.values()),'late analytical-open same-dt import differs: '+str(checks))
        out.result['import_validation']=dict(passed=True,checks=checks);out.save()
        out.native('ordinary-changed-dt',0,.005,parent=parent,mismatch=True)
        initial=out.native('fine-initial',0,.005,ivp_parent=parent)
        require(state_hash(initial['state'])==state_hash(seed['coupled']) and initial['initialization_lift_calls_iterations']==[0,0],
                'fine physical seed changed')
        require(initial['pimple_attempts']==[],'initialization advanced a physical step')
        cp=Path(initial['checkpoint']);out.result.update(status='running_fine_local_control',accepted_fine_steps=0);out.save()
        for end in [32,64,68]:
            start=out.result['accepted_fine_steps'];remaining=wall_limit(a.root)-(time.time()-load(a.root/'origin.json')['started_epoch'])
            require(remaining>(end-start)*27+90,'remaining cumulative budget cannot support next native block')
            report=out.native('fine-step'+str(end),end-start,.005,parent=cp);cp=Path(report['checkpoint'])
            out.result.update(accepted_fine_steps=end,actual_time=report['time'],latest_checkpoint=str(cp),latest_checkpoint_sha256=sha(cp/'State.bin'));out.save()
        out.result.update(status='completed_local_control',local_elapsed_time=report['time']-seed['time'],
            interpretation='Common-state local sensitivity only. Original global perturbation recovery and WP2 closure remain OPEN.')
    except Exception as error:
        out.result.update(status='terminal_stop',failure=str(error));raise
    finally:
        out.result.update(parent_preserved=sha(parent/'State.bin')==before,driver_seconds=time.monotonic()-started);out.save()


if __name__=='__main__':main()
