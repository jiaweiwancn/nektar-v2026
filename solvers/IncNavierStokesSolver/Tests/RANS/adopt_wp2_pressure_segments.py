#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Revalidate complete native segments after a timing-wrapper interruption.

The original unknown process exit status is retained. Only a complete native
report, all accepted histories/checkpoints, and a successful exact zero-step
restart can supply a successor prefix. No physical step is replayed here.
"""
import argparse,concurrent.futures,json,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import require,sha,state_hash
from wp2_native import Run
from wp2_qualification import verify_report
from run_wp2_preservation import US,K0,TREQ


def adopt(a,phase,model):
    root=a.root;old_path=root/('pressure-'+phase+'-'+model)/'progress.json';old=load(old_path)
    require(old['solver_sha256']==sha(a.solver) and old['case_definition']==dict(level=4,order=16,projection='RadialEnergyV1',pressure_space='inner-p2-v1'),'interrupted segment runtime/case differs')
    require(old['status']=='terminal_stop' and 'timed out' in old['failure'],'not the identified wrapper interruption')
    file=old_path.parent/'step30/case.evolution.json';report=load(file);report_sha=sha(file)
    require(not report['failure'] and report['implicit_callback_count']==30 and report['steps']==30 and
            abs(report['time']-.3)<2e-14 and report['dt']==.01,'native segment did not finish its declared endpoint')
    require(len(report['accepted_states'])==len(report['state_audits'])==len(report['wp2_per_state_monitoring'])==31 and
            len(report['pimple_attempts'])==30,'incomplete native accepted-state inventory')
    for row in report['pimple_attempts']:
        for key in ['outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls']:
            require(row[key]==1,'completed segment governing count '+key)
        require(row['solved_epsilon_retained_bitwise'],'completed segment reset solved epsilon')
    checks=verify_report(report,US,K0[model],allow_recovery=phase=='recovery')
    parent=Path(report['checkpoint']);before=sha(parent/'State.bin');cp=checkpoint(parent)
    require(cp['signature']==report['signature'] and cp['steps']==30 and cp['time']==report['time'] and cp['dt']==.01,'parent metadata differs')
    require(state_hash(cp['coupled'])==state_hash(report['state']),'parent accepted coefficients differ')
    for key in ['reference','defect']:require(state_hash(cp[key])==state_hash(report['fixed_'+key]),'parent fixed data differs')
    timing=[json.loads(line) for line in (file.parent/'timing.jsonl').read_text().splitlines()]
    total=[row for row in timing if row['name']=='process'];require(len(total)==1,'original native process timing is incomplete')
    require(not Path('/proc/'+str(total[0]['pid'])).exists(),'original native PID is still present; no duplicate launch')
    output=root/('pressure-adopted-prefix-'+phase+'-'+model)
    if output.exists():
        saved=load(output/'progress.json');require(saved['status']=='completed_bounded_trajectory' and
            saved['solver_sha256']==sha(a.solver) and saved['results']['completed_original_segment']['report_sha256']==report_sha,
            'existing prefix is not the already validated exact native segment')
        item=saved['results']['exact-checkpoint-revalidation'];require(item['code']==0 and sha(Path(item['report']))==item['report_sha256'],'saved exact restart changed')
        archive=output/'before-interface-scope-correction.json';require(not archive.exists(),'prefix already reassessed');atomic(archive,saved)
        x=Run.__new__(Run);x.result=saved;x.output=output;native=load(Path(item['report']))
    else:
        x=Run(root,a.solver,a.h1,a.preflight,output,model,4,16,'RadialEnergyV1',False,'inner-p2-v1')
        native=x.native('exact-checkpoint-revalidation',0,.01,parent=parent)
    require(same_history(report,native),'zero-step restart changed accepted state/history/time')
    for field in ['fixed_reference','fixed_defect']:require(state_hash(report[field])==state_hash(native[field]),'zero-step restart changed '+field)
    require(sha(parent/'State.bin')==before and sha(file)==report_sha,'revalidation modified its immutable parent data')
    x.result['results']['completed_original_segment']=dict(report=str(file),report_sha256=report_sha,
        code=None,native_process_exit_status='not observed after timing-wrapper SIGKILL; never relabelled zero',
        observed_physical_attempts=30,accepted_step_executions=30,requested_attempt_work=30,
        native_process_timing=total[0],original_driver_progress=str(old_path),original_driver_progress_sha256=sha(old_path),
        original_driver_failure=old['failure'],complete_native_report_independently_checked=True)
    x.result.update(status='completed_bounded_trajectory',mode='recover' if phase=='recovery' else 'preserve',dt=.01,Treq=TREQ,
        accepted_steps=30,actual_time=report['time'],physical_endpoint=.3,reference_defect_preserved=True,
        final_physical_budgets_passed=checks[-1]['passed'],resource_recovery=True,
        qualification='complete native accepted segment, unknown original process return code, successful exact zero-step native restart; original failed driver remains immutable',
        observations=[dict(time=report['time'],accepted=30,checkpoint=str(parent),checkpoint_sha256=before,report=str(file),
            physical_budgets_passed=all(v['passed'] for v in checks),final_physical_checks=checks[-1],seed_sha256=state_hash(report['accepted_states'][0]))])
    x.save();return dict(path=str(output/'progress.json'),sha256=sha(output/'progress.json'),original_report_sha256=report_sha,
                         original_exit_status_unknown=True,exact_restart_return_code=0,physical_steps_added=0)


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','h1','preflight']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();record=a.root/'pressure-adoption-v2.json';require(not record.exists(),'adoption already has an owner')
    result=dict(version='WP2CompleteSegmentRevalidationV2',status='validating',prefixes={},started_epoch=time.time())
    atomic(record,result)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            jobs={(phase,model):pool.submit(adopt,a,phase,model) for phase in ['recovery','fine'] for model in ['CLS','standard']}
            for (phase,model),job in jobs.items():
                result['prefixes'].setdefault(model,{})[phase]=job.result();atomic(record,result)
        result.update(status='passed_complete_segment_revalidation',source_sha256=sha(Path(__file__)),solver_sha256=sha(a.solver),
            no_physical_replay=True,original_failures_preserved=True)
    except Exception as e:result.update(status='terminal_stop',failure=str(e));raise
    finally:atomic(record,result)


if __name__=='__main__':main()
