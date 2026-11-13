#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh CT5 regression and actual same-method CT6 restart/rollback controls."""
import argparse,json,os,resource,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
from ct_accepted6_cases import make
from ct_accepted_open import run,require
from ct_accepted_pilot import artifact_hash
from run_cls_epawf_picard import set_parameter,set_info,write
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','ct5','member','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    reference=a.ct5/'runs/A95/0.002/T0.008/native'
    old=json.loads((reference/'case.evolution.json').read_text())
    folder=a.member/'trajectory';mms=json.loads((folder/'case.evolution.json').read_text())
    require(mms['dt']==.002 and mms['time']>=.004 and not mms['failure'],'required actual MMS restart source unavailable')
    manifest=RequestManifest(a.output,dict(schema='CT6FlatAndNativeHistoryControlsV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        old_CT5_report_sha256=sha(reference/'case.evolution.json'),MMS_report_sha256=sha(folder/'case.evolution.json'),source_sha256=sha(Path(__file__))))
    checks={}
    flat=a.output/'flat-physical';flat.mkdir();root=ET.parse(reference/'case.xml').getroot()
    set_parameter(root,'NumSteps',2);set_info(root,'RANSWallRestart','None')
    set_info(root,'RANSWallCheckpoint',str(flat/'checkpoint'));set_info(root,'RANSWallPimpleSnapshotDirectory',str(flat/'pimple-states'))
    write(root,flat/'case.xml');os.environ['NEKTAR_RANS_TIMING_FILE']=str(flat/'timing.jsonl')
    e=run([a.runtime/'IncNavierStokesSolver','case.xml'],flat,'solver.log');f=json.loads((flat/'case.evolution.json').read_text())
    expected={round(t,12):s for t,s in zip(old['audit_times'],old['accepted_states'],strict=True)}
    states=all(s==expected[round(t,12)] for t,s in zip(f['audit_times'],f['accepted_states'],strict=True))
    historical=reference/'checkpoint.1/State.bin'
    history=bool(f['checkpoint']) and sha(Path(f['checkpoint'])/'State.bin')==sha(historical)
    checks['fresh_CT5_A95_two_steps']=dict(status='pass' if e['exit_code']==0 and not f['failure'] and f['time']==.004 and states and history else 'fail',
        execution=e,accepted_states_bitwise_identical=states,history_payload_bitwise_identical=history,final_time=f['time'],
        current_report_sha256=sha(flat/'case.evolution.json'),historical_report_sha256=sha(reference/'case.evolution.json'),
        historical_history_sha256=sha(historical),fresh_history_sha256=sha(Path(f['checkpoint'])/'State.bin') if f['checkpoint'] else None)
    write_json(a.output/'progress.json',checks)
    first=folder/'checkpoint';second=folder/'checkpoint.1'
    for name,reject,amplitude in (('unchanged_restart',0,.005),('deliberate_rejection_rollback',2,.005),('incompatible_geometry_history',0,0.)):
        # The manufactured boundary identity deliberately binds all non-work
        # parameters, including the deliberate rejection fixture. Give that
        # fixture its own fresh history; never edit or relabel a signature.
        sub=a.output/name;make(sub,a.runtime,model='A95',degree=4,dt=.002,
            steps=2 if reject else 1,restart=None if reject else first,amplitude=amplitude)
        if reject:
            root=ET.parse(sub/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',reject);write(root,sub/'case.xml')
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(sub/'timing.jsonl')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],sub,'solver.log')
        path=sub/'case.evolution.json';r=json.loads(path.read_text()) if path.exists() else {};log=(sub/'solver.log').read_text()
        if name=='incompatible_geometry_history':
            passed=e['exit_code']!=0 and ('signature' in log.lower() or 'physics' in log.lower())
            row=dict(status='pass' if passed else 'fail',intended_rejection='checkpoint physics/signature mismatch',execution=e,log_sha256=sha(sub/'solver.log'),log_tail=log[-1200:])
        else:
            t=.004 if not reject else .002
            expected_state=next(s for tt,s in zip(mms['audit_times'],mms['accepted_states'],strict=True) if abs(tt-t)<1e-13)
            identical=r.get('state')==expected_state
            if reject:
                candidate=next(s for tt,s in zip(mms['audit_times'],mms['accepted_states'],strict=True) if abs(tt-.004)<1e-13)
                retained=r.get('rejected_state')==candidate
                passed=e['exit_code']!=0 and identical and retained and abs(r.get('time',-1)-.002)<1e-13 and 'deliberate' in r.get('failure','')
                row=dict(status='pass' if passed else 'fail',execution=e,last_accepted_state_bitwise_identical=identical,
                    rejected_candidate_equals_uninterrupted_next_state=retained,failure=r.get('failure'),
                    initial_state_equals_uninterrupted_initial=r.get('initial')==mms['initial'],
                    time=r.get('time'),native_history_times=r.get('native_history_times'),report_sha256=sha(path))
            else:
                same_history=bool(r.get('checkpoint')) and sha(Path(r['checkpoint'])/'State.bin')==sha(second/'State.bin')
                row=dict(status='pass' if e['exit_code']==0 and not r.get('failure') and identical and same_history else 'fail',
                    execution=e,accepted_state_bitwise_identical=identical,history_payload_bitwise_identical=same_history,
                    time=r.get('time'),report_sha256=sha(path),history_sha256=sha(Path(r['checkpoint'])/'State.bin') if r.get('checkpoint') else None)
        checks[name]=row;write_json(a.output/'progress.json',checks);print('CT6 regression',name,row['status'],flush=True)
    os.environ.pop('NEKTAR_RANS_TIMING_FILE',None)
    result=dict(schema='CT6FlatAndNativeHistoryControlsV1',status='pass' if all(q['status']=='pass' for q in checks.values()) else 'fail',checks=checks)
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
