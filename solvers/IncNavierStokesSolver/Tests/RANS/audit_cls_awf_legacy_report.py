#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only audit of an older terminal report with periodic checkpoint lineage.

This does not fabricate the missing final checkpoint or offer a restart path.
The exact-checkpoint spatial-readback contract remains unchanged.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import shutil
from diagnose_cls_accepted_covariance import common_packet,run_probe,independent_witnesses,STAGES
from run_cls_epawf_picard import load,sha,state_hash,require
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',required=True,type=Path)
    p.add_argument('--report',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--contract-tests',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    native=load(a.report);attempt=native['pimple_attempts'][-1];source=Path(attempt['k_snapshot'])
    require(sha(source)==attempt['k_snapshot_sha256'],'native terminal packet differs')
    shutil.copy2(a.report,a.output/'native-report.json');shutil.copy2(source,a.output/'archived.json')
    cp=a.output/'parent-checkpoint';cp.mkdir();original_cp=Path(native['checkpoint'])/'State.bin'
    shutil.copy2(original_cp,cp/'State.bin');q=load(a.output/'archived.json')
    require(state_hash(q['committed'])==state_hash(native['state']),'packet/report committed state differs')
    audit=native['physical_range_events'][-1]['audit'];witnesses=[]
    for name in ('accepted_awf_inner_minimum_point','accepted_awf_outer_minimum_point'):
        if name in audit:witnesses.append(audit[name][:3])
    packet=dict(common_packet(cp),operation='archived_report_stages',native_report=str(a.output/'native-report.json'),
        native_report_sha256=sha(a.output/'native-report.json'),archived_packet=str(a.output/'archived.json'),
        archived_packet_sha256=sha(a.output/'archived.json'),stage_hashes={n:state_hash(q[n]) for n in STAGES},
        native_history_hashes_json=json.dumps([state_hash(v) for v in native['native_history']],separators=(',',':')),
        native_history_times_sha256=state_hash(native['native_history_times']),inner_stress_model=native['inner_stress_model'],witnesses=witnesses)
    d,rc=run_probe(a.probe,a.report.parent/'case.xml',packet,a.output/'stages')
    checks=independent_witnesses(d,q)
    result=dict(status='passed',operation=d['operation'],audit_revision=d['audit_revision'],receipt=rc,
        native_report=str(a.report),native_report_sha256=sha(a.report),archived_packet=str(source),archived_packet_sha256=sha(source),
        periodic_parent_checkpoint=str(original_cp.parent),periodic_parent_checkpoint_sha256=sha(original_cp),
        parent_step=d['parent_checkpoint_step'],last_accepted_step=d['last_native_accepted_step'],
        stage_hashes=packet['stage_hashes'],complete_audits=d['complete_audits'],independent_tensor_checks=checks,
        governing_solves=d['governing_solves'],physical_commits=d['physical_commits'],native_history_advanced=d['native_history_advanced'],
        contract_rejections={})
    if a.contract_tests:
        for name,change in [('report_hash',dict(native_report_sha256='0'*64)),('history',dict(native_history_times_sha256='0'*64)),
                            ('signature',dict(native_signature='incorrect'))]:
            bad=deepcopy(packet);bad.update(change);reason,receipt=run_probe(a.probe,a.report.parent/'case.xml',bad,a.output/('reject-'+name),False)
            result['contract_rejections'][name]=dict(reason=reason[-1800:],receipt=receipt)
    require(sha(a.report)==result['native_report_sha256'] and sha(source)==result['archived_packet_sha256'] and
        sha(original_cp)==result['periodic_parent_checkpoint_sha256'],'readback modified old evidence')
    atomic(a.output/'summary.json',result);print(a.output/'summary.json')

if __name__=='__main__':main()
