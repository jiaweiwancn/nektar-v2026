#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only accepted-103 / rejected-104 branch audit. Never advances old tail."""
import argparse
import shutil
import sys
from pathlib import Path
from continue_cls_epawf_pimple import header
from run_cls_epawf_picard import load,write_json,require,run,sha,state_hash
from run_cls_epawf_pimple import identity

CP_SHA='11d1c23373ed98a3f79c7ce6f242c23883398a3e369e1210e34a81a1c42445e9'
REPORT_SHA='32867b868a613d90b0288dabbafc0df5dec65aad7c79d919dd307ca664f0f43d'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--segment',type=Path,required=True);p.add_argument('--boundary-report',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();require(sys.platform=='linux','Linux authoritative')
    a.output.mkdir(parents=True,exist_ok=False)
    report=a.archive/'upstream-terminal.evolution.json';cp=a.archive/'upstream-checkpoint'
    require(sha(report)==REPORT_SHA and sha(cp/'State.bin')==CP_SHA,'immutable terminal artifact differs')
    e=load(report);checkpoint=header(cp,True);require(checkpoint['steps']==103 and e['steps']==103,'wrong historical frontier')
    accepted=checkpoint['coupled_state'];rejected=e['rejected_state'];offset=e['offsets'];k=int(offset[2]);end=int(offset[3])
    require(accepted[k:end]==rejected[k:end],'flow changed k coefficients: investigate a different mechanism')
    attempt=e['pimple_attempts'][-1]
    # last_stage_history was transactionally rolled back. Ask the native GLM
    # to form the callback from the exact saved slots, then abort before any
    # solve and verify the archived callback hash (including signed zeros).
    shutil.copytree(cp,a.output/'checkpoint103');shutil.copy2(a.segment/'case.xml',a.output/'case.xml')
    boundary=load(a.boundary_report)['fixed_external_coefficients']
    request=dict(operation='branch_coverage',state=accepted,compare_state=rejected,committed=accepted,
        form_native_callback_history=True,expected_history_sha256=attempt['history_sha256'],time=attempt['time'],alpha=attempt['alpha'],
        native_checkpoint=str((a.output/'checkpoint103').resolve()),checkpoint_sha256=CP_SHA,native_signature=e['signature'],
        fixed_external_coefficients=boundary,projection_only=False,require_history_h0=False,
        inner_stress_model='SaturatedR2',stations=[-3.5,-2.625,-1.,1.6])
    write_json(a.output/'input.json',request)
    receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],a.output,'probe.log',1800)
    manifest=dict(status='read_only_replay',checkpoint_sha256=CP_SHA,terminal_report_sha256=REPORT_SHA,
        report_path=str(report),boundary_path=str(a.boundary_report),boundary_sha256=sha(a.boundary_report),
        accepted_state_sha256=state_hash(accepted),rejected_state_sha256=state_hash(rejected),
        history_sha256=attempt['history_sha256'],k_coefficient_difference=0,probe=identity(a.probe.resolve()),
        receipt=receipt,governing_solves=0,old_tail_evolved=False,checkpoint_migration=False)
    write_json(a.output/'manifest.json',manifest);require(receipt['exit_code']==0,'read-only replay failed')
    r=load(a.output/'result.json');require(r['k_coefficient_difference']==0,'replayed flow changed k')
    before=r['state']['anchors'][0];after=r['compare']['anchors'][0]
    manifest.update(before=before,after=after,result_sha256=sha(a.output/'result.json'))
    require(before['chi']==0 and after['chi']>0,'archived gate switch was not reproduced')
    branches={v['branch']:v for v in before['profiles']}
    require(branches['current']['minimum']>=0 and branches['prospective_active']['minimum']<0,
        'prospective branch did not expose the saved old-acceptance gap')
    write_json(a.output/'manifest.json',manifest);print(a.output/'manifest.json')

if __name__=='__main__':main()
