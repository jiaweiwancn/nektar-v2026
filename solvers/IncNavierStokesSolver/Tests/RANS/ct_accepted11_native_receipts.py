#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Expose retained native trace/splitting receipts without re-evaluating a state."""
import argparse,json,os
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha
from ct_accepted_feasibility import require,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--qualified',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();read=lambda f:json.loads(f.read_text())
    seal=read(a.qualified/'completed-artifacts.json')
    require(seal['request_sha256']==read(a.qualified/'request-manifest.json')['sha256'] and
            all(sha(a.qualified/k)==v for k,v in seal['files'].items()),'changed CT11 collection')
    target=read(a.qualified/'target.json');bindings={};members={}
    keys=('time','alpha','imposed_trace_defect','retained_trace_column_max',
          'trace_mass_column_max','trace_physical_storage_contribution_max',
          'refreshed_trace_defect','epsilon_frozen_linear_residual',
          'k_frozen_linear_residual','epsilon_same_state_source_storage_defect',
          'k_same_state_source_storage_defect','refreshed_original_block_norms',
          'corrected_continuity_gauge','frozen_momentum_defect',
          'governing_solves_complete','scheduled_complete','failure',
          'outer_sweeps','bound_calls','additional_governing_equation_solves',
          'additional_correction_or_restoration_solves')
    for dt,member in target['members'].items():
        f=Path(member['directory'])/'case.evolution.json'
        require(sha(f)==member['native_report_sha256'],'changed native report:'+dt)
        report=read(f);rows=[]
        for attempt in report['pimple_attempts']:
            row={k:attempt[k] for k in keys}
            for name in ('flow_snapshot','k_snapshot'):
                path=Path(attempt[name]);digest=attempt[name+'_sha256']
                require(sha(path)==digest,'changed actual native stage:'+str(path))
                row[name]=dict(path=str(path),sha256=digest)
            rows.append(row)
        require(len(rows)==report['steps'] and report['failure']=='','incomplete native accepted receipts')
        bindings[dt]=dict(path=str(f),sha256=member['native_report_sha256'])
        members[dt]=dict(dt=report['dt'],time=report['time'],state_sha256=member['state_sha256'],
            history_sha256=member['history_sha256'],configuration_sha256=member['configuration_sha256'],
            raw_traction=member['traction_history'],stages=rows,
            effective_inputs={k:v['actual_coefficient_stages'] for k,v in member['budgets']['budgets'].items()})
    m=RequestManifest(a.output,dict(schema='CT11RetainedNativeReceiptsV1',
        qualified_seal=sha(a.qualified/'completed-artifacts.json'),native_reports=bindings,
        source=sha(Path(__file__)),scope='direct extraction of sealed native fields; no new solve,state evaluation or changed qualification'))
    write_json(a.output/'results.json',dict(schema='CT11RetainedNativeReceiptsV1',members=members,
        semantics='raw accepted-state traction remains separate from refreshed dependent-epsilon trace lag; neither is repaired or substituted for the other'))
    m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
