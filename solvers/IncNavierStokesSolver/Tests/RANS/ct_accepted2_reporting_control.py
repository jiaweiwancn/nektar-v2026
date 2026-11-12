#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reporting-only rebuild: same frozen inputs, fields, failure and native history."""
import argparse,json,os,copy
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import require,run,probe
from ct_accepted2_diagnosis import make_case,WITNESS
from ct_accepted_feasibility import write_json
from compare_cls_epawf_pimple import same_history
from rans_request_manifest import RequestManifest,runtime_identity,sha
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();old=r/'runtime';new=r/'reporting/runtime';os.environ['LD_LIBRARY_PATH']=str(old/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT2ProjectionMetadataOnlyV1',old_runtime=runtime_identity(old/'IncNavierStokesSolver'),new_runtime=runtime_identity(new/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__)),arithmetic='no numerical changes; incorrect open Functions projection description corrected'))
    source=r/'matrix-v1/work/B-0.002';c=json.loads((source/'construction.json').read_text());case=a.output/'case';make_case(case,source/'case.xml',c,'B',.002,2)
    packets=[]
    for label,rt in [('old',old),('reporting',new)]:
        os.environ['LD_LIBRARY_PATH']=str(rt/'lib');packet,e=probe(rt,case,dict(operation='functions',grid=[WITNESS],grid_sides=[0],export_basis=False),'initial-'+label);require(e['exit_code']==0,'reporting_probe_execution');packets.append(packet)
    x,y=[copy.deepcopy(q['observation']) for q in packets];before=x['actual_quadrature'].pop('reference_projection');after=y['actual_quadrature'].pop('reference_projection')
    require(x==y and 'Gauss32' in before and 'Functions:' in after and 'max(8,p+3)' in after,'reporting_metadata_or_numeric_observation_mismatch')
    e=run([new/'IncNavierStokesSolver','case.xml'],case,'reporting-solver.log');actual=json.loads((case/'case.evolution.json').read_text());baseline=json.loads((source/'case.evolution.json').read_text())
    require(e['exit_code']!=0 and same_history(actual,baseline) and actual['rejected_state']==baseline['rejected_state'] and actual['failure']==baseline['failure'],'reporting_only_rebuild_changed_trajectory')
    result=dict(status='pass',numerical_observation='bitwise identical apart from corrected projection description',projection_before=before,projection_after=after,accepted_rejected_states_native_history_and_failure='bitwise identical',governing_cohorts_runtime='original CT2 runtime only; reporting rebuild not substituted into sealed trajectories',execution=e,report_sha256=sha(case/'case.evolution.json'))
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print('metadata regression pass');return 0
if __name__=='__main__':raise SystemExit(main())
