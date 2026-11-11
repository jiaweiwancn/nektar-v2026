#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh native open histories, pressure offsets, rollback and closed controls."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import fixture,probe,sampling_request,run,require
from run_cls_epawf_pimple import prepare
from run_cls_epawf_picard import set_info,write
from compare_cls_epawf_pimple import same_history
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--old-runtime',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    manifest=RequestManifest(a.output,dict(schema='CTAcceptedOpenNativeControlsV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        old_runtime=runtime_identity(a.old_runtime/'IncNavierStokesSolver'),dt=1e-7,steps=2,pressure_offset=.125,
        pressure_offset_limit=2e-10,restart='bitwise state and history',closed='bitwise state/history at unchanged identities'))
    result=dict(status='fail',checks={},executions=[])
    def execute(name,steps,**options):
        folder=a.output/name;fixture(folder,a.runtime,steps=steps,**options)
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(dict(name=name,**e))
        path=folder/'case.evolution.json';r=json.loads(path.read_text()) if path.exists() else None
        write_json(a.output/'progress.json',result)
        return r,e,(folder/'solver.log').read_text()
    try:
        full,e,log=execute('full',2);require(e['exit_code']==0 and not full['failure'],'open_full_execution:'+log[-1500:])
        first,e,log=execute('first',1);require(e['exit_code']==0,'open_first_execution')
        continued,e,log=execute('continued',1,restart=Path(first['checkpoint']));require(e['exit_code']==0 and same_history(full,continued),'same_method_restart_mismatch')
        result['checks']['open_restart']=dict(status='pass',state_sha256=exact_vector_hash(full['state']),history=full['native_history_times'])
        shifted,e,log=execute('shifted',2,pressure=.125);require(e['exit_code']==0,'pressure_offset_execution')
        state=np.array(full['state']);other=np.array(shifted['state']);offset=int(full['differential_size']) if 'differential_size' in full else len(state)
        # Obtain the actual field offsets, constants and pressure layout.
        packet,e=probe(a.runtime,a.output/'full',dict(operation='functions',grid=[],export_basis=True),'initial-layout')
        o=packet['observation'];off=np.array(o['offsets'],int);expected=np.zeros(len(state))
        for b in (4,5):
            for q in o['pressure_mass']:
                if int(q['region'])==b-4:
                    ids=np.array(q['pressure']['ids'],int)+off[b];expected[ids]=.125*np.array(q['pressure']['constant'])
        defect=float(np.max(abs(other-state-expected)));require(defect<=2e-10,'pressure_offset_evolved_velocity_scalar_or_level')
        result['checks']['pressure_offset_evolution']=dict(status='pass',maximum_coefficient_defect=defect,offset=.125)
        rejected,e,log=execute('rejected',2,reject=2)
        require(e['exit_code']!=0 and 'deliberate PIMPLE rejection' in rejected['failure'] and same_history(first,rejected),'rejected_step_did_not_fully_rollback')
        result['checks']['open_rollback']=dict(status='pass',reason=rejected['failure'],last_accepted_time=rejected['time'])
        wrong,e,log=execute('changed_traction_restart',0,restart=Path(first['checkpoint']),traction_shift=.01)
        require(e['exit_code']!=0 and ('signature' in log.lower() or 'fingerprint' in log.lower()),'changed_boundary_restart_wrong_rejection')
        result['checks']['changed_boundary_restart']=dict(status='pass',expected_rejection='physics signature',witness=log[-2000:])
        folder=a.output/'internal_reverse';fixture(folder,a.runtime,internal_reverse=True)
        request=dict(operation='functions',grid=[[.5,0,1.],[1.,0,1.]],grid_sides=[0,0],export_basis=False)
        packet,e=probe(a.runtime,folder,request);require(packet['open_boundary']['outlet_status']['admissible'] and
            packet['observation']['grid_values'][0][0]<0 and packet['observation']['grid_values'][1][0]>0,'internal_reversal_misclassified_as_outlet_backflow')
        result['checks']['internal_reverse']=dict(status='pass',native_physical_status=packet['status'],outlet=packet['open_boundary']['outlet_status'],velocities=packet['observation']['grid_values'])
        reports={}
        for name,runtime,cache in [('closed_old',a.old_runtime,'Off'),('closed_new',a.runtime,'Off'),('closed_cached',a.runtime,'ExactStepV1')]:
            os.environ['LD_LIBRARY_PATH']=str(runtime/'lib');folder=a.output/name
            xml=prepare(folder,1e-4,2,diagnostic=True,checkpoint=True);root=ET.parse(xml).getroot();set_info(root,'RANSWallEvaluationReuse',cache);write(root,xml)
            e=run([runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(dict(name=name,**e));require(e['exit_code']==0,'closed_CLS_execution:'+name)
            reports[name]=json.loads((folder/'case.evolution.json').read_text())
        require(same_history(reports['closed_old'],reports['closed_new']) and same_history(reports['closed_new'],reports['closed_cached']),'closed_CLS_trajectory_history_changed')
        result['checks']['closed_CLS_compatibility_cache']=dict(status='pass',state_sha256=exact_vector_hash(reports['closed_new']['state']),scope='two ordinary CLS diagnostic steps; no renewed physical qualification')
        os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
        e=run([a.runtime/'RANSNativeDiagnosticChecks'],a.output,'native-diagnostics.log');result['executions'].append(e);require(e['exit_code']==0,'existing_native_diagnostic_negative_controls')
        result['checks']['existing_diagnostic_controls']=dict(status='pass')
        result['status']='pass'
    except Exception as error:result['failure']=repr(error)
    write_json(a.output/'results.json',result)
    if result['status']=='pass':manifest.seal(a.output,['results.json']+[str(f.relative_to(a.output)) for f in a.output.rglob('case.xml')]+[str(f.relative_to(a.output)) for f in a.output.rglob('solver.log')])
    print(result['status'],result.get('failure',''));return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
