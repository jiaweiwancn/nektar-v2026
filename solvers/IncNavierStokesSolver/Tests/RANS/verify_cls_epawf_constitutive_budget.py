#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent frozen accounting for immutable V2 projection-only repairs."""
import argparse
import math
import shutil
from pathlib import Path
import numpy as np
from run_cls_epawf_picard import load,require,write_json,run,state_hash


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--legacy',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);a=p.parse_args()
    mass_source=load(a.legacy/'analysis-final/contact0025/result.json')
    mass=np.array(mass_source['free_k_mass']);free=np.array(mass_source['free_k_rows'],dtype=int)
    result=dict(schema=1,production_governing_solves=0,diagnostic_transport_solves=0,cases={})
    for name,report_path in [('primary','main/steps-000006-000021/case.evolution.json'),('finer','finer-from-h0/case.evolution.json')]:
        report=load(a.legacy/report_path);packet=load(Path(report['pimple_attempts'][-1]['k_snapshot']))
        repaired=load(a.fixtures/(name+'-v2/result.json'));require(repaired['k_bound']['success'],'missing successful repair')
        audit_folder=a.output.parent/(name+'-repaired-audit');audit_folder.mkdir()
        shutil.copy2(a.fixtures/(name+'-v2/case.xml'),audit_folder/'case.xml')
        request=load(a.fixtures/(name+'-v2/input.json'));request.update(operation='audit_witness',state=repaired['corrected'])
        write_json(audit_folder/'input.json',request)
        receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],audit_folder,'probe.log',3600)
        require(receipt['exit_code']==0,'hard-acceptance replay failed');audit_result=load(audit_folder/'result.json')
        require(not audit_result['failure'],'repaired audit raised a constitutive failure')
        audit=audit_result['audit'];K0=repaired['k_bound']['K0']
        require(audit['admissible']==[1] and min(audit['r2_covariance'])>=-1e-10*K0 and audit['r2_scalar_domain'][8]>=.5,
                'repaired fixture failed hard scalar/covariance/normalization acceptance')
        if not repaired['pre_k_inventory']['admissible']:
            guard=a.output.parent/(name+'-production-entry-refused');guard.mkdir()
            shutil.copy2(a.fixtures/(name+'-v2/case.xml'),guard/'case.xml')
            request=load(a.fixtures/(name+'-v2/input.json'));request['projection_only']=False
            write_json(guard/'input.json',request)
            receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],guard,'probe.log',3600)
            require(receipt['exit_code']==0,'guard diagnostic failed');refused=load(guard/'result.json')
            require(not refused['k_bound']['success'] and 'pre-k state is inadmissible' in refused['k_bound']['failure'],
                    'invalid legacy pre-k state bypassed production guard')
            require(state_hash(refused['raw'])==state_hash(refused['corrected']),'refused entry modified the candidate')

        rows=np.array(packet['transport_rows'],dtype=int);require(np.array_equal(rows,free),'native free mass order differs')
        n=len(rows);matrix=np.frombuffer(Path(packet['matrix_file']).read_bytes(),dtype='<f8').reshape(n,n)
        raw_delta=np.array(packet['raw_increment']);raw_res=np.array(packet['transport_residual'])+matrix@raw_delta
        correction=(np.array(repaired['corrected'])-np.array(repaired['raw']))[rows]
        added=matrix@correction;storage=packet['alpha']*(mass@correction);corrected=raw_res+added
        mass_norm=math.sqrt(correction@(mass@correction))
        require(abs(mass_norm-repaired['k_bound']['mass_norm_correction'])<1e-12,'independent correction metric differs')
        require(np.max(np.abs(raw_res))<1e-12,'raw transport residual not roundoff')
        require(np.max(np.abs(corrected-added))<1e-12,'corrected residual identity differs')
        raw_full=np.array(packet['raw']);new_full=np.array(repaired['corrected']);before=np.array(packet['pre_k'])
        net=(new_full-before)[rows];raw=(raw_full-before)[rows]
        require(np.max(np.abs(net-raw-correction))<1e-15,'increment cancellation identity differs')
        offsets=list(map(int,repaired['offsets']));trace=np.array(repaired['scalar_domain']['trace_target_state'])
        trace_defect=float(np.max(np.abs(trace[offsets[3]:offsets[4]]-new_full[offsets[3]:offsets[4]])))
        result['cases'][name]=dict(hard_acceptance_passed=True,constitutive_qualified=audit['r2_qualified']==[1],
            complete_covariance=audit['r2_covariance'],lambda_min=audit['r2_scalar_domain'][8],
            legacy_pre_k_admissible_under_v2=repaired['pre_k_inventory']['admissible'],
            raw_frozen_residual=float(np.max(np.abs(raw_res))),
            corrected_frozen_residual=float(np.max(np.abs(corrected))),added_A_delta=float(np.max(np.abs(added))),
            storage_alpha_M_delta=float(np.max(np.abs(storage))),mass_norm=mass_norm,
            corrected_step_mass_norm=math.sqrt(net@(mass@net)),raw_step_mass_norm=math.sqrt(raw@(mass@raw)),
            retained_epsilon_matching_coefficient_defect=trace_defect,
            A_delta_minus_storage_max=float(np.max(np.abs(added-storage))),
            cancellation_defect=float(np.max(np.abs(net-raw-correction))))
        write_json(a.output.with_name(a.output.stem+'-'+name+'-vectors.json'),dict(A_delta=added.tolist(),alpha_M_delta=storage.tolist()))
    result['status']='passed';write_json(a.output,result);print(a.output)


if __name__=='__main__':main()
