#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Recover only reporting from completed, identified native temporal runs.

No executable is launched, no physical state/history is changed. A saved replay
packet must match the native final coefficients, time and initial reference.
"""
import argparse
import math
from pathlib import Path
import sys
from rans_qualification import write_json,exact_vector_hash,audited_extrema,coverage,compare_samples
from qualify_rans_preservation import load,MODELS,POLICIES,WALLS,residual_metrics
from rans_preservation_oracle import Reference,representation

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--model',choices=MODELS,required=True);a=p.parse_args()
    folder=a.root/('preservation-temporal-'+a.label);grid=load(a.root/'preservation-spatial/evaluation-grid.json');weights=grid['weights']
    result=dict(schema=1,scope='read-only reporting recovery; original failed receipts retained',checks=[],expected_check_names=[],temporal_tables=[])
    for policy in POLICIES:
        for wall in WALLS:
            ref=Reference(a.model,policy,wall);tag='-'.join((a.model,policy,wall));runs=[]
            for j,dt in enumerate((1e-4,5e-5,2.5e-5)):
                name=tag+'-dt'+str(j);result['expected_check_names'].append(name);f=folder/name
                try:
                    row=load(f/'result.json');report=load(f/'case.evolution.json');full=load(f/'state-packet.json');packet=full['observation']
                    assert row['status']=='pass' and row['exit_code']==0 and not report['failure'],'native execution failed'
                    assert report['steps']==round(.0016/dt) and math.isclose(report['time'],.0016,abs_tol=1e-14),'native terminal time/steps differ'
                    assert exact_vector_hash(packet['state'])==exact_vector_hash(report['state']),'replay coefficient identity differs'
                    assert exact_vector_hash(full['reference'])==exact_vector_hash(report['accepted_states'][0]),'initial replay identity differs'
                    assert packet['grid']==grid['grid'],'evaluation grid identity differs'
                    assert full['status']=='pass' and packet['read_only'],'invalid or mutating replay'
                    packet['grid_weights']=weights
                    metrics=dict(status='pass',execution_status='completed',time=report['time'],steps=report['steps'],dt=dt,
                        initial_sha256=exact_vector_hash(report['accepted_states'][0]),reference_sha256=exact_vector_hash(full['reference']),
                        final_spatial_alpha_zero=residual_metrics(packet['audit_alpha_zero']),final_stage_residual=report['stage_audits'][-1],
                        refreshed_trace_defect_max=max(s['refreshed_trace_defect'] for s in report['pimple_attempts']),
                        source_splitting_defects={k:max(abs(s.get(k,0.)) for s in report['pimple_attempts']) for k in ('momentum_splitting_identity_defect','k_same_state_source_storage_defect','epsilon_same_state_source_storage_defect')},
                        extrema=audited_extrema(report),drift=representation(ref,packet,weights),recovered_reporting=True)
                    write_json(f/'temporal-metrics.json',metrics);write_json(f/'endpoint-samples.json',packet)
                    result['checks'].append(dict(name=name,status='pass'));runs.append((packet,metrics))
                except Exception as error:result['checks'].append(dict(name=name,status='fail',failure=repr(error)))
                write_json(folder/'reanalysis.json',result)
            table=dict(model=a.model,kappa_policy=policy,wall=wall,status='fail',fields={})
            if len(runs)==3:
                assert len({m['initial_sha256'] for _,m in runs})==1,'triplet initial states differ'
                ab=compare_samples(runs[0][0],runs[1][0],weights);bc=compare_samples(runs[1][0],runs[2][0],weights)
                for field in ('u','v','k','epsilon_outer','Pi'):
                    d1=ab['fields'][field]['weighted_RMS'];d2=bc['fields'][field]['weighted_RMS'];uncertainty=1e-10
                    measured=min(d1,d2)>5*uncertainty;order=math.log2(d1/d2) if measured else None
                    table['fields'][field]=dict(differences_weighted_RMS=[d1,d2],differences_L2=[ab['fields'][field]['L2'],bc['fields'][field]['L2']],order=order,
                        uncertainty_floor=uncertainty,status='unresolved_floor' if not measured else 'pass' if .8<=order<=1.2 else 'fail')
                table['status']='fail' if any(x['status']=='fail' for x in table['fields'].values()) else 'pass';table['initial_state_sha256']=runs[0][1]['initial_sha256'];table['terminal_time']=.0016
            result['temporal_tables'].append(table)
    result['coverage']=coverage(result['checks'],set(result['expected_check_names']),key=lambda x:x['name'])
    result['status']='pass' if result['coverage']['status']=='pass' and all(x['status']=='pass' for x in result['checks']) and all(x['status']=='pass' for x in result['temporal_tables']) else 'fail'
    write_json(folder/'reanalysis.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
