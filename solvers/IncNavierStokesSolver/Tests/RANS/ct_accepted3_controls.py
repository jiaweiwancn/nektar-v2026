#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT3 missing-coverage negatives and qualified P4/P2 open pressure controls."""
import argparse, copy, json, os, subprocess, sys, xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS, exact_vector_hash, write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted_budget_gate import physical_budgets, paired_coverage
from ct_accepted_open import fixture, probe, run, require, sampling_request, classify
from cls_epawf_m4_cases import tensor_mesh
from ct_accepted3_pressure import pressure_action
from run_cls_epawf_picard import write, set_parameter
from rans_request_manifest import RequestManifest, runtime_identity, sha


def negatives(ct2):
    data=json.loads((ct2/'continue2-audits-v1/results.json').read_text())
    rows=copy.deepcopy(data['cases'])
    good=rows['Zminus-0.002']
    require(physical_budgets(good)['status']=='pass','complete_CT2_budget_not_preserved')
    tests={}
    for name,change in [
        ('missing_outer_epsilon',lambda r:r['budgets']['native']['equations'].pop('outer_epsilon')),
        ('missing_reconstruction',lambda r:r['budgets']['q48f16']['reconstruction'].pop('k')),
        ('missing_control',lambda r:r['budgets'].pop('q64f32')),
        ('nonfinite_primitive',lambda r:r['budgets']['native']['equations']['inner_k'].update(stage_source=float('nan'))),
        ('nonpositive_scale',lambda r:r['budgets']['native']['equations']['outer_k'].update(physical_scale=0.)),
        ('pass_label_failed_accuracy',lambda r:r['budgets']['native']['equations']['outer_k'].update(stage_source=1.,frozen_stage_physical_defect=-1.))]:
        r=copy.deepcopy(good);change(r);v=physical_budgets(r)
        require(v['status']=='fail','missing_coverage_false_pass:'+name)
        tests[name]=dict(status='pass',observed_rejection=v)
    identity={dt:dict(state='actual-'+dt,configuration='cfg-'+dt,runtime='pinned',model='A95') for dt in ('0.002','0.001')}
    members={dt:dict(rows['Zminus-'+dt],time=.2,failure='',identity=identity[dt],initial_state_sha256='same') for dt in identity}
    require(paired_coverage(members,.2,identity)['status']=='pass','complete_CT2_pair_not_preserved')
    for name,change in [
        ('missing_member',lambda r:r.pop('0.001')),
        ('incomplete_endpoint',lambda r:r['0.001'].update(time=.02)),
        ('foreign_identity',lambda r:r['0.002'].update(identity=dict(identity['0.002'],model='A93')))]:
        r=copy.deepcopy(members);change(r);v=paired_coverage(r,.2,identity)
        require(v['status']=='fail','paired_coverage_false_pass:'+name)
        tests[name]=dict(status='pass',observed_rejection=v)
    return dict(status='pass',historical_complete_CT2='both final members retain primitive pass; no historical reclassification',
                ct2_receipt_sha256=sha(ct2/'continue2-audits-v1/results.json'),tests=tests)


def degree_fixture(folder,runtime,pressure=0.,steps=0,roughness=0.):
    xml=fixture(folder,runtime,pressure=pressure,steps=steps)
    root=ET.parse(xml).getroot()
    tensor_mesh(root,[0.,2.,4.],4,2,1.,1.,points=[8,8],xknots=[0.,.5,1.])
    set_parameter(root,'RANSWallRoughness',roughness)
    write(root,xml)
    return xml


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ct2',type=Path,required=True);p.add_argument('--runtime',type=Path,required=True)
    p.add_argument('--output',type=Path);p.add_argument('--intentional-fail',action='store_true')
    a=p.parse_args()
    if a.intentional_fail:
        q=json.loads((a.ct2/'continue2-audits-v1/results.json').read_text())['cases']['Zminus-0.002']
        q['budgets']['native']['equations'].pop('outer_epsilon')
        v=physical_budgets(q);print(json.dumps(v));return 1 if v['status']=='fail' else 0
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT3AdmissionAndP2PressureV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        ct2=sha(a.ct2/'continue2-audits-v1/results.json'),pressure_degree=2,velocity_degree=4,source_sha256=sha(Path(__file__))))
    result=dict(status='fail',coverage=negatives(a.ct2),checks={},executions=[])
    try:
        z=subprocess.run([sys.executable,'-B',__file__,'--ct2',str(a.ct2),'--runtime',str(a.runtime),'--intentional-fail'],capture_output=True,text=True)
        require(z.returncode==1 and 'missing_or_foreign_budget_equation' in z.stdout,'required_failure_exit_not_propagated')
        result['negative_command_exit']=dict(status='pass',exit_code=z.returncode,stdout=z.stdout)
        packets={}
        for name,pressure in [('base',0.),('shifted',.125)]:
            folder=a.output/name;degree_fixture(folder,a.runtime,pressure=pressure)
            packet,e=probe(a.runtime,folder,sampling_request(),'sample');result['executions'].append(e)
            require(e['exit_code']==0,'P4_P2_native_construction_failed')
            result['checks'][name]=classify(packet);packets[name]=packet
        b,s=packets['base'],packets['shifted'];off=b['observation']['offsets']
        external={int(i) for i,_ in b['open_frozen_systems']['external_constraints']}
        free=sorted(set(range(len(b['observation']['state'])))-external)
        d=np.array(s['open_frozen_systems']['raw_residual'])-b['open_frozen_systems']['raw_residual']
        require(np.max(abs(d[free]))<2e-10,'P2_absolute_traction_pressure_offset_action')
        result['checks']['absolute_traction_action']=dict(status='pass',maximum_defect=float(np.max(abs(d[free]))),shift=.125)
        folder=a.output/'rough';degree_fixture(folder,a.runtime,roughness=.1)
        packet,e=probe(a.runtime,folder,dict(operation='functions',grid=[],export_basis=True,export_open_verification=True),'sample');result['executions'].append(e)
        require(e['exit_code']==0,'P2_rough_fixture_failed')
        o=packet['observation'];f=packet['open_frozen_systems'];direction=np.zeros(len(o['state']));off=o['offsets']
        for block in (4,5):
            mass=np.array(f['masses'][block]);load=np.zeros(len(mass))
            for q in o['pressure_mass']:
                if q['region']==block-4:
                    v=q['pressure'];load[np.array(v['ids'],int)]+=q['weight']*np.array(v['phi'])*q['coordinate'][0]**2
            direction[off[block]:off[block+1]]=np.linalg.solve(mass,load)
        expected,rough,provenance=pressure_action(o,direction,.1,2)
        flow=next(x for x in f['blocks'] if x['block']==0);ids=np.array(flow['rows'],int);nv=sum(ids<off[4])
        actual=np.array(flow['matrix'])[:nv]@direction[ids]
        defect=float(np.max(abs(actual-expected[ids[:nv]])))
        deriv=max(abs(q['dpdx']-2*q['coordinate'][0]) for q in provenance)
        # Reusing one face quadrature row cannot reproduce derivative 2*x.
        wrong=max(abs(provenance[0]['dpdx']-q['dpdx']) for q in provenance)
        require(defect<2e-10 and deriv<2e-10 and wrong>.1,'P2_actual_coordinate_rough_pressure_action')
        result['checks']['rough_p_x_squared']=dict(status='pass',maximum_action_defect=defect,maximum_derivative_defect=deriv,
            single_row_derivative_failure=wrong,rough_load_maximum=float(np.max(abs(rough))),provenance=provenance,unknowns=len(direction))
        reports={}
        for name,pressure in [('evolved_base',0.),('evolved_shifted',.125)]:
            folder=a.output/name;degree_fixture(folder,a.runtime,pressure=pressure,steps=2)
            e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(e)
            report=json.loads((folder/'case.evolution.json').read_text());require(e['exit_code']==0 and not report['failure'],'P4_P2_absolute_traction_evolution_failed')
            reports[name]=report
        expected=np.array(s['observation']['state'])-b['observation']['state']
        d=np.array(reports['evolved_shifted']['state'])-reports['evolved_base']['state']-expected
        require(np.max(abs(d))<2e-10,'P2_evolved_pressure_offset')
        result['checks']['evolved_absolute_traction']=dict(status='pass',maximum_defect=float(np.max(abs(d))),base_state_sha256=exact_vector_hash(reports['evolved_base']['state']))
        result['status']='pass'
    except Exception as e:result['failure']=repr(e)
    write_json(a.output/'results.json',result)
    manifest.seal(a.output,['results.json'])
    print(result['status'],result.get('failure',''));return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
