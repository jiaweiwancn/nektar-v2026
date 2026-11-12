#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed CT3 cases, fresh histories, primitive-gated same-method milestones.

Invoke through ct_accepted3_command.py; one family at a time. A failed native
trajectory is sealed as a failure, never repaired or silently replaced.
"""
import argparse, hashlib, json, os, shutil, xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS, exact_vector_hash, write_json, audited_extrema
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import require, run, probe, row_value
from ct_accepted2_diagnosis import function_identity, WITNESS
from ct_accepted3_metrics import grid, domain, point, inventory
from ct_accepted3_pressure import projection, lagrange, pressure_action
from ct_accepted2_storage import artifact
from ct_accepted2_audit import mass, values
from ct_accepted_pilot import replay, compare, budget_classification, constitutive_comparison, artifact_hash
from ct_accepted_budget_gate import physical_budgets, paired_coverage
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_picard import set_parameter, set_info, write
from rans_request_manifest import RequestManifest, sha, runtime_identity, canonical

CASES={'A95-XZminus-P3':('A95',3,[.1875],1669),'A95-XZminus-P4':('A95',4,[.1875],2981),
       'A93-Zminus-P3':('A93',3,[],1172),'E94E95-Zminus-P3':('E94E95',3,[],1172)}
CONTROLS=[('native',[32,8,10,10,10],16),('q48f16',[48,12,12,12,12],16),
          ('q64f16',[64,16,16,16,16],16),('q64f32',[64,16,16,16,16],32)]


def make_case(folder, source, c, degree, dt, steps, restart=None):
    folder.mkdir(parents=True,exist_ok=False)
    root=ET.parse(source).getroot();functions=function_identity(root)
    tensor_mesh(root,c['knots'],degree,2,c['x'][-1]-c['x'][0],1.,points=[8,8],xknots=c['x'])
    set_parameter(root,'TimeStep',dt);set_parameter(root,'NumSteps',steps)
    set_info(root,'RANSWallRestart',str(restart) if restart else 'None')
    set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'))
    set_info(root,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'))
    require(function_identity(root)==functions,'mesh_builder_changed_physical_functions')
    write(root,folder/'case.xml')
    actual=dict(c,dt=dt,steps=steps,target_time=dt*steps,segment_duration=dt*steps,
        target_time_semantics='segment duration; absolute endpoint belongs to native report',xml_sha256=sha(folder/'case.xml'))
    write_json(folder/'construction.json',actual)
    return folder/'case.xml'


def signature_inputs(c, xml, runtime):
    return dict(model=c['model'],policy='Natural',wall='Accepted',degree=c['degree'],pressure_degree=c['pressure_degree'],
                x=c['x'],knots=c['knots'],functions=function_identity(ET.parse(xml).getroot()),
                configuration_sha256=sha(xml),runtime=runtime_identity(runtime/'IncNavierStokesSolver'))


def initial(a):
    model,degree,xadd,predicted=CASES[a.case]
    source=a.ct1/'pilot-v3'/(model+'-dt-'+('0.002' if model=='A95' else '0.001')+'-zero')/'case.xml'
    old=a.ct1/'pilot-v3'/(model+'-dt-'+('0.002' if model=='A95' else '0.001'))/'construction.json'
    c=json.loads(old.read_text())
    c['source_template_execution_metadata']={key:c.pop(key) for key in ('dt','steps','target_time','xml_sha256') if key in c}
    c.update(x=sorted(set(c['x']+xadd)),knots=sorted(set(c['knots']+[5.125])),degree=degree,pressure_degree=degree-2)
    points,sides,weights,n=grid(c)
    manifest=RequestManifest(a.output,dict(schema='CT3FreshInitializationV1',case=a.case,contract_sha256=sha(a.contract),
        source_XML_sha256=sha(source),source_construction_sha256=sha(old),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        construction=c,grid_sha256=hashlib.sha256(canonical([points,sides,weights])).hexdigest(),source_sha256=sha(Path(__file__))))
    result=dict(schema='CT3InitialV1',status='fail',case=a.case,construction=c,global_count=n,global_weights=weights,
                grid=points,sides=sides,members={})
    for dt in (.002,.001):
        folder=a.output/str(dt);xml=make_case(folder,source,c,degree,dt,0)
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'zero.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        require(e['exit_code']==0,'fresh_zero_step_native_failure:'+str(dt))
        report=json.loads((folder/'case.evolution.json').read_text())
        packet,pe=probe(a.runtime,folder,dict(operation='functions',grid=points,grid_sides=sides,export_basis=True,export_test_basis=True),'initial')
        require(pe['exit_code']==0 and packet['observation']['state']==report['state'],'fresh_effective_initial_state_mismatch')
        require(len(report['state'])==predicted and predicted<=3000,'actual_native_count_not_confirmed')
        d=domain(packet,report,c);require(d['status']=='pass','initial_physical_domain_failure')
        pp=projection(ET.parse(xml).getroot(),c,degree-2,packet['observation'],packet['test_basis'][-1])
        row=dict(status='pass',slots=len(report['state']),block_counts=packet['observation']['sizes'],
                 identity=signature_inputs(c,xml,a.runtime),state_sha256=exact_vector_hash(report['state']),
                 initial_history_sha256=artifact_hash(report['checkpoint']),checkpoint=report['checkpoint'],
                 native_report_sha256=sha(folder/'case.evolution.json'),packet_sha256=sha(folder/'initial.json'),
                 initialization_receipt=report['initialization_receipt'],domain=d,pressure_projection=pp,executions=[e,pe])
        result['members'][str(dt)]=row
        if dt==.001:require(row['state_sha256']==result['members']['0.002']['state_sha256'],'same_space_initial_coefficients_differ')
        write_json(a.output/'progress.json',result)
    result['status']='pass';write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    return result


def outlet(initial,state,xml,model):
    o=initial['observation'];off=o['offsets'];s=np.array(state);root=ET.parse(xml).getroot()
    f=root.find("./CONDITIONS/FUNCTION[@NAME='OutletTotalTractionOuter']")
    expressions=[f.find("E[@VAR='"+v+"']").get('VALUE').replace('^','**') for v in ('u','v')]
    maximum_x=max(e['upper'][0] for e in o['elements']);rows=[];square=0.;weight=0.
    for q in o['faces']:
        if q['interface'] or not q['region'] or q['coordinate'][0]!=maximum_x:continue
        v=point(q,s,off,model);x,_,y=q['coordinate']
        given=np.array([eval(e,{'__builtins__':{}},dict(x=x,y=y,z=0.,t=0.,sqrt=np.sqrt,abs=np.abs,exp=np.exp,log=np.log)) for e in expressions],float)
        if 'total_stress' in v:
            actual=np.array(v['total_stress'])[[0,2],0];diff=actual-given
            square+=q['weight']*float(diff@diff);weight+=q['weight']
            rows.append(dict(coordinate=q['coordinate'],geometry=q['velocity_k']['geometry'],actual=actual.tolist(),prescribed=given.tolist(),defect=diff.tolist(),Pi=v['Pi'],Tzz=v['total_stress'][2][2]))
        else:rows.append(dict(coordinate=q['coordinate'],status='undefined',invariants=v))
    return dict(weighted_RMS=float(np.sqrt(square/weight)) if weight else None,rows=rows,
                scope='pointwise reconstructed total outlet traction, prescribed weakly, not forced pointwise')


def attempt_diagnosis(initial,report,c):
    attempts=report['pimple_attempts'];selected=[x for i,x in enumerate(attempts) if i in {0,len(attempts)-1}]
    out=[];o=initial['observation'];off=o['offsets']
    for a in selected:
        if 'flow_snapshot' not in a:continue
        f=json.loads(artifact(a['flow_snapshot'],a['flow_snapshot_sha256']).read_text())
        row=dict(time=f['time'],snapshot_sha256=a['flow_snapshot_sha256'],stages={},predictor_is_diagnostic_only=True)
        for name in ('input','momentum_predictor','corrected_flow'):
            state=np.array(f[name]);row['stages'][name]=dict(native=inventory(o['volume'],[q['coordinate'] for q in o['volume']],state,off,c['model']),
                fixed=inventory(initial['test_basis'],o['grid'],state,off,c['model']))
        direction=np.array(f['corrected_flow'])-f['momentum_predictor']
        expected,rough,owners=pressure_action(o,direction,c['roughness_input'],c['pressure_degree'])
        ids=np.array(f['rows'],int);nv=f['velocity_size'];actual=np.array(f['full_pressure_action'])
        defect=float(np.max(abs(expected[ids[:nv]]-actual[ids[:nv]])))
        require(defect<2e-9,'actual_coordinate_pressure_action_mismatch')
        row['full_pressure_action']=dict(maximum_defect=defect,rough_load_maximum=float(np.max(abs(rough))),trace_owners=sorted({x['trace_geometry'] for x in owners}),actual_coordinate_count=len(owners))
        for key in ('predictor_linear_residual','pressure_linear_residual','splitting_identity_defect'):
            row[key]=float(np.max(abs(np.array(f[key]))))
        out.append(row)
    return out


def audits(a,initial_packet,report,c,folder,audit_folder=None):
    result=dict(budgets={},native_report_sha256=sha(folder/'case.evolution.json'),state_sha256=exact_vector_hash(report['state']),time=report['time'])
    result['budget_binding']=dict(state_sha256=result['state_sha256'],native_report_sha256=result['native_report_sha256'],
        time=report['time'],dt=report['dt'],verification_identity=initial_packet['observation']['identity'],model=c['model'],
        configuration_sha256=sha(folder/'case.xml'),functions_sha256=hashlib.sha256(canonical(function_identity(ET.parse(folder/'case.xml').getroot()))).hexdigest())
    result['domain']=domain(initial_packet,report,c);result['mass']=mass(initial_packet,report['state'],c)
    result['outlet_traction']=outlet(initial_packet,report['state'],folder/'case.xml',c['model'])
    n=a.initial_receipt['global_count'];weights=a.initial_receipt['global_weights']
    q=values(initial_packet,report['state'],n);base=values(initial_packet,initial_packet['observation']['state'],n)
    q['observation']['time']=base['observation']['time']=report['time']
    result['actual_initial_drift']=compare(q,base,weights,c)
    g=np.array(q['observation']['grid_gradients']);w=np.array(weights);d=g[:,0,0]+g[:,1,2]
    result['divergence']=dict(weighted_RMS=float(np.sqrt(w@(d*d)/w.sum())),maximum=float(np.max(abs(d))),signed_mean=float(w@d/w.sum()))
    audit_folder=audit_folder or folder
    if audit_folder!=folder:
        audit_folder.mkdir(parents=True,exist_ok=False)
        shutil.copy2(folder/'case.xml',audit_folder/'case.xml')
    request=replay(initial_packet,report['state'],report['time'],initial_packet['observation']['grid'],export_constitutive=True)
    request['grid_sides']=initial_packet['observation']['grid_sides']
    endpoint,e=probe(a.runtime,audit_folder,request,'constitutive')
    result['constitutive_execution']=e
    if e['exit_code']==0:result['independent_constitutive']=constitutive_comparison(endpoint,c)
    else:result['independent_constitutive']=dict(status='failed_native_actual_domain',packet=endpoint)
    if report['steps']>0:
        attempt=next(x for x in reversed(report['pimple_attempts']) if x['scheduled_complete'])
        snapshot=json.loads(artifact(attempt['k_snapshot'],attempt['k_snapshot_sha256']).read_text())
        i=next(i for i,t in enumerate(report['audit_times']) if abs(t-report['time'])<1e-13)
        require(i>0,'missing_consecutive_accepted_state')
        previous=report['accepted_states'][i-1]
        for label,orders,face in CONTROLS:
            sub=audit_folder/label;sub.mkdir();root=ET.parse(folder/'case.xml').getroot()
            for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
            write(root,sub/'case.xml')
            packet,e=probe(a.runtime,sub,replay(initial_packet,report['state'],report['time'],[],export_physical_scalar_budget=True,
                previous_state=previous,stage_history=snapshot['history'],export_frozen_scalar_reconciliation=True,physical_face_order=face),'budget')
            if e['exit_code']==0:
                require(packet['observation']['state']==report['state'] and packet['observation']['identity']==initial_packet['observation']['identity'],'budget_consumed_state_or_configuration_mismatch')
                result['budgets'][label]=dict(budget_classification(packet,c),execution=e,packet_sha256=sha(sub/'budget.json'),binding=result['budget_binding'])
            else:result['budgets'][label]=dict(status='failed_probe',execution=e,packet=packet)
    result['budget_gate']=physical_budgets(result)
    result['attempt_diagnosis']=attempt_diagnosis(initial_packet,report,c)
    return result


def advance(a):
    a.initial_receipt=json.loads((a.initial/'results.json').read_text());ir=a.initial_receipt
    require(ir['status']=='pass' and ir['case']==a.case,'foreign_or_failed_initialization')
    sealed=json.loads((a.initial/'completed-artifacts.json').read_text())
    require(all(sha(a.initial/f)==h for f,h in sealed['files'].items()),'changed_initialization_receipt')
    for name,q in ir['members'].items():
        zero=a.initial/name
        require(sha(zero/'case.evolution.json')==q['native_report_sha256'] and sha(zero/'initial.json')==q['packet_sha256'],'changed_zero_step_artifact')
        require(q['identity']['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign_initial_runtime')
    c=ir['construction'];degree=c['degree']
    if a.target==.008:require(a.parent is None,'first_segment_must_use_fresh_zero_history');start=0.
    else:
        require(a.parent is not None,'continuation_requires_admitted_parent')
        parent=json.loads((a.parent/'results.json').read_text())
        start=.008 if a.target==.02 else .02
        require(parent['case']==a.case and parent['qualification']['status']=='pass' and abs(parent['target']-start)<1e-13,'parent_primitive_gate_not_admitted')
        # Sealed bytes, not an enclosing label, bind the parent primitive receipt.
        seal=json.loads((a.parent/'completed-artifacts.json').read_text())
        require(all(sha(a.parent/f)==h for f,h in seal['files'].items()),'changed_parent_evidence')
        require(parent['qualification']==qualify(parent,start),'parent_primitives_do_not_reclassify')
        for dt,q in parent['members'].items():
            f=Path(q['directory']);native=json.loads((f/'case.evolution.json').read_text())
            require(sha(f/'case.evolution.json')==q['native_report_sha256'] and exact_vector_hash(native['state'])==q['state_sha256'],'parent_native_state_binding_mismatch')
            require(sha(f/'case.xml')==q['identity']['configuration_sha256'] and q['identity']['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'parent_configuration_runtime_binding_mismatch')
            require(artifact_hash(native['checkpoint'])==q['history_sha256'],'parent_native_history_changed')
    recovered={}
    if a.recover_from:
        previous=json.loads((a.recover_from/'request-manifest.json').read_text())['request']
        require(previous['case']==a.case and previous['target']==a.target and previous['initialization_sha256']==sha(a.initial/'results.json')
            and previous['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign_recovery_request')
        recovered={str(dt):sha(a.recover_from/str(dt)/'case.evolution.json') for dt in (.002,.001) if (a.recover_from/str(dt)/'case.evolution.json').exists()}
    manifest=RequestManifest(a.output,dict(schema='CT3MatchedSegmentV1',case=a.case,target=a.target,start=start,dt=[.002,.001],
        initialization_sha256=sha(a.initial/'results.json'),parent_sha256=sha(a.parent/'results.json') if a.parent else None,
        contract_sha256=sha(a.contract),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__)),
        recovery_request_sha256=sha(a.recover_from/'request-manifest.json') if a.recover_from else None,reused_native_reports=recovered))
    result=dict(schema='CT3SegmentV1',case=a.case,target=a.target,start=start,diagnostic_execution='incomplete',members={},temporal={})
    endpoints={}
    for dt in (.002,.001):
        name=str(dt);zero=a.initial/name;initial_packet=json.loads((zero/'initial.json').read_text())
        source=(Path(parent['members'][name]['directory']) if a.parent else zero)
        p=json.loads((source/'case.evolution.json').read_text())
        require(abs(p['time']-start)<1e-13 and abs(p['dt']-dt)<1e-15,'parent_history_time_or_dt_mismatch')
        reused=name in recovered
        folder=a.recover_from/name if reused else a.output/name
        xml=folder/'case.xml' if reused else make_case(folder,source/'case.xml',c,degree,dt,round((a.target-start)/dt),Path(p['checkpoint']))
        require(function_identity(ET.parse(xml).getroot())==function_identity(ET.parse(zero/'case.xml').getroot()),'changed_physical_functions')
        if reused:e=dict(reused_native_report_sha256=recovered[name],execution='original immutable trajectory; diagnostic collector previously failed; not rerun')
        else:
            os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl')
            e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        report=json.loads((folder/'case.evolution.json').read_text())
        require(report['accepted_states'][0]==p['state'],'native_history_consumed_different_parent_state')
        row=dict(time=report['time'],failure=report['failure'],execution=e,steps=report['steps'],directory=str(folder),
            state_sha256=exact_vector_hash(report['state']),initial_state_sha256=ir['members'][name]['state_sha256'],
            native_report_sha256=sha(folder/'case.evolution.json'),identity=signature_inputs(c,xml,a.runtime),
            checkpoint=report['checkpoint'],history_sha256=artifact_hash(report['checkpoint']),native_history_times=report['native_history_times'],
            initialization_receipt=report['initialization_receipt'],audit_extrema=audited_extrema(report),
            maximum_refreshed_trace_lag=max((x.get('refreshed_trace_defect',0.) for x in report['pimple_attempts']),default=0.),
            corrections=sum(x.get('bound_calls',0) for x in report['pimple_attempts']),
            additional_governing_solves=sum(x.get('additional_governing_equation_solves',0) for x in report['pimple_attempts']))
        # Independent audits of a last accepted prefix remain useful after failure.
        row.update(audits(a,initial_packet,report,c,folder,a.output/(name+'-recovered-audits') if reused else None))
        result['members'][name]=row;endpoints[name]=(initial_packet,report)
        write_json(a.output/'progress.json',result)
    aa,bb=endpoints['0.002'],endpoints['0.001']
    common=sorted(set(round(t,12) for t in aa[1]['audit_times']) & set(round(t,12) for t in bb[1]['audit_times']))
    n=ir['global_count'];weights=ir['global_weights'];rows=[]
    # All common accepted times; no unequal-terminal-time spatial comparison.
    for t in common:
        packets=[]
        for ini,r in (aa,bb):
            i=next(i for i,v in enumerate(r['audit_times']) if abs(v-t)<1e-13)
            q=values(ini,r['accepted_states'][i],n);q['observation']['time']=t;packets.append(q)
        rows.append(compare(*packets,weights,c))
    result['temporal']=dict(comparisons=rows,status='pass' if rows and all(x['status']=='pass' for x in rows) else 'fail')
    result['diagnostic_execution']='pass';result['qualification']=qualify(result,a.target)
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']+[str(p.relative_to(a.output)) for p in a.output.glob('*/case.evolution.json')])
    return result


def qualify(result,target):
    rows=result['members'];expected={dt:q['identity'] for dt,q in rows.items()}
    coverage=paired_coverage(rows,target,expected);reasons=list(coverage['reasons'])
    for dt,q in rows.items():
        if q['domain']['status']!='pass':reasons.append(dt+':independent_domain')
        if q['mass']['scaled_defect']>1e-8:reasons.append(dt+':mass')
        if q['independent_constitutive']['status']!='pass':reasons.append(dt+':constitutive_oracle')
        if q['corrections'] or q['additional_governing_solves']:reasons.append(dt+':forbidden_correction_or_solve')
    temporal=result['temporal']['comparisons']
    if not any(abs(q['time']-target)<1e-13 for q in temporal):reasons.append('missing_equal_time_endpoint_comparison')
    if any(q['status']!='pass' for q in temporal):reasons.append('temporal_accuracy')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,coverage=coverage,
                physical_completion='pass' if all(not x['failure'] and abs(x['time']-target)<1e-13 for x in rows.values()) and len(rows)==2 else 'fail',
                spatial_uncertainty='separate_cross_space_receipt_required',resource_compliance='separate_owned_launcher_receipt_required')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['initial','advance'],required=True);p.add_argument('--case',choices=CASES,required=True)
    for x in ('runtime','ct1','contract','output'):p.add_argument('--'+x,type=Path,required=True)
    p.add_argument('--initial',type=Path);p.add_argument('--parent',type=Path);p.add_argument('--recover-from',type=Path);p.add_argument('--target',type=float,choices=[.008,.02,.2])
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    r=initial(a) if a.mode=='initial' else advance(a)
    status=r.get('qualification',r)['status'];print(a.case,a.mode,a.target,status)
    return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
