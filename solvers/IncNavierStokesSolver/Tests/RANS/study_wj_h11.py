#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H11 selected immutable imports, every-step audit and strict common-time views."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from run_wj_h11 import H10,load,save,sha,update_journal,record_coverage_halt
from study_wj_h9 import identity,check,serial,norm_metrics
from wj_h11_analysis import MATRIX,VIEWS,FACE18,semantics,admit_xml,selected,checkpoint_match
from wj_h10_analysis import gates,arithmetic,matrix
from wj_h7_analysis import Rows,context
from wj_h8_analysis import extend_witnesses,FACE,local_rows
from wj_h9_local import same,phase_bindings,MomentAudit
from wj_h9_weak import ProjectedTests,WeakAudit,compare,compare_anchor
from wj_h4_points import vector_hash
from wj3_measure import checkpoint_payload
from study_wj_h6 import configuration
from cls_epawf_m3_cases import set_info,set_parameter,write
from validate_wj1 import Probe


def halt(root,status,reason,**kwargs):
    update_journal(root/'accounting.json',lambda b:b.setdefault('h11_halt',dict(status=status,reason=reason,**kwargs)))


def imports(root):
    out=root/'science';out.mkdir();frozen=root/'source/docs/plans/wj-h10-linux';bindings=load(frozen/'publication-v1/state-bindings.json');verified=[];cases=[]
    for name in ('publication-v1/seal.json','closure/seal.json','closure/accounting-final.json','closure/compiled-inputs.json','closure/runtime-observed.json'):
        a=identity(frozen/name);b=identity(H10/name);assert a['sha256']==b['sha256'];verified.append(b)
    assert load(H10/'accounting.json')['closed']
    for v in load(H10/'closure/compiled-inputs.json')['dependencies']:check(v)
    runtime=load(H10/'science/runtime-build.json')
    for v in runtime['executables']:
        for k in ('ELF','wrapper'):verified.append(check(v[k]))
    for v in runtime['private_libraries']:verified.append(check(v))
    for v in load(H10/'closure/runtime-observed.json')['libraries']:verified.append(check(v))
    for k in ('geometry','kernel'):verified.append(check(runtime[k]))
    inv=next(v for v in load(H10/'publication-v1/manifest.json')['heavy']if v['path']==str(H10/'science/final-common/inventory.json'));check(inv)
    save(out/'inventory-parent.json',load(inv['path']));save(out/'projected-tests.json',load(H10/'science/projected-tests.json'))
    expected={'coarse':'317e2cb182d0a7b93208ed8697926c12df866709a83eece004a1b451893749ba','fine':'896203983504129e507e467a52a6ea1e34c393256bc1ae056b94455fe885f362','quarter':'8f71a0b6c6fd742ac7a6e4d68548965600b98f604424b289f8e9d3b070e23854'}
    for branch,dt,n,target in MATRIX:
        old=bindings[branch+'-run'];check(old['producer']);check(old['checkpoint']);assert old['checkpoint']['sha256']==expected[branch]
        r=load(old['producer']['path']);cp=checkpoint_payload(Path(old['checkpoint']['path']).parent,include_reference=True)
        assert checkpoint_match(cp,r,1)and cp['dt']==dt and cp['order']==1 and cp['eta']==0
        assert same(cp['fixed_reference'],r['fixed_reference'])and same(cp['fixed_defect'],r['fixed_defect'])
        assert same(cp['native_history'],r['native_history'])and same(cp['native_history_times'],r['native_history_times'])
        assert same(cp['native_differential_state'],r['state'][:r['offsets'][4]])
        assert vector_hash(r['accepted_states'][0])==bindings[branch+'-cp0']['invariants']['post_sha256']
        original=Path(old['producer']['path']).with_name('case.xml');xml=ET.parse(original).getroot();folder=out/'cases'/branch;folder.mkdir(parents=True)
        set_parameter(xml,'NumSteps',n);set_info(xml,'RANSWallRestart',str(Path(old['checkpoint']['path']).parent));set_info(xml,'RANSWallCheckpoint',str(folder/'CP'))
        set_info(xml,'RANSWallPimpleSnapshotDirectory',str(folder/'snapshots'));set_info(xml,'RANSWallPimpleArithmetic','On');set_info(xml,'RANSWallPimpleArithmeticTime','.01')
        path=folder/'case.xml';write(xml,path);admit_xml(path,branch,original)
        cx=context(r,1);cases.append(dict(branch=branch,dt=dt,additional=n,target=target,xml=identity(path),original_xml=identity(original),
            report=str(folder/'case.evolution.json'),parent=old,source_context=dict(state_sha256=vector_hash(cx['state']),history_sha256=vector_hash(cx['history']),time=cx['time'],alpha=cx['alpha']),
            imported_history=dict(slots=[vector_hash(x)for x in cp['native_history']],times=cp['native_history_times'],primary_sha256=vector_hash(cp['native_differential_state']))))
    save(out/'input-intent.json',dict(cases=cases,inventory=inv,runtime=runtime,runtime_files=load(H10/'closure/runtime-observed.json')['libraries']+[runtime['kernel']],verified_imports=verified,
        historical_strict='reuse sealed H10 compact bindings/weak packets; no old heavy strict readback or regeneration',scope='25 additional attempts, three restart blocks, no CP0/build'))
    update_journal(root/'accounting.json',lambda b:b.setdefault('intent_sha256',sha(out/'input-intent.json')))
    print('H11 imports PASS, exact runtime/history/reference and immutable additional-step inputs',flush=True)


def parent(root,branch):
    return next(v for v in load(root/'science/input-intent.json')['cases']if v['branch']==branch)


def checkpoints(root,branch,r):
    result=[]
    for directory in sorted((root/'science/cases'/branch).glob('CP*')):
        if not(directory/'State.bin').is_file():continue
        p=checkpoint_payload(directory,include_reference=True);s=p['steps'];assert checkpoint_match(p,r,s)
        assert same(p['fixed_reference'],r['fixed_reference'])and same(p['fixed_defect'],r['fixed_defect'])
        result.append(dict(step=s,time=p['time'],dt=p['dt'],checkpoint=identity(directory/'State.bin'),state_sha256=vector_hash(p['state']),
            native_history=[vector_hash(x)for x in p['native_history']],native_history_times=p['native_history_times'],signature_sha256=__import__('hashlib').sha256(p['signature'].encode()).hexdigest()))
    return result


def scan(root,branch,version='v1'):
    rec=parent(root,branch);r=load(rec['report']);out=root/f'science/audits/{branch}-{version}';out.mkdir(parents=True)
    old=load(rec['parent']['producer']['path']);assert r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
    assert r['signature']==old['signature']and same(r['accepted_states'][0],old['state'])
    for k in ('fixed_reference','fixed_defect'):assert same(r[k],old[k])
    first=r['steps']+1-len(r['accepted_states']);assert first==1
    index=checkpoints(root,branch,r);save(out/'checkpoints.json',index)
    cp=checkpoint_payload(r['checkpoint'],include_reference=True)
    assert same(cp['native_history'],r['native_history'])and same(cp['native_history_times'],r['native_history_times'])
    inv=extend_witnesses(r,load(root/'science/inventory-parent.json'),configuration(6.));save(out/'inventory.json',inv)
    rows=Rows(r,[v['physical']for v in inv['points']],configuration(6.));initial=rows.evaluate(old['accepted_states'][0]);mom=MomentAudit(rows);records=[];bindings=[];projections=[];kernel=Probe(str(H10/'runtime/WJKernelProbe'));queries=0
    try:
        for step in range(first+1,r['steps']+1):
            cx,a,local=selected(r,step);check(dict(path=a['k_snapshot'],bytes=Path(a['k_snapshot']).stat().st_size,sha256=a['k_snapshot_sha256']))
            packet=load(a['k_snapshot']);states,binding=phase_bindings(r,packet,step);bindings.append(dict(branch=branch,**binding))
            for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):assert a[k]==1
            assert a['bound_activations']==0 and a['additional_correction_or_restoration_solves']==0 and a['additional_governing_equation_solves']==0
            data=rows.evaluate(cx['state']);summary=rows.summary(data,initial);summary['divergence_scaled_max']=float(np.max(abs(data['divergence'])/np.maximum(1.,np.linalg.norm(data['G'],axis=(1,2)))))
            pp=mom.projection.measure(data);moments={k:dict(moments=v['pressure_test_moments'],projection_L2=v['continuous_Q1_projection_L2'])for k,v in pp['regions'].items()};summary['moments']=moments
            status=gates(summary,moments,startup=False);ii=np.flatnonzero(rows.outer);j=ii[np.argmin(data['margin'][ii])]
            if data['margin'][j]<.002068763923454:status='QUALITY_LOST'if data['margin'][j]>2e-11 else 'UNRESOLVED'if data['margin'][j]>=-2e-11 else 'OBSTRUCTED'
            j=np.flatnonzero(rows.outer[:rows.point_count])[np.argmin(data['margin'][:rows.point_count][rows.outer[:rows.point_count]])]
            from wj_h3_analysis import principal
            independent=principal(data['k'][j],data['epsilon'][j],data['G'][j]);angle=data['angle'][j];xi=np.array([np.cos(angle),0,np.sin(angle)]);H=np.outer([-xi[2],0,xi[0]],xi)
            value=kernel.closure(data['k'][j],data['epsilon'][j],data['G'][j]);derivative=kernel.closure(data['k'][j],data['epsilon'][j],data['G'][j],(0,0,H));queries+=2
            compiled=-np.sum(H*derivative['R'])/value['scalar'][1];assert max(abs(compiled-data['margin'][j]),abs(independent['in_plane']-data['margin'][j]))<2e-11
            records.append(dict(branch=branch,step=step,time=cx['time'],alpha=cx['alpha'],role='accepted',native='PASS',status=status,
                producer=identity(Path(rec['report'])),state_sha256=vector_hash(cx['state']),history_sha256=vector_hash(cx['history']),summary=serial(summary),compiled_replay=float(compiled),
                frozen_subsolves={k:v for k,v in a.items()if isinstance(v,(float,int,bool))or k in ('corrected_continuity_gauge','refreshed_original_block_norms')}))
            projections.append(dict(branch=branch,step=step,time=cx['time'],**pp))
            if status!='PASS':
                if status in ('OBSTRUCTED','UNRESOLVED'):record_coverage_halt(root/'accounting.json',dict(status=status,trajectory=branch,branch=branch,time=cx['time'],step=step,role='accepted',witness=serial(summary['witness'])))
                halt(root,status,'retrospective every-step qualification failed',branch=branch,step=step,time=cx['time'])
    finally:kernel.close()
    native_status='OBSTRUCTED'if r['failure']else 'PASS'
    if r['failure']:halt(root,'OBSTRUCTED',r['failure'],branch=branch,producer=identity(Path(rec['report'])))
    successful=not r['failure']and r['steps']==rec['target']and all(v['status']=='PASS'for v in records)
    save(out/'states.json',serial(records));save(out/'phases.json',serial(bindings));save(out/'pressure-projections.json',serial(projections))
    arithmetic_steps=[first+i+1 for i,a in enumerate(r['pimple_attempts'])if 'momentum_arithmetic'in a]
    assert arithmetic_steps==([]if branch=='coarse'else [2]if branch=='fine'else [4])
    save(out/'result.json',dict(status='PASS'if successful else 'OBSTRUCTED',native=native_status,producer=identity(Path(rec['report'])),
        reserved=rec['additional'],attempted=len(r['pimple_attempts']),accepted=len(records),qualified=sum(v['status']=='PASS'for v in records),
        native_failure=r['failure'],index=identity(out/'checkpoints.json'),states=identity(out/'states.json'),phases=identity(out/'phases.json'),inventory=identity(out/'inventory.json'),
        arithmetic_steps=arithmetic_steps,diagnostic_mass_factorizations=mom.projection.factorizations,diagnostic_mass_solves=sum(v['mass_solves']for v in projections),kernel_queries=queries))
    if r.get('rejected_state'):save(out/'rejected-candidate.json',dict(producer=identity(Path(rec['report'])),role='rejected; no checkpoint/strict eligibility',state_sha256=vector_hash(r['rejected_state']),original_failure=r['failure']))
    selection=root/'science/audit-selection.json';v=load(selection)if selection.exists()else{};v[branch]=identity(out/'result.json');save(selection,v)
    def reconcile(b):b.setdefault('actual',{})[branch]=dict(attempted=len(r['pimple_attempts']),accepted=len(records),qualified=sum(v['status']=='PASS'for v in records))
    update_journal(root/'accounting.json',reconcile);print(branch,'every-step',successful,'accepted',len(records),flush=True)


def audit(root,branch):
    v=load(root/'science/audit-selection.json')[branch];check(v);return load(v['path'])


def view_request(root,branch,step):
    rec=parent(root,branch);r=load(rec['report']);cx,attempt,local=selected(r,step);assert step in VIEWS[branch]
    a=audit(root,branch);assert a['status']=='PASS';candidates=[v for v in load(a['index']['path'])if v['step']==step];assert candidates
    cp=candidates[0];out=root/f'science/views/{branch}-{step}';out.mkdir(parents=True)
    req=dict(operation='open_checks',inner_stress_model='Linear',native_report=rec['report'],native_report_sha256=sha(rec['report']),native_checkpoint=str(Path(cp['checkpoint']['path']).parent),
        checkpoint_sha256=cp['checkpoint']['sha256'],selected_absolute_step=step,physical_observations=[[8,0,FACE],[8,0,FACE18]],weak=True,weak_outer_order=32,assembly=False,reuse=False)
    save(out/'request.json',req);save(out/'context.json',dict(report=identity(Path(rec['report'])),xml=rec['xml'],checkpoint=cp['checkpoint'],step=step,
        state_sha256=vector_hash(cx['state']),history_sha256=vector_hash(cx['history']),alpha=cx['alpha'],time=cx['time']))


def view_finish(root,branch,step):
    out=root/f'science/views/{branch}-{step}';cx=load(out/'context.json');r=load(cx['report']['path']);n=load(out/'native.json');check(cx['report']);check(cx['checkpoint'])
    assert n['selected_absolute_step']==step and n['report_sha256']==cx['report']['sha256']and n['checkpoint_sha256']==cx['checkpoint']['sha256']
    assert n['state_sha256']==cx['state_sha256']and n['callback_history_sha256']==cx['history_sha256']and n['weak_checks']['alpha_mass']==cx['alpha']and n['open_boundary']['time']==cx['time']
    assert n['initialization_work']==[0,0]and n['initial_velocity_calls']==0 and n['governing_solves']==0
    tests=ProjectedTests.from_export(r,load(root/'science/projected-tests.json'));assert tests.validate_boundary(n['weak_checks'])<2e-11
    c=context(r,step);wk=WeakAudit(r,cx['xml']['path'],tests);values={str(q):wk.evaluate(c['state'],c['history'],c['alpha'],q)for q in (24,32)}
    path=out/'signed-weak.json';save(path,serial(values));anchor=compare_anchor(values['32'],dict(components=n['weak_checks']['native_terms'],regional=[n['weak_checks']['native_inner'],n['weak_checks']['native_outer']]))
    quadrature=compare(values['24'],values['32']);native_moments={str(b):n['spatial_residual'][r['offsets'][b]:r['offsets'][b+1]]for b in(4,5)}
    assert [len(v)for v in native_moments.values()]==[34,70]and max(abs(x)for v in native_moments.values()for x in v)<=2e-15
    status='PASS'if anchor['status']==quadrature['status']=='PASS'else 'UNRESOLVED'
    save(out/'result.json',dict(status=status,context=cx,native=identity(out/'native.json'),signed_weak=identity(path),anchor=anchor,quadrature=quadrature,
        native_Q1_moments=native_moments,norms=n['physical_mass_dual_norms'],open_boundary={k:v for k,v in n['open_boundary'].items()if k not in ('boundary_rows','stations')}))
    if status!='PASS':halt(root,status,'common-time weak anchor/quadrature failed',branch=branch,step=step)
    print(branch,step,'weak',status,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('imports','scan','view-request','view-finish'),required=True);p.add_argument('--branch',choices=('coarse','fine','quarter'));p.add_argument('--step',type=int);p.add_argument('--version',default='v1');a=p.parse_args()
    if a.mode=='imports':imports(a.root)
    elif a.mode=='scan':scan(a.root,a.branch,a.version)
    else:globals()[a.mode.replace('-','_')](a.root,a.branch,a.step)
