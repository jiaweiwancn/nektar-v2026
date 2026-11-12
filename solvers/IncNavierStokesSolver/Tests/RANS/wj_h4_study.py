#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H4 import, native controls and retrospective physical inventories."""
import argparse
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h4 import load,save,sha,THREADS
from wj_h3_analysis import actual,identity,state_context,outer_loads
from wj_h4_points import inventory,scan,owners,evaluate,owner_point,vector_hash,aligned,extend_scan
from validate_wj1 import Probe

PARENT=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h3-20261112')
RUNTIME=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/runtime-v4/IncNavierStokesSolver')
OBSERVER=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/observer-v6/RANSWallResolutionProbe')
KERNEL=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/build-v1/WJKernelProbe')
SOURCE=Path(__file__).resolve().parents[4]
CASES={'base-x1-half':dict(refined=False,nx=1,dt=.005,steps=16,size=1062),
       'refined-x1':dict(refined=True,nx=1,dt=.01,steps=8,size=1966),
       'refined-x1-half':dict(refined=True,nx=1,dt=.005,steps=16,size=1966),
       'base-x3-half':dict(refined=False,nx=3,dt=.005,steps=16,size=2630),
       'local-x3':dict(refined=False,nx=3,dt=.01,steps=8,size=2878,extra_physical_z=[.8818157002137355,1.1947077409162703]),
       'local-x3-half':dict(refined=False,nx=3,dt=.005,steps=16,size=2878,extra_physical_z=[.8818157002137355,1.1947077409162703])}


def runtime_check(root):
    inputs=load(SOURCE/'docs/plans/wj-h2-linux/native-source-inputs.json')
    for f in inputs:assert sha(SOURCE/f['path'])==f['sha256'],f['path']
    libs=load(SOURCE/'docs/plans/wj-h2-linux/loaded-runtime.json')
    for f in libs['loaded_files']:assert sha(actual(f['path']))==f['sha256'],f['path']
    assert sha(RUNTIME)=='9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe'
    assert sha(OBSERVER)=='6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d'
    save(root/'receipts/runtime-equivalence.json',dict(status='PASS',native_inputs=len(inputs),closure_files=len(libs['loaded_files']),
        solver=identity(RUNTIME),observer=identity(OBSERVER),kernel=identity(KERNEL),source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=SOURCE,text=True).strip()))


def bind(report,step=None,rejected=False):
    state,history,alpha,t,p=state_context(report,step,rejected)
    assert vector_hash(state)==p['state_sha256']
    if p['callback_history']:assert vector_hash(p['callback_history'])==p['stage_history_sha256']
    if not rejected:aligned(report['steps'] if step is None else step,report['dt'],t)
    return state,dict(time=t,alpha=alpha,packet=p,state_sha256=p['state_sha256'],callback_history_sha256=p['stage_history_sha256'],
        native_status='rejected candidate' if rejected else 'accepted under original native inventory',step=step)


def assess_report(root,label,rp,inv,first_step=None):
    """Every accepted state, including odd half steps; no callback substitution."""
    folder=root/'scans'/label;folder.mkdir(parents=True,exist_ok=False);report=load(rp);result=[]
    start=report['steps']+1-len(report['accepted_states']);earliest=None
    selections=[(step,False) for step in range(start,report['steps']+1) if first_step is None or step>=first_step]
    if report['failure'] and report['physical_range_events']:selections.append((None,True))
    probe=Probe(str(KERNEL))
    for step,rejected in selections:
        state,context=bind(report,step,rejected);packet=scan(report,state,inv)
        # Independently replay every state's limiting outer action, including
        # accepted states newly found inadmissible by the expanded inventory.
        if packet['minimum'] is not None:
            w=packet['minimum']['row'];G=np.array(w['G']);h=np.outer(w['principal']['polarization'],w['principal']['xi'])
            c=probe.closure(w['k'],w['epsilon'],G);d=probe.closure(w['k'],w['epsilon'],G,(0,0,h))
            replay=float(-np.sum(h*d['R'])/c['scalar'][1]);assert_allclose(replay,w['principal']['in_plane'],rtol=2e-9,atol=2e-11)
            packet['compiled_in_plane_replay']=replay
        if packet['status']!='PASS' and not rejected and earliest is None:earliest=context['time']
        packet.update(context,producer=identity(rp),prefix_qualified=earliest is None,
            qualification='PASS' if packet['status']=='PASS' and earliest is None else 'FAIL' if packet['status']!='PASS' else 'UNQUALIFIED_AFTER_EARLIER_EXPANDED_FAILURE')
        path=folder/(f'step-{step}.json' if not rejected else 'rejected.json');save(path,packet)
        result.append(dict(path=str(path),sha256=sha(path),step=step,time=context['time'],native_status=context['native_status'],
            expanded_status=packet['status'],prefix_qualified=packet['prefix_qualified'],minimum=packet['minimum'],failures=len(packet['failures']),
            state_sha256=context['state_sha256'],callback_history_sha256=context['callback_history_sha256']))
    probe.close();summary=dict(producer=identity(rp),inventory_sha256=inv['sha256'],states=result,earliest_extra_failure=earliest,
        retrospective=True,scope='no claim to an earlier native precommit stop')
    save(folder/'index.json',summary);return summary


def preflight(root):
    runtime_check(root)
    published=load(SOURCE/'docs/plans/wj-h3-linux/runs.json');wanted=('B-zero-x1','B-zero-x3','A-dt02','A-dt01')
    selected={k:next(v['producer'] for v in published if v['case']==k) for k in wanted}
    inv_reports=[]
    for name in (*wanted,'B-frozen-x1','B-frozen-x3'):
        record=next(v['producer'] for v in published if v['case']==name);rp=actual(record['path']);assert sha(rp)==record['sha256'];inv_reports.append(load(rp))
    inv=inventory(inv_reports);save(root/'inventory-v0.json',inv)
    result={}
    for name,record in selected.items():
        rp=actual(record['path']);assert sha(rp)==record['sha256'];result[name]=assess_report(root,'retained-'+name,rp,inv)
    # Check representative existing native samples on their exact declared owner.
    entry=load(PARENT/'native-observations-v1/index.json')['B-zero-x3'];p=actual(entry['native']['path']);assert sha(p)==entry['native']['sha256']
    native=load(p)['weak_checks'];report=load(selected['B-zero-x3']['path'])
    errors=[]
    for row in native['samples'][::191]:
        el=next(v for v in report['fields'][0]['mesh'] if v['geometry']==row['geometry']);o=owner_point(report,el,row['physical'])
        v=evaluate(report,report['state'],row['physical'],o);assert_allclose(v['G'],row['G'],rtol=3e-10,atol=2e-12)
        assert_allclose(v['k'],row['k'],rtol=3e-12,atol=2e-13);errors.append(float(np.max(abs(np.array(v['G'])-row['G']))))
    missed=any(result[k]['earliest_extra_failure'] is not None for k in ('B-zero-x1','B-zero-x3'))
    # A .01 is also an affected native-accepted trajectory if expanded points
    # reveal an omitted negative; A's already-rejected candidate is not a gap.
    missed=missed or result['A-dt01']['earliest_extra_failure'] is not None or result['A-dt02']['earliest_extra_failure'] is not None
    save(root/'preflight.json',dict(status='EXPANDED_GAP' if missed else 'PASS',strict_path_required=missed,
        cases=result,inventory=identity(root/'inventory-v0.json'),native_sample_max_error=max(errors),native_sample_count=len(errors)))
    # The two fixed splits are admitted only as a retained-data normal question.
    # A sign change at its normal face and a small native B margin provide that
    # question; resource admission remains separate and is updated after pilots.
    question=any(s['minimum']['margin']<.01 for k in wanted for s in result[k]['states'])
    save(root/'optional-activation.json',dict(activated=question,resource_admission='PENDING_MEASURED_CORE_COST',
        physical_splits=[.8818157002137355,1.1947077409162703],
        evidence='retained limiting B neighborhood and A normal-face sign/timestep sensitivity',
        preflight_sha256=sha(root/'preflight.json'),no_fresh_evolution_before_choice=True))
    print('Retained expanded preflight',missed,'points',len(inv['points']),flush=True)
    catalog=[]
    for name,record in result.items():
        catalog.extend(dict(case='retained-'+name,**v)for v in record['states'])
    save(root/'catalog.json',catalog);save(root/'active-inventory.json',identity(root/'inventory-v0.json'))
    save(root/'execution-mode.json',dict(strict_required=missed,path_qualified=False,source='retained expanded preflight'))


COUNTERS=('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')
def attempt_counts(report):
    records=report.get('pimple_attempts',[]);rows=[[v.get(k,0) for k in COUNTERS]for v in records]
    # A partial hard exit retains its ORIGINAL failure. Only a falsely successful
    # report requires the all-solves invariant; do not mask the native witness.
    return dict(attempts=len(rows),counts=[sum(v[i]for v in rows)for i in range(5)],
        complete=sum(v==[1]*5 for v in rows),partial=[v for v in rows if v!=[1]*5],
        no_extra_solves=all(all(n in (0,1)for n in v)for v in rows),
        success_counts_valid=bool(report.get('failure')) or all(v==[1]*5 for v in rows),
        native_failure=report.get('failure',''),corrections=sum(v.get('bound_activations',0)for v in records))


def active_inventory(root):return load(load(root/'active-inventory.json')['path'])


def initialize_catalog(root):
    files=[root/v for v in ('catalog.json','active-inventory.json','execution-mode.json')]
    if all(p.exists()for p in files):return
    if any(p.exists()for p in files):raise RuntimeError('partial catalog metadata: inspect, do not reset')
    p=load(root/'preflight.json');inv=p['inventory'];assert sha(inv['path'])==inv['sha256']
    catalog=[]
    for name,record in p['cases'].items():
        for v in record['states']:
            assert sha(v['path'])==v['sha256'];catalog.append(dict(case='retained-'+name,**v))
    save(files[0],catalog);save(files[1],inv)
    save(files[2],dict(strict_required=p['strict_path_required'],path_qualified=False,source='verified retained preflight'))
    save(root/'catalog-bootstrap.json',dict(preflight=identity(root/'preflight.json'),bindings=len(catalog),governing_attempts=0))


def refresh(root,new,tag):
    old=active_inventory(root)
    if new['sha256']==old['sha256']:return
    folder=root/'inventories';folder.mkdir(exist_ok=True);dest=folder/(tag+'.json');save(dest,new)
    catalog=load(root/'catalog.json');updated=[];mode=load(root/'execution-mode.json')
    out=root/'scan-extensions'/tag;out.mkdir(parents=True)
    for i,entry in enumerate(catalog):
        prev=load(entry['path']);rp=prev['producer']['path'];r=load(rp);rejected=entry['native_status']=='rejected candidate'
        state=state_context(r,entry['step'],rejected)[0];v=extend_scan(r,state,prev,old,new)
        if not rejected and v['status']!='PASS':mode.update(strict_required=True,source='expanded inventory extension '+tag)
        path=out/(str(i)+'.json');save(path,v)
        updated.append(dict(entry,path=str(path),sha256=sha(path),expanded_status=v['status'],minimum=v['minimum'],failures=len(v['failures'])))
    earliest={}
    def trajectory(entry):
        return next((name for name in sorted(CASES,key=len,reverse=True) if entry['case'].startswith(name+'-')),entry['case'])
    for entry in updated:
        if entry['native_status']!='rejected candidate' and entry['expanded_status']!='PASS':
            key=trajectory(entry);earliest[key]=min(earliest.get(key,float('inf')),entry['time'])
    for entry in updated:entry['prefix_qualified']=entry['time']<earliest.get(trajectory(entry),float('inf'))
    save(root/'catalog.json',updated);save(root/'active-inventory.json',identity(dest));save(root/'execution-mode.json',mode)


def register_scan(root,label,rp,first_step=None):
    r=load(rp);new=inventory([r],active_inventory(root));refresh(root,new,label)
    result=assess_report(root,label,rp,new,first_step)
    catalog=load(root/'catalog.json');catalog.extend(dict(case=label,**v)for v in result['states']);save(root/'catalog.json',catalog)
    if result['earliest_extra_failure'] is not None:
        mode=load(root/'execution-mode.json');mode.update(strict_required=True,source=label);save(root/'execution-mode.json',mode)
    return result


def native_call(root,name,label,cfg,steps,restart=None):
    from wj_h2_cases import case_xml
    from cls_epawf_m3_cases import set_info,set_parameter,write
    folder=root/'cases'/name/label;folder.parent.mkdir(parents=True,exist_ok=True)
    if restart is None:
        xml=case_xml(folder,'physical',steps,cfg['dt'],cfg['refined'],'ExactStepV1',pressure_seed='VerticalBalanceV1',
            outlet='ZeroTractionV1',xknots=np.linspace(0,8,cfg['nx']+1).tolist(),tangential_order=12,
            extra_physical_z=cfg.get('extra_physical_z'))
    else:
        folder.mkdir();template=root/'cases'/name/'cp0/case.xml';r=ET.parse(template).getroot()
        set_info(r,'RANSWallRestart',str(Path(restart).resolve()));set_info(r,'RANSWallCheckpoint',str(folder/'CP'))
        set_parameter(r,'NumSteps',steps);xml=folder/'case.xml';write(r,xml)
    verify=ET.parse(xml).getroot().find('./CONDITIONS/SOLVERINFO/I[@PROPERTY="RANSWallPimpleVerify"]')
    assert verify.get('VALUE')=='Off'
    request=dict(argv=[str(RUNTIME),str(xml)],runtime=identity(RUNTIME),xml=identity(xml),steps=steps,dt=cfg['dt'],
        restart=identity(Path(restart)/'State.bin') if restart else None,threads=THREADS,role='fresh native CP0' if restart is None else 'governing attempt block')
    save(folder/'request.json',request);start=time.monotonic();env=dict(os.environ,**THREADS,NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'))
    with (folder/'command.log').open('x') as f:
        child=subprocess.run([str(RUNTIME),str(xml)],cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600)
    rp=xml.with_suffix('.evolution.json');r=load(rp) if rp.exists() else None
    record=dict(request,seconds=time.monotonic()-start,exit=child.returncode,log=identity(folder/'command.log'),
        producer=identity(rp) if r else None,accounting=attempt_counts(r) if r else None)
    record['raw_corrected_equal']=True
    if r:
        for a in r['pimple_attempts']:
            if a.get('k_snapshot'):
                snapshot=load(a['k_snapshot'])
                record['raw_corrected_equal'] &= np.asarray(snapshot['raw'],dtype='<f8').tobytes()==np.asarray(snapshot['corrected'],dtype='<f8').tobytes()
    save(folder/'receipt.json',record)
    if r is None:raise RuntimeError('native report absent; preserve command receipt and stop this configuration')
    assert len(r['state'])==cfg['size'] and len(r['state'])<=3000
    return rp,r,record


def outer_check(r):
    result={}
    for rejected in ([False,True] if r['failure'] and r['physical_range_events'] else [False]):
        state,h,alpha,t,p=state_context(r,rejected=rejected)
        v={q:outer_loads(r,state,h,alpha,q)for q in (12,24,32)};comparisons={}
        for q in (12,24):
            errors={name:float(np.max(abs(np.array(a)-v[32]['components'][name])))for name,a in v[q]['components'].items()}
            reg=float(np.max(abs(np.array(v[q]['regional'])-v[32]['regional'])))
            # An absolute-floor-only test is a STRONGER sufficient condition for
            # the original mixed tolerances with any identical exterior load.
            comparisons[str(q)+'-32']=dict(status='PASS' if max(errors.values())<=2e-9 and reg<=4e-9 else 'UNRESOLVED',
                errors=errors,regional_error=reg,original_component_floor=2e-9,original_regional_floor=4e-9)
        result['rejected' if rejected else 'accepted']=dict(time=t,views=v,comparisons=comparisons)
    return result


def prepare(root,name):
    initialize_catalog(root)
    cfg=CASES[name]
    if name.startswith('local-') and not load(root/'optional-activation.json')['activated']:raise ValueError('optional design inactive')
    completed=root/'cases'/name/'cp0/receipt.json'
    if completed.exists():
        # Recover only an actually completed zero-step initialization after a
        # harness failure. This is NOT a zero-launch claim or a governing retry.
        record=load(completed);rp=Path(record['producer']['path']);assert sha(rp)==record['producer']['sha256']
        assert sha(record['xml']['path'])==record['xml']['sha256'] and record['steps']==0 and record['exit']==0
        r=load(rp);assert not r['pimple_attempts'] and r['steps']==0 and r['dt']==cfg['dt'] and len(r['state'])==cfg['size']
        saved=load(rp.with_name('case-design.json'))
        assert saved['refined']==cfg['refined'] and saved['tangential_order']==12 and len(saved['x'])==cfg['nx']+1
        assert saved.get('extra_physical_z',[])==cfg.get('extra_physical_z',[])
        save(root/'cases'/name/'completed-cp0-recovery.json',dict(original=identity(completed),role='completed startup recovered; no new native attempt'))
    else:rp,r,record=native_call(root,name,'cp0',cfg,0)
    result=register_scan(root,name+'-cp0',rp)
    save(root/'cases'/name/'prepared.json',dict(configuration=cfg,cp0=record,expanded=result,
        ready=not r['failure'] and result['earliest_extra_failure'] is None))
    print(name,'CP0',result['states'][-1]['expanded_status'],flush=True)


def reuse(root):
    from measure_wj_h2 import checkpoint_for
    from wj3_measure import checkpoint_payload
    published=load(SOURCE/'docs/plans/wj-h3-linux/runs.json');obs=load(PARENT/'native-observations-v1/index.json');result={}
    for name in ('B-zero-x1','B-zero-x3'):
        rec=next(v for v in published if v['case']==name);rp=actual(rec['producer']['path']);assert sha(rp)==rec['producer']['sha256'];r=load(rp)
        cp,payload=checkpoint_for(dict(xml=str(rp.with_name('case.xml')),report=str(rp)),r['steps'])
        assert sha(Path(cp)/'State.bin')==obs[name]['checkpoint']['sha256']
        for key in ('native_history','native_history_times'):
            assert np.asarray(payload[key],dtype='<f8').tobytes()==np.asarray(r[key],dtype='<f8').tobytes()
        assert payload['reference_sha256']==vector_hash(r['fixed_reference']) and payload['time']==r['time'] and payload['dt']==.01
        initial=next(v for v in published if v['case']==name+'-cp0');ir=load(initial['producer']['path'])
        icp=Path(ir['checkpoint']);ip=checkpoint_payload(icp)
        assert vector_hash(ip['state'])==vector_hash(r['initial']) and ip['reference_sha256']==payload['reference_sha256']
        result[name]=dict(report=identity(rp),checkpoint=identity(Path(cp)/'State.bin'),original_cp0=identity(icp/'State.bin'),
            time=r['time'],dt=r['dt'],steps=r['steps'],signature=r['signature'],reference_sha256=payload['reference_sha256'],
            genuine_completed_callback_sha256=r['wj_audit_packets'][-1]['stage_history_sha256'],
            reused_native_readback=obs[name]['native'],native_readback_sha_verified=sha(obs[name]['native']['path'])==obs[name]['native']['sha256'])
        assert result[name]['native_readback_sha_verified']
    save(root/'receipts/reused-native-bindings.json',dict(status='PASS',trajectories=result,new_governing_attempts=0))


def evolve(root,name):
    prepared=load(root/'cases'/name/'prepared.json');assert prepared['ready'],'CP0 not qualified'
    mode=load(root/'execution-mode.json')
    if mode['strict_required'] and not mode['path_qualified']:raise RuntimeError('strict execution path requires bounded 2 vs1+1 qualification before launch')
    cfg=CASES[name];cp0=load(prepared['cp0']['producer']['path']);cp=cp0['checkpoint'];calls=[];stop=None;remaining=cfg['steps']
    while remaining:
        count=1 if mode['strict_required'] else remaining
        start=cfg['steps']-remaining;rp,r,record=native_call(root,name,f'run-{start:02d}',cfg,count,restart=cp);calls.append(record)
        assessment=register_scan(root,name+f'-run-{start:02d}',rp,first_step=start+1)
        if r['failure'] or record['accounting']['corrections'] or not record['raw_corrected_equal'] or not record['accounting']['success_counts_valid'] or not record['accounting']['no_extra_solves']:
            stop='native_hard_or_correction_exit';break
        if assessment['earliest_extra_failure'] is not None:
            stop='expanded_retrospective_failure';break
        remaining-=count;cp=r['checkpoint']
    result=dict(configuration=cfg,execution_path='strict_one_step' if mode['strict_required'] else 'whole_trajectory',calls=calls,
        accepted_steps=r['steps'],accepted_time=r['time'],requested_target=.08,stop=stop,
        accounting=attempt_counts(r) if len(calls)==1 else dict(attempts=sum(v['accounting']['attempts']for v in calls)),
        outer_checks=outer_check(r),status='PASS_SHORT_CONTROL' if stop is None else 'OBSTRUCTED')
    save(root/'cases'/name/'trajectory.json',result);print(name,result['status'],r['time'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('preflight','prepare','evolve','reuse'),required=True)
    p.add_argument('--case',choices=tuple(CASES));a=p.parse_args()
    if a.mode in ('preflight','reuse'):globals()[a.mode](a.root)
    else:globals()[a.mode](a.root,a.case)
