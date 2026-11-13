#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Narrow fixed-six-millimetre H14 qualification and immutable input adapters.

Historical states below are operator fixtures, never target restart states.
"""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import xml.etree.ElementTree as ET

from run_wj_h4 import load,save,sha
from wj_h14_execution import clone,file_id,validate_pair
from wj_h7_analysis import context
from wj_h4_points import vector_hash
from cls_epawf_m3_cases import set_info,set_parameter,write

H10=Path('/home/jay/build/wj-h10-flat-inlet-20261113')
H13=Path('/home/jay/build/wj-h13-fullheight-paired-domain-20261113')

def verify(v):
    p=Path(v['path'])
    assert p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],str(p)
    return p

def runtime(root):
    rt=root/'runtime';rt.mkdir(exist_ok=True);(rt/'glibc').mkdir(exist_ok=True);(rt/'geometry').mkdir(exist_ok=True)
    def identical_copy(source,dest):
        if dest.exists():assert sha(source)==sha(dest),'partial preparation differs: '+str(dest)
        else:shutil.copyfile(source,dest)
        dest.chmod(0o755)
    old=load(H10/'setup.json');aliases=[]
    for rec in old['private_libraries']+[old['geometry']]:
        src=verify(rec);dest=rt/('geometry'if rec==old['geometry']else'glibc')/src.name
        identical_copy(src,dest)
        assert sha(dest)==rec['sha256'];aliases.append(dict(original=rec,copy=file_id(dest)))
    executables=[]
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe'):
        source=root/'build'/name;dest=rt/(name+'.elf');identical_copy(source,dest)
        wrapper=rt/name
        text='#!/bin/sh\nexec '+shlex.quote(str(rt/'glibc/ld-linux-x86-64.so.2'))+' --library-path '+\
            shlex.quote(str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64')+' '+shlex.quote(str(dest))+' "$@"\n'
        if wrapper.exists():assert wrapper.read_text()==text
        else:wrapper.write_text(text)
        wrapper.chmod(0o755);executables.append(dict(build=file_id(source),ELF=file_id(dest),wrapper=file_id(wrapper)))
    commands=load(root/'build/compile_commands.json')
    assert all('-ffp-contract=off'in x['command']and'-fno-fast-math'in x['command']for x in commands)
    save(root/'science/runtime-build.json',dict(executables=executables,private_libraries=aliases,
        compile_commands=file_id(root/'build/compile_commands.json'),compiled_inputs=[dict(path=x['file'],sha256=sha(x['file']),command=x['command'])for x in commands],
        source_capsules=[file_id(p)for p in sorted(root.glob('*source*.json'))],
        compiler='/usr/bin/c++ GNU11.4.0; Release, serial native rank, compile parallel8',
        scope='new H14 governor/reader, exact H10 original glibc/geometry bytes; actual native-exec maps required separately'))

def inputs(root):
    out=root/'science/small';out.mkdir();bindings=load(H10/'publication-v1/state-bindings.json')
    cp0rec=bindings['coarse-cp0']['producer'];runrec=bindings['coarse-run']['producer']
    cp0=load(verify(cp0rec));r=load(verify(runrec));ctx=context(r,1)
    assert ctx['step']==1 and ctx['alpha']==100 and ctx['time']==.01
    assert vector_hash(cp0['state'])==vector_hash(r['accepted_states'][0])
    xml=H10/'science/cases/coarse/cp0/case.xml'
    assert sha(xml)=='f33d721a2f8bd3ceb6cba40f769b15d10bf6a648a1f9e922e30a02c7dadc3eb3'
    case=clone(xml,out/'operator',dt=.01)
    packet=dict(operation='mapped_backend_checks',inner_stress_model='Linear',state=cp0['state'],history=ctx['history'],
        fixed_reference=cp0['fixed_reference'],fixed_defect=cp0['fixed_defect'],alpha=100.,time=.01,
        state_sha256=vector_hash(cp0['state']),history_sha256=vector_hash(ctx['history']))
    save(case.parent/'request.json',packet)
    original=H10/'science/fixture-v2/case.xml'
    assert sha(original)=='08abd8398c2e682e783a1c434c1d96050b0693fb4997981f38f5d418776bf84c'
    fixture=clone(original,out/'initializer',dt=.01)
    save(fixture.parent/'request.json',dict(operation='initializer_checks',inner_stress_model='Linear'))
    cp=clone(xml,out/'cp0',dt=.01)
    runs={}
    for name,steps,restart in (('two',2,out/'cp0/CP'),('one',1,out/'cp0/CP'),('restart',1,out/'one/CP.1'),('reject',1,out/'one/CP.1')):
        path=clone(xml,out/name,dt=.01,steps=steps,restart=restart)
        if name=='reject':
            x=ET.parse(path).getroot();set_parameter(x,'RANSWallPimpleHardAuditRejectStep',2);write(x,path)
        runs[name]=file_id(path)
    save(out/'intent.json',dict(scope='seven maximum small verification attempts: two parity sweeps, four integration, one deliberate rejection; no target',
        original_cp0_report=cp0rec,original_step1_report=runrec,original_xml=file_id(xml),operator_xml=file_id(case),
        operator_packet=file_id(case.parent/'request.json'),initializer_original=file_id(original),initializer_xml=file_id(fixture),cp0=file_id(cp),runs=runs,
        height_mm=6,steps=7,state=2328,operator_state_role='original CP0 as test vector; original actual step1 callback history; never a checkpoint relabeling'))

def small_checks(root):
    """Compare real new histories; the rejected attempt stays rejected."""
    import numpy as np
    from wj_h9_local import same,phase_bindings
    from wj3_measure import checkpoint_payload
    out=root/'science/small';reports={k:load(out/k/'case.evolution.json')for k in ('cp0','two','one','restart','reject')}
    cp0=reports['cp0'];two=reports['two'];one=reports['one'];restart=reports['restart'];reject=reports['reject']
    assert cp0['initial_velocity_calls']==1 and cp0['initial_velocity_mode']=='FlatInletTraceV1'
    assert cp0['steps']==0 and two['steps']==restart['steps']==2 and one['steps']==1
    err=float(np.max(abs(np.array(two['state'])-restart['state'])))
    assert err<=1e-9*max(1.,max(abs(x)for x in two['state']))
    counts={};phase=[]
    for name,r in reports.items():
        assert r['initial_velocity_mode']=='FlatInletTraceV1'
        for key in ('fixed_reference','fixed_defect'):assert same(r[key],cp0[key])
        if name!='cp0':
            assert r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
        for attempt in r['pimple_attempts']:
            assert attempt['scheduled_complete']
            p=load(attempt['k_snapshot']);assert sha(attempt['k_snapshot'])==attempt['k_snapshot_sha256']
            assert same(p['raw'],p['corrected'])and p['k_bound']['activations']==0 and p['k_bound']['mass_norm_correction']==0
            assert p['matrix_encoding'].startswith('COO native-endian int32 row,int32 column') and p['schema']==2
        counts[name]=dict(steps=r['steps'],time=r['time'],failure=r['failure'],attempts=len(r['pimple_attempts']),
            initialization=r['initialization_lift_calls_iterations'],initializer_calls=r['initial_velocity_calls'],
            producer=file_id(out/name/'case.evolution.json'),checkpoint=r['checkpoint'])
    for key in ('state','native_history_times'):assert same(two[key],restart[key])
    assert len(two['native_history'])==len(restart['native_history'])
    assert all(same(x,y)for x,y in zip(two['native_history'],restart['native_history']))
    assert reject['failure']and reject['steps']==1 and reject['time']==.01 and same(reject['state'],one['state'])
    assert len(reject['physical_range_events'])>0
    for p in(out/'reject').glob('CP*/State.bin'):
        c=checkpoint_payload(p.parent);assert c['steps']<=1
    for name in ('two','one','restart'):
        r=reports[name]
        first=r['steps']+1-len(r['accepted_states'])
        for step in range(first+1,r['steps']+1):
            a=r['pimple_attempts'][step-first-1];p=load(a['k_snapshot']);_,v=phase_bindings(r,p,step);assert v['P_C_bitwise'];phase.append(dict(case=name,**v))
    save(out/'restart-rollback-checks.json',dict(status='PASS',state_max_error=err,state_bitwise=same(two['state'],restart['state']),
        native_history_bitwise=True,reference_defect_bitwise=True,rejected_attempt_rollback=True,counts=counts,phase_bindings=phase))

def view(root,case,reuse=False):
    folder=root/'science/small'/case;rp=folder/'case.evolution.json';r=load(rp);out=folder/'view';out.mkdir()
    save(out/'request.json',dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
        native_checkpoint=r['checkpoint'],checkpoint_sha256=sha(Path(r['checkpoint'])/'State.bin'),
        physical_observations=[[0.,0.,.2375],[4.,0.,.2375],[8.,0.,.2375],[8.,0.,.25607517888672315]],
        weak=False,assembly=False,reuse=reuse))

def stored_reference(manifest=None):
    """Restore ONLY the settled H13 six-mm reference; never run an inverse/fit."""
    from wj_h5_neutral import NeutralReference
    from wj_h2_profiles import wall
    v=load(manifest or H13/'publication-v1/paired-case/paired-domain.json')['primary_reference']
    assert v['H_w']==.15 and v['h_s_mm']==4.341801417014184 and v['physical_z_H_mm']==9.5
    # Deserialization of already validated immutable data, not construction of
    # another candidate through NeutralReference.__post_init__/roughness inverse.
    c=object.__new__(NeutralReference)
    for k,value in dict(H=v['H_w'],h=v['h_s'],z0_e=v['z0_e'],u_tau=v['u_tau'],nu=v['nu'],d0=v['d0'],top=v['top'],
                        d_nu=v['d_nu'],fit_identity=v['fit_identity']).items():object.__setattr__(c,k,value)
    w=wall(v['K'],c.H,c.h,nu=c.nu);assert abs(w['rc']-v['rc'])<2e-14
    object.__setattr__(c,'w',w);assert c.K==v['K']and c.A==v['A'];return c

def assessment(root,relative,version='v1'):
    import math
    import numpy as np
    from wj_h4_points import physical
    from wj_h7_analysis import Rows,classify
    from wj_h9_local import same
    from study_wj_h9 import serial
    from wj3_measure import checkpoint_payload
    folder=root/'science'/relative;rp=folder/'case.evolution.json';r=load(rp);c=stored_reference()
    out=folder/('assessment'if version=='v1'else'assessment-'+version);out.mkdir();before=file_id(rp)
    if r['failure']:
        save(out/'result.json',dict(status='OBSTRUCTED',native_failure=r['failure'],producer=before));return
    cp=checkpoint_payload(r['checkpoint'],include_reference=True)
    for key in ('state','fixed_reference','fixed_defect'):assert same(cp[key],r[key])
    assert cp['signature']==r['signature']and cp['steps']==r['steps']and cp['time']==r['time']and cp['dt']==r['dt']
    pts={};xmin=min(e['lower'][0]for e in r['fields'][0]['mesh']);xmax=max(e['upper'][0]for e in r['fields'][0]['mesh'])
    def add(p,tag):
        p=list(map(float,p))
        if not xmin<=p[0]<=xmax or not c.d0<=p[2]<=c.top:return
        key=tuple(x.hex()for x in p);v=pts.setdefault(key,dict(physical=p,tags=[]))
        if tag not in v['tags']:v['tags'].append(tag)
    g=(-1.,-math.sqrt(3/7),0.,math.sqrt(3/7),1.)
    for e in r['fields'][0]['mesh']:
        lo,hi=e['lower'],e['upper']
        for u in g:
            for v in g:
                x=lo[0]if u==-1 else hi[0]if u==1 else (lo[0]+hi[0])/2+(hi[0]-lo[0])*u/2
                z=lo[2]if v==-1 else hi[2]if v==1 else (lo[2]+hi[2])/2+(hi[2]-lo[2])*v/2
                add(physical(r,[x,0.,z]),'every stored element5x5 Lobatto; all legitimate owners')
    for x in sorted(set([xmin,xmax,(xmin+xmax)/2]+[x for x in (0.,4.,8.,10.,-2.5,2.5,6.5,7.5)if xmin<=x<=xmax])):
        for z in (c.d0,c.d0+c.H,.175,.25,.25607517888672315,.31946961248954786,.375,c.top):add([x,0.,z],'fixed interface/measurement/prior-witness physical heights')
    for packet in r['wj_audit_packets']:
        for w in packet['principal_witnesses'].values():
            p=w['position'];add(p,'genuine native witness')
            for a in (0,2):
                for s in (-1,1):q=p.copy();q[a]+=s*1e-6;add(q,'fixed physical inward/witness offset1e-6')
    inventory=dict(rule='unchanged native hard inventory plus every5x5 stored-element Lobatto/all owners, fixed physical observations/native witnesses+1e-6; Gauss8 physical volume support split at settled wall breakpoints',points=list(pts.values()))
    save(out/'inventory.json',inventory);rows=Rows(r,[v['physical']for v in pts.values()],c)
    if relative.startswith('target/'):
        assert len(r['state'])==24688 and [f['size']for f in r['fields']]==[6493,6493,6493,4429,255,525]
    initial_report=folder.parent/'cp0/case.evolution.json'
    initial_state=load(initial_report)['state']if initial_report.exists()else r['accepted_states'][0]
    initial=rows.evaluate(initial_state,principal=False);results=[]
    first=r['steps']+1-len(r['accepted_states'])
    for j,state in enumerate(r['accepted_states']):
        d=rows.evaluate(state);s=rows.summary(d,initial);active=np.flatnonzero(rows.outer);i=active[np.argmin(d['margin'][active])]
        s.update(min_in_plane=float(d['margin'][i]),min_transverse=float(np.min(d['transverse'][active])),
            min_covariance=float(np.min(d['covariance'])),min_scalar=min(float(np.min(d['k'])),float(np.min(d['epsilon'][active]))),
            witness=rows.witness(d,i),quality_scope='all explicit owners AND physical Gauss8 support',
            divergence_scaled_max=float(np.max(abs(d['divergence'])/np.maximum(1.,np.linalg.norm(d['G'],axis=(1,2))))))
        s['status']=classify(s['min_in_plane'],s['min_transverse'],s['min_covariance'],s['min_scalar']);moments={}
        groups={}
        for n,(p,o)in enumerate(rows.rows[:rows.point_count]):
            if not o['interface_role'].startswith('matching_'):continue
            key=(float(p[0]).hex(),float(o['lower'][0]).hex(),float(o['upper'][0]).hex())
            groups.setdefault(key,{})[o['region']]=n
        interfaces=[]
        for pair in groups.values():
            if set(pair)!={'inner','outer'}:continue
            a,b=pair['inner'],pair['outer'];K=float(d['K'][b]);target=(K/(81/(2*math.sqrt(143))))**1.5/(c.H*.3819178635467331)
            interfaces.append(dict(physical=rows.rows[b][0],inner_owner=rows.rows[a][1],outer_owner=rows.rows[b][1],K=K,rc=float(d['rc'][b]),Lambda=float(d['Lambda'][b]),
                shared_U=(d['U'][a]-d['U'][b]).tolist(),shared_k=float(d['k'][a]-d['k'][b]),epsilon=float(d['epsilon'][b]),epsilon_target=target,
                epsilon_trace_defect=float(d['epsilon'][b]-target),total_traction_jump=(d['T'][a,:,2]-d['T'][b,:,2]).tolist(),
                k_flux_jump=float(d['nu_d'][a]*d['grad_k'][a,2]-d['nu_d'][b]*d['grad_k'][b,2]),reduced_pressure_jump=float(d['pi'][a]-d['pi'][b])))
        s['owner_interfaces']=interfaces
        for b in (4,5):
            ii,(ids,phi,_,_)=rows.maps[b];rhs=np.zeros(r['fields'][b]['size']);np.add.at(rhs,(ids-r['offsets'][b]).ravel(),(phi*(rows.weights[ii]*d['divergence'][ii])[:,None]).ravel())
            moments[str(b)]=dict(moments=rhs.tolist(),maximum=float(np.max(abs(rhs))),signed_integral=float(np.sum(rhs)))
        step=first+j
        if max(v['maximum']for v in moments.values())>2e-15:s['status']='OBSTRUCTED'
        if step==0 and s['divergence_scaled_max']>1e-12:s['status']='OBSTRUCTED'
        results.append(dict(absolute_step=step,time=r['wj_audit_packets'][j]['time'],state_sha256=vector_hash(state),
            summary=s,all_regional_Q1_moments=moments,native_audit_state_sha256=r['wj_audit_packets'][j]['state_sha256']))
    invariants={}
    if r['initial_velocity_calls']:
        diag=r['initial_velocity_diagnostic'];a=diag['pre_state'];b=diag['post_state'];off=r['offsets']
        assert same(a[off[2]:],b[off[2]:])and same(a,r['fixed_reference'])and vector_hash(b)==vector_hash(r['accepted_states'][0])
        essential=set()
        for block in (0,1):
            for e in r['fields'][block]['mesh']:
                ids=np.array(e['ids']).reshape(4,4)+off[block]
                if e['lower'][0]==xmin:essential.update(ids[:,0])
                if e['lower'][2]==0:essential.update(ids[0,:])
                if block==1 and e['upper'][2]==c.top:essential.update(ids[1,:])
        essential=sorted(essential);assert same(np.array(a)[essential],np.array(b)[essential])
        pre=rows.evaluate(a,principal=False);post=rows.evaluate(b,principal=False)
        for k in ('K','rc','Lambda'):assert same(pre[k],post[k])
        invariants=dict(essential_velocity_bits=True,essential_count=len(essential),nonvelocity_bits=True,fixed_reference_defect_bits=True,
            K_cutoff_trace_inputs_unchanged=True,initializer_calls=r['initial_velocity_calls'],velocity_change=rows.norms(post,pre))
    else:
        assert r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
        invariants=dict(initializer_calls=0,project_lift_calls=[0,0])
    status='PASS'if all(v['summary']['status']=='PASS'for v in results)else next(v['summary']['status']for v in results if v['summary']['status']!='PASS')
    save(out/'result.json',serial(dict(status=status,producer=before,checkpoint=file_id(Path(r['checkpoint'])/'State.bin'),
        inventory=file_id(out/'inventory.json'),physical_row_count=len(rows.rows),owner_row_count=rows.point_count,
        states=results,invariants=invariants,diagnostic_mass_solves=0,raw_corrected_k_scope='actual native snapshots/counters verified separately')))
    if status!='PASS':
        from run_wj_h7 import update_journal
        failed=next(v for v in results if v['summary']['status']!='PASS')
        def halt(b):b.setdefault('coverage_halt',dict(status=status,producer=before,absolute_step=failed['absolute_step'],time=failed['time'],witness=serial(failed['summary']['witness'])))
        update_journal(root/'accounting.json',halt)
    assert file_id(rp)==before;print(relative,status,'minima',[v['summary']['min_in_plane']for v in results],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('runtime','inputs','small-checks','view','assessment'),required=True)
    p.add_argument('--case',default='two');p.add_argument('--reuse',action='store_true');p.add_argument('--version',default='v1');a=p.parse_args()
    if a.mode=='view':view(a.root,a.case,a.reuse)
    elif a.mode=='assessment':assessment(a.root,a.case,a.version)
    else:globals()[a.mode.replace('-','_')](a.root)
