#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed H10 imports, native inputs, owned CP0/endpoint gates and passive audits."""
import argparse
from pathlib import Path
import shlex
import shutil
import xml.etree.ElementTree as ET
import numpy as np
from run_wj_h10 import load,save,sha
from study_wj_h9 import identity,check,serial,norm_metrics
from wj_h9_local import MomentAudit,actions,DIRECTION
from wj_h9_weak import WeakAudit,ProjectedTests,compare,compare_anchor
from wj_h10_analysis import MATRIX,QUALITY,preservation,gates,arithmetic,matrix
from wj_h7_analysis import Rows,context
from wj_h8_analysis import extend_witnesses,FACE
from wj_h4_points import vector_hash
from study_wj_h6 import configuration
from cls_epawf_m3_cases import set_info,set_parameter,write
from cls_epawf_m4_cases import tensor_mesh
from wj3_measure import checkpoint_payload
from validate_wj1 import Probe

H9=Path('/home/jay/build/wj-h9-startup-compatibility-20261113')

def inputs(root):
    out=root/'science/imports';out.mkdir();old=load(root/'source/docs/plans/wj-h9-linux/publication-v3/imports.json');records=[]
    def copy(v,name):
        check(v);dest=out/name;shutil.copyfile(v['path'],dest);assert sha(dest)==v['sha256'];records.append(dict(original=v,copy=identity(dest)));return dest
    copy(old['inventory'],'inventory.json')
    for b,dt in MATRIX:
        rec=old['reports'][b];copy(rec['xml'],b+'.xml');copy(rec['report'],b+'.report.json')
        a=next(v for v in old['attempts']if v['branch']==b and v['step']==1);copy(a['packet'],b+'.k1.json')
        expected=load(root/'source/docs/plans/wj-h9-linux/publication-v3/manifest.json')['heavy'];v=next(v for v in expected if v['path']==str(H9/f'science-v2/impulse-{b}.json'));copy(v,b+'.impulse.json')
        tree=ET.parse(out/(b+'.xml')).getroot()
        for phase,steps in (('cp0',0),('run',1)):
            folder=root/'science/cases'/b/phase;folder.mkdir(parents=True);x=ET.fromstring(ET.tostring(tree))
            set_info(x,'RANSWallInitialVelocity','FlatInletTraceV1');set_info(x,'RANSWallCheckpoint',str(folder/'CP'))
            set_info(x,'RANSWallRestart','None'if not steps else str(root/'science/cases'/b/'cp0/CP'))
            set_info(x,'RANSWallPimpleSnapshotDirectory',str(folder/'snapshots'));set_info(x,'RANSWallPimpleArithmetic','On');set_info(x,'RANSWallPimpleArithmeticTime',format(dt,'.17g'))
            set_parameter(x,'NumSteps',steps);set_parameter(x,'TimeStep',dt);set_parameter(x,'Time',0.)
            assert [ET.tostring(f)for f in x.findall('./CONDITIONS/FUNCTION')]==[ET.tostring(f)for f in tree.findall('./CONDITIONS/FUNCTION')]
            assert ET.tostring(x.find('GEOMETRY'))==ET.tostring(tree.find('GEOMETRY'))and ET.tostring(x.find('EXPANSIONS'))==ET.tostring(tree.find('EXPANSIONS'))
            write(x,folder/'case.xml')
    make_fixture(root,'fixture')
    fixture=root/'science/fixture'
    save(root/'science/input-intent.json',dict(records=records,cases=[dict(branch=b,dt=dt,cp0=identity(root/f'science/cases/{b}/cp0/case.xml'),run=identity(root/f'science/cases/{b}/run/case.xml'))for b,dt in MATRIX],fixture=identity(fixture/'case.xml'),scope='only added initializer/passive observation/output and authorized dt/step/restart controls; mesh/Functions fixed'))


def small_fixture_tree(base):
    x=ET.fromstring(ET.tostring(base));tensor_mesh(x,[0.,.075,.15,10.075,20.],3,2,8.,1.,points=[8,8],xknots=[0.,4.,8.])
    for e in x.findall('./EXPANSIONS/E'):e.set('FIELDS',e.get('FIELDS').replace('eps_tilde','epsilon'))
    for v in x.findall('./CONDITIONS/VARIABLES/V'):
        if v.text=='eps_tilde':v.text='epsilon'
    return x


def make_fixture(root,name='fixture-v2'):
    x=small_fixture_tree(ET.parse(root/'science/cases/coarse/cp0/case.xml').getroot())
    out=root/'science'/name;out.mkdir();write(x,out/'case.xml');save(out/'request.json',dict(operation='initializer_checks',inner_stress_model='Linear'))


def runtime(root):
    rt=root/'runtime';records=[]
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe'):
        p=root/'build'/name;dest=rt/(name+'.elf');shutil.copyfile(p,dest);dest.chmod(0o755)
        wrapper=rt/name;wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(rt/'glibc/ld-linux-x86-64.so.2'))+' --library-path '+shlex.quote(str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64')+' '+shlex.quote(str(dest))+' "$@"\n');wrapper.chmod(0o755);records.append(dict(ELF=identity(dest),wrapper=identity(wrapper)))
    wrapper=rt/'WJKernelProbe';wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(rt/'glibc/ld-linux-x86-64.so.2'))+' --library-path '+shlex.quote(str(rt/'glibc')+':'+str(rt/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64')+' '+shlex.quote(str(rt/'WJKernelProbe.elf'))+' "$@"\n');wrapper.chmod(0o755)
    compile_commands=load(root/'build/compile_commands.json');flags=sorted(set(v['command'].split(' -o ')[0]for v in compile_commands))
    assert all('-ffp-contract=off'in v and '-fno-fast-math'in v for v in flags)
    save(root/'science/runtime-build.json',dict(executables=records,kernel=identity(rt/'WJKernelProbe.elf'),private_libraries=load(root/'setup.json')['private_libraries'],geometry=load(root/'setup.json')['geometry'],compiler=load(root/'build/CMakeFiles/3.22.1/CMakeCXXCompiler.cmake')if False else '/usr/bin/c++ GNU11.4.0',compiled_inputs=[dict(path=v['file'],sha256=sha(v['file']),command=v['command'])for v in compile_commands if (Path(v['directory'])/shlex.split(v['command'])[shlex.split(v['command']).index('-o')+1]).exists()],configured_commands=identity(root/'build/compile_commands.json'),scope='new governor/reader ELF identities; complete copied original private glibc closure qualified by actual fixture/maps, no host-library substitution'))


def halt(root,status,reason,**details):
    from run_wj_h7 import update_journal
    def change(book):book.setdefault('h10_halt',dict(status=status,reason=reason,**details))
    update_journal(root/'accounting.json',change)


def scan(root,branch,phase,version="v1"):
    folder=root/f'science/cases/{branch}/{phase}';rp=folder/'case.evolution.json';r=load(rp);out=root/('science/analysis/'+branch+'-'+phase+(''if version=='v1'else '-'+version));out.mkdir(parents=True)
    before=identity(rp);assert r['initial_velocity_mode']=='FlatInletTraceV1'and len(r['state'])==2328
    if r['failure']:
        save(out/'result.json',dict(status='OBSTRUCTED',failure=r['failure'],producer=before));halt(root,'OBSTRUCTED',r['failure']);return
    cp=checkpoint_payload(r['checkpoint']);assert vector_hash(cp['state'])==vector_hash(r['state'])and cp['signature']==r['signature']and cp['dt']==r['dt']
    if phase=='cp0':
        invariants=preservation(r)
        if branch!='coarse':
            baseline=load(root/'science/cases/coarse/cp0/case.evolution.json')
            from wj_h9_local import same
            for key in ('state','fixed_reference','fixed_defect'):assert same(r[key],baseline[key])
            other=checkpoint_payload(baseline['checkpoint'])
            assert len(cp['native_history'])==len(other['native_history'])and all(same(a,b)for a,b in zip(cp['native_history'],other['native_history']))
    else:
        invariants=dict(status='PASS',initializer_calls=r['initial_velocity_calls']);assert r['initial_velocity_calls']==0 and r['initialization_lift_calls_iterations']==[0,0]
        cp0=load(root/f'science/cases/{branch}/cp0/case.evolution.json');assert vector_hash(r['accepted_states'][0])==vector_hash(cp0['state'])
        from wj_h9_local import same
        for key in ('fixed_reference','fixed_defect'):assert same(r[key],cp0[key])
    inv=extend_witnesses(r,load(root/'science/imports/inventory.json'),configuration(6.));save(out/'inventory.json',inv)
    cache=Rows(r,[v['physical']for v in inv['points']],configuration(6.));d=cache.evaluate(r['state']);initial=cache.evaluate(r['accepted_states'][0]);summary=cache.summary(d,initial)
    summary['divergence_scaled_max']=float(np.max(abs(d['divergence'])/np.maximum(1.,np.linalg.norm(d['G'],axis=(1,2)))))
    local_keys={tuple(v)for v in inv['local_keys']};local_idx=[j for j,(p,o)in enumerate(cache.rows[:cache.point_count])if tuple(float(v).hex()for v in p)in local_keys]
    oo=[j for j in local_idx if cache.outer[j]];ad=actions(d['k'][oo],d['epsilon'][oo],d['G'][oo]);aa={j:n for n,j in enumerate(oo)};local=[]
    for j in local_idx:
        v=dict(physical=cache.rows[j][0],owner=cache.rows[j][1],U=d['U'][j],G=d['G'][j],k=d['k'][j],epsilon=d['epsilon'][j]if cache.outer[j]else None,pi=d['pi'][j],mechanical_pressure=d['mechanical_pressure'][j],R=d['R'][j],covariance=d['covariance'][j],K=d['K'][j],rc=d['rc'][j],Lambda=d['Lambda'][j],raw_work=d['production'][j])
        if j in aa:v.update({k:np.asarray(a[aa[j]])for k,a in ad.items()})
        local.append(v)
    save(out/'local.json',serial(local))
    mom=MomentAudit(cache);mm=mom.moments(r['state']);summary['moments']={k:dict(moments=v['moments'].tolist(),projection_L2=v['projection_L2'],signed_integral=v['integral'])for k,v in mm.items()}
    state_summary=dict(status=gates(summary,summary['moments'],startup=phase=='cp0'),producer=before,checkpoint=identity(Path(r['checkpoint'])/'State.bin'),invariants=invariants,summary=serial(summary),diagnostic_mass_factorizations=mom.projection.factorizations,diagnostic_mass_solves=mom.solves)
    if phase=='cp0':
        pre=cache.evaluate(r['initial_velocity_diagnostic']['pre_state'],principal=False);state_summary['velocity_change']=cache.norms(d,pre)['full']['U'];state_summary['reference_error']=cache.norms(d,cache.exact)
        assert np.array_equal(pre['K'],d['K'])and np.array_equal(pre['rc'],d['rc'])and np.array_equal(pre['Lambda'],d['Lambda'])
    # One most adverse actual owner with the shared compiled action/root derivative.
    i=np.flatnonzero(cache.outer[:cache.point_count])[np.argmin(d['margin'][:cache.point_count][cache.outer[:cache.point_count]])];from wj_h3_analysis import principal
    independent=principal(d['k'][i],d['epsilon'][i],d['G'][i]);assert abs(independent['in_plane']-d['margin'][i])<2e-11
    angle=d['angle'][i];xi=np.array([np.cos(angle),0.,np.sin(angle)]);a=np.array([-np.sin(angle),0.,np.cos(angle)]);H=np.outer(a,xi)
    probe=Probe(str(root/'runtime/WJKernelProbe'))
    try:
        v=probe.closure(d['k'][i],d['epsilon'][i],d['G'][i]);dv=probe.closure(d['k'][i],d['epsilon'][i],d['G'][i],(0,0,H));value=-float(np.sum(H*dv['R']))/v['scalar'][1]
        assert abs(value-d['margin'][i])<2e-11;state_summary['compiled_replay']=dict(normalized=value,independent=float(d['margin'][i]),point=cache.rows[i],queries=2)
    finally:probe.close()
    if phase=='run':
        attempt=r['pimple_attempts'][0];p=load(attempt['k_snapshot']);assert sha(attempt['k_snapshot'])==attempt['k_snapshot_sha256'];m=load(attempt['momentum_arithmetic']['path']);pp=load(attempt['pressure_arithmetic']['path'])
        for desc in (attempt['momentum_arithmetic'],attempt['pressure_arithmetic']):assert sha(desc['path'])==desc['sha256']
        for v in (m,pp):assert v['history_sha256']==vector_hash(p['history'])and v['time']==r['time']and v['alpha_mass']==1/r['dt']
        result=arithmetic(r,p,m,pp,matrix(m['mass_or_correction_map']['path'],m['mass_or_correction_map']),matrix(pp['mass_or_correction_map']['path'],pp['mass_or_correction_map']),matrix(pp['matrix']['path'],pp['matrix']),matrix(m['matrix']['path'],m['matrix']),mom)
        flow=np.asarray(p['post_flow'])-p['committed'];dpi=cache.field(flow,4)[0]+cache.field(flow,5)[0];result['pressure_physical']=norm_metrics(cache,dpi);result['impulse_physical']=norm_metrics(cache,r['dt']*dpi)
        profiles=[]
        for field,values in [('predictor',result['predictor']),('pressure_correction',result['pressure_correction_velocity']),('total',result['total_velocity_increment']),('predictor_per_dt',result['predictor_per_dt'])]:
            vec=np.zeros(len(flow));vec[np.array(result['velocity_rows'],int)]=values;u=cache.field(vec,0)[0];w=cache.field(vec,1)[0];result[field+'_physical_L2']=float(np.sqrt(cache.weights@(u*u+w*w)))
        for j,(pt,owner)in enumerate(cache.rows[:cache.point_count]):
            if pt[0]in(0.,4.,8.):profiles.append(dict(physical=pt,owner=owner,delta_pi=float(dpi[j]),j_dt=float(r['dt']*dpi[j])))
        result['pressure_profiles']=profiles;result['new_vs_H9_first_step']=dict(old=load(root/f'science/imports/{branch}.impulse.json')['physical'],new=result['impulse_physical'],comparison='different first-step times; same pressure level, no regauge')
        save(out/'arithmetic.json',serial(result));state_summary['arithmetic']=identity(out/'arithmetic.json')
    save(out/'result.json',serial(state_summary));assert identity(rp)==before
    selection=root/'science/analysis-selection.json';index=load(selection)if selection.exists()else {};index[branch+'-'+phase]=identity(out/'result.json');save(selection,index)
    if state_summary['status']!='PASS':halt(root,state_summary['status'],'native state failed expanded/quality/compatibility gate',producer=before,time=r['time'],witness=serial(summary['witness']))
    print(branch,phase,state_summary['status'],'margin',summary['min_in_plane'],'scaled div',summary['divergence_scaled_max'],flush=True)


def selected_analysis(root,branch,phase):
    p=root/'science/analysis-selection.json';index=load(p)if p.exists()else {}
    v=index.get(branch+'-'+phase)
    if v:check(v);return load(v['path'])
    return load(root/f'science/analysis/{branch}-{phase}/result.json')


def view_request(root,branch,phase):
    rp=root/f'science/cases/{branch}/{phase}/case.evolution.json';r=load(rp);ctx=context(r);out=root/f'science/views/{branch}-{phase}';out.mkdir(parents=True)
    cp=Path(r['checkpoint']);save(out/'request.json',dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),native_checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),physical_observations=[[8.,0.,FACE]],weak=True,weak_outer_order=32,assembly=False,reuse=False))
    save(out/'context.json',dict(report=identity(rp),xml=identity(rp.with_name('case.xml')),checkpoint=identity(cp/'State.bin'),state_sha256=vector_hash(ctx['state']),history_sha256=vector_hash(ctx['history']),alpha=ctx['alpha'],time=ctx['time']))


def view_finish(root,branch,phase):
    out=root/f'science/views/{branch}-{phase}';ctx=load(out/'context.json');r=load(ctx['report']['path']);p=load(out/'native.json');check(ctx['report']);check(ctx['checkpoint'])
    assert p['state_sha256']==ctx['state_sha256']and p['callback_history_sha256']==ctx['history_sha256']and p['initialization_work']==[0,0]and p['initial_velocity_calls']==0 and p['governing_solves']==0
    export=root/'science/projected-tests.json'
    if export.exists():tests=ProjectedTests.from_export(r,load(export));assert tests.validate_boundary(p['weak_checks'])<2e-11
    else:tests=ProjectedTests(r,p['weak_checks']);save(export,tests.export())
    wk=WeakAudit(r,ctx['xml']['path'],tests);cx=context(r);key=vector_hash(cx['state'])+'-'+vector_hash(cx['history'])+'-'+str(cx['alpha']);cache=root/'science/weak-cache';cache.mkdir(exist_ok=True);path=cache/(key+'.json')
    if path.exists():values=load(path);evaluations=0
    else:
        values={str(q):wk.evaluate(cx['state'],cx['history'],cx['alpha'],q)for q in(24,32)};save(path,serial(values));evaluations=2
    native_moments={str(b):p['spatial_residual'][r['offsets'][b]:r['offsets'][b+1]]for b in(4,5)}
    assert [len(native_moments[str(b)])for b in(4,5)]==[34,70]
    assert max(abs(x)for v in native_moments.values()for x in v)<=2e-15
    observed=selected_analysis(root,branch,phase)['summary']['moments']
    assert max(abs(x-y)for k,v in native_moments.items()for x,y in zip(v,observed[k]['moments']))<=2e-15
    anchor=compare_anchor(values['32'],dict(components=p['weak_checks']['native_terms'],regional=[p['weak_checks']['native_inner'],p['weak_checks']['native_outer']]))
    quadrature=compare(values['24'],values['32']);status='PASS'if anchor['status']==quadrature['status']=='PASS'else 'UNRESOLVED'
    save(out/'result.json',dict(status=status,context=ctx,native=identity(out/'native.json'),anchor=anchor,quadrature=quadrature,signed_weak=identity(path),native_Q1_moments=native_moments,weak_evaluations=evaluations,norms=p['physical_mass_dual_norms'],open_boundary={k:v for k,v in p['open_boundary'].items()if k not in ('boundary_rows','stations')}))
    if status!='PASS':halt(root,status,'strict/native-independent weak anchor or q24/32 failed')
    print(branch,phase,'strict',status,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('inputs','runtime','fixture','scan','view-request','view-finish'),required=True);p.add_argument('--version',default='v1');p.add_argument('--branch',choices=('coarse','fine','quarter'));p.add_argument('--phase',choices=('cp0','run'));a=p.parse_args()
    if a.mode=='fixture':make_fixture(a.root)
    elif a.mode in ('inputs','runtime'):globals()[a.mode](a.root)
    elif a.mode=='scan':scan(a.root,a.branch,a.phase,a.version)
    else:globals()[a.mode.replace('-','_')](a.root,a.branch,a.phase)
