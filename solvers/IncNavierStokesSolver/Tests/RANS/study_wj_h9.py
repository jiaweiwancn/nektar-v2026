#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seven saved H7/H8 attempts, one bounded passive batch; never a PDE entrypoint."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import shutil
import time
import numpy as np
from run_wj_h9 import load,save,sha
from wj_h4_points import vector_hash,aligned
from wj_h7_analysis import Rows,context
from wj_h9_local import MATRIX,DIRECTION,XI,POL,same,phase_bindings,external_ids,delta_B,telescope,safe_ratio,MomentAudit,actions,matrix_check
from wj_h9_weak import ProjectedTests,WeakAudit,compare
from study_wj_h6 import configuration
from validate_wj1 import Probe

H8=Path('/home/jay/build/wj-h8-local-timestep-20261113')

def identity(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())

def check(v):
    x=identity(v['path']);assert x['sha256']==v['sha256']and x['bytes']==v['bytes'],v['path'];return x

def copied(v,p):
    check(v);p.parent.mkdir(parents=True,exist_ok=True);assert not p.exists();shutil.copyfile(v['path'],p);assert identity(p)['sha256']==v['sha256'];return identity(p)

def stage(root,expected):
    e=load(expected);aliases=[];reports={};attempts=[]
    for branch in ('coarse','fine','quarter'):
        rec=e['reports'][branch];dst=root/'inputs'/branch;rp=copied(rec['report'],dst/'report.json');xml=copied(rec['xml'],dst/'case.xml');r=load(rp['path']);reports[branch]=dict(report=rp,xml=xml,original_report=rec['report'])
        for b,dt,step in MATRIX:
            if b!=branch:continue
            a=r['pimple_attempts'][step-1];v=dict(path=a['k_snapshot'],sha256=a['k_snapshot_sha256'],bytes=Path(a['k_snapshot']).stat().st_size)
            if branch!='quarter':
                known=rec['snapshots'][step-1];assert v==known
            packet=copied(v,dst/f'k-{step}.json');s=load(packet['path']);mp=Path(s['matrix_file']);mat=copied(dict(path=str(mp),bytes=mp.stat().st_size,sha256=s['matrix_sha256']),dst/f'k-{step}.matrix.bin')
            ctx=context(r,step);assert a['scalar_order']=='k_then_epsilon'and not r['failure'];assert len(ctx['state'])==2328;assert r['dt']==dt
            states,binding=phase_bindings(r,s,step)
            changes=[i for i in range(len(states['P']))if not same(states['P'][i:i+1],states['C'][i:i+1])]
            assert set(changes)<=external_ids(r),'C-P includes a non-exterior-lift change'
            binding['external_lift_coefficients']=[dict(row=i,P=float(states['P'][i]),C=float(states['C'][i]),delta=float(states['C'][i]-states['P'][i]))for i in changes]
            assert a['history_sha256']==vector_hash(s['history'])and a['incoming_sha256']==binding['states']['C']
            for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):assert a[key]==1
            assert a['bound_activations']==0
            # Hash both accepted checkpoint endpoints already imported/retained
            # by H8, selected by native absolute step and producer identity.
            cp=e['checkpoints'][branch][str(step)];check(cp)
            attempts.append(dict(branch=branch,dt=dt,step=step,packet=packet,matrix=mat,original_packet=v,original_matrix=identity(mp),binding=binding,checkpoint=cp))
            aliases.extend([dict(original=v,copy=packet),dict(original=identity(mp),copy=mat)])
    anchors={}
    for branch,v in e['anchors'].items():check(v['native']);anchors[branch]=v
    inv=copied(e['inventory'],root/'inputs/H8-inventory.json');records=copied(e['records'],root/'inputs/H8-records.json');check(e['local']);check(e['H8_accounting'])
    runtime=e['runtime'];check(runtime['ELF']);check(runtime['wrapper'])
    for v in runtime['private_libraries']:check(v)
    save(root/'imports.json',dict(reports=reports,attempts=attempts,anchors=anchors,inventory=inv,records=records,local=e['local'],runtime=runtime,H8_accounting=e['H8_accounting'],aliases=aliases,
        source=identity(__file__),new_governing_attempts=0,new_CP0=0,new_strict_readbacks=0,new_builds=0))
    print('PASS exact seven saved packets/matrices, three reports and existing anchors',flush=True)


def serial(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:serial(x)for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x)for x in v]
    return v

def norm_metrics(cache,a):
    out={}
    for side,mask in [('full',np.ones(len(cache.rows),bool)),('inner',~cache.outer),('outer',cache.outer)]:
        w=cache.weights[mask];v=a[mask];out[side]=dict(L2=float(np.sqrt(w@(v*v))),mean=float(w@v/sum(w)),minimum=float(min(v)),maximum=float(max(v)))
    return out


def run(root):
    start=time.monotonic();imp=load(root/'imports.json');out=root/'science';out.mkdir();reps={b:load(v['report']['path'])for b,v in imp['reports'].items()};r0=reps['coarse'];c=configuration(6.)
    inventory=load(imp['inventory']['path']);keys={tuple(v)for v in inventory['local_keys']};points=[v['physical']for v in inventory['points']if tuple(float(x).hex()for x in v['physical'])in keys]
    cache=Rows(r0,points,c);mom=MomentAudit(cache);owned=list(range(cache.point_count));outer_idx=np.flatnonzero(cache.outer[:cache.point_count]);fixed=next(i for i,(p,o)in enumerate(cache.rows[:cache.point_count])if p[0]==8 and abs(p[2]-.25607517888672315)<1e-14 and o['geometry']==17)
    old_records=load(imp['records']['path']);local=[];bindings=[];moments=[];bt=[];tel=[];matrices=[];weak_results={};anchors={};weak_cache={};data_cache={};action_cache={};impulses={};impulse_fields={};phase_ids={};probe_checks=[];compiled_cache={};history_states={};interface_phases=[];flow_increments=[];kernel_calls=0
    save(out/'inventory.json',dict(parent=imp['inventory'],points=points,rows=[dict(physical=p,owner=o)for p,o in cache.rows[:cache.point_count]],rule='H8 local_keys only; all owners; no new point search',fixed_direction=dict(xi=XI.tolist(),a=POL.tolist(),H=DIRECTION.tolist())))
    # Only one full anchor payload is resident at a time; immutable tests are
    # recovered once. No copies of its large scalar rows are published.
    anchor0=load(imp['anchors']['coarse']['native']['path'])
    def validate_anchor(branch,p,step):
        ctx=context(reps[branch],step);assert p['state_sha256']==ctx['packet']['state_sha256']and p['callback_history_sha256']==vector_hash(ctx['history'])
        assert p['weak_checks']['alpha_mass']==ctx['alpha']and p['open_boundary']['time']==ctx['time']
        assert p['report_sha256']==imp['reports'][branch]['report']['sha256']and p['checkpoint_sha256']==imp['anchors'][branch]['context']['checkpoint']['sha256']
    validate_anchor('coarse',anchor0,1);tests=ProjectedTests(r0,anchor0['weak_checks']);weak=WeakAudit(r0,imp['reports']['coarse']['xml']['path'],tests)
    def compact_anchor(p):return dict(components=p['weak_checks']['native_terms'],regional=[p['weak_checks']['native_inner'],p['weak_checks']['native_outer']])
    anchors['coarse']=compact_anchor(anchor0);save(out/'projected-tests.json',tests.export());del anchor0;gc.collect()
    for branch in ('fine','quarter'):
        p=load(imp['anchors'][branch]['native']['path']);validate_anchor(branch,p,2 if branch=='fine'else 4);anchors[branch]=compact_anchor(p);assert p['state_sha256']==context(reps[branch],2 if branch=='fine'else 4)['packet']['state_sha256']
        assert tests.validate_boundary(p['weak_checks'])<2e-11;del p;gc.collect()
    import ast,xml.etree.ElementTree as ET
    first=[load(v['packet']['path'])for v in imp['attempts']if v['step']==1]
    assert len({vector_hash(v['committed'])for v in first})==1
    static=[]
    for branch,item in imp['reports'].items():
        for f in ET.parse(item['xml']['path']).findall('./CONDITIONS/FUNCTION'):
            if f.get('NAME').startswith(('Boundary','Outlet','Top')):
                for equation in f.findall('E'):
                    tree=ast.parse(equation.get('VALUE').replace('^','**'),mode='eval');assert not any(isinstance(n,ast.Name)and n.id=='t'for n in ast.walk(tree))
        static.append(branch)
    save(out/'first-step-premises.json',dict(common_committed_state_sha256=vector_hash(first[0]['committed']),static_boundary_branches=static,
        incoming_differential_storage_max=[float(max(abs(np.array(v['committed'])[:r0['offsets'][4]]-np.array(v['history'])[:r0['offsets'][4]])))for v in first],
        pressure_seed_scope='pressure is not stored in the differential mass term; no governing Schur/predictor assembled',S0_limit='conditional source asymptotics, not computed predictor/correction vectors'))
    setup_seconds=time.monotonic()-start;probe=Probe(imp['runtime']['wrapper']['path'])
    def data(state):
        key=vector_hash(state)
        if key not in data_cache:data_cache[key]=cache.evaluate(state,principal=False)
        return data_cache[key]
    def weak_value(state,history,alpha,q):
        key=(vector_hash(state),vector_hash(history),alpha,q)
        if key not in weak_cache:weak_cache[key]=weak.evaluate(state,history,alpha,q)
        return weak_cache[key]
    def own_fields(d):return {k:v[outer_idx]for k,v in d.items()}
    try:
        for number,item in enumerate(imp['attempts']):
            timer=time.monotonic();branch,step,dt=item['branch'],item['step'],item['dt'];r=reps[branch];packet=load(item['packet']['path']);states,binding=phase_bindings(r,packet,step);binding['external_lift_coefficients']=item['binding']['external_lift_coefficients'];bindings.append(dict(branch=branch,**binding));history=np.array(packet['history']);alpha=packet['alpha']
            ids={key:vector_hash(value)for key,value in states.items()};phase_ids[(branch,step)]=ids;ds={key:data(value)for key,value in states.items()}
            matrices.append(dict(branch=branch,step=step,matrix=item['matrix'],**matrix_check(packet,item['matrix']['path'])))
            mm=mom.identity(states['C'],states['F']);old=next(v for v in old_records if v['branch']==branch and v['step']==step and v['role']=='accepted');assert old['state_sha256']==ids['A']and old['callback_history_sha256']==binding['stage_history_sha256']and old['alpha']==alpha
            for region in mm:
                assert mm[region]['linear_identity_max']<2e-15 and mm[region]['post_flow_max']<2e-15 and mm[region]['repair_max']<2e-15
                prior=old['pressure_projection']['regions'][region]['pressure_test_moments'];assert max(abs(np.array(prior)-mm[region]['F']))<2e-15
            moments.append(dict(branch=branch,step=step,time=packet['time'],regions=mm,F_reuses_A_only_after_bitwise_velocity_equality=True))
            actual={}
            for phase in ('P','C','F','K','A'):
                d=ds[phase]
                if ids[phase]not in action_cache:action_cache[ids[phase]]=actions(d['k'][outer_idx],d['epsilon'][outer_idx],d['G'][outer_idx])
                ad=action_cache[ids[phase]];actual[phase]=ad;assert max(abs(ad['implicit_derivative_error']))<2e-11
                for j,i in enumerate(outer_idx):
                    p,o=cache.rows[i];v={k:np.asarray(x[j]).tolist()for k,x in ad.items()};local.append(dict(branch=branch,step=step,time=packet['time'],phase=phase,state_sha256=ids[phase],physical=p,owner=o,
                        U=d['U'][i].tolist(),G=d['G'][i].tolist(),k=float(d['k'][i]),epsilon=float(d['epsilon'][i]),pi=float(d['pi'][i]),mechanical_pressure=float(d['mechanical_pressure'][i]),**v))
                
                for i in owned:
                    if not cache.outer[i]:
                        local.append(dict(branch=branch,step=step,time=packet['time'],phase=phase,state_sha256=ids[phase],physical=cache.rows[i][0],owner=cache.rows[i][1],
                            U=d['U'][i].tolist(),G=d['G'][i].tolist(),k=float(d['k'][i]),epsilon=None,epsilon_role='no inner PDE',pi=float(d['pi'][i]),R=d['R'][i].tolist(),covariance=np.linalg.eigvalsh(d['R'][i]).tolist(),raw_production=float(d['production'][i]),nu_eff=float(d['nu_eff'][i])))
                trace=cache.interface_history(d);interface_phases.append(dict(branch=branch,step=step,time=packet['time'],phase=phase,state_sha256=ids[phase],**trace));history_states[(branch,step,phase)]=dict(state=states[phase],data=d,action=ad,interface=trace)
            for a,b in (('P','C'),('C','F'),('F','K'),('K','A')):
                diff=delta_B(own_fields(ds[a]),own_fields(ds[b]));assert np.max(abs(diff['residual']))<2e-13
                for j,i in enumerate(outer_idx):
                    bt.append(dict(kind='stage',branch=branch,step=step,time=packet['time'],transition=a+'->'+b,physical=cache.rows[i][0],owner=cache.rows[i][1],**{k:v[j].tolist()for k,v in diff.items()},
                        gradient_norm=float(np.linalg.norm(diff['gradient_term'][j])),timescale_norm=float(np.linalg.norm(diff['timescale_term'][j]))))
            for field in ('dimensional_fixed_action','normalized_fixed_action','in_plane','transverse','N','nu_eff'):
                vv=np.array([actual[p][field]for p in ('C','F','K','A')]);pieces,error=telescope(vv);assert np.max(abs(error))<2e-13
                tel.append(dict(kind='local',branch=branch,step=step,functional=field,maximum_identity_defect=float(np.max(abs(error))),phases=vv.tolist(),increments=pieces.tolist(),lift=(actual['C'][field]-actual['P'][field]).tolist()))
            delta_pressure=states['F']-states['C'];dpi=cache.field(delta_pressure,4)[0]+cache.field(delta_pressure,5)[0]
            flow_increments.append(dict(branch=branch,step=step,time=packet['time'],dt=dt,pressure=norm_metrics(cache,dpi),impulse=norm_metrics(cache,dt*dpi),velocity_increment_L2=float(np.sqrt(cache.weights@np.sum((ds['F']['U']-ds['C']['U'])**2,axis=1)))))
            if step==1:
                pressure=(states['F']-states['C'])*dt;pi=cache.field(pressure,4)[0]+cache.field(pressure,5)[0];impulse_fields[branch]=pi;impulses[branch]=dict(dt=dt,time=packet['time'],coefficient_sha256=vector_hash(pressure[r['offsets'][4]:]),coefficients=pressure[r['offsets'][4]:].tolist(),physical=norm_metrics(cache,pi),
                    pressure_increment=norm_metrics(cache,pi/dt),velocity_increment_L2=float(np.sqrt(cache.weights@np.sum((ds['F']['U']-ds['C']['U'])**2,axis=1))),
                    profiles=[dict(physical=p,owner=o,value=float(pi[i]))for i,(p,o)in enumerate(cache.rows[:cache.point_count])])
                save(out/f'impulse-{branch}.json',impulses[branch])
                pa={q:weak.pressure_action(pressure,q)for q in (24,32)};save(out/f'pressure-impulse-loads-{branch}.json',pa)
            if abs(packet['time']-.01)<1e-14:
                values={phase:{q:weak_value(states[phase],history,alpha,q)for q in (24,32)}for phase in ('C','F','K','A')}
                if not binding['P_C_bitwise']:values['P']={q:weak_value(states['P'],history,alpha,q)for q in (24,32)}
                checks={phase:compare(v[24],v[32])for phase,v in values.items()};anchor=compare(values['A'][32],anchors[branch]);weak_results[branch]=dict(checks=checks,anchor=anchor,values=values)
                # Same history/time/test for every phase: no synthetic BE packet.
                for q in (24,32):
                    vv=[np.array(values[p][q]['regional'])for p in ('C','F','K','A')];increments,error=telescope(vv);assert np.max(abs(error))<2e-13
                    tel.append(dict(kind='weak',branch=branch,step=step,order=q,phases=vv,increments=increments,identity_max=float(np.max(abs(error)))))
                    dp=states['F']-states['C'];act=weak.pressure_action(dp,q);save(out/f'pressure-increment-loads-{branch}-q{q}.json',act)
                    def pressure(v):return sum(np.array(a)for name,a in v['splits'].items()if name=='reduced_pressure_volume'or name=='roughness_pressure'or name.startswith('pressure_reaction_'))
                    err=float(np.max(abs(pressure(values['F'][q])-pressure(values['C'][q])-np.array(act['regional']))));assert err<2e-13
                    tel.append(dict(kind='pressure_linear',branch=branch,step=step,order=q,identity_max=err))
                # Bounded actual-phase shared compiled action replay: one point,
                # four roles. Equal states remain aliases, not new solutions.
                j=list(outer_idx).index(fixed)
                for phase in ('C','F','K','A'):
                    d=ds[phase];k,e,G=d['k'][fixed],d['epsilon'][fixed],d['G'][fixed];key=np.r_[k,e,G.ravel()].astype('<f8').tobytes()
                    if key not in compiled_cache:compiled_cache[key]=(probe.closure(k,e,G),probe.closure(k,e,G,(0,0,DIRECTION)));kernel_calls+=2
                    v,dv=compiled_cache[key]
                    numerical=-float(np.sum(DIRECTION*dv['R']));err=abs(numerical-actual[phase]['dimensional_fixed_action'][j]);assert err<2e-13
                    assert abs(dv['scalar'][12]-actual[phase]['implicit_dN'][j])<2e-11
                    probe_checks.append(dict(branch=branch,step=step,phase=phase,dimensional_action=numerical,independent=actual[phase]['dimensional_fixed_action'][j],error=err,
                        nu_eff=float(v['scalar'][1]),normalized=numerical/v['scalar'][1],implicit_dN_error=float(dv['scalar'][12]-actual[phase]['implicit_dN'][j])))
            save(out/f'attempt-{branch}-{step}.json',serial(dict(binding=binding,continuity=mm,matrix=matrices[-1],local=[v for v in local if v['branch']==branch and v['step']==step],
                delta_B=[v for v in bt if v['branch']==branch and v['step']==step],weak=weak_results.get(branch)if abs(packet['time']-.01)<1e-14 else None,seconds=time.monotonic()-timer)))
            save(out/'progress.json',dict(completed_attempts=number+1,states=len(data_cache),weak_evaluations=weak.evaluations,elapsed=time.monotonic()-start))
            if number==0:
                pilot=time.monotonic()-timer;remaining=(2*pilot+6*max(1.,pilot*.15)+2*setup_seconds+20)*1.25
                ledger=load(root/'accounting.json');pre=ledger['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for k,v in ledger['entries'].items()if k!='science-batch')
                projection=dict(pilot_attempt='coarse1 retained',setup_seconds=setup_seconds,pilot_seconds=pilot,remaining_with_25percent_headroom=remaining,
                    projected_batch_seconds=time.monotonic()-start+remaining,batch_ceiling=360,projected_package_seconds=pre+time.monotonic()-start+remaining+90,package_ceiling=600,closure_reserved=90)
                save(out/'admission.json',projection)
                if projection['projected_batch_seconds']>360 or projection['projected_package_seconds']>600:raise RuntimeError('RESOURCE_LIMIT after retained pilot; remaining batch NOT_RUN')
        same_time=[]
        for a,b in [('coarse','fine'),('fine','quarter')]:
            sa=1 if a=='coarse'else 2;sb=2 if b=='fine'else 4;aa=history_states[(a,sa,'A')]['data'];bb=history_states[(b,sb,'A')]['data'];d=delta_B(own_fields(aa),own_fields(bb));j=list(outer_idx).index(fixed)
            ratios=[]
            for jj,i in enumerate(outer_idx):
                gn=float(np.linalg.norm(d['gradient_term'][jj]));tn=float(np.linalg.norm(d['timescale_term'][jj]));ratios.append(dict(physical=cache.rows[i][0],owner=cache.rows[i][1],gradient_norm=gn,timescale_norm=tn,ratio=safe_ratio(tn,gn),**{k:v[jj]for k,v in d.items()}))
            same_time.append(dict(pair=a+'-'+b,time=.01,records=ratios,fixed_upper_ratio=safe_ratio(float(np.linalg.norm(d['timescale_term'][j])),float(np.linalg.norm(d['gradient_term'][j])))))
        impulse_comparisons=[]
        for a,b in [('coarse','fine'),('fine','quarter')]:
            diff=impulse_fields[b]-impulse_fields[a];metric=norm_metrics(cache,diff);impulse_comparisons.append(dict(pair=a+'-'+b,physical_difference=metric,relative_full_L2=safe_ratio(metric['full']['L2'],impulses[a]['physical']['full']['L2']),
                coefficient_difference=(np.array(impulses[b]['coefficients'])-impulses[a]['coefficients']).tolist()))
        save(out/'interface-phases.json',serial(interface_phases));save(out/'flow-increments.json',serial(flow_increments));save(out/'phase-bindings.json',serial(bindings));save(out/'continuity.json',serial(moments));save(out/'local-phases.json',serial(local));save(out/'delta-B-stage.json',serial(bt));save(out/'delta-B-common-time.json',serial(same_time));save(out/'telescopes.json',serial(tel));save(out/'matrix-checks.json',serial(matrices));save(out/'compiled-checks.json',serial(probe_checks));save(out/'weak-budgets.json',serial(weak_results));save(out/'impulse-comparisons.json',serial(impulse_comparisons))
        native_ok=all(v['anchor']['status']=='PASS'and all(x['status']=='PASS'for x in v['checks'].values())for v in weak_results.values())
        save(out/'summary.json',dict(status='PASS'if native_ok else 'UNRESOLVED',attempts_audited=len(bindings),new_governing=0,new_CP0=0,new_strict=0,new_builds=0,k_matvecs=sum(v['matvecs']for v in matrices),k_factorizations=0,k_solves=0,
            diagnostic_pressure_mass_factorizations=mom.projection.factorizations,diagnostic_pressure_mass_solves=mom.solves,projected_test_fits=tests.fits,kernel_calls=kernel_calls,weak_evaluations=weak.evaluations,unique_actual_states=len(data_cache),
            no_saved_predictor_or_pressure_correction=True,causal_split='UNRESOLVED',seconds=time.monotonic()-start,source=identity(__file__)))
        print('H9 audit',load(out/'summary.json'),flush=True)
    finally:probe.close()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--mode',choices=('stage','science'),required=True);p.add_argument('--expected',type=Path);a=p.parse_args()
    if a.mode=='stage':stage(a.root,a.expected)
    else:run(a.root)
