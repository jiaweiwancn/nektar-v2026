#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent current equations and complete histories for native A/B/A replays."""
import argparse,json,math,os,resource,sys,time
from pathlib import Path
from types import SimpleNamespace
sys.dont_write_bytecode=True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import coupled_oracle as core
import history_recovery_check as hr
import physical_cache_history_check as history
import physical_transient_oracle as transient
import implicit_algorithm_check as fields
import cutoff_owner_transition_check as transition
require,read=core.require,history.read
ACCEPTED=('warm-a-first','warm-b','warm-a-return','retry-a','cold-a','cold-b')
SOURCE='56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
TABLE='b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
RESIDUAL_KEYS = ('momentum', 'continuity', 'gauge', 'k', 'epsilon', 'trace')


def paths(manifest,folder):
    require(manifest['schema']=='WJH27CoupledReturnManifestV1','wrong return manifest')
    root=hr.resolve(folder,manifest['root'])
    return {k:history.input_path(v,root)[0]for k,v in manifest['runs'].items()},root


def current_case(manifest,folder,table,name,deadline,out):
    require(name in ACCEPTED,'not an accepted replay');runs,_=paths(manifest,folder);packet=read(runs[name]);d=packet['state']
    require(packet['schema']=='WJH27CoupledReturnCaseV1'and packet['normal_main']is False and
        packet['actual_native_Coupling']is True and packet['accepted']is True and packet['case']==name,
        'input must honestly identify the native diagnostic replay')
    require(d['source_identity']==SOURCE and d['reference_identity']==TABLE,'native source/table differs')
    args=SimpleNamespace(nx=20,nr=24,refine_nx=24,refine_nr=32,direct=False,direct_seconds=40.)
    out.update(native=hr.identity(runs[name]),native_source_identity=d['source_identity'],native_signature=d['signature'],
        origin='Separate public-Coupling probe linked to exact source19 libraries; not normal-main application')
    core.audit(d,read(table),args,out,deadline)
    require(out['current_coupled_step_gate_pass']is True,'independent complete current equation gate failed')
    require(out['negative_controls']['frozen_prior_scalar_and_source']['status']=='DETECTED','frozen scalar/source mutation not detected')


def fixed_time_equations(data, table, state, when, face, deadline):
    """Independent zero-mass current actions at the exact saved state and clock."""
    mms = core.coupled_mms_check.IndependentMMS(
        table, data['mms_body_amplitude'], data['mms_variant'])
    assemblies = []
    for nx, nr in ((20, 24), (24, 32)):
        oracle = core.Oracle(data, mms, nx, nr, deadline)
        load = oracle.forcing(state, when)
        advection = oracle.advection(state, when)[0]
        roots = []
        for element in oracle.elements:
            oracle.clock()
            x, r, weight = oracle.rule(element, state, when)
            xl, xr, rl, rr = element['bounds']
            area = (xr-xl)*(rr-rl)
            require(np.all(weight > 0) and abs(np.sum(weight)-area) <= 1e-12*area,
                    'owner rule has invalid weights or area')
            exact = oracle.exact(x, r, when, bool(element['epsilon']))
            bu = core.basis(element, 'u', x, r)[0]
            uids = np.asarray(element['u'])
            load['flow'][uids] -= bu.T@(weight*exact['u'][3])
            load['flow'][oracle.nu+uids] -= bu.T@(weight*exact['w'][3])
            for key in ('k', 'epsilon'):
                if element[key]:
                    basis = core.basis(element, key, x, r)[0]
                    load[key][element[key]] -= basis.T@(weight*exact[key][3])
            if min(abs(rl-face), abs(rr-face)) > 1e-12:
                continue
            for point in oracle.cutoff_x_splits(element, state, None, (xl, xr, face, face)):
                k_trace = oracle.trace(element, 'k', np.array([point]), state)
                error = abs(float(oracle.cutoff(k_trace)[0])-face)
                if error <= 1e-12*mms.H:
                    roots.append(dict(owner=element['geometry'], x=point,
                                      x_column=[xl, xr], cutoff_residual=error))
        require(len(roots) == 2 and abs(roots[0]['x']-roots[1]['x']) <= 1e-10
                and roots[0]['x_column'] == roots[1]['x_column'],
                'exact saved state does not retain both cutoff-face owners')
        base = core.combine(load, 1., advection, -1.)
        _, norms = oracle.residual(state, when, 0., base)
        trace = oracle.trace_error(state, when)
        require(all(value <= 1e-10 for value in trace.values()),
                'independent exact-state essential/dependent trace differs')
        norms['trace'] = max(trace.values())
        assemblies.append(dict(rule=[nx, nr], residuals=norms, trace=trace,
                               roots=roots))
    difference = {key: abs(assemblies[0]['residuals'][key]-assemblies[1]['residuals'][key])
                  for key in RESIDUAL_KEYS}
    require(all(value <= 1e-10 for value in difference.values()),
            'independent exact-state quadrature difference exceeds 1e-10')
    return dict(assemblies=assemblies, quadrature_difference=difference)


def compare(manifest,folder,table,deadline,out):
    runs,root=paths(manifest,folder);packets={k:read(p)for k,p in runs.items()};data={k:p['state']for k,p in packets.items()};a=data['warm-a-first']
    require(list(runs)==list(ACCEPTED)+['failed-a'],'replay order differs');ref=fields.po.Reference(read(table));owners,_=fields.inventory(a,ref)
    original={}
    for key in('A','B'):
        p,pid=history.input_path(manifest['inputs'][key],root);s,meta=history.binary_history(p,a,.02,deadline,pid['sha256']);original[key]=dict(path=p,states=s,meta=meta)
    stale_path, stale_identity = history.input_path(manifest['inputs']['stale_A'], root)
    _, stale_meta = history.binary_history(stale_path, a, .02, deadline, stale_identity['sha256'])
    require(stale_meta['step'] == 10 and all(
        stale_meta['hashes'][role] == original['B' if role == 'advection' else 'A']['meta']['hashes'][role]
        for role in history.HISTORY_NAMES), 'stale mutation changed other history bytes')
    require(stale_meta['hashes']['advection'] != original['A']['meta']['hashes']['advection'],
            'stale explicit load is identical')
    rejection_path, rejection_identity = history.input_path(manifest['stale_rejection'], root)
    rejection = read(rejection_path)
    require(rejection['source_path'] == str(stale_path)
            and 'restored explicit history differs' in rejection['reason'],
            'stale history was not rejected by its native reconstruction guard')
    geometry_path, _ = history.input_path(manifest['crossing_geometry'], root)
    cutoff_face = read(geometry_path)['cutoff_face']
    table_data = read(table)
    fixed = {key: fixed_time_equations(a, table_data, value['states']['current'],
                                       value['meta']['time'], cutoff_face, deadline)
             for key, value in original.items()}
    roots_a, roots_b = (fixed[key]['assemblies'][-1]['roots'] for key in ('A', 'B'))
    require(roots_a[0]['x'] > .5 > roots_b[0]['x'], 'A/B does not cross the x owner edge')
    states={};binaries={};equations={};events={};last_age=0;previous_calls=0
    exact_inputs = {}
    for index,name in enumerate(ACCEPTED+('failed-a',)):
        fields.clock(deadline);p=packets[name];d=data[name];source='B'if name in('warm-b','cold-b')else'A';origin=original[source]
        require(p['schema']=='WJH27CoupledReturnCaseV1'and p['normal_main']is False and p['actual_native_Coupling']is True and
            p['case']==name and p['sequence']==index,'native replay identity/order differs')
        for key in fields.EXACT_KEYS:
            if key not in('step','time'):require(d[key]==a[key],'native configuration differs: '+key)
        require(d['source_identity']==SOURCE and d['reference_identity']==TABLE and d['all_row_action_crosschecks']is True,'source/current-action checks differ')
        checks=d['current_action_checks'];require(checks['all_row_calls']>0 and
            all(math.isfinite(v)and 0<=v<=1e-10 for v in checks['maximum_scaled_differences_flow_k_epsilon']),
            'current all-row direct/assembled action check failed')
        start=10 if source=='A'else 18;when=start*.02;expect=name!='failed-a'
        require(p['accepted']is expect and p['instance']==(0 if index<4 else 1 if index<6 else 2)and
            p['explicit_reset']is(name.startswith('cold-'))and p['injected_retry']is(name=='retry-a'),
            'native instance/reset/retry policy differs')
        h=d['input_history'];require(h['step']==start and h['time']==when and h['dt']==.02,'input scientific clock differs')
        audit = p['exact_input_current_audit']
        require(audit['time'] == when and audit['step'] == start and
                all(math.isfinite(audit[key]) and audit[key] > 0 for key in history.GUARDS),
                'exact saved state clock/admissibility differs')
        reference = fixed[source]['assemblies'][-1]['residuals']
        delta = {key: abs(audit['residuals'][i]-reference[key])
                 for i, key in enumerate(RESIDUAL_KEYS)}
        require(all(value <= 1e-10 for value in delta.values()),
                'exact-state native/independent residual mismatch')
        exact_inputs[name] = dict(source=source, native=audit, differences=delta)
        prefix=runs[name].with_name(name);incoming=Path(str(prefix)+'.input.bin');final=Path(str(prefix)+'.final.bin')
        before,bid=history.binary_history(incoming,d,.02,deadline);after,fid=history.binary_history(final,d,.02,deadline)
        for role in history.HISTORY_NAMES:
            require(transient.state_hash(before[role])==origin['meta']['hashes'][role]and
                transient.state_hash(core.state(h[role],d['sizes']))==origin['meta']['hashes'][role],
                'complete replay input differs from actual bound history')
        if expect:
            require(incoming.read_bytes()==origin['path'].read_bytes()and d['signature']==a['signature'],
                'accepted replay changed native input bytes/signature')
            states[name]=transition.states(d)
            for role,key in zip(history.HISTORY_NAMES,history.HASH_NAMES):require(fid['hashes'][role]==d[key],'final binary/report history differs')
            require(d['step']==start+1 and d['time']==(start+1)*.02 and d['dt']==.02 and
                max(d['residuals'])<=3e-12 and all(math.isfinite(d[k])and d[k]>0 for k in history.GUARDS),
                'accepted clock/residual/admissibility failed')
            report_path,rid=history.input_path(manifest['current_equations'][name],root);j=read(report_path)
            require(j['native']['sha256']==hr.identity(runs[name])['sha256']and j['current_coupled_step_gate_pass']is True and
                j['native_source_identity']==SOURCE and j['negative_controls']['frozen_prior_scalar_and_source']['status']=='DETECTED',
                'separate current equation/scalar-staleness proof failed')
            equations[name]=dict(file=rid,status=j['status'],residuals=j['independent_assemblies'][-1]['current_residual'])
            stats=fields.cache_policy_statistics(d,'CurrentDataDriftV1');es=d['coupled_newton']['preconditioner_events'];events[name]=es
            require(es and len(es)==stats['refreshes']+stats['reuses']and d['rebuilds']==3*stats['refreshes'],'inverse events/counts incomplete')
            if index<4:
                require(checks['all_row_calls']>previous_calls,'live instance did not execute fresh current actions');previous_calls=checks['all_row_calls']
            for i,event in enumerate(es):
                require(event['alpha']==75 and event['age_used']==(1 if event['refreshed']else event['age_before']+1),'wrong BDF alpha or inverse age')
                if name in('warm-b','warm-a-return')and i==0:
                    require(event['reason']=='reuse'and event['refreshed']is False and event['age_before']==last_age,
                        'return/changed-state base did not retain a live inverse')
            if name in('warm-a-first','cold-a','cold-b'):require(es[0]['reason']=='cold'and es[0]['age_before']==0,'cold control was warm')
            if name=='retry-a':
                require(d['retries']==1 and es[0]['attempt']==0 and es[0]['reason']=='injected'and
                    any(v['attempt']==1 for v in es)and all(v['reason']=='retry'and v['refreshed']for v in es if v['attempt']==1)and
                    d['coupled_newton']['counts_complete']is False,'injected retry evidence differs')
            else:require(d['retries']==0 and d['coupled_newton']['counts_complete']is True,'unexpected retry/incomplete statistics')
            if index<4:last_age=es[-1]['age_used']
        else:
            require(p['history_unchanged']is True and'NEWTON_NOT_CONVERGED'in p['failure_reason']and
                incoming.read_bytes()==final.read_bytes()and d['step']==start and d['time']==when and d['retries']==0 and
                d['signature']!=a['signature']and d['coupled_newton']['refreshes']==1,
                'one-iteration rejection changed history or failed for another reason')
            require(all(bid['hashes'][role]==fid['hashes'][role]==origin['meta']['hashes'][role]for role in history.HISTORY_NAMES),'rejected binary histories changed')
        binaries[name]=dict(input=bid,output=fid)
    results={}
    for name,control in [('warm-a-first','cold-a'),('warm-b','cold-b'),('warm-a-return','cold-a'),('retry-a','cold-a')]:
        d=data[name];require(all(d[k]==data[control][k]for k in fields.EXACT_KEYS),'matched control configuration differs')
        result=transition.compare(states[name],states[control],d,ref,owners,deadline);lifts={}
        left = exact_inputs[name]['native']['residuals']
        right = exact_inputs[control]['native']['residuals']
        result['exact_input_cold_residual_MAX'] = max(abs(x-y) for x, y in zip(left, right))
        require(result['exact_input_cold_residual_MAX'] <= 1e-10,
                'exact-input native residual differs from a cold evaluator')
        for key,scale in [('essential_flow_lift',1.),('essential_k_lift',ref.Kin),('essential_epsilon_lift',ref.A/ref.H)]:
            v,w=np.asarray(d[key]),np.asarray(data[control][key]);require(v.shape==w.shape,'lift inventory differs');lifts[key]=float(np.max(abs(v-w))/scale)
        require(all(v<=1e-8 for v in lifts.values()),'matched essential lift differs');result['lifts']=lifts;results[name]=result
    changed=dict(K_over_Kin=0.,cutoff_over_H=0.,outer_epsilon_over_local_reference=0.,gradient_over_Uref_per_Lref=0.)
    changed.update(U_over_Uref=0., transported_k_over_Kin=0., dependent_Ew_over_local_reference=0.)
    sa, sb = (original[key]['states']['current'] for key in ('A', 'B'))
    for e in a['elements']:
        x,r,_=fields.rule(e,17,True);va=fields.fields(sa,e,a['sizes'],x,r);vb=fields.fields(sb,e,a['sizes'],x,r)
        changed['U_over_Uref'] = max(changed['U_over_Uref'], float(np.max(np.linalg.norm(va['U']-vb['U'], axis=1))))
        changed['transported_k_over_Kin'] = max(changed['transported_k_over_Kin'], float(np.max(abs(va['k']-vb['k'])))/ref.Kin)
        changed['gradient_over_Uref_per_Lref']=max(changed['gradient_over_Uref_per_Lref'],float(np.max(abs(va['gradU']-vb['gradU']))))
        if e['epsilon']:changed['outer_epsilon_over_local_reference']=max(changed['outer_epsilon_over_local_reference'],float(np.max(abs(va['epsilon']-vb['epsilon'])/(ref.A/r))))
        else:
            owner = owners[e['geometry']]
            bk = fields.po.basis(owner, 'k', x, np.full(len(x), ref.H))[0]
            ka, kb = (bk@s['k'][owner['k']] for s in (sa, sb))
            ea, eb = (K**1.5/ref.cl/np.maximum(r, ref.cutoff(K)) for K in (ka, kb))
            scale = ref.A/np.maximum(r, float(ref.cutoff(ref.Kin)))
            changed['dependent_Ew_over_local_reference'] = max(
                changed['dependent_Ew_over_local_reference'], float(np.max(abs(ea-eb)/scale)))
    seen=set()
    for e in owners.values():
        if e['geometry']in seen:continue
        seen.add(e['geometry']);x=np.linspace(e['bounds'][0],e['bounds'][1],17);r=np.full(len(x),ref.H);b=fields.po.basis(e,'k',x,r)[0]
        ka,kb=b@sa['k'][e['k']],b@sb['k'][e['k']]
        changed['K_over_Kin']=max(changed['K_over_Kin'],float(np.max(abs(ka-kb)))/ref.Kin)
        changed['cutoff_over_H']=max(changed['cutoff_over_H'],float(np.max(abs(ref.cutoff(ka)-ref.cutoff(kb))))/ref.H)
    changed['essential_velocity_over_Uref']=float(np.max(abs(np.asarray(a['essential_flow_lift'])-data['warm-b']['essential_flow_lift'])))
    changed['essential_k_over_Kin']=float(np.max(abs(np.asarray(a['essential_k_lift'])-data['warm-b']['essential_k_lift'])))/ref.Kin
    limits=dict(K_over_Kin=1e-3,cutoff_over_H=1e-5,outer_epsilon_over_local_reference=1e-3,gradient_over_Uref_per_Lref=1e-5,essential_velocity_over_Uref=1e-7,essential_k_over_Kin=1e-4)
    limits.update(U_over_Uref=1e-7, transported_k_over_Kin=1e-3, dependent_Ew_over_local_reference=1e-3)
    require(all(changed[k]>=v for k,v in limits.items()),'A/B state/boundary excursion is insufficient')
    out.update(status='PASS_COMPLETE_COUPLED_RETURN_FOUNDATION_ONLY',return_gate_pass=True,comparisons=results,current_equations=equations,
        binaries=binaries,cache_events=events,observed_A_B_changes=changed,minimum_excursion=limits,
        exact_input_equations=fixed, exact_input_comparisons=exact_inputs, stale_history_rejection=dict(
            file=rejection_identity, checkpoint=stale_meta, reason=rejection['reason']),
        scope='Actual public native Coupling A/B/A and typed-retry/history controls on source19; probe origin is explicit. No physical pilot, selected production recipe, monitoring/storage/cost or Task3 closure.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in('manifest','table','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--case',choices=ACCEPTED);p.add_argument('--seconds',type=float,required=True);a=p.parse_args();began=time.monotonic()
    require(not a.out.exists(),'immutable return report exists')
    out=dict(schema='WJH27CoupledReturnCurrentV1'if a.case else'WJH27CoupledReturnCheckV1',status='FAIL',return_gate_pass=False,
        full_H27_qualification=False,TASK3_READY=False,manifest=hr.identity(a.manifest),table=hr.identity(a.table),
        sources={k:hr.identity(v)for k,v in dict(checker=__file__,oracle=core.__file__,fields=fields.__file__,binary=history.__file__).items()})
    try:
        require(sys.platform.startswith('linux')and 0<a.seconds<=120,'bounded authoritative Linux assessment required')
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in(soft,hard)if v!=resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
        require(hr.identity(a.table)['sha256']==TABLE,'reference table differs')
        if a.case:current_case(read(a.manifest),a.manifest.resolve().parent,a.table,a.case,began+a.seconds,out)
        else:compare(read(a.manifest),a.manifest.resolve().parent,a.table,began+a.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out.update(elapsed_seconds=time.monotonic()-began,maximum_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    with a.out.open('x')as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:out[k]for k in('status','reason','elapsed_seconds')if k in out}))
    passed = out.get('current_coupled_step_gate_pass', False) if a.case else out['return_gate_pass']
    if 'exception' in out or not passed:
        raise SystemExit(2)


if __name__=='__main__':main()
