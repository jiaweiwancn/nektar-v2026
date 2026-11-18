#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent cache/history checks across an actual moving-cutoff owner edge."""
import argparse,json,math,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import coupled_oracle as core
import coupled_mms_check as mms
import history_recovery_check as hr
import physical_cache_history_check as history
import implicit_algorithm_check as fields
import physical_transient_oracle as transient
require,read = core.require,history.read
TABLE = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'


def states(data):
    h=data['input_history'];sizes=data['sizes']
    raw=dict(accepted=data['coefficients'],input_current=h['current'],input_previous=h['previous'],
        current_advection=data['current_advection'],previous_advection=data['previous_advection'],
        input_advection=h['advection'],input_previous_advection=h['previous_advection'])
    out={key:transient.state(value,sizes) for key,value in raw.items()}
    require(transient.state_hash(out['accepted'])==data['state_sha256'] and
        transient.state_hash(out['input_current'])==data['previous_sha256'] and
        transient.state_hash(out['current_advection'])==data['advection_sha256'] and
        transient.state_hash(out['previous_advection'])==transient.state_hash(out['input_advection'])==data['previous_advection_sha256'],
        'scientific frame histories/hash shift differ')
    return out


def compare(a,b,data,ref,owners,deadline):
    result=dict(coefficients={},owners={},dependent={})
    for key in fields.STATE_NAMES+fields.CACHE_NAMES:
        delta=fields.coefficient_delta(a[key],b[key],data['sizes'],ref)
        require(all(v['scaled_MAX']<=1e-8 for v in delta.values()),'complete history coefficient limit exceeded')
        result['coefficients'][key]=delta
    for key in fields.STATE_NAMES:
        values=[fields.physical_delta(a[key],b[key],data,ref,q,grid,deadline)for q,grid in((8,9),(12,17))]
        require(all(v.get('status','MET')=='MET'for row in values for v in row['norms'].values()),'owner field limit exceeded')
        result['owners'][key]=values;delta=fields.dependent_delta(a[key],b[key],data,owners,ref,deadline)
        require(delta['gate_pass'],'dependent K/Gamma/cutoff/Ew limit exceeded');result['dependent'][key]=delta
    duals=fields.cache_duals(a,b,data,ref,deadline)
    require(all(v<=1e-8 for row in duals.values()for v in row.values()),'explicit all-row mass-dual limit exceeded')
    result['explicit_mass_duals']=duals
    return result


def run(manifest,folder,table,deadline,out):
    require(manifest['schema']=='WJH27CutoffOwnerTransitionV1','wrong transition manifest');root=hr.resolve(folder,manifest['root'])
    require(hr.identity(table)['sha256']==TABLE,'reference changed');ref=history.po.Reference(read(table));paths={}
    for name in ('fresh','cached','prefix','cold','reload15','reload18'):
        paths[name],_=history.input_path(manifest['runs'][name],root)
    loaded={name:fields.load(paths[name],'CoupledJFNKV1')for name in('fresh','cached','prefix','cold')}
    cached=loaded['cached'];d=cached['data'];require(d['source_identity']==manifest['source_identity']and
        d['manufactured']is True and d['mms_variant']=='Historical'and d['mms_body_amplitude']==1 and
        d['scalar_quadrature_policy']=='CubicCutoffFacePanelsV1'and d['step']==20 and d['dt']==.02 and d['time']==.4,
        'transition source/science/clock differs')
    owners,_=fields.inventory(d,ref);target=mms.IndependentMMS(read(table),1,'Historical');oracle=core.Oracle(d,target,20,24,deadline)
    gp,_=history.input_path(manifest['geometry'],root);face=read(gp)['cutoff_face'];frames={};trajectories={};native_costs={}
    for name,item in loaded.items():
        data=item['data'];begin=10 if name=='cold'else 0;end=10 if name=='prefix'else 20
        for key in fields.EXACT_KEYS:
            if key not in ('step','time'):require(data[key]==d[key],'transition configuration differs: '+key)
        require(data['step']==end and data['time']==end*.02,'transition endpoint clock differs')
        policy='RebuildEveryNewtonV1'if name=='fresh'else'CurrentDataDriftV1'
        packet=read(paths[name]);native_costs[name]={key:packet[key]for key in('construction_seconds','initialisation_seconds','complete_advancing_seconds','complete_step_seconds','final_export_seconds','recovery_seconds')}
        require(len(packet['complete_step_seconds'])==end-begin,'wrong executed-step count')
        index_path=history.replace_suffix(paths[name],'.result.json','.science-index.json');index=read(index_path)
        require(index['run_complete']is True and index['initial_step']==begin and index['final_step']==end and
            index['coefficients_included']is True and [v['step']for v in index['frames']]==list(range(begin+1,end+1)),
            'scientific output interval incomplete')
        for key in('source_identity','reference_identity','signature','dt'):require(index[key]==data[key],'scientific identity differs')
        rows=[hr.native_json(line)for line in history.replace_suffix(paths[name],'.result.json','.steps.jsonl').read_text().splitlines()]
        require([v['step']for v in rows]==list(range(begin+1,end+1)),'accepted log interval differs')
        require(fields.cache_log(paths[name],data,policy,deadline)['complete_logged_interval'],'cache log incomplete')
        frames[name]={};trajectory=[]
        for entry,row in zip(index['frames'],rows):
            fields.clock(deadline);frame_path,ident=history.input_path(entry,root);record=read(frame_path);step=row['step'];when=step*.02
            require(entry['origin']=='current_process_accepted_step'and entry['reason_mask']&1 and
                record['step']==step and record['time']==when and all(record[key]==row[key]for key in history.HASH_NAMES),
                'scientific frame clock/origin/history differs')
            require(all(record[key]==data[key]for key in('source_identity','reference_identity','signature','dt')),'frame scientific identity differs')
            require(all(math.isfinite(v)and 0<=v<=3e-12 for v in record['residuals'])and
                all(math.isfinite(record[key])and record[key]>0 for key in('minimum_k','minimum_epsilon','minimum_covariance_over_k','minimum_principal_over_mu')),
                'accepted native residual/admissibility failed')
            stats=fields.cache_policy_statistics(record,policy)
            require(stats['preconditioner_events']and all(event['alpha']==(1 if step==1 else 1.5)/.02 for event in stats['preconditioner_events']),
                'cache event used wrong BE/BDF2 alpha')
            value=states(record);frames[name][step]=dict(record=record,states=value,file=ident);roots=[]
            for e in d['elements']:
                xl,xr,rl,rr=e['bounds'];x,r,w=oracle.rule(e,value['accepted'],when);area=(xr-xl)*(rr-rl)
                require(np.all(w>0)and abs(float(np.sum(w))-area)<=1e-12*area,'current positive quadrature area does not close')
                if min(abs(rl-face),abs(rr-face))>1e-12:continue
                for point in oracle.cutoff_x_splits(e,value['accepted'],None,(xl,xr,face,face)):
                    error=abs(float(oracle.cutoff(oracle.trace(e,'k',np.asarray([point]),value['accepted']))[0])-face)
                    if error<=1e-12*ref.H:roots.append(dict(owner=e['geometry'],x=point,bounds=e['bounds'],cutoff_residual=error))
            require(len(roots)==2 and abs(roots[0]['x']-roots[1]['x'])<=1e-10,'current root does not retain two adjacent original owners')
            column=tuple(roots[0]['bounds'][:2]);require(all(tuple(v['bounds'][:2])==column for v in roots),'crossing radial owners disagree on column')
            columns={v['geometry']:v for v in oracle.owner.values()};values=[]
            for point in(.25,.75):
                candidates=[v for v in columns.values()if v['bounds'][0]<point<v['bounds'][1]]
                require(len(candidates)==1,'interface sample has no unique original owner')
                values.append(float(oracle.trace(candidates[0],'k',np.asarray([point]),value['accepted'])[0]))
            K=np.asarray(values);cut=oracle.cutoff(K)
            trajectory.append(dict(step=step,time=when,column=column,roots=roots,K_samples=K.tolist(),cutoff_samples=cut.tolist(),refreshes=stats['refreshes'],reuses=stats['reuses'],
                maximum_age=stats['maximum_used_age'],newton=record['nonlinear_iterations'],krylov=stats['krylov_iterations']))
        trajectories[name]=trajectory
    require(trajectories['cached'][0]['roots'][0]['x']>.55 and trajectories['cached'][-1]['roots'][0]['x']<.45,'actual current crossing did not traverse x=.5')
    changes=[(a,b)for a,b in zip(trajectories['cached'],trajectories['cached'][1:])if a['column']!=b['column']]
    require(changes,'no actual crossing owner transition observed');event_steps=sorted({v['step']for pair in changes for v in pair})
    require(any(v['reuses']>0 for v in trajectories['cached']if any(abs(v['step']-step)<=1 for step in event_steps)),
        'no inverse reuse near the actual owner transition')
    Kchange=float(np.max(abs(np.asarray(trajectories['cached'][-1]['K_samples'])-trajectories['cached'][0]['K_samples'])))/ref.Kin
    cchange=float(np.max(abs(np.asarray(trajectories['cached'][-1]['cutoff_samples'])-trajectories['cached'][0]['cutoff_samples'])))/ref.H
    require(Kchange>=1e-3 and cchange>=1e-5,'trajectory has insufficient current K/cutoff change')
    frame_comparisons={};full_comparisons={}
    for name,steps in [('fresh',range(1,21)),('cold',range(11,21))]:
        maxima={}
        for step in steps:
            a,b=frames['cached'][step]['states'],frames[name][step]['states'];summary={}
            for key in fields.STATE_NAMES+fields.CACHE_NAMES:
                delta=fields.coefficient_delta(a[key],b[key],d['sizes'],ref);require(all(v['scaled_MAX']<=1e-8 for v in delta.values()),'intermediate frame/cache history disagreement')
                summary[key]=delta
            maxima[str(step)]=summary
        frame_comparisons[name]=maxima
        for step in sorted(set(event_steps+[20])&set(steps)):
            full_comparisons[name+':'+str(step)]=compare(frames['cached'][step]['states'],frames[name][step]['states'],d,ref,owners,deadline)
    for step in range(1,11):
        require(all(frames['prefix'][step]['record'][key]==frames['cached'][step]['record'][key]for key in history.HASH_NAMES),'prefix is not the uninterrupted accepted history')
    incoming,input_id=history.input_path(manifest['cold_input'],root);require(incoming==history.replace_suffix(paths['prefix'],'.result.json','.final.bin'),'wrong prefix checkpoint')
    binary,pid=history.binary_history(incoming,loaded['prefix']['data'],.02,deadline,input_id['sha256'])
    for role,key in zip(history.HISTORY_NAMES,history.HASH_NAMES):require(pid['hashes'][role]==frames['cached'][10]['record'][key],'cold binary history differs from uninterrupted step10')
    first=frames['cold'][11]['record'];event=first['coupled_newton']['preconditioner_events'][0]
    require(event['reason']=='cold'and event['refreshed']is True and event['age_before']==0 and
        first['previous_sha256']==pid['hashes']['current']and first['previous_advection_sha256']==pid['hashes']['advection'],'cold first event/history differs')
    cp=read(paths['cached']);require(cp['rotated_recovery_count']==4 and len(cp['recovery_seconds'])==6,'six writes/four actual rotations not observed')
    prefix=history.replace_suffix(paths['cached'],'.result.json','');retained={}
    expected=[Path(str(prefix)+'.recovery-'+str(step)+'.bin')for step in(15,18)]
    require([Path(p)for p in cp['retained_recovery']]==expected and all(not Path(str(prefix)+'.recovery-'+str(step)+'.bin').exists()for step in(3,6,9,12)),'retained/removed recovery inventory differs')
    for step,path in zip((15,18),expected):
        source,source_id=history.input_path(manifest['reload_inputs'][str(step)],root);require(source==path,'wrong retained prelaunch input')
        record=frames['cached'][step]['record'];_,bid=history.binary_history(path,record,.02,deadline,source_id['sha256'])
        require(all(bid['hashes'][role]==record[key]for role,key in zip(history.HISTORY_NAMES,history.HASH_NAMES)),'retained binary/science binding differs')
        zp=paths['reload'+str(step)];packet=read(zp);z=packet['state'];final=history.replace_suffix(zp,'.result.json','.final.bin')
        require(packet['restart']is True and z['step']==step and z['time']==step*.02 and final.read_bytes()==path.read_bytes()and
            packet['complete_step_seconds']==[]and packet['complete_advancing_seconds']==0 and z['rebuilds']==0 and
            z['coupled_newton']['refreshes']==z['coupled_newton']['reuses']==0 and not history.replace_suffix(zp,'.result.json','.steps.jsonl').read_bytes(),
            'retained zero import changed binary histories or advanced')
        retained[str(step)]=dict(input=bid,output=hr.identity(final),result=hr.identity(zp),exact_four_histories=True,advancing_evidence_claimed=False)
    equations={}
    for name in('fresh','cached','cold'):
        rp,rid=history.input_path(manifest['current_equations'][name],root);report=read(rp)
        require(report['native_sha256']==hr.identity(paths[name])['sha256']and report['table_sha256']==TABLE and
            report['native_source_identity']==manifest['source_identity']and report['current_coupled_step_gate_pass']is True,
            'separate current MMS equation gate failed')
        equations[name]=dict(file=rid,status=report['status'],residuals=report['independent_assemblies'][-1]['current_residual'])
    compare_path,compare_id=history.input_path(manifest['cache_comparison'],root);comparison=read(compare_path)
    require(comparison['comparison_gate_pass']is True and comparison['source_identity']==manifest['source_identity'],'complete fresh/cache endpoint comparison failed')
    out.update(status='PASS_CUTOFF_OWNER_CACHE_HISTORY_ONLY',foundation_gate_pass=True,source_identity=manifest['source_identity'],
        trajectories=trajectories,owner_transitions=[dict(before=a,after=b)for a,b in changes],sampled_K_change_over_Kin=Kchange,
        sampled_cutoff_change_over_H=cchange,frame_coefficient_comparisons=frame_comparisons,transition_and_endpoint_comparisons=full_comparisons,
        cold_input=pid,first_cold_inverse_event=event,retained_recovery_imports=retained,observed_writes=[3,6,9,12,15,18],observed_removed=[3,6,9,12],removed_binary_byte_audit_claimed=False,
        independent_equations=equations,complete_cache_comparison=compare_id,native_costs=native_costs,
        scope='Actual twenty-step moving-cutoff MMS, inverse-cache and cold-history/retention foundation. No physical pilot, selected recipe, monitoring/storage/cost or Task3 closure.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in('manifest','table','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True);a=p.parse_args();began=time.monotonic();require(not a.out.exists(),'transition assessment is immutable')
    out=dict(schema='WJH27CutoffOwnerTransitionCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,manifest=hr.identity(a.manifest),
        sources={key:hr.identity(path)for key,path in dict(checker=__file__,fields=fields.__file__,oracle=core.__file__,history=history.__file__,continuous_mms=mms.__file__).items()},workers=1)
    try:
        require(sys.platform.startswith('linux')and math.isfinite(a.seconds)and 0<a.seconds<=120,'bounded Linux assessment required')
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in(soft,hard)if v!=resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
        run(read(a.manifest),a.manifest.resolve().parent,a.table.resolve(),began+a.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds']=time.monotonic()-began;out['maximum_resident_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with a.out.open('x')as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({key:out[key]for key in('status','reason','elapsed_seconds')if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__=='__main__':main()
