#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One selected CT5 split, paired traction gates and conditional model transfer."""
import argparse,copy,hashlib,json,os
from pathlib import Path
from types import SimpleNamespace
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted4_run as base
from ct_accepted5_traction import profiles,difference
from ct_accepted_pilot import artifact_hash,compare,replay
from ct_accepted2_audit import values
from ct_accepted_budget_gate import physical_budgets
from ct_accepted4_reference import write_json
from ct_accepted_open import require,run,probe
from rans_request_manifest import RequestManifest,runtime_identity,sha,canonical

def read(p):return json.loads(Path(p).read_text())

def physical(c):
    keys=('model','H','Z','nu','ustar','B','K0','roughness_input','provider','reference_policy','coordinate','pressure')
    return {k:c[k] for k in keys}

def traction_gate(members,c,target):
    names=list(members);require(len(names)==2,'traction_pair_coverage')
    if any(not np.isfinite(r['time']) or abs(r['time']-target)>1e-13 or r['failure'] for r in members.values()):return dict(status='not_evaluated',maximum_normalized_spread=None,reason='incomplete_or_rejected_paired_endpoint')
    ps=[members[n]['traction_history'][-1]['profile'] for n in names];d=difference(*ps,c)
    spread=max(abs(d['signed_mean_difference']),*(q['normalized_weighted_RMS'] for q in d['profiles'].values()))
    limits={};states=[]
    for key in ('normalized_wall_traction_error','normalized_interface_traction_RMS'):
        vals=[q[key] for q in ps];require(all(v is not None and np.isfinite(v) for v in vals),'nonfinite_traction_gate_primitive')
        classification='pass' if max(vals)+2*spread<=.001 else 'resolved_failure' if min(vals)-2*spread>.001 else 'ambiguous'
        limits[key]=dict(member_errors=vals,timestep_uncertainty=spread,limit=.001,upper=max(vals)+2*spread,lower=min(vals)-2*spread,status=classification);states.append(classification)
    return dict(status='pass' if all(q=='pass' for q in states) else 'ambiguous' if 'ambiguous' in states and 'resolved_failure' not in states else 'fail',
        signed_profile_difference=d,maximum_normalized_spread=spread,quantities=limits,definition='signed means and actual two-sided full traction profiles, not a difference of absolute error norms')

def basic_gate(rows,temporal,target,dts,initial_hash):
    reasons=[]
    if set(rows)!={str(x) for x in dts}:reasons.append('missing_timestep_member')
    for dt,r in rows.items():
        if r['failure'] or not np.isfinite(r['time']) or abs(r['time']-target)>1e-13:reasons.append(dt+':incomplete_or_rejected_endpoint')
        if r['initial_state_sha256']!=initial_hash:reasons.append(dt+':foreign_actual_initial_state')
        for key in ('domain','traction','reference_error','constitutive','mass'):
            if r[key]['status']!='pass':reasons.append(dt+':'+key)
        fields=r['reference_error']['fields']
        if set(fields)!={'u','w','k','outer_epsilon','Pi'}:reasons.append(dt+':field_coverage')
        for name,v in fields.items():
            if not np.isfinite(v['weighted_RMS']) or v['weighted_RMS']<0 or not np.isfinite(v['scale']) or v['scale']<=0 or v['weighted_RMS']/v['scale']>.001:reasons.append(dt+':field_primitive:'+name)
        mass=r['mass']
        if not np.isfinite(mass['outward_net']) or mass['physical_scale']<=0 or abs(mass['outward_net'])/mass['physical_scale']>1e-8:reasons.append(dt+':mass_primitive')
        source=r['constitutive'].get('maximum_scaled_defects',{})
        if set(source)!={'R','Q','P','nu_t','epsilon_production','epsilon_destruction'} or any(not np.isfinite(v) or v>2e-10 for v in source.values()):reasons.append(dt+':model_source_primitives')
        if physical_budgets(r)['status']!='pass':reasons.append(dt+':physical_budget_primitives')
        if r['corrections'] or r['additional_solves']:reasons.append(dt+':forbidden_correction')
        if r['reference_error']['fields']['w']['weighted_RMS']>1e-4:reasons.append(dt+':vertical_velocity_accuracy')
    if not temporal or abs(temporal[-1]['time']-target)>1e-13 or any(q['status']!='pass' for q in temporal):reasons.append('temporal_fields_or_coverage')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,physical_completion='pass' if len(rows)==2 and all(not r['failure'] and abs(r['time']-target)<1e-13 for r in rows.values()) else 'fail')

class Campaign:
    def __init__(self,a):
        self.a=a;self.root=a.root;self.rt=a.runtime;self.old=a.ct4;self.contract=a.contract;self.model=None;self.c=None
        self.diag=read(a.diagnosis/'results.json');require(self.diag['status']=='pass','diagnosis_not_admitted')
        self.inits={};self.members={};self.initial_hash=None;self.baseline_packet=None;self.quarter=False
    def prepare(self,model,dts=(.002,.001),quarter=False):
        self.model=model
        source=self.old/('initial-v3/'+model+'-6.0-P'+('4' if model=='A95' else '3'))/'results.json';old=read(source);c=copy.deepcopy(old['case']);before=physical(c)
        c['degree']=4;c['pressure_degree']=2;point=self.diag['selection']['point']/self.diag['selection']['case']['H']*c['H']
        require(point not in c['knots'],'selected_split_already_present');c['knots']=sorted(c['knots']+[point]);require(physical(c)==before,'continuous_input_changed_by_split')
        self.c=c;folder=self.root/('initial-quarter' if quarter else 'initial')/model
        prep=base.initial(SimpleNamespace(references=self.old/'reference-v1/results.json',model=model,height=6.,degree=4,runtime=self.rt,contract=self.contract,output=folder),case=c,timesteps=dts)
        require(prep['status']=='pass','new_space_initial_accuracy_or_domain_failure')
        for dt,row in prep['members'].items():
            require(row['slots']<=3000,'native_unknown_guard')
            if self.diag['selection']['proposed_first_slab_selected']:require(row['slots']==2957,'predicted_native_count_wrong')
            if self.initial_hash is None:self.initial_hash=row['state_sha256']
            require(row['state_sha256']==self.initial_hash,'fresh_timestep_initial_coefficients_differ')
            self.inits[dt]=(folder,prep)
        return dict(source_initial_receipt_sha256=sha(source),new_initial_receipt_sha256=sha(folder/'results.json'),initialization=prep,
            physical_input_before=before,physical_input_after=physical(c),continuous_input_sha256=hashlib.sha256(canonical(before)).hexdigest(),selected_point=point)
    def member(self,dt,target):
        name=str(dt);folder=self.root/'runs'/self.model/name/('T'+str(target));ip,prep=self.inits[name];zero=ip/name;ini=read(zero/'initial.json')
        previous=self.members.get((name,.008)) if target==.02 else None
        require(target==.008 or previous is not None,'missing_matching_parent_history')
        parent=Path(previous['directory']) if previous else zero;old=read(parent/'case.evolution.json');start=old['time']
        require(abs(start-(.008 if target==.02 else 0.))<1e-13 and old['dt']==dt,'changed_native_timestep_or_time')
        manifest=RequestManifest(folder,dict(schema='CT5NativeMemberV1',case=self.c,dt=dt,start=start,target=target,contract_sha256=sha(self.contract),
            initial_receipt_sha256=sha(ip/'results.json'),parent_report_sha256=sha(parent/'case.evolution.json'),parent_history_sha256=artifact_hash(old['checkpoint']),
            runtime=runtime_identity(self.rt/'IncNavierStokesSolver'),sources={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted5_run.py','ct_accepted4_run.py','ct_accepted5_traction.py')}),resume=folder.exists())
        if manifest.completed(folder):
            row=read(folder/'member.json');nf=Path(row['directory']);require(row['native_report_sha256']==sha(nf/'case.evolution.json') and row['history_sha256']==artifact_hash(read(nf/'case.evolution.json')['checkpoint']),'changed_completed_member');return row
        # base.make creates a new directory; the immutable manifest already
        # owns this folder, so generate in a child and keep a precise link.
        native_folder=folder/'native';base.make(native_folder,self.rt,self.c,dt,round((target-start)/dt),Path(old['checkpoint']),parent/'case.xml')
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(native_folder/'timing.jsonl');e=run([self.rt/'IncNavierStokesSolver','case.xml'],native_folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        n=read(native_folder/'case.evolution.json');require(n['accepted_states'][0]==old['state'],'native_parent_consumption_changed')
        a=SimpleNamespace(prep=prep,runtime=self.rt);row=dict(directory=str(native_folder),execution=e,identity=base.identity(self.c,native_folder,self.rt),native_report_sha256=sha(native_folder/'case.evolution.json'),
            initial_state_sha256=self.initial_hash,state_sha256=exact_vector_hash(n['state']),history_sha256=artifact_hash(n['checkpoint']),initialization_receipt=n['initialization_receipt'],
            corrections=sum(x.get('bound_calls',0) for x in n['pimple_attempts']),additional_solves=sum(x.get('additional_governing_equation_solves',0) for x in n['pimple_attempts']),
            maximum_refreshed_epsilon_lag=max((x.get('refreshed_trace_defect',0) for x in n['pimple_attempts']),default=0))
        row.update(base.endpoint(a,ini,n,self.c,native_folder));row['traction_history']=[dict(time=t,profile=profiles(ini,s,self.c)) for t,s in zip(n['audit_times'],n['accepted_states'],strict=True)]
        row['first_accuracy_crossing']=next((dict(time=q['time'],wall=q['profile']['normalized_wall_traction_error'],interface=q['profile']['normalized_interface_traction_RMS']) for q in row['traction_history'] if q['profile']['normalized_wall_traction_error']>.001 or q['profile']['normalized_interface_traction_RMS'] is None or q['profile']['normalized_interface_traction_RMS']>.001),None)
        write_json(folder/'member.json',row);manifest.seal(folder,['member.json']);return row
    def pair(self,target,dts=(.002,.001)):
        rows={}
        for dt in dts:
            row=self.member(dt,target);self.members[(str(dt),target)]=row;rows[str(dt)]=row
        temporal=[];packets=[]
        for dt,row in rows.items():
            ip,prep=self.inits[dt];packet=read(ip/dt/'initial.json');n=read(Path(row['directory'])/'case.evolution.json');packets.append((packet,n))
        common=sorted(set.intersection(*(set(round(t,12) for t in n['audit_times']) for _,n in packets)))
        for t in common:
            obs=[]
            for packet,n in packets:
                idx=next(i for i,x in enumerate(n['audit_times']) if abs(t-x)<1e-13);v=values(packet,n['accepted_states'][idx],prep['global_count']);v['observation']['time']=t;obs.append(v)
            temporal.append(compare(*obs,prep['weights'],self.c))
        basic=basic_gate(rows,temporal,target,dts,self.initial_hash);traction=traction_gate(rows,self.c,target)
        folder=self.root/'pairs'/self.model/('T'+str(target)+'-'+str(dts[-1]));m=RequestManifest(folder,dict(schema='CT5PairedMilestoneV1',case=self.c,target=target,dts=dts,
            member_receipts={dt:sha(Path(row['directory']).parent/'member.json') for dt,row in rows.items()},contract_sha256=sha(self.contract),source_sha256=sha(Path(__file__))))
        q=dict(case=self.c,target=target,dts=list(dts),members=rows,temporal=temporal,basic=basic,traction_uncertainty=traction,
            qualification=dict(status='pass' if basic['status']==traction['status']=='pass' else 'fail',physical_completion=basic['physical_completion']))
        write_json(folder/'results.json',q);m.seal(folder,['results.json']);q['receipt']=str(folder/'results.json');return q
    def spatial(self,q):
        require(self.model=='A95','spatial_control_only_primary')
        target=q['target'];require(target==.008,'no_unequal_baseline_horizon_comparison')
        oldinit=self.old/'initial-v3/A95-6.0-P4';oldrun=self.old/'evolution-v1/A95-6.0-P4/T0.008';ini=read(oldinit/'0.001/initial.json');basec=read(oldinit/'results.json')['case']
        require(physical(basec)==physical(self.c),'foreign_continuum_spatial_comparator')
        newip,prep=self.inits['0.001'];pts,sides,weights,n=base.grid(self.c);folder=self.root/'spatial'/str(q['dts'][-1]);manifest=RequestManifest(folder,dict(schema='CT5PhysicalSpatialContrastV1',case=self.c,endpoint=q['target'],dts=q['dts'],source_sha256=sha(Path(__file__))))
        shutil=__import__('shutil');shutil.copy2(oldinit/'0.001/case.xml',folder/'case.xml')
        bp,e=probe(self.rt,folder,replay(ini,ini['observation']['state'],0.,pts,grid_sides=sides,export_test_basis=True),'common-basis');require(e['exit_code']==0,'baseline_common_grid_export_failed')
        oldprofiles=[];differences={};initial=[]
        for dt in ('.002','.001'):
            name=str(float(dt));nold=read(oldrun/name/'case.evolution.json');oldprofiles.append(profiles(ini,nold['state'],basec))
            if name not in q['members']:continue
            rnew=q['members'][name];nnew=read(Path(rnew['directory'])/'case.evolution.json');npkt=read(newip/name/'initial.json')
            require(abs(nold['time']-target)<1e-13 and abs(nnew['time']-target)<1e-13,'unequal_spatial_endpoints')
            aa=values(bp,nold['state'],n);bb=values(npkt,nnew['state'],n);aa['observation']['time']=bb['observation']['time']=target;differences[name]=compare(aa,bb,weights,self.c)
            a0=values(bp,ini['observation']['state'],n);b0=values(npkt,npkt['observation']['state'],n);initial.append(compare(a0,b0,weights,self.c))
        uncertainty=difference(*oldprofiles,basec);oldspread=max(abs(uncertainty['signed_mean_difference']),*(x['normalized_weighted_RMS'] for x in uncertainty['profiles'].values()))
        newprofiles=[q['members'][str(x)]['traction_history'][-1]['profile'] for x in q['dts']];cross=difference(oldprofiles[-1],newprofiles[-1],self.c)
        improvement=abs(cross['signed_mean_difference']);limit=4*max(oldspread,q['traction_uncertainty']['maximum_normalized_spread'],2e-7)
        r=dict(status='resolved_improvement' if improvement>limit and newprofiles[-1]['normalized_wall_traction_error']<oldprofiles[-1]['normalized_wall_traction_error'] else 'unresolved',
            physical_identity=physical(self.c),common_grid_sha256=hashlib.sha256(canonical([pts,sides,weights])).hexdigest(),baseline_common_packet_sha256=sha(folder/'common-basis.json'),
            endpoint_differences=differences,initial_projection_differences=initial,baseline_timestep_traction=uncertainty,spatial_traction=cross,normalized_spatial_signal=improvement,uncertainty_limit=limit,
            historical_pressure_interpretation='CT4 P3/P4 pressure uncertainty remains unresolved; this new contrast does not relabel it',scope='same continuous functions, union physical inventory, compatible P4/P2 spaces; no convergence order')
        write_json(folder/'results.json',r);manifest.seal(folder,['results.json']);return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct4','runtime','diagnosis','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');m=RequestManifest(a.output,dict(schema='CT5FixedCampaignV1',contract_sha256=sha(a.contract),diagnosis_sha256=sha(a.diagnosis/'results.json'),
        runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__))))
    result=dict(status='incomplete',models={},scope='one split, heights5.5/6.5 deferred, no .2/curved/hill work')
    for model in ('A95','A93','E94E95'):
        if model!='A95' and result['models']['A95'].get('status')!='pass':result['models'][model]=dict(status='not_admitted_primary_failed');continue
        camp=Campaign(a);row=dict(initial=camp.prepare(model),milestones={});result['models'][model]=row;dts=(.002,.001)
        for target in (.008,.02):
            if target==.02 and row['milestones']['0.008']['qualification']['status']!='pass':row['continuation']='not_admitted';break
            q=camp.pair(target,dts);row['milestones'][str(target)]=q
            if model=='A95' and target==.008:row['spatial']=camp.spatial(q) if q['basic']['physical_completion']=='pass' else dict(status='not_evaluated_incomplete_endpoint')
            unresolved_spatial=model=='A95' and row['spatial']['status']!='resolved_improvement'
            allowed_half=all(r.endswith(':traction') or r=='temporal_fields_or_coverage' for r in q['basic']['reasons'])
            if (q['traction_uncertainty']['status']=='ambiguous' or unresolved_spatial) and allowed_half and not camp.quarter and q['basic']['physical_completion']=='pass' and all(x['domain']['status']=='pass' for x in q['members'].values()):
                camp.quarter=True;row['quarter_initial']=camp.prepare(model,(.0005,),True);dts=(.001,.0005)
                for t in (.008,.02) if target==.02 else (.008,):
                    fine=camp.pair(t,dts);row['milestones'][str(t)+'-quarter']=fine
                    if fine['qualification']['status']!='pass':break
                if fine['target']!=target:row['quarter_continuation']='not_admitted_finer_short_gate_failed';break
                q=fine;row['milestones'][str(target)]=q
                if model=='A95' and target==.008:
                    row['spatial']=camp.spatial(q)
            if q['qualification']['status']!='pass' or (model=='A95' and row['spatial']['status']!='resolved_improvement'):break
            write_json(a.output/'progress.json',result)
        row['status']='pass' if all(t in row['milestones'] and row['milestones'][t]['qualification']['status']=='pass' for t in ('0.008','0.02')) and (model!='A95' or row['spatial']['status']=='resolved_improvement') else 'failed_or_unresolved'
        write_json(a.output/'progress.json',result)
    result['status']='pass' if all(q['status']=='pass' for q in result['models'].values()) else 'failed_or_unresolved';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT5 campaign', {k:v['status'] for k,v in result['models'].items()});return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
