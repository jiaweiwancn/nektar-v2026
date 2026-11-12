#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT4 continuous-reference case generation and primitive-gated short pilots."""
import argparse,copy,hashlib,json,math,os,shutil,xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted4_reference import Column,parameters,reference,write_json,require
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info,set_parameter,write
from ct_accepted_open import run,probe,row_value
from ct_accepted_pilot import replay,compare,budget_classification,artifact_hash,constitutive_comparison
from ct_accepted2_audit import values,mass
from ct_accepted2_storage import artifact
from ct_accepted3_metrics import domain,point
from ct_accepted3_run import CONTROLS,attempt_diagnosis
from ct_accepted_budget_gate import paired_coverage,physical_budgets
from rans_request_manifest import RequestManifest,sha,runtime_identity,canonical


def construction(ref,model,mm,degree=3):
    q=ref['cases'][model+'-'+str(float(mm))];H=q['H'];hr=q['hr'];cut=q['cutoff']
    inner=[0,H/32,H/16,H/8,H/4,H/2,H,hr,cut]
    if 0<q['zv']<H:inner.append(q['zv'])
    # A pre-run log-series estimate rules out H->.25 as a first P3 cell
    # for the 1e-3 gradient/traction target. Freeze one 1.25 radial grading
    # before constructing/evolving any native case; no outcome-driven search.
    outer=[H];z=H
    while z*1.25<2.:z*=1.25;outer.append(z)
    knots=sorted(set(inner+outer+[2.]))
    return dict(model=model,H_mm=mm,H=H,Z=2.,nu=q['nu'],ustar=q['utau'],B=q['B'],K0=q['K'],
        roughness_input=hr,reference_column=q,x=[0.,.125,.25],knots=knots,degree=degree,pressure_degree=degree-2,
        provider='NeutralCompositeV1',reference_policy='Natural',coordinate='physical z, no displacement shift',pressure='absolute Pi_w=0;Pi_o=-Qzz')


def make(folder,rt,c,dt,steps,parent=None,template=None):
    if template is None:
        xml=case_xml(folder,c['model'],'Natural','Accepted',steps=steps,dt=dt,cache='Off',
            descriptor_probe=rt/'RANSNewModelsProbe',radial_knots=c['knots'],degree=c['degree'],xknots=c['x'],quadrature=[32,8,10,10,10])
        root=ET.parse(xml).getroot();conditions=root.find('CONDITIONS')
        for f in list(conditions.findall('FUNCTION')):conditions.remove(f)
        for key,value in dict(RANSWallBoundaryPolicy='AcceptedOpenFlatV1',RANSOpenTop='DirichletCompatibleV1',RANSOpenBackflow='RejectV1',
            RANSWallReference='Functions',RANSOpenReference='NeutralCompositeV1',RANSWallFunctionCoordinates='Physical',RANSWallPimpleSnapshots='On').items():set_info(root,key,value)
        for key,value in dict(RANSWallHeight=c['H'],RANSWallSampleHeight=c['H'],RANSWallTailStart=c['H'],RANSWallRoughness=c['roughness_input'],
            Kinvis=c['nu'],RANSWallFrictionVelocity=c['ustar'],RANSOpenReferenceLength=1.,RANSOpenMeanIntercept=c['B'],
            RANSOpenEquilibriumQ=parameters(c['model'])['q'],RANSOpenNaturalKappa=parameters(c['model'])['kappa'],
            K0=c['K0'],RANSKappa=parameters(c['model'])['kappa'],RANSOpenVelocityScale=1.).items():set_parameter(root,key,value)
    else:folder.mkdir(parents=True,exist_ok=False);root=ET.parse(template).getroot();xml=folder/'case.xml'
    set_parameter(root,'NumSteps',steps);set_parameter(root,'TimeStep',dt)
    set_info(root,'RANSWallRestart',str(parent) if parent else 'None');set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'))
    set_info(root,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'));write(root,xml)
    write_json(folder/'construction.json',dict(c,dt=dt,segment_steps=steps,configuration_sha256=sha(xml)))
    return xml


def identity(c,folder,rt):
    return dict(model=c['model'],wall='AcceptedBudgetAWF',policy='Natural',provider=c['provider'],H=c['H'],h_s=c['roughness_input'],
        physical_parameters_sha256=hashlib.sha256(canonical(c)).hexdigest(),configuration_sha256=sha(folder/'case.xml'),runtime=runtime_identity(rt/'IncNavierStokesSolver'))


def grid(c):
    g,w=leggauss(6);points=[];weights=[]
    for a,b in zip(c['x'][:-1],c['x'][1:]):
        for d,e in zip(c['knots'][:-1],c['knots'][1:]):
            for i in range(6):
                for j in range(6):points.append([(a+b)/2+(b-a)*g[i]/2,0.,(d+e)/2+(e-d)*g[j]/2]);weights.append((b-a)*(e-d)*w[i]*w[j]/4)
    n=len(points);sides=[0]*n
    xs=sorted(set([a+(b-a)*j/p for a,b in zip(c['x'][:-1],c['x'][1:]) for p in (3,4) for j in range(p+1)]))
    for x in xs:
        for z in (0.,c['reference_column']['cutoff']*(1-1e-8),c['reference_column']['cutoff']*(1+1e-8),c['roughness_input'],c['H']/2,
                  c['H']*(1-1e-8),c['H']*(1+1e-8),1.,2.):points.append([x,0.,z]);sides.append(0)
        for side in (-1,1):points.append([x,0.,c['H']]);sides.append(side)
    return points,sides,weights,n


def expected(packet,c):
    o=packet['observation'];col=Column(c['model'],c['K0'],c['H'],c['roughness_input'],c['nu'],order=96)
    at=lru_cache(None)(lambda z,outer:reference(col,z,outer))
    vs=[];gs=[]
    for p,r in zip(o['grid'],o['grid_regions'],strict=True):
        v=at(p[2],bool(r));vs.append([v[0],0,v[1],v[2] if r else 0,v[3]])
        gs.append([[0,0,v[4]],[0,0,0],[0,0,v[5]],[0,0,v[6] if r else 0],[0,0,0]])
    return dict(observation=dict(time=o['time'],grid=o['grid'],grid_sides=o['grid_sides'],grid_regions=o['grid_regions'],grid_values=vs,grid_gradients=gs)),at


def reference_agreement(packet,c,at):
    defects=np.zeros(9);k_nonconstant=0.;skipped=0
    for row in packet['neutral_composite_reference']['samples']:
        z=row['coordinate'][2];r=bool(row['region']);a=np.array(row['values_U_k_epsilon_Pi_gradients_Q0_D']);b=at(z,r)
        delta=abs(a-b);scales=np.array([1,c['K0'],c['reference_column']['counterfactual']['epsilon_amplitude']/c['H'],1,1,c['K0']/c['H'],c['reference_column']['counterfactual']['epsilon_amplitude']/c['H']**2,c['reference_column']['counterfactual']['epsilon_amplitude'],1])
        if not r and abs(z-c['reference_column']['cutoff'])<4e-13*c['H']:delta[6]=0;skipped+=1
        defects=np.maximum(defects,delta/scales)
        if not r:k_nonconstant=max(k_nonconstant,abs(b[1]-c['K0'])/c['K0'])
    return dict(status='pass' if np.max(defects)<=2e-9 else 'fail',maximum_scaled_defects=defects.tolist(),
        cutoff_derivative_boundary_exemptions=skipped,nonconstant_inner_k_signal=k_nonconstant,
        scope='independent finite-nu split quadrature and original stress oracle vs same continuous provider consumed by driver/live data/probe')


def traction(initial,state,c):
    o=initial['observation'];s=np.array(state);off=o['offsets'];H=c['H'];p=c['degree'];col=parameters(c['model'])
    # Exact polynomial restriction of the current shared endpoint trace. This
    # is a field evaluation, not reuse of constitutive validity or fixed K.
    traces=[]
    for a,b in zip(c['x'][:-1],c['x'][1:]):
        xs=[a+(b-a)*j/p for j in range(p+1)];ks=[]
        for x in xs:
            i=next(i for i,(q,side) in enumerate(zip(o['grid'],o['grid_sides'],strict=True)) if q==[x,0.,H] and side==1)
            ks.append(row_value(initial['test_basis'][i]['velocity_k'],s,off[2])[0])
        traces.append((a,b,np.polynomial.polynomial.polyfit((np.array(xs)-a)/(b-a),ks,p)))
    def nt(x,z):
        a,b,coef=next(t for t in traces if t[0]<=x<=t[1]);K=float(np.polynomial.polynomial.polyval((x-a)/(b-a),coef))
        require(K>0,'invalid_actual_endpoint_in_traction')
        q=c['roughness_input']*math.sqrt(K)/(70*c['nu']);B=10.7 if q==0 else 10.7*(1-q**max(.5-.4*q**.7,1-.79*q**(-.28)))
        return col['alpha']*math.sqrt(K)*max(0.,z-c['nu']*B/math.sqrt(K))
    wall=np.zeros(3);wall_weight=0.;outlet=0.;outlet_weight=0.;interfaces={};rows=[]
    for f in o['faces']:
        x,_,z=f['coordinate'];region=f['region'];G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):G[axis]=row_value(f['velocity_k'],s,off[b])[1]
        pi=row_value(f['pressure'],s,off[4+region])[0]
        if region:
            v=point(f,s,off,c['model'])
            if 'total_stress' not in v:return dict(status='undefined_actual_domain',witness=v)
            T=np.array(v['total_stress'])
        else:T=(c['nu']+nt(x,z))*(G+G.T)-np.eye(3)*pi
        if f['interface']:
            key=(round(x,14),round(z,14));interfaces.setdefault(key,[]).append((T[:,2],f['weight'],f['velocity_k']['geometry']))
        elif z==0:
            wall+=f['weight']*T[:,2];wall_weight+=f['weight']
        elif x==c['x'][-1]:
            given=np.array([c['reference_column']['outer_equilibrium']['outlet_traction'][0] if region else 0.,0.,c['ustar']**2])
            err=T[:,0]-given;outlet+=f['weight']*float(err@err);outlet_weight+=f['weight']
    jump=0.;iw=0.
    for key,ts in interfaces.items():
        require(len(ts)==2,'one_sided_interface_traction_pair_missing')
        d=ts[0][0]-ts[1][0];jump+=ts[0][1]*float(d@d);iw+=ts[0][1]
        rows.append(dict(coordinate=key,traction_jump=d.tolist(),geometries=[ts[0][2],ts[1][2]]))
    t=c['ustar']**2;wall_error=abs(wall[0]/wall_weight-t)/t;interface=math.sqrt(jump/iw)/t
    return dict(status='pass' if max(wall_error,interface)<=1e-3 else 'fail',mean_wall_shear=wall[0]/wall_weight,
        normalized_wall_traction_error=wall_error,normalized_interface_traction_RMS=interface,
        normalized_outlet_reconstruction_RMS=math.sqrt(outlet/outlet_weight)/t,interface_rows=rows,
        scope='actual native one-sided gradients and current polynomial endpoint K; outlet datum imposed weakly, reconstructed outlet mismatch reported separately')


def initial(a):
    ref=json.loads(a.references.read_text());c=construction(ref,a.model,a.height,a.degree);pts,sides,weights,n=grid(c)
    m=RequestManifest(a.output,dict(schema='CT4FreshInitialV1',contract_sha256=sha(a.contract),reference_sha256=sha(a.references),case=c,
        runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__))))
    result=dict(schema='CT4InitialV1',case=c,global_count=n,weights=weights,members={},status='fail')
    for dt in (.002,.001):
        folder=a.output/str(dt);make(folder,a.runtime,c,dt,0)
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'zero.log');row=dict(status='fail',execution=e,directory=str(folder));result['members'][str(dt)]=row
        if not (folder/'case.evolution.json').exists():row['reason']='native_zero_initialization_failed';continue
        native=json.loads((folder/'case.evolution.json').read_text())
        if 'state' not in native:
            row.update(reason='native_initialization_exception',native_failure=native,native_report_sha256=sha(folder/'case.evolution.json'))
            write_json(a.output/'progress.json',result);continue
        row.update(slots=len(native['state']),native_report_sha256=sha(folder/'case.evolution.json'),failure=native['failure'])
        packet,pe=probe(a.runtime,folder,dict(operation='functions',grid=pts,grid_sides=sides,export_basis=True,export_test_basis=True,
            export_constitutive=True,export_neutral_composite_reference=True),'initial');row['probe_execution']=pe
        if e['exit_code'] or pe['exit_code']:row['reason']='native_or_independent_initial_probe_rejection';continue
        require(packet['observation']['state']==native['state'],'driver_and_probe_effective_initial_differ')
        require(len(native['state'])<=3000,'native_unknown_guard')
        truth,at=expected(packet,c);small=values(packet,native['state'],n);truth['observation']={k:v[:n] if isinstance(v,list) else v for k,v in truth['observation'].items()}
        err=compare(small,truth,weights,c);d=domain(packet,native,c);tr=traction(packet,native['state'],c)
        row.update(domain=d,traction=tr,reference_error=err,continuous_oracle=reference_agreement(packet,c,at),
            constitutive=constitutive_comparison(packet,c),block_counts=packet['observation']['sizes'],state_sha256=exact_vector_hash(native['state']),
            initial_history_sha256=artifact_hash(native['checkpoint']),initialization_receipt=native['initialization_receipt'],
            raw_unforced_alpha_zero_norms=native['raw_reference_norms'],actual_quadrature=packet['observation']['actual_quadrature'],
            alpha_zero_native=packet['observation']['audit_alpha_zero'],identity=identity(c,folder,a.runtime),packet_sha256=sha(folder/'initial.json'))
        row['status']='pass' if all(x['status']=='pass' for x in (d,tr,err,row['continuous_oracle'],row['constitutive'])) else 'fail'
        write_json(a.output/'progress.json',result)
    rows=list(result['members'].values());result['identical_pair_initial_coefficients']=len(rows)==2 and all('state_sha256' in q for q in rows) and rows[0]['state_sha256']==rows[1]['state_sha256']
    result['status']='pass' if result['identical_pair_initial_coefficients'] and all(q['status']=='pass' for q in rows) else 'fail'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result


def endpoint(a,initial,native,c,folder):
    n=a.prep['global_count'];weights=a.prep['weights'];result=dict(budgets={},time=native['time'],failure=native['failure'])
    result.update(domain=domain(initial,native,c),mass=mass(initial,native['state'],c),traction=traction(initial,native['state'],c))
    obs=values(initial,native['state'],n);base=values(initial,initial['observation']['state'],n);truth,_=expected(initial,c)
    truth['observation']={k:v[:n] if isinstance(v,list) else v for k,v in truth['observation'].items()}
    for q in (obs,base,truth):q['observation']['time']=native['time']
    result['reference_error']=compare(obs,truth,weights,c);result['drift_from_actual_initial']=compare(obs,base,weights,c)
    request=replay(initial,native['state'],native['time'],initial['observation']['grid'],export_constitutive=True);request['grid_sides']=initial['observation']['grid_sides']
    packet,e=probe(a.runtime,folder,request,'endpoint');result['constitutive_execution']=e
    result['constitutive']=constitutive_comparison(packet,c) if e['exit_code']==0 else dict(status='fail',reason='actual_domain_probe_rejected')
    scheduled=[x for x in native['pimple_attempts'] if x['scheduled_complete']]
    if native['steps'] and scheduled:
        attempt=scheduled[-1]
        snap=json.loads(artifact(attempt['k_snapshot'],attempt['k_snapshot_sha256']).read_text());idx=next(i for i,t in enumerate(native['audit_times']) if abs(t-native['time'])<1e-13)
        if idx:
            previous=native['accepted_states'][idx-1]
            for label,orders,face in CONTROLS:
                sub=folder/label;sub.mkdir();root=ET.parse(folder/'case.xml').getroot()
                for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
                write(root,sub/'case.xml');req=replay(initial,native['state'],native['time'],[],export_physical_scalar_budget=True,
                    previous_state=previous,stage_history=snap['history'],export_frozen_scalar_reconciliation=True,physical_face_order=face)
                b,e=probe(a.runtime,sub,req,'budget')
                result['budgets'][label]=dict(budget_classification(b,c),execution=e,packet_sha256=sha(sub/'budget.json')) if e['exit_code']==0 else dict(status='failed_probe',execution=e,raw_packet_sha256=sha(sub/'budget.json'))
    result['budget_gate']=physical_budgets(result);result['alpha_zero_native']=packet.get('observation',{}).get('audit_alpha_zero')
    result['stage_diagnostics']=attempt_diagnosis(initial,native,c)
    result['all_accepted_native_audits']=audited_extrema(native)
    return result


def qualify(q):
    rows=q['members'];expected={dt:r['identity'] for dt,r in rows.items()};cov=paired_coverage(rows,q['target'],expected);reasons=list(cov['reasons'])
    for dt,r in rows.items():
        for key in ('domain','traction','reference_error','constitutive','mass'):
            if r[key]['status']!='pass':reasons.append(dt+':'+key)
        if r['corrections'] or r['additional_solves']:reasons.append(dt+':forbidden_correction')
        if r['reference_error']['fields']['w']['weighted_RMS']>1e-4:reasons.append(dt+':vertical_velocity_accuracy')
    if not q['temporal'] or not any(abs(r['time']-q['target'])<1e-13 for r in q['temporal']):reasons.append('missing_paired_endpoint')
    if any(r['status']!='pass' for r in q['temporal']):reasons.append('temporal_accuracy')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,coverage=cov,physical_completion='pass' if all(not r['failure'] and abs(r['time']-q['target'])<1e-13 for r in rows.values()) and len(rows)==2 else 'fail')


def advance(a):
    prep=json.loads((a.initial/'results.json').read_text());a.prep=prep;c=prep['case'];require(prep['status']=='pass','initialization_not_admitted')
    seal=json.loads((a.initial/'completed-artifacts.json').read_text());require(all(sha(a.initial/f)==h for f,h in seal['files'].items()),'initialization_changed')
    start=0.;parent=None
    if a.target==.02:
        require(a.parent is not None,'missing_parent');parent=json.loads((a.parent/'results.json').read_text())
        require(parent['qualification']==qualify(parent) and parent['qualification']['status']=='pass' and parent['target']==.008,'parent_gate_not_admitted')
        start=.008
    m=RequestManifest(a.output,dict(schema='CT4MatchedSegmentV1',case=c,start=start,target=a.target,dt=[.002,.001],initial_sha256=sha(a.initial/'results.json'),
        parent_sha256=sha(a.parent/'results.json') if a.parent else None,runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__))))
    result=dict(schema='CT4SegmentV1',case=c,start=start,target=a.target,members={},temporal=[],diagnostic_execution='incomplete');endpoints={}
    for dt in (.002,.001):
        name=str(dt);zero=a.initial/name;ini=json.loads((zero/'initial.json').read_text());source=Path(parent['members'][name]['directory']) if parent else zero
        old=json.loads((source/'case.evolution.json').read_text());require(abs(old['time']-start)<1e-13 and old['dt']==dt,'history_timestep_or_time_mismatch')
        if parent:require(sha(source/'case.evolution.json')==parent['members'][name]['native_report_sha256'] and artifact_hash(old['checkpoint'])==parent['members'][name]['history_sha256'],'parent_state_or_history_changed')
        folder=a.output/name;make(folder,a.runtime,c,dt,round((a.target-start)/dt),Path(old['checkpoint']),source/'case.xml')
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        native=json.loads((folder/'case.evolution.json').read_text());require(native['accepted_states'][0]==old['state'],'consumed_parent_differs')
        row=dict(execution=e,directory=str(folder),identity=identity(c,folder,a.runtime),native_report_sha256=sha(folder/'case.evolution.json'),
            state_sha256=exact_vector_hash(native['state']),initial_state_sha256=prep['members'][name]['state_sha256'],history_sha256=artifact_hash(native['checkpoint']),
            history_times=native['native_history_times'],initialization_receipt=native['initialization_receipt'],
            corrections=sum(r.get('bound_calls',0) for r in native['pimple_attempts']),additional_solves=sum(r.get('additional_governing_equation_solves',0) for r in native['pimple_attempts']),
            maximum_refreshed_epsilon_lag=max((r.get('refreshed_trace_defect',0) for r in native['pimple_attempts']),default=0))
        row.update(endpoint(a,ini,native,c,folder));result['members'][name]=row;endpoints[name]=(ini,native);write_json(a.output/'progress.json',result)
    common=sorted(set(round(t,12) for t in endpoints['0.002'][1]['audit_times'])&set(round(t,12) for t in endpoints['0.001'][1]['audit_times']))
    for t in common:
        ps=[]
        for ini,native in endpoints.values():
            idx=next(i for i,v in enumerate(native['audit_times']) if abs(v-t)<1e-13);packet=values(ini,native['accepted_states'][idx],prep['global_count']);packet['observation']['time']=t;ps.append(packet)
        result['temporal'].append(compare(*ps,prep['weights'],c))
    result['diagnostic_execution']='pass';result['qualification']=qualify(result)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['initial','advance'],required=True)
    for k in ('runtime','references','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--model',choices=['A93','A95','E94E95']);p.add_argument('--height',type=float,choices=[5.5,6.,6.5]);p.add_argument('--degree',type=int,choices=[3,4],default=3)
    p.add_argument('--initial',type=Path);p.add_argument('--parent',type=Path);p.add_argument('--target',type=float,choices=[.008,.02])
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');r=initial(a) if a.mode=='initial' else advance(a)
    status=r.get('qualification',r)['status'];print('CT4',a.model,a.height,a.mode,a.target,status);return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
