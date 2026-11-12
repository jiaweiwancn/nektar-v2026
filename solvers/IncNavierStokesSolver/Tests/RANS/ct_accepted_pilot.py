#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded, unforced CT rough-floor pilot. A stopped case is never repaired."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted_open import pilot_case,probe,run,require
from rans_request_manifest import RequestManifest,runtime_identity,sha
from run_cls_epawf_picard import set_parameter,write
from validate_a93_a95_e94e95 import oracle
import hashlib


def artifact_hash(path):
    path=Path(path)
    if path.is_file():return sha(path)
    require(path.is_dir(),'missing_native_history_artifact')
    items={str(p.relative_to(path)):sha(p) for p in sorted(path.rglob('*')) if p.is_file()}
    require(items,'empty_native_history_artifact')
    return hashlib.sha256(json.dumps(items,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def grid(construction):
    x,w=leggauss(6);points=[];weights=[]
    for a,b in zip(construction['x'][:-1],construction['x'][1:]):
        for c,d in zip(construction['knots'][:-1],construction['knots'][1:]):
            for i in range(6):
                for j in range(6):
                    points.append([(a+b)/2+(b-a)*x[i]/2,0.,(c+d)/2+(d-c)*x[j]/2]);weights.append((b-a)*(d-c)*w[i]*w[j]/4)
    return points,weights


def replay(initial,state,time,points,**extra):
    o=initial['observation']
    return dict(operation='replay',identity=o['identity'],descriptor=o['resolved_descriptor'],
        state=state,reference=initial['reference'],fixed_boundary=initial['fixed_boundary'],time=time,
        grid=points,grid_sides=[0]*len(points),**extra)


def physical_observations(packet,weights,construction):
    o=packet['observation'];v=np.array(o['grid_values']);g=np.array(o['grid_gradients']);w=np.array(weights);outer=np.array(o['grid_regions'])==1
    divergence=g[:,0,0]+g[:,1,2];faces=o['faces'];mass=0.;by_face={}
    off=np.array(o['offsets'],int);s=np.array(o['state'])
    for f in faces:
        if f['interface']:continue
        q=f['velocity_k'];ids=np.array(q['ids'],int);phi=np.array(q['phi']);u=np.array([phi@s[ids+off[b]] for b in (0,1)])
        value=f['weight']*(u[0]*f['area'][0]+u[1]*f['area'][2]);mass+=value
        axis=0 if f['area'][0] else 2;side=-1 if f['area'][axis]<0 else 1;name=str((axis,side));by_face[name]=by_face.get(name,0.)+value
    isotropic=-2*v[outer,2]*divergence[outer]/3
    return dict(absolute_mass_outward=mass,scaled_mass_defect=abs(mass)/construction['Z'],mass_faces=by_face,
        mass_scale='U(Hh)*Z, U(Hh)=1',pointwise_divergence_RMS=float(np.sqrt(w@(divergence**2)/w.sum())),
        pointwise_divergence_maximum=float(np.max(abs(divergence))),divergence_scale='U(Hh)/Hh=1',
        outer_isotropic_production_RMS=float(np.sqrt(w[outer]@(isotropic**2)/w[outer].sum())),
        outer_isotropic_production_maximum=float(np.max(abs(isotropic))),
        pressure_level='absolute total traction, no gauge shift applied',audit=o['audit_alpha_zero'])


def compare(a,b,weights,construction):
    oa,ob=a['observation'],b['observation'];require(abs(oa['time']-ob['time'])<1e-13,'endpoint_time_mismatch')
    require(oa['grid']==ob['grid'] and oa['grid_sides']==ob['grid_sides'] and oa['grid_regions']==ob['grid_regions'],'sample_identity_mismatch')
    va,vb=np.array(oa['grid_values']),np.array(ob['grid_values']);require(va.shape==vb.shape==(len(weights),5),'sample_array_mismatch')
    w=np.array(weights);outer=np.array(oa['grid_regions'])==1
    c=construction;scale=[1.,1.,c['K0'],c['reference_column']['counterfactual']['epsilon_amplitude']/c['H'],1.]
    rows={}
    for j,name in enumerate(('u','w','k','outer_epsilon','Pi')):
        m=outer if j==3 else np.ones(len(w),bool);d=va[m,j]-vb[m,j];rms=float(np.sqrt(w[m]@(d*d)/w[m].sum()))
        rows[name]=dict(weighted_RMS=rms,scaled_weighted_RMS=rms/scale[j],scale=scale[j],limit=1e-3,
            status='pass' if rms/scale[j]<=1e-3 else 'fail')
    return dict(status='pass' if all(x['status']=='pass' for x in rows.values()) else 'fail',fields=rows,
        metric='physical weighted RMS, absolute Pi level from traction, no pressure gauge subtraction',time=oa['time'])


def constitutive_comparison(packet,construction):
    defects={k:0. for k in ('R','Q','P','nu_t','epsilon_production','epsilon_destruction')};ratios=[];finiteRe=[]
    samples=packet['native_constitutive']['samples'];require(samples,'missing_model_source_observations')
    for row in samples:
        expected=oracle(construction['model'],row['k'],row['epsilon'],np.array(row['raw_gradient']))
        for name,key in [('R','R'),('Q','Q'),('P','production'),('nu_t','nut'),('epsilon_production','eps_production'),('epsilon_destruction','eps_destruction')]:
            d=float(np.max(abs(np.array(row[name])-np.array(expected[key]))))/max(1.,float(np.max(abs(expected[key]))));defects[name]=max(defects[name],d)
        if row['nu_t']>0:ratios.append(construction['nu']/row['nu_t'])
        finiteRe.append(np.sqrt(construction['nu']*row['epsilon'])/row['k'])
    require(max(defects.values())<=2e-10,'independent_model_stress_or_source_mismatch')
    return dict(status='pass',maximum_scaled_defects=defects,samples=len(samples),
        outer_nu_over_nut_range=[min(ratios),max(ratios)],finite_Re_denominator_indicator_range=[float(min(finiteRe)),float(max(finiteRe))],
        scope='actual raw gradients, full signed isotropic production work; finite-Re indicators assessed without changing high-Re equations')


def budget_classification(packet,construction):
    b=packet['physical_scalar_budget'];c=construction;K=c['K0'];eH=c['reference_column']['counterfactual']['epsilon_amplitude']/c['H']
    equations={};total_k={k:sum(b['equations'][r][k] for r in ('inner_k','outer_k')) for k in b['equations']['inner_k']}
    interfaces={0:0.,1:0.};current_interfaces={0:0.,1:0.}
    for f in b['physical_faces']:
        if f['axis']==2 and f['coordinate'][2]==c['H']:
            interfaces[int(f['region'])]+=f['weight']*f['k_frozen_flux']
            current_interfaces[int(f['region'])]+=f['weight']*f['k_current_flux']
    total_k['stage_physical_boundary_flux']-=sum(interfaces.values())
    total_k['instantaneous_physical_boundary_flux']-=sum(current_interfaces.values())
    total_k['frozen_stage_physical_defect']=total_k['storage_backward_difference']-total_k['stage_source']-total_k['stage_physical_boundary_flux']
    total_k['instantaneous_physical_defect']=total_k['storage_backward_difference']-total_k['instantaneous_source']-total_k['instantaneous_physical_boundary_flux']
    total_k['frozen_operator_splitting']=total_k['stage_source']+total_k['stage_physical_boundary_flux']-total_k['instantaneous_source']-total_k['instantaneous_physical_boundary_flux']
    for name,row in dict(**b['equations'],whole_k=total_k).items():
        # Sources/storage have units scalar*volume/time. Keep the same fixed
        # Hh/U(Hh) scales, rather than a local source that may vanish.
        volume=(c['x'][-1]-c['x'][0])*(c['knots'][-1]-c['knots'][0])
        require(volume>0,'nonpositive_physical_budget_geometry')
        scale=(eH if name=='outer_epsilon' else K)*volume
        equations[name]=dict(**row,physical_scale=scale,scaled_frozen_defect=abs(row['frozen_stage_physical_defect'])/scale,
            accuracy_limit=.02,status='pass' if abs(row['frozen_stage_physical_defect'])/scale<=.02 else 'fail')
    rec=b['reconciliation'];reconstruction={}
    if rec:
        for name,row in rec['assembled'].items():
            defect=max(abs(x) for x in row['reconstruction_defect']);normal=max(1.,max(abs(x) for x in row['native_raw_frozen_load']))
            reconstruction[name]=dict(maximum_defect=defect,normalization=normal,scaled_defect=defect/normal,status='pass' if defect/normal<=2e-10 else 'fail')
    return dict(reconstruction_status='pass' if reconstruction and all(r['status']=='pass' for r in reconstruction.values()) else 'fail',
        accuracy_status='pass' if all(r['status']=='pass' for r in equations.values()) else 'fail',equations=equations,
        reconstruction=reconstruction,interface_frozen_k_fluxes=interfaces,interface_current_k_fluxes=current_interfaces,
        face_order=b['face_order'],actual_quadrature=packet['observation']['actual_quadrature'],physical_volume=volume,
        scope='independent physical gradient fluxes, dependent-epsilon storage retained; prescribed-row reactions not substituted; inner/outer interface flux cancels only when physically matching')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--feasibility',type=Path,required=True);p.add_argument('--contract',type=Path,required=True);p.add_argument('--verification',type=Path,required=True);p.add_argument('--controls',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    for f in (a.verification,a.controls):
        receipt=json.loads(f.read_text());require(receipt['status']=='pass' and all(r['status']=='pass' for r in receipt['checks'].values()),'boundary_capability_not_admitted:'+str(f))
        request=json.loads((f.parent/'request-manifest.json').read_text())['request'];bound=request['runtime']
        solver_identity=bound['IncNavierStokesSolver'] if 'IncNavierStokesSolver' in bound else bound
        require(solver_identity==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign_boundary_verification_runtime')
        sealed=json.loads((f.parent/'completed-artifacts.json').read_text())
        require(all(sha(f.parent/name)==value for name,value in sealed['files'].items()),'changed_boundary_verification_artifact')
    manifest=RequestManifest(a.output,dict(schema='CTAcceptedOpenPilotV1',contract_sha256=sha(a.contract),feasibility_sha256=sha(a.feasibility),
        code_sha256={name:sha(Path(__file__).with_name(name)) for name in ('ct_accepted_open.py','ct_accepted_pilot.py','ct_accepted_feasibility.py')},
        required_verification={str(f):sha(f) for f in (a.verification,a.controls)},
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        horizons={'A95':.2,'A93':.002,'E94E95':.002},dt={'A95':[.002,.001],'A93':[.001],'E94E95':[.001]},
        initialization='CTAcceptedFiniteDepthInitialV1, fresh zero-step native histories for every dt',
        limits=dict(temporal_scaled_RMS=1e-3,mass=1e-8,physical_budget=.02),
        quadrature_controls=[[48,12,12,12,12,16],[64,16,16,16,16,16],[64,16,16,16,16,32]]))
    result=dict(status='fail',cases={},executions=[],experimental_prediction_validation='not_attempted')
    endpoints={};initials={};observed={};constructions={}
    def persist():write_json(a.output/'progress.json',result)
    for model,dts,T in [('A95',[.002,.001],.2),('A93',[.001],.002),('E94E95',[.001],.002)]:
        for dt in dts:
            name=model+'-dt-'+str(dt);row=dict(status='fail');result['cases'][name]=row
            try:
                zero=a.output/(name+'-zero');pilot_case(zero,a.runtime,a.feasibility,model,0,dt)
                c=json.loads((zero/'construction.json').read_text());constructions[name]=c
                e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero,'solver.log');result['executions'].append(dict(case=name,phase='native_zero',**e))
                require((zero/'case.evolution.json').exists(),'native_initialization_produced_no_report: '+(zero/'solver.log').read_text()[-2000:])
                native=json.loads((zero/'case.evolution.json').read_text());require(e['exit_code']==0 and not native['failure'],'native_initialization_failed:'+str(native.get('failure')))
                points,weights=grid(c);initial,e=probe(a.runtime,zero,dict(operation='functions',grid=points,grid_sides=[0]*len(points),export_basis=True),'initial')
                result['executions'].append(e);require(e['exit_code']==0 and initial['observation']['state']==native['state'],'initial_effective_native_input_mismatch')
                initials[name]=initial;row.update(unknowns=len(native['state']),physical_unknowns=len(native['state'])-1,initial_state_sha256=exact_vector_hash(native['state']),
                    initial_history_sha256=artifact_hash(native['checkpoint']),boundary=initial['initialization_receipt'],initial_native_report_sha256=sha(zero/'case.evolution.json'),
                    initial_observations=physical_observations(initial,weights,c),construction=c)
                if model=='A95' and dt==.001:require(row['initial_state_sha256']==result['cases']['A95-dt-0.002']['initial_state_sha256'],'same_space_timestep_initial_coefficients_differ')
                folder=a.output/name;pilot_case(folder,a.runtime,a.feasibility,model,round(T/dt),dt,Path(native['checkpoint']))
                e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(dict(case=name,phase='evolution',**e))
                report=json.loads((folder/'case.evolution.json').read_text());row.update(exit_code=e['exit_code'],time=report['time'],steps=report['steps'],
                    failure=report['failure'],native_failure=report['failure'],report_sha256=sha(folder/'case.evolution.json'),state_sha256=exact_vector_hash(report['state']),
                    extrema=audited_extrema(report),initialization_receipt=report['initialization_receipt'],
                    checkpoint=report['checkpoint'],history_sha256=artifact_hash(report['checkpoint']) if report['checkpoint'] else None,
                    pimple_stage_count=len(report['pimple_attempts']),zero_correction_calls=all(x.get('bound_calls',0)==0 for x in report['pimple_attempts']))
                require(e['exit_code']==0 and not report['failure'] and abs(report['time']-T)<1e-13 and report['steps']==round(T/dt),'trajectory_stopped_or_incomplete')
                request=replay(initial,report['state'],report['time'],points,export_basis=True,export_constitutive=True)
                packet,e=probe(a.runtime,folder,request,'endpoint');result['executions'].append(e);require(e['exit_code']==0,'endpoint_actual_domain_rejected')
                row['endpoint_observations']=physical_observations(packet,weights,c);require(row['endpoint_observations']['scaled_mass_defect']<=1e-8,'physical_mass_defect')
                row['independent_constitutive']=constitutive_comparison(packet,c)
                previous=report['accepted_states'][-2];snapshot=json.loads(Path(report['pimple_attempts'][-1]['k_snapshot']).read_text())
                budgets={}
                for label,orders,face in [('native',[32,8,10,10,10],16),('q48f16',[48,12,12,12,12],16),('q64f16',[64,16,16,16,16],16),('q64f32',[64,16,16,16,16],32)]:
                    controls=folder/label;controls.mkdir();root=ET.parse(folder/'case.xml').getroot()
                    for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,key,value)
                    write(root,controls/'case.xml')
                    req=replay(initial,report['state'],report['time'],[],export_basis=False,export_physical_scalar_budget=True,
                        previous_state=previous,stage_history=snapshot['history'],export_frozen_scalar_reconciliation=True,physical_face_order=face)
                    b,e=probe(a.runtime,controls,req,'budget');result['executions'].append(e);require(e['exit_code']==0,'frozen_budget_replay_failed:'+label)
                    budgets[label]=budget_classification(b,c)
                row['physical_budgets']=budgets;row['status']='pass';endpoints[name]=packet;observed[name]=report
            except Exception as error:row['failure']=repr(error)
            persist()
    if all(n in endpoints for n in ('A95-dt-0.002','A95-dt-0.001')):
        name='A95-dt-0.002';_,weights=grid(constructions[name]);result['timestep_comparison']=compare(endpoints[name],endpoints['A95-dt-0.001'],weights,constructions[name])
    else:result['timestep_comparison']=dict(status='unattempted',reason='one or both intended A95 endpoints unavailable')
    result['execution_status']='pass' if all(r['status']=='pass' for r in result['cases'].values()) else 'fail'
    result['physical_budget_accuracy_status']='pass' if result['execution_status']=='pass' and all(b['accuracy_status']=='pass' and b['reconstruction_status']=='pass' for r in result['cases'].values() for b in r['physical_budgets'].values()) else 'failed_or_unresolved'
    result['status']='pass' if result['execution_status']=='pass' and result['timestep_comparison']['status']=='pass' and result['physical_budget_accuracy_status']=='pass' else 'fail'
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']+[str(f.relative_to(a.output)) for f in a.output.rglob('case.xml')])
    print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
