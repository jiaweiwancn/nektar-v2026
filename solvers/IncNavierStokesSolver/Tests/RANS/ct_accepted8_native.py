#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh CT8 native admission and finite unforced strips, primitive gated.

Use the owned command launcher. New immutable destinations are required.
"""
import argparse,copy,json,math,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted_feasibility import write_json,require,parameters
from ct_accepted8_reference import column,reference,POLICY,PROVIDER,FIT
from ct_accepted4_run import make as old_make,traction
from ct_accepted7_top import classify as operator_classify
from ct_accepted4_analysis import weak
from ct_accepted3_metrics import domain
from ct_accepted_pilot import replay,budget_classification,constitutive_comparison,artifact_hash
from ct_accepted_budget_gate import physical_budgets
from ct_accepted2_audit import values
from ct_accepted_open import run,probe,row_value
from ct_accepted6_audit import captured_k_matrix
from ct_accepted8_gate import pair as pair_gate
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha

def read(p):return json.loads(Path(p).read_text())
def col(c):return column(c['K0'],c['H'],c['roughness_input'],c['nu'],c['displacement'])

def make(folder,rt,c,dt,steps,parent=None,template=None):
    xml=old_make(folder,rt,c,dt,steps,parent,template);r=ET.parse(xml).getroot()
    for k,v in dict(RANSWallBoundaryPolicy='AcceptedRichardsNorrisFlatV1',RANSOpenTop='RichardsNorris2011V1',
        RANSOpenReference=PROVIDER,RANSWallDisplacementPolicy=POLICY,RANSOpenFitPolicy=FIT).items():set_info(r,k,v)
    for k,v in dict(RANSOpenDisplacement=c['displacement'],RANSWallMassOrder=12,RANSWallVolumeOrder=12,RANSWallInterfaceOrder=64).items():set_parameter(r,k,v)
    write(r,xml);write_json(folder/'construction.json',dict(c,dt=dt,steps=steps,configuration_sha256=sha(xml)))
    return xml

def inventory(c):
    x,w=leggauss(6);grid=[];weights=[]
    cuts=sorted(set(c['knots']+[c['H'],c['roughness_input'],c['reference_column']['cutoff']]))
    for a,b in zip(c['x'][:-1],c['x'][1:]):
        for z,t in zip(cuts[:-1],cuts[1:]):
            for i in range(6):
                for j in range(6):grid.append([(a+b)/2+(b-a)*x[i]/2,0.,(z+t)/2+(t-z)*x[j]/2]);weights.append((b-a)*(t-z)*w[i]*w[j]/4)
    n=len(grid);sides=[0]*n
    xs=sorted({a+(b-a)*j/p for a,b in zip(c['x'][:-1],c['x'][1:]) for p in (3,4) for j in range(p+1)})
    for xx in xs:
        for zz in sorted(set([0.,.175,c['displacement'],c['H']/2,c['reference_column']['cutoff'],c['roughness_input'],c['Z']])):
            grid.append([xx,0.,zz]);sides.append(0)
        for side in (-1,1):grid.append([xx,0.,c['H']]);sides.append(side)
    return grid,sides,weights,n

def request(c,operator=False,tensor=True):
    g,s,_,_=inventory(c)
    return dict(operation='functions',grid=g,grid_sides=s,export_basis=True,export_test_basis=True,export_constitutive=True,
        export_neutral_composite_reference=True,export_resolution=operator,export_open_verification=operator,
        export_mixed_space=operator,export_tensor_validation=operator and tensor)

def binding(p,c):
    o=p['observation'];require(o['resolved_descriptor'].startswith('A95:OriginalHighRe-v1'),'foreign actual descriptor')
    require(POLICY in o['identity'] and PROVIDER in o['identity'] and FIT in o['identity'],'missing physical/provider/fit signature')
    require(p['neutral_composite_reference']['provider']==PROVIDER,'foreign reference provider')
    r=p['initialization_receipt'];require(r['open_explicit_implicit_agreement'] and r['explicit_boundary_sha256']==r['implicit_boundary_sha256'],'live boundary mismatch')
    return r

def oracle_agreement(p,c):
    at=lru_cache(None)(lambda z,r:reference(col(c),z,r));defect=np.zeros(9)
    A=col(c).A;H=c['H'];K=c['K0'];scales=np.array([c['ustar'],K,A/(H-c['displacement']),c['ustar']**2,c['ustar']/H,K/H,A/H**2,A,col(c).D(H)])
    nonconstant=0.
    for row in p['neutral_composite_reference']['samples']:
        z=row['coordinate'][2];r=bool(row['region']);a=np.asarray(row['values_U_k_epsilon_Pi_gradients_Q0_D']);b=at(z,r)
        delta=abs(a-b)
        if not r and abs(z-col(c).cutoff)<4e-13*H:delta[6]=0.
        defect=np.maximum(defect,delta/scales)
        if not r:nonconstant=max(nonconstant,abs(b[1]/K-1))
    return dict(status='pass' if max(defect)<=2e-9 and nonconstant>1e-3 else 'fail',maximum_scaled_defects=defect.tolist(),nonconstant_inner_k_signal=nonconstant)

def observations(p,state,n):return np.asarray(values(p,state,n)['observation']['grid_values'])

def norms(a,b,weights,regions,c):
    require(a.shape==b.shape and len(a)==len(weights)==len(regions),'mismatched physical observations')
    scales=[c['ustar'],c['ustar'],c['K0'],col(c).A/(c['H']-c['displacement']),c['ustar']**2];out={}
    for j,name in enumerate(('U','w','k','outer_epsilon','Pi')):
        keep=np.asarray(regions,dtype=bool) if j==3 else np.ones(len(a),dtype=bool);w=np.asarray(weights)[keep];d=(a-b)[keep,j]
        rms=float(np.sqrt(w@(d*d)/w.sum()));out[name]=dict(weighted_RMS=rms,true_L2=float(np.sqrt(w@(d*d))),scale=scales[j],scaled=rms/scales[j])
    return out

def physical_metrics(p,state,c,weights,n):
    g=p['observation']['grid'];r=p['observation']['grid_regions'];s=np.asarray(state);off=p['observation']['offsets']
    v=observations(p,state,n);base=observations(p,p['observation']['state'],n)
    truth=[]
    for q,region in zip(g[:n],r[:n],strict=True):
        a=reference(col(c),q[2],bool(region));truth.append([a[0],0.,a[1],a[2] if region else 0.,a[3]])
    drift=norms(v,base,weights,r[:n],c);error=norms(v,np.asarray(truth),weights,r[:n],c)
    net=0.;ground=[]
    for f in p['observation']['faces']:
        if f['interface']:continue
        vel=np.array([row_value(f['velocity_k'],s,off[b])[0] for b in (0,1)])
        net+=f['weight']*float(vel@np.asarray(f['area'])[[0,2]])
    tr=traction(p,state,c)
    reasons=[name for name,v in drift.items() if v['scaled']>(1e-4 if name=='w' else 1e-3)]
    if abs(net)/(c['ustar']*c['Z'])>1e-8:reasons.append('mass')
    if tr['status']!='pass':reasons.append('traction')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,drift=drift,absolute_reference_error=error,traction=tr,
        mass=dict(outward_net=net,physical_scale=c['ustar']*c['Z'],scaled=abs(net)/(c['ustar']*c['Z'])),
        scope='weighted physical RMS; true L2 separate; absolute pressure, fixed reference and actual initial baseline')

CONTROLS=[('native',(32,8,12,12,64),16),('q48f16',(48,12,16,16,16),16),('q64f16',(64,16,20,20,20),16),('q64f32',(64,16,20,20,20),32)]

def budgets(rt,folder,ini,native,c):
    body=dict(budgets={});attempt=next(a for a in reversed(native['pimple_attempts']) if a['scheduled_complete'])
    snap_path=Path(attempt['k_snapshot']);require(sha(snap_path)==attempt['k_snapshot_sha256'],'changed actual stage');snap=read(snap_path)
    for label,orders,face in CONTROLS:
        dest=folder/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
        for key,val in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,key,val)
        write(root,dest/'case.xml');req=replay(ini,native['state'],native['time'],[],export_basis=True,export_physical_scalar_budget=True,
            export_frozen_scalar_reconciliation=True,previous_state=native['accepted_states'][-2],stage_history=snap['history'],actual_flow_stage=snap['post_flow'],actual_pre_k_stage=snap['pre_k'],physical_face_order=face)
        packet,e=probe(rt,dest,req,'budget');require(e['exit_code']==0,'physical budget rejected:'+str(packet.get('failure')))
        require(packet['observation']['state']==native['state'] and packet['observation']['identity']==ini['observation']['identity'],'budget identity mismatch')
        b=budget_classification(packet,c);b.update(execution=e,packet_sha256=sha(dest/'budget.json'),state_sha256=exact_vector_hash(native['state']))
        body['budgets'][label]=b
        if label=='native':body['captured_k_matrix']=captured_k_matrix(packet,snap)
    body['gate']=physical_budgets(body);return body

def preflight(a):
    c=read(a.reference)['case'];m=RequestManifest(a.output,dict(schema='CT8NativePreflightV1',case=c,contract=sha(a.contract),reference=sha(a.reference),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__))))
    proof_path=a.basis_proof;proof=read(proof_path);proof_result=read(proof_path.parent.parent/'results.json')
    require(proof_result['status']=='pass' and proof_result['checks']['A95-p4']['packet_sha256']==sha(proof_path),'unqualified or changed all-mode basis proof')
    require(read(proof_path.parent.parent/'request-manifest.json')['request']['runtime']==runtime_identity(a.runtime/'RANSPreservationProbe'),'foreign basis verification runtime')
    result=dict(schema='CT8NativePreflightV1',case=c,counts={},small_operator={},initial={},status='fail',basis_proof=dict(path=str(proof_path),sha256=sha(proof_path),
        scope='same rebuilt runtime, native Modified_A P4/P2 transform/derivative/face and mass checks; separate small d0 basis fixture, no replacement of displaced governing rows'))
    # Count both declared degree layouts without allocating a governing matrix.
    for degree in (3,4):
        q=copy.deepcopy(c);q.update(degree=degree,pressure_degree=degree-2);folder=a.output/('count-p'+str(degree));make(folder,a.runtime,q,.002,0)
        packet,e=probe(a.runtime,folder,dict(operation='native_layout_inventory'),'count');require(e['exit_code']==0,'native count failed')
        result['counts'][str(degree)]=dict(packet,packet_sha256=sha(folder/'count.json'))
    require(result['counts']['4']['stored_slots']<=3000,'CT8 declared P4 count exceeds cap; no trajectory admitted')
    # Small physical-data fixture for every free native mode and full matrices.
    # Resolve the actual analytical breakpoints. Tensor transform checks use
    # their separately bound small fixture; the displaced full operator and
    # pressure rank use THESE actual rows/matrices, never substituted bases.
    q=copy.deepcopy(c);q.update(Z=.5,knots=sorted(set([0.,c['H']/4,c['H']/2,col(c).cutoff,c['roughness_input'],c['H'],.3,.5])))
    for degree in (3,4):
        cc=copy.deepcopy(q);cc.update(degree=degree,pressure_degree=degree-2);folder=a.output/('operator-p'+str(degree));make(folder,a.runtime,cc,.002,0)
        packet,e=probe(a.runtime,folder,request(cc,True,tensor=False),'operator');require(e['exit_code']==0,'native displaced operator rejected:'+str(packet.get('failure')))
        binding(packet,cc);out=operator_classify(packet,cc,column=col(cc),reference_function=reference,tensor_proof=proof['tensor_validation'])
        out['continuous_provider']=oracle_agreement(packet,cc)
        require(out['continuous_provider']['status']=='pass','continuous provider mismatch')
        for name,kwargs in [('pressure',dict(wrong_pressure=True)),('constant_inner_k',dict(constant_inner_k=True))]:
            wrong=weak(packet,cc,column=col(cc),reference_function=reference,**kwargs);require(wrong['status']=='fail','nondiscriminating '+name+' negative');out[name+'_negative']=wrong
        out.update(execution=e,packet_sha256=sha(folder/'operator.json'));result['small_operator'][str(degree)]=out
        write_json(a.output/'progress.json',result)
    grid,sides,weights,n=inventory(c);result.update(weights=weights,global_count=n)
    hashes=[]
    for dt in (.002,.001):
        folder=a.output/('initial-'+str(dt));make(folder,a.runtime,c,dt,0)
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');native=read(folder/'case.evolution.json');require(e['exit_code']==0 and not native['failure'],'zero native rejection:'+str(native['failure']))
        packet,pe=probe(a.runtime,folder,request(c),'initial');require(pe['exit_code']==0,'initial physical probe rejected')
        require(packet['observation']['state']==native['state'],'driver/probe initialization differ');bind=binding(packet,c)
        mat=physical_metrics(packet,native['state'],c,weights,n);dom=domain(packet,native,c);ref=oracle_agreement(packet,c)
        out=dict(directory=str(folder),state_sha256=exact_vector_hash(native['state']),history_sha256=artifact_hash(native['checkpoint']),initialization_receipt=bind,
            metrics=mat,domain=dom,reference_agreement=ref,native_report_sha256=sha(folder/'case.evolution.json'),packet_sha256=sha(folder/'initial.json'),raw_unforced_spatial_norms=native['raw_reference_norms'],slots=len(native['state']))
        out['status']='pass' if mat['status']==dom['status']==ref['status']=='pass' else 'fail';result['initial'][str(dt)]=out;hashes.append(out['state_sha256'])
        write_json(a.output/'progress.json',result)
    result['same_space_actual_initial_equal']=len(set(hashes))==1
    # One declared read-only degree contrast on the identical physical grid.
    # Different native coefficients are expected; compare physical values.
    cc=copy.deepcopy(c);cc.update(degree=3,pressure_degree=1);dest=a.output/'physical-p3';make(dest,a.runtime,cc,.002,0)
    pp,pe=probe(a.runtime,dest,request(cc),'projection');require(pe['exit_code']==0,'P3 represented reference domain failed')
    pm=physical_metrics(pp,pp['observation']['state'],cc,weights,n);fine=read(a.output/'initial-0.002/initial.json')
    delta=norms(observations(pp,pp['observation']['state'],n),observations(fine,fine['observation']['state'],n),weights,pp['observation']['grid_regions'][:n],c)
    result['projection_contrast']=dict(coarse_P3=pm,fine_P4=result['initial']['0.002']['metrics'],physical_difference=delta,
        coarse_unanchored_spatial_audit=pp['observation'].get('audit_alpha_zero'),fine_unanchored_spatial_audit=fine['observation'].get('audit_alpha_zero'),
        continuous_inputs_identical=True,grid_and_sides_identical=pp['observation']['grid']==fine['observation']['grid'] and pp['observation']['grid_sides']==fine['observation']['grid_sides'],
        scope='static reference representations, not evolved spatial convergence or a measured degree order')
    result['status']='pass' if result['same_space_actual_initial_equal'] and all(r['status']=='pass' for r in result['initial'].values()) else 'failed_initial_accuracy_or_domain'
    result['flow_through_estimate']=dict(mean_speed=col(c).t0*col(c).J0,time=(c['x'][-1]-c['x'][0])/(col(c).t0*col(c).J0),definition='interface velocity scale; full-height logarithm has no single uniform transit')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result

def advance(a):
    pre=read(a.preflight/'results.json');require(pre['status']=='pass','preflight not admitted');c=pre['case'];parent=read(a.parent/'results.json') if a.parent else None
    if parent:require(parent['qualification']['status']=='pass' and parent['target']==.008,'parent failed/incompatible')
    start=.008 if parent else 0.;target=a.target;require(target in (.008,.02) and target>start,'wrong milestone')
    m=RequestManifest(a.output,dict(schema='CT8BoundedPairV1',case=c,start=start,target=target,dt=[.002,.001],preflight=sha(a.preflight/'results.json'),parent=sha(a.parent/'results.json') if parent else None,contract=sha(a.contract),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source=sha(Path(__file__))))
    result=dict(schema='CT8BoundedPairV1',case=c,target=target,members={},status='incomplete');packets=[]
    for dt in (.002,.001):
        name=str(dt);zero=Path(pre['initial'][name]['directory']);ini=read(zero/'initial.json');source=Path(parent['members'][name]['directory']) if parent else zero;old=read(source/'case.evolution.json')
        require(old['time']==start and old['dt']==dt,'changed history time/method');folder=a.output/name;make(folder,a.runtime,c,dt,round((target-start)/dt),Path(old['checkpoint']),source/'case.xml')
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
        native=read(folder/'case.evolution.json');require(native['accepted_states'][0]==old['state'],'actual consumed history differs')
        row=dict(directory=str(folder),execution=e,time=native['time'],failure=native['failure'],state_sha256=exact_vector_hash(native['state']),native_report_sha256=sha(folder/'case.evolution.json'),history_sha256=artifact_hash(native['checkpoint']),signature=native['signature'],initialization_receipt=native['initialization_receipt'],history_times=native['native_history_times'])
        require(':resolved=A95:OriginalHighRe-v1' in native['signature'] and POLICY in native['signature'],'foreign history descriptor')
        row['accepted_metrics']=[dict(time=t,**physical_metrics(ini,s,c,pre['weights'],pre['global_count'])) for t,s in zip(native['audit_times'],native['accepted_states'],strict=True)]
        row['first_accuracy_crossing']=next((v for v in row['accepted_metrics'] if v['status']!='pass'),None)
        row['domain']=domain(ini,native,c);row['budgets']=budgets(a.runtime,folder,ini,native,c) if native['steps'] else dict(gate=dict(status='missing'))
        g,sides,_,_=inventory(c);req=replay(ini,native['state'],native['time'],g,export_constitutive=True);req['grid_sides']=sides
        packet,pe=probe(a.runtime,folder,req,'endpoint');row['constitutive']=constitutive_comparison(packet,c) if pe['exit_code']==0 else dict(status='fail',reason=packet.get('failure'));row['corrections']=sum(v.get('bound_calls',0) for v in native['pimple_attempts']);row['extra_solves']=sum(v.get('additional_governing_equation_solves',0) for v in native['pimple_attempts'])
        row['refreshed_epsilon_trace_lag']=max((v.get('refreshed_trace_defect',0) for v in native['pimple_attempts']),default=0.)
        result['members'][name]=row;packets.append((ini,native));write_json(a.output/'progress.json',result)
    reasons=[];temporal={};times=[n['time'] for _,n in packets]
    if any(abs(t-target)>1e-13 for t in times):reasons.append('incomplete_endpoints')
    else:
        obs=[observations(p,n['state'],pre['global_count']) for p,n in packets];regions=packets[0][0]['observation']['grid_regions'][:pre['global_count']]
        delta=norms(*obs,pre['weights'],regions,c)
        for name,diff in delta.items():
            drift=max(row['accepted_metrics'][-1]['drift'][name]['scaled'] for row in result['members'].values());limit=max(1e-7,.1*drift)
            temporal[name]=dict(difference=diff,measured_drift=drift,limit=limit,status='pass' if diff['scaled']<=limit else 'unresolved')
        if any(v['status']!='pass' for v in temporal.values()):reasons.append('timestep_uncertainty')
    for name,row in result['members'].items():
        if row['failure'] or row['accepted_metrics'][-1]['status']!='pass' or row['domain']['status']!='pass':reasons.append(name+':physical_or_accuracy')
        if row['budgets']['gate']['status']!='pass' or row['constitutive']['status']!='pass':reasons.append(name+':budget_or_constitutive')
        if row['corrections'] or row['extra_solves']:reasons.append(name+':forbidden_correction')
    result.update(temporal=temporal,status='complete')
    result['qualification']=pair_gate(result)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=('preflight','advance'),required=True)
    for name in ('runtime','reference','contract','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--preflight',type=Path);p.add_argument('--parent',type=Path);p.add_argument('--target',type=float)
    p.add_argument('--basis-proof',type=Path)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    try:r=preflight(a) if a.mode=='preflight' else advance(a)
    except Exception as error:
        a.output.mkdir(exist_ok=True,parents=True);write_json(a.output/'failure.json',dict(schema='CT8ExecutionFailureV1',status='fail',failure=repr(error)));raise
    status=r.get('qualification',r)['status'];print('CT8',a.mode,status,flush=True);return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
