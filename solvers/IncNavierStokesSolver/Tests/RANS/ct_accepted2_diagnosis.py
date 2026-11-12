#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT2 bounded mesh matrix and independent computed-flow diagnosis.
Invoke ONLY through ct_accepted_command.py. Heavy packets remain in --output.
"""
import argparse,json,os,copy,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted_feasibility import write_json,parameters
from ct_accepted_open import run,probe,require,row_value,fixture
from ct_accepted_pilot import replay,compare,physical_observations,budget_classification
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_picard import set_parameter,set_info,write
from rans_request_manifest import RequestManifest,runtime_identity,sha
from compare_cls_epawf_pimple import same_history
from validate_a93_a95_e94e95 import oracle

WITNESS=[.24836915803232323,0.,6.2989252590303035]
LAYOUTS={'B':([],[],1085),'X':([.1875],[],1545),'Zplus':([],[6.5],1172),'XZplus':([.1875],[6.5],1669),'Zminus':([],[5.125],1172)}

def function_identity(root):
    return {n.get('NAME'):{e.get('VAR'):e.get('VALUE') for e in n.findall('E')} for n in root.findall('./CONDITIONS/FUNCTION')}
def analytic(root,suffix,variable,y):
    f=root.find("./CONDITIONS/FUNCTION[@NAME='InitialConditions"+suffix+"']/E[@VAR='"+variable+"']").get('VALUE')
    return np.asarray(eval(f.replace('^','**'),{'__builtins__':{}},{'x':0.,'y':np.asarray(y),'z':0.,'t':0.,'sqrt':np.sqrt,'abs':np.abs,'exp':np.exp,'log':np.log}),float)+np.zeros_like(y)
def make_case(folder,source,c,layout,dt,steps,restart=None):
    folder.mkdir(parents=True,exist_ok=False);root=ET.parse(source).getroot();before=function_identity(root)
    xs,zs,_=LAYOUTS[layout]
    if xs or zs:tensor_mesh(root,sorted(set(c['knots']+zs)),3,2,.25,1.,points=[8,8],xknots=sorted(set(c['x']+xs)))
    set_parameter(root,'TimeStep',dt);set_parameter(root,'NumSteps',steps)
    set_info(root,'RANSWallRestart',str(restart) if restart else 'None')
    set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'))
    set_info(root,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'))
    require(function_identity(root)==before,'physical_FUNCTION_expressions_changed')
    write(root,folder/'case.xml')
    nc=dict(c,x=sorted(set(c['x']+xs)),knots=sorted(set(c['knots']+zs)),dt=dt,steps=steps)
    write_json(folder/'construction.json',nc)
    return nc

def grids(c):
    # One physical grid for ALL layouts; union partition and physical breaks.
    xs=[0.,.125,.1875,.25];zs=sorted(set(c['knots']+[5.125,6.5]));x,w=leggauss(4)
    points=[];weights=[]
    for a,b in zip(xs[:-1],xs[1:]):
        for d,e in zip(zs[:-1],zs[1:]):
            for i in range(4):
                for j in range(4):points.append([(a+b)/2+(b-a)*x[i]/2,0.,(d+e)/2+(e-d)*x[j]/2]);weights.append((b-a)*(e-d)*w[i]*w[j]/4)
    n_global=len(points)
    zz=sorted(set(np.linspace(5.8,6.8,81).tolist()+np.linspace(6.251,6.35,41).tolist()+[WITNESS[2],6.25-1e-8,6.25+1e-8,6.5-1e-8,6.5+1e-8]))
    for xx in (.125+1e-8,.17,.1875-1e-8,.1875+1e-8,.22,.24,WITNESS[0],.2499,.25):
        for z in zz:points.append([xx,0.,z])
    points.append(WITNESS)
    return points,weights,n_global

def invariant(G,k,e):
    Sd=(G+G.T)/2-np.eye(3)*np.trace(G)/3;O=(G-G.T)/2
    ev=np.linalg.eigvalsh(Sd);sn=float(np.linalg.norm(Sd));on=float(np.linalg.norm(O));tau=k/e
    if sn==0 and on!=0:return dict(status='undefined_pure_rotation',raw_gradient=G.tolist(),k=k,epsilon=e,raw_divergence=float(np.trace(G)))
    C=1/6.5 if sn==0 else 1/(6.5+3*ev[-1]/sn*np.hypot(sn,on)*tau)
    rad=float(1-9*C*C*sn*sn*tau*tau)
    return dict(status='admissible' if k>0 and e>0 and rad>=0 else 'constitutive_domain_failure',raw_gradient=G.tolist(),raw_divergence=float(np.trace(G)),strain_eigenvalues=ev.tolist(),strain_norm=sn,rotation_norm=on,k=float(k),epsilon=float(e),tau=float(tau),Cmu=float(C),radicand=rad,
        isotropic_production_work=float(-2*k*np.trace(G)/3),full_stress_and_production='defined' if rad>=0 else 'undefined: original negative A95 radicand')
def at(row,s,off):
    G=np.zeros((3,3));G[0]=row_value(row['velocity_k'],s,off[0])[1];G[2]=row_value(row['velocity_k'],s,off[1])[1]
    k=row_value(row['velocity_k'],s,off[2])[0];e=row_value(row['epsilon'],s,off[3])[0]
    d=invariant(G,k,e);pr,dp=row_value(row['pressure'],s,off[5]);d.update(Pi=float(pr),pressure_gradient=dp.tolist())
    planar=G.copy();planar[0,0]-=np.trace(G)/2;planar[2,2]-=np.trace(G)/2
    d['diagnostic_planar_trace_removed']=invariant(planar,k,e)
    if d['status']=='admissible':
        law=oracle('A95',k,e,G);T=-law['R']+np.eye(3)*(2*k/3-pr)
        d.update(R=law['R'].tolist(),Q=law['Q'].tolist(),production=float(law['production']),total_stress=T.tolist())
    return d

def inventory(rows,coords,s,off,region=None):
    minimum=None;negative=0;div2=0.;divmax=0.;count=0
    for i,row in enumerate(rows):
        if row['epsilon'] is None:continue
        v=at(row,s,off);count+=1;div2+=v['raw_divergence']**2;divmax=max(divmax,abs(v['raw_divergence']))
        if v['status']!='admissible':negative+=1
        if minimum is None or v.get('radicand',-1e300)<minimum.get('radicand',-1e300):
            minimum=dict(v,coordinate=coords[i],geometry=row['velocity_k']['geometry'],element=row['velocity_k'].get('element'),side=0,sample_index=i)
    return dict(minimum=minimum,domain_failure_count=negative,count=count,divergence_unweighted_RMS=float(np.sqrt(div2/count)),divergence_maximum=divmax)

def independent_pressure_action(o,direction,hr):
    off=np.array(o['offsets'],int);s=np.array(direction);result=np.zeros(len(s));rough=np.zeros(len(s))
    elements={x['geometry']:x for x in o['elements']}
    outerH=[f for f in o['faces'] if f['interface'] and f['region']==1]
    for q in o['volume']:
        row=q['velocity_k'];ix=np.array(row['ids'],int);p=row_value(q['pressure'],s,off[4+q['region']])[0]
        for b,axis in enumerate((0,2)):result[ix+off[b]]-=q['weight']*p*np.array(row['derivative'][axis])
        if q['region']==0 and q['coordinate'][2]<=hr:
            x=q['coordinate'][0]
            face=next(f for f in outerH if elements[f['velocity_k']['geometry']]['lower'][0]<=x<=elements[f['velocity_k']['geometry']]['upper'][0])
            dp=row_value(face['pressure'],s,off[5])[1][0]
            rough[ix+off[0]]-=q['weight']*dp*np.array(row['phi'])
    return result+rough,rough

def rough_control(runtime,output):
    output.mkdir();xml=fixture(output/'case',runtime);root=ET.parse(xml).getroot();set_parameter(root,'RANSWallRoughness',.1);write(root,xml)
    packet,e=probe(runtime,xml.parent,dict(operation='functions',grid=[],export_basis=True,export_open_verification=True),'rough')
    require(e['exit_code']==0 and packet['status']=='pass','rough_admissible_fixture_failed')
    o=packet['observation'];f=packet['open_frozen_systems'];off=np.array(o['offsets'],int);s=np.zeros(len(o['state']))
    for b in (4,5):
        M=np.array(f['masses'][b]);load=np.zeros(len(M))
        for q in o['pressure_mass']:
            if q['region']==b-4:
                row=q['pressure'];load[np.array(row['ids'],int)]+=q['weight']*np.array(row['phi'])*q['coordinate'][0]
        s[off[b]:off[b+1]]=np.linalg.solve(M,load)
    expected,rough=independent_pressure_action(o,s,.1)
    flow=next(x for x in f['blocks'] if x['block']==0);ids=np.array(flow['rows'],int);nv=sum(ids<off[4])
    actual=np.array(flow['matrix'])[:nv]@s[ids]
    defect=float(np.max(abs(actual-expected[ids[:nv]])));omitted=float(np.max(abs(rough[ids[:nv]])))
    require(defect<=2e-10 and omitted>1e-6,'roughness_pressure_volume_action_unverified')
    return dict(status='pass',configuration_sha256=sha(xml),report_sha256=sha(xml.parent/'rough.json'),hr=.1,actual_unknowns=len(s),pressure_direction='independent mass projection p=x in both regions',maximum_action_defect=defect,omitted_roughness_defect=omitted,expected_roughness_volume_load=rough.tolist(),execution=e)

def pressure_projection(xml,c,o,basis):
    root=ET.parse(xml).getroot();z=np.array([x for x in c['knots'] if x>=c['H']]);n=len(z);M=np.zeros((n,n));load=np.zeros(n);g,w=leggauss(8)
    for i,(a,b) in enumerate(zip(z[:-1],z[1:])):
        points=(a+b)/2+(b-a)*g/2;phi=np.array([(1-g)/2,(1+g)/2]);p=analytic(root,'Outer','p',points)
        M[i:i+2,i:i+2]+=(b-a)/2*(phi*w)@phi.T;load[i:i+2]+=(b-a)/2*(phi@(w*p))
    coeff=np.linalg.solve(M,load);j=int(np.searchsorted(z,WITNESS[2])-1);gradient=(coeff[j+1]-coeff[j])/(z[j+1]-z[j])
    s=np.array(o['state']);off=np.array(o['offsets'],int);native=row_value(basis[-1]['pressure'],s,off[5])[1][2]
    require(abs(gradient-native)<2e-11,'independent_initial_pressure_projection_mismatch')
    pressure_load=np.zeros(off[4]);stress_load=np.zeros(off[4]);common=np.zeros(3)
    for q in o['volume']:
        if not q['region']:continue
        row=q['velocity_k'];ix=np.array(row['ids'],int)+off[1];p=row_value(q['pressure'],s,off[5])[0]
        law=at(q,s,off);Q=np.array(law['Q']);G=np.array(law['raw_gradient'])
        pressure_load[ix]-=q['weight']*np.array(row['derivative'][2])*p
        stress_load[ix]-=q['weight']*(np.array(row['derivative'][2])*Q[2,2]+np.array(row['derivative'][0])*Q[2,0])
        x,_,h=q['coordinate'];v=np.exp(-((h-6.25)/1.5)**2)*x*(.25-x)
        dx=np.exp(-((h-6.25)/1.5)**2)*(.25-2*x);dz=-2*(h-6.25)/1.5**2*v
        common+=q['weight']*np.array([-p*dz,-Q[2,2]*dz-Q[2,0]*dx,0.])
    return dict(status='pass',independent_gradient_at_witness=float(gradient),native_gradient_at_witness=float(native),gradient_defect=float(abs(native-gradient)),nodes=z.tolist(),independent_nodal_projection=coeff.tolist(),unmasked_w_pressure_load_l2=float(np.linalg.norm(pressure_load)),unmasked_w_Q_load_l2=float(np.linalg.norm(stress_load)),signed_common_test=dict(pressure=float(common[0]),nonlinear_stress=float(common[1]),sum=float(common[:2].sum())),common_test='phi=x*(.25-x)*exp(-((z-6.25)/1.5)^2), full volume without residual subtraction',analytical_Pi_above_delta=0.,projection_order=8)

def initialize(runtime,previous,output,layouts):
    output.mkdir();construction=json.loads((previous/'pilot-v3/A95-dt-0.002/construction.json').read_text());source=previous/'pilot-v3/A95-dt-0.002-zero/case.xml'
    points,weights,ng=grids(construction);result=dict(schema='CT2InitialMeshControlsV1',status='pass',layouts={},runtime=runtime_identity(runtime/'IncNavierStokesSolver'),fixed_grid=points,global_weights=weights,global_count=ng)
    for layout in layouts:
        row={};result['layouts'][layout]=row
        for dt in (.002,.001):
            name=layout+'-'+str(dt);folder=output/name;c=make_case(folder,source,construction,layout,dt,0)
            e=run([runtime/'IncNavierStokesSolver','case.xml'],folder,'zero.log');require(e['exit_code']==0,'zero_step_failed:'+name)
            native=json.loads((folder/'case.evolution.json').read_text())
            packet,pe=probe(runtime,folder,dict(operation='functions',grid=points,grid_sides=[0]*len(points),export_basis=True,export_test_basis=True),'initial')
            require(pe['exit_code']==0 and packet['status']=='pass' and packet['observation']['state']==native['state'],'initial_actual_state_mismatch:'+name)
            require(len(native['state'])==LAYOUTS[layout][2],'predicted_count_not_confirmed:'+layout)
            o=packet['observation'];basis=packet['test_basis'];off=np.array(o['offsets'],int)
            entry=dict(status='pass',slots=len(native['state']),state_sha256=exact_vector_hash(native['state']),report_sha256=sha(folder/'case.evolution.json'),packet_sha256=sha(folder/'initial.json'),configuration_sha256=sha(folder/'case.xml'),initialization_receipt=native['initialization_receipt'],functions=function_identity(ET.parse(folder/'case.xml').getroot()),executions=[e,pe],native_domain=inventory(o['volume'],[q['coordinate'] for q in o['volume']],np.array(native['state']),off),fixed_grid_domain=inventory(basis,points,np.array(native['state']),off),pressure_projection=pressure_projection(folder/'case.xml',c,o,basis))
            row[str(dt)]=entry
            if dt==.001:require(entry['state_sha256']==row['0.002']['state_sha256'],'same_space_initial_mismatch')
            write_json(output/'results.json',result)
    return result

def analyze_attempt(initial,report,folder):
    o=initial['observation'];off=np.array(o['offsets'],int);attempts=[]
    for attempt in report['pimple_attempts']:
        require('flow_snapshot' in attempt,'missing_computed_flow_snapshot')
        path=Path(attempt['flow_snapshot']);f=json.loads(path.read_text());row=dict(time=f['time'],snapshot_sha256=sha(path),stages={})
        for name in ('input','momentum_predictor','corrected_flow'):
            s=np.array(f[name]);row['stages'][name]=dict(native=inventory(o['volume'],[q['coordinate'] for q in o['volume']],s,off),fixed=inventory(initial['test_basis'],o['grid'],s,off))
        row['earliest_sampled_domain_loss']=next((x for x in ('input','momentum_predictor','corrected_flow') if any(row['stages'][x][kind]['domain_failure_count'] for kind in ('native','fixed'))),'none')
        witness=row['stages']['corrected_flow']['native']['minimum'];coord=witness['coordinate']
        # Use the SAME original native row for every contribution; do not mix minimizers.
        ix=witness['sample_index'];q=o['volume'][ix];parts={}
        for name in ('input','momentum_predictor','corrected_flow'):parts[name]=at(q,np.array(f[name]),off)
        g0=np.array(parts['input']['raw_gradient']);gp=np.array(parts['momentum_predictor']['raw_gradient']);gc=np.array(parts['corrected_flow']['raw_gradient'])
        row['same_native_point']=dict(coordinate=coord,geometry=q['generating_geometry'],stages=parts,predictor_gradient_increment=(gp-g0).tolist(),pressure_mass_gradient_increment=(gc-gp).tolist())
        expected,rough=independent_pressure_action(o,np.array(f['corrected_flow'])-np.array(f['momentum_predictor']),folder['roughness_input'])
        ids=np.array(f['rows'],int);nv=f['velocity_size'];actual=np.array(f['full_pressure_action'])
        action_defect=float(np.max(abs(expected[ids[:nv]]-actual[ids[:nv]])))
        row['independent_full_pressure_action_defect']=action_defect;row['roughness_pressure_action_maximum']=float(np.max(abs(rough)))
        require(action_defect<2e-9,'independent_flow_pressure_action_mismatch')
        for key in ('predictor_linear_residual','pressure_linear_residual','splitting_identity_defect'):
            row[key]=float(np.max(abs(np.array(f[key]))))
        row['continuity']=float(np.max(abs(np.array(f['frozen_momentum_and_continuity_defect'])[nv:])))
        attempts.append(row)
    return attempts

def advance(runtime,previous,preparation,output,layouts,T=.008):
    output.mkdir();result=dict(schema='CT2MatchedMeshTrajectoriesV1',diagnostic_execution='pass',physical_completion='fail',cases={},comparisons={})
    endpoints={};constructions={}
    for layout in layouts:
        for dt in (.002,.001):
            name=layout+'-'+str(dt);zero=preparation/name;ini=json.loads((zero/'initial.json').read_text());c=json.loads((zero/'construction.json').read_text());native0=json.loads((zero/'case.evolution.json').read_text())
            folder=output/name;make_case(folder,zero/'case.xml',c,'B',dt,round(T/dt),Path(native0['checkpoint']))
            e=run([runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');report=json.loads((folder/'case.evolution.json').read_text())
            row=dict(execution=e,requested_time=T,time=report['time'],steps=report['steps'],failure=report['failure'],native_report_sha256=sha(folder/'case.evolution.json'),state_sha256=exact_vector_hash(report['state']),history=report['native_history_times'],checkpoint=report['checkpoint'],status='pass' if not report['failure'] and abs(report['time']-T)<1e-13 else 'physical_failure',attempts=analyze_attempt(ini,report,c))
            result['cases'][name]=row;constructions[name]=c
            if layout=='B':
                old=json.loads((previous/('pilot-v3/A95-dt-'+str(dt))/'case.evolution.json').read_text())
                require(same_history(old,report) and old['rejected_state']==report['rejected_state'] and old['failure']==report['failure'],'instrumentation_changed_sealed_failure')
                row['sealed_baseline_bitwise_equality']='accepted state, rejected state, native history and failure unchanged'
            # Preserve all accepted common-time observations even for failed endpoints.
            endpoints[name]={}
            for t,s in zip(report['audit_times'],report['accepted_states'],strict=True):
                endpoints[name][round(t,12)]=(s,ini,c,report)
            write_json(output/'results.json',result)
    # All layouts share the same actual physical grid and weights.
    pre=json.loads((preparation/'results.json').read_text());ng=pre['global_count'];weights=np.array(pre['global_weights'])
    for layout in layouts:
        a,b=layout+'-0.002',layout+'-0.001';common=sorted(set(endpoints[a])&set(endpoints[b]));rows=[]
        for t in common:
            packets=[]
            for name in (a,b):
                s,ini,c,r=endpoints[name][t];s=np.asarray(s);o=ini['observation'];off=np.array(o['offsets'],int);values=[]
                for i,q in enumerate(ini['test_basis'][:ng]):
                    region=0 if o['grid'][i][2]<c['H'] else 1
                    values.append([row_value(q['velocity_k'],s,off[j])[0] for j in (0,1,2)]+[row_value(q['epsilon'],s,off[3])[0] if region else 0.,row_value(q['pressure'],s,off[4+region])[0]])
                ob=dict(time=t,grid=o['grid'][:ng],grid_sides=o['grid_sides'][:ng],grid_regions=o['grid_regions'][:ng],grid_values=values)
                packets.append(dict(observation=ob))
            rows.append(compare(*packets,weights,c))
        result['comparisons'][layout]=rows
    result['physical_completion']='pass' if all(r['status']=='pass' for r in result['cases'].values()) else 'fail'
    write_json(output/'results.json',result)
    return result

def recover(runtime,preparation,trajectories,output):
    # Read-only completion after the V1 collector's list/array error. No native
    # process, stored result, request, coefficient or history is overwritten.
    from ct_accepted2_audit import values
    request=json.loads((trajectories/'request-manifest.json').read_text())['request']
    require(request['runtime']['IncNavierStokesSolver']==runtime_identity(runtime/'IncNavierStokesSolver'),'foreign_recovery_runtime')
    require(request['preparation_sha256']==sha(preparation/'results.json'),'foreign_recovery_initialization')
    result=json.loads((trajectories/'work/results.json').read_text())
    require(set(result['cases'])=={x+'-'+str(dt) for x in request['layouts'] for dt in (.002,.001)},'incomplete_native_matrix_coverage')
    pr=json.loads((preparation/'results.json').read_text());n=pr['global_count'];weights=pr['global_weights'];reports={};initials={}
    for name,row in result['cases'].items():
        f=trajectories/'work'/name/'case.evolution.json'
        require(sha(f)==row['native_report_sha256'],'changed_native_matrix_report')
        reports[name]=json.loads(f.read_text());initials[name]=json.loads((preparation/'work'/name/'initial.json').read_text())
    for layout in request['layouts']:
        a,b=layout+'-0.002',layout+'-0.001';times=[];observed={}
        for name in (a,b):
            report=reports[name];observed[name]={}
            for t,s in zip(report['audit_times'],report['accepted_states'],strict=True):
                q=values(initials[name],s,n);q['observation']['time']=t;observed[name][round(t,12)]=q
        c=json.loads((trajectories/'work'/a/'construction.json').read_text())
        result['comparisons'][layout]=[compare(observed[a][t],observed[b][t],weights,c) for t in sorted(set(observed[a])&set(observed[b]))]
    result['physical_completion']='pass' if all(r['status']=='pass' for r in result['cases'].values()) else 'fail'
    result['trajectory_directory']=str(trajectories/'work')
    result['recovery']=dict(original_request_sha256=sha(trajectories/'request-manifest.json'),original_progress_sha256=sha(trajectories/'work/results.json'),fix='convert stored list to array in comparison only; all eight native trajectories reused without mutation',unmasked_load_label='V1 initial receipts named raw unmasked coefficient loads free_w; corrected source and final supplement use unmasked_w')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=['initial','advance','rough','recover'],required=True);p.add_argument('--preparation',type=Path);p.add_argument('--trajectories',type=Path);p.add_argument('--layouts',nargs='+',default=['B','X','Zplus','XZplus']);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    req=dict(schema='CTAccepted2DiagnosisV1',mode=a.mode,layouts=a.layouts,runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe')},previous=str(a.previous),source_sha256=sha(Path(__file__)),preparation_sha256=sha(a.preparation/'results.json') if a.preparation else None,recovered_request_sha256=sha(a.trajectories/'request-manifest.json') if a.trajectories else None)
    manifest=RequestManifest(a.output,req)
    # RequestManifest creates output, so subdirectories hold immutable work.
    result=recover(a.runtime,a.preparation,a.trajectories,a.output) if a.mode=='recover' else rough_control(a.runtime,a.output/'work') if a.mode=='rough' else initialize(a.runtime,a.previous,a.output/'work',a.layouts) if a.mode=='initial' else advance(a.runtime,a.previous,a.preparation/'work',a.output/'work',a.layouts)
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    print('diagnostics completed;',result.get('physical_completion',result.get('status')))
    return 0 if result.get('diagnostic_execution',result.get('status'))=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
