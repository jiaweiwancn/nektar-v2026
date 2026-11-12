#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT2 physical comparisons/budgets on accepted states, never repair a failure."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import require,run,probe,row_value
from ct_accepted_pilot import replay,compare,budget_classification
from ct_accepted_feasibility import write_json
from ct_accepted2_diagnosis import at,grids,make_case
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

def values(initial,state,n):
    o=initial['observation'];off=np.array(o['offsets'],int);s=np.array(state);vs=[];gs=[]
    for i,q in enumerate(initial['test_basis'][:n]):
        v=[];g=[]
        for b,key in ((0,'velocity_k'),(1,'velocity_k'),(2,'velocity_k'),(3,'epsilon'),(4+o['grid_regions'][i],'pressure')):
            a,d=row_value(q[key],s,off[b]) if q[key] is not None else (0.,np.zeros(3));v.append(a);g.append(d.tolist())
        vs.append(v);gs.append(g)
    return dict(observation=dict(time=0.,grid=o['grid'][:n],grid_sides=o['grid_sides'][:n],grid_regions=o['grid_regions'][:n],grid_values=vs,grid_gradients=gs))

def mass(initial,s,c):
    o=initial['observation'];off=np.array(o['offsets'],int);s=np.array(s);total=0.;faces={}
    for f in o['faces']:
        if f['interface']:continue
        u=np.array([row_value(f['velocity_k'],s,off[b])[0] for b in (0,1)]);area=np.array(f['area'])[[0,2]]
        axis=0 if area[0] else 2;side=-1 if (area[0] if axis==0 else area[1])<0 else 1
        value=f['weight']*float(u@area);key=str((axis,side));faces[key]=faces.get(key,0.)+value;total+=value
    return dict(outward_net=total,by_face=faces,physical_scale=c['Z'],scaled_defect=abs(total)/c['Z'],limit=1e-8,status='pass' if abs(total)/c['Z']<=1e-8 else 'fail')

def outlet_tractions(initial,state,xml):
    o=initial['observation'];s=np.array(state);off=np.array(o['offsets'],int);root=ET.parse(xml).getroot();rows=[];square=0.;weight=0.;undefined=0
    f=root.find("./CONDITIONS/FUNCTION[@NAME='OutletTotalTractionOuter']");expressions=[f.find("E[@VAR='"+v+"']").get('VALUE').replace('^','**') for v in ('u','v')]
    for q in o['faces']:
        if q['interface'] or not q['region'] or q['coordinate'][0] != max(x['upper'][0] for x in o['elements']):continue
        d=at(q,s,off);y=q['coordinate'][2];given=[float(eval(e,{'__builtins__':{}},dict(x=q['coordinate'][0],y=y,z=0.,t=0.,sqrt=np.sqrt,abs=np.abs,exp=np.exp,log=np.log))) for e in expressions]
        if d['status']=='admissible':
            T=np.array(d['total_stress']);actual=T[[0,2],0];defect=actual-np.array(given);square+=q['weight']*float(defect@defect);weight+=q['weight']
            rows.append(dict(coordinate=q['coordinate'],geometry=q['velocity_k']['geometry'],actual=actual.tolist(),prescribed=given,defect=defect.tolist(),Pi=d['Pi'],Qzz=d['Q'][2][2],Tzz=d['total_stress'][2][2]))
        else:undefined+=1;rows.append(dict(coordinate=q['coordinate'],geometry=q['velocity_k']['geometry'],status='undefined_actual_stress',invariants=d))
    return dict(status='all_evaluable' if not undefined else 'domain_limited',undefined_points=undefined,vector_traction_weighted_RMS=float(np.sqrt(square/weight)) if weight else None,rows=rows,scope='pointwise reconstructed outlet traction; prescribed total traction is enforced weakly, no pointwise equality imposed')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--preparation',type=Path,required=True);p.add_argument('--trajectories',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    prep=a.preparation/'work';results=json.loads((a.trajectories/'results.json').read_text());tra=Path(results.get('trajectory_directory',a.trajectories/'work'));pr=json.loads((a.preparation/'results.json').read_text())
    manifest=RequestManifest(a.output,dict(schema='CT2AcceptedStateAuditsV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),preparation_sha256=sha(a.preparation/'results.json'),trajectories_sha256=sha(a.trajectories/'results.json'),source_sha256=sha(Path(__file__)),limits=dict(temporal=1e-3,mass=1e-8,scalar=.02),geometry='actual x extent and native radial extent',frozen_controls=[[48,12,12,12,12,16],[64,16,16,16,16,16],[64,16,16,16,16,32]]))
    n=pr['global_count'];weights=pr['global_weights'];res=dict(schema='CT2PhysicalAuditsV1',diagnostic_execution='pass',cases={},representation_comparisons={});observed={}
    for name in results['cases']:
        folder=tra/name;report=json.loads((folder/'case.evolution.json').read_text());initial=json.loads((prep/name/'initial.json').read_text());c=json.loads((folder/'construction.json').read_text());row=dict(physical_completion=results['cases'][name]['status'],budgets={},native_report_sha256=sha(folder/'case.evolution.json'));res['cases'][name]=row
        state=report['state'];row['mass']=mass(initial,state,c);row['outlet_traction']=outlet_tractions(initial,state,folder/'case.xml')
        observed[name]={}
        for t,s in zip(report['audit_times'],report['accepted_states'],strict=True):
            packet=values(initial,s,n);packet['observation']['time']=t;observed[name][round(t,12)]=packet
        # Actual initialized state is distinct from the inferred continuum.
        start=values(initial,initial['observation']['state'],n);end=values(initial,state,n);start['observation']['time']=end['observation']['time']=report['time']
        row['drift_from_actual_initial']=compare(start,end,weights,c)
        g=np.array(end['observation']['grid_gradients']);w=np.array(weights);d=g[:,0,0]+g[:,1,2]
        row['pointwise_divergence']=dict(weighted_RMS=float(np.sqrt(w@(d*d)/w.sum())),maximum=float(np.max(abs(d))),signed_mean=float(w@d/w.sum()))
        if report['steps']>0:
            attempt=next(x for x in reversed(report['pimple_attempts']) if x['scheduled_complete']);snapshot=json.loads(Path(attempt['k_snapshot']).read_text())
            index=next(i for i,t in enumerate(report['audit_times']) if abs(t-report['time'])<1e-13);previous=report['accepted_states'][index-1]
            for label,orders,face in [('native',[32,8,10,10,10],16),('q48f16',[48,12,12,12,12],16),('q64f16',[64,16,16,16,16],16),('q64f32',[64,16,16,16,16],32)]:
                sub=a.output/name/label;sub.mkdir(parents=True);root=ET.parse(folder/'case.xml').getroot()
                for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
                write(root,sub/'case.xml');req=replay(initial,state,report['time'],[],export_physical_scalar_budget=True,previous_state=previous,stage_history=snapshot['history'],export_frozen_scalar_reconciliation=True,physical_face_order=face)
                b,e=probe(a.runtime,sub,req,'budget')
                row['budgets'][label]=dict(budget_classification(b,c),execution=e) if e['exit_code']==0 else dict(status='domain_or_replay_failure',packet=b,execution=e)
        write_json(a.output/'progress.json',res)
    for name in observed:
        if name.startswith('B-'):continue
        baseline='B-'+name.split('-',1)[1]
        if baseline not in observed:continue
        c=json.loads((tra/name/'construction.json').read_text())
        res['representation_comparisons'][name+'-vs-'+baseline]=[compare(observed[name][t],observed[baseline][t],weights,c) for t in sorted(set(observed[name])&set(observed[baseline]))]
    res['mass_status']='pass' if all(r['mass']['status']=='pass' for r in res['cases'].values()) else 'fail'
    res['scalar_accuracy_status']='pass' if all(r['budgets'] and all(b.get('accuracy_status')=='pass' and b.get('reconstruction_status')=='pass' for b in r['budgets'].values()) for r in res['cases'].values()) else 'failed_or_unresolved'
    write_json(a.output/'results.json',res);manifest.seal(a.output,['results.json']);print('audits',res['mass_status'],res['scalar_accuracy_status']);return 0
if __name__=='__main__':raise SystemExit(main())
