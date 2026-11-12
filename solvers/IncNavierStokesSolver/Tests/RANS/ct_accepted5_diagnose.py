#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reconcile CT4 native loads and attribute tractions, without advancing a field."""
import argparse,json,os,shutil
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from ct_accepted4_reference import Column,reference,write_json
from ct_accepted4_run import expected
from ct_accepted4_analysis import weak
from ct_accepted5_traction import profiles,delta
from ct_accepted_pilot import replay,artifact_hash
from ct_accepted_open import probe,require,row_value
from ct_accepted2_storage import artifact
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha,runtime_identity

def add(a,row,block,off,value):a[np.array(row['ids'],int)+off[block]]+=value

def reference_loads(o,c):
    """Independent continuous stress/advection on the *same* native rows."""
    n=len(o['state']);off=o['offsets'];col=Column(c['model'],c['K0'],c['H'],c['roughness_input'],c['nu'],order=96)
    at=lru_cache(None)(lambda z,outer:reference(col,z,outer));elements={};pressures={};external={}
    for q in o['volume']:
        region=bool(q['region']);z=q['coordinate'][2];v=at(z,region);G=np.zeros((3,3));G[0,2]=v[4];key=('outer:' if region else 'inner:')+str(q['velocity_k']['geometry'])
        a=elements.setdefault(key,np.zeros(n));p=pressures.setdefault(key,np.zeros(n));d=np.array(q['velocity_k']['derivative']);w=q['weight']
        if region:T=-oracle(c['model'],v[1],v[2],G)['R']+np.eye(3)*(2*v[1]/3-v[3])
        else:T=v[8]*(G+G.T)-np.eye(3)*v[3]
        for b,axis in enumerate((0,2)):
            flux=T[axis].copy();flux[0]-=v[0]*(v[0] if axis==0 else 0.)
            add(a,q['velocity_k'],b,off,w*(flux@d));add(p,q['velocity_k'],b,off,-w*v[3]*d[axis])
    for q in o['faces']:
        if q['interface'] or q['coordinate'][0]!=c['x'][-1]:continue
        reg=bool(q['region']);v=at(q['coordinate'][2],reg);key=('outer:' if reg else 'inner:')+str(q['velocity_k']['geometry'])
        a=elements.setdefault(key,np.zeros(n));ext=external.setdefault(key,np.zeros(n));phi=np.array(q['velocity_k']['phi']);tr=c['reference_column']['outer_equilibrium']['outlet_traction'];w=q['weight']
        for b,value in enumerate((v[0]**2-(tr[0] if reg else 0.),-tr[1])):
            add(a,q['velocity_k'],b,off,w*value*phi);add(ext,q['velocity_k'],b,off,w*value*phi)
    return elements,pressures,external

def reconcile(packet,c):
    o=packet['observation'];r=packet['resolution'];off=o['offsets'];s=np.array(o['state']);n=len(s);names=r['physical_term_names'];catalog={q['geometry']:q for q in o['elements']}
    volume={key:sum((np.array(v) for k,v in terms.items() if k in names),np.zeros(n)) for key,terms in r['element_terms'].items()}
    ext={};open_terms={k:np.zeros(n) for k in names};pressure={key:np.zeros(n) for key in volume}
    for q in o['volume']:
        region=q['region'];key=('outer:' if region else 'inner:')+str(q['velocity_k']['geometry']);pi=row_value(q['pressure'],s,off[4+region])[0]
        for b,axis in enumerate((0,2)):add(pressure[key],q['velocity_k'],b,off,-q['weight']*pi*np.array(q['velocity_k']['derivative'][axis]))
    for q in o['faces']:
        if q['interface'] or q['coordinate'][0]!=c['x'][-1]:continue
        region=q['region'];row=q['velocity_k'];key=('outer:' if region else 'inner:')+str(row['geometry']);a=ext.setdefault(key,np.zeros(n));w=q['weight'];phi=np.array(row['phi']);u=[row_value(row,s,off[b])[0] for b in (0,1)]
        tr=c['reference_column']['outer_equilibrium']['outlet_traction']
        for b in range(4):
            test=row if b<3 else q['epsilon']
            if test is None:continue
            value=row_value(test,s,off[b])[0];v=w*u[0]*value*np.array(test['phi']);add(open_terms['conservative_advection'],test,b,off,v);add(a,test,b,off,v)
        name='outer_momentum_stress' if region else 'inner_momentum_stress'
        for b,t in enumerate((tr[0] if region else 0.,tr[1])):
            v=-w*t*phi;add(open_terms[name],row,b,off,v);add(a,row,b,off,v)
    defects={}
    for name in names:
        vv=sum((np.array(q.get(name,np.zeros(n))) for q in r['element_terms'].values()),np.zeros(n))
        defects[name]=float(np.max(abs(vv+open_terms[name]-r['terms'][name])))
    total=sum(volume.values(),np.zeros(n))+sum(ext.values(),np.zeros(n));native=np.array(r['native']['spatial_unmasked'])
    require(max(defects.values())<2e-10 and np.max(abs(total-native))<2e-10,'volume_plus_open_face_native_reconstruction')
    ref,refp,refext=reference_loads(o,c);masked=set(int(i) for i,_ in r['constraints']);free={b:[j for j in range(o['sizes'][b]) if off[b]+j not in masked] for b in (0,1)}
    actual=total.copy();expected=sum(ref.values(),np.zeros(n));direction=actual-expected;dual={};energy=0.
    for b in (0,1):
        ids=np.array(free[b],int);mat=np.array(r['mass'][b]);factor=cho_factor(mat[np.ix_(ids,ids)]);dual[b]=cho_solve(factor,direction[ids+off[b]]);energy+=float(direction[ids+off[b]]@dual[b])
    def signed(v):return sum(float(v[np.array(free[b])+off[b]]@dual[b]) for b in (0,1))
    slabs={};element_rows=[]
    for key,v in volume.items():
        geom=int(key.split(':')[1]);e=catalog[geom];lo,hi=e['lower'][2],e['upper'][2];a=v+ext.get(key,np.zeros(n))-ref[key]
        block={b:float(np.linalg.norm(a[np.array(free[b])+off[b]])) for b in (0,1)}
        label=str(lo)+':'+str(hi);q=slabs.setdefault(label,dict(lower=lo,upper=hi,signed_share=0.,geometries=[]));q['signed_share']+=signed(a);q['geometries'].append(geom)
        element_rows.append(dict(geometry=geom,region=key.split(':')[0],lower=e['lower'],upper=e['upper'],signed_dual_projection=signed(a),
            pressure_difference_projection=signed(pressure[key]-refp[key]),open_difference_projection=signed(ext.get(key,np.zeros(n))-refext.get(key,np.zeros(n))),
            free_row_euclidean_diagnostic=block,scope='volume+open-face represented-minus-continuous, not isolated raw stress norm'))
    ordered=sorted(slabs.values(),key=lambda q:q['signed_share'],reverse=True)
    for q in ordered:q['fraction_of_global_dual_squared']=q['signed_share']/energy
    require(abs(sum(q['signed_share'] for q in ordered)-energy)<2e-10*max(1.,energy),'signed_localization_sum')
    pressure_total=sum(pressure.values(),np.zeros(n));ref_pressure=sum(refp.values(),np.zeros(n));terms={}
    for name in ('inner_momentum_stress','outer_momentum_stress','conservative_advection','roughness_force'):
        v=np.array(r['terms'][name]);terms[name]=dict(signed_projection=signed(v),native_free_maximum=float(max(np.max(abs(v[np.array(free[b])+off[b]])) for b in (0,1))))
    tops=[]
    for b in (0,1):
        ids=sorted(free[b],key=lambda j:abs(direction[off[b]+j]),reverse=True)[:8]
        for j in ids:
            support=sorted({q['velocity_k']['geometry'] for q in o['volume'] if j in q['velocity_k']['ids']})
            tops.append(dict(block=b,global_row=off[b]+j,local_row=j,signed_defect=direction[off[b]+j],native_load=native[off[b]+j],continuous_load=expected[off[b]+j],support_geometries=support,
                support_bounds=[dict(geometry=g,lower=catalog[g]['lower'],upper=catalog[g]['upper']) for g in support]))
    return dict(status='pass',physical_reconstruction_maximum=float(np.max(abs(total-native))),term_reconstruction=defects,global_momentum_dual_squared=energy,
        radial_slabs=ordered,elements=element_rows,decisive_free_rows=tops,native_norms=r['norms']['native_stage'],momentum_terms=terms,
        pressure_volume_projection=signed(pressure_total-ref_pressure),continuous_free_weak_maximum=float(max(np.max(abs(expected[np.array(free[b])+off[b]])) for b in (0,1))),
        prescribed_row_reaction_maximum=float(np.max(abs(native[list(masked)]))),
        shared_interface_rows=[dict(block=b,row=off[b]+j,load=native[off[b]+j]) for b in (0,1) for j in free[b] if any(q['interface'] and j in q['velocity_k']['ids'] for q in o['faces'])],
        owner_mismatch_count=sum(q['velocity_k']['geometry']!=q['pressure']['geometry'] or (q['epsilon'] is not None and q['epsilon']['geometry']!=q['velocity_k']['geometry']) for q in o['volume']),
        scope='unmasked volume+independent open flux reconstructed before constraints; shared interface rows once, actual pointwise traction and prescribed reactions separate; mass-dual signed localization is diagnostic only')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct4','runtime','output','contract'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');root=a.ct4
    iz=root/'initial-v3/A95-6.0-P4/0.001';run=root/'evolution-v1/A95-6.0-P4/T0.008/0.001'
    prep=json.loads((iz.parent/'results.json').read_text());c=prep['case'];receipt=json.loads((run.parent/'results.json').read_text());native=json.loads((run/'case.evolution.json').read_text());ini=json.loads((iz/'initial.json').read_text())
    require(runtime_identity(a.runtime/'IncNavierStokesSolver')==prep['members']['0.001']['identity']['runtime'],'foreign_governing_runtime')
    require(sha(run/'case.evolution.json')==receipt['members']['0.001']['native_report_sha256'] and artifact_hash(native['checkpoint'])==receipt['members']['0.001']['history_sha256'],'changed_CT4_report_or_history')
    m=RequestManifest(a.output,dict(schema='CT5NativeDiagnosisV1',contract_sha256=sha(a.contract),baseline_sha256=sha(run.parent/'results.json'),initial_sha256=sha(iz/'initial.json'),
        runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source_sha256={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted5_diagnose.py','ct_accepted5_traction.py')}))
    result=dict(status='fail',accepted={},stages=[],governing_evolution_runs=0,baseline_native_sha256=sha(run/'case.evolution.json'))
    for t in (0.,.001,.002,.004,.008):
        idx=next(i for i,v in enumerate(native['audit_times']) if abs(v-t)<1e-13);state=native['accepted_states'][idx];folder=a.output/str(t);folder.mkdir();shutil.copy2(iz/'case.xml',folder/'case.xml')
        packet,e=probe(a.runtime,folder,replay(ini,state,t,[],export_basis=True,export_resolution=True),'resolution')
        require(e['exit_code']==0 and packet['observation']['state']==state,'read_only_replay_changed_input')
        result['accepted'][str(t)]=dict(state_sha256=exact_vector_hash(state),packet_sha256=sha(folder/'resolution.json'),execution=e,reconciliation=reconcile(packet,c),
            analytical_weak=weak(packet,c),traction=profiles(ini,state,c))
        write_json(a.output/'progress.json',result)
    for attempt in native['pimple_attempts']:
        if attempt['time'] not in (.001,.002,.004,.008):continue
        fpath=artifact(attempt['flow_snapshot'],attempt['flow_snapshot_sha256']);kpath=artifact(attempt['k_snapshot'],attempt['k_snapshot_sha256']);f=json.loads(fpath.read_text());k=json.loads(kpath.read_text())
        states=[('input',f['input']),('momentum_predictor',f['momentum_predictor']),('corrected_flow',f['corrected_flow']),('post_flow',k['post_flow']),('post_epsilon',k['pre_k']),('post_k',k['corrected'])]
        entries={name:profiles(ini,s,c,f['input']) for name,s in states};parts={}
        for (left,_),(right,_) in zip(states[:-1],states[1:]):parts[left+'->'+right]=delta(entries[left],entries[right],c)
        parts['input->post_k']=delta(entries['input'],entries['post_k'],c)
        result['stages'].append(dict(time=attempt['time'],incoming_sha256=attempt['incoming_sha256'],history_sha256=attempt['history_sha256'],flow_snapshot_sha256=sha(fpath),k_snapshot_sha256=sha(kpath),
            traction=entries,decomposition=parts,pressure_action_maximum=float(np.max(abs(f['full_pressure_action']))),predictor_residual_maximum=float(np.max(abs(f['predictor_linear_residual']))),pressure_residual_maximum=float(np.max(abs(f['pressure_linear_residual']))),
            scope='no stress-domain rejection added at diagnostic predictor; each actual D uses that stage endpoint state, frozen-input D shown only as counterfactual'))
    require(all(v['reconciliation']['owner_mismatch_count']==0 for v in result['accepted'].values()),'volume_owner_inconsistency')
    require(result['accepted']['0.0']['analytical_weak']['status']=='pass','initial_continuous_weak_identity')
    initial=result['accepted']['0.0']['reconciliation'];winner=initial['radial_slabs'][0]
    target=(winner['lower']+winner['upper'])/2
    gradient=sum(abs(q['decomposition']['momentum_predictor->corrected_flow']['gradient'])+abs(q['decomposition']['input->momentum_predictor']['gradient']) for q in result['stages'])
    coefficient=sum(abs(q['decomposition']['input->post_k']['coefficient'])+abs(q['decomposition']['input->post_k']['cross']) for q in result['stages'])
    require(winner['signed_share']>0 and gradient>coefficient,'no_supported_gradient_split_selection')
    result.update(status='pass',selection=dict(selected_interval=[winner['lower'],winner['upper']],point=target,point_mm=40*target,
        proposed_first_slab_selected=winner['lower']==0 and winner['upper']==c['H']/32,signed_share=winner['signed_share'],fraction=winner['fraction_of_global_dual_squared'],
        gradient_attribution_magnitude=gradient,coefficient_attribution_magnitude=coefficient,case=c,
        scope='one maximal positive assembled defect share; selected before constructing a new space/history; no mesh search'))
    require(sha(run/'case.evolution.json')==result['baseline_native_sha256'],'postmortem_changed_CT4_native')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT5 selected interval',result['selection']['selected_interval'],'fraction',result['selection']['fraction']);return 0
if __name__=='__main__':raise SystemExit(main())
