#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only Stokes/element-ownership witness on frozen native quadrature.

The counterfactual changes exported test-basis ownership only. It does not
change production assembly, a solution coefficient, a boundary, or a history.
An explained native failure still returns nonzero and remains unqualified.
"""
import argparse
import gzip
import json
import os
import shutil
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from diagnose_rans_balanced_spatial import run_probe
from qualify_rans_preservation import replay_request
from rans_duration_contract import require
from rans_preservation_oracle import Reference, weak_balance
from rans_qualification import exact_vector_hash, write_json
from rans_request_manifest import RequestManifest, runtime_identity, sha


def load(path):
    with gzip.open(path,'rt') as f:return json.load(f)


def stokes(observation,ref,resolution):
    n=observation['sizes'][0];s=np.zeros((3,n),dtype=np.longdouble)
    for q in observation['volume']:
        row=q['velocity_k'];ids=np.asarray(row['ids'],int)
        s[:,ids]+=np.longdouble(q['weight'])*np.asarray(row['derivative'],np.longdouble)
    for q in observation['faces']:
        row=q['velocity_k'];ids=np.asarray(row['ids'],int)
        s[:,ids]-=np.longdouble(q['weight'])*np.asarray(q['area'],np.longdouble)[:,None]*np.asarray(row['phi'],np.longdouble)
    flux=np.asarray(ref.t0*s[0]-ref.pi[0]*s[2],float)
    offset=observation['offsets'][1];fixed={int(i)-offset for i,v in resolution['constraints'] if offset<=i<offset+n}
    free=np.asarray([i for i in range(n) if i not in fixed],int);mass=np.asarray(resolution['mass'][1])[np.ix_(free,free)]
    dual=float(np.sqrt(max(0.,flux[free]@cho_solve(cho_factor(mass,lower=True),flux[free])))/(ref.us**2/ref.H*np.sqrt(8.)))
    return dict(maximum_stokes_defect_by_axis=np.max(np.abs(s),axis=1).astype(float).tolist(),
        vertical_constant_flux_maximum=float(np.max(np.abs(flux))),
        vertical_constant_flux_scaled_native_free_mass_dual=dual,
        accumulation='long double on the exported double weights/basis',
        dominant_rows=[dict(row=int(i),constant_flux_load=float(flux[i]))
            for i in np.argsort(np.abs(flux))[-20:][::-1]])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--runtime',type=Path,required=True);p.add_argument('--model',choices=('A93','A95','E94E95'),required=True)
    p.add_argument('--policy',choices=('Natural','Prescribed'),default='Natural');p.add_argument('--wall',choices=('Accepted','EP'),required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    runtime=a.runtime.resolve();os.environ['LD_LIBRARY_PATH']=str(runtime/'lib')
    packet=load(a.packet);o=packet['observation'];ref=Reference(a.model,a.policy,a.wall)
    require(o['resolved_descriptor'].startswith(a.model+':'),'ownership_model_mismatch')
    require(':policy='+a.policy+':' in o['resolved_descriptor'],'ownership_policy_mismatch')
    require(':wall='+('Accepted' if a.wall=='Accepted' else 'EP') in o['resolved_descriptor'],'ownership_wall_mismatch')
    request=dict(schema='FrozenVolumeOwnershipDiagnosisV1',source_packet_sha256=sha(a.packet),
        source_configuration_sha256=sha(a.packet.parent/'case.xml'),model=a.model,policy=a.policy,wall=a.wall,
        descriptor=o['resolved_descriptor'],state_sha256=exact_vector_hash(o['state']),time=o['time'],
        runtime=runtime_identity(runtime/'RANSPreservationProbe'),diagnostic_source_sha256=sha(Path(__file__)),
        scope='same points/weights/coefficients; explicit containing-element test basis only; no production or trajectory modification')
    seal=RequestManifest(a.output,request,a.resume)
    if seal.completed(a.output):return 0 if json.loads((a.output/'results.json').read_text())['native_qualification_status']=='pass' else 1
    elements={q['geometry']:q for q in o['elements']};owners=[];wrong=[]
    for i,q in enumerate(o['volume']):
        current=q['velocity_k']['geometry'];e=elements[current];x=q['coordinate']
        outside=max(max(e['lower'][axis]-x[axis],x[axis]-e['upper'][axis],0.) for axis in (0,2))
        if outside:
            matches=[e['geometry'] for e in elements.values() if all(e['lower'][axis]<=x[axis]<=e['upper'][axis] for axis in (0,2))]
            require(len(matches)==1,'physical_point_has_no_unique_containing_element');owner=matches[0]
            wrong.append(dict(index=i,coordinate=x,selected_geometry=current,containing_geometry=owner,
                selected_interval=[e['lower'][2],e['upper'][2]],containing_interval=[elements[owner]['lower'][2],elements[owner]['upper'][2]],
                outside_distance=outside,weight=q['weight']))
        else:owner=current
        owners.append(owner)
    native_weak=weak_balance(ref,o);native_stokes=stokes(o,ref,packet['resolution'])
    target=a.output/'explicit-containing-elements';target.mkdir();shutil.copy2(a.packet.parent/'case.xml',target/'case.xml')
    probe=replay_request(packet,[q['coordinate'] for q in o['volume']],state=o['state'],time=o['time'])
    probe.update(export_test_basis=True,grid_test_geometries=owners,
        grid_sides=[(-1 if elements[g]['upper'][2]<=ref.H else 1) if abs(q['coordinate'][2]-ref.H)<1e-12 else 0 for q,g in zip(o['volume'],owners,strict=True)])
    output,seconds=run_probe(runtime,target,probe)
    require(output['performed_time_steps']==0,'ownership_diagnostic_advanced_state')
    require(exact_vector_hash(output['observation']['state'])==request['state_sha256'],'ownership_diagnostic_changed_state')
    require(output['initialization_receipt']['explicit_boundary_sha256']==packet['initialization_receipt']['explicit_boundary_sha256'],'ownership_diagnostic_changed_boundary')
    rows=output['test_basis'];require(len(rows)==len(owners),'ownership_basis_coverage')
    volume=[]
    for q,row,owner in zip(o['volume'],rows,owners,strict=True):
        require(row['velocity_k']['geometry']==owner,'explicit_owner_not_consumed')
        volume.append(dict(q,velocity_k=row['velocity_k'],pressure=row['pressure'],epsilon=row['epsilon']))
    corrected=dict(o,volume=volume);control_weak=weak_balance(ref,corrected);control_stokes=stokes(corrected,ref,packet['resolution'])
    result=dict(status='diagnosed_native_quadrature_failure' if native_weak['status']!='pass' and control_weak['status']=='pass' else 'unresolved' if native_weak['status']!='pass' else 'pass',
        execution_status='completed',native_qualification_status=native_weak['status'],
        configuration={k:request[k] for k in ('model','policy','wall','descriptor')},time=o['time'],source_packet_sha256=sha(a.packet),
        state_sha256=request['state_sha256'],boundary_sha256=packet['initialization_receipt']['explicit_boundary_sha256'],
        wrong_selected_points=len(wrong),maximum_outside_distance=max((q['outside_distance'] for q in wrong),default=0.),
        total_wrong_point_weight=sum(q['weight'] for q in wrong),ownership_witnesses=wrong,
        original_native=dict(stokes=native_stokes,analytical_weak_balance=native_weak),
        explicit_containing_element_counterfactual=dict(stokes=control_stokes,analytical_weak_balance=control_weak,
            source_packet_sha256=sha(target/'packet.json.gz'),seconds=seconds),
        trajectory_requalified=False,production_geometry_policy_unchanged=True,scope=request['scope'])
    write_json(a.output/'results.json',result);seal.seal(a.output,['results.json','explicit-containing-elements/packet.json.gz'])
    print(result['status'],result['wrong_selected_points'],native_stokes['vertical_constant_flux_maximum'],control_stokes['vertical_constant_flux_maximum'])
    return 0 if result['native_qualification_status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
