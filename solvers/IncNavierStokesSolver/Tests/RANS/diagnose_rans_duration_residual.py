#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only localization and timestep sensitivity of a duration residual gate.

The original gate is reapplied without changing its floor. This supplement
compares residual *loads*, using one fixed common mass, not differences of norms.
It never admits a failed duration cohort or advances a state.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from rans_duration_contract import require, identity, residual_reduction, MEMBERS
from rans_preservation_oracle import Reference
from rans_request_manifest import RequestManifest, sha
from rans_spatial_resolution import basis_matrix, scales
from rans_qualification import write_json


def load(path):
    if str(path).endswith('.gz'):
        with gzip.open(path,'rt') as f:return json.load(f)
    return json.loads(path.read_text())


def common_vertical(common, packet, cross, analysis, ref):
    co=common['observation'];o=packet['observation'];r=packet['resolution']
    n=o['sizes'][1];offset=o['offsets'][1];rows=o['volume']
    fine=basis_matrix([q['velocity_k'] for q in rows],n)
    test=basis_matrix([q['velocity_k'] for q in cross['test_basis']],co['sizes'][1])
    weights=np.asarray([q['weight'] for q in rows])
    P=cho_solve(cho_factor(np.asarray(r['mass'][1]),lower=True),
                (fine.T@test.multiply(weights[:,None])).toarray())
    original=analysis['common_test']['fields']['w'];free=np.asarray(original['common_free_ids'],int)
    M=np.asarray(common['resolution']['mass'][1])[np.ix_(free,free)]
    require(sha_bytes(M)==original['common_mass_sha256'],'fixed_mass_identity_mismatch')
    p=P[:,free];scale=scales(ref)[1]*np.sqrt(8.);factor=cho_factor(M,lower=True)
    norm=lambda v:float(np.sqrt(max(0.,v@cho_solve(factor,v)))/scale)
    project=lambda v:p.T@np.asarray(v)[offset:offset+n]
    terms={k:project(r['terms'][k]) for k in r['physical_term_names']}
    total=project(r['native']['spatial_unmasked'])
    require(np.max(np.abs(total-original['restricted_load']))<2e-11,'stored_common_load_changed')
    require(np.max(np.abs(sum(terms.values())-total))<2e-11,'common_term_reconstruction_failed')
    local=[];elements={e['geometry']:e for e in o['elements']}
    for key,et in r['element_terms'].items():
        geometry=int(key.split(':')[1]);e=elements[geometry]
        vectors={k:project(v) for k,v in et.items() if k in terms}
        value=sum(vectors.values(),np.zeros(len(free)))
        local.append(dict(element=key,interval=[e['lower'][2],e['upper'][2]],
            physical_bounds=[e['lower'],e['upper']],total_common_mass_dual=norm(value),
            physical_terms={k:norm(v) for k,v in vectors.items()}))
    support={i:[] for i in free}
    for e in co['elements']:
        ids=set()
        for q in co['volume']:
            if q['velocity_k']['geometry']==e['geometry']:ids.update(q['velocity_k']['ids'])
        for i in ids.intersection(support):support[i].append(dict(geometry=e['geometry'],
            interval=[e['lower'][2],e['upper'][2]]))
    ranking=sorted(range(len(free)),key=lambda j:abs(total[j])/np.sqrt(M[j,j]),reverse=True)[:20]
    body=dict(total=norm(total),common_free_ids=free.tolist(),common_mass_sha256=sha_bytes(M),
        physical_term_norms={k:norm(v) for k,v in terms.items()},
        term_sum_reconstruction_maximum=float(np.max(np.abs(sum(terms.values())-total))),
        term_cancellation_ratio=norm(total)/sum(norm(v) for v in terms.values()),
        dominant_elements=sorted(local,key=lambda q:q['total_common_mass_dual'],reverse=True)[:16],
        dominant_rows=[dict(common_id=int(free[j]),load=float(total[j]),
            diagonal_mass_scaled=float(total[j]/np.sqrt(M[j,j])/scale),support=support[free[j]]) for j in ranking])
    return body,total,norm


def sha_bytes(array):
    import hashlib
    return hashlib.sha256(np.asarray(array,dtype='<f8').tobytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--admission',type=Path,required=True)
    p.add_argument('--cohort',required=True);p.add_argument('--time',type=float,default=.5)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true')
    a=p.parse_args();ad=load(a.admission);cfg=ad['configuration'];ref=Reference(cfg['model'],cfg['policy'],cfg['wall'])
    cp_path=Path(ad['common_directory'])/'packet.json.gz';require(sha(cp_path)==ad['common_packet_sha256'],'common_packet_changed')
    common=load(cp_path);paths={};analyses={}
    for member in MEMBERS:
        folder=a.root/'evolved-audits'/(a.cohort+'-'+str(a.time)+'-'+member)
        packet=folder/('T-'+format(a.time,'.4f'))/'packet.json.gz'
        cross=folder/('T-'+format(a.time,'.4f')+'-common-cross')/'packet.json.gz'
        analysis=packet.parent/'analysis.json';paths[member]=dict(packet=packet,cross=cross,analysis=analysis)
        row=load(analysis);identity(row['configuration'],cfg);require(abs(row['time']-a.time)<1e-13,'wrong_time')
        require(row['status']=='pass','source_analysis_failed');analyses[member]=row
    request=dict(schema='DurationResidualFailureDiagnosisV1',configuration=cfg,time=a.time,
        admission_sha256=sha(a.admission),source_sha256={m:{k:sha(v) for k,v in q.items()} for m,q in paths.items()},
        common_sha256=sha(cp_path),scope='read-only; retain original gates; compare residual loads on fixed tests',
        diagnostic_uncertainty=dict(residual_floor=1e-8,accumulation=1e-10,maximum_ratio=.85,
            temporal_load_floor=1e-10,temporal_fraction_of_spatial_load_difference=.1))
    seal=RequestManifest(a.output,request,a.resume)
    if seal.completed(a.output):return 0 if load(a.output/'results.json')['original_duration_gate']=='pass' else 1
    rows={};loads={};norm=None
    for member in MEMBERS:
        body,total,norm=common_vertical(common,load(paths[member]['packet']),load(paths[member]['cross']),analyses[member],ref)
        body.update(source_native_sha256=analyses[member]['native_report_sha256'],history=analyses[member]['history'],
            native_full=analyses[member]['spatial']['native_full_mass_dual'],
            reconstruction=analyses[member]['spatial']['reconstruction_defects'])
        rows[member]=body;loads[member]=total
    reductions={dt:residual_reduction(analyses[c],analyses[f],cfg['model'])
        for dt,c,f in (('primary','B5','B11'),('half','B5-dt2','B11-dt2'))}
    spatial=norm(loads['B5']-loads['B11']);temporal={m:norm(loads[m]-loads[m+'-dt2']) for m in MEMBERS[:2]}
    limit=max(1e-10,.1*spatial)
    result=dict(status='pass' if reductions['primary']['status']=='pass' else 'failed_or_unresolved',
        execution_status='completed',configuration=cfg,time=a.time,original_duration_gate=reductions['primary']['status'],
        original_gate_unchanged=True,reductions=reductions,vertical_common_load_comparison=dict(
            spatial_difference=spatial,timestep_differences=temporal,uncertainty_limit=limit,
            status='resolved' if max(temporal.values())<=limit else 'temporally_unresolved',
            interpretation='diagnostic load sensitivity; does not replace the original residual-reduction gate'),
        members=rows,scope=request['scope'])
    write_json(a.output/'results.json',result);seal.seal(a.output,['results.json'])
    print(result['status'],result['vertical_common_load_comparison'])
    return 0 if result['original_duration_gate']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
