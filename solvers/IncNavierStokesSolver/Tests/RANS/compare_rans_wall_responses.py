#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Observe both wall responses on a common grid; wall differences are predictions."""
import argparse
import os
from pathlib import Path
import shutil
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from rans_duration_contract import require, MEMBERS, identity
from rans_request_manifest import RequestManifest, sha, runtime_identity
from rans_preservation_oracle import Reference, physical_grid
from qualify_rans_preservation import load, replay_request
from diagnose_rans_balanced_spatial import run_probe
from rans_qualification import write_json, exact_vector_hash


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--model',choices=('A93','A95'),required=True)
    p.add_argument('--accepted-prefix',required=True)
    p.add_argument('--ep-prefix',required=True)
    p.add_argument('--dt',type=float,required=True)
    p.add_argument('--time',type=float,default=.016)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=a.root.resolve();runtime=r/'runtime-resolution-v3'
    os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    refs={wall:Reference(a.model,'Natural',wall) for wall in ('Accepted','EP')}
    cases={};knots=set()
    for wall,name in (('Accepted',a.accepted_prefix),('EP',a.ep_prefix)):
        for member in MEMBERS:
            labels=(name+'-'+member,a.model+'-N-'+wall+'-duration-'+str(a.dt)+'-'+member)
            for label in labels:
                folder=r/'campaigns'/label;q=load(folder/'request-manifest.json')['request']
                identity(q,dict(model=a.model,policy='Natural',wall=wall))
                endpoint=folder/('T-'+format(a.time,'.4f'))
                require(RequestManifest(folder,q,True).completed(endpoint),'wall_response_endpoint_unsealed')
                require(abs(load(endpoint/'case.evolution.json')['time']-a.time)<1e-13,
                        'wall_response_endpoint_time_mismatch')
                knots.update(q['knots']);cases[label]=(q,folder,endpoint)
        knots.update(refs[wall].knots())
    grid,weights=physical_grid(sorted(knots),order=12,xorder=8)
    sides=[0]*len(grid)
    # Zero-weight trace observations retain explicit region/derivative provenance.
    for side in (-1,1):
        for x in (-.75,0.,.75):
            grid.append([x,0.,2.]);weights.append(0.);sides.append(side)
    request=dict(schema='CommonGridWallResponseV1',model=a.model,policy='Natural',time=a.time,
        dt=a.dt,prefixes=dict(Accepted=a.accepted_prefix,EP=a.ep_prefix),grid=grid,weights=weights,
        sides=sides,knots=sorted(knots),orders=[12,8],runtime=runtime_identity(runtime/'RANSPreservationProbe'),
        campaigns={label:sha(folder/'request-manifest.json') for label,(q,folder,endpoint) in cases.items()},
        interpretation='same-wall uncertainty is gated; differences between wall models are measured without an equality requirement')
    seal=RequestManifest(a.output,request,False);observations={};provenance=[]
    for label,(q,folder,endpoint) in cases.items():
        target=a.output/label;target.mkdir();shutil.copy2(endpoint/'case.xml',target/'case.xml')
        native=load(endpoint/'case.evolution.json');seed=load(folder/'initial-packet.json')
        replay=replay_request(seed,grid,state=native['state'],time=native['time'])
        replay.update(fixed_boundary=native['initialization_receipt']['live_fixed_boundary'],
                      grid_sides=sides,export_constitutive=True)
        packet,seconds=run_probe(runtime,target,replay);o=packet['observation']
        require(o['grid']==grid and o['grid_sides']==sides,'wall_response_grid_or_side_mismatch')
        require(exact_vector_hash(o['state'])==exact_vector_hash(native['state']),
                'wall_response_replay_changed_state')
        observations[label]=packet
        provenance.append(dict(label=label,request_sha256=sha(folder/'request-manifest.json'),
            native_report_sha256=sha(endpoint/'case.evolution.json'),packet_sha256=sha(target/'packet.json.gz'),
            state_sha256=exact_vector_hash(native['state']),seconds=seconds))
    weights=np.asarray(weights);outer=np.asarray(grid)[:,2]>2
    outer |= (np.asarray(grid)[:,2]==2)&(np.asarray(sides)==1)
    def rms(x,w):return float(np.sqrt(np.sum(w*x*x)/np.sum(w)))
    responses={};fields={};all_pass=True
    for wall,name in (('Accepted',a.accepted_prefix),('EP',a.ep_prefix)):
        ref=refs[wall];arrays={}
        for member in MEMBERS:
            pp=observations[name+'-'+member]
            bb=observations[a.model+'-N-'+wall+'-duration-'+str(a.dt)+'-'+member]
            delta=np.asarray(pp['observation']['grid_values'])-np.asarray(bb['observation']['grid_values'])
            require(delta.shape==(len(grid),5),'wall_response_field_shape')
            delta[:,4]-=np.dot(weights,delta[:,4])/weights.sum()
            value={key:(delta[outer,j],weights[outer]) if j==3 else (delta[:,j],weights)
                for j,key in enumerate(('u','w','k','epsilon_outer','Pi'))}
            scales=dict(u=ref.us,w=ref.us,k=ref.K,epsilon_outer=ref.A/ref.H,Pi=ref.us**2)
            for key,(v,w) in list(value.items()):value[key]=(v/scales[key],w)
            ps,bs=pp['native_constitutive']['samples'],bb['native_constitutive']['samples']
            require(len(ps)==len(bs),'wall_response_tensor_shape')
            stress=[];sw=[]
            for first,last in zip(ps,bs,strict=True):
                require(first['sample_index']==last['sample_index'] and first['side']==last['side'],
                        'wall_response_tensor_identity')
                R=np.asarray(first['R'])-np.asarray(last['R']);Q=np.asarray(first['Q'])-np.asarray(last['Q'])
                stress.append([R[0,0],R[1,1],R[2,2],R[0,2],Q[0,2],first['P']-last['P']])
                sw.append(weights[first['sample_index']])
            stress=np.asarray(stress);sw=np.asarray(sw)
            for j,key in enumerate(('Rxx','Ryy','Rzz','Rxz','Qxz','P')):
                value[key]=(stress[:,j]/(ref.A/ref.H if key=='P' else ref.us**2),sw)
            arrays[member]=value
        responses[wall]=arrays;wallfields={}
        for key in arrays['B5']:
            v=[arrays[m][key][0] for m in MEMBERS];w=arrays['B5'][key][1]
            require(all(np.array_equal(w,arrays[m][key][1]) for m in MEMBERS),'wall_response_weights_mismatch')
            signal=max(rms(x,w) for x in v);dh=rms(v[0]-v[1],w);dt=max(rms(v[0]-v[2],w),rms(v[1]-v[3],w))
            hlimit=max(5e-7,.1*signal);tlimit=max(1e-7,.05*signal)
            passed=dh<=hlimit and dt<=tlimit;all_pass &= passed
            wallfields[key]=dict(response_RMS=[rms(x,w) for x in v],spatial_difference=dh,
                temporal_difference=dt,spatial_limit=hlimit,temporal_limit=tlimit,
                uncertainty_status='pass' if passed else 'unresolved',
                signal_status='resolved' if signal>5e-7 else 'unresolved_floor')
        fields[wall]=wallfields
    wall_differences={}
    for key in responses['Accepted']['B5']:
        values=[]
        for member in MEMBERS:
            x,w=responses['Accepted'][member][key];y,v=responses['EP'][member][key]
            require(np.array_equal(w,v) and x.shape==y.shape,'cross_wall_sample_mismatch')
            values.append(rms(x-y,w))
        wall_differences[key]=dict(scaled_weighted_RMS=values,
            interpretation='difference between baseline-subtracted wall-model responses; no equality or convergence requirement across wall laws')
    result=dict(status='pass' if all_pass else 'unresolved',configuration=dict(model=a.model,policy='Natural',wall='EP-vs-Accepted'),
        time=a.time,grid_request_sha256=sha(a.output/'request-manifest.json'),fields=fields,
        wall_differences=wall_differences,provenance=provenance,
        scope='common physical grid, explicit two-sided interface observations, one global pressure gauge; prescribed-reference patch only')
    write_json(a.output/'results.json',result)
    seal.seal(a.output,['results.json',*map(str,a.output.glob('*/packet.json.gz'))])
    return 0 if all_pass else 1


if __name__=='__main__':raise SystemExit(main())
