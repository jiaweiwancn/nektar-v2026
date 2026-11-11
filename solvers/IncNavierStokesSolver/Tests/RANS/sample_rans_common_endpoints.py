# SPDX-License-Identifier: MIT
"""Read-only common physical sampling, with sealed native states as the input."""
import json,os,shutil
from pathlib import Path
from qualify_rans_preservation import load,replay_request
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import write_json,exact_vector_hash,compare_samples
from rans_duration_contract import require
from rans_preservation_oracle import Reference,physical_grid
from diagnose_rans_balanced_spatial import run_probe


def sample(root,labels,time,output,runtime):
    folders=[root/'campaigns'/s for s in labels];requests=[load(p/'request-manifest.json') for p in folders]
    ref=Reference(*[requests[0]['request'][k] for k in ('model','policy','wall')])
    knots=sorted(set(ref.knots()+[z for q in requests for z in q['request']['knots']]))
    grid,weights=physical_grid(knots,10,8)
    native=[load(p/('T-'+format(time,'.4f'))/'case.evolution.json') for p in folders]
    intent=dict(schema='SealedCommonEndpointSamplingV1',requests={label:sha(f/'request-manifest.json') for label,f in zip(labels,folders,strict=True)},
        time=time,source_reports={label:sha(f/('T-'+format(time,'.4f'))/'case.evolution.json') for label,f in zip(labels,folders,strict=True)},
        grid=grid,weights=weights,sides=[0]*len(grid),runtime=runtime_identity(runtime/'RANSPreservationProbe'),
        scope='unchanged physical coordinates; union of native partitions and analytical active profile breaks; actual initial fields, one common pressure gauge; no trajectory or stored sample overwrite')
    seal=RequestManifest(output,intent,output.exists())
    if seal.completed(output):return load(output/'results.json')['sampling']
    os.environ['LD_LIBRARY_PATH']=str(runtime/'lib');sampling=[];files=[]
    for label,folder,req,n in zip(labels,folders,requests,native,strict=True):
        own=RequestManifest(folder,req['request'],True);endpoint=folder/('T-'+format(time,'.4f'))
        require(own.completed(endpoint),'common_sampling_unsealed_endpoint')
        first=folder/('T-'+format(req['request']['horizons'][0],'.4f'));initial=load(first/'case.evolution.json')['accepted_states'][0]
        seed=load(folder/'initial-packet.json');values={}
        for kind,state,t in [('initial',initial,0.),('final',n['state'],n['time'])]:
            dest=output/(label+'-'+kind);dest.mkdir();shutil.copy2(endpoint/'case.xml',dest/'case.xml')
            request=replay_request(seed,grid,state=state,time=t)
            request.update(fixed_boundary=n['initialization_receipt']['live_fixed_boundary'],grid_sides=[0]*len(grid))
            packet,_=run_probe(runtime,dest,request);o=dict(packet['observation'],grid_weights=weights)
            require(exact_vector_hash(o['state'])==exact_vector_hash(state),'common_sampling_changed_state')
            require(o['grid']==grid and o['grid_sides']==[0]*len(grid),'common_sampling_changed_coordinates_or_side')
            values[kind]=o;files += [dest/'request.json',dest/'packet.json.gz',dest/'case.xml',dest/'probe.log']
        provenance=dict(native_report_sha256=sha(endpoint/'case.evolution.json'),
                        initial_report_sha256=sha(first/'case.evolution.json'),sampling_request_sha256=seal.payload['sha256'])
        sampling.append(dict(final=values['final'],initial=values['initial'],native_report_sha256=provenance['native_report_sha256'],provenance=provenance))
    for row in sampling:compare_samples(sampling[0]['final'],row['final'],weights)
    write_json(output/'results.json',dict(status='pass',sampling=sampling));seal.seal(output,['results.json',*files]);return sampling
