#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Admit one model/policy/wall from sealed, identity-checked short evidence."""
import argparse,gzip,json,os
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import THREADS,write_json
from rans_duration_contract import identity,require,required_times,residual_reduction,scalar_coverage,prefix,SHORT_TIMES


def runtime_bytes(value):
    return dict(sha256=value['sha256'],libraries={Path(p).name:h for p,h in value['libraries'].items()})


def validate_audit(path,campaign,config,times,common):
    data=load(path);identity(load(campaign/'request-manifest.json')['request'],config)
    require(data['status']=='pass','spatial_audit_failed')
    required_times(data['checks'],times)
    intent=load(path.parent/'immutable-request.json')
    require(intent['common_packet_sha256']==sha(common/'packet.json.gz'),'wrong_common_test_packet')
    req=load(campaign/'request-manifest.json');manifest=RequestManifest(campaign,req['request'],True)
    for row in data['checks']:
        endpoint=campaign/('T-'+format(row['time'],'.4f'))
        require(manifest.completed(endpoint),'unsealed_spatial_endpoint')
        require(row['trajectory_request_sha256']==req['sha256'],'spatial_request_mismatch')
        require(row['native_report_sha256']==sha(endpoint/'case.evolution.json'),'spatial_native_report_mismatch')
        require(row['state_sha256']==load(endpoint/'metrics.json')['final_state_sha256'],'spatial_state_mismatch')
        if 'configuration' in row:identity(row['configuration'],config)
    return dict(data,configuration=config,source=str(path),source_sha256=sha(path),common_packet_sha256=intent['common_packet_sha256'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--previous-root',type=Path,required=True);p.add_argument('--model',choices=('A93','A95','E94E95'),required=True);p.add_argument('--policy',choices=('Natural','Prescribed'),default='Natural');p.add_argument('--wall',choices=('EP','Accepted'),required=True);a=p.parse_args();r=a.root.resolve();old=a.previous_root.resolve();os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
    tag=prefix(a.model,a.policy,a.wall);config=dict(model=a.model,policy=a.policy,wall=a.wall)
    first=old/'campaigns'/(tag+'-B5');req=load(first/'request-manifest.json')['request'];identity(req,config);config['descriptor']=req['descriptor']
    common=old/'diagnosis'/((a.model+'-'+a.policy+'-'+a.wall) if a.model!='E94E95' else 'initial-'+a.wall)/'common'
    with gzip.open(common/'packet.json.gz','rt') as f:cp=json.load(f)
    require(cp['observation']['resolved_descriptor']==config['descriptor'],'common_descriptor_mismatch')
    initial=old/('initial-gate-'+(a.model+'-'+a.policy+'-' if a.model!='E94E95' else '')+a.wall+'.json')
    require(load(initial)['status']=='pass','initial_spatial_not_qualified')
    runtimes={n:runtime_identity(r/'runtime'/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')}
    bridges={}
    for name,value in runtimes.items():
        require(runtime_bytes(value)==runtime_bytes(req['runtime'][name]),'reused_runtime_bytes_differ:'+name)
        bridges[name]=dict(current=value,previous=req['runtime'][name],bytes_and_loaded_dependencies_equal=True)
    # Preserve content provenance even if a future source checkout has advanced.
    files=load(old/'mathematical-source-reuse.json')['governing_files'];source_proof=[]
    for row in files:
        current=r/'source'/row['file'];require(sha(current)==row['sha256'],'governing_source_changed:'+row['file']);source_proof.append(row)
    meshes={};audits={};companions={};comparisons={};proof=[]
    for member in ('B5','B11','B5-dt2','B11-dt2'):
        folder=old/'campaigns'/(tag+'-'+member);m=load(folder/'request-manifest.json');q=m['request'];identity(q,config);manifest=RequestManifest(folder,q,True)
        require(runtime_bytes(q['runtime']['IncNavierStokesSolver'])==runtime_bytes(runtimes['IncNavierStokesSolver']),'companion_runtime_mismatch')
        reports=[]
        for t in SHORT_TIMES:
            endpoint=folder/('T-'+format(t,'.4f'));require(manifest.completed(endpoint),'missing_short_companion:'+member)
            metric=load(endpoint/'metrics.json');require(metric['status']=='pass' and abs(metric['time']-t)<1e-13,'short_endpoint_failed')
            reports.append(load(endpoint/'case.evolution.json'))
        obs_path=folder/'scalar-observations/results.json';obs=load(obs_path)
        with gzip.open(obs_path.parent/'native-basis.json.gz','rb') as stream:basis_bytes=stream.read()
        require(__import__('hashlib').sha256(basis_bytes).hexdigest()==obs['source_packet_sha256'],'scalar_basis_identity_mismatch')
        require(json.loads(basis_bytes)['observation']['resolved_descriptor']==config['descriptor'],'scalar_configuration_mismatch')
        coverage=scalar_coverage(obs,reports,config,max(SHORT_TIMES))
        companions[member]=dict(status='pass',configuration=config,scalar_coverage=coverage,observer_source=str(obs_path),observer_sha256=sha(obs_path),request_sha256=m['sha256'])
        if member in ('B5','B11'):
            suffix='Balanced5V1' if member=='B5' else 'Balanced11V1';mesh=r/'meshes'/tag/(suffix+'.json');spec=load(mesh)
            require(spec['knots']==q['knots'] and spec['degree']==q['degree'],'mesh_request_mismatch')
            meshes[member]=dict(path=str(mesh),sha256=sha(mesh),specification=spec)
            audit_path=old/'evolved-audits'/('balanced-'+tag+'-'+member)/'results.json'
            audits[member]=validate_audit(audit_path,folder,config,SHORT_TIMES,common)
    reductions=[]
    for t in SHORT_TIMES:
        rows=[next(c for c in audits[m]['checks'] if abs(c['time']-t)<1e-13) for m in ('B5','B11')]
        gate=residual_reduction(*rows,a.model);require(gate['status']=='pass','evolved_reduction_failed:'+str(t));reductions.append(dict(time=t,**gate))
        path=old/'comparisons'/(a.model+'-'+a.policy+'-'+a.wall+'-'+str(t)+'.json');c=load(path);require(c['status']=='pass','short_comparison_failed');require(abs(c['time']-t)<1e-13,'short_comparison_time_mismatch')
        require(len(c['provenance'])==4,'missing_short_companion')
        require([v['label'] for v in c['provenance']]==[tag+'-'+m for m in ('B5','B11','B5-dt2','B11-dt2')],'short_configuration_mismatch')
        comparisons[str(t)]=dict(c,source=str(path),source_sha256=sha(path))
    quadpath=r/('frozen-quadrature-'+tag)/'results.json';quad=load(quadpath)
    require(quad['status']=='pass' and len(quad['checks'])==2,'model_quadrature_incomplete')
    for check in quad['checks']:identity(check['configuration'],config);require(abs(check['time']-.016)<1e-13,'quadrature_time_mismatch')
    output=dict(schema='ConfigurationDurationAdmissionV2',status='pass',configuration=config,
        required_short_times=list(SHORT_TIMES),milestones=[.016,.1,.5],timestep_trials=[[.002,.001],[.001,.0005]],
        meshes=meshes,common_directory=str(common),common_packet_sha256=sha(common/'packet.json.gz'),
        runtime=runtimes,diagnostic_runtime=runtime_identity(r/'runtime-resolution-v2/RANSPreservationProbe'),
        runtime_reuse_bridges=bridges,governing_source_proof=source_proof,initial_gate=dict(source=str(initial),sha256=sha(initial),body=load(initial)),
        short_spatial_audits=audits,short_companions=companions,short_comparisons=comparisons,evolved_reductions=reductions,
        frozen_quadrature=dict(quad,source=str(quadpath),sha256=sha(quadpath)),previous_root=str(old),
        neutral_gates=dict(field_RMS=1e-3,vertical_RMS=1e-4,wall_traction=1e-3,interface_traction=1e-3,temporal_floor=1e-7,temporal_fraction=.1,residual_floor=1e-8,residual_ratio=.85,residual_uncertainty=1e-10),
        history='independently initialized at each dt; unchanged-method restart within trajectory')
    dest=r/'admissions'/(tag+'.json');dest.parent.mkdir(exist_ok=True);require(not dest.exists(),'admission_already_exists');write_json(dest,output);print(str(dest));return 0


if __name__=='__main__':raise SystemExit(main())
