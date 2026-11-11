#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Admit fresh GeneratedVolumeOwnerV1 cohorts, never bridge old evolution."""
import argparse
import gzip
import json
import os
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest, sha, runtime_identity
from rans_duration_contract import MEMBERS, SHORT_TIMES, identity, require, prefix, required_times, aggregate
from admit_rans_configuration import validate_audit
from rans_qualification import write_json

OPERATOR = 'GeneratedVolumeOwnerV1'


def pinned_runtime(path):
    previous=os.environ.get('LD_LIBRARY_PATH')
    try:
        os.environ['LD_LIBRARY_PATH']=str(path.parent/'lib')
        return runtime_identity(path)
    finally:
        if previous is None:os.environ.pop('LD_LIBRARY_PATH',None)
        else:os.environ['LD_LIBRARY_PATH']=previous


def admit(root, model, wall, output=None):
    r=root.resolve(); cfg=dict(model=model,policy='Natural',wall=wall);tag=prefix(**cfg)
    geometry=r/'geometric-regressions-v1/results.json';g=load(geometry)
    require(g['status']=='pass' and g['contract']['numerical_operator']==OPERATOR,'native_geometry_not_qualified')
    first=r/'campaigns'/(tag+'-B5');q=load(first/'request-manifest.json')['request']
    identity(q,cfg);cfg['descriptor']=q['descriptor']
    runtimes={n:pinned_runtime(r/'runtime'/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')}
    build=r/'runtime/source-build.json';require(load(build)['numerical_operator']==OPERATOR,'legacy_compiled_operator')
    common=r/'diagnosis'/('initial-'+wall if model=='E94E95' else model+'-Natural-'+wall)/'common'
    with gzip.open(common/'packet.json.gz','rt') as stream:cp=json.load(stream)
    require(cp['observation']['resolved_descriptor']==cfg['descriptor'],'common_descriptor_mismatch')
    require(cp['observation']['numerical_operator']==OPERATOR,'old_common_test_operator')
    initial=r/('initial-gate-'+(model+'-Natural-' if model!='E94E95' else '')+wall+'.json')
    require(load(initial)['status']=='pass','initial_spatial_not_qualified')
    receipts={};audits={m:dict(status='pass',configuration=cfg,checks=[]) for m in MEMBERS[:2]}
    companions={};comparisons={};quadratures={}
    for t in SHORT_TIMES:
        p=r/'short-audits'/(tag+'-'+str(t))/'results.json';a=load(p)
        simple={k:cfg[k] for k in ('model','policy','wall')}
        identity(a['configuration'],simple);require(a['numerical_operator']==OPERATOR,'old_short_operator')
        require(a['status']=='pass' and aggregate(dict(simple,audit_times=[t]),a)['status']=='pass','short_required_gate_failed:'+str(t))
        require(a['fixed_common_inner_k']['reduction']['status'] in ('pass','unresolved_floor'),'common_inner_k_failed')
        require(a['quadrature']['status']=='pass' and len(a['quadrature']['checks'])==2,'incomplete_frozen_quadrature')
        require(len(a['constitutive'])==2 and all(s['status']=='pass' for s in a['constitutive'].values()),'constitutive_gate_failed')
        for m in MEMBERS[:2]:
            source=Path(a['spatial_sources'][m]['path']);require(sha(source)==a['spatial_sources'][m]['sha256'],'short_audit_changed')
            validated=validate_audit(source,r/'campaigns'/(tag+'-'+m),cfg,[t],common,
                                     expected_runtime=pinned_runtime(r/'runtime/RANSPreservationProbe'))
            audits[m]['checks'].extend(validated['checks'])
        receipts[str(t)]=dict(path=str(p),sha256=sha(p),body=a)
        comparisons[str(t)]=a['comparison'];quadratures[str(t)]=a['quadrature']
        if t==max(SHORT_TIMES):companions=a['companions']
    for a in audits.values():required_times(a['checks'],SHORT_TIMES)
    meshes={};bindings={}
    for m in MEMBERS:
        folder=r/'campaigns'/(tag+'-'+m);rm=load(folder/'request-manifest.json');request=rm['request']
        identity(request,cfg);require(request['numerical_operator']==OPERATOR,'old_trajectory_operator')
        require(request['runtime']==runtimes,'trajectory_runtime_changed')
        seal=RequestManifest(folder,request,True);reports={}
        for t in SHORT_TIMES:
            endpoint=folder/('T-'+format(t,'.4f'));require(seal.completed(endpoint),'unsealed_short_endpoint')
            report=load(endpoint/'case.evolution.json');require(report['numerical_operator']==OPERATOR,'old_native_report')
            reports[str(t)]={n:sha(endpoint/n) for n in ('case.evolution.json','metrics.json')}
        bindings[m]=dict(request_sha256=rm['sha256'],native_reports=reports,
                         consumed_inputs=load(folder/'T-0.0016/metrics.json').get('initial_identity',{}))
        if m in MEMBERS[:2]:
            path=r/'meshes'/tag/('Balanced5V1.json' if m=='B5' else 'Balanced11V1.json');spec=load(path)
            require(spec['knots']==request['knots'] and spec['degree']==request['degree'],'admitted_mesh_mismatch')
            meshes[m]=dict(path=str(path),sha256=sha(path),specification=spec)
    value=dict(schema='NativeOwnershipDurationAdmissionV1',status='pass',configuration=cfg,
        numerical_operator=OPERATOR,corrective_contract=dict(schema='NativeOwnershipRequiredScopeV1',
            response_models=['A95'],response_max_time=.016,A93_response='outside_scope'),
        required_short_times=list(SHORT_TIMES),milestones=[.016,.1,.5],timestep_trials=[[.002,.001],[.001,.0005]],
        meshes=meshes,common_directory=str(common),common_packet_sha256=sha(common/'packet.json.gz'),runtime=runtimes,
        diagnostic_runtime=pinned_runtime(r/'runtime-resolution-v2/RANSPreservationProbe'),
        duration_diagnostic_runtime=pinned_runtime(r/'runtime-resolution-v3/RANSPreservationProbe'),
        compiled_source_build=dict(path=str(build),sha256=sha(build),body=load(build)),
        geometric_qualification=dict(path=str(geometry),sha256=sha(geometry)),
        initial_gate=dict(source=str(initial),sha256=sha(initial),body=load(initial)),
        short_spatial_audits=audits,short_companions=companions,short_comparisons=comparisons,
        short_complete_receipts=receipts,short_bindings=bindings,frozen_quadrature=quadratures,previous_root=str(r),
        neutral_gates=dict(field_RMS=1e-3,vertical_RMS=1e-4,wall_traction=1e-3,interface_traction=1e-3,
            temporal_floor=1e-7,temporal_fraction=.1,residual_floor=1e-8,residual_ratio=.85,residual_uncertainty=1e-10),
        history='fresh corrected initialization at every timestep; unchanged-method native continuation within trajectory',
        reuse_scope='old coefficients only in separately labelled read-only diagnostics; no old-operator evolution admission')
    dest=output or r/'admissions'/(tag+'.json');dest.parent.mkdir(exist_ok=True)
    require(not dest.exists(),'admission_already_exists');write_json(dest,value);return value


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--model',choices=('A93','A95','E94E95'),required=True)
    p.add_argument('--wall',choices=('Accepted','EP'),required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();admit(a.root,a.model,a.wall,a.output);print('pass');return 0


if __name__=='__main__':raise SystemExit(main())
