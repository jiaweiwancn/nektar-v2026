#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh corrected-operator short matrix, with independent configuration gates."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest, sha, runtime_identity
from rans_duration_contract import MEMBERS, prefix, aggregate, scalar_coverage, require
from rans_preservation_oracle import Reference
from make_rans_balanced_meshes import specifications
from audit_rans_common_scalar_partitions import partitions
from qualify_rans_duration_inner_k import reduction
from analyze_rans_duration_constitutive import analyze as constitutive
from rans_worker_budget import spawn

OPERATOR='GeneratedVolumeOwnerV1'
TIMES=(.0016,.016)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--geometry-receipt',type=Path,required=True)
    p.add_argument('--concurrency',type=int,choices=(1,2),default=2)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();r=a.root.resolve();source=Path(__file__).parent
    geometry=load(a.geometry_receipt);require(geometry['status']=='pass','native_geometry_not_qualified')
    require(geometry['contract']['numerical_operator']==OPERATOR,'foreign_geometry_operator')
    runtime=r/'runtime';os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    native_manifest=load(runtime/'source-build.json')
    require(native_manifest['numerical_operator']==OPERATOR,'legacy_runtime_not_qualifying')
    runtime_ids={n:runtime_identity(runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')}
    config_list=[dict(model=m,policy='Natural',wall=w) for m in ('E94E95','A95','A93') for w in ('Accepted','EP')]
    critical=['run_rans_neutral_pilot.py','run_a93_a95_e94e95_matrix.py','rans_preservation_oracle.py',
              'rans_qualification.py','rans_request_manifest.py','qualify_rans_preservation.py',
              'compare_rans_balanced_evolution.py','diagnose_rans_balanced_spatial.py','rans_spatial_resolution.py',
              'make_rans_balanced_meshes.py','classify_rans_balanced_spatial.py','audit_rans_balanced_evolved.py',
              'audit_rans_common_scalar_partitions.py','audit_rans_pilot_samples.py',
              'check_rans_balanced_state_quadrature.py','rans_duration_contract.py',
              'analyze_rans_duration_constitutive.py','validate_a93_a95_e94e95.py',
              'a93_a95_e94e95_note_wall_oracle.py','run_cls_epawf_picard.py','run_cls_epawf_pimple.py',
              'continue_cls_epawf_pimple.py','cls_epawf_m4_cases.py']
    critical.append('qualify_rans_duration_inner_k.py')
    source_ids={n:sha(source/n) for n in critical}
    intent=dict(schema='NativeOwnershipShortMatrixV1',numerical_operator=OPERATOR,
        geometry_receipt_sha256=sha(a.geometry_receipt),configurations=config_list,
        spaces=['Balanced5V1','Balanced11V1'],timesteps=[.0004,.0002],required_times=TIMES,
        runtime=runtime_ids,source_sha256=source_ids,controller_sha256=sha(Path(__file__)),
        gates='unchanged physical field/traction/domain, .85 dominant native/common and inner-only k; original temporal floors',
        history='fresh corrected initializations, native same-method continuation only')
    out=r/'ownership-short-qualification';seal=RequestManifest(out,intent,a.resume)
    if seal.completed(out):return 0 if load(out/'results.json')['status']=='pass' else 1
    lock=threading.Lock();events=[];results={}
    env=dict(os.environ,**THREADS,NEKTAR_RANS_PROBE_TIMEOUT='14400')

    def call(name,arguments,artifacts):
        # Snapshot exact scientific scripts before every subprocess launch.
        require(all(sha(source/n)==h for n,h in source_ids.items()),'short_scientific_source_changed')
        cmd=[sys.executable,'-B',*map(str,arguments)]
        job=out/'jobs'/name
        request=dict(schema='OwnedShortJobV1',parent_request_sha256=seal.payload['sha256'],command=cmd,
                     numerical_operator=OPERATOR,runtime=runtime_ids)
        manifest=RequestManifest(job,request,a.resume and job.exists())
        if manifest.completed(job):return load(job/'results.json')['exit_code']
        with (job/'execution.log').open('x') as f:
            proc=None
            while proc is None:
                proc=spawn(r,cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
                if proc is None:time.sleep(2)
            code=proc.wait()
        value=dict(exit_code=code,status='completed' if code==0 else 'failed_or_unresolved',command=cmd)
        write_json(job/'results.json',value)
        files=['results.json','execution.log']
        for path in artifacts:
            if path.is_dir():files.extend(str(p) for p in path.rglob('*') if p.is_file())
            elif path.exists():files.append(str(path))
        manifest.seal(job,files)
        with lock:
            events.append(dict(name=name,**value));write_json(out/'execution.json',dict(events=events))
        return code

    def run_batch(jobs):
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(call,*job) for job in jobs]
            return [f.result() for f in futures]

    def configuration(cfg):
        model,policy,wall=(cfg[k] for k in ('model','policy','wall'))
        tag=prefix(model,policy,wall);ref=Reference(model,policy,wall)
        screen='initial-'+wall if model=='E94E95' else model+'-Natural-'+wall
        mesh=r/'meshes'/tag
        if not mesh.exists():
            mesh.mkdir(parents=True);spaces=specifications(ref)
            for n,spec in spaces.items():write_json(mesh/(n+'.json'),spec)
            spaces['BaseP5']=dict(degree=5,knots=ref.knots())
            write_json(mesh/'screen-spaces.json',spaces)
        initial=r/'diagnosis'/screen
        code=call(tag+'-initial',[source/'diagnose_rans_balanced_spatial.py','--root',r,
            '--analysis-runtime',runtime,'--model',model,'--wall',wall,'--label',screen,
            '--specifications',mesh/'screen-spaces.json'],[initial])
        initial_gate=r/('initial-gate-'+(wall if model=='E94E95' else model+'-Natural-'+wall)+'.json')
        gate_code=call(tag+'-initial-gate',[source/'classify_rans_balanced_spatial.py','--root',r,
            '--model',model,'--wall',wall,'--screen',screen,'--output',initial_gate],[initial_gate]) if code==0 else 1
        row=dict(configuration=cfg,numerical_operator=OPERATOR,initial_status='pass' if gate_code==0 else 'fail',milestones=[])
        if gate_code:
            row['status']='failed_initial_spatial_gate'
            return row
        common=initial/'common'
        with gzip.open(common/'packet.json.gz','rt') as f:cp=json.load(f)
        require(cp['observation']['numerical_operator']==OPERATOR,'old_common_test_packet')
        for T in TIMES:
            jobs=[]
            for member in MEMBERS:
                label=tag+'-'+member;folder=r/'campaigns'/label
                mesh_name='Balanced5V1' if member.startswith('B5') else 'Balanced11V1'
                cmd=[source/'run_rans_neutral_pilot.py','--root',r,'--label',label,'--model',model,
                     '--wall',wall,'--mesh-file',mesh/(mesh_name+'.json'),'--dt',
                     str(.0002 if member.endswith('dt2') else .0004),'--milestones','.0016','.016',
                     '--max-time',str(T),'--timeout','14400']
                if folder.exists():cmd+=['--resume']
                jobs.append((tag+'-evolve-'+member+'-'+str(T),cmd,
                             [folder/'request-manifest.json',folder/('T-'+format(T,'.4f'))]))
            executions=run_batch(jobs)
            comparison=r/'comparisons'/(model+'-'+policy+'-'+wall+'-'+str(T)+'.json')
            compare_code=call(tag+'-compare-'+str(T),[source/'compare_rans_balanced_evolution.py',
                '--root',r,'--labels',*[tag+'-'+m for m in MEMBERS],'--time',str(T),'--output',comparison],
                [comparison]) if all(v==0 for v in executions) else 1
            if not all(v==0 for v in executions):
                row['milestones'].append(dict(time=T,status='failed_or_unresolved',comparison=str(comparison),
                                              process_exit_codes=executions));break
            audit_jobs=[];spatial_paths={};scalar_paths={}
            for member in MEMBERS:
                label=tag+'-'+member
                if member in MEMBERS[:2]:
                    name='owned-short-'+tag+'-'+str(T)+'-'+member;dest=r/'evolved-audits'/name
                    spatial_paths[member]=dest/'results.json'
                    audit_jobs.append((name,[source/'audit_rans_balanced_evolved.py','--root',r,
                        '--analysis-runtime',runtime,'--campaign-root',r,'--label',label,'--output-label',name,
                        '--times',str(T),'--common-directory',common,'--constitutive-samples'],[dest]))
                observer='owned-short-scalars-'+str(T);dest=r/'campaigns'/label/observer
                scalar_paths[member]=dest/'results.json'
                audit_jobs.append((tag+'-scalars-'+member+'-'+str(T),[source/'audit_rans_pilot_samples.py',
                    '--root',r,'--label',label,'--through-time',str(T),'--output-label',observer],[dest]))
            qlabel='owned-short-quadrature-'+tag+'-'+str(T);qpath=r/qlabel/'results.json'
            audit_jobs.append((qlabel,[source/'check_rans_balanced_state_quadrature.py','--root',r,
                '--analysis-runtime',runtime,'--model',model,'--walls',wall,'--time',str(T),
                '--output-label',qlabel],[r/qlabel]))
            audits_execution=run_batch(audit_jobs)
            audits={m:load(p) for m,p in spatial_paths.items() if p.exists()}
            companions={};native_reports={};stress={}
            for member in MEMBERS:
                folder=r/'campaigns'/(tag+'-'+member);request=load(folder/'request-manifest.json')['request']
                require(request['numerical_operator']==OPERATOR,'old_trajectory_request')
                reports=[load(folder/('T-'+format(t,'.4f'))/'case.evolution.json') for t in TIMES if t<=T]
                require(all(v['numerical_operator']==OPERATOR for v in reports),'old_native_report')
                native_reports[member]=[sha(folder/('T-'+format(t,'.4f'))/'case.evolution.json') for t in TIMES if t<=T]
                obs=load(scalar_paths[member]);coverage=scalar_coverage(obs,reports,cfg,T)
                companions[member]=dict(status=obs['status'],configuration=cfg,scalar_coverage=coverage,
                                        observer_source=str(scalar_paths[member]),observer_sha256=sha(scalar_paths[member]))
                if member in MEMBERS[:2] and member in audits:
                    packet=load_gzip(spatial_paths[member].parent/('T-'+format(T,'.4f'))/'packet.json.gz')
                    require(packet['numerical_operator']==OPERATOR,'old_spatial_packet')
                    stress[member]=constitutive(packet,dict(cfg,descriptor=request['descriptor']),request['evaluation']['weights'])
            groups={m:partitions(cp,audits[m]['checks'][0],ref) for m in MEMBERS[:2]}
            inner=reduction(*[groups[m]['inner']['mass_dual'] for m in MEMBERS[:2]])
            receipt=dict(configuration=cfg,numerical_operator=OPERATOR,time=T,companions=companions,
                spatial_audits={m:dict(audits[m],configuration=cfg) for m in audits},
                constitutive=stress,quadrature=load(qpath),native_report_sha256=native_reports,
                spatial_sources={m:dict(path=str(p),sha256=sha(p)) for m,p in spatial_paths.items()},
                comparison=dict(load(comparison),source=str(comparison),sha256=sha(comparison)),
                fixed_common_inner_k=dict(reduction=inner,groups=groups),execution_exit_codes=audits_execution)
            gate=aggregate(dict(cfg,audit_times=[T]),receipt)
            receipt.update(aggregate=gate,status='pass' if compare_code==0 and gate['status']=='pass' and inner['status'] in ('pass','unresolved_floor') and
                           all(v==0 for v in audits_execution) and all(v['status']=='pass' for v in stress.values()) else 'fail')
            target=r/'short-audits'/(tag+'-'+str(T));target.mkdir(parents=True,exist_ok=False)
            write_json(target/'results.json',receipt)
            row['milestones'].append(dict(time=T,status=receipt['status'],receipt=str(target/'results.json'),sha256=sha(target/'results.json')))
            if receipt['status']!='pass':break
        row['status']='pass' if len(row['milestones'])==2 and all(v['status']=='pass' for v in row['milestones']) else 'failed_or_unresolved'
        return row

    def completed(cfg):
        tag=prefix(**cfg)
        try:row=configuration(cfg)
        except Exception as error:row=dict(configuration=cfg,numerical_operator=OPERATOR,status='fail',failure=repr(error))
        with lock:
            results[tag]=row;write_json(out/'progress.json',dict(status='running',configurations=results))
        print(tag,row.get('status','fail'),flush=True)
        return row
    with ThreadPoolExecutor(max_workers=a.concurrency) as pool:
        rows=list(pool.map(completed,config_list))
    result=dict(status='pass' if all(q.get('status')=='pass' for q in rows) else 'failed_unresolved_or_incomplete',
                numerical_operator=OPERATOR,configurations=results,execution=events,
                scope='six Natural configurations; corrected runtime and fresh histories only; no old evolution reuse')
    write_json(out/'results.json',result)
    files=['results.json']+[str(p) for p in (r/'short-audits').glob('*/results.json')]
    seal.seal(out,files)
    return 0 if result['status']=='pass' else 1


def load_gzip(path):
    with gzip.open(path,'rt') as f:return json.load(f)


if __name__=='__main__':raise SystemExit(main())
