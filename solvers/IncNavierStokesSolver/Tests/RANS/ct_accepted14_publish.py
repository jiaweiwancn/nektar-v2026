#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish sealed CT14 evidence; add this command's closed receipt afterward."""
import argparse,json,os,resource,subprocess,tarfile,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted14_compact import compact
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def read(path):return json.loads(path.read_text())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','capacity','component','static','pair','collection','runtime','self-control','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--admission',type=Path)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    bodies={key:sealed(getattr(a,key)) for key in ('capacity','component','static','pair','collection')}
    pair_request=read(a.pair/'request-manifest.json')['request']
    if pair_request['dt']!=[.001,.0005]:
        require(a.admission is not None and sha(a.admission/'results.json')==pair_request['timestep_admission'],'missing publication timestep admission')
        bodies['admission']=sealed(a.admission);require(bodies['admission']['status']=='pass','unqualified publication timestep admission')
    require(bodies['collection']['closure'] in ('PASS','UNRESOLVED'),'physical collection incomplete')
    source=a.root/'source';git=['git','-c','index.threads=1','-c','core.preloadIndex=false','-c','core.filemode=false','-C',str(source)]
    revision=subprocess.check_output(git+['rev-parse','HEAD'],text=True).strip()
    require(not subprocess.check_output(git+['status','--porcelain'],text=True).strip(),'dirty publication source')
    manifest=source/'docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json'
    require(sha(manifest)=='004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296','changed physical manifest')
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    runtime={key:runtime_identity(a.runtime/name) for key,name in (('governing','IncNavierStokesSolver'),('probe','RANSPreservationProbe'))}
    require(runtime['governing']==bodies['collection']['runtime'] and runtime['probe']==bodies['collection']['probe'],'foreign publication runtime')
    request=RequestManifest(a.output,dict(schema='CT14PublicationV1',source=sha(Path(__file__)),source_revision=revision,inputs={key:sha(getattr(a,key)/'results.json') for key in bodies},runtime=runtime))
    for key in ('pair','collection'):
        original=getattr(a,key)/'results.json';body,omitted=compact(bodies[key])
        body['compact_publication']=dict(schema='CT14CompactPublicationV1',complete_external_body=str(original),complete_body_sha256=sha(original),omitted_payloads=omitted,rule='aggregate primitives/classifications unchanged; large vectors,per-point rows and native identity strings retained externally; not a classifier/replay input')
        write_json(a.output/(key+'.json'),body)
    families={};roots=set();cpu=wall=0.;aggregate=0
    for folder in sorted(a.root.glob('resources-*')):
        if folder==a.self_control:continue
        f=folder/'execution.json';require(f.is_file(),'unclosed owned family:'+str(folder));e=read(f);r=e['resource']
        require(r['project_status']==r['aggregate_status']=='pass','worker violation:'+str(folder))
        families[folder.name]=dict(path=str(f),sha256=sha(f),execution=e)
        roots.update(e[k] for k in ('control_parent_pid','guard_pid','owned_family_pid'))
        cpu+=e['process_cpu_seconds'];wall+=max(0.,e['wall_seconds']-e['allocation_wait_seconds']);aggregate=max(aggregate,r['maximum_aggregate'])
    peak=0;witness=None
    for name in families:
        high=a.root/name/'worker-guard-highwater.jsonl'
        if not high.exists():continue
        for line in high.read_text().splitlines():
            row=json.loads(line);owned=set(roots)
            while True:
                add={v['pid'] for v in row['processes'] if v['parent'] in owned}-owned
                if not add:break
                owned.update(add)
            count=sum(v['threads'] for v in row['processes'] if v['pid'] in owned)
            if count>peak:peak=count;witness=dict(time_unix=row['time_unix'],threads=count,processes=[v for v in row['processes'] if v['pid'] in owned])
    require(peak<=16 and aggregate<=72,'combined worker allocation')
    output=sum(f.lstat().st_size for f in a.root.rglob('*') if f.is_file() or f.is_symlink());whole=time.time()-1794751358
    require(cpu+120<=24000 and max(wall,whole)+240<=28800 and output+256*1024**2<=16*1024**3,'bounded publication allocation')
    write_json(a.output/'prepublication-accounting.json',dict(schema='CT14PrepublicationAccountingV1',status='pass',families=families,all_previous_families_closed=True,closed_CPU_seconds=cpu,active_wall_upper_seconds=max(wall,whole),whole_elapsed_package_seconds=whole,retained_new_bytes=output,maximum_combined_project=peak,maximum_aggregate=aggregate,combined_project_witness=witness,CPU_reserve_seconds=120,wall_reserve_seconds=240,output_reserve_bytes=256*1024**2,scope='add this publication command closed receipt afterward; source/transport/inspection reserve remains explicit; no future exit attested'))
    index={}
    def record(path,body=False):
        require(path.is_file(),'missing external artifact:'+str(path));r=dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
        if body:r['body']=read(path)
        index[str(path.relative_to(a.root))]=r
    for folder in sorted(a.root.glob('resources-*')):
        watch=folder/'package-budget-watch.json'
        if watch.exists():record(watch,True)
        for filename in ('live-memory-observations.jsonl','live-memory-summary.json'):
            memory=folder/filename
            if memory.exists():record(memory,filename.endswith('summary.json'))
    guard_controls=a.root/'package-guard-controls-v1/results.json'
    if guard_controls.exists():record(guard_controls,True)
    for key in bodies:
        folder=getattr(a,key)
        for name in ('request-manifest.json','completed-artifacts.json'):record(folder/name,True)
        record(folder/'results.json')
    for name,row in bodies['pair']['members'].items():
        folder=Path(row['directory'])
        for filename in ('case.xml','case.evolution.json','sparse-receipts.jsonl'):record(folder/filename)
        native=read(folder/'case.evolution.json')
        for stage in native['pimple_attempts']:
            for key in ('flow_snapshot','k_snapshot'):
                if stage.get(key):record(Path(stage[key]))
        zero=a.pair/('zero-'+name)
        for filename in ('case.xml','case.evolution.json','sparse-receipts.jsonl'):record(zero/filename)
    mount=subprocess.check_output(['findmnt','-T',str(a.root),'-n','-o','TARGET,SOURCE,FSTYPE,OPTIONS'],text=True).strip()
    require('/media/jay/新加卷1 ' in mount,'unexpected authoritative mount')
    write_json(a.output/'external-artifact-index.json',dict(schema='CT14ExternalArtifactIndexV1',files=index,scope='complete input/runtime/state/history/matrix bindings remain in the unchanged request/result/completion bodies; large payloads remain external'))
    write_json(a.output/'provenance.json',dict(schema='CT14FinalPublicationProvenanceV1',source_revision=revision,source_clean=True,archive=str(a.root),mount=mount,runtime=runtime,physical_manifest_sha256=sha(manifest),component_identity='CT14-VectorCOLAMDV1',native_tested_source='6477d9d988e969686a418c8e333fadfb85bb3f9d',controller_scopes='component at6477d9d9,static at087fa720,first pair ateca8f3ab,changed-alpha admission atae18441f,retry atb4867df7; individual source hashes in immutable requests',closure=bodies['collection']['closure'],threads=THREADS,history='historical results,failed attempts and shared installations preserved'))
    files=[f.name for f in a.output.iterdir() if f.is_file() and f.name!='request-manifest.json'];request.seal(a.output,files)
    archive=a.output/'ct14-compact-v1.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=6) as tar:
        for f in sorted(a.output.iterdir()):
            if f!=archive:tar.add(f,arcname=f.name)
    write_json(a.output/'transfer-sha256.json',dict(file=archive.name,bytes=archive.stat().st_size,sha256=sha(archive)))
    print(json.dumps(dict(status='pass',closure=bodies['collection']['closure'],closed_CPU=cpu,whole_elapsed=whole,retained_bytes=output,transfer=read(a.output/'transfer-sha256.json'))),flush=True)

if __name__=='__main__':main()
