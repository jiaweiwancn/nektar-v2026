#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Close CT13 publication accounting without modifying any earlier seal.

Run under ct_accepted3_command.py. The final launcher's closed execution receipt
must be added afterward; this program cannot attest to its own future exit.
"""
import argparse, json, os, shutil, subprocess, tarfile
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_request_manifest import sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','ct9','ct11','ct12','self-control','output'):
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);ev=a.root/'evaluated-v1'
    completion=json.loads((ev/'completed-artifacts.json').read_text())
    require(all(sha(ev/name)==digest for name,digest in completion['files'].items()),'changed evaluated artifact')
    families={};roots=set();cpu=wall=0.;aggregate=0
    for d in sorted(a.root.glob('resources-*')):
        if d==a.self_control:continue
        e=json.loads((d/'execution.json').read_text());r=e['resource']
        require(r['project_status']==r['aggregate_status']=='pass','resource violation:'+d.name)
        high=d/'worker-guard-highwater.jsonl';events=d/'worker-guard-events.jsonl'
        families[d.name]=dict(execution=e,sha256=sha(d/'execution.json'),
            highwater=[json.loads(x) for x in high.read_text().splitlines()] if high.exists() else [],
            events=[json.loads(x) for x in events.read_text().splitlines()] if events.exists() else [])
        roots.update(e[k] for k in ('control_parent_pid','guard_pid','owned_family_pid'))
        cpu+=e['process_cpu_seconds'];wall+=max(0.,e['wall_seconds']-e['allocation_wait_seconds'])
        aggregate=max(aggregate,r['maximum_aggregate'])
    peak=0;witness=None
    for v in families.values():
        for row in v['highwater']:
            owned=set(roots)
            while True:
                add={x['pid'] for x in row['processes'] if x['parent'] in owned}-owned
                if not add:break
                owned.update(add)
            n=sum(x['threads'] for x in row['processes'] if x['pid'] in owned)
            if n>peak:
                peak=n;witness=dict(time_unix=row['time_unix'],threads=n,
                    processes=[x for x in row['processes'] if x['pid'] in owned])
    output=sum(f.stat().st_size for f in a.root.rglob('*') if f.is_file()
        and f.relative_to(a.root).parts[0] not in ('source','tmp'))
    require(peak<=16 and aggregate<=72,'combined worker limit')
    require(cpu+60<24000 and wall+120<28800 and output+64*1024**2<16*1024**3,'bounded publication reserve')
    write_json(a.output/'prepublication-accounting.json',dict(
        schema='CT13ClosedPublicationSupplementV1',status='pass',families=families,
        all_previous_families_closed=True,family_count=len(families),
        process_CPU_seconds=cpu,active_wall_upper_seconds=wall,
        output_bytes_before_packaging=output,maximum_combined_project=peak,
        maximum_aggregate=aggregate,combined_project_witness=witness,
        publication_CPU_reserve=60,publication_wall_reserve=120,
        publication_output_reserve_bytes=64*1024**2,
        scope='add the closed publication receipt separately; summed elapsed minus allocation wait is a conservative active-wall upper bound'))
    source=a.root/'source'
    rev=subprocess.check_output(['git','-c','core.preloadIndex=false','rev-parse','HEAD'],cwd=source,text=True).strip()
    require(not subprocess.check_output(['git','-c','core.preloadIndex=false','status','--porcelain'],cwd=source,text=True).strip(),'dirty authoritative source')
    manifest=source/'docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json'
    require(sha(manifest)=='004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296','changed physical manifest')
    runtime={}
    for key,path,digest in (
        ('governing',a.ct9/'build-v8/runtime/IncNavierStokesSolver','5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7'),
        ('legacy_probe',a.ct11/'build-v2/runtime/RANSPreservationProbe','20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8')):
        require(sha(path)==digest,'changed pinned runtime:'+key)
        runtime[key]=dict(path=str(path),sha256=digest)
    for key,path in (
        ('diagnostic_probe',a.root/'build-v6/runtime/RANSPreservationProbe'),
        ('ordering_probe',a.root/'build-v2/runtime/CT13SparseOrderingProbe')):
        runtime[key]=dict(path=str(path),sha256=sha(path))
    mount=subprocess.check_output(['findmnt','-T',str(a.root),'-n','-o','TARGET,SOURCE,FSTYPE,OPTIONS'],text=True).strip()
    require('/media/jay/新加卷1 ' in mount,'unexpected experiment mount')
    changed=subprocess.check_output(['git','-c','core.preloadIndex=false','diff','--name-only','5799b240917e6ca7dcf4ac3d44f2b520f5629955','HEAD'],cwd=source,text=True).splitlines()
    require(not any(name.startswith('docs/plans/cao-tamura-accepted-') and not name.startswith(('docs/plans/cao-tamura-accepted-13.','docs/plans/cao-tamura-accepted-13-evidence/')) for name in changed),'historical plan/evidence changed')
    write_json(a.output/'publication-provenance.json',dict(
        schema='CT13PublicationProvenanceV1',source_revision=rev,source_clean=True,
        archive=str(a.root),mount=mount,runtime=runtime,
        physical_manifest_sha256=sha(manifest),evaluated_seal_sha256=sha(ev/'completed-artifacts.json'),
        changed_files={name:sha(source/name) for name in changed if (source/name).is_file()},
        driver_sha256=sha(Path(__file__)),threads=THREADS,
        governing='pinned CT9,not rebuilt or advanced; no new component enabled',
        history='CT4-CT12 seals,failed classifications and installations unchanged'))
    # Preserve complete bodies and selected failure/guard witnesses. Heavy rows,
    # CSR,solutions,states and histories remain at their hashed archive paths.
    records={}
    for relative in ('build-v1','build-v2','build-v3','build-v4','build-v5','build-v6',
                     'exports-v1','exports-v2','orderings-v1','orderings-v3','orderings-v4',
                     'memory-v1','memory-v2','initial-image-v1'):
        directory=a.root/relative
        for name in ('request-manifest.json','results.json','progress.json','live-work-witness.json'):
            path=directory/name
            if path.is_file():
                record=dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))
                if name=='results.json' and relative.startswith('build-'):
                    record['body']=json.loads(path.read_text())
                records[relative+'/'+name]=record
    old=json.loads((source/'docs/plans/cao-tamura-accepted-12-evidence/closed-publication-accounting.json').read_text())
    write_json(a.output/'forecast-before-final-closure.json',dict(
        schema='CT13IndependentLongCampaignForecastV1',endpoint=.407,
        campaign=old['campaign'],
        CT13_spent_CPU_before_publication=cpu,CT13_active_wall_upper_before_publication=wall,
        conservative_total_CPU_before_final=old['campaign']['CPU_forecast_spent_plus_remaining']+cpu,
        conservative_total_wall_before_final=old['campaign']['active_wall_upper_spent_plus_remaining']+wall,
        conservative_output_planning_upper_before_final=old['campaign']['output_planning_upper_bytes']+output,
        CT12_closed_accounting_sha256=sha(source/'docs/plans/cao-tamura-accepted-12-evidence/closed-publication-accounting.json'),
        governing_speedup_demonstrated=False,physical_continuation=False,
        scope='add closed publication costs/reserve separately; preserve25percent campaign headroom,no assumed ordering speedup or online output savings'))
    write_json(a.output/'external-artifact-index.json',dict(
        schema='CT13SelectedExternalArtifactsV1',files=records,
        complete_bindings='evaluated exports/orderings/memory/initial-image and their immutable request/completion seals; no heavy payload copied into Git'))
    shutil.copy2(a.ct12/'static-v1/paired-definition.json',a.output/'selected-paired-definition.json')
    shutil.copy2(Path(__file__),a.output/'publication-driver.py')
    archive=a.output/'ct13-compact-v1.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=6) as tar:
        tar.add(ev,arcname='evaluated-v1')
        for f in sorted(a.output.iterdir()):
            if f!=archive:tar.add(f,arcname=f.name)
    write_json(a.output/'transfer-sha256.json',dict(file=archive.name,bytes=archive.stat().st_size,sha256=sha(archive)))
    print(json.dumps(dict(status='pass',source=rev,CPU_before=cpu,active_wall_before=wall,
        project=peak,aggregate=aggregate,transfer=json.loads((a.output/'transfer-sha256.json').read_text()))),flush=True)

if __name__=='__main__':main()

