#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only final provenance, accepted-state extrema and resource supplement."""
import argparse,json,os,platform,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
import scipy
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha


def read(path):return json.loads(path.read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct2','compact','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();source=a.root/'source';git=['git','-c','index.threads=1','-c','core.preloadIndex=false','-c','pack.threads=1']
    revision=subprocess.check_output(git+['rev-parse','HEAD'],cwd=source,text=True).strip()
    seals=read(a.compact/'completed-artifacts.json')
    require(all(sha(a.compact/f)==h for f,h in seals['files'].items()),'changed_compact_evidence')
    m=RequestManifest(a.output,dict(schema='CT3FinalReadOnlySupplementV1',compact_seal_sha256=sha(a.compact/'completed-artifacts.json'),
        source_revision=revision,source_sha256=sha(Path(__file__))))
    previous=read(source/'docs/plans/cao-tamura-accepted-2-evidence/final-provenance.json')
    specs={f:sha(source/f) for f in previous['specifications']}
    data={f:sha(source/'docs/plans/cao-tamura-validation-data'/f) for f in previous['frozen_dataset']['hashes']}
    kernels={f:sha(source/f) for f in previous['unchanged_reused_kernel_sources']}
    require(specs==previous['specifications'],'changed_mathematical_specification')
    require(data==previous['frozen_dataset']['hashes'],'changed_frozen_observations_or_roles')
    require(kernels==previous['unchanged_reused_kernel_sources'],'changed_reused_kernel')
    cache=a.root/'profile-build-v1/build/CMakeCache.txt'
    compiler={line.split(':',1)[0]:line.split('=',1)[1] for line in cache.read_text().splitlines()
        if '=' in line and line.startswith(('CMAKE_CXX_','CMAKE_BUILD_TYPE:','RANS_','Nektar++_DIR:'))}
    flags=a.root/'profile-build-v1/build/CMakeFiles/IncNavierStokesSolver.dir/flags.make'
    provenance=dict(schema='CT3FinalProvenanceV1',starting_revision='db37c6b11882e7a06413a29eb6623ab212c4d4a5',
        CT2_final_source=previous.get('source_revision_before_evidence_commit'),source_revision=revision,
        source_status=subprocess.check_output(git+['status','--porcelain'],cwd=source,text=True),
        specifications=specs,frozen_dataset=data,unchanged_kernels=kernels,
        environment=dict(python=sys.version,python_executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,
            platform=platform.platform(),threads={k:os.environ.get(k) for k in THREADS},
            compiler=subprocess.check_output([compiler['CMAKE_CXX_COMPILER'],'--version'],text=True).splitlines()[0],
            CMake_cache=compiler,CMake_cache_sha256=sha(cache),actual_target_flags=flags.read_text(),actual_target_flags_sha256=sha(flags),
            actual_mount=json.loads(subprocess.check_output(['findmnt','-T',str(a.root),'-J'],text=True)),
            MPI='MPI-enabled build; all supported native tests one serial rank; library/helper threads counted'),
        runtime_scope=dict(primary_source='a897896524b8761d3a954d491f96d1b33941d94f',
            primary_runtime_bridge_sha256=sha(a.ct2/'relocation-runtime-bridge.json'),
            primary_runtime_build_identity_sha256=sha(source/'docs/plans/cao-tamura-accepted-2-evidence/runtime-build-identity.json'),
            timing_build_request=read(a.root/'profile-build-v1/request-manifest.json'),
            primary_and_profile_loaded_libraries='exact per-process loaded-image hashes in execution objects; linked-image hashes in source-scope/timing-build; transient MPI mappings are retained as observations, not persistent libraries'),
        source_files={str(f.relative_to(source)):sha(f) for f in (source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('ct_accepted3*.py')})
    write_json(a.output/'final-provenance.json',provenance)
    extrema={};summaries={}
    keys=('accepted_awf_covariance','accepted_awf_scalar_domain','physical_sample_K','physical_sample_outer',
        'native_A95_radicand_range','native_raw_divergence_range','native_outer_isotropic_production_work_maximum',
        'common_velocity_k_dependent_epsilon_trace_defects','wall_cutoff_required_budget',
        'strong_continuity','momentum_k_load_and_boundary_reaction')
    witness_keys=('native_A95_minimum_radicand_point','native_A95_minimum_radicand_geometry_element_side_block',
        'native_raw_divergence_maximum_witness','native_outer_isotropic_production_work_witness',
        'accepted_awf_inner_minimum_point','accepted_awf_outer_minimum_point',
        'native_inner_minimum_geometry_element_side_block','native_outer_minimum_geometry_element_side_block')
    for case in read(a.compact/'decision.json')['cases']:
        for target in (.008,.02,.2):
            file=a.compact/(case+'-T'+str(target)+'.json');receipt=read(file)
            for dt,row in receipt['members'].items():
                folder=Path(row['directory']);native=read(folder/'case.evolution.json');label=case+'/'+str(target)+'/'+dt
                require(sha(folder/'case.evolution.json')==row['native_report_sha256'],'changed_native_endpoint')
                values={}
                for key in keys:
                    bounds=[]
                    for t,audit in zip(native['audit_times'],native['state_audits'],strict=True):
                        v=audit.get(key)
                        if not isinstance(v,list):continue
                        while len(bounds)<len(v):bounds.append({})
                        for j,x in enumerate(v):
                            require(isinstance(x,(float,int)) and np.isfinite(x),'nonfinite_required_native_audit')
                            b=bounds[j]
                            for mode in ('min','max'):
                                if mode not in b or (x<b[mode]['value'] if mode=='min' else x>b[mode]['value']):
                                    b[mode]=dict(value=x,time=t,witnesses={k:audit[k] for k in witness_keys if k in audit})
                    values[key]=bounds
                extrema[label]=dict(native_report_sha256=sha(folder/'case.evolution.json'),
                    scope='all accepted-state native inventories in this sealed segment; component order is the native audit schema; rejected-state and staged inventories separate',
                    audited_times=native['audit_times'],component_extrema=values,final_audit=native['state_audits'][-1])
                if row['failure']:continue
                defects={n:max(abs(q['equations'][n]['frozen_stage_physical_defect'])/q['equations'][n]['physical_scale'] for q in row['budgets'].values())
                    for n in ('inner_k','outer_k','outer_epsilon','whole_k')}
                summaries[label]=dict(time=row['time'],mass=row['mass'],maximum_absolute_scaled_budget_defect=defects,
                    endpoint_domain=row['domain'],epsilon_trace_and_cutoff=native['state_audits'][-1]['common_velocity_k_dependent_epsilon_trace_defects'],
                    budget_status=row['budget_gate']['status'])
    write_json(a.output/'accepted-state-extrema.json',extrema);write_json(a.output/'physical-summary.json',summaries)
    resources={};pending=[];events=[];highwater=[]
    for folder in sorted(a.root.glob('launch-*')):
        final=folder/'resource-final.json'
        if final.exists():
            resources[folder.name]=dict(read(final),execution=read(folder/'execution.json'))
        else:pending.append(str(folder))
        for name,target in [('worker-guard-events.jsonl',events),('worker-guard-highwater.jsonl',highwater),('project-highwater.jsonl',highwater)]:
            if (folder/name).exists():target.extend(dict(control=folder.name,kind=name,**json.loads(line)) for line in (folder/name).read_text().splitlines())
    result=dict(schema='CT3CompletedFamilyResourcesV1',status='pass' if resources and all(q['aggregate_status']==q['project_status']=='pass' for q in resources.values()) else 'fail',
        maximum_aggregate=max(q['maximum_aggregate'] for q in resources.values()),maximum_project=max(q['maximum_project'] for q in resources.values()),
        project_limit=16,aggregate_limit=72,pending_finalization_monitors=pending,
        scope='completed families, including failed collection attempts; copy this finalizer own monitor receipt after it exits',
        resources=resources,events=events,highwater=highwater)
    write_json(a.output/'completed-family-resources.json',result)
    m.seal(a.output,['final-provenance.json','accepted-state-extrema.json','physical-summary.json','completed-family-resources.json'])
    print(json.dumps(dict(resources=result['status'],maximum_aggregate=result['maximum_aggregate'],maximum_project=result['maximum_project'],
        pending_finalization_monitors=pending,unchanged_specifications_data_kernels=True,successful_endpoint_members=len(summaries))))
    return 0
if __name__=='__main__':raise SystemExit(main())
