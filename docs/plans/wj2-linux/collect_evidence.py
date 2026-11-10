#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only source/runtime audit; write a compact WJ-2 evidence index."""
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

root=Path(sys.argv[1]).resolve();source=root/'source';build=root/'build'
reviewed='7a06e1c1f7e4e81d42f53f5c5d94bea33668bb85'
required='''WJKernelContracts WJIndependentOracle RANSKernelUnitTests RANSCLSEPAWFM1Oracle
RANSCLSEPAWFM2Oracle RANSAcceptedAWFKernel RANSEquilibriumProfilesNative RANSProfileBranchChecks
RANSKBoundQPChecks RANSHistoryRejectionChecks RANSAcceptedAWFCases RANSStressTangentOracle
RANSReturnedFieldsChecks RANSObservationKernels RANSWallTimingChecks RANSEquilibriumTimeNative
RANSReducedModelChecks RANSStandardTransportSmoke RANSCLSTransportSmoke RANSAWFOriginalSmoke
RANSAWFEquilibriumSmoke RANSStandardAWFEquilibriumSmoke RANSEquilibriumWallEvolutionSmoke
RANSEquilibriumWallStageSmoke WJNativeEP WJNativeBudget RANSWJSharedLegacyTransaction'''.split()

def command(*args):return subprocess.check_output(args,text=True).strip()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def identity(path):return dict(path=str(path),sha256=sha(path),bytes=Path(path).stat().st_size)
def vector_hash(values):return hashlib.sha256(struct.pack(f'<{len(values)}d',*values)).hexdigest()
def object_hash(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def transaction_receipt(path):
    data=json.loads(path.read_text(),parse_int=float);offsets=[int(i) for i in data['offsets']]
    attempts=[]
    for a in data['pimple_attempts']:
        p=Path(a['k_snapshot']);stages=json.loads(p.read_text(),parse_int=float)
        attempts.append(dict(time=a['time'],alpha=a['alpha'],scalar_order=a.get('scalar_order'),
            counts={k:a[k] for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves',
                'k_transport_solves','bound_calls','additional_governing_equation_solves','bound_activations')},
            epsilon_frozen_defect=a['epsilon_frozen_linear_residual'],k_frozen_defect=a['k_frozen_linear_residual'],
            epsilon_final_trace_defect=a['refreshed_trace_defect'],k_bound=a['k_bound'],
            correction_budgets={k:a[k] for k in ('k_corrected_frozen_linear_residual','k_added_frozen_residual','k_added_storage_residual')},
            stage_snapshot=identity(p),stage_history_hashes={k:vector_hash(stages[k]) for k in
                ('committed','post_flow','pre_k','raw','corrected','history')},
            transport_matrix=identity(stages['matrix_file'])))
    checkpoint=data['checkpoint'];role='published accepted checkpoint'
    if not checkpoint:
        info=ET.parse(path.with_name('case.xml')).getroot().findall('./CONDITIONS/SOLVERINFO/I')
        checkpoint=next(i.get('VALUE') for i in info if i.get('PROPERTY')=='RANSWallRestart')
        role='unchanged restart checkpoint; rejected attempt did not publish a candidate'
    cp=Path(checkpoint)/'State.bin'
    return dict(report=identity(path),time=data['time'],accepted_steps=data['steps'],
        state=vector_hash(data['state']),pressure_and_gauge=vector_hash(data['state'][offsets[4]:]),
        native_history=[vector_hash(v) for v in data['native_history']],native_history_times=vector_hash(data['native_history_times']),
        accepted_audits=object_hash(data['state_audits']),physical_samples=object_hash(data['physical_samples']),
        last_stage_history=vector_hash(data['last_stage_history']),last_stage_time=data['last_stage_time'],
        rejected_state=vector_hash(data['rejected_state']),failure=data['failure'],
        checkpoint=dict(**identity(cp),role=role),attempts=attempts)

entries=subprocess.check_output(['git','-C',str(source),'ls-tree','-rz','HEAD']).split(b'\0')
content=[];mismatches=[]
for entry in entries:
    if not entry:continue
    header,name=entry.split(b'\t',1);mode,kind,oid=header.split();name=os.fsdecode(name);path=source/name
    raw=os.fsencode(os.readlink(path)) if mode==b'120000' and path.is_symlink() else path.read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if actual!=oid.decode():mismatches.append(name)
    content.append(f'{hashlib.sha256(raw).hexdigest()}  {name}\n')
assert not mismatches,mismatches
(root/'source-content.sha256').write_text(''.join(sorted(content)))
changed=command('git','-C',str(source),'diff','--name-only',reviewed,'HEAD').splitlines()
(root/'changed-source.sha256').write_text(''.join(f'{sha(source/p)}  {p}\n' for p in changed))

inventory=json.loads((root/'test-inventory.json').read_text());registered={t['name']:t for t in inventory['tests']}
assert set(required)<=registered.keys(),sorted(set(required)-registered.keys())
tests=[];ctest_logs=[]
for name,selected in (('ctest-native-final-2.log',required[-3:]),('ctest-regressions-final-2.log',required[:-3])):
    path=root/name;log=path.read_text()
    rows=re.findall(r'Test\s+#\d+:\s+(\w+)\s+\.+\s+(Passed|\*\*\*Failed)\s+([\d.]+) sec',log)
    assert rows and '100% tests passed, 0 tests failed' in log,log[-2000:]
    tests.extend(dict(name=t,status=status,seconds=float(seconds)) for t,status,seconds in rows)
    ctest_logs.append(dict(**identity(path),exit_code=0,
        command=['ctest','--test-dir',str(build),'-j','1','--output-on-failure','-R','^('+'|'.join(selected)+')$'],
        seconds=float(re.search(r'Total Test time \(real\) =\s+([\d.]+) sec',log)[1]),
        receipt_basis='completed serial CTest output and successful enclosing shell',tests=[r[0] for r in rows]))
assert set(required)=={t['name'] for t in tests}

executables='''IncNavierStokesSolver RANSWallResolutionProbe WJKernelProbe RANSKernelProbe
RANSConstitutiveProbe RANSAcceptedAWFProbe RANSEquilibriumProfilesProbe RANSProfileBranchChecks
RANSKBoundQPChecks RANSHistoryRejectionChecks RANSStressTangentAudit RANSReturnedFieldsProbe
RANSObservationKernelProbe RANSWallTimingChecks RANSEquilibriumTimeProbe RANSReducedModelChecks'''.split()
binaries={};libraries={}
for name in executables:
    path=build/name;ldd=command('ldd',str(path));binaries[name]=dict(**identity(path),ldd=ldd)
    for lib in re.findall(r'(/\S+)\s+\(0x[0-9a-fA-F]+\)',ldd):
        resolved=Path(lib).resolve();libraries[str(resolved)]=sha(resolved)
for name in ('IncNavierStokesSolver','RANSWallResolutionProbe'):
    assert str(build/'geometry-lib') in binaries[name]['ldd'],'isolated geometry dependency not selected'
    symbols=command('nm','-C','--defined-only',str(build/name))
    selected=[line for line in symbols.splitlines() if re.search(
        r'\bT Nektar::SolverUtils::(EvaluateRANSClosure\(|EvaluateWJClosure\(|ApplyWJClosureDirection\(|EvaluateRANSRobustClosure\(|RANSNeutralReference::WJ\()',line)]
    assert any('EvaluateWJClosure(' in line for line in selected)
    assert any('EvaluateRANSClosure(' in line for line in selected)
    binaries[name]['executable_defined_closure_symbols']=selected

native={};artifacts={};earlier=[]
for policy in ('EP','Budget'):
    summaries=list((build/f'wj2-native-{policy}').glob('run-*/summary.json'))
    passing=[p for p in summaries if json.loads(p.read_text())['status']=='PASS']
    latest=max(passing,key=lambda p:p.stat().st_mtime);report=json.loads(latest.read_text());folder=latest.parent
    assert report['source_revision']==command('git','-C',str(source),'rev-parse','HEAD')
    native[policy]=dict(**identity(latest),report=report,
        transactions={name:transaction_receipt(folder/name/'case.evolution.json') for name in
            ('whole2','split1','restart1','accepted-zero-step','reject-after-hard-audit')})
    # Bind all generated XML, JSON, binary checkpoints and transport matrices.
    # Large content stays under the immutable Linux run directory.
    for p in sorted(folder.rglob('*')):
        if p.is_file() and (p.suffix in ('.xml','.json','.bin','.f64') or p.name=='State.bin'):
            artifacts[str(p)]=dict(sha256=sha(p),bytes=p.stat().st_size)
    for p in summaries:
        if p!=latest:earlier.append(dict(**identity(p),status=json.loads(p.read_text())['status']))
legacy_reports=list((build/'wj2-legacy-transaction').glob('run-*/summary.json'))
legacy_path=max(legacy_reports,key=lambda p:p.stat().st_mtime);legacy=json.loads(legacy_path.read_text())
assert legacy['status']=='passed' and legacy['restart_exact'] and legacy['rollback_exact']
for p in legacy_path.parent.rglob('*'):
    if p.is_file() and p.suffix in ('.json','.xml','.bin','.f64'):
        artifacts[str(p)]=dict(sha256=sha(p),bytes=p.stat().st_size)
write(root/'native-artifact-index.json',artifacts)

flags={name:(build/'CMakeFiles'/f'{name}.dir'/'flags.make').read_text()
    for name in ('RANSKernels','IncNavierStokesSolver','RANSWallResolutionProbe','RANSSpatialDomains','WJKernelProbe')}
assert all('-ffp-contract=off' in flag and '-fno-fast-math' in flag for flag in flags.values())
snapshot=root/'processes-pre-final-build.txt'
worker_rows=[];worker_snapshots={}
for p in (snapshot,root/'processes-during-final-tests.txt'):
    workers=0
    for row in p.read_text().splitlines()[1:]:
        pid,nlwp,comm=row.split(None,2)
        if re.fullmatch(r'IncNavierStokes.*|cc1plus|cc1|cmake|ctest|python3|orted|mpirun|RANS.*',comm):
            workers+=int(nlwp)
            if p==snapshot:worker_rows.append(row)
    worker_snapshots[str(p)]=dict(**identity(p),relevant_threads=workers)
assert max(q['relevant_threads'] for q in worker_snapshots.values())+4<=72
manifest=dict(schema='WJ2EvidenceV1',reviewed_parent=reviewed,
    tested_source_revision=command('git','-C',str(source),'rev-parse','HEAD'),
    tested_git_tree=command('git','-C',str(source),'rev-parse','HEAD^{tree}'),
    source_files=len(content),source_content=identity(root/'source-content.sha256'),
    changed_source=identity(root/'changed-source.sha256'),source_content_mismatches=mismatches,
    working_tree_ignoring_ntfs_modes=command('git','-c','core.filemode=false','-C',str(source),'status','--porcelain'),
    specs={p:sha(source/'docs/literature'/p) for p in
        ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')},
    run_root=str(root),host=platform.node(),platform=platform.platform(),compiler=command('/usr/bin/c++','--version'),
    cmake=command('cmake','--version'),python=sys.version,
    floating_point_policy='Release -O3 -DNDEBUG -ffp-contract=off -fno-fast-math; C++17',flags=flags,
    configure_command=['cmake','-S',str(source/'solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(build),
        '-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON',
        '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake',
        '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math'],
    build_command=['cmake','--build',str(build),'--parallel','1','--target']+executables,
    allocation=dict(build_workers=1,ctest_parallel=1,native_ranks=1,threads={k:os.environ.get(k) for k in
        ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        process_snapshots=worker_snapshots,relevant_rows_before_build=worker_rows,
        ompi_tmpdir_base=os.environ.get('OMPI_MCA_orte_tmpdir_base')),
    inventory=dict(**identity(root/'test-inventory.json'),registered=len(registered),missing=[],
        selected_properties={k:registered[k].get('properties',[]) for k in required}),
    ctest_receipts=ctest_logs,tests=tests,executables=binaries,linked_libraries_sha256=libraries,
    command_receipts=dict(**identity(root/'final-command-receipts.json'),
        receipts=[dict(**q,log_sha256=sha(q['log'])) for q in json.loads((root/'final-command-receipts.json').read_text())],
        native_artifact_index_sha256=sha(root/'native-artifact-index.json'),
        independent_oracle_sha256=sha(build/'wj1-oracle.json')),
    static_kernels=identity(build/'libRANSKernels.a'),
    independent_oracle=identity(build/'wj1-oracle.json'),
    native=native,legacy_transaction=dict(**identity(legacy_path),report=legacy),
    artifact_index=identity(root/'native-artifact-index.json'),earlier_native_attempts=earlier,
    logs={p.name:identity(p) for p in sorted(root.glob('*.log')) if p.name!='collection.log'},
    scope='native serial flat 2D short steps and transactions; no long-time/global stability or 3D/MPI claim')
write(root/'manifest.json',manifest)
print(json.dumps({k:manifest[k] for k in ('tested_source_revision','source_files','source_content_mismatches','working_tree_ignoring_ntfs_modes')},indent=2))
