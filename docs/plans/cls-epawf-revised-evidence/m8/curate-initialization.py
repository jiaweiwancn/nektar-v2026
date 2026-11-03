"""Copy small, source-pinned M8 verification evidence; never runs a solver."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess

HOST = 'jay@100.105.247.110'
REMOTE = '/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/initialization'
DEST = Path('/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m8')
FILES = {
    'mms-initialization-attempt-001.json': 'mms-runner-001/receipt.json',
    'mms-initialization-prefix-001.json': 'mms-001/validation.json',
    'mms-initialization-launcher-failure.log': 'mms-001/flat-mpi2/run.log',
    'mms-initialization-runner.json': 'mms-runner-002/receipt.json',
    'mms-initialization-qualification.json': 'mms-002/validation.json',
    'mms-checkpoint-compatibility.json': 'mms-runner-002/checkpoint-compatibility-signatures-003.json',
    'mms-initialization-mpi-bindings.log': 'mms-002/flat-mpi2/run.log',
    'mms-stale-checkpoint-rejection.log': 'mms-runner-002/reject-old-functions-capped.log',
    'run-initialization-mms-001.py': 'run-corrected-mms-verification.py',
    'run-initialization-mms-002.py': 'run-corrected-mms-verification-002.py',
    'audit-initialization-compatibility-001.py': 'audit-checkpoint-compatibility-signatures.py',
    'audit-initialization-compatibility-002.py': 'audit-checkpoint-compatibility-signatures-002.py',
    'audit-initialization-compatibility-003.py': 'audit-checkpoint-compatibility-signatures-003.py',
}
FP = (
    'mms-runner-001/runtime-before.json', 'mms-runner-001/runtime-after.json',
    'mms-runner-002/runtime-before.json', 'mms-runner-002/runtime-after.json',
    'mms-runner-002/signature-audit-runtime-before.json',
    'mms-runner-002/signature-audit-002-runtime-before.json',
    'mms-runner-002/signature-audit-002-runtime-after.json',
    'mms-runner-002/signature-audit-003-runtime-before.json',
    'mms-runner-002/signature-audit-003-runtime-after.json',
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def remote_manifest():
    code = '''from pathlib import Path
import hashlib,json,subprocess
root=Path(ROOT)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
source=root/'source'
head=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
status=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
assert head=='ed34052fe4dd967e963c22272a0ec0888cf35b70' and not status
files={name:dict(original=str(root/path),bytes=(root/path).stat().st_size,sha256=sha(root/path)) for name,path in FILES.items()}
fps={}
for path in FP:
 p=root/path;raw=json.loads(p.read_text())
 fps[path]=dict(original=str(p),sha256=sha(p),success=raw['success'],source_commit=raw['source_commit'],
  source_input_count=len(raw['source_inputs']),runtime_image_count=len(raw['runtime_artifacts']),
  runtime_environment=raw['runtime_environment'],failures=raw.get('failures',[]),
  identical_qualification_inputs_and_runtime=raw.get('identical_qualification_inputs_and_runtime'),
  differences=raw.get('differences'))
print(json.dumps(dict(source_commit=head,source_clean=True,files=files,fingerprints=fps)))
'''.replace('ROOT', repr(REMOTE)).replace('FILES', repr(FILES)).replace('FP', repr(FP))
    result = subprocess.run(['ssh', HOST, 'python3 -'], input=code, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return json.loads(result.stdout)


before = remote_manifest()


def copy(item):
    name, row = item
    target = DEST/name
    if target.exists():
        assert sha(target) == row['sha256'], name
    else:
        subprocess.run(['scp', HOST+':'+row['original'], str(target)], check=True)
    assert target.stat().st_size == row['bytes'] and sha(target) == row['sha256'], name
    return dict(file=name, transformation='byte-for-byte copy', **row)


with ThreadPoolExecutor(max_workers=4) as pool:
    copied = list(pool.map(copy, before['files'].items()))
after = remote_manifest()
assert after == before
index = dict(schema=1, source_commit=before['source_commit'], source_clean=True,
             scope='Compact corrected-initialization MMS and scoped checkpoint compatibility evidence',
             connection='Tailscale 100.105.247.110', files=copied,
             remote_inputs_and_fingerprints_rechecked_after_copy=True,
             curation_recipe_sha256=sha(__file__))
(DEST/'mms-curation-index.json').write_text(json.dumps(index, indent=2)+'\n')
metadata = dict(schema=1, fingerprints=before['fingerprints'],
    qualification=dict(status='passed', fresh_cases=7, all_source_runtime_old_inputs_unchanged=True,
        note='mms-002 executes every solver case afresh after the source fix; no report replay'),
    preserved_attempts=[
        dict(id='mms-001', stage='MPI launcher', status='stopped_before_MPI_solver_initialization',
            reason='Explicit --map-by core conflicts with OMPI_MCA_hwloc_base_cpu_list=40-43; rank-file mapping policy conflict',
            scientific_scope='Four serial MMS cases passed; no MPI solver state was initialized. All seven cases were rerun fresh in mms-002.'),
        dict(id='signature-audit-001', stage='read-only metadata audit', status='assertion_guard_stopped',
            reason='Initial LD-path comparison did not drop the actual trailing empty component, whereas the fingerprint tool normalizes it away.',
            output_receipt_emitted=False),
        dict(id='signature-audit-002', stage='read-only metadata audit', status='post-fingerprint_guard_stopped',
            reason='Importing cls_epawf_m8_oracle changes VECLIB_MAXIMUM_THREADS from unset to 1; full post-fingerprint correctly rejected this environment difference.',
            scientific_scope='Exact signatures, field/history norms and all nine actual-case loader checks passed before the metadata guard; no CFD rerun was performed.',
            output_receipt_emitted=False),
        dict(id='signature-audit-003', stage='read-only metadata audit', status='passed',
            reason='VECLIB_MAXIMUM_THREADS is initialized explicitly before the audit pre-pin. Its environment difference from the outer solver runner is recorded, and its own complete pre/post fingerprints are identical. Actual duplicate/empty LD paths are retained and loader bindings are checked from every actual solver case directory.')])
(DEST/'mms-initialization-environment-audits.json').write_text(json.dumps(metadata, indent=2)+'\n')
affinity = dict(schema=1, classification='nonnumerical MPI launcher affinity diagnostic',
    evidence_origin='Exact command and observed output from this task tool transcript; saved after execution',
    command=['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'OMPI_MCA_hwloc_base_cpu_list=40-43',
             'mpirun', '--report-bindings', '--bind-to', 'core', '-np', '2', 'python3', '-c',
             'import os; print(os.sched_getaffinity(0))'],
    explicit_environment_overrides=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                                      OMPI_MCA_hwloc_base_cpu_list='40-43'),
    inherited_environment_snapshot_captured=False,
    observed_stdout=['{40, 41, 42, 43}', '{40, 41, 42, 43}'],
    observed_stderr=['[jay-Super-Server:34869] MCW rank 0 is not bound (or bound to all available processors)',
                     '[jay-Super-Server:34869] MCW rank 1 is not bound (or bound to all available processors)'],
    interpretation='Two launcher-control processes imported only os and printed affinity. No model/solver or numerical experiment was executed. The proposed adapter was not created or executed; the final CFD matrix uses the unchanged driver after unsetting the conflicting CPU-list variable.')
(DEST/'mms-affinity-diagnostic.json').write_text(json.dumps(affinity, indent=2)+'\n')
print(json.dumps(dict(copied_files=len(copied), copied_bytes=sum(row['bytes'] for row in copied),
                     source_commit=before['source_commit'], destination=str(DEST))))
