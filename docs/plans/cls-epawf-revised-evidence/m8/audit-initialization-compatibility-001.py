"""Read-only exact checkpoint/signature and native-history compatibility audit."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
from datetime import datetime, timezone

M8 = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
ROOT = M8/'initialization'
SOURCE = ROOT/'source'
OUTPUT = ROOT/'mms-runner-002/checkpoint-compatibility-signatures.json'
REVISION = 'ed34052fe4dd967e963c22272a0ec0888cf35b70'
MARKER = ':capped_initial_trace=1'
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
os.sched_setaffinity(0, {40, 41, 42, 43})
lock = (ROOT/'validation.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(SOURCE), *args], text=True).strip()


assert not OUTPUT.exists()
assert git('rev-parse', 'HEAD') == REVISION and not git('status', '--porcelain')
producer_path = ROOT/'mms-runner-002/receipt.json'
producer = json.loads(producer_path.read_text())
assert producer['status'] == 'passed' and producer['full_runtime_and_source_unchanged']
assert sha(ROOT/'build/IncNavierStokesSolver') == producer['program']['sha256']
fingerprint_path = ROOT/'mms-runner-002/runtime-before.json'
fingerprint = json.loads(fingerprint_path.read_text())
for key, value in fingerprint['runtime_environment'].items():
    if value is None:
        os.environ.pop(key, None)
    else:
        os.environ[key] = value
os.environ['TMPDIR'] = producer['environment']['TMPDIR']
fingerprint_tool = M8.parent/'transfers/fingerprint-runtime.py'
audit_before = OUTPUT.with_name('signature-audit-runtime-before.json')
audit_after = OUTPUT.with_name('signature-audit-runtime-after.json')
assert not audit_before.exists() and not audit_after.exists()
subprocess.run(['python3', '-B', str(fingerprint_tool), '--root', str(ROOT),
                '--output', str(audit_before), '--compare', str(fingerprint_path)], check=True)
audit_before_data = json.loads(audit_before.read_text())
assert audit_before_data['success'] and audit_before_data['identical_qualification_inputs_and_runtime']
launcher_pin = fingerprint['tool_images']['mpiexec']
launcher_aliases = []
for name in ('mpirun', 'mpiexec'):
    path = Path(shutil.which(name))
    assert sha(path) == launcher_pin['sha256'] and str(path.resolve()) == launcher_pin['resolved']
    launcher_aliases.append(dict(name=name, path=str(path), resolved=str(path.resolve()), sha256=sha(path)))
scripts = SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS'
sys.path.insert(0, str(scripts))
import numpy as np
from validate_cls_epawf_m4 import state_difference, native_history_difference
from cls_epawf_m8_mms import verify as verify_MMS

pins = {}


def load(path):
    path = Path(path)
    pins[str(path)] = sha(path)
    return json.loads(path.read_text())


def signature_record(report_path):
    raw = load(report_path)
    text = raw['signature']
    xml_path = Path(report_path).parent/'case.xml'
    pins[str(xml_path)] = sha(xml_path)
    result = dict(report=str(report_path), report_sha256=sha(report_path),
                  session=str(xml_path), session_sha256=sha(xml_path),
                  signature=text, signature_utf8_sha256=hashlib.sha256(text.encode('utf-8')).hexdigest(),
                  capped_initial_trace_marker_present=MARKER in text,
                  time=raw['time'], order=raw['order'], dt=raw['dt'])
    return raw, result


def checkpoint_record(report_path):
    raw, result = signature_record(report_path)
    payload_path = Path(raw['checkpoint'])/'State.bin'
    payload = payload_path.read_bytes()
    pins[str(payload_path)] = sha(payload_path)
    magic, endian, schema, count = struct.unpack_from('=QQQQ', payload)
    assert (magic, endian, schema) == (0x4550415746544931, 0x0102030405060708, 1)
    assert 32 + count <= len(payload) - 8
    text = payload[32:32+count].decode('utf-8')
    assert text == raw['signature']
    fnv = 14695981039346656037
    for value in payload[:-8]:
        fnv = ((fnv ^ value) * 1099511628211) & ((1 << 64) - 1)
    checksum, = struct.unpack_from('=Q', payload, len(payload)-8)
    assert checksum == fnv
    result.update(checkpoint_payload=str(payload_path), checkpoint_payload_sha256=sha(payload_path),
                  checkpoint_payload_bytes=len(payload), binary_signature_matches_report=True,
                  independently_verified_FNV64_checksum=True,
                  checkpoint_format=dict(magic=hex(magic), endian=hex(endian), schema=schema,
                                         signature_bytes=count))
    return raw, result


def compare(a, b, name):
    assert a['dimension'] == b['dimension'] and a['order'] == b['order'] and a['dt'] == b['dt']
    assert a['offsets'] == b['offsets']
    assert len(a['native_history']) == len(b['native_history'])
    assert all(len(x) == len(y) for x, y in zip(a['native_history'], b['native_history']))
    fields = state_difference(a, b, name+'_full_fields', tolerance=1e-8)
    histories = native_history_difference(a, b, name+'_complete_native_history', tolerance=1e-8)
    auxiliary = {}
    for key in ('native_history_times', 'fixed_reference', 'fixed_defect'):
        first, second = np.asarray(a[key]), np.asarray(b[key])
        assert first.shape == second.shape
        maximum = float(np.max(np.abs(first-second))) if first.size else 0.
        tolerance = 1e-8 if key == 'native_history_times' else 1e-10
        assert maximum <= tolerance, (name, key, maximum)
        auxiliary[key+'_maximum_absolute_difference'] = maximum
    return dict(field_mass_norms=fields, maximum_field_mass_norm=max(fields),
                history_mass_norms=histories, maximum_history_mass_norm=max(histories),
                native_history_slots=len(a['native_history']),
                history_slot_vector_sizes=[len(slot) for slot in a['native_history']],
                field_and_history_tolerance=1e-8, auxiliary=auxiliary,
                all_shapes_checked=True)


old_functions, old_f = checkpoint_record(M8/'mms-completion-001/flat-checkpoint/case.evolution.json')
new_functions, new_f = checkpoint_record(ROOT/'mms-002/flat-checkpoint/case.evolution.json')
fresh_restart, new_restart = signature_record(ROOT/'mms-002/flat-restart/case.evolution.json')
fresh_full, new_full = signature_record(ROOT/'mms-002/flat/case.evolution.json')
fresh_MPI, new_MPI = signature_record(ROOT/'mms-002/flat-mpi2/case.evolution.json')
old_neutral, old_n = checkpoint_record(M8/'native-completion-001/scaled-native/checkpoint/case.evolution.json')
neutral_full, neutral_reference = signature_record(M8/'native-completion-001/scaled-native/joint-four/case.evolution.json')
neutral_restart, neutral_resumed = signature_record(ROOT/'compatibility-002/old-neutral-compatible/case.evolution.json')
negative_path = ROOT/'compatibility-002/old-functions-capped/case.evolution.json'
negative = load(negative_path)
negative_xml = negative_path.parent/'case.xml'
negative_log = ROOT/'mms-runner-002/reject-old-functions-capped.log'
pins[str(negative_xml)], pins[str(negative_log)] = sha(negative_xml), sha(negative_log)
assert negative['accepted_states'] == [] and negative['initialization_complete'] is False
assert 'checkpoint signature mismatch' in negative['failure']
assert MARKER not in old_f['signature'] and MARKER in new_f['signature']
assert new_f['signature'].count(MARKER) == 1
assert new_f['signature'].replace(MARKER, '') == old_f['signature']
assert new_f['signature'] == new_restart['signature'] == new_full['signature']
assert old_n['signature'] == neutral_reference['signature'] == neutral_resumed['signature']
assert MARKER not in old_n['signature']
validation_path = ROOT/'mms-002/validation.json'
validation = load(validation_path)
build_pin_path = ROOT/'build-002/runtime-fingerprint.json'
build_pin = json.loads(build_pin_path.read_text())
image_keys = ('source_inputs', 'runtime_artifacts', 'runtime_bindings', 'python_packages', 'tool_images')
assert all(fingerprint[key] == build_pin[key] for key in image_keys)
environment_changes = {key: dict(build=build_pin['runtime_environment'].get(key),
                                 qualification=fingerprint['runtime_environment'].get(key))
                       for key in set(build_pin['runtime_environment']) | set(fingerprint['runtime_environment'])
                       if build_pin['runtime_environment'].get(key) != fingerprint['runtime_environment'].get(key)}
actual_solver_ld_path = validation['environment']['LD_LIBRARY_PATH']
normalized_ld_path = ':'.join(dict.fromkeys(actual_solver_ld_path.split(':')))
assert normalized_ld_path == fingerprint['runtime_environment']['LD_LIBRARY_PATH']
commands = []
for origin, rows in ((producer_path, producer['commands']), (validation_path, validation['commands'])):
    for row in rows:
        canonical = json.dumps(dict(command=row['command'], cwd=row['cwd']), sort_keys=True,
                               separators=(',', ':'), ensure_ascii=False).encode('utf-8')
        commands.append(dict(receipt=str(origin), name=row['name'], command=row['command'], cwd=row['cwd'],
                             command_and_cwd_canonical_json_sha256=hashlib.sha256(canonical).hexdigest(),
                             exit_code=row['exit_code'], log_sha256=row['log_sha256']))
result = dict(schema=1, status='passed', generated_utc=datetime.now(timezone.utc).isoformat(),
              scope='Read-only exact signature, checkpoint checksum/hash and complete physical native-history audit',
              source_commit=REVISION, solver=producer['program'],
              producer_receipt=str(producer_path), producer_receipt_sha256=sha(producer_path),
              audit_recipe_sha256=sha(__file__),
              norm_oracle_sha256=sha(scripts/'validate_cls_epawf_m4.py'),
              command_hash_encoding='UTF-8 JSON of {command: argv array, cwd: absolute path}; sort_keys=True, separators=(comma,colon), ensure_ascii=False',
              commands=commands, actual_MPI_launcher_aliases=launcher_aliases,
              actual_MPI_launcher_and_pinned_loader_image_identical=True,
              pinned_MPI_loader_bindings=fingerprint['runtime_bindings']['tool/mpiexec'],
              build_pin_comparison=dict(path=str(build_pin_path), sha256=sha(build_pin_path),
                  identical_fields=list(image_keys), full_runtime_environment_identical=not environment_changes,
                  runtime_environment_changes=environment_changes),
              actual_solver_environment=dict(validation['environment'],
                  **{key: value for key, value in fingerprint['runtime_environment'].items()
                     if key.startswith(('OMPI_MCA_', 'PMIX_MCA_', 'UCX_', 'I_MPI_', 'MPICH_'))}),
              fingerprint_LD_LIBRARY_PATH_normalization=dict(actual_solver=actual_solver_ld_path,
                  fingerprint=fingerprint['runtime_environment']['LD_LIBRARY_PATH'],
                  only_identical_duplicate_directories_removed=True),
              old_functions_checkpoint=old_f, new_functions_checkpoint=new_f,
              new_functions_restart=new_restart,
              stale_functions_control=dict(status='expected_rejection_verified',
                  report=str(negative_path), report_sha256=sha(negative_path),
                  XML_sha256=sha(negative_xml), failure=negative['failure'],
                  log_sha256=sha(negative_log), accepted_states=0, initialization_complete=False,
                  exact_signature_change_only_scoped_marker=True,
                  interpretation='Expected scoped checkpoint versioning rejection; not a physical robustness failure'),
              old_neutral_checkpoint=old_n, old_neutral_reference=neutral_reference,
              new_neutral_restart=neutral_resumed,
              unaffected_neutral_signature_exactly_preserved=True,
              fresh_MMS_checkpoint_analytic_verification=verify_MMS(new_functions),
              fresh_MMS_complete_MPI=compare(fresh_full, fresh_MPI, 'fresh_MMS_MPI'),
              fresh_MMS_complete_restart=compare(fresh_full, fresh_restart, 'fresh_MMS'),
              old_neutral_complete_restart=compare(neutral_full, neutral_restart, 'old_Neutral'),
              source_and_runtime_fingerprint_before_sha256=sha(ROOT/'mms-runner-002/runtime-before.json'),
              source_and_runtime_fingerprint_after_sha256=sha(ROOT/'mms-runner-002/runtime-after.json'),
              input_sha256=pins)
assert all(sha(path) == expected for path, expected in pins.items())
assert git('rev-parse', 'HEAD') == REVISION and not git('status', '--porcelain')
assert sha(ROOT/'build/IncNavierStokesSolver') == producer['program']['sha256']
subprocess.run(['python3', '-B', str(fingerprint_tool), '--root', str(ROOT),
                '--output', str(audit_after), '--compare', str(audit_before)], check=True)
audit_after_data = json.loads(audit_after.read_text())
assert audit_after_data['success'] and audit_after_data['identical_qualification_inputs_and_runtime']
result['readonly_audit_runtime'] = dict(before=str(audit_before), before_sha256=sha(audit_before),
    after=str(audit_after), after_sha256=sha(audit_after),
    identical_to_solver_cohort_before=True, identical_before_after=True,
    command=['python3', '-B', str(Path(__file__).resolve())],
    numerical_threads={name: os.environ[name] for name in
        ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')},
    affinity=sorted(os.sched_getaffinity(0)))
result['all_input_hashes_rechecked_at_end'] = True
OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'], output=str(OUTPUT),
    fresh_MMS_history_max=result['fresh_MMS_complete_restart']['maximum_history_mass_norm'],
    old_neutral_history_max=result['old_neutral_complete_restart']['maximum_history_mass_norm'],
    scoped_marker_only=True, unaffected_neutral_exact_signature=True)), flush=True)
