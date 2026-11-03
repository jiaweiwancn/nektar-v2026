#!/usr/bin/env python3
"""Strict metadata-only sealing/reuse of completed R2 auxiliary analyses.

No numerical library, native executable or scientific oracle is imported.
Sealing observes completed historical files; it does not invent an original
output hash that the earlier wrapper did not record. Reuse requires the seal
and verifies the actual files again. Missing/incompatible inputs are errors.
"""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

HELPER = 'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
METADATA = 'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_reuse.py'
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'runtime_environment', 'python_packages', 'tool_images')
IGNORED_ORCHESTRATION = frozenset(('main', 'verify'))


def require(condition, message):
    if not condition:
        raise ValueError('R2 solve reuse: '+message)


def digest(path):
    output = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024), b''):
            output.update(block)
    return output.hexdigest()


def canonical(value):
    output = hashlib.sha256()
    encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), allow_nan=False)
    for piece in encoder.iterencode(value):
        output.update(piece.encode())
    return output.hexdigest()


def vector_digest(values):
    output = hashlib.sha256()
    for value in values:
        require(not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value),
            'nonfinite/nonnumeric coefficient in identity check')
        output.update(struct.pack('<d', float(value)))
    return output.hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    require(path.is_file(), 'missing file '+str(path))
    return dict(bytes=path.stat().st_size, sha256=digest(path))


def check_files(files):
    for name, expected in files.items():
        require(file_identity(name) == {key:expected[key] for key in ('bytes', 'sha256')},
            'changed file '+name)


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1',
        '-c', 'checkout.workers=1', '-C', str(source), *arguments], text=True).strip()


def option(receipt, name):
    arguments = receipt['command_argv']
    require(arguments.count('--'+name) == 1, 'missing/ambiguous historical --'+name)
    return Path(arguments[arguments.index('--'+name)+1]).resolve()


def prefix_digest(path, count):
    output = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while count:
            block = stream.read(min(1024*1024, count))
            require(bool(block), 'native trace prefix was truncated')
            output.update(block); count -= len(block)
    return output.hexdigest()


def scientific_core(path):
    """Everything except the two explicitly reviewed orchestration functions.

    Analyzer, existing scientific utilities, imports, constants and module
    statements must therefore remain AST-identical. No whitespace exception
    permits a changed scientific expression or an added global binding.
    """
    tree = ast.parse(Path(path).read_text())
    tree.body = [node for node in tree.body if not (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name in IGNORED_ORCHESTRATION)]
    return hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest()


def scientific_imports(imported, source, entrypoint):
    result, source = {}, Path(source).resolve()
    entrypoint = Path(entrypoint).resolve()
    for name, value in imported.items():
        path = Path(name).resolve()
        try:
            relative = str(path.relative_to(source))
        except ValueError:
            continue
        if relative == METADATA or path == entrypoint and relative != HELPER:
            continue
        result[relative] = value
    require(HELPER in result, 'historical/current helper import was not hash-bound')
    return result


def runtime_identity(pin):
    require(pin.get('success') and pin.get('clean'), 'failed/dirty producing runtime pin')
    return {key:canonical(pin[key]) for key in RUNTIME_KEYS}


def seal(directory, expected_commit=None):
    directory = Path(directory).resolve()
    target = directory/'reuse-manifest.json'
    require(not target.exists(), 'refusing to overwrite an existing reuse manifest')
    receipt_path, result_path = directory/'receipt.json', directory/'validation.json'
    receipt, result = (json.loads(path.read_text()) for path in (receipt_path, result_path))
    require(receipt.get('success') and not receipt.get('failures') and result.get('analysis_success')
        and result.get('success') and not result.get('failures'), 'auxiliary analysis was not successfully completed')
    require(receipt.get('native_calls') == 0 and result.get('native_calls') == 0 and
        result.get('only_one_completed_native_record_analyzed') is True, 'not a read-only single-record analysis')
    require(receipt['source_clean_start'] and receipt['source_clean_end'] and
        receipt['source_commit'] == receipt['source_commit_end'], 'historical source changed/was dirty')
    old_source = option(receipt, 'source')
    require(git(old_source, 'rev-parse', 'HEAD') == receipt['source_commit'] and not git(old_source, 'status', '--porcelain'),
        'historical immutable analysis source no longer matches')
    files = dict(receipt['inputs']); check_files(files)
    for name, value in receipt['imported_source_files_start'].items():
        require(digest(name) == value, 'historical imported source changed '+name)
        files[name] = file_identity(name)
    wrapper = Path(receipt['command_argv'][0]).resolve()
    require(digest(wrapper) == receipt['wrapper_sha256'], 'historical wrapper changed')
    files[str(wrapper)] = file_identity(wrapper)
    prefix = receipt['active_trace_prefix']
    require(prefix.get('prefix_sha256_end') == prefix['prefix_sha256'], 'historical native prefix verification failed')
    require(prefix_digest(prefix['path'], prefix['prefix_bytes']) == prefix['prefix_sha256'], 'native trace prefix changed')
    snapshot = Path(prefix['snapshot']).resolve()
    raw = snapshot.read_bytes()
    require(len(raw) == prefix['selected_bytes'] and hashlib.sha256(raw).hexdigest() == prefix['selected_sha256']
        and raw.endswith(b'\n'), 'selected native record snapshot changed/incomplete')
    record = json.loads(raw)
    require(record['alpha'] == result['selected_alpha'] == result['solved_record']['alpha'] and
        record['time'] == result['physical_stage_time'], 'cached solve metadata does not match its native record')
    preflight = json.loads(option(receipt, 'preflight').read_text())
    require(vector_digest(record['history']) == vector_digest(preflight['initial']) == result['source_bound_h0_float64_LE_sha256'],
        'historical h0 binding differs')
    require(canonical(preflight['fields']) == result['native_field_descriptors_sha256'], 'historical field-descriptor hash differs')
    pin = json.loads(option(receipt, 'runtime-pin').read_text())
    runtime = runtime_identity(pin)
    require(runtime == receipt['producing_runtime_sections'] and pin['source_commit'] == receipt['producing_source_commit'],
        'historical producing runtime sections differ')
    require(digest(receipt['producing_solver']['path']) == receipt['producing_solver']['sha256'], 'producing native solver changed')
    require(set(result['dense_matrix_files']).issubset(files), 'matrix pins were not included in the successful receipt')
    for name, entry in result['dense_matrix_files'].items():
        require(entry['sha256'] == files[name]['sha256'] and entry['bytes'] == files[name]['bytes'], 'matrix result/receipt mismatch')
    for path in (receipt_path, result_path, snapshot):
        files[str(path)] = file_identity(path)
    current_source = Path(__file__).resolve().parents[4]
    current_commit = git(current_source, 'rev-parse', 'HEAD')
    require(not git(current_source, 'status', '--porcelain') and (expected_commit is None or current_commit == expected_commit),
        'seal command is not from the requested clean source')
    imports = scientific_imports(receipt['imported_source_files_start'], old_source, wrapper)
    manifest = dict(schema=1, kind='M8-R2 completed auxiliary solve reuse seal',
        sealed_utc=datetime.now(timezone.utc).isoformat(),
        seal_definition='Metadata-only observation after successful completion; original wrapper did not record a validation output hash',
        sealing_source=dict(root=str(current_source), commit=current_commit, tool=str(Path(__file__).resolve()), sha256=digest(__file__)),
        result=str(result_path), receipt=str(receipt_path), files=files,
        old_analysis_source=dict(root=str(old_source), commit=receipt['source_commit'],
            helper=str(old_source/HELPER), helper_sha256=digest(old_source/HELPER),
            scientific_core_ast_sha256=scientific_core(old_source/HELPER), scientific_imports=imports),
        runtime_sections=runtime, producing_source_commit=pin['source_commit'],
        identity=dict(record_index=receipt['record_index'], raw_line_sha256=hashlib.sha256(raw).hexdigest(),
            canonical_record_sha256=canonical(record), stage_xml_sha256=digest(option(receipt, 'stage-session')),
            protected_physics_sha256=result['protected_physics_sha256'],
            native_field_descriptors_sha256=result['native_field_descriptors_sha256'],
            h0_float64_LE_sha256=result['source_bound_h0_float64_LE_sha256'],
            alpha=record['alpha'], time=record['time'], inner_stress_model=record['inner_stress_model'],
            policy=record['policy'], frozen=record['frozen_velocity_pressure_diagnostic']),
        no_gate_reclassification=True, original_failed_gates_must_remain_false=True)
    temporary = target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False)+'\n')
    temporary.replace(target)
    return target


def context(trace, records, raw, session, source, pin, imported, protected_physics_sha256, matrix_root=None):
    identities = []
    with Path(trace).open('rb') as stream:
        for line in stream:
            require(line.strip() and line.endswith(b'\n'), 'final native trace has a blank/incomplete record')
            identities.append(hashlib.sha256(line).hexdigest())
    require(len(identities) == len(records), 'native record/line identity counts differ')
    return dict(records=records, raw_line_sha256=identities, source=Path(source).resolve(),
        session=Path(session).resolve(), stage_xml_sha256=digest(session), runtime_sections=runtime_identity(pin),
        producing_source_commit=pin['source_commit'], fields_sha256=canonical(raw['fields']),
        h0_sha256=vector_digest(raw['initial']), imported=imported,
        protected_physics_sha256=protected_physics_sha256,
        matrix_root=Path(matrix_root or Path(trace).parent).resolve())


def load_reuse(paths, current):
    """Validate every explicitly requested cache before any fresh solve work."""
    output = {}
    for supplied in paths:
        path = Path(supplied).resolve()
        if path.is_dir():
            path = path/'reuse-manifest.json'
        require(path.is_file(), 'missing mandatory reuse manifest '+str(path))
        manifest = json.loads(path.read_text())
        require(manifest.get('schema') == 1 and manifest.get('kind') == 'M8-R2 completed auxiliary solve reuse seal',
            'unsupported reuse manifest')
        check_files(manifest['files'])
        old = manifest['old_analysis_source']
        require(git(old['root'], 'rev-parse', 'HEAD') == old['commit'] and not git(old['root'], 'status', '--porcelain'),
            'old immutable analysis checkout changed')
        require(digest(old['helper']) == old['helper_sha256'], 'old helper bytes changed before AST compatibility check')
        require(scientific_core(old['helper']) == old['scientific_core_ast_sha256'] == scientific_core(current['source']/HELPER),
            'Analyzer/scientific utilities/constants/imports changed')
        sealing = manifest['sealing_source']
        require(digest(sealing['tool']) == sealing['sha256'], 'sealing metadata tool changed')
        require(git(sealing['root'], 'rev-parse', 'HEAD') == sealing['commit'] and not git(sealing['root'], 'status', '--porcelain'),
            'immutable sealing checkout changed')
        current_imports = scientific_imports(current['imported'], current['source'], Path(sys.argv[0]))
        require(set(current_imports) == set(old['scientific_imports']), 'imported scientific module set changed')
        for relative, expected in old['scientific_imports'].items():
            if relative != HELPER:
                require(current_imports[relative] == expected and digest(current['source']/relative) == expected,
                    'imported scientific oracle bytes changed: '+relative)
        identity = manifest['identity']; index = identity['record_index']
        require(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(current['records']) and index not in output,
            'missing/duplicate native record index')
        record = current['records'][index]
        require(current['raw_line_sha256'][index] == identity['raw_line_sha256'] and
            canonical(record) == identity['canonical_record_sha256'], 'selected native raw line/canonical record differs')
        require(current['stage_xml_sha256'] == identity['stage_xml_sha256'] and
            current['protected_physics_sha256'] == identity['protected_physics_sha256'], 'actual XML/protected physics differs')
        require(current['fields_sha256'] == identity['native_field_descriptors_sha256'] and
            current['h0_sha256'] == identity['h0_float64_LE_sha256'] == vector_digest(record['history']),
            'field descriptors or actual h0 differ')
        require(current['runtime_sections'] == manifest['runtime_sections'] and
            current['producing_source_commit'] == manifest['producing_source_commit'], 'producing source/runtime sections differ')
        for name in ('alpha', 'time', 'inner_stress_model', 'policy'):
            require(record[name] == identity[name], 'selected solve '+name+' differs')
        require(record['frozen_velocity_pressure_diagnostic'] == identity['frozen'], 'selected solve frozen flag differs')
        result = json.loads(Path(manifest['result']).read_text())
        require(result['analysis_success'] and result['success'] and not result['failures'], 'sealed result was not successfully analyzed')
        expected_matrices = {}
        for iteration in record['iterations']:
            matrix = Path(iteration['jacobian_binary_file'])
            matrix = (current['matrix_root']/matrix.name).resolve()
            size = len(record['state'])
            require(iteration['jacobian_binary_shape'] == [size, size], 'native current matrix shape differs')
            expected_matrices[str(matrix)] = file_identity(matrix)
        require(set(expected_matrices) == set(result['dense_matrix_files']), 'retained current matrix paths/count differ')
        for name, expected in expected_matrices.items():
            require(expected == {key:result['dense_matrix_files'][name][key] for key in ('bytes', 'sha256')}, 'current matrix bytes differ')
            matrix = result['dense_matrix_files'][name]
            require(matrix['shape'] == [len(record['state']), len(record['state'])] and
                matrix['bytes'] == 8*len(record['state'])**2 and matrix['dtype'] == 'float64' and
                matrix['byteorder'] == sys.byteorder, 'current matrix shape/type/byteorder differs')
        verified_files = dict(manifest['files'])
        verified_files[sealing['tool']] = file_identity(sealing['tool'])
        output[index] = dict(result=result, verified_files=verified_files,
            provenance=dict(manifest=str(path), manifest_file=file_identity(path),
            result=manifest['result'], result_file=file_identity(manifest['result']),
            receipt=manifest['receipt'], receipt_file=file_identity(manifest['receipt']),
            original_analysis_source=old['commit'], original_analysis_root=old['root'], identity=identity,
            scientific_core_ast_sha256=old['scientific_core_ast_sha256'],
            ignored_orchestration_functions=sorted(IGNORED_ORCHESTRATION),
            runtime_sections=manifest['runtime_sections'], seal_definition=manifest['seal_definition'],
            original_checks=result['checks'], original_errors=result['errors'], gate_values_reused_without_reclassification=True))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seal', type=Path, required=True)
    parser.add_argument('--expected-commit')
    args = parser.parse_args()
    target = seal(args.seal, args.expected_commit)
    print(json.dumps(dict(success=True, manifest=str(target), sha256=digest(target), native_calls=0, scientific_calls=0)))


if __name__ == '__main__':
    main()
