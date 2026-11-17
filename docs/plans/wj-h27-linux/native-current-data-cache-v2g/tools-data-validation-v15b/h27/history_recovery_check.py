#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only tiny MMS history/retry/rotation FOUNDATION checker.

Reads WJH27HistoryFoundationManifestV1; launches nothing. Paths resolve under
manifest.root. runs is a dictionary control,prefix,cold,retry,recovery4,recovery6,
whose values are result paths or {result,steps}. The .steps.jsonl path otherwise
follows the .result.json prefix. input_checkpoints maps cold/recovery4/recovery6
to {path,sha256}; control_prefix and negative_controls are required. table
defaults to root/science/reference-tables-v1.json. All outputs are create-only.

This does not qualify physical cold restart, terminal handling, production
monitoring/storage, H27-D, or the full coupled model. Native residual/guard
records are checked against their fixed limits, not recomputed as a PDE oracle.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
import time

DT = .02
FINAL_STEP = 7
COMPARISON_LIMIT = 1e-8
RESIDUAL_LIMIT = 3e-12
STARTS = dict(control=0, prefix=0, cold=1, retry=0, recovery4=4, recovery6=6)
HASH_KEYS = ('state_sha256', 'previous_sha256', 'advection_sha256', 'previous_advection_sha256')
FIELDS = ('flow', 'k', 'epsilon')


def require(condition, message):
    if not condition: raise ValueError(message)


def identity(path):
    data = Path(path).read_bytes()
    return dict(path=str(Path(path).resolve()), sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))


def native_json(text):
    # std::ostream may emit a negative IEEE zero as the integer-looking token
    # -0. Preserve its sign for binary history SHA binding, while retaining
    # integer step/map tokens as integers.
    return json.loads(text, parse_int=lambda token: -0. if token == '-0' else int(token))


def clock(deadline):
    if time.monotonic() >= deadline: raise TimeoutError('history foundation checker reached its admitted cap')


def resolve(root, path):
    path = Path(path)
    return (path if path.is_absolute() else root/path).resolve()


def state(packet, sizes):
    result = {}
    for key, size in zip(FIELDS, (2*sizes[0]+sizes[1], sizes[2], sizes[3])):
        value = packet[key]
        require(isinstance(value, list) and len(value) == size and all(math.isfinite(v) for v in value),
                'invalid native '+key+' history inventory')
        result[key] = list(map(float, value))
    return result


def state_hash(value):
    digest = hashlib.sha256()
    for key in FIELDS:
        for v in value[key]: digest.update(struct.pack('<d', v))
    return digest.hexdigest()


def same_time(value, expected, label):
    require(math.isfinite(value) and abs(value-expected) <= 1e-13, label+' time differs')


def configuration(record, control):
    for key in ('schema', 'model', 'signature', 'source_identity', 'reference_identity', 'sizes', 'gauge', 'basis',
                'quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations',
                'anderson_depth', 'stress_split', 'constitutive_backend', 'mms_variant', 'mms_body_amplitude'):
        require(record[key] == control[key], 'scientific/configuration binding differs: '+key)
    require(record['schema'] == 'WJH27NativeCoupledStateV2' and record['manufactured'] is True and record['gauge'] is True,
            'fixture is not the tiny gauged MMS state')
    require(record['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and record['basis'] == 'Modified_A', 'wrong model/basis')
    require(record['mms_variant'] == 'PolynomialScalarsV1' and record['mms_body_amplitude'] == 1., 'wrong verification fixture')
    require(record['dt'] == DT and record['linear_tolerance'] == 3e-13 and record['nonlinear_tolerance'] == RESIDUAL_LIMIT and
            record['anderson_depth'] == 0 and record['maximum_iterations'] == 60, 'prospective time/iteration policy differs')
    require(len(record['source_identity']) == 64 and len(record['signature']) == 64, 'missing scientific identity')


def inventory(record, control):
    for key in ('geometry', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids'):
        require(record[key] == control[key], 'native geometry/essential map differs: '+key)
    elements = lambda d: sorted(d['elements'], key=lambda e: e['geometry'])
    require(elements(record) == elements(control), 'native signed element maps differ')
    require(len(record['elements']) == 8, 'foundation must use nx2/nr2 (eight cells)')
    geometry = record['geometry']
    for key, expected in dict(left=0., right=1., top=.5, H=.15, d0=.0875).items():
        require(abs(geometry[key]-expected) <= 1e-12, 'tiny MMS geometry differs: '+key)
    axes = ((0., .5, 1.), (0., .075, .15, .325, .5))
    def axis_index(value, axis):
        # Native affine geometry may report the same planned edge a few ulps
        # apart. Classify bounds solely for this inventory check: retain every
        # original geometry owner, signed coefficient map and physical bound.
        tolerance = 64*sys.float_info.epsilon*(axis[-1]-axis[0])
        matches = [i for i, expected in enumerate(axis)
                   if math.isfinite(value) and abs(value-expected) <= tolerance]
        require(len(matches) == 1, 'unexpected tiny mesh edge outside affine-roundoff allowance')
        return matches[0]
    cells = []
    for element in record['elements']:
        bounds = element['bounds']; require(len(bounds) == 4, 'invalid tiny element bounds')
        xl, xr = (axis_index(bounds[i], axes[0]) for i in (0, 1))
        rl, rr = (axis_index(bounds[i], axes[1]) for i in (2, 3))
        require(xr == xl+1 and rr == rl+1, 'tiny element does not span adjacent prospective edges')
        cells.append((xl, rl))
    require(len(set(cells)) == 8 and set(cells) == {(x, r) for x in range(2) for r in range(4)},
            'unexpected tiny mesh cell coverage')


def guards(record):
    residual = record['residuals']
    require(len(residual) == 6 and all(math.isfinite(v) and 0 <= v <= RESIDUAL_LIMIT for v in residual),
            'native fixed-norm current residual exceeds 3e-12')
    for key in ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k', 'minimum_principal_over_mu'):
        require(math.isfinite(record[key]) and record[key] > 0, 'nonpositive/nonfinite accepted-state guard: '+key)


def exported_histories(record):
    sizes = record['sizes']; step = record['step']; h = record['input_history']
    require(isinstance(step, int) and step >= 1 and h['step'] == step-1 and h['dt'] == DT,
            'input history is not the complete pre-final accepted history')
    same_time(record['time'], step*DT, 'accepted')
    same_time(h['time'], (step-1)*DT, 'input history')
    histories = dict(current=state(record['coefficients'], sizes), previous=state(h['current'], sizes),
        advection=state(record['current_advection'], sizes), previous_advection=state(record['previous_advection'], sizes),
        older=state(h['previous'], sizes), older_advection=state(h['previous_advection'], sizes))
    for name, key in zip(('current', 'previous', 'advection', 'previous_advection'), HASH_KEYS):
        require(state_hash(histories[name]) == record[key], 'exported '+name+' hash differs')
    require(state_hash(state(h['advection'], sizes)) == record['previous_advection_sha256'], 'EX2 current history is stale/mislabeled')
    for name in ('advection', 'previous_advection', 'older_advection'):
        require(all(v == 0. for v in histories[name]['flow'][2*sizes[0]:]), 'advection history has a pressure load')
    if step == 1:
        require(state_hash(histories['previous']) == state_hash(histories['older']) and
                state_hash(histories['previous_advection']) == state_hash(histories['older_advection']),
                'BE startup does not retain duplicate step-zero histories')
    labels = dict(current=dict(step=step, time=step*DT), previous=dict(step=step-1, time=(step-1)*DT),
        older=dict(step=max(0, step-2), time=max(0, step-2)*DT),
        advection=dict(step=step, time=step*DT), previous_advection=dict(step=step-1, time=(step-1)*DT),
        older_advection=dict(step=max(0, step-2), time=max(0, step-2)*DT),
        final_equation='BE_STARTUP' if step == 1 else 'BDF2_EX2',
        final_EX_coefficients=[1., 0.] if step == 1 else [2., -1.])
    return histories, labels


def read_checkpoint(path, expected_identity, control, deadline):
    clock(deadline); data = path.read_bytes(); require(64 < len(data) < 128*1024*1024, 'invalid checkpoint size')
    file_identity = identity(path)
    require(file_identity['sha256'] == expected_identity, 'input checkpoint file SHA differs')
    payload, checksum = data[:-64], data[-64:]
    require(hashlib.sha256(payload).hexdigest().encode('ascii') == checksum, 'independent checkpoint checksum mismatch')
    magic, signature, body = payload.split(b'\n', 2)
    require(magic == b'WJ_H27_NATIVE_COUPLED_V2' and signature.decode('ascii') == control['signature'],
            'independent checkpoint scientific identity differs')
    cursor = 0
    def unpack(fmt):
        nonlocal cursor
        length = struct.calcsize(fmt); require(cursor+length <= len(body), 'truncated checkpoint')
        value = struct.unpack_from(fmt, body, cursor); cursor += length; return value
    step, when, dt = unpack('<Qdd'); same_time(when, step*DT, 'checkpoint'); require(dt == DT, 'checkpoint dt differs')
    histories = {}; sizes = control['sizes']
    for name in ('current', 'previous', 'advection', 'previous_advection'):
        item = {}
        for key, size in zip(FIELDS, (2*sizes[0]+sizes[1], sizes[2], sizes[3])):
            length, = unpack('<Q'); require(length == size, 'checkpoint vector size differs')
            item[key] = list(unpack('<'+str(length)+'d'))
            require(all(math.isfinite(v) for v in item[key]), 'nonfinite checkpoint value')
        histories[name] = item
    require(cursor == len(body), 'checkpoint contains trailing bytes')
    return histories, dict(**file_identity, payload_sha256=checksum.decode('ascii'), signature=signature.decode('ascii'),
        step=step, time=when, dt=dt, history_hashes={key: state_hash(v) for key, v in histories.items()})


def compare(first, second, sizes, scales):
    n = sizes[0]
    difference = dict(U=max((abs(a-b) for a, b in zip(first['flow'][:2*n], second['flow'][:2*n])), default=0.),
        Pi=max((abs(a-b) for a, b in zip(first['flow'][2*n:], second['flow'][2*n:])), default=0.)/scales['Pi'],
        k=max((abs(a-b) for a, b in zip(first['k'], second['k'])), default=0.)/scales['k'],
        epsilon=max((abs(a-b) for a, b in zip(first['epsilon'], second['epsilon'])), default=0.)/scales['epsilon'])
    require(max(difference.values()) <= COMPARISON_LIMIT, 'fixed-scale coefficient/history disagreement exceeds 1e-8')
    return dict(coefficient_MAX_reference_scaled=difference, exact_hash_equal=state_hash(first) == state_hash(second),
                first_sha256=state_hash(first), second_sha256=state_hash(second))


def step_log(path, label, final, control, input_checkpoint, deadline):
    rows = [native_json(line) for line in path.read_text().splitlines() if line.strip()]
    expected_end = 1 if label == 'prefix' else FINAL_STEP
    require([row['step'] for row in rows] == list(range(STARTS[label]+1, expected_end+1)),
            'missing/duplicate accepted step log; IO_InfoSteps=1 is required')
    previous = None; retry_steps = []
    for row in rows:
        clock(deadline); configuration(row, control); guards(row); same_time(row['time'], row['step']*DT, 'step log')
        require(type(row['retries']) is int and row['retries'] in (0, 1), 'invalid retry count')
        if row['retries']:
            retry_steps.append(row['step'])
            require('qualification injected recoverable Krylov failure' in row['retry_reason'], 'retry lacks injected-failure marker')
        if previous is not None:
            require(row['previous_sha256'] == previous['state_sha256'] and
                    row['previous_advection_sha256'] == previous['advection_sha256'], 'accepted history chain skipped/changed a state')
        elif input_checkpoint is not None:
            require(row['previous_sha256'] == input_checkpoint['history_hashes']['current'] and
                    row['previous_advection_sha256'] == input_checkpoint['history_hashes']['advection'], 'cold start did not use the declared checkpoint histories')
        previous = row
    require(retry_steps == ([2] if label == 'retry' else []), 'retry must occur exactly once at step2 in retry run only')
    require(all(rows[-1][key] == final[key] for key in HASH_KEYS), 'final output and final accepted log disagree')
    return rows, dict(file=identity(path), accepted_steps=[v['step'] for v in rows], retry_steps=retry_steps,
        maximum_native_reported_residual=max(v for row in rows for v in row['residuals']),
        minimum_native_reported_guards={key: min(row[key] for row in rows) for key in
            ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k', 'minimum_principal_over_mu')})


def run_check(manifest, folder, deadline, out):
    root = resolve(folder, manifest['root']); runs = manifest['runs']; require(set(runs) == set(STARTS), 'run inventory differs')
    paths = {}; packets = {}; histories = {}; labels = {}
    for name, item in runs.items():
        result = item if isinstance(item, str) else item['result']; path = resolve(root, result)
        require(path.name.endswith('.result.json'), 'result filename must retain the native suffix')
        steps = (resolve(root, item['steps']) if isinstance(item, dict) and 'steps' in item else
                 path.with_name(path.name[:-len('.result.json')]+'.steps.jsonl'))
        paths[name] = (path, steps); packets[name] = native_json(path.read_text())
    control = packets['control']['state']; table_path = resolve(root, manifest.get('table', 'science/reference-tables-v1.json'))
    table = json.loads(table_path.read_text()); require(table['status'] == 'PASS' and identity(table_path)['sha256'] == control['reference_identity'], 'reference table identity differs')
    Kin = table['parameters']['K']; Ain = Kin**1.5/table['constants']['c_length']
    scales = dict(U=1., Pi=Kin/table['constants']['q_e'], k=Kin, epsilon=Ain/.5)
    out.update(root=str(root), source_identity=control['source_identity'], reference=identity(table_path), signature=control['signature'], scales=scales)
    if 'source_identity' in manifest: require(manifest['source_identity'] == control['source_identity'], 'manifest source binding differs')
    out['runs'] = {}
    for name, packet in packets.items():
        clock(deadline)
        require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and packet.get('steady_mode', False) is False,
                'not a completed normal-main advancing run')
        require(packet['restart'] is (name in ('cold', 'recovery4', 'recovery6')), 'cold/warm run label differs')
        record = packet['state']; configuration(record, control); inventory(record, control); guards(record)
        expected = 1 if name == 'prefix' else FINAL_STEP
        require(record['step'] == expected, 'wrong completed endpoint step'); same_time(record['time'], expected*DT, name)
        require(len(packet['complete_step_seconds']) == expected-STARTS[name], 'completed-call timing inventory differs')
        histories[name], labels[name] = exported_histories(record)
        out['runs'][name] = dict(result=identity(paths[name][0]), endpoint_step=expected, history_labels=labels[name],
            exact_history_hashes={key: state_hash(v) for key, v in histories[name].items()})
    checkpoints = {}; checkpoint_states = {}; out['input_checkpoints'] = {}
    require(set(manifest['input_checkpoints']) == {'cold', 'recovery4', 'recovery6'}, 'input checkpoint inventory differs')
    for name, spec in manifest['input_checkpoints'].items():
        h, record = read_checkpoint(resolve(root, spec['path']), spec['sha256'], control, deadline)
        require(record['step'] == STARTS[name], 'cold checkpoint step label differs')
        checkpoints[name] = record; checkpoint_states[name] = h; out['input_checkpoints'][name] = record
    for key in ('current', 'previous', 'advection', 'previous_advection'):
        require(state_hash(checkpoint_states['cold'][key]) == state_hash(histories['prefix'][key]), 'prefix final checkpoint is not the exported prefix history')
    logs = {}
    for name in runs:
        logs[name], evidence = step_log(paths[name][1], name, packets[name]['state'], control, checkpoints.get(name), deadline)
        out['runs'][name]['accepted_log'] = evidence
    out['final_comparisons'] = {}
    for name in ('cold', 'retry', 'recovery4', 'recovery6'):
        out['final_comparisons'][name] = {key: compare(histories['control'][key], histories[name][key], control['sizes'], scales)
            for key in ('current', 'previous', 'advection', 'previous_advection', 'older', 'older_advection')}
    # The previously accepted state and its EX2 source remain exactly linked
    # within each run; cross-run differences may be roundoff, never relabeling.
    prefix = resolve(root, manifest['control_prefix'])
    require(paths['control'][0] == Path(str(prefix)+'.result.json'), 'control prefix/result mismatch')
    old = Path(str(prefix)+'.recovery-2.bin'); retained = [Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (4, 6)]
    require(not old.exists() and all(p.is_file() for p in retained), 'actual recovery rotation/retention differs')
    actual_recovery = {p.resolve() for p in prefix.parent.glob(prefix.name+'.recovery-*.bin')}
    reported = {resolve(root, p) for p in packets['control']['retained_recovery']}
    require(actual_recovery == reported == set(retained), 'reported/actual retained recovery inventory differs')
    require(len(packets['control']['recovery_seconds']) == 3 and packets['control']['rotated_recovery_count'] == 1,
            'recovery must report three writes and one rotation')
    require(not list(prefix.parent.glob(prefix.name+'.recovery-*.bin.tmp')), 'unfinished recovery temporary file')
    out['rotation'] = dict(native_reported_write_count=3, native_reported_rotated_count=1, removed_path=str(old),
        removed_path_absent=True, retained=[identity(p) for p in retained], science=[])
    for step in (2, 4, 6):
        path = Path(str(prefix)+'.science-'+str(step)+'.json'); record = native_json(path.read_text())
        configuration(record, control); inventory(record, control); guards(record)
        require(record['step'] == step, 'science step label differs')
        science_histories, science_labels = exported_histories(record)
        log = next(v for v in logs['control'] if v['step'] == step)
        require(all(record[key] == log[key] for key in HASH_KEYS), 'science and accepted-step hashes differ')
        entry = dict(file=identity(path), history_labels=science_labels, hashes={key: record[key] for key in HASH_KEYS})
        if step in (4, 6):
            name = 'recovery'+str(step); require(Path(checkpoints[name]['path']) == retained[(step-4)//2], 'recovery run did not use the retained control file')
            require(checkpoints[name]['sha256'] == identity(path.with_name(prefix.name+'.recovery-'+str(step)+'.bin'))['sha256'], 'retained checkpoint identity changed')
            for key in ('current', 'previous', 'advection', 'previous_advection'):
                require(state_hash(checkpoint_states[name][key]) == state_hash(science_histories[key]), 'recovery checkpoint/science history mismatch')
            entry['checkpoint_history_binding'] = 'EXACT'
        out['rotation']['science'].append(entry)
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.science-*.json')} ==
            {Path(str(prefix)+'.science-'+str(s)+'.json') for s in (2, 4, 6)}, 'science retention inventory differs')
    controls = manifest['negative_controls']; require(len(controls) == 2, 'two restart rejection controls required')
    out['negative_controls'] = []; categories = set()
    for item in controls:
        clock(deadline); marker = item['expected_marker']; text = resolve(root, item['log']).read_text(errors='replace')
        category = ('dt_incompatible' if 'scientific identity differs' in marker else
                    'checksum_corrupted' if 'checksum mismatch' in marker else None)
        require(category is not None and category not in categories, 'negative-control markers must independently exercise dt and checksum rejection')
        categories.add(category)
        require(type(item['exit']) is int and item['exit'] != 0 and marker in text, 'invalid checkpoint was accepted or lacked expected rejection marker')
        out['negative_controls'].append(dict(name=item['name'], category=category, exit=item['exit'],
            expected_marker=marker, marker_found=True, log=identity(resolve(root, item['log']))))
    out.update(status='PASS_MMS_HISTORY_RECOVERY_FOUNDATION_ONLY', foundation_gate_pass=True,
        comparison_norm='MAX over identical signed native coefficient inventories; fixed physical/reference scales, no RHS normalization',
        residual_evidence='Every accepted-step native fixed-normalized residual/guard record checked; this checker does not recompute the nonlinear PDE.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True); args = parser.parse_args(); start = time.monotonic()
    out = dict(schema='WJH27IndependentHistoryFoundationCheckV1', status='FAIL', foundation_gate_pass=False,
        full_coupled_gate_pass=False, H27_D_pass=False, physical_cold_restart_qualified=False,
        terminal_handling_qualified=False, production_monitoring_storage_qualified=False,
        scope='Tiny PolynomialScalarsV1 body1 MMS/history/retry/recovery-rotation FOUNDATION only.',
        script=identity(__file__), manifest=identity(args.manifest),
        prospective=dict(dt=DT, final_step=FINAL_STEP, final_time=.14, coefficient_comparison_limit=COMPARISON_LIMIT,
            native_fixed_residual_limit=RESIDUAL_LIMIT, expected_retry_step=2, recovery_writes=[2, 4, 6], retained_recovery=[4, 6]))
    failure = None
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little', 'authoritative little-endian Linux analysis required')
        require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite positive admitted analysis cap required')
        manifest = json.loads(args.manifest.read_text()); require(manifest['schema'] == 'WJH27HistoryFoundationManifestV1', 'wrong manifest schema')
        run_check(manifest, args.manifest.resolve().parent, start+args.seconds, out)
    except Exception as error:
        failure = error; out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-start; args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as file: file.write(json.dumps(out, indent=2, allow_nan=False)+'\n')
    if failure is not None: raise failure
    print(json.dumps(dict(status=out['status'], foundation_gate_pass=out['foundation_gate_pass'], elapsed_seconds=out['elapsed_seconds'])))


if __name__ == '__main__': main()
