#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only CoupledJFNKV1 MMS history/retry/recovery FOUNDATION audit.

No solver is launched. Manifest schema WJH27ImplicitHistoryFoundationV1:
runs maps control,prefix,cold,retry,recovery3 to completed .result.json paths;
input_checkpoints maps cold,recovery3 to {path,sha256}; negatives contains
{name,category,exit,log,expected_marker,session,failure?,checkpoint?}.
root defaults to the manifest directory. Relative paths resolve under root.
Failure/checkpoint are required for iteration_bound and deadline. A checkpoint
may be a path or {path,sha256}; positive input checkpoints require the SHA.
All output is create-only. Native PDE residual records are checked, not solved
again; cache equality/lineage is not an independent advection-equation audit.
"""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
import history_recovery_check as h

START = dict(control=0, prefix=0, cold=1, retry=0, recovery3=3)
END = 4
DT = .02
HISTORY_NAMES = ('current', 'previous', 'advection', 'previous_advection')
TIMERS = ('constitutive', 'wall', 'assembly', 'flow', 'scalar', 'rebuild',
          'guards', 'forcing', 'residual', 'nonlinear', 'io', 'other')
GUARDS = ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k', 'minimum_principal_over_mu')
require = h.require


def read(path):
    return h.native_json(Path(path).read_text())


def suffix(path, old, new):
    require(path.name.endswith(old), 'unexpected native artifact suffix: '+str(path))
    return path.with_name(path.name[:-len(old)]+new)


def timing(record):
    seconds = record['seconds']; categories = record['disjoint_seconds']
    require(math.isfinite(seconds) and seconds > 0 and set(categories) == set(TIMERS), 'missing attempt/step timing inventory')
    reserve = 128*sys.float_info.epsilon*seconds
    require(all(math.isfinite(v) and v >= -reserve for v in categories.values()), 'nonfinite/negative disjoint timing')
    total = math.fsum(categories.values())
    require(abs(total-seconds) <= reserve, 'disjoint timing sum differs from actual attempt/step elapsed time')
    return dict(seconds=seconds, categories=categories, sum=total,
                attributed_fraction=1-categories['other']/seconds,
                attribution_gate='NOT_EVALUATED: tiny history fixture, not the warm-cost gate')


def integer(value, label):
    require(type(value) is int and value >= 0, 'invalid nonnegative count: '+label)
    return value


def configuration(record, control, negative=None):
    if negative is None:
        h.configuration(record, control)
        require(record['coupled_newton']['step_seconds_cap'] == 5., 'positive step cap differs from prospective 5 seconds')
    else:
        # A deliberate numerical-policy mutation has a different signature.
        # Never require that changed signature to equal the positive control.
        keys = ('schema', 'model', 'source_identity', 'reference_identity', 'sizes', 'gauge', 'basis',
                'manufactured', 'quadratureX', 'quadratureR', 'linear_tolerance', 'maximum_iterations',
                'anderson_depth', 'stress_split', 'constitutive_backend', 'mms_variant', 'mms_body_amplitude', 'dt')
        for key in keys: require(record[key] == control[key], 'negative control changed physical/source inventory: '+key)
        require(record['signature'] != control['signature'], 'intentional negative options did not change signature')
        require(record['nonlinear_tolerance'] == (1e-14 if negative == 'iteration_bound' else 3e-12),
                'negative tolerance differs from declared mutation')
        require(record['coupled_newton']['step_seconds_cap'] == (1e-9 if negative == 'deadline' else 5.),
                'negative step cap differs from declared mutation')
    require(record['implicit_solver'] == 'CoupledJFNKV1', 'history fixture is not the opt-in Newton path')
    for key in ('source_identity', 'reference_identity', 'signature'):
        require(len(record[key]) == 64 and all(c in '0123456789abcdef' for c in record[key]), 'invalid SHA256 '+key)
    require(record['steady_velocity_preconditioner'] == control['steady_velocity_preconditioner'], 'preconditioner policy differs')


def accepted(record, control, retry=False):
    configuration(record, control); h.guards(record)
    stats = record['coupled_newton']
    for key in ('krylov_iterations', 'jacobian_actions', 'residual_evaluations', 'preconditioner_applications', 'rejected_trials'):
        integer(stats[key], key)
    require(type(stats['counts_complete']) is bool and stats['counts_complete'] is (not retry),
            'interrupted retry statistics must be explicitly incomplete; other accepted steps complete')
    require(integer(record['retries'], 'retries') == int(retry), 'unexpected accepted retry count')
    require(integer(record['rebuilds'], 'rebuilds') >= (6 if retry else 3), 'missing current block-factor rebuilds')
    require(stats['residual_evaluations'] > 0 and stats['preconditioner_applications'] > 0,
            'accepted Newton step lacks live residual/preconditioner callback work')
    if retry:
        require('qualification injected recoverable Krylov failure' in record['retry_reason'], 'missing injected retry reason')
    else:
        require(not record['retry_reason'], 'unexpected retry reason on clean accepted step')
    return dict(native_reported_residuals=record['residuals'], guards={key: record[key] for key in GUARDS},
                newton=stats, retries=record['retries'], rebuilds=record['rebuilds'], timing=timing(record))


def checkpoint(spec, root, reference, deadline):
    path = h.resolve(root, spec if isinstance(spec, str) else spec['path'])
    expected = h.identity(path)['sha256'] if isinstance(spec, str) else spec['sha256']
    return h.read_checkpoint(path, expected, reference, deadline)


def bind_checkpoint(states, native, hashes):
    for key, header in zip(HISTORY_NAMES, h.HASH_KEYS):
        require(h.state_hash(states[key]) == native[header] == hashes[key], 'checkpoint/native history mismatch: '+key)


def session(path):
    root = ET.parse(path).getroot(); parameters = {}; solver = {}
    for item in root.findall('.//PARAMETERS/P'):
        text = ''.join(item.itertext()).strip(); key, value = text.split('=', 1)
        require(key.strip() not in parameters, 'duplicate session parameter')
        parameters[key.strip()] = value.strip()
    for item in root.findall('.//SOLVERINFO/I'):
        key = item.attrib['PROPERTY']; require(key not in solver, 'duplicate session solver option')
        solver[key] = item.attrib['VALUE']
    return parameters, solver


def failed_initial(record):
    require(record['step'] == 0 and record['time'] == 0 and record['dt'] == DT, 'failed step advanced the accepted clock')
    sizes = record['sizes']; before = record['input_history']
    require(before['step'] == 0 and before['time'] == 0 and before['dt'] == DT, 'failed attempt input history is not fresh step0')
    states = dict(current=h.state(record['coefficients'], sizes), previous=h.state(before['previous'], sizes),
                  advection=h.state(record['current_advection'], sizes),
                  previous_advection=h.state(record['previous_advection'], sizes))
    require(h.state_hash(states['current']) == h.state_hash(h.state(before['current'], sizes)) == h.state_hash(states['previous']),
            'failed attempt changed fresh physical history')
    require(h.state_hash(states['advection']) == h.state_hash(h.state(before['advection'], sizes)) ==
            h.state_hash(h.state(before['previous_advection'], sizes)) == h.state_hash(states['previous_advection']),
            'failed attempt changed fresh explicit history')
    for key, header in zip(HISTORY_NAMES, h.HASH_KEYS):
        require(h.state_hash(states[key]) == record[header], 'failed accepted state history hash differs')
    for key in ('advection', 'previous_advection'):
        require(all(v == 0 for v in states[key]['flow'][2*sizes[0]:]), 'failed explicit cache contains pressure loads')
    return states


def run(manifest, folder, table_path, deadline, out):
    root = h.resolve(folder, manifest.get('root', '.'))
    require(set(manifest['runs']) == set(START), 'wrong implicit-history run inventory')
    paths = {name: h.resolve(root, value) for name, value in manifest['runs'].items()}
    packets = {name: read(path) for name, path in paths.items()}
    control = packets['control']['state']; table = read(table_path)
    require(table['status'] == 'PASS' and h.identity(table_path)['sha256'] == control['reference_identity'], 'reference identity differs')
    Kin = table['parameters']['K']; Ain = Kin**1.5/table['constants']['c_length']
    scales = dict(U=1., Pi=Kin/table['constants']['q_e'], k=Kin, epsilon=Ain/.5)
    require(all(math.isfinite(v) and v > 0 for v in scales.values()), 'invalid fixed reference scales')
    if 'source_identity' in manifest: require(manifest['source_identity'] == control['source_identity'], 'manifest source identity differs')
    out.update(root=str(root), source_identity=control['source_identity'], signature=control['signature'], scales=scales,
               runs={}, input_checkpoints={}, comparisons={}, negatives=[])
    histories = {}; logs = {}; checks = {}; checkpoint_states = {}
    require(set(manifest['input_checkpoints']) == {'cold', 'recovery3'}, 'missing pre-launch input checkpoint identities')
    for name, spec in manifest['input_checkpoints'].items():
        require(isinstance(spec, dict) and set(spec) >= {'path', 'sha256'}, 'input checkpoint must pin path and SHA256')
        states, evidence = checkpoint(spec, root, control, deadline)
        require(evidence['step'] == START[name], 'input checkpoint starts at wrong step')
        checks[name] = evidence; checkpoint_states[name] = states; out['input_checkpoints'][name] = evidence
    for name, packet in packets.items():
        h.clock(deadline); path = paths[name]; step = 1 if name == 'prefix' else END
        require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and
                packet['factory'] == 'WJNativeCoupledV2' and packet['steady_mode'] is False, 'not a completed advancing run')
        require(packet['restart'] is (name in ('cold', 'recovery3')), 'restart label differs')
        record = packet['state']; h.inventory(record, control); accepted(record, control)
        require(record['step'] == step, 'wrong accepted endpoint'); h.same_time(record['time'], step*DT, name)
        require(len(packet['complete_step_seconds']) == step-START[name], 'incomplete native call timing inventory')
        require(all(math.isfinite(v) and v > 0 for v in packet['complete_step_seconds']), 'invalid complete step timing')
        histories[name], labels = h.exported_histories(record)
        log_path = suffix(path, '.result.json', '.steps.jsonl')
        rows = [h.native_json(line) for line in log_path.read_text().splitlines() if line.strip()]
        require([row['step'] for row in rows] == list(range(START[name]+1, step+1)), 'missing/duplicate accepted step log')
        previous = None; diagnostics = []
        for row in rows:
            h.clock(deadline); h.same_time(row['time'], row['step']*DT, 'step log')
            diagnostics.append(accepted(row, control, name == 'retry' and row['step'] == 2))
            if previous is not None:
                require(row['previous_sha256'] == previous['state_sha256'] and
                        row['previous_advection_sha256'] == previous['advection_sha256'], 'accepted history/cache chain differs')
            elif name in checks:
                require(row['previous_sha256'] == checks[name]['history_hashes']['current'] and
                        row['previous_advection_sha256'] == checks[name]['history_hashes']['advection'], 'restart did not use declared checkpoint')
            previous = row
        require(all(rows[-1][key] == record[key] for key in h.HASH_KEYS), 'final output differs from accepted log')
        logs[name] = rows
        final_path = suffix(path, '.result.json', '.final.bin')
        final_states, final = checkpoint(str(final_path), root, control, deadline)
        bind_checkpoint(final_states, record, {key: h.state_hash(histories[name][key]) for key in HISTORY_NAMES})
        require(final['step'] == step, 'final checkpoint clock differs')
        out['runs'][name] = dict(result=h.identity(path), log=h.identity(log_path), final_checkpoint=final,
            history_labels=labels, exact_hashes={key: h.state_hash(v) for key, v in histories[name].items()}, accepted_steps=diagnostics)
    for key in HISTORY_NAMES:
        require(h.state_hash(checkpoint_states['cold'][key]) == h.state_hash(histories['prefix'][key]), 'cold input is not prefix final history')
    require(checks['cold']['path'] == str(suffix(paths['prefix'], '.result.json', '.final.bin')), 'cold restart used a different prefix file')
    for name in ('cold', 'retry', 'recovery3'):
        out['comparisons'][name] = {key: h.compare(histories['control'][key], histories[name][key], control['sizes'], scales)
            for key in ('current', 'previous', 'advection', 'previous_advection', 'older', 'older_advection')}
        out['comparisons'][name]['essential_lifts'] = {}
        for key, scale in (('essential_flow_lift', 1.), ('essential_k_lift', scales['k']), ('essential_epsilon_lift', scales['epsilon'])):
            a, b = control[key], packets[name]['state'][key]
            require(len(a) == len(b) and all(math.isfinite(v) for v in a+b), 'essential lift inventory differs')
            error = max((abs(x-y) for x, y in zip(a, b)), default=0.)/scale
            require(error <= h.COMPARISON_LIMIT, 'essential/Gamma lift disagreement')
            out['comparisons'][name]['essential_lifts'][key] = error

    prefix = suffix(paths['control'], '.result.json', '')
    retained = {Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (3, 4)}
    removed = [Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (1, 2)]
    require(not any(path.exists() for path in removed), 'old recovery1/2 was not rotated')
    require({path.resolve() for path in prefix.parent.glob(prefix.name+'.recovery-*.bin')} == retained ==
            {h.resolve(root, path) for path in packets['control']['retained_recovery']}, 'retained recovery inventory differs')
    require(len(packets['control']['recovery_seconds']) == 4 and packets['control']['rotated_recovery_count'] == 2,
            'control must report four recovery writes and two rotations')
    require(not list(prefix.parent.glob(prefix.name+'.recovery-*.bin.tmp')), 'unfinished recovery temporary file')
    out['rotation'] = dict(removed_absent=[str(p) for p in removed], retained=[h.identity(p) for p in sorted(retained)], science=[])
    science = {}; science_paths = {Path(str(prefix)+'.science-'+str(step)+'.json') for step in range(1, 5)}
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.science-*.json')} == science_paths, 'science1..4 retention differs')
    for step in range(1, 5):
        h.clock(deadline); path = Path(str(prefix)+'.science-'+str(step)+'.json'); record = read(path)
        h.inventory(record, control); accepted(record, control); require(record['step'] == step, 'science clock differs')
        values, labels = h.exported_histories(record); science[step] = values
        require(all(record[key] == logs['control'][step-1][key] for key in h.HASH_KEYS), 'science/log history hashes differ')
        item = dict(file=h.identity(path), history_labels=labels, hashes={key: h.state_hash(v) for key, v in values.items()})
        if step in (3, 4):
            states, evidence = checkpoint(str(Path(str(prefix)+'.recovery-'+str(step)+'.bin')), root, control, deadline)
            bind_checkpoint(states, record, {key: h.state_hash(values[key]) for key in HISTORY_NAMES})
            item['checkpoint'] = evidence
            if step == 3: require(evidence == checks['recovery3'], 'recovery3 import identity differs from retained control checkpoint')
        out['rotation']['science'].append(item)
    out['comparisons']['prefix_at_step1'] = {key: h.compare(science[1][key], histories['prefix'][key], control['sizes'], scales)
        for key in ('current', 'previous', 'advection', 'previous_advection', 'older', 'older_advection')}
    initial = dict(current=science[1]['previous'], previous=science[1]['older'],
                   advection=science[1]['previous_advection'], previous_advection=science[1]['older_advection'])

    negative = manifest['negatives']
    require(len(negative) == 3 and {v['category'] for v in negative} == {'iteration_bound', 'deadline', 'algorithm_identity'}, 'wrong negative-control inventory')
    for spec in negative:
        h.clock(deadline); category = spec['category']; log_path = h.resolve(root, spec['log']); text = log_path.read_text(errors='replace')
        require(type(spec['exit']) is int and spec['exit'] != 0 and spec['expected_marker'] in text, 'negative did not fail with expected log marker')
        session_path = h.resolve(root, spec['session']); parameters, solver = session(session_path)
        require(float(parameters['TimeStep']) == DT and solver['WJNativeMode'] == 'Manufactured' and
                solver['WJMMSVariant'] == 'PolynomialScalarsV1' and float(parameters['WJMMSBodyAmplitude']) == 1.,
                'negative session changed physical fixture')
        require(h.identity(h.resolve(root, solver['WJReferenceTable']))['sha256'] == control['reference_identity'],
                'negative session used a different reference table')
        item = dict(name=spec['name'], category=category, exit=spec['exit'], marker=spec['expected_marker'],
                    log=h.identity(log_path), session=h.identity(session_path))
        if category == 'algorithm_identity':
            require(solver['WJImplicitSolver'] == 'PicardV1' and 'scientific identity differs' in text,
                    'algorithm-identity control did not reject Picard reading Newton history')
            require(h.resolve(root, solver['WJRestart']) == Path(checks['cold']['path']), 'algorithm control used wrong Newton checkpoint')
            item['input_checkpoint'] = checks['cold']; item['initialisation_rejection_only'] = True
        else:
            require(solver['WJImplicitSolver'] == 'CoupledJFNKV1', 'negative is not Newton')
            if category == 'iteration_bound':
                require(float(parameters['WJStepNewtonIterations']) == 1 and float(parameters['WJCoupledTolerance']) == 1e-14,
                        'iteration-bound session does not pin one Newton iteration and 1e-14')
            else:
                require(float(parameters['WJStepNewtonSeconds']) == 1e-9, 'deadline session does not pin 1e-9')
            failure_path = h.resolve(root, spec['failure']); packet = read(failure_path); record = packet['state']
            require(packet['schema'] == 'WJH27CoupledFailureV1' and packet['status'] == 'FAILED_WITH_ACCEPTED_HISTORY_PRESERVED', 'missing failed accepted-history journal')
            require(spec['expected_marker'] in packet['reason'], 'failure journal/log reason differs')
            require(('NEWTON_NOT_CONVERGED' if category == 'iteration_bound' else 'cooperative deadline') in packet['reason'], 'wrong numerical failure category')
            configuration(record, control, category); h.inventory(record, control)
            states = failed_initial(record)
            for key in GUARDS: require(math.isfinite(record[key]) and record[key] > 0, 'failed accepted seed is scientifically inadmissible')
            compare = {key: h.compare(initial[key], states[key], control['sizes'], scales) for key in HISTORY_NAMES}
            binary, evidence = checkpoint(spec['checkpoint'], root, record, deadline)
            require(evidence['step'] == 0, 'failed checkpoint advanced the accepted clock')
            bind_checkpoint(binary, record, {key: h.state_hash(states[key]) for key in HISTORY_NAMES})
            attempt = packet['failed_attempt']; times = timing(attempt)
            require(type(attempt['coupled_statistics_complete']) is bool, 'missing failed statistics completeness flag')
            require(integer(attempt['retries'], 'failed retries') == 0, 'nonlinear/deadline failure triggered a retry')
            for key in ('nonlinear_iterations', 'coupled_krylov_iterations', 'coupled_residual_evaluations', 'coupled_preconditioner_applications', 'rebuilds'):
                integer(attempt[key], key)
            if category == 'iteration_bound':
                require(attempt['nonlinear_iterations'] <= 1 and attempt['coupled_statistics_complete'] is True,
                        'iteration-bound failure lacks complete returned statistics')
                require(attempt['rebuilds'] >= 3 and attempt['coupled_residual_evaluations'] > 0 and
                        attempt['coupled_preconditioner_applications'] > 0, 'iteration-bound failure lacks actual Newton work')
            require(math.isfinite(packet['elapsed_solve_seconds']) and packet['elapsed_solve_seconds'] >= attempt['seconds'], 'failed-attempt wall clock exceeds whole failure elapsed')
            stem = suffix(failure_path, '.failure.json', '')
            require(h.resolve(root, solver['WJOutputPrefix']) == stem, 'negative session/failure prefix differs')
            require(not Path(str(stem)+'.result.json').exists() and not Path(str(stem)+'.final.bin').exists(), 'failed run also exported a successful final artifact')
            log_steps = Path(str(stem)+'.steps.jsonl')
            require(log_steps.is_file() and not log_steps.read_text().strip(), 'failed first step wrote accepted step records')
            item.update(failure=h.identity(failure_path), checkpoint=evidence, source_identity=record['source_identity'],
                altered_signature=record['signature'], accepted_seed_comparisons=compare, failed_attempt=attempt,
                timing=times, accepted_residual_kind=packet['accepted_state_residual_kind'])
        out['negatives'].append(item)
    out.update(status='PASS_IMPLICIT_HISTORY_FOUNDATION_ONLY', foundation_gate_pass=True,
               cache_scope='Complete explicit vectors, binary hashes, shifts, duplicate BE seeds and restart equivalence; no independent advection recomputation here.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('manifest', 'table', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True); args = parser.parse_args()
    began = time.monotonic(); out = dict(schema='WJH27ImplicitHistoryFoundationCheckV1', status='FAIL', foundation_gate_pass=False,
        full_coupled_gate_pass=False, H27_D_pass=False, physical_cold_restart_qualified=False,
        terminal_handling_qualified=False, production_monitoring_storage_qualified=False,
        scope='Tiny PolynomialScalarsV1/body1 implicit-Newton history/retry/rotation FOUNDATION only.',
        script=h.identity(__file__), history_reader=h.identity(h.__file__), manifest=h.identity(args.manifest), table=h.identity(args.table),
        prospective=dict(dt=DT, final_step=END, final_time=.08, coefficient_limit=1e-8, native_residual_limit=3e-12,
            retry_step=2, recovery_writes=[1, 2, 3, 4], retained_recovery=[3, 4]), workers=1)
    failure = None
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little', 'authoritative little-endian Linux analysis required')
        require(math.isfinite(args.seconds) and 0 < args.seconds <= 180, 'analysis cap must be in (0,180] seconds')
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS, (cap, hard)); out['memory_cap_bytes'] = cap
        manifest = read(args.manifest); require(manifest['schema'] == 'WJH27ImplicitHistoryFoundationV1', 'wrong manifest schema')
        run(manifest, args.manifest.resolve().parent, args.table.resolve(), began+args.seconds, out)
    except Exception as error:
        failure = error; out.update(exception=type(error).__name__, reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream: json.dump(out, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: out[key] for key in ('status', 'reason', 'elapsed_seconds') if key in out}))
    if failure is not None: raise SystemExit(2)


if __name__ == '__main__': main()
