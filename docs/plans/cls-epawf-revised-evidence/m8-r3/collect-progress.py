#!/usr/bin/env python3
"""Collect durable R3 completed states, including a timeout/incomplete Jacobian.

Stdlib metadata/hashes only. Convergence observations are not independent
operator, scalar-sign, constitutive or physical-trajectory qualification.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import struct

H0 = '4ad26ed0167bb1765bd0be2c150928280b69dbd82ace8a27ea5f9c8c7f4d0409'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def identity(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def load(path):
    return json.loads(Path(path).read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def collect(progress, expected_h0=H0, expected_unknowns=2165):
    progress = Path(progress).resolve(); files = {}
    def bind(pin):
        require(isinstance(pin, dict) and Path(pin['path']).resolve().parent == progress.parent and
            Path(pin['path']).name.startswith(progress.name+'.solve1.'), 'Payload escapes its R3 solve')
        current = identity(pin['path']); require(current == pin, 'Durable payload changed: '+pin['path'])
        old = files.setdefault(pin['path'], pin); require(old == pin, 'Conflicting payload identities')
        return pin
    def vector(pin, size):
        bind(pin); require(pin['bytes'] == 8*size, 'Wrong durable vector/constraint shape')
        with Path(pin['path']).open('rb') as stream:
            require(all(math.isfinite(value[0]) for value in struct.iter_unpack('<d', stream.read())), 'Nonfinite durable vector')
    snapshots = {}; jacobians = {}; trials = []; terminal = None; protocol = None
    events = []; prefix = hashlib.sha256(); prefix_bytes = 0; tail = None
    with progress.open('rb') as stream:
        for index, line in enumerate(stream, 1):
            if not line.endswith(b'\n'):
                tail = dict(bytes=len(line), sha256=hashlib.sha256(line).hexdigest()); break
            event = json.loads(line)
            require(event.get('schema') == 1 and event.get('policy') == 'ScalarFeasibleR3:v1' and
                event.get('action_policy') == 'VerifiedR3Action:v1' and event.get('solve_number') == 1 and
                event.get('alpha') == 400. and event.get('time') == .0025 and event.get('sequence') == index,
                'Changed R3 policy/target/history or event ordering')
            kind = event['kind']
            if kind == 'solve_started':
                require(index == 1 and event['history_sha256'] == expected_h0 and not event['history_advanced'],
                    'R3 solve did not retain original h0')
                vector(event['history'], expected_unknowns)
                require(event['history']['sha256'] == expected_h0, 'History payload differs from original h0')
            elif kind == 'protocol':
                require(protocol is None and event['tolerance'] == 2e-12 and
                    0 < event['maximum_Newton_updates'] <= 16 and event['no_auxiliary_alpha_sequence'] is True,
                    'Changed R3 tolerance/update/alpha protocol')
                protocol = event
            elif kind == 'iterate_completed':
                pin = bind(event['snapshot']); value = load(pin['path'])
                require(value['iteration'] == event['iteration'] and value['iteration'] not in snapshots,
                    'Duplicate or mismatched completed iterate')
                snapshots[value['iteration']] = dict(pin=pin, value=value)
            elif kind in ('jacobian_started_incomplete', 'jacobian_completed'):
                iteration = event['iteration']
                require(iteration in snapshots, 'A Jacobian started before its completed state was durable')
                if kind == 'jacobian_started_incomplete':
                    require(iteration not in jacobians, 'Repeated current Jacobian')
                    jacobians[iteration] = dict(status='incomplete', event=event)
                else:
                    require(iteration in jacobians and jacobians[iteration]['status'] == 'incomplete', 'No preceding Jacobian start')
                    matrix = bind(event['matrix']); require(matrix['bytes'] == 8*expected_unknowns**2, 'Wrong complete matrix shape')
                    jacobians[iteration] = dict(status='complete', event=event)
            elif kind == 'line_search_trial':
                vector(event['state'], expected_unknowns)
                if 'residual' in event: vector(event['residual'], expected_unknowns)
                trials.append(event)
            elif kind == 'solve_terminal':
                require(terminal is None, 'Repeated terminal R3 record'); terminal = event
            events.append(event); prefix.update(line); prefix_bytes += len(line)
    require(protocol is not None and snapshots, 'No complete R3 protocol/iterate snapshot')
    last_path = Path(str(progress)+'.last.json')
    last = load(last_path); last_index = last['iteration']
    expected_last = Path(str(progress)+'.solve1.iterate'+str(last_index)+'.snapshot.json')
    require(last == load(expected_last), 'Atomic last snapshot differs from its immutable completed snapshot')
    last_pin = bind(identity(expected_last))
    last_precedes_event = last_index not in snapshots
    if last_precedes_event:
        require(last_index == max(snapshots)+1, 'Unpublished last snapshot skipped an internal iterate')
        snapshots[last_index] = dict(pin=last_pin, value=last)
    require(sorted(snapshots) == list(range(max(snapshots)+1)) and last_index == max(snapshots),
        'Completed state sequence or atomic last state is inconsistent')
    compact = []
    for iteration, item in sorted(snapshots.items()):
        value = item['value']
        require(value['schema'] == 1 and value['policy'] == 'ScalarFeasibleR3:v1' and
            value['action_policy'] == 'VerifiedR3Action:v1' and value['solve_number'] == 1 and
            value['alpha'] == 400. and value['time'] == .0025 and value['is_complete_internal_iterate'] is True and
            value['is_physical_time_sample'] is False and value['history_advanced'] is False,
            'Completed snapshot changed its equation/history scope')
        vector(value['state'], expected_unknowns); vector(value['residual'], expected_unknowns)
        vector(value['history'], expected_unknowns)
        require(value['history']['sha256'] == expected_h0, 'Completed snapshot changed h0')
        vector(value['constraints'], 2*value['constraint_pairs'])
        require(value['scalar_domain']['admissible'] is True, 'Completed iterate was outside its scalar evaluation domain')
        compact.append(dict(snapshot=item['pin'], **value))
    # A status marker can be fully published immediately before its matching
    # compact event. Preserve it, without claiming a completed factorization.
    statuses = []
    for path in sorted(progress.parent.glob(progress.name+'.solve1.J*.status.json')):
        value = load(path); iteration = value['iteration']
        require(iteration in snapshots and value['status'] in ('incomplete', 'complete') and
            value['action_policy'] == 'VerifiedR3Action:v1', 'Invalid durable Jacobian marker')
        if value['status'] == 'complete':
            matrix = bind(value['matrix']); require(matrix['bytes'] == 8*expected_unknowns**2, 'Wrong completed matrix marker shape')
        if iteration in jacobians and jacobians[iteration]['status'] == 'complete':
            event = jacobians[iteration]['event']
            require(value['status'] == 'complete' and value['matrix'] == event['matrix'] and value['details'] == event['details'],
                'Jacobian status marker differs from its published completion event')
        statuses.append(dict(file=identity(path), **value))
    if terminal and terminal.get('converged') is True:
        require(max(last['scaled_block_norms']) < protocol['tolerance'], 'Native convergence flag exceeds unchanged residual tolerance')
    benchmark = bool(protocol['benchmark_only'])
    factor_events = [event for event in events if event['kind'] == 'factorization_completed']
    benchmark_complete = bool(benchmark and terminal and terminal.get('benchmark_completed') is True and
        len(factor_events) == 1 and len(jacobians) == 1 and all(row['status'] == 'complete' for row in jacobians.values()))
    for pin in files.values(): require(identity(pin['path']) == pin, 'Durable payload changed while collecting')
    frozen = hashlib.sha256()
    with progress.open('rb') as stream:
        remaining = prefix_bytes
        while remaining:
            block = stream.read(min(1024*1024, remaining)); require(bool(block), 'Published progress prefix was truncated')
            frozen.update(block); remaining -= len(block)
    require(frozen.hexdigest() == prefix.hexdigest(), 'Published progress prefix changed during collection')
    return dict(schema=1, success=True, method='R3 durable metadata/binary collection; no scientific reconstruction',
        native_calls=0, operator_reconstructions=0, progress=dict(path=str(progress),
            complete_prefix_bytes=prefix_bytes, complete_prefix_sha256=prefix.hexdigest(), unpublished_tail=tail),
        protocol=protocol, original_h0_sha256=expected_h0, unknowns=expected_unknowns,
        completed_iterates=compact, completed_iterate_count=len(compact), accepted_update_count=max(snapshots),
        last_complete_snapshot=identity(last_path), last_snapshot=compact[-1],
        atomic_last_snapshot_precedes_completed_event=last_precedes_event,
        jacobians=[dict(iteration=index, **row) for index, row in sorted(jacobians.items())],
        jacobian_status_markers=statuses, factorization_events=factor_events,
        trial_count=len(trials), trial_classifications=dict(Counter(row['classification'] for row in trials)),
        line_search_trials=trials, terminal=terminal, terminal_record_available=terminal is not None,
        complete_J_benchmark=benchmark_complete, native_converged_observation=terminal.get('converged') if terminal else None,
        current_action_accuracy_independently_assessed=False, physical_trajectory_qualification=False,
        root_absence_proven=False, M8_closed=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--progress', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); require(not args.output.exists(), 'Use a fresh collection output')
    value = collect(args.progress); value['collected_utc'] = datetime.now(timezone.utc).isoformat()
    value['collector'] = identity(__file__)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(success=True, output=identity(args.output), completed_iterates=value['completed_iterate_count'],
        terminal_available=value['terminal_record_available'])))


if __name__ == '__main__':
    main()
