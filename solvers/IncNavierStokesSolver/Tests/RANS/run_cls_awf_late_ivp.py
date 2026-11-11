#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded late-state IMEX1 IVPs; parent dt validation is never bypassed.

Each level imports exactly the same physical state/reference, then forms its
own native history. Subsequent blocks use ordinary SAME-dt restart. The first
same-dt step is compared against an independent ordinary parent continuation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
import json

from continue_cls_epawf_pimple import CampaignLock, atomic, process_identity
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_accepted_awf import case_xml, check_report
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, set_info, write


def prepare(folder, steps, dt, parent=None, restart=None, refinement='none', bad_hash=False,inner_k=None,hard_reject=0):
    xml = case_xml('Accepted-AWF-05G', folder, steps, dt, restart=restart,
        companion='regularized', companion_reason='same retained physical state diagnostic IVP; faithful covariance blocker retained',
        inner_refinement=refinement, checkpoint_every=1,hard_reject=hard_reject)
    tree = ET.parse(xml).getroot()
    set_info(tree, 'RANSWallEvaluationReuse', 'ExactStepV1')
    if parent:
        require(not restart, 'diagnostic import is not restart')
        q = checkpoint(parent)
        for k, v in dict(RANSWallDiagnosticIVP='LateStateIMEX1V1', RANSWallIVPParent=str(parent),
            RANSWallIVPParentSHA256='0'*64 if bad_hash else sha(parent/'State.bin'),
            RANSWallIVPParentTimeStep=format(q['dt'], '.17g')).items(): set_info(tree, k, v)
        if inner_k is not None:
            set_info(tree,'RANSWallDiagnosticIVP','LateStateInnerKIMEX1V1')
            set_info(tree,'RANSWallIVPInnerKAmplitude',format(inner_k,'.17g'))
    require(inner_k is None or parent is not None,'inner-k seed belongs only to explicit diagnostic IVP initialization')
    write(tree, xml)
    definition=load(folder/'definition.json');definition['xml_sha256']=sha(xml)
    if parent:
        definition.update(fresh_native_initialization=False,initialization_contract='explicit retained-state IMEX1 IVP; no Lift; own native histories',
            diagnostic_ivp=dict(parent_checkpoint=str(parent),parent_sha256=sha(parent/'State.bin'),physical_t0=q['time'],
                parent_dt=q['dt'],target_dt=dt,inner_k_amplitude=inner_k,refinement=refinement))
    else:definition['initialization_contract']='ordinary compatible native restart' if restart else 'fresh analytical-family native initialization'
    atomic(folder/'definition.json',definition)
    return xml


def native(solver, xml, timeout, lock, expected_success=True):
    folder = xml.parent
    env = dict(os.environ, **THREADS, LD_LIBRARY_PATH=str(solver.parent/'lib'), TMPDIR=str(folder))
    start = time.monotonic()
    with (folder/'solver.log').open('w') as log:
        child = subprocess.Popen(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(folder/'usage.txt'), str(solver), str(xml)],
            cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, pass_fds=(lock.stream.fileno(),))
        atomic(folder/'process.json', dict(identity=process_identity(child.pid), command=[str(solver), str(xml)]))
        try: code = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            # The native child inherits the campaign lock. Do not launch a
            # replacement on timeout while any child could still be advancing.
            atomic(folder/'interruption.json', dict(reason='block timeout', pid=child.pid))
            raise RuntimeError('block timeout: inspect recorded native child before any recovery')
    file = folder/'case.evolution.json'
    report = load(file) if file.exists() else None
    receipt = dict(exit_code=code, seconds=time.monotonic()-start, report=str(file),
        report_sha256=sha(file) if file.exists() else None, xml_sha256=sha(xml), log_sha256=sha(folder/'solver.log'),
        usage=(folder/'usage.txt').read_text().strip())
    atomic(folder/'receipt.json', receipt)
    if expected_success:
        require(code == 0 and report is not None and not report.get('failure'), 'native IVP failed: '+str(folder))
        check_report(report, bool(report.get('steps')), load(folder/'definition.json'))
        require(report['evaluation_reuse']['mode'] == 'ExactStepV1' and report['audit_revision'] == 'NativeElementSideV2', 'wrong reuse/audit')
    else:
        require(code != 0, 'incompatible reader request unexpectedly accepted')
    return report, receipt


def packet(report):
    attempts = report['pimple_attempts']
    require(len(attempts) == 1, 'comparison requires one attempted step')
    return load(Path(attempts[0]['k_snapshot']))


def validate(args, parent, lock):
    root = args.output/'import-validation'; root.mkdir()
    q = checkpoint(parent)
    ordinary, ro = native(args.solver, prepare(root/'ordinary', 1, q['dt'], restart=parent), args.timeout, lock)
    imported, ri = native(args.solver, prepare(root/'imported', 1, q['dt'], parent=parent), args.timeout, lock)
    a, b = packet(ordinary), packet(imported)
    fields = ('committed', 'post_flow', 'pre_k', 'raw', 'corrected', 'history', 'fixed_reference')
    checks = {k: state_hash(a[k]) == state_hash(b[k]) for k in fields}
    checks.update(seed_exact=state_hash(imported['accepted_states'][0]) == state_hash(q['coupled']),
        fixed_reference_exact=state_hash(imported['fixed_reference']) == state_hash(q['reference']),
        fixed_defect_exact=state_hash(imported['fixed_defect']) == state_hash(q['defect']),
        first_accepted_exact=state_hash(ordinary['state']) == state_hash(imported['state']),
        histories_exact=[state_hash(v) for v in ordinary['native_history']] == [state_hash(v) for v in imported['native_history']],
        time_labels_exact=state_hash(ordinary['native_history_times']) == state_hash(imported['native_history_times']),
        physical_time_exact=ordinary['time'] == imported['time'],
        initialization_no_lift=imported['initialization_lift_calls_iterations'] == [0, 0],
        startup_bound_noop=imported['initial_k_bound_check']['activations'] == 0,
        counters_distinct=ordinary['steps'] == q['steps']+1 and imported['steps'] == 1)
    require(all(checks.values()), 'same-dt IVP contract differs: '+str(checks))
    rejected = {}
    tests = {'ordinary_changed_dt': prepare(root/'ordinary-changed-dt', 0, q['dt']/2, restart=parent),
        'parent_hash': prepare(root/'bad-parent-hash', 0, q['dt'], parent=parent, bad_hash=True),
        'physical_mesh': prepare(root/'wrong-mesh', 0, q['dt'], parent=parent, refinement='cutoff-v1')}
    corrupted = root/'corrupted-test-copy'; shutil.copytree(parent, corrupted)
    data = bytearray((corrupted/'State.bin').read_bytes()); data[-1] ^= 1; (corrupted/'State.bin').write_bytes(data)
    tests['checksum'] = prepare(root/'bad-checksum', 0, q['dt'], parent=corrupted)
    for name, xml in tests.items():
        _, rc = native(args.solver, xml, args.timeout, lock, False)
        rejected[name] = dict(receipt=rc, reason=(xml.parent/'solver.log').read_text()[-2000:])
    result = dict(status='passed', checks=checks, ordinary_receipt=ro, imported_receipt=ri,
        parent_sha256=sha(parent/'State.bin'), rejected_contracts=rejected,
        imported_identity=imported['diagnostic_ivp'], ordinary_first_step=ordinary['pimple_attempts'][0],
        imported_first_step=imported['pimple_attempts'][0])
    atomic(root/'summary.json', result)
    return result


def level(label, dt, args, parent, lock):
    root = args.output/label; root.mkdir()
    started = time.monotonic()
    out = dict(status='initializing', dt=dt, parent_sha256=sha(parent/'State.bin'), segments=[], observations={},
        elapsed_time=0., accepted_steps=0, physical_t0=checkpoint(parent)['time'])
    atomic(root/'progress.json', out)
    initial, receipt = native(args.solver, prepare(root/'initial', 0, dt, parent=parent), args.timeout, lock)
    q = checkpoint(parent)
    require(state_hash(initial['state']) == state_hash(q['coupled']) and state_hash(initial['fixed_reference']) == state_hash(q['reference']) and
        state_hash(initial['fixed_defect']) == state_hash(q['defect']), 'IVP seed/reference changed')
    require(initial['initialization_lift_calls_iterations'] == [0, 0] and initial['initial_k_bound_check']['activations'] == 0, 'IVP startup correction')
    cp = Path(initial['checkpoint'])
    out.update(status='running', initialization_receipt=receipt, initialization_identity=initial['diagnostic_ivp'],
        native_seed_sha256=state_hash(initial['state']), initial_history_sha256=[state_hash(v) for v in initial['native_history']],
        initial_history_times_sha256=state_hash(initial['native_history_times']))
    atomic(root/'progress.json', out)
    for elapsed in args.times:
        target = round(elapsed/dt); require(abs(target*dt-elapsed) < 1e-13, 'nonintegral elapsed target')
        while out['accepted_steps'] < target:
            require(time.monotonic()-started < args.max_hours*3600, 'declared IVP wall-cost limit reached')
            start = out['accepted_steps']; end = min(start+args.block_steps, target)
            folder = root/f'steps-{start+1:06d}-{end:06d}'
            report, rc = native(args.solver, prepare(folder, end-start, dt, restart=cp), args.timeout, lock)
            require(state_hash(report['accepted_states'][0]) == state_hash(checkpoint(cp)['coupled']), 'same-dt block changed state')
            cp = Path(report['checkpoint'])
            require(report['steps'] == end and abs(report['time']-(out['physical_t0']+end*dt)) < 5e-11, 'IVP physical clock/count differs')
            out.update(accepted_steps=end, physical_time=report['time'], elapsed_time=end*dt,
                checkpoint=str(cp), checkpoint_sha256=sha(cp/'State.bin'))
            out['segments'].append(dict(start_step=start, end_step=end, receipt=rc, bound_activations=sum(t['bound_activations'] for t in report['pimple_attempts']),
                correction_integral=sum(t['k_bound']['integral_total'] for t in report['pimple_attempts'])))
            if end == target:
                out['observations'][str(elapsed)] = dict(checkpoint=str(cp), xml=str(folder/'case.xml'),
                    report=str(folder/'case.evolution.json.gz'), elapsed_time=end*dt, physical_time=report['time'])
            with (folder/'case.evolution.json').open('rb') as src, gzip.open(folder/'case.evolution.json.gz', 'wb') as dst: shutil.copyfileobj(src, dst)
            (folder/'case.evolution.json').unlink()
            # Hash and discard only accepted unselected dense replay matrices.
            # Every native checkpoint and the complete failed packet survive.
            with (root/'discarded-packets.jsonl').open('a') as tomb:
                for t in report['pimple_attempts'][1:-1]:
                    p = Path(t['k_snapshot']); saved = load(p)
                    tomb.write(json.dumps(dict(packet=str(p), sha256=sha(p), matrix_sha256=saved['matrix_sha256']))+'\n')
                    Path(saved['matrix_file']).unlink(); p.unlink()
            atomic(root/'progress.json', out)
    out.update(status='targets_completed', seconds=time.monotonic()-started)
    atomic(root/'progress.json', out)
    return label, out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('solver', 'parent', 'policy', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--times', nargs='+', type=float, default=[.05, .1]); p.add_argument('--block-steps', type=int, default=32)
    p.add_argument('--timeout', type=float, default=7200); p.add_argument('--max-hours', type=float, default=8)
    p.add_argument('--validate-only', action='store_true')
    a = p.parse_args(); require(platform.system() == 'Linux', 'Linux native validation required')
    require(1 <= a.block_steps <= 32 and a.times == sorted(set(a.times)) and 0 < max(a.times) <= .5, 'bounded chronological IVP targets required')
    a.output = a.output.resolve(); a.output.mkdir(parents=True, exist_ok=False); lock = CampaignLock(a.output)
    parent = a.output/'parent'; shutil.copytree(a.parent, parent)
    require(sha(parent/'State.bin') == sha(a.parent/'State.bin'), 'parent copy differs')
    shutil.copy2(a.policy, a.output/'declared-policy.json')
    out = dict(status='validating_import', started_UTC=datetime.now(timezone.utc).isoformat(), solver_sha256=sha(a.solver),
        parent=str(a.parent), parent_sha256=sha(parent/'State.bin'), script_sha256=sha(Path(__file__)),
        policy_sha256=sha(a.policy), threads=THREADS, jobs=2, cases={})
    atomic(a.output/'summary.json', out)
    out['import_validation'] = validate(a, parent, lock); atomic(a.output/'summary.json', out)
    if not a.validate_only:
        out['status'] = 'running'; atomic(a.output/'summary.json', out)
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(level, n, dt, a, parent, lock) for n, dt in [('coarse', .000625), ('middle', .0003125), ('fine', .00015625)]]
            for label, future in zip(('coarse', 'middle', 'fine'), futures):
                try:
                    name, result = future.result(); out['cases'][name] = result
                except Exception as error:
                    progress = load(a.output/label/'progress.json') if (a.output/label/'progress.json').exists() else {}
                    failures = list((a.output/label).glob('*/case.evolution.json'))
                    rejected = [f for f in failures if load(f).get('failure')]
                    progress.update(status='numerically_rejected' if rejected else 'infrastructure_interrupted',
                        failure=str(error), failure_reports=[str(f) for f in rejected])
                    atomic(a.output/label/'progress.json', progress); out['cases'][label] = progress
                atomic(a.output/'summary.json', out)
    out.update(status='completed' if all(v['status']=='targets_completed' for v in out['cases'].values()) else 'recorded_terminal_outcomes',
        completed_UTC=datetime.now(timezone.utc).isoformat()); atomic(a.output/'summary.json', out)


if __name__ == '__main__': main()
