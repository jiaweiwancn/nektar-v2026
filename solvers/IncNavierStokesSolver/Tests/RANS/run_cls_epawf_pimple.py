#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native IMEX1 single-sweep outer6 evolution. Linux only; independent fixed dt.

A numerical completion outside R2LocalV1 is explicitly labelled diagnostic.
The optional exact-flow work is verification-only, never a production sweep.
"""
import argparse
import json
import platform
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import (ROOT, ORIGINAL, H0, THREADS, session, set_info,
    set_parameter, write, run, load, sha, state_hash, write_json, require, readback)


def prepare(output, dt, steps, diagnostic=False, verify=False, restart=None,
            checkpoint=False, reject_step=0):
    output.mkdir(parents=True, exist_ok=False)
    xml = session(ORIGINAL, output/'case.xml', 'outer6')
    root = ET.parse(xml).getroot()
    root.find('./CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text = '1'
    for key, value in {
        'RANSWallTemporalCoupling': 'PimpleSingleSweep',
        'RANSWallNonlinearPolicy': 'Legacy',
        'RANSWallPhysicalRangeGuard': 'R2DiagnosticV1' if diagnostic else 'R2LocalV1',
        'RANSWallSolveDiagnostics': 'Off', 'RANSWallPimpleVerify': 'On' if verify else 'Off',
        'RANSWallRestart': str(restart) if restart else 'None',
        'RANSWallCheckpoint': str(output/'checkpoint') if checkpoint else 'None',
    }.items():
        set_info(root, key, value)
    for key, value in {'NumSteps': steps, 'TimeStep': dt, 'RANSWallEta': 0,
        'RANSWallDiagnostics': 1, 'RANSWallPimpleRejectStep': reject_step,
        'RANSWallPicardFlowRelaxation': 1, 'RANSWallPicardKRelaxation': 1,
        'RANSWallPicardEpsilonRelaxation': 1, 'RANSWallPicardPseudoAlpha': 0}.items():
        set_parameter(root, key, value)
    write(root, xml)
    return xml


def identity(executable):
    linked = subprocess.check_output(['ldd', str(executable)], text=True)
    require('not found' not in linked, 'unresolved shared library')
    libraries = {Path(w).resolve() for w in linked.split() if w.startswith('/')}
    return dict(path=str(executable), sha256=sha(executable), ldd=linked,
                libraries={str(p): sha(p) for p in sorted(libraries)})


def check_report(report):
    require(report['temporal_coupling'] == 'PimpleSingleSweep' and report['order'] == 1,
            'wrong temporal mode')
    require(state_hash(report['initial']) == H0, 'changed represented outer6 initial state')
    require(report['eta'] == 0 and report['inner_stress_model'] == 'SaturatedR2', 'changed physics')
    attempts = report['pimple_attempts']
    require(report['implicit_callback_count'] == len(attempts), 'callback/attempt mismatch')
    for row in attempts:
        if not row['scheduled_complete']:
            continue
        for key in ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
                    'epsilon_transport_solves', 'k_transport_solves'):
            require(row[key] == 1, 'extra or missing governing solve: '+key)
        require(row['mass_inverse_applications'] == row['pressure_gauge_size'] and
                row['additional_correction_or_restoration_solves'] == 0, 'hidden correction')
        require(max(row['corrected_continuity_gauge']) < 1e-10 and
                row['boundary_defect'] < 1e-13 and row['momentum_splitting_identity_defect'] < 1e-9,
                'flow constraints/identity failed')
        require(row['imposed_trace_defect'] < 1e-12 and row['trace_mass_column_max'] > 0 and
                row['retained_trace_column_max'] > 0, 'epsilon trace rows/columns failed')
        require(row['epsilon_frozen_linear_residual'] < 1e-9 and
                row['k_frozen_linear_residual'] < 1e-9, 'scalar equation residual')
    return attempts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--solver', type=Path)
    p.add_argument('--resolution-probe', type=Path)
    p.add_argument('--dt', type=float, default=.0025)
    p.add_argument('--steps', type=int, default=1)
    p.add_argument('--diagnostic', action='store_true')
    p.add_argument('--verify', action='store_true')
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--checkpoint', action='store_true')
    p.add_argument('--restart', type=Path)
    p.add_argument('--reject-step', type=int, default=0)
    p.add_argument('--timeout', type=int, default=3600)
    a = p.parse_args()
    require(a.dt > 0 and a.steps >= 0 and a.reject_step >= 0, 'invalid time configuration')
    a.output = a.output.resolve()
    if not a.prepare_only:
        require(platform.system() == 'Linux' and a.solver, 'authoritative runs require Linux')
        require(a.output.is_relative_to(Path('/media/jay/新加卷/nektar-v2026-runs')), 'heavy run storage')
    xml = prepare(a.output, a.dt, a.steps, a.diagnostic, a.verify,
                  a.restart.resolve() if a.restart else None, a.checkpoint, a.reject_step)
    git = lambda *args: subprocess.check_output(['git', '-C', str(ROOT), *args], text=True).strip()
    manifest = dict(schema=1, source=git('rev-parse', 'HEAD'), source_status=git('status', '--porcelain'),
        archived_xml_sha256=sha(ORIGINAL), xml_sha256=sha(xml), expected_h0_sha256=H0,
        ranks=1, threads=THREADS, dt=a.dt, requested_steps=a.steps,
        qualification='R2DiagnosticV1; outside-box states are not physically qualified M8' if a.diagnostic else 'R2LocalV1 strict',
        status='prepared')
    if not a.prepare_only:
        a.solver = a.solver.resolve(strict=True)
        manifest['executable'] = identity(a.solver)
        manifest['execution'] = run([str(a.solver), 'case.xml'], a.output, 'solver.log', a.timeout)
        report = load(a.output/'case.evolution.json')
        attempts = check_report(report)
        manifest.update(status='numerical_trajectory_completed' if manifest['execution']['exit_code'] == 0 else 'attempt_rejected',
            accepted_steps=report['steps'], final_time=report['time'], failure=report['failure'],
            callback_count=report['implicit_callback_count'], attempts=attempts,
            fully_converged_M8_stage=False, qualification_events=report.get('physical_range_events', []))
        # Large states/audits stay in the native report, rather than duplicating them.
        manifest['qualification_events'] = [{k:v for k,v in e.items() if k not in ('state','history','audit')}
            for e in manifest['qualification_events']]
        if a.resolution_probe and report['steps'] > 0:
            a.resolution_probe = a.resolution_probe.resolve(strict=True)
            manifest['readback_executable'] = identity(a.resolution_probe)
            norms = readback(a, a.output, report['state'], report['last_stage_history'],
                             report['last_stage_time'], report['last_stage_alpha'], manifest)
            require(max(abs(x-y) for x,y in zip(norms, report['stage_audits'][-1][6:])) < 1e-8,
                    'independent original residual differs')
        manifest['runtime_identity_unchanged'] = identity(a.solver) == manifest['executable']
        require(manifest['runtime_identity_unchanged'], 'runtime binary/library changed')
    write_json(a.output/'manifest.json', manifest)
    print(json.dumps({k:manifest[k] for k in ('status','accepted_steps','final_time','failure') if k in manifest}))
    # A recorded numerical rejection is a result; structural/check failures throw.
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
