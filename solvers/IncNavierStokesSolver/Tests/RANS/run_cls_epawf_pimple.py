#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native IMEX1 single-sweep outer6 evolution. Linux only; independent fixed dt.

A numerical completion outside R2LocalV1 is explicitly labelled diagnostic.
The optional exact-flow work is verification-only, never a production sweep.
"""
import argparse
import json
import platform
import struct
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import (ROOT, ORIGINAL, H0, THREADS, session, set_info,
    set_parameter, write, run, load, sha, state_hash, write_json, require, readback)


BOUND_WORK = dict(add='RANSWallPimpleKBoundMaxAdd', total='RANSWallPimpleKBoundMaxTotal',
                  separations='RANSWallPimpleKBoundMaxSeparations', seconds='RANSWallPimpleKBoundMaxSeconds')


def add_bound_work_arguments(parser):
    for key in BOUND_WORK:
        parser.add_argument('--bound-max-'+key, type=float if key=='seconds' else int)


def bound_work_arguments(args):
    return {key:getattr(args,'bound_max_'+key) for key in BOUND_WORK
            if getattr(args,'bound_max_'+key) is not None}


def bound_work_command(controls):
    return [word for key,value in controls.items() for word in ('--bound-max-'+key,str(value))]


def prepare(output, dt, steps, diagnostic=False, verify=False, restart=None,
            checkpoint=False, reject_step=0, k_bound="Off", floor_relative=None, snapshots=False, budget_time=-1., hard_reject_step=0, inventory="LegacyV1", wall_policy="Upstream", bound_work=None):
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
    require(wall_policy in ('Upstream','UniformEPFlat'),'unknown flat wall policy')
    set_info(root, 'RANSWallTransition', wall_policy)
    set_info(root, 'RANSWallPimpleKBound', k_bound)
    set_info(root, 'RANSWallPimpleKInventory', inventory)
    set_info(root, 'RANSWallPimpleSnapshots', 'On' if snapshots else 'Off')
    if floor_relative is not None:set_parameter(root, 'RANSWallPimpleKFloorRelative', floor_relative)
    if budget_time>=0:set_parameter(root, 'RANSWallPimpleBudgetTime', budget_time)
    if hard_reject_step:set_parameter(root,'RANSWallPimpleHardAuditRejectStep',hard_reject_step)
    for key,value in (bound_work or {}).items():
        require(key in BOUND_WORK and value>0,'invalid projection work control')
        set_parameter(root,BOUND_WORK[key],value)
    write(root, xml)
    return xml


def identity(executable):
    linked = subprocess.check_output(['ldd', str(executable)], text=True)
    require('not found' not in linked, 'unresolved shared library')
    libraries = {Path(w).resolve() for w in linked.split() if w.startswith('/')}
    return dict(path=str(executable), sha256=sha(executable), ldd=linked,
                libraries={str(p): sha(p) for p in sorted(libraries)})


def check_report(report, restarting=False):
    require(report['temporal_coupling'] == 'PimpleSingleSweep' and report['order'] == 1,
            'wrong temporal mode')
    if not restarting and report.get('wall_policy','Upstream')=='Upstream':
        require(state_hash(report['accepted_states'][0]) == H0, 'changed represented outer6 initial state')
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
                row.get('additional_governing_equation_solves',0) == 0, 'hidden PDE correction')
        if report.get('k_bound_mode','Off')=='Off':
            require(row['additional_correction_or_restoration_solves']==0 and row.get('bound_calls',0)==0,'Off correction work')
        else:
            require(row['bound_calls']==1 and row['additional_correction_or_restoration_solves']==row['bound_kkt_solves'],'bound accounting')
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
    p.add_argument('--steps', type=int)
    p.add_argument('--final-time', type=float)
    p.add_argument('--k-bound',choices=('Off','MassProjection'),default='Off')
    p.add_argument('--wall-policy',choices=('Upstream','UniformEPFlat'),default='Upstream')
    p.add_argument('--inventory',choices=('LegacyV1','ConstitutiveV2','ConstitutiveV3'),default='LegacyV1')
    p.add_argument('--floor-relative',type=float)
    p.add_argument('--snapshots',action='store_true')
    p.add_argument('--budget-time',type=float,default=-1.)
    p.add_argument('--diagnostic', action='store_true')
    p.add_argument('--verify', action='store_true')
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--checkpoint', action='store_true')
    p.add_argument('--checkpoint-interval', type=int, default=0)
    p.add_argument('--restart', type=Path)
    p.add_argument('--reject-step', type=int, default=0)
    p.add_argument('--hard-reject-step',type=int,default=0)
    p.add_argument('--timeout', type=int, default=3600)
    add_bound_work_arguments(p)
    a = p.parse_args()
    start_time=0.
    if a.restart:
        payload=(a.restart/'State.bin').read_bytes();n=struct.unpack_from('<Q',payload,24)[0]
        _,_,start_time,_=struct.unpack_from('<QQdd',payload,32+n)
    if a.final_time is not None:
        count=round((a.final_time-start_time)/a.dt)
        require(count>=0 and abs(start_time+count*a.dt-a.final_time)<1e-13,'final time must lie on fixed dt grid')
        require(a.steps is None or a.steps==count,'steps conflict with fixed final time');a.steps=count
    if a.steps is None:a.steps=1
    require(a.dt > 0 and a.steps >= 0 and a.reject_step >= 0, 'invalid time configuration')
    a.output = a.output.resolve()
    if not a.prepare_only:
        require(platform.system() == 'Linux' and a.solver, 'authoritative runs require Linux')
        require(a.output.is_relative_to(Path('/media/jay/新加卷/nektar-v2026-runs')), 'heavy run storage')
    xml = prepare(a.output, a.dt, a.steps, a.diagnostic, a.verify,
                  a.restart.resolve() if a.restart else None, a.checkpoint, a.reject_step,
                  a.k_bound,a.floor_relative,a.snapshots,a.budget_time,a.hard_reject_step,a.inventory,a.wall_policy,bound_work_arguments(a))
    require(a.checkpoint_interval>=0,'negative checkpoint interval')
    if a.checkpoint_interval:
        require(a.checkpoint,'checkpoint interval requires checkpoint output')
        root=ET.parse(xml).getroot();set_parameter(root,'IO_CheckSteps',a.checkpoint_interval);write(root,xml)
    git = lambda *args: subprocess.check_output(['git', '-C', str(ROOT), *args], text=True).strip()
    manifest = dict(schema=1, source=git('rev-parse', 'HEAD'), source_status=git('status', '--porcelain'),
        wall_policy=a.wall_policy, actual_initial_sha256=None, archived_xml_sha256=sha(ORIGINAL), xml_sha256=sha(xml), expected_h0_sha256=H0 if a.wall_policy=='Upstream' else None,
        ranks=1, threads=THREADS, dt=a.dt, requested_steps=a.steps,
        qualification='R2DiagnosticV1; outside-box states are not physically qualified M8' if a.diagnostic else 'R2LocalV1 strict',
        status='prepared', projection_work_controls=bound_work_arguments(a))
    if not a.prepare_only:
        a.solver = a.solver.resolve(strict=True)
        manifest['executable'] = identity(a.solver)
        manifest['execution'] = run([str(a.solver), 'case.xml'], a.output, 'solver.log', a.timeout)
        report = load(a.output/'case.evolution.json')
        attempts = check_report(report, bool(a.restart))
        manifest['actual_initial_sha256']=state_hash(report['accepted_states'][0])
        accepted_times={x[0] for x in report['stage_audits']}
        cumulative=0.
        for x in attempts:
            x['numerically_accepted']=x['time'] in accepted_times
            if x['numerically_accepted']:cumulative+=x.get('k_bound',{}).get('integral_total',0.)
            x['cumulative_bound_integral_current_run']=cumulative
        if manifest['execution']['exit_code']==0:
            require(len(accepted_times)==a.steps and abs(report['time']-(start_time+a.steps*a.dt))<1e-12,'incomplete fixed horizon')
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
        final_identity = identity(a.solver)
        manifest['runtime_identity_unchanged'] = all(final_identity[key] == manifest['executable'][key]
            for key in ('path','sha256','libraries'))
        require(manifest['runtime_identity_unchanged'], 'runtime binary/library changed')
    write_json(a.output/'manifest.json', manifest)
    print(json.dumps({k:manifest[k] for k in ('status','accepted_steps','final_time','failure') if k in manifest}))
    # A recorded numerical rejection is a result; structural/check failures throw.
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
