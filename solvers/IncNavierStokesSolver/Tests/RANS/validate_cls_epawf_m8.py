#!/usr/bin/env python3
"""Staged Linux M8 native qualification, with explicit failed physical gates.

Every invocation freezes its model/observation protocol before launching any
solver.  Failed/unfinished native trajectories retain their reports and cannot
be promoted to a robust or developed terrain result.  A successful accounting
test and a physically qualified candidate are separate receipt fields.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

from cls_epawf_m8_cases import (
    VARIANTS, depletion_session, flow_session, manufactured_session,
    perturbation_session,
)
from cls_epawf_m8_native_oracle import verify as verify_energy
from cls_epawf_m8_mms import flat_session as regularized_mms_session, verify as verify_regularized_mms
from cls_epawf_m8_protocol import VERSION, contract, sha
from validate_cls_epawf_m4 import spaces, state_difference, native_history_difference
import validate_cls_epawf_m7 as m7
from validate_cls_epawf_m7_transition_flow import (
    Audit, failure_observation, paired, time_comparison,
)


SOURCE = Path(__file__).resolve().parents[4]
TEMPLATE = Path(__file__).with_name('EquilibriumWallEvolution.xml')


def parameters(args):
    return dict(variant=args.variant, beta=args.beta, gamma=args.gamma,
                gamma_scale=args.gamma_scale, theta=args.theta,
                delta=args.delta, smooth=args.smooth)


def physical_gate(raw):
    """Finite evaluated constitutive gates; never infer all-state bounds."""
    K0 = raw['parameters'][5]
    minima = {}
    minimum_k = minimum_epsilon = math.inf
    minimum_tangent = math.inf
    for record in raw.get('state_audits', []):
        for name in ('source_volume_inner', 'source_volume_outer', 'physical_sample_K',
                     'physical_sample_outer', 'interface_inner', 'interface_outer'):
            row = record.get(name, [])
            if len(row) != 6 or not row[0]:
                continue
            minima[name] = min(minima.get(name, math.inf), row[3])
            minimum_k = min(minimum_k, row[1])
            if 'outer' in name:
                minimum_epsilon = min(minimum_epsilon, row[2])
                minimum_tangent = min(minimum_tangent, row[4])
    covariance_pass = bool(minima) and min(minima.values()) >= -1e-10*K0
    return dict(covariance_pass=covariance_pass,
                covariance_minima=minima,
                finite_sample_turbulence_positive=minimum_k >= 0 and minimum_epsilon > 0,
                minimum_k=minimum_k if math.isfinite(minimum_k) else None,
                minimum_outer_epsilon=minimum_epsilon if math.isfinite(minimum_epsilon) else None,
                minimum_yx_tangent=minimum_tangent if math.isfinite(minimum_tangent) else None,
                scope='Native finite constitutive evaluation sets; no unsampled polynomial or all-direction/state-envelope theorem')


class Runner:
    def __init__(self, args, output, receipt):
        self.args, self.output, self.receipt = args, output, receipt
        self.build = args.build.resolve()
        self.environment = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH'] = str(self.build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH', '')
        self.raws = {}
        self.audit = Audit()
        receipt['environment'] = dict(THREADS, TMPDIR=self.environment.get('TMPDIR'),
            LD_LIBRARY_PATH=self.environment['LD_LIBRARY_PATH'], cpu_affinity=sorted(os.sched_getaffinity(0)))

    def save(self):
        self.receipt['checks'] = self.audit.checks
        self.receipt['errors'] = self.audit.errors
        (self.output/'validation.json').write_text(json.dumps(self.receipt, indent=2, allow_nan=False)+'\n')

    def command(self, name, argv, cwd, ranks=1, timeout=1800):
        argv = list(map(str, argv))
        if ranks > 1:
            argv = ['mpirun', '--report-bindings', '--bind-to', 'core', '--map-by', 'core',
                    '-np', str(ranks), *argv]
        log = cwd/'run.log'
        start = time.monotonic(); code, error = 127, None
        try:
            with log.open('w') as stream:
                result = subprocess.run(argv, cwd=cwd, env=self.environment,
                    stdout=stream, stderr=subprocess.STDOUT, timeout=timeout, check=False)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code, error = 124, 'Hard command timeout'
        except OSError as exception:
            error = str(exception)
        if error:
            with log.open('a') as stream:
                stream.write('\n'+error+'\n')
        record = dict(name=name, command=argv, cwd=str(cwd), ranks=ranks,
                      numerical_threads=THREADS, exit_code=code,
                      seconds=time.monotonic()-start, timeout_seconds=timeout,
                      log=str(log), log_sha256=sha(log))
        self.receipt['commands'].append(record)
        self.save()
        print(f'{name}: exit {code} ({record["seconds"]:.2f}s)', flush=True)
        return code, record

    def case(self, name, generator=perturbation_session, ranks=1,
             source_mechanism=False, full_energy=True, **kwargs):
        directory = self.output/name; directory.mkdir()
        selected = dict(parameters(self.args), **kwargs)
        xml = generator(TEMPLATE, directory/'case.xml', **selected)
        record = dict(name=name, parameters=selected, requested_ranks=ranks,
                      xml=str(xml), xml_sha256=sha(xml), status='running',
                      source_mechanism_only=source_mechanism)
        self.receipt['cases'].append(record); self.save()
        code, command = self.command(name, [self.build/'IncNavierStokesSolver', xml], directory,
                                     ranks, self.args.timeout)
        path = directory/'case.evolution.json'
        raw = {}
        if path.is_file():
            record.update(native_report=str(path), native_report_sha256=sha(path))
            try:
                raw = json.loads(path.read_text())
            except (ValueError, OSError) as error:
                record['native_report_error'] = str(error)
        record.update(exit_code=code, seconds=command['seconds'])
        if code or not raw or raw.get('failure'):
            record.update(status='stopped', failure_observation=failure_observation(raw, code))
            if raw.get('state_audits'):
                record['accepted_prefix_physical_gates'] = physical_gate(raw)
            self.receipt['failures'].append(dict(case=name, **record['failure_observation']))
            self.save()
            return None
        try:
            self.audit.require(raw['native_time_map'] and raw['ranks'] == ranks, name+': missing native collective time map')
            self.audit.require(raw['model_variant'] == selected['variant'], name+': model selector was not applied')
            self.audit.near(raw['robust_parameters'][:6],
                [selected['beta'], selected['gamma'], selected['theta'], selected['delta'],
                 int(selected['smooth']), selected['gamma_scale']], name+': explicit model parameters', 1e-13)
            self.audit.require(bool(raw['robust_parameters'][6]) == (selected['variant'] != 'Original'),
                               name+': wrong capped source selection')
            record['native_audit'] = m7.audits(raw, name)
            self.audit.checks += 1
            record['independent_energy'] = verify_energy(raw, full_volume=full_energy)
            self.audit.checks += record['independent_energy']['checks']
            self.audit.require(len(raw['native_history']) == 2*raw['order'], name+': complete native history')
            for stage in raw['stage_audits']:
                self.audit.require(stage[4] <= 1e-9 and max(stage[6:]) <= 1e-9, name+': unconverged native stage')
            self.audit.near(raw['native_history'][0], raw['state'][:raw['offsets'][raw['dimension']+2]],
                            name+': complete differential native output', 1e-10)
            if selected.get('restart') is None:
                self.audit.require(abs(raw['time']-selected.get('steps', 4)*selected.get('dt', .0025)) < 1e-11,
                                   name+': incomplete declared duration')
            record['physical_gates'] = physical_gate(raw)
            if source_mechanism:
                sources = [row['robust_source'] for row in raw['state_audits']]
                self.audit.require(any(row[4] > 0 for row in sources), name+': depletion did not activate source limiter')
                self.audit.require(all(abs(row[3]) < 1e-10*max(1., abs(row[1]), abs(row[2])) for row in sources),
                                   name+': redistribution does not retain zero metric integral')
                record['source_repair_verified'] = True
                record['interpretation'] = 'Source mechanism and accounting verified; covariance remains an independent recorded gate'
            record['status'] = 'verified'
            self.raws[name] = raw
        except Exception as error:
            record.update(status='audit_failed', audit_failure=str(error))
            self.receipt['failures'].append(dict(case=name, classification='verification_audit_failure', error=str(error)))
            self.save()
            return None
        self.save()
        return raw


def pilot(run):
    for name, generator, options in [
        ('reference', perturbation_session, dict(kind='none', eta=1, steps=1)),
        ('joint', perturbation_session, dict(kind='joint', eta=1, steps=1)),
        ('depletion', depletion_session, dict(eta=1, steps=1, source_mechanism=True)),
    ]:
        if run.case(name, generator, **options) is None:
            return False
    return True


def native(run):
    if not pilot(run):
        return False
    matrix = [
        ('joint-four', dict(kind='joint', eta=1, steps=4)),
        ('joint-dt2', dict(kind='joint', eta=1, dt=.00125, steps=8)),
        ('joint-dt4', dict(kind='joint', eta=1, dt=.000625, steps=16)),
        ('joint-p4', dict(kind='joint', eta=1, degree=4, steps=4)),
        ('joint-h2', dict(kind='joint', eta=1, refine=2, steps=4)),
        ('velocity', dict(kind='velocity', eta=0, steps=4)),
        ('k', dict(kind='k', eta=0, steps=4)),
        ('epsilon', dict(kind='epsilon', eta=0, steps=4)),
        ('interface', dict(kind='interface', eta=0, steps=4)),
        ('joint-3d', dict(kind='joint', eta=1, dimension=3, steps=1)),
    ]
    for name, options in matrix:
        if run.case(name, **options) is None:
            return False
    differences = [m7.temporal_difference(run.raws[a], run.raws[b]) for a, b in
                   [('joint-four', 'joint-dt2'), ('joint-dt2', 'joint-dt4')]]
    run.audit.require(differences[1]['differential'] < differences[0]['differential'],
                      'M8 native dt refinement did not reduce field difference')
    run.receipt['temporal_refinement'] = dict(comparisons=differences,
        ratio=differences[0]['differential']/differences[1]['differential'])
    parallel = run.case('joint-mpi2', ranks=2, kind='joint', eta=1, steps=4)
    if parallel is None:
        return False
    run.receipt['mpi_field_norms'] = state_difference(run.raws['joint-four'], parallel, 'M8 MPI all-field parity')
    first = run.case('checkpoint', kind='joint', eta=1, steps=2, checkpoint=run.output/'history')
    if first is None:
        return False
    resumed = run.case('restart', kind='joint', eta=1, steps=2, restart=first['checkpoint'])
    if resumed is None:
        return False
    run.receipt['restart_field_norms'] = state_difference(run.raws['joint-four'], resumed, 'M8 full-field restart')
    native_history_difference(run.raws['joint-four'], resumed, 'M8 complete native histories')
    # The original manufactured source is valid only where the persistent
    # stress increment is absent.  A nonzero-gamma MMS must be derived anew.
    if run.args.gamma == 0 and run.args.gamma_scale == 0:
        for family in ('flat', 'feet'):
            if run.case('mms-'+family, manufactured_session, family=family,
                        steps=2, dt=.005, time=0., full_energy=False) is None:
                return False
    else:
        flat = run.case('regularized-flat-mms', regularized_mms_session, steps=2,
                        dt=.005, full_energy=False)
        if flat is None:
            return False
        run.receipt['regularized_mms'] = verify_regularized_mms(flat)
        run.audit.checks += run.receipt['regularized_mms']['checks']
        run.receipt['regularized_curved_mms'] = dict(status='required_separate_gate',
            reason='Flat nonzero-gamma MMS does not establish curved traction/body consistency')
    run.save()
    return True


def mms(run):
    matrix = [('flat', dict()), ('flat-p4', dict(degree=4)),
              ('flat-h2', dict(refine=2)), ('flat-dt2', dict(dt=.0025, steps=4))]
    answers = {}
    for name, options in matrix:
        raw = run.case(name, regularized_mms_session, full_energy=False,
                       **dict(dict(steps=2, dt=.005), **options))
        if raw is None:
            return False
        answers[name] = verify_regularized_mms(raw)
        run.audit.checks += answers[name]['checks']
    parallel = run.case('flat-mpi2', regularized_mms_session, ranks=2,
                        steps=2, dt=.005, full_energy=False)
    if parallel is None:
        return False
    answers['flat-mpi2'] = verify_regularized_mms(parallel)
    run.audit.checks += answers['flat-mpi2']['checks']
    run.receipt['mms_mpi_field_norms'] = state_difference(run.raws['flat'], parallel, 'M8 genuine flat MMS MPI parity')
    first = run.case('flat-checkpoint', regularized_mms_session, steps=1, dt=.005,
                     checkpoint=run.output/'mms-history', full_energy=False)
    if first is None:
        return False
    resumed = run.case('flat-restart', regularized_mms_session, steps=1, dt=.005,
                       restart=first['checkpoint'], full_energy=False)
    if resumed is None:
        return False
    answers['flat-restart'] = verify_regularized_mms(resumed)
    run.audit.checks += answers['flat-restart']['checks']
    run.receipt['mms_restart_field_norms'] = state_difference(run.raws['flat'], resumed, 'M8 genuine flat MMS restart')
    native_history_difference(run.raws['flat'], resumed, 'M8 genuine flat MMS complete histories')
    run.receipt['regularized_mms'] = answers
    run.receipt['regularized_curved_mms'] = dict(status='required_separate_gate',
        reason='Flat all-field verification does not qualify the persistent regularization over curved terrain')
    run.save()
    return True


def flow(run, complete=False):
    results = {}
    matrix = [(dt, amplitude) for dt in ((.005, .0025) if complete else (.005,))
              for amplitude in (0., .001)]
    for dt, amplitude in matrix:
        name = ('control' if amplitude == 0 else 'pulse')+('-coarse' if dt == .005 else '-fine')
        steps = round(.3/dt) if complete else 1
        raw = run.case(name, flow_session, amplitude=amplitude, dt=dt, steps=steps)
        if raw is None:
            return False
        run.audit.require(raw['eta'] == 0 and raw['upstream_transition'] and raw['neutral_reference'],
                          'Unforced flow changed required transition/anchor/reference mode')
        run.audit.near(raw['transition_controls'][:7], [-3., -2.25, 2.5, 3.5, -2., .02, .1],
                       'M8 same physical transition/sensor controls')
        results[name] = raw
    pairs = {}
    for suffix in (('coarse', 'fine') if complete else ('coarse',)):
        pairs[suffix] = paired(results['control-'+suffix], results['pulse-'+suffix], run.audit, complete)
    run.receipt['paired_flow'] = pairs
    if complete:
        run.receipt['paired_dt_refinement'] = time_comparison(
            results['control-coarse'], results['pulse-coarse'],
            results['control-fine'], results['pulse-fine'], pairs['coarse'], pairs['fine'], run.audit)
    run.receipt['transition_passage_established'] = bool(complete)
    run.receipt['long_time_qualification'] = False
    run.receipt['duration_limitation'] = 'The .3s passage fixture is not the required multi-turnover/multi-domain-flow-time robust qualification'
    run.save()
    return True


def prerequisite(path, run, terrain=False):
    run.audit.require(path is not None and path.is_file(), 'A prior source-bound candidate qualification receipt is required')
    candidate = json.loads(path.read_text())
    run.audit.require(candidate.get('candidate_qualified') is True and candidate.get('implementation_verified') is True,
                      'Prior candidate has not passed all robust qualification gates')
    run.audit.require(candidate.get('parameters') == parameters(run.args), 'Prior qualification selected different model parameters')
    run.audit.require(candidate.get('source_commit') == run.receipt['source_commit'], 'Prior qualification belongs to another source revision')
    if terrain:
        run.audit.require(candidate.get('long_time_qualification') is True, 'Physical terrain validation requires complete long-time qualification')
    run.receipt['qualification_prerequisite'] = dict(path=str(path), sha256=sha(path))
    return candidate


def terrain(run):
    prerequisite(run.args.qualification, run, terrain=True)
    run.audit.require(run.args.terrain_sessions is not None and run.args.terrain_sessions.is_file(),
                      'A versioned revised-native Cao--Tamura floor/hill session manifest is required')
    manifest = json.loads(run.args.terrain_sessions.read_text())
    run.audit.require(manifest.get('protocol_version') == VERSION, 'Terrain manifest does not use versioned M8 data split')
    run.audit.require(manifest.get('parameters') == parameters(run.args), 'Terrain sessions select a different model')
    run.audit.require(manifest.get('calibration_data_role') == 'flat_only' and
                      manifest.get('reserved_hill_used_for_selection') is False,
                      'Terrain model selection/calibration consumed reserved validation data')
    # This entry point cannot fall back to historical CV-strip cases.  A
    # future physical adapter supplies complete revised-native session paths,
    # development windows and identical floor/hill physical input hashes.
    run.audit.require(manifest.get('backend') == 'RANSEquilibriumWallEvolution', 'Wrong wall representation in physical manifest')
    groups = manifest.get('pairs', [])
    run.audit.require(bool(groups), 'No paired physical floor/hill sessions supplied')
    run.receipt['terrain_manifest'] = dict(path=str(run.args.terrain_sessions), sha256=sha(run.args.terrain_sessions))
    for group in groups:
        run.audit.require(group['surface'] in ('S', 'R') and group['route'] in ('EPTransition', 'AllBaseline'), 'Invalid physical surface/route')
        run.audit.require(group['floor_inputs_sha256'] == group['hill_inputs_sha256'], 'Physical floor/hill inputs differ')
        run.audit.require(group['unforced_flow_times'] >= 5 and group['developed'] is True,
                          'Undeveloped physical trajectory cannot enter a validation score')
        for kind in ('floor', 'hill'):
            path = Path(group[kind]['report'])
            run.audit.require(path.is_file() and sha(path) == group[kind]['sha256'], 'Physical native report provenance mismatch')
            raw = json.loads(path.read_text())
            run.audit.require(not raw.get('failure') and raw['eta'] == 0 and raw['model_variant'] == run.args.variant,
                              'Stopped, anchored or wrong-model physical report')
            m7.audits(raw, group['surface']+'-'+group['route']+'-'+kind)
            run.audit.require(physical_gate(raw)['covariance_pass'], 'Physical native covariance gate failed')
    run.receipt['physical_validation'] = dict(status='developed_reports_verified_scoring_separate',
        actual_pairs=len(groups), experimental_accuracy_claimed=False,
        limitation='Observable extraction, sampling uncertainty and every reserved score must be supplied by the versioned native adapter')
    return True


def provenance():
    revision = subprocess.check_output(['git', '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip()
    status = subprocess.check_output(['git', '-C', str(SOURCE), 'status', '--porcelain'], text=True)
    if status:
        raise RuntimeError('Numerical campaign requires a clean committed isolated source')
    for name, expected in m7.FROZEN.items():
        if sha(SOURCE/name) != expected:
            raise RuntimeError('Frozen original model/contract changed: '+name)
    return revision


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=('protocol', 'pilot', 'native', 'mms', 'flow-pilot', 'flow', 'terrain'), default='pilot')
    p.add_argument('--variant', choices=VARIANTS, default='CubicScaled')
    p.add_argument('--beta', type=float, default=.4)
    p.add_argument('--gamma', type=float, default=0.)
    p.add_argument('--gamma-scale', type=float, default=0.)
    p.add_argument('--theta', type=float, default=.2)
    p.add_argument('--delta', type=float, default=.1)
    p.add_argument('--smooth', action='store_true')
    p.add_argument('--qualification', type=Path)
    p.add_argument('--terrain-sessions', type=Path)
    p.add_argument('--diagnostic', action='store_true', help='Allow unqualified flow observations, with no robust/physical promotion')
    p.add_argument('--timeout', type=float, default=3600.)
    args = p.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Numerical qualification must run on the authoritative Linux server')
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    revision = provenance()
    protocol = contract(SOURCE, args.mode, parameters(args))
    protocol_path = output/'preregistration.json'
    protocol_path.write_text(json.dumps(protocol, indent=2, allow_nan=False)+'\n')
    solver = args.build.resolve()/'IncNavierStokesSolver'
    receipt = dict(schema=1, milestone='M8', protocol_version=VERSION, mode=args.mode,
        source_commit=revision, source_clean_start=True, parameters=parameters(args),
        program=dict(path=str(solver), sha256=sha(solver)),
        preregistration=dict(path=str(protocol_path), sha256=sha(protocol_path)),
        started_utc=datetime.now(timezone.utc).isoformat(), commands=[], cases=[], failures=[],
        implementation_verified=False, candidate_qualified=False,
        long_time_qualification=False, experimental_accuracy_claimed=False,
        connection='tailscale:100.105.247.110', numerical_threads=THREADS,
        driver_sha256=sha(__file__))
    run = Runner(args, output, receipt)
    run.save(); success = False
    try:
        if args.mode == 'protocol':
            success = True
        elif args.mode == 'pilot':
            success = pilot(run)
        elif args.mode == 'native':
            success = native(run)
        elif args.mode == 'mms':
            success = mms(run)
        elif args.mode in ('flow-pilot', 'flow'):
            if not args.diagnostic:
                prerequisite(args.qualification, run)
            receipt['diagnostic_unqualified_flow'] = args.diagnostic
            success = flow(run, complete=args.mode == 'flow')
        else:
            success = terrain(run)
    except Exception as error:
        receipt['failures'].append(dict(classification='qualification_or_prerequisite_failure', error=str(error)))
    receipt['implementation_verified'] = bool(success)
    physical_failures = [dict(case=row['name'], gates=row['physical_gates']) for row in receipt['cases']
                         if row.get('physical_gates') and not row['physical_gates']['covariance_pass']
                         and not row.get('source_mechanism_only')]
    receipt['finite_covariance_failures'] = physical_failures
    # A pilot/native short fixture cannot self-promote to a robust terrain
    # candidate without state-envelope and complete long-time qualification.
    receipt['candidate_qualified'] = False
    receipt['qualification_status'] = ('implementation_verified_candidate_unqualified' if success else 'failed_or_stopped')
    if physical_failures:
        receipt['qualification_status'] += '_finite_covariance_gate_failed'
    receipt['remaining_gates'] = ['certified admissible constitutive state envelope',
        'outer turbulence positivity and high-order discrete admissibility',
        'complete long-time perturbation/refinement qualification',
        'qualified versioned native paired terrain observations and reserved scores']
    try:
        run.audit.require(provenance() == revision, 'M8 source revision changed during execution')
        run.audit.require(sha(solver) == receipt['program']['sha256'], 'M8 solver binary changed during execution')
        receipt['source_clean_end'] = True
    except Exception as error:
        receipt['implementation_verified'] = False
        receipt['failures'].append(dict(classification='provenance_failure', error=str(error)))
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    run.save()
    print(json.dumps(dict(output=str(output), implementation_verified=receipt['implementation_verified'],
        candidate_qualified=receipt['candidate_qualified'], qualification_status=receipt['qualification_status'],
        cases=len(receipt['cases']), checks=run.audit.checks, failures=receipt['failures']), allow_nan=False), flush=True)
    return 0 if receipt['implementation_verified'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
