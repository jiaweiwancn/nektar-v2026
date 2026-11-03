#!/usr/bin/env python3
"""Preregistered paired unforced finite-flow passage through the M6 mask.

No manufactured source, body force, changed gate or modified model completes
the pulse. Paired physical polynomial norms distinguish the perturbation from
the unanchored control's own evolution. An incomplete Newton solve is not a
physical-instability finding. Explicitly rejected negative sampled states are
retained as physical admissibility observations, without claiming passage.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from cls_epawf_m6_cases import transition_session
from cls_epawf_m3_cases import set_info, set_parameter, write
from validate_cls_epawf_m4 import spaces
import validate_cls_epawf_m7 as native_audit
from cls_epawf_m7_energy_oracle import verify as verify_energy

PHYSICS_REVISION = 'bf202be'
XMIN, XMAX, ZMAX, H, CENTER, WIDTH, US, AMPLITUDE, DURATION = -4., 2., 4., 2., -3.5, .3, .4, .001, .3
START, END, BASELINE = -3., -2.25, -2.


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def plan(mode):
    return dict(schema=1, milestone='M7', experiment='paired unforced finite-flow flat Upstream passage',
        preregistered_before_execution=True, physics_revision=PHYSICS_REVISION, mode=mode,
        domain=dict(x=[XMIN, XMAX], z=[0., ZMAX], H=H, degree=3, dimension=2),
        mask=dict(start=START, end=END, baseline_faces=BASELINE, sensor_distances=[2.5, 3.5],
                  state_on=.02, state_off=.1, eta=0., chi='prescribed Upstream product'),
        reference='Neutral native constrained projection; identical frozen reference and ordinary essential data in each pair',
        forcing='none: all pre-existing FUNCTION blocks removed; only InitialPerturbation functions added',
        nonlinear_solver=dict(maximum_iterations=80,tolerance=2e-12,coupled_acceptance=1e-9,
            rationale='iteration budget extended after preserved 35-iteration stop; physical operator and stopping tolerance unchanged'),
        pulse=dict(amplitude=AMPLITUDE, center=CENTER, nominal_width=WIDTH,
            definition='psi=A*us*sigma*f(x)*sin(pi*z/4)^2; u=partial_z psi,w=-partial_x psi; k/epsilon localized differential fields',
            envelope='essential-boundary polynomial bubble times Gaussian with first log-derivative cancelled at center; f(center)=1,fprime(center)=0',
            epsilon='outer normal bubble vanishes at H/top; actual current trace is supplied by coupled Lift'),
        matrix=dict(preflight=[dict(amplitude=a, dt=.005, steps=1) for a in (0., AMPLITUDE)],
                    complete=[dict(amplitude=a, dt=dt, steps=round(DURATION/dt)) for dt in (.005, .0025) for a in (0., AMPLITUDE)]),
        measurement='physical Jacobi/GL squared norm of paired deltaU/us,deltaK/K0,deltaEpsilon/EH; moment in physical x',
        acceptance=dict(all_accepted_state_native_audits=True, independent_upstream_mass_work_audit=True,
            nominal_matching_velocity_travel_exceeds_mask_width=True, initial_centroid_before=START,
            initial_upstream_norm_fraction_minimum=.75, final_centroid_after=END,
            final_downstream_norm_fraction_minimum=.5, fine_vs_coarse_paired_field_relative_maximum=.2,
            fine_vs_coarse_centroid_maximum=.15, fixed_time_step_refinement='two steps; bounded agreement, no order fit'),
        failure_policy='always preserve command/native diagnostic paths; stop before dependent full cases if a preflight fails; stop at first full-case failure; no automatic physical pass',
        interpretation='a coupled perturbation-norm drift witness in finite flow, not Lagrangian particle tracking or proof that advection alone causes the redistribution',
        resources=dict(ranks=1, numerical_threads=1, sequential_commands=True))


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, tolerance=1e-10):
        actual, expected = np.asarray(actual), np.asarray(expected)
        self.require(actual.shape == expected.shape and np.all(np.isfinite(actual)) and
                     np.all(np.isfinite(expected)), 'Invalid '+name)
        error = float(np.max(np.abs(actual-expected), initial=0.))/max(1., float(np.max(np.abs(expected), initial=0.)))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: {error} > {tolerance}')


def session(template, destination, amplitude, dt, steps):
    transition_session(template, destination, family='flat', degree=3, dimension=2,
                       steps=steps, dt=dt, order=2, time=0., probe=True,
                       quadrature=8, normal_order=16)
    root = ET.parse(destination).getroot()
    conditions = root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):
        conditions.remove(node)
    set_info(root, 'RANSWallReference', 'Neutral')
    for name, value in {'RANSWallDiagnostics': 1, 'RANSWallEta': 0, 'RANSWallChi': 1,
                        'RANSWallSensorOn': .02, 'RANSWallSensorOff': .1,
                        'RANSWallMaxIterations': 80}.items():
        set_parameter(root, name, value)
    denominator = (CENTER-XMIN)**2*(XMAX-CENTER)**2
    beta = 2/(CENTER-XMIN)-2/(XMAX-CENTER)
    polynomial = f'((x-({XMIN}))^2*(({XMAX})-x)^2/({denominator:.17g}))'
    polynomial_derivative = f'(2*(x-({XMIN}))*(({XMAX})-x)*(({XMIN+XMAX})-2*x)/({denominator:.17g}))'
    displacement = f'(x-({CENTER}))'
    gaussian = f'exp(-({displacement}/({WIDTH}))^2-({beta:.17g})*{displacement})'
    f = f'({polynomial}*{gaussian})'
    fp = f'(({polynomial_derivative}-{polynomial}*(2*{displacement}/({WIDTH*WIDTH:.17g})+({beta:.17g})))*{gaussian})'
    g = f'(sin(PI*y/({ZMAX})))^2'
    ux = f'({amplitude*US*WIDTH*math.pi/ZMAX:.17g})*{f}*sin(2*PI*y/({ZMAX}))'
    uz = f'({-amplitude*US*WIDTH:.17g})*{fp}*{g}'
    k = f'({amplitude:.17g})*Qstar*0.16*{f}*{g}'
    epsilon = f'({amplitude*US**3/H:.17g})/Kappa*{f}*(sin(PI*(y-({H}))/({ZMAX-H})))^2'
    for region in ('Inner', 'Outer'):
        node = ET.SubElement(conditions, 'FUNCTION', {'NAME': 'InitialPerturbation'+region})
        values = dict(u=ux, v=uz, k=k, eps_tilde=epsilon if region == 'Outer' else '0', p='0')
        for variable, value in values.items():
            ET.SubElement(node, 'E', {'VAR': variable, 'VALUE': value})
    names = [node.attrib['NAME'] for node in conditions.findall('FUNCTION')]
    if names != ['InitialPerturbationInner', 'InitialPerturbationOuter']:
        raise RuntimeError('Unforced flow session retained unexpected functions')
    return write(root, destination)


def failure_observation(raw, code):
    actual = []
    audits = raw.get('state_audits', [])
    for index, record in enumerate(audits):
        for name, row in record.items():
            if name in ('admissible', 'strong_continuity', 'interface', 'coupled_norms', 'energy',
                        'roughness', 'wall_sources', 'momentum_k_load_and_boundary_reaction') or len(row) != 6:
                continue
            if row[0] <= 0:
                continue
            negative_k = row[1] < 0
            bad_outer_epsilon = 'outer' in name and row[5] > 0 and row[2] <= 0
            if negative_k or bad_outer_epsilon:
                actual.append(dict(audit=index, locations=name, count=row[0], minimum_k=row[1],
                                   minimum_epsilon=row[2], invalid_count=row[5]))
    accepted = len(raw.get('accepted_states', []))
    times = raw.get('audit_times', [])
    last = times[accepted-1] if accepted and len(times) >= accepted else None
    return dict(classification='explicit_sampled_physical_admissibility_failure' if actual else
                               'unqualified_numerical_solver_or_initialization_failure',
                explicit_invalid_samples=actual, native_failure=raw.get('failure'), exit_code=code,
                accepted_states=accepted, last_accepted_time=last,
                last_reported_time=raw.get('time'), rejected_state_present=bool(raw.get('rejected_state')),
                pulse_passage_established=False,
                limitation='Newton nonconvergence alone is not a positivity proof or a continuum-instability proof')


class Runner:
    def __init__(self, build, output, audit, receipt):
        self.build, self.output, self.audit, self.receipt = build, output, audit, receipt
        self.environment = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH', '')
        receipt['environment'] = dict(THREADS, TMPDIR=self.environment.get('TMPDIR'),
                                      LD_LIBRARY_PATH=self.environment['LD_LIBRARY_PATH'],
                                      cpu_affinity=sorted(os.sched_getaffinity(0)))

    def save(self):
        self.receipt.update(checks=self.audit.checks, errors=self.audit.errors)
        (self.output/'validation.json').write_text(json.dumps(self.receipt, indent=2, allow_nan=False)+'\n')

    def command(self, name, argv, cwd, timeout=60):
        argv, log, start = list(map(str, argv)), self.output/(name+'.log'), time.monotonic()
        code, error = 127, None
        try:
            with log.open('w') as stream:
                completed = subprocess.run(argv, cwd=cwd, env=self.environment, stdout=stream,
                                           stderr=subprocess.STDOUT, timeout=timeout, check=False)
                code = completed.returncode
        except subprocess.TimeoutExpired:
            code, error = 124, 'Command timed out'
        except OSError as exception:
            error = str(exception)
        if error:
            with log.open('a') as stream:
                stream.write('\n'+error+'\n')
        record = dict(name=name, argv=argv, cwd=str(cwd), exit_code=code, ranks=1, threads=THREADS,
                      seconds=time.monotonic()-start, log=str(log), log_sha256=digest(log))
        self.receipt['commands'].append(record)
        self.save()
        print(f'{name}: exit {code}, {record["seconds"]:.2f}s', flush=True)
        return code, log.read_text().strip(), record

    def git(self, source, name, *arguments):
        code, text, _ = self.command(name, ['git', '-c', 'gc.auto=0', '-C', source, *arguments], source)
        self.audit.require(code == 0, 'Git metadata command failed: '+name)
        return text

    def case(self, name, template, amplitude, dt, steps):
        directory = self.output/name
        directory.mkdir()
        xml = session(template, directory/'case.xml', amplitude, dt, steps)
        record = dict(name=name, amplitude=amplitude, dt=dt, requested_steps=steps,
                      xml=str(xml), xml_sha256=digest(xml), status='running')
        self.receipt['cases'].append(record)
        code, _, command = self.command(name, [self.build/'IncNavierStokesSolver', xml], directory,
                                        timeout=1800 if steps == 1 else 10800)
        path = directory/'case.evolution.json'
        raw = {}
        if path.is_file():
            record.update(native_report=str(path), native_report_sha256=digest(path))
            try:
                raw = json.loads(path.read_text())
            except (ValueError, OSError) as error:
                record['native_report_error'] = str(error)
        record['exit_code'] = code
        record['seconds'] = command['seconds']
        if code or not raw or raw.get('failure'):
            record.update(status='stopped', observation=failure_observation(raw, code))
            self.receipt['failures'].append(dict(case=name, **record['observation']))
            self.save()
            return None
        try:
            self.audit.require(raw['eta'] == 0 and raw['upstream_transition'] and raw['neutral_reference'], 'Wrong unforced flow mode')
            self.audit.require(raw['chi'] == 1 and raw['transition_controls'][:5] == [START, END, 2.5, 3.5, BASELINE], 'Changed physical mask/sampler controls')
            self.audit.near(raw['transition_controls'][5:7], [.02, .1], 'unchanged_state_thresholds')
            for sample in raw['transition_samples']:
                if sample['input'][0] >= END:
                    self.audit.near(sample['gates'][0], 0., 'no_downstream_reactivation', 1e-15)
            record['native_audit'] = native_audit.audits(raw, name)
            record['independent_energy'] = verify_energy(raw)
            self.audit.require(abs(raw['time']-steps*dt) < 1e-11, 'Native flow did not complete its requested physical duration')
            record['status'] = 'passed'
        except Exception as error:
            record.update(status='audit_failed', audit_failure=str(error))
            self.receipt['failures'].append(dict(case=name, classification='verification_audit_failure', error=str(error)))
            self.save()
            return None
        self.save()
        return raw


def norm_matrices(field):
    """Independent physical polynomial integrals, split at every mask cut."""
    matrices = {name: np.zeros((field.size, field.size)) for name in ('mass', 'x_moment', 'upstream', 'mask', 'downstream', 'baseline')}
    z, w = np.polynomial.legendre.leggauss(16)
    for element in field.elements:
        cuts = sorted({element['lower'][0], element['upper'][0]}|
                      {cut for cut in (START, END, BASELINE) if element['lower'][0] < cut < element['upper'][0]})
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            points, weights = [], []
            for i, j in itertools.product(range(16), repeat=2):
                points.append([lo+(z[i]+1)*(hi-lo)/2, 0., element['lower'][2]+(z[j]+1)*(element['upper'][2]-element['lower'][2])/2])
                weights.append(w[i]*w[j]*(hi-lo)*(element['upper'][2]-element['lower'][2])/4)
            points, weights = np.asarray(points), np.asarray(weights)
            ids, phi, _ = field.rows(element, points)
            measures = dict(mass=weights, x_moment=weights*points[:, 0],
                            upstream=weights*(points[:, 0] < START),
                            mask=weights*((points[:, 0] >= START)&(points[:, 0] < END)),
                            downstream=weights*(points[:, 0] >= END), baseline=weights*(points[:, 0] >= BASELINE))
            for name, measure in measures.items():
                matrices[name][np.ix_(ids, ids)] += phi.T@(measure[:, None]*phi)
    return matrices


def paired(control, perturbed, audit, passage):
    audit.require(control['dimension'] == perturbed['dimension'] == 2 and
                  control['offsets'] == perturbed['offsets'] and control['fields'] == perturbed['fields'], 'Paired native spaces/orderings differ')
    audit.near(control['fixed_reference'], perturbed['fixed_reference'], 'identical_paired_fixed_reference', 1e-12)
    times = np.asarray(control['audit_times'])
    audit.near(times, perturbed['audit_times'], 'paired_physical_times', 1e-12)
    audit.require(len(control['accepted_states']) == len(perturbed['accepted_states']) == len(times), 'Missing paired accepted state')
    fields, offsets = spaces(control), control['offsets']
    scales = [US, US, control['parameters'][5], US**3/(control['parameters'][7]*H)]
    matrices = [norm_matrices(field) for field in fields[:4]]
    records = []
    for index, (c, p) in enumerate(zip(control['accepted_states'], perturbed['accepted_states'])):
        delta = np.asarray(p)-c
        values = {name: 0. for name in matrices[0]}
        kinetic = 0.
        for block in range(4):
            part = delta[offsets[block]:offsets[block+1]]/scales[block]
            for name, matrix in matrices[block].items():
                values[name] += float(part@matrix@part)
            if block < 2:
                physical = delta[offsets[block]:offsets[block+1]]
                kinetic += .5*float(physical@matrices[block]['mass']@physical)
        audit.require(values['mass'] > 1e-25, 'Paired perturbation norm vanished or is unresolved')
        centroid = values['x_moment']/values['mass']
        fractions = {name: values[name]/values['mass'] for name in ('upstream', 'mask', 'downstream', 'baseline')}
        audit.near(fractions['upstream']+fractions['mask']+fractions['downstream'], 1., 'physical_partition_norm_sum', 2e-12)
        audit.require(XMIN-1e-10 <= centroid <= XMAX+1e-10, 'Unphysical centroid')
        records.append(dict(time=float(times[index]), physical_centroid=centroid, scaled_squared_norm=values['mass'],
                            region_fractions=fractions, physical_velocity_perturbation_energy=kinetic,
                            resolved_plus_k_energy_difference=perturbed['state_audits'][index]['energy'][8]-control['state_audits'][index]['energy'][8]))
    UH = fields[0].evaluate(control['fixed_reference'], offsets[0], [CENTER, 0., H], 1)[0]
    audit.require(UH > 0 and UH*DURATION > END-START, 'Nominal finite flow cannot span the prescribed mask')
    audit.require(records[0]['physical_centroid'] < START and records[0]['region_fractions']['upstream'] > .75,
                  'Initial native projected pulse is not localized on the upstream plateau')
    if passage:
        audit.require(records[-1]['physical_centroid'] > END and records[-1]['region_fractions']['downstream'] > .5,
                      'Completed paired trajectory did not establish the preregistered downstream passage')
    return dict(matching_velocity=UH, nominal_travel_at_T_03=UH*DURATION,
                duration=float(times[-1]-times[0]), physical_centroid_drift=records[-1]['physical_centroid']-records[0]['physical_centroid'],
                passage_established=bool(passage), observations=records,
                norm='integral [(deltaU/us)^2+(deltaK/K0)^2+(deltaEpsilon/EH)^2] with outer-only epsilon',
                interpretation='paired coupled perturbation-norm redistribution; not Lagrangian parcel tracking')


def time_comparison(coarse_c, coarse_p, fine_c, fine_p, coarse, fine, audit):
    audit.require(coarse_c['fields'] == fine_c['fields'] and coarse_c['offsets'] == fine_c['offsets'], 'dt refinement changed native space')
    coarse_t = np.asarray([record['time'] for record in coarse['observations']])
    fine_t = np.asarray([record['time'] for record in fine['observations']])
    audit.near(coarse_t, fine_t[::2], 'same_time_dt_refinement_epochs', 1e-12)
    delta = (np.asarray(coarse_p['state'])-coarse_c['state'])-(np.asarray(fine_p['state'])-fine_c['state'])
    reference = np.asarray(fine_p['state'])-fine_c['state']
    fields, offsets = spaces(fine_c), fine_c['offsets']
    scales = [US, US, fine_c['parameters'][5], US**3/(fine_c['parameters'][7]*H)]
    square = norm = 0.
    for block in range(4):
        mass = fields[block].mass()
        a = delta[offsets[block]:offsets[block+1]]/scales[block]
        b = reference[offsets[block]:offsets[block+1]]/scales[block]
        square += float(a@mass@a)
        norm += float(b@mass@b)
    error = math.sqrt(max(0., square)/norm)
    centroid = abs(coarse['observations'][-1]['physical_centroid']-fine['observations'][-1]['physical_centroid'])
    audit.require(error < .2 and centroid < .15, 'Paired physical passage not retained with the preregistered dt-half agreement')
    return dict(paired_final_field_relative_error=error, final_centroid_difference=centroid,
                coarse_dt=.005, fine_dt=.0025, conclusion='two-step bounded agreement; no fitted temporal order')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('pilot', 'all'), default='pilot')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific flow passage computations belong on Linux')
    output, build = args.output.resolve(), args.build.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output/'experiment-plan.json').write_text(json.dumps(plan(args.mode), indent=2)+'\n')
    audit = Audit()
    receipt = dict(schema=1, milestone='M7', success=False, verification_status='failed', mode=args.mode,
                   started_utc=datetime.now(timezone.utc).isoformat(), commands=[], cases=[], pairs={}, failures=[],
                   experiment_plan=str(output/'experiment-plan.json'), experiment_plan_sha256=digest(output/'experiment-plan.json'),
                   passage_established=False, qualification_scope='preflight only' if args.mode == 'pilot' else 'paired physical passage with dt-half control')
    run = Runner(build, output, audit, receipt)
    source = Path(__file__).resolve().parents[4]
    template = Path(__file__).with_name('EquilibriumWallFlat.xml')
    try:
        receipt['source_commit'] = run.git(source, 'git-start-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-start-status', 'status', '--porcelain=v1', '--untracked-files=all')
        audit.require(not status, 'Source is not clean before physical-flow evidence')
        receipt['source_clean_start'] = True
        receipt['compiled_physics_revision'] = run.git(source, 'git-physics-revision', 'rev-parse', PHYSICS_REVISION)
        changed = run.git(source, 'git-physics-diff', 'diff', '--name-only', PHYSICS_REVISION, '--', 'library', 'solvers/IncNavierStokesSolver')
        compiled_changes = [name for name in changed.splitlines() if Path(name).suffix in ('.cpp', '.h', '.hpp', '.cc', '.cxx', '.c')]
        audit.require(not compiled_changes, 'Compiled physics changed relative to locked '+PHYSICS_REVISION+': '+str(compiled_changes))
        program = build/'IncNavierStokesSolver'
        receipt['program'] = dict(path=str(program), sha256=digest(program))
        preflight_control = run.case('preflight-control', template, 0., .005, 1)
        preflight_pulse = run.case('preflight-pulse', template, AMPLITUDE, .005, 1)
        audit.require(preflight_control is not None and preflight_pulse is not None, 'Preflight stopped; full physical trajectory not attempted')
        receipt['pairs']['preflight'] = paired(preflight_control, preflight_pulse, audit, False)
        if args.mode == 'all':
            cases = {}
            for dt in (.005, .0025):
                tag = str(dt).replace('.', 'p')
                control = run.case('control-dt-'+tag, template, 0., dt, round(DURATION/dt))
                audit.require(control is not None, 'Unforced control stopped before completion; no paired passage claim')
                pulse = run.case('pulse-dt-'+tag, template, AMPLITUDE, dt, round(DURATION/dt))
                audit.require(pulse is not None, 'Perturbed unforced trajectory stopped before completion; no passage claim')
                pair = paired(control, pulse, audit, True)
                cases[dt] = control, pulse, pair
                receipt['pairs'][str(dt)] = pair
                run.save()
            c, p, pair = cases[.005]
            fc, fp, fpair = cases[.0025]
            receipt['time_refinement'] = time_comparison(c, p, fc, fp, pair, fpair, audit)
            receipt['passage_established'] = True
        receipt['success'] = True
    except Exception as error:
        receipt['failures'].append(dict(classification='incomplete_or_failed_qualification', error=str(error)))
    try:
        end = run.git(source, 'git-end-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-end-status', 'status', '--porcelain=v1', '--untracked-files=all')
        receipt.update(source_commit_end=end, source_clean_end=not status)
        audit.require(not status and end == receipt.get('source_commit'), 'Source changed during physical-flow evidence')
        if 'program' in receipt:
            receipt['program']['sha256_end'] = digest(receipt['program']['path'])
            audit.require(receipt['program']['sha256_end'] == receipt['program']['sha256'], 'Native solver binary changed during flow evidence')
    except Exception as error:
        receipt['success'] = False
        receipt['failures'].append(dict(classification='provenance_failure', error=str(error)))
    receipt.update(finished_utc=datetime.now(timezone.utc).isoformat(), verification_status='passed' if receipt['success'] else 'failed')
    run.save()
    print(json.dumps(dict(success=receipt['success'], passage_established=receipt['passage_established'], checks=audit.checks,
                          failures=receipt['failures'], output=str(output))), flush=True)
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
