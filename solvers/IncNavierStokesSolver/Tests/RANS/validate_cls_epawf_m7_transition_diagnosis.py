#!/usr/bin/env python3
"""Collect the preregistered unforced-flow amplitude/dt diagnosis after stops.

A successful collection verifies provenance and accepted prefixes; it does not
promote failed trajectories to completed passage. All original reports, logs,
rejected states and failure messages remain intact. Newton nonconvergence is
not a positivity proof. No source, gate, constitutive law or limiter is added.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import numpy as np
import validate_cls_epawf_m7_transition_flow as flow
import validate_cls_epawf_m7 as native
import cls_epawf_m7_energy_oracle as energy
from cls_epawf_m3_cases import set_parameter, write


def preregistration(previous, transverse):
    physical = flow.plan('all')
    physical.pop('failure_policy', None)
    return dict(schema=1, milestone='M7', experiment='complete unforced transition amplitude/dt diagnosis',
        registered_before_new_commands=True, prior_attempt=str(previous), strict_transverse=str(transverse),
        compiled_physics_revision=flow.PHYSICS_REVISION, final_requested_time=flow.DURATION,
        reused_case=dict(dt=.005, amplitude=0., steps=60, nonlinear_iterations=80, tolerance=2e-12),
        new_cases=[dict(name='coarse-pulse', dt=.005, amplitude=flow.AMPLITUDE, steps=60),
                   dict(name='fine-control', dt=.0025, amplitude=0., steps=120),
                   dict(name='fine-pulse', dt=.0025, amplitude=flow.AMPLITUDE, steps=120)],
        fixed_solver=dict(maximum_iterations=80, tolerance=2e-12, ranks=1, threads=flow.THREADS),
        physical_inputs=physical,
        analysis='every accepted prefix audited; paired physical norms only at common accepted epochs; no extrapolation or order claim',
        failure_policy='run all three remaining cases sequentially even after solver stops; preserve full native failure reports',
        collection_success='all available observations/provenance verified, independent of trajectory completion',
        passage_success='all four complete and the original full-time passage and dt agreement criteria pass',
        physical_failure_policy='explicit negative sampled values only; Newton limit/line-search rejection without such data is an unqualified numerical stop')


def read_json(path):
    return json.loads(Path(path).read_text())


def resolve(path, base):
    path = Path(path)
    return path.resolve() if path.is_absolute() else (base/path).resolve()


def check_session(path, audit, dt, amplitude):
    root = ET.parse(path).getroot()
    parameters = {}
    for item in root.findall('./CONDITIONS/PARAMETERS/P'):
        if item.text and '=' in item.text:
            name, value = item.text.split('=', 1)
            parameters[name.strip().upper()] = value.strip()
    expected = dict(RANSWALLMAXITERATIONS=80., RANSWALLNONLINEARTOLERANCE=2e-12,
                    RANSWALLDIAGNOSTICS=1., RANSWALLETA=0., RANSWALLCHI=1.,
                    RANSWALLSENSORON=.02, RANSWALLSENSOROFF=.1,
                    RANSWALLTRANSITIONSTART=flow.START, RANSWALLTRANSITIONEND=flow.END,
                    RANSWALLBASELINESTART=flow.BASELINE, RANSWALLSENSORDISTANCEA=2.5,
                    RANSWALLSENSORDISTANCEB=3.5, RANSWALLHEIGHT=flow.H, TIMESTEP=dt,
                    NUMSTEPS=round(flow.DURATION/dt))
    for name, value in expected.items():
        audit.require(name in parameters, 'Missing locked session parameter '+name)
        audit.near(float(parameters[name]), value, 'locked_xml_'+name, 1e-13)
    infos = {item.attrib['PROPERTY']: item.attrib['VALUE'] for item in root.findall('./CONDITIONS/SOLVERINFO/I')}
    audit.require(infos.get('RANSWallReference') == 'Neutral' and infos.get('RANSWallTransition') == 'Upstream' and
                  infos.get('RANSWallMap') == 'Flat', 'Changed unforced reference/geometry/transition')
    functions = root.findall('./CONDITIONS/FUNCTION')
    audit.require({item.attrib['NAME'] for item in functions} == {'InitialPerturbationInner', 'InitialPerturbationOuter'},
                  'Unforced diagnosis retained a source/body/MMS function')
    return dict(path=str(Path(path).resolve()), sha256=flow.digest(path), dt=dt, amplitude=amplitude,
                maximum_iterations=80, nonlinear_tolerance=2e-12, function_names=[item.attrib['NAME'] for item in functions])


def accepted_prefix(raw, audit, name):
    audit.require(raw.get('diagnostics_schema') == 1, name+': missing full native diagnostic report')
    states, records, times = raw.get('accepted_states', []), raw.get('state_audits', []), raw.get('audit_times', [])
    count = len(states)
    if not count:
        audit.require(bool(raw.get('failure')), name+': no accepted states and no failure explanation')
        return None, dict(accepted_states=0, scope='no accepted prefix exists to validate')
    audit.require(len(records) >= count and len(times) >= count, name+': accepted prefix data missing')
    audit.require(len(raw.get('stage_audits', [])) == count-1, name+': accepted stage count mismatch')
    audit.require(len(raw.get('step_audits', [])) == count-1, name+': completed step count mismatch')
    audit.near(raw['state'], states[-1], name+'_report_state_is_last_accepted', 1e-12)
    audit.near(raw['time'], times[count-1], name+'_report_time_is_last_accepted', 1e-12)
    if 'steps' in raw:
        audit.require(raw['steps'] == count-1, name+': native history step count differs from accepted prefix')
    view = dict(raw, failure='', state_audits=records[:count], audit_times=times[:count],
                accepted_states=states[:count], accepted_prefix_view_only=True)
    before = native.CHECKS
    summary = native.audits(view, name+' accepted-prefix-only')
    summary.update(scope='only the actual accepted prefix; failed trajectory not promoted',
                   native_failure_preserved_in_original=raw.get('failure'), accepted_prefix_states=count,
                   accepted_prefix_last_time=times[count-1], additional_audit_checks=native.CHECKS-before)
    return view, summary


def last_accepted_energy(raw, view, audit, name):
    """Direct volume/mass check, without inventing absent coefficient loads."""
    if view is None:
        return dict(verified=False, reason='no accepted state')
    state = np.asarray(view['state'])
    time = float(view['time'])
    count = len(view['accepted_states'])
    if count == 1:
        history, alpha = state, 0.
    else:
        history = np.asarray(raw.get('last_stage_history', []))
        audit.require(history.shape == state.shape, name+': actual last accepted history unavailable')
        audit.near(raw['last_stage_time'], time, name+'_actual_last_accepted_stage_epoch', 1e-12)
        alpha = float(raw['last_stage_alpha'])
        audit.require(alpha > 0, name+': missing actual accepted-stage alpha')
    fields, local = energy._fields(view), energy.Audit()
    storage, stored = energy._mass(view, fields, state, history, alpha)
    direct = energy._direct(view, fields, state, local)
    observed = np.asarray(view['state_audits'][-1]['energy'])
    total = storage+direct['stress_work_minus_k_source']+direct['resolved_advection']+direct['roughness_force_work']
    # observed[0] is explicitly a captured native scalar contraction. After
    # a failed solve there may be no exported coefficient load to reconstruct.
    expected = [storage, direct['stress_work_minus_k_source'], direct['resolved_advection'],
                direct['roughness_force_work'], total, observed[0]-total, stored]
    audit.near(observed[[1, 2, 3, 4, 6, 7, 8]], expected, name+'_last_accepted_direct_energy', 1e-10)
    audit.near(view['state_audits'][-1]['source_volume_outer'][3],
               direct['outer_covariance_minimum_on_declared_Gauss_rule'],
               name+'_independent_outer_covariance_eigenvalue',1e-10)
    return dict(verified=True, scope='last accepted actual history, physical mass and independent direct-volume components',
                time=time, alpha=alpha, storage=storage, stored_energy=stored, direct=direct,
                independently_reconstructed_total=total, native_captured_unmasked_work=float(observed[0]),
                direct_to_native_scalar_defect=float(observed[0]-total),
                coefficient_load_contraction_reconstructed=False,
                essential_reaction_reconstructed=False,
                limitation='failed report may omit final unmasked coefficient loads; native captured work/reaction are not fabricated as independent data',
                independent_internal_checks=local.checks, independent_internal_errors=local.errors)


def analyze(raw, record, audit):
    name = record['name']
    observation = flow.failure_observation(raw, record['exit_code'])
    completed = not raw.get('failure') and record['exit_code'] == 0 and abs(raw.get('time', -1)-flow.DURATION) < 1e-11
    record.update(trajectory_completed=completed, native_failure=raw.get('failure'), observation=observation)
    if completed:
        record['observation']['classification'] = 'requested_trajectory_completed'
    else:
        audit.require(record['exit_code'] != 0 and bool(raw.get('failure')),
                      name+': incomplete run has inconsistent exit/failure evidence')
    view, summary = accepted_prefix(raw, audit, name)
    record['accepted_prefix'] = summary
    record['independent_last_accepted_energy'] = last_accepted_energy(raw, view, audit, name)
    record['observation_verified'] = True
    record['status']='verified_completed' if completed else 'verified_stop'
    return view


def reuse_coarse(run, previous, receipt, audit, current_program):
    prior_path = previous/'validation.json'
    prior = read_json(prior_path)
    audit.require(prior.get('source_clean_start') and prior.get('source_clean_end') and
                  prior.get('source_commit') == prior.get('source_commit_end'), 'Previous attempt source provenance did not close cleanly')
    pin = prior['program']
    audit.require(pin['sha256'] == pin.get('sha256_end') == current_program,
                  'Previous coarse-control binary differs from the locked current solver')
    cases = [record for record in prior['cases'] if record.get('amplitude') == 0 and record.get('dt') == .005 and
             record.get('requested_steps') == 60]
    audit.require(len(cases) == 1, 'Previous attempt does not uniquely contain coarse control80')
    original = cases[0]
    xml, report = resolve(original['xml'], previous), resolve(original['native_report'], previous)
    audit.require(flow.digest(xml) == original['xml_sha256'] and flow.digest(report) == original['native_report_sha256'],
                  'Previous coarse-control XML/report changed')
    command = next(item for item in prior['commands'] if item['name'] == original['name'])
    log = resolve(command['log'], previous)
    audit.require(flow.digest(log) == command['log_sha256'], 'Previous coarse-control log changed')
    audit.require(command['exit_code'] == original['exit_code'] and command['ranks'] == 1 and
                  all(str(value) == '1' for value in command['threads'].values()), 'Previous coarse-control command allocation differs')
    record = dict(name='coarse-control', reused=True, exit_code=original['exit_code'],
                  amplitude=0., dt=.005, requested_steps=60, xml=check_session(xml, audit, .005, 0.),
                  native_report=str(report), native_report_sha256=flow.digest(report),
                  original_command=command, original_case=original)
    receipt['previous_attempt'] = dict(path=str(prior_path), sha256=flow.digest(prior_path),
                                      source_commit=prior['source_commit'], failures=prior.get('failures', []),
                                      success=prior.get('success'), program=pin)
    receipt['cases'].append(record)
    return read_json(report), record


def run_case(run, template, name, amplitude, dt, receipt, audit):
    directory = run.output/name
    directory.mkdir()
    xml = flow.session(template, directory/'case.xml', amplitude, dt, round(flow.DURATION/dt))
    root = ET.parse(xml).getroot()
    set_parameter(root, 'RANSWallMaxIterations', 80)
    set_parameter(root, 'RANSWallNonlinearTolerance', '2e-12')
    write(root, xml)
    record = dict(name=name, reused=False, amplitude=amplitude, dt=dt,
                  requested_steps=round(flow.DURATION/dt), xml=check_session(xml, audit, dt, amplitude), status='running')
    receipt['cases'].append(record)
    code, _, command = run.command(name, [run.build/'IncNavierStokesSolver', xml], directory, timeout=10800)
    report = directory/'case.evolution.json'
    record.update(exit_code=code, command=command, native_report=str(report))
    audit.require(report.is_file(), name+': no native diagnostic report; command/log retained')
    record['native_report_sha256'] = flow.digest(report)
    raw = read_json(report)
    return raw, record


def slice_view(view, count):
    result = dict(view)
    for name in ('accepted_states', 'state_audits', 'audit_times'):
        result[name] = view[name][:count]
    for name in ('stage_audits', 'step_audits'):
        result[name] = view[name][:count-1]
    result.update(state=result['accepted_states'][-1], time=result['audit_times'][-1], steps=count-1,
                  view_scope='common accepted epoch prefix only; no missing state extrapolated')
    # Those original fields belong to a different final epoch after cropping.
    for name in ('last_stage_history', 'last_stage_alpha', 'last_stage_time', 'energy_budget',
                 'final_stage_loads', 'final_stage_constraints'):
        result.pop(name, None)
    return result


def compare_prefixes(views, audit):
    paired_views, result = {}, {}
    for dt, control, pulse in ((.005, 'coarse-control', 'coarse-pulse'), (.0025, 'fine-control', 'fine-pulse')):
        if views.get(control) is None or views.get(pulse) is None:
            result[str(dt)] = dict(available=False, reason='a trajectory has no verified accepted prefix')
            continue
        count = min(len(views[control]['accepted_states']), len(views[pulse]['accepted_states']))
        c, p = slice_view(views[control], count), slice_view(views[pulse], count)
        audit.near(c['audit_times'], p['audit_times'], 'common_same_dt_prefix_epochs', 1e-12)
        pair = flow.paired(c, p, audit, False)
        pair.update(scope='common accepted epochs only', requested_final_time=flow.DURATION,
                    common_last_time=c['time'], full_passage_not_asserted=True)
        result[str(dt)] = pair
        paired_views[dt] = c, p
    differences = []
    if .005 in paired_views and .0025 in paired_views:
        cc, cp = paired_views[.005]
        fc, fp = paired_views[.0025]
        coarse, fine = result['0.005']['observations'], result['0.0025']['observations']
        audit.require(cc['fields'] == fc['fields'] and cc['offsets'] == fc['offsets'], 'Prefix dt comparison changed native space')
        fields, offsets = flow.spaces(cc), cc['offsets']
        scales = [flow.US, flow.US, cc['parameters'][5], flow.US**3/(cc['parameters'][7]*flow.H)]
        for i, row in enumerate(coarse):
            matches = [j for j, other in enumerate(fine) if abs(other['time']-row['time']) < 1e-11]
            if not matches:
                continue
            j = matches[0]
            a = np.asarray(cp['accepted_states'][i])-cc['accepted_states'][i]
            b = np.asarray(fp['accepted_states'][j])-fc['accepted_states'][j]
            square = norm = 0.
            for block in range(4):
                mass = fields[block].mass()
                aa = (a-b)[offsets[block]:offsets[block+1]]/scales[block]
                bb = b[offsets[block]:offsets[block+1]]/scales[block]
                square += float(aa@mass@aa)
                norm += float(bb@mass@bb)
            differences.append(dict(time=row['time'], paired_field_relative_difference=math.sqrt(max(0., square)/norm),
                                    centroid_difference=abs(row['physical_centroid']-fine[j]['physical_centroid'])))
    result['time_comparison'] = dict(scope='actual common accepted epochs only; no temporal-order or requested-T convergence claim',
                                     observations=differences)
    return result, paired_views


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--transverse', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific transition diagnosis belongs on Linux')
    output, build, previous, transverse = (value.resolve() for value in (args.output, args.build, args.previous, args.transverse))
    if not (previous/'validation.json').is_file() and (previous/'qualification'/'validation.json').is_file():
        previous = previous/'qualification'
    output.mkdir(parents=True, exist_ok=False)
    (output/'experiment-plan.json').write_text(json.dumps(preregistration(previous, transverse), indent=2)+'\n')
    audit = flow.Audit()
    receipt = dict(schema=1, milestone='M7', status='collecting', success=False, passage_qualified=False,
                   started_utc=datetime.now(timezone.utc).isoformat(), cases=[], commands=[], failures=[], pairs={},
                   collection_errors=[], physical_stability_gate='failed_as_specified',
                   experiment_plan_sha256=flow.digest(output/'experiment-plan.json'))
    run = flow.Runner(build, output, audit, receipt)
    source, template = Path(__file__).resolve().parents[4], Path(__file__).with_name('EquilibriumWallFlat.xml')
    views, raws = {}, {}
    try:
        receipt['source_commit'] = run.git(source, 'git-start-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-start-status', 'status', '--porcelain=v1', '--untracked-files=all')
        audit.require(not status, 'Audit source is not clean')
        receipt['source_clean_start'] = True
        receipt['compiled_physics_revision'] = run.git(source, 'git-physics-revision', 'rev-parse', flow.PHYSICS_REVISION)
        changed = run.git(source, 'git-physics-diff', 'diff', '--name-only', flow.PHYSICS_REVISION, '--', 'library', 'solvers/IncNavierStokesSolver')
        audit.require(not [name for name in changed.splitlines() if Path(name).suffix in ('.cpp', '.h', '.hpp', '.cc', '.cxx', '.c')],
                      'Compiled source differs from locked physics')
        program = build/'IncNavierStokesSolver'
        receipt['program'] = dict(path=str(program), sha256=flow.digest(program))
        strict = read_json(transverse)
        audit.require(strict.get('success') and strict.get('verification_status') == 'passed' and
                      strict.get('physical_stability_gate') == 'failed_as_specified', 'Missing passed strict transverse failure witness')
        receipt['strict_transverse_evidence'] = dict(path=str(transverse), sha256=flow.digest(transverse),
            checks=strict.get('checks'), source_commit=strict.get('source_commit'),
            physical_stability_gate=strict['physical_stability_gate'], relationship='separate definitive original-law failure witness; not inferred from these nonlinear stops')
        try:
            raw, record = reuse_coarse(run, previous, receipt, audit, receipt['program']['sha256'])
            raws[record['name']] = raw
            views[record['name']] = analyze(raw, record, audit)
        except Exception as error:
            receipt['collection_errors'].append(dict(case='coarse-control-reuse', error=str(error)))
        run.save()
        for name, amplitude, dt in (('coarse-pulse', flow.AMPLITUDE, .005), ('fine-control', 0., .0025), ('fine-pulse', flow.AMPLITUDE, .0025)):
            try:
                raw, record = run_case(run, template, name, amplitude, dt, receipt, audit)
                raws[name] = raw
                views[name] = analyze(raw, record, audit)
            except Exception as error:
                receipt['collection_errors'].append(dict(case=name, error=str(error)))
            run.save()
        pairs, common = compare_prefixes(views, audit)
        receipt['pairs'] = pairs
        receipt['all_trajectories_completed'] = len(receipt['cases']) == 4 and all(record.get('trajectory_completed') for record in receipt['cases'])
        if receipt['all_trajectories_completed'] and not receipt['collection_errors']:
            try:
                c, p = common[.005]
                fc, fp = common[.0025]
                coarse, fine = flow.paired(c, p, audit, True), flow.paired(fc, fp, audit, True)
                receipt['full_passage_metrics'] = dict(coarse=coarse, fine=fine,
                    time_refinement=flow.time_comparison(c, p, fc, fp, coarse, fine, audit))
                receipt['passage_qualified'] = True
            except Exception as error:
                receipt['passage_not_qualified_reason'] = str(error)
        else:
            receipt['passage_not_qualified_reason'] = 'At least one requested trajectory stopped or could not be verified; common prefixes are not completed passage'
        audit.require(len(receipt['cases']) == 4 and all(record.get('observation_verified') for record in receipt['cases']),
                      'Not all four requested observations were verified')
        audit.require(not receipt['collection_errors'], 'Diagnosis contains observation/provenance/accepted-prefix verification errors')
        receipt['success'] = True
    except Exception as error:
        receipt['collection_errors'].append(dict(case='diagnosis', error=str(error)))
    try:
        end = run.git(source, 'git-end-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-end-status', 'status', '--porcelain=v1', '--untracked-files=all')
        receipt.update(source_commit_end=end, source_clean_end=not status)
        audit.require(not status and end == receipt.get('source_commit'), 'Audit source changed during diagnosis')
        if 'program' in receipt:
            receipt['program']['sha256_end'] = flow.digest(receipt['program']['path'])
            audit.require(receipt['program']['sha256_end'] == receipt['program']['sha256'], 'Locked native solver changed during diagnosis')
    except Exception as error:
        receipt['success'] = False
        receipt['collection_errors'].append(dict(case='final-provenance', error=str(error)))
    receipt.update(status='diagnosis_collected' if receipt['success'] else 'diagnosis_incomplete',
                   verification_status='observations_verified' if receipt['success'] else 'failed',
                   finished_utc=datetime.now(timezone.utc).isoformat(),
                   interpretation='collection success is not trajectory success; no positivity proof is inferred from a solver stop')
    run.save()
    print(json.dumps(dict(status=receipt['status'], checks=audit.checks, passage_qualified=receipt['passage_qualified'],
                          collection_errors=receipt['collection_errors'], output=str(output))), flush=True)
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
