#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H2b operator, representation, mapped-startup and exact-cache checks."""
import argparse
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import time
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic, freeze, CampaignLock
from run_cao_tamura_budget_h2a import prepare, physical
from cls_epawf_m4_cases import tensor_mesh
from compare_cls_epawf_pimple import same_history
from diagnose_cls_accepted_covariance import checkpoint, common_packet, run_probe
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, set_info, set_parameter, write

X = [-.2, -.1, -.05, 0., .05, .1, .6]
Z = [0., .00025, .001, .003, .015, .03, .06, .1, .2, .25, .3, .4, .5, .8]


def case(folder, h1, steps=0, dt=5e-6, restart=None, geometry='hill',
         audit=True, split=None, pressure=0., old_mesh=False):
    xml = prepare(folder, h1, steps, dt, restart)
    t = ET.parse(xml).getroot()
    x = [-.2, -.1, 0., .1, .2] if old_mesh else X
    z = list(Z)
    if split:
        z.remove(.5)  # unchanged constant freestream; retain .25/.4 transitions
        z.append(split)
        z.sort()
    points = [8, 8] if old_mesh else [16, 16]
    if not (old_mesh and geometry == 'affine' and not split):
        tensor_mesh(t, z, 3, 2, x[-1]-x[0], 1., points=points, xknots=x)
    for k, v in dict(RANSWallGeometry='LegacyV1' if geometry == 'affine' else 'CompositeHillV1',
                     RANSWallMap='Flat' if geometry == 'affine' else 'CosineHill',
                     RANSWallAuditReuse='ExactObservationV1' if audit else 'Off').items():
        set_info(t, k, v)
    parameters = dict(RANSWallOutletTotalPressure=pressure)
    if not old_mesh:
        parameters.update(RANSWallMapAmplitude=.04 if geometry == 'hill' else 0.,
                          RANSWallMapLengthX=.1, RANSWallTangentialOrder=16, RANSWallNormalOrder=16)
    for k, v in parameters.items():
        set_parameter(t, k, v)
    write(t, xml)
    ni = sum(a <= .015 for a in z[1:]); no = len(z)-1-ni; nx = len(x)-1
    sizes = [(3*nx+1)*(3*(ni+no)+1)]*3 + [(3*nx+1)*(3*no+1), (nx+1)*(ni+1), (nx+1)*(no+1), 1]
    require(sum(sizes) <= 3000, 'H2b dense guard exceeded')
    definition = load(folder/'definition.json')
    definition.update(version='CaoTamuraBudgetH2bV1', case='cao-tamura-budget-h2b-'+geometry+('-near-wall-split' if split else ''),
                      geometry=geometry,
                      chart='inner physical r; outer eta=(zHat-Hw)/(Z-Hw)',
                      x=x, z=z, physical_steps=steps, dt=dt, size=sum(sizes),
                      field_sizes=sizes, points=points,
                      audit_reuse=audit, split=split,
                      native_size_estimate=sum(sizes), native_field_sizes=sizes,
                      mesh=dict(x=x, z=z, velocity_scalar_degree=3, pressure_degree=1, points=points,
                                revision='H2bV1', feet=[-.1,.1], crest=0., near_wall_split=split),
                      physical_map='affine flat' if geometry=='affine' else 'CompositeHillV1',
                      initialization='native finite/AWF radial setup' if geometry!='hill' else 'physical streamfunction lift, constrained physical-mass projection and native continuity',
                      hill_amplitude=.04 if geometry=='hill' else 0.,
                      requested_local_elapsed=steps*dt,
                      endpoint_time=None if restart else steps*dt,
                      xml_sha256=sha(xml), interpretation='operator/startup mesh; no hill validation score')
    atomic(folder/'definition.json', definition)
    atomic(folder/'h2b-definition.json', definition)
    (folder/'h2a-definition.json').rename(folder/'unexecuted-h2a-template-definition.json')
    return xml


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for n in ('solver', 'probe', 'h1', 'output'):
        p.add_argument('--'+n, type=Path, required=True)
    p.add_argument('--phase', choices=['operators', 'curved-operators', 'initialize', 'flat-limit-steps', 'startup', 'timing', 'representation'], default='operators')
    p.add_argument('--audit-mode', choices=['Off', 'ExactObservationV1'], default='ExactObservationV1')
    p.add_argument('--outlet-datum', type=float, default=0.)
    p.add_argument('--operators-root', type=Path)
    a = p.parse_args()
    require(platform.system() == 'Linux', 'authoritative Linux required')
    a.output.mkdir(parents=True, exist_ok=False)
    lock = CampaignLock(a.output)
    os.environ.update(THREADS)
    solver, image = freeze(a.solver.resolve(), a.output/'runtime')
    probe, pimage = freeze(a.probe.resolve(), a.output/'probe-runtime')
    result = dict(status='started', phase=a.phase, runtime=image, probe=pimage, tests={}, observations={})
    def save(**kw):
        result.update(kw); atomic(a.output/'progress.json', result)
    def operators(name, **kw):
        kw.setdefault('audit', a.audit_mode != 'Off')
        xml = case(a.output/name, a.h1, **kw)
        packet = dict(operation='h2b_operators', inner_stress_model='Linear', xml_sha256=sha(xml))
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        d, r = run_probe(probe, xml, packet, a.output/(name+'-probe'))
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        r['native_children_CPU_delta'] = dict(user=after.ru_utime-before.ru_utime, system=after.ru_stime-before.ru_stime)
        r['root_cumulative_children_peak_RSS_KiB'] = after.ru_maxrss
        result['observations'][name] = dict(receipt=r, checks=d['checks'], audit=d['audit'],
            state_sha256=d['state_sha256'], setup=d['seed_setup'], minimum_inner=d['seed_constitutive']['minimum_inner'],
            minimum_outer=d['seed_constitutive']['minimum_outer'])
        save(); return d
    def native(name, steps=0, dt=5e-6, cp=None, mutate=None, expected=True, **kw):
        xml = case(a.output/name, a.h1, steps, dt, cp, **kw)
        if mutate:
            t = ET.parse(xml).getroot(); mutate(t); write(t, xml)
        definition=load(xml.parent/'definition.json')
        definition.update(executed_xml_sha256=sha(xml), diagnostic_mutation=bool(mutate))
        atomic(xml.parent/'definition.json',definition)
        env = dict(os.environ, LD_LIBRARY_PATH=str(solver.parent/'lib'), TMPDIR=str(xml.parent),
                   NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        cmd = ['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(xml.parent/'usage.txt'), str(solver), str(xml)]
        before = time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:
            child = subprocess.Popen(cmd, cwd=xml.parent, env=env, stdout=log, stderr=subprocess.STDOUT,
                                     pass_fds=(lock.stream.fileno(),))
            atomic(xml.parent/'launch.json', dict(pid=child.pid, command=cmd))
            code = child.wait(timeout=1800)
        f = xml.with_suffix('.evolution.json'); q = load(f) if f.exists() else None
        receipt = dict(code=code, seconds=time.monotonic()-before, usage=(xml.parent/'usage.txt').read_text(),
                       report=str(f), report_sha256=sha(f) if f.exists() else None, xml_sha256=sha(xml),
                       failure=q.get('failure') if q else (xml.parent/'solver.log').read_text()[-3000:])
        atomic(xml.parent/'receipt.json', receipt)
        require((code == 0 and q and not q['failure']) == expected, 'unexpected native outcome '+name+': '+str(receipt))
        if expected:
            for row in q['pimple_attempts']:
                for k in ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
                          'epsilon_transport_solves', 'k_transport_solves', 'bound_calls'):
                    require(row[k] == 1, 'wrong physical schedule '+k)
                require(row['additional_governing_equation_solves'] == 0 and row['solved_epsilon_retained_bitwise'],
                        'extra governing solve or epsilon replacement')
        return q, receipt
    def observation(name, q, finite=False):
        cp = Path(q['checkpoint']); s = checkpoint(cp); report = a.output/name/'case.evolution.json'
        packet = dict(common_packet(cp), operation='finite_inflow_h1' if finite else 'h2b_checkpoint', inner_stress_model='Linear',
                      state=s['coupled'], history=s['coupled'], alpha=0., time=s['time'],
                      native_report_file=str(report), native_report_sha256=sha(report))
        d, r = run_probe(probe, report.parent/'case.xml', packet, a.output/(name+'-readback'))
        result['observations'][name] = dict(receipt=r, audit=d.get('audit', d.get('complete_audit')), spatial=d['selected_spatial_block_norms'],
            boundary=d['developing_boundary'], budget=d['native_turbulence_budget'])
        save(); return d
    save()
    try:
        if a.phase == 'curved-operators':
            hill = operators('composite-hill', pressure=a.outlet_datum)
            result['tests']['curved_operators'] = hill['checks']['status'] == 'passed'
            result['tests']['mapped_seed_audit'] = hill['audit']['accepted_awf_hard_admissible'] == [1]
            save(status='completed_operators' if result['tests']['mapped_seed_audit'] else 'operators_passed_seed_rejected')
        elif a.phase == 'operators':
            affine = operators('affine', geometry='affine', pressure=.013)
            flat = operators('composite-flat', geometry='flat', pressure=.013)
            require(affine['state_sha256'] == flat['state_sha256'], 'composite flat seed differs')
            differences = {}
            for key in ('physical_mass_action', 'native_pressure_action', 'native_spatial_residual',
                        'scalar_frozen_load', 'epsilon_frozen_load'):
                errors = [abs(x-y) for x, y in zip(affine['checks'][key], flat['checks'][key])]
                differences[key] = max(errors, default=0.)
                require(differences[key] < 1e-11, 'flat limit action differs '+key)
            comparable = {k: v for k, v in flat['audit'].items() if not k.endswith('_physical_point')}
            require(affine['audit'].keys() == comparable.keys(), 'flat-limit audit membership differs')
            audit_errors = {k: max((abs(x-y) for x, y in zip(v, comparable[k])), default=0.)
                            for k, v in affine['audit'].items()}
            require(max(audit_errors.values()) < 2e-12, 'flat-limit complete inventory/audit differs '+str(audit_errors))
            result['tests']['composite_flat_limit'] = dict(passed=True, maximum_errors=differences,
                tolerance=1e-11, audit_errors=audit_errors,
                reason='identical flat physical quadrature; explicit composite provenance retained; no approximate mass')
            save()
            hill = operators('composite-hill')
            result['tests']['curved_operators'] = hill['checks']['status'] == 'passed'
            result['tests']['mapped_seed_audit'] = hill['audit']['accepted_awf_hard_admissible'] == [1]
            save(status='completed_operators' if result['tests']['mapped_seed_audit'] else 'operators_passed_seed_rejected')
        elif a.phase == 'flat-limit-steps':
            affine0, _ = native('affine-h0', geometry='affine')
            flat0, _ = native('composite-flat-h0', geometry='flat')
            require(state_hash(affine0['state']) == state_hash(flat0['state']), 'flat-limit physical seeds differ')
            affine, ra = native('affine2', 2, cp=Path(affine0['checkpoint']), geometry='affine')
            whole, rw = native('composite-flat2', 2, cp=Path(flat0['checkpoint']), geometry='flat')
            difference = max(abs(x-y) for x, y in zip(affine['state'], whole['state']))
            require(difference < 1e-10, 'composite flat-limit stepping exceeds declared coefficient tolerance')
            first, _ = native('composite-flat1', 1, cp=Path(flat0['checkpoint']), geometry='flat')
            second, _ = native('composite-flat1-restart', 1, cp=Path(first['checkpoint']), geometry='flat')
            require(same_history(whole, second), 'composite split restart changed state/history bits')
            before = sha(Path(flat0['checkpoint'])/'State.bin')
            rejected, _ = native('composite-flat-assigned-reject', 1, cp=Path(flat0['checkpoint']), geometry='flat', expected=False,
                mutate=lambda t: set_parameter(t, 'RANSWallPimplePostAssignRejectStep', 1))
            require(same_history(flat0, rejected) and sha(Path(flat0['checkpoint'])/'State.bin') == before,
                'composite post-assignment rejection changed parent/accepted state/history')
            incompatible, _ = native('composite-geometry-mismatch', cp=Path(flat0['checkpoint']), geometry='flat', expected=False,
                mutate=lambda t: set_parameter(t, 'RANSWallMapAmplitude', .04))
            require(incompatible and 'signature mismatch' in incompatible['failure'], 'composite geometry restart guard bypassed')
            result['tests'].update(flat_limit_steps=True, split_restart_bitwise=True, post_assignment_rollback=True,
                incompatible_geometry_rejected=True)
            result['comparison'] = dict(affine=ra, composite=rw, maximum_coefficient_difference=difference,
                declared_tolerance=1e-10, physical_differences=physical(affine, whole), hill_physical_steps=0)
            save(status='completed_flat_limit_steps')
        elif a.phase == 'initialize':
            zero, receipt = native('hill-h0')
            d = observation('hill-h0', zero)
            result['tests']['native_seed_hard_admissible'] = d['audit']['accepted_awf_hard_admissible'] == [1]
            result['initialization'] = receipt
            save(status='completed_initialization')
        elif a.phase == 'startup':
            require(a.operators_root, 'operator qualification root required')
            gate = load(a.operators_root/'progress.json')
            require(gate['status'] == 'completed_operators' and gate['tests']['mapped_seed_audit'], 'mapped seed/operator gate not passed')
            zero, _ = native('hill-h0')
            observation('hill-h0', zero)
            cp = Path(zero['checkpoint']); original = sha(cp/'State.bin')
            whole, receipt = native('hill8', 8, cp=cp)
            first, _ = native('hill4', 4, cp=cp)
            second, _ = native('hill4-restart', 4, cp=Path(first['checkpoint']))
            require(same_history(whole, second), 'mapped split restart differs')
            bad, _ = native('hill-reject', 1, cp=cp, expected=False,
                mutate=lambda t: set_parameter(t, 'RANSWallPimplePostAssignRejectStep', 1))
            require(same_history(zero, bad) and sha(cp/'State.bin') == original, 'mapped rollback changed parent/accepted history')
            fine0, _ = native('hill-half-h0', dt=2.5e-6)
            require(state_hash(zero['state']) == state_hash(fine0['state']), 'mapped dt changed h0')
            fine, rr = native('hill-half16', 16, dt=2.5e-6, cp=Path(fine0['checkpoint']))
            observation('hill8', whole); observation('hill-half16', fine)
            result['tests'].update(split_restart=True, rejection_rollback=True, same_h0=True)
            result['temporal'] = dict(adjacent_physical_differences=physical(whole, fine), order=None,
                reason='two levels only; startup capability, not hill validation', coarse=receipt, fine=rr)
            for label, mutate in [('geometry', lambda t: set_parameter(t, 'RANSWallMapAmplitude', .03)),
                                  ('datum', lambda t: set_parameter(t, 'RANSWallOutletTotalPressure', .01))]:
                q, r = native('reject-'+label, cp=cp, mutate=mutate, expected=False)
                require(q and 'signature mismatch' in q['failure'], 'restart identity not enforced '+label)
            save(status='completed_startup')
        elif a.phase == 'representation':
            base, rb = native('retained-mesh-h0', geometry='affine', old_mesh=True)
            before = observation('retained-mesh-h0', base, finite=True)
            scales = before['provider']['wall_reference_scales']
            split = scales['positive_viscosity_intercept']
            require(0 < split < .00025, 'targeted scale is not in first retained panel')
            changed, rr = native('split-mesh-h0', geometry='affine', old_mesh=True, split=split)
            after = observation('split-mesh-h0', changed, finite=True)
            result['comparison'] = dict(scales=scales, before_receipt=rb, after_receipt=rr,
                before_projection=before['physical_projection_metrics'], after_projection=after['physical_projection_metrics'],
                before_budget=before['native_turbulence_budget'], after_budget=after['native_turbulence_budget'],
                before_profiles=before['physical_profiles'], after_profiles=after['physical_profiles'],
                steps=0, altered_physics=False, selection='no automatic mesh adoption from one traction number')
            save(status='completed_representation')
        else:
            # Alternation avoids mistaking cold-cache or machine-load effects
            # for speedup. Every run has identical complete hard checks/output.
            zero, _ = native('timing-h0', geometry='affine', old_mesh=True)
            rows = []
            for pair in range(3):
                states = []
                for on in ((False, True) if pair % 2 == 0 else (True, False)):
                    name = 'pair'+str(pair)+('-on' if on else '-off')
                    q, receipt = native(name, 8, cp=Path(zero['checkpoint']), geometry='affine', old_mesh=True, audit=on)
                    states.append(q); rows.append(dict(pair=pair, on=on, receipt=receipt, cache=q['evaluation_reuse']))
                require(same_history(*states), 'audit-cache changed state or native histories')
                for b, o in zip(states[0]['pimple_attempts'], states[1]['pimple_attempts']):
                    for k in ('history_sha256', 'incoming_sha256', 'post_flow_sha256', 'pre_k_sha256', 'raw_k_sha256', 'corrected_k_sha256'):
                        require(b[k] == o[k], 'audit-cache changed stage '+k)
                result['timings'] = rows; save(status='timings_running')
            long, rr = native('continued32', 32, cp=Path(zero['checkpoint']), geometry='affine', old_mesh=True)
            result['timings'] = rows; result['continued32'] = dict(receipt=rr, cache=long['evaluation_reuse'])
            result['tests']['cache_state_history_stages_bitwise'] = True
            save(status='completed_timings')
    except Exception as e:
        save(status='gate_or_numerical_stop', failure=str(e)); raise
    finally:
        lock.stream.close()


if __name__ == '__main__':
    main()
