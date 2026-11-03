"""Read-only component diagnosis of the retained initial outer-k J0 action.

The original failed absolute gate is preserved. Exact affine roughness forcing,
complex-k indexed outer tensors and affine scalar transport are compared with
the stored matrix and the already recorded whole-oracle differences. No native
calls, field changes, coefficient search or nonlinear-policy changes occur.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
INITIAL_REPORT_SHA = '52491ec3974757cf4b177b0502360171e5715448c68ab1ceb19b8bf564e38977'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-c', 'pack.threads=1',
        '-c', 'checkout.workers=1', '-C', str(source), *arguments], text=True).strip()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'initial-report', 'initial-receipt', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-commit', required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Scientific component diagnosis belongs to Linux')
    args.source = args.source.resolve(); args.output.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema=1, success=False, started_utc=stamp(), native_calls=0,
        source_commit=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(args.source, 'status', '--porcelain')),
        cpu_affinity=sorted(os.sched_getaffinity(0)), numerical_threads={key:os.environ.get(key) for key in THREADS},
        driver_sha256=sha(__file__), command_argv=sys.argv, inputs={}, failures=[],
        solved_path_assessed=False, target_stage_status_assessed=False, M8_closed=False)
    for path in (args.initial_report, args.initial_receipt, Path(__file__)):
        receipt['inputs'][str(path.resolve())] = dict(bytes=path.stat().st_size, sha256=sha(path))
    save(args.output/'receipt.json', receipt); started = time.monotonic()
    try:
        if receipt['source_commit'] != args.expected_commit or not receipt['source_clean_start']:
            raise ValueError('Unpinned analysis source')
        if receipt['cpu_affinity'] != [48] or any(receipt['numerical_threads'][key] != '1' for key in THREADS):
            raise ValueError('Diagnosis requires CPU48 and six numerical thread settings equal to one')
        if sha(args.initial_report) != INITIAL_REPORT_SHA:
            raise ValueError('This bounded diagnosis requires the exact retained failed initial-J0 report')
        previous = json.loads(args.initial_receipt.read_text())
        if not previous['success'] or previous['producing_source_commit'] != '9f20b3053c303e8805f5183f668a8ec468d9656b':
            raise ValueError('Original initial-J0 analysis/provenance did not complete')
        for path, item in previous['inputs'].items():
            if sha(path) != item['sha256'] or Path(path).stat().st_size != item['bytes']:
                raise ValueError('Original bound input changed: '+path)
            receipt['inputs'][path] = item
        if sha(previous['producing_solver']['path']) != previous['producing_solver']['sha256']:
            raise ValueError('Producing native solver changed')
        argv = previous['command_argv']
        option = lambda name:Path(argv[argv.index('--'+name)+1])
        sys.dont_write_bytecode = True
        tests = args.source/'solvers/IncNavierStokesSolver/Tests/RANS'; sys.path.insert(0, str(tests))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as stage
        from cls_epawf_m8_r2_oracle import indexed_outer, selected_outer
        from cls_epawf_m8_r1_preflight import add
        from cls_epawf_m8_r2_preflight import Inventory
        from cls_epawf_m7_energy_oracle import _horizontal, _velocity, _sample
        from cls_epawf_m6_operator_oracle import _diagnostic
        from validate_cls_epawf_m5 import physical_rows
        imported = {str(Path(module.__file__).resolve()):sha(module.__file__)
            for module in list(sys.modules.values()) if getattr(module, '__file__', None)
            and str(Path(module.__file__).resolve()).startswith(str(args.source))}
        receipt['imported_source_files_start'] = imported
        old = json.loads(args.initial_report.read_text()); raw = json.loads(option('preflight').read_text())
        analyzer = stage.Analyzer(raw, option('stage-session'), args.source, option('matrix').parent)
        state = np.asarray(raw['initial']); direction = np.asarray(old['outer_k_action']['direction'])
        d, off, H, nu = raw['dimension'], raw['offsets'], raw['parameters'][1], raw['parameters'][0]
        if raw['inner_stress_model'] != 'Linear' or raw['model_variant'] != 'Regularized':
            raise ValueError('This derivation is restricted to the retained Linear/Regularized initial state')
        if stage.vector_hash(state) != old['initial_h0_float64_LE_sha256']:
            raise ValueError('Initial coefficients changed')
        if np.any(direction[:off[d]]) or np.any(direction[off[d+1]:]):
            raise ValueError('Stored direction is not pure transported k')
        matrix = analyzer.matrix(dict(jacobian_binary_file=str(option('matrix')), jacobian_binary_shape=[2165, 2165]))
        current = matrix@direction
        analyzer.audit.near(current, old['outer_k_action']['actual_current_J_action'], 'same_retained_J0_and_k_direction', 0.)
        constraints = analyzer.constraints(state, raw['reference_constraints'])
        inventory = Inventory(raw, analyzer.fields, state, option('stage-session'))
        nonzero = {int(index) for index in np.flatnonzero(direction[off[d]:off[d+1]])}
        sample_height = raw['parameters'][3]
        sample_trace = analyzer.fields[d].boundary(lambda axis, side, element:
            axis == 2 and side == 1 and abs(element['upper'][2]-sample_height) < 1e-12)
        if not sample_trace or nonzero.intersection(sample_trace):
            raise ValueError('Physical K-sampler trace has no structural zero direction proof')
        constant = analyzer.fields[d].constant()
        pivot = int(np.flatnonzero(constant)[0]); Kconstant = state[off[d]+pivot]/constant[pivot]
        constant_error = float(np.max(np.abs(state[off[d]:off[d+1]]-Kconstant*constant)))
        if any(element['degree'] != 3 for element in analyzer.fields[d].elements):
            raise ValueError('Initial limiter proof is restricted to the declared p3 Modified_A tensor basis')
        # Every p3 Modified_A tensor factor has absolute value at most one.
        # This coefficient bound checks the entire represented initial field,
        # including any last-bit constant-projection coefficients; it changes
        # no field value and is not a positivity repair.
        deviation = state[off[d]:off[d+1]]-Kconstant*constant
        deviation_bound = max(float(np.sum(np.abs(deviation[np.asarray(element['ids'], int)])))
            for element in analyzer.fields[d].elements)*(1+64*np.finfo(float).eps)
        limiter_margin = (Kconstant-deviation_bound)/(raw['robust_parameters'][2]*(Kconstant+deviation_bound))-1
        if limiter_margin <= 0.:
            raise ValueError('Initial limiter plateau has no positive coefficient-based margin')
        anchors = sorted({tuple(point) for point, _ in analyzer.source_anchors}.union(tuple(point) for point in analyzer.trace_anchors))
        view = dict(raw, transition=dict(controls=raw['transition_controls'], weights=raw['transition_weights'],
            reference=raw['transition_reference']))
        state_gate_margin, sample_action = float('inf'), 0.
        for anchor in anchors:
            point = np.asarray(anchor); diagnostic = _diagnostic(view, state, analyzer.fields, point, zero_body=True)
            if diagnostic['gates'][1] != 0.:
                state_gate_margin = min(state_gate_margin, raw['transition_controls'][5]**2-float(diagnostic['gates'][-1]))
            sample = [point[0], point[1], sample_height]
            value, gradient = analyzer.fields[d].evaluate(direction, off[d], sample)
            sample_action = max(sample_action, abs(value), abs(gradient[0]))
        if sample_action != 0. or state_gate_margin <= 0.:
            raise ValueError('Sampler or transition first-derivative remainder is not verified zero')
        assumptions = dict(pure_k_direction=True, physical_K_sampler_trace_DOF_overlap=0,
            physical_K_sampler_value_and_tangent_action=sample_action,
            constant_initial_k_coefficient_error=constant_error,
            initial_k_deviation_from_constant_coefficient_bound=deviation_bound,
            limiter_distance_above_plateau_lower_bound=limiter_margin,
            source_and_trace_state_gate_squared_distance_to_on=state_gate_margin,
            chi_action_zero_from_open_plateau_and_zero_sampler=True,
            normalization_lambda_action_zero_from_limiter_plateau=True,
            projected_epsilon_H_action_zero=True,
            gamma_action_zero='GammaScale*nu_t0 is prescribed in physical r, independent of current k',
            scope='Only the retained initial direction/state; no conclusion at later Newton iterates')
        # Exact affine counterpart of -F_A, where delta F_A=-2 grad_Gamma(delta k_H)/3.
        rough = np.zeros(analyzer.size); rough_rows = set(); outer_rows = set()
        for element in analyzer.fields[0].elements:
            if element['lower'][2] >= H-1e-12:
                outer_rows.update(off[b]+int(index) for b in range(d) for index in element['ids'])
            lo, hi = max(0., element['lower'][2]), min(raw['parameters'][2], element['upper'][2])
            if hi <= lo:
                continue
            z, w = np.polynomial.legendre.leggauss(8)
            points, weights = [], []
            for i in range(8):
                for j in range(8):
                    points.append([element['lower'][0]+(z[i]+1)*(element['upper'][0]-element['lower'][0])/2,
                        0., lo+(z[j]+1)*(hi-lo)/2])
                    weights.append(w[i]*w[j]*(element['upper'][0]-element['lower'][0])*(hi-lo)/4)
            points, weights = np.asarray(points), np.asarray(weights)
            ids, phi, gradients = analyzer.fields[0].rows(element, points)
            trace_points = points.copy(); trace_points[:, 2] = H
            derivative = np.array([analyzer.fields[d].evaluate(direction, off[d], point, 1)[1][0] for point in trace_points])
            add(rough, off[0], ids, phi, gradients, weights, 2*derivative/3, np.zeros((len(points), 3)))
            rough_rows.update(off[0]+int(index) for index in ids)
        rough_rows.difference_update(constraints); outer_rows.difference_update(constraints)
        if rough_rows.intersection(outer_rows):
            raise ValueError('Roughness and outer momentum support sets unexpectedly overlap')
        scalar_inner = np.zeros(analyzer.size)
        for chart, horizontal in _horizontal(analyzer.fields, d, H, int(raw['quadrature_orders'][0])):
            ray = inventory.ray(chart)
            groups = {}
            for index, point in enumerate(ray['points']):
                groups.setdefault(analyzer.fields[0].element(point)['geometry'], []).append(index)
            for identity, indices in groups.items():
                indices = np.asarray(indices); points = ray['points'][indices]
                ids, phi, gradient, delta, grad_delta = physical_rows(analyzer.fields[d], direction, off[d], points, identity)
                U, G = _velocity(analyzer.fields, state, off, d, points, identity)
                nt = float(ray['profile'].a)*np.maximum(ray['r'][indices]-float(ray['profile'].zv), 0.)
                weights = horizontal*ray['metric']['baseArea']*ray['measure'][indices]
                add(scalar_inner, off[d], ids, phi, gradient, weights,
                    2*delta*np.trace(G, axis1=1, axis2=2)/3, (nu+nt)[:, None]*grad_delta-U*delta[:, None])
        beta = raw['robust_parameters'][0]
        def outer_action(complex_step, branch):
            output = np.zeros(analyzer.size); minimum_yap_gap = float('inf')
            for element in analyzer.fields[d+3].elements:
                points, weights = analyzer.fields[d+3].quadrature(element, int(raw['quadrature_orders'][0]))
                identity = element['geometry']; U, G = _velocity(analyzer.fields, state, off, d, points, identity)
                ki, kp, kg, k, gk = physical_rows(analyzer.fields[d], state, off[d], points, identity)
                _, _, _, delta, gd = physical_rows(analyzer.fields[d], direction, off[d], points, identity)
                ei, ep, eg, epsilon, ge = physical_rows(analyzer.fields[d+1], state, off[d+1], points, identity)
                dT, dnt, dPk, dSe = [], [], [], []
                for index, point in enumerate(points):
                    complex_k = k[index]+1j*complex_step*delta[index]
                    value = indexed_outer(G[index], complex_k, epsilon[index], branch)
                    S = G[index]+G[index].T
                    derivative = (value['nut'].imag*S-value['quadratic'].imag-beta*value['cubic'].imag)/complex_step
                    nt_derivative = float(value['nut'].imag/complex_step)
                    production_derivative = float(np.sum(derivative*G[index])-2*delta[index]*np.trace(G[index])/3)
                    R, _ = selected_outer(raw, G[index], k[index], epsilon[index], point)
                    production = -float(np.sum(R*G[index]))
                    ratio = complex_k**1.5/(2.5*epsilon[index]*point[2])
                    minimum_yap_gap = min(minimum_yap_gap, abs(float(ratio.real)-1))
                    yap = (ratio-1)*ratio**2 if ratio.real > 1 else 0.
                    eps_source = (1.44*epsilon[index]/complex_k*(production+1j*complex_step*production_derivative)
                        -1.92*epsilon[index]**2/complex_k+.83*epsilon[index]**2/complex_k*yap)
                    dT.append(derivative); dnt.append(nt_derivative); dPk.append(production_derivative)
                    dSe.append(float(np.imag(eps_source)/complex_step))
                dT, dnt = np.asarray(dT), np.asarray(dnt)
                for block, axis in enumerate((0, 2)):
                    vi, vp, vg, _, _ = physical_rows(analyzer.fields[0], state, off[block], points, identity)
                    add(output, off[block], vi, vp, vg, weights, np.zeros(len(points)), dT[:, axis])
                nt = np.array([selected_outer(raw, g, kval, eps, point)[1] for g, kval, eps, point in zip(G, k, epsilon, points)])
                add(output, off[d], ki, kp, kg, weights, -np.asarray(dPk), dnt[:, None]*gk+nt[:, None]*gd-U*delta[:, None])
                add(output, off[d+1], ei, ep, eg, weights, -np.asarray(dSe), dnt[:, None]/1.3*ge)
            return output, minimum_yap_gap
        components = []
        mass = analyzer.time_mass(direction, np.zeros(analyzer.size), 3200., constraints)
        for branch in ('strain', 'rotation'):
            for step in (1e-20, 1e-30):
                outer, yap_gap = outer_action(step, branch)
                full = rough+scalar_inner+outer+mass
                for row in constraints:
                    full[row] = direction[row]  # projected target action was independently verified zero above.
                components.append(dict(branch=branch, complex_step=step, minimum_Yap_switch_gap=yap_gap,
                    outer_action=outer, complete_action=full))
        reference = components[0]['complete_action']
        controls = [dict(branch=row['branch'], complex_step=row['complex_step'],
            minimum_Yap_switch_gap=row['minimum_Yap_switch_gap'],
            action_difference_from_first=analyzer.audit.norm_record(row['complete_action']-reference, constraints)) for row in components]
        stable = all(row['action_difference_from_first']['maximum'] < 1e-10 and row['minimum_Yap_switch_gap'] > 1e-10 for row in controls)
        error = current-reference
        def masked(vector, rows):
            value = np.zeros_like(vector); indices = sorted(rows); value[indices] = vector[indices]; return value
        mask_records = dict(roughness_velocity_load_rows=sorted(rough_rows), outer_velocity_load_rows=sorted(outer_rows),
            disjoint=True, interpretation='Coefficient load-support masks, not orthogonal energy regions; their norms are not added')
        levels = []
        for level in old['outer_k_action']['levels']:
            h = level['step']
            if not all(level[name]['evaluated'] for name in ('plus', 'minus')):
                levels.append(dict(step=h, available=False)); continue
            difference = (np.asarray(level['plus']['spatial_residual'])-np.asarray(level['minus']['spatial_residual']))/(2*h)+mass
            fd_error = difference-reference
            closure = (current-difference)-(error-fd_error)
            analyzer.audit.near(closure, np.zeros(analyzer.size), 'J_minus_FD_exact_error_decomposition', 2e-10)
            levels.append(dict(step=h, available=True,
                retained_J_minus_recorded_FD=analyzer.audit.norm_record(current-difference, constraints),
                recorded_oracle_FD_minus_analytic=analyzer.audit.norm_record(fd_error, constraints),
                recorded_FD_roughness_support_error=analyzer.audit.norm_record(masked(fd_error, rough_rows), constraints),
                recorded_FD_outer_support_error=analyzer.audit.norm_record(masked(fd_error, outer_rows), constraints)))
        result = dict(schema=1, original_failed_absolute_gate_preserved=True, original_report_sha256=INITIAL_REPORT_SHA,
            assumptions=assumptions, analytic_controls=controls, analytic_controls_met=bool(stable),
            stored_J_minus_analytic=analyzer.audit.norm_record(error, constraints),
            stored_J_roughness_support_error=analyzer.audit.norm_record(masked(error, rough_rows), constraints),
            stored_J_outer_support_error=analyzer.audit.norm_record(masked(error, outer_rows), constraints),
            exact_roughness_action=rough.tolist(), affine_inner_scalar_action=scalar_inner.tolist(),
            indexed_outer_action=components[0]['outer_action'].tolist(), exact_mass_action=mass.tolist(),
            complete_initial_analytic_action=reference.tolist(), masks=mask_records, retained_FD_decomposition=levels,
            fixed_absolute_action_tolerance=1e-7,
            fixed_absolute_gate_met=bool(analyzer.audit.norm_record(error, constraints)['maximum'] <= 1e-7),
            no_interval_floating_point_error_certificate=True,
            interpretation='A stable analytic component comparison localizes differentiation/subtraction error numerically; it is not a rigorous floating-point proof, a nonlinear root test or a physical validation',
            solved_path_assessed=False, target_stage_status_assessed=False, M8_closed=False,
            checks=analyzer.audit.checks, errors=analyzer.audit.errors)
        result['localization_classification'] = ('unresolved_analytic_reference_controls' if not stable else
            'retained_J_fails_against_stable_analytic_reference' if not result['fixed_absolute_gate_met'] else
            'retained_J_matches_analytic_reference_original_FD_gate_still_preserved')
        result['rigorous_derivative_defect_proven'] = False
        result['native_intervention_or_retry_requested'] = False
        save(args.output/'validation.json', result)
        receipt['success'] = True; receipt['analytic_controls_met'] = bool(stable)
        receipt['fixed_absolute_gate_met'] = result['fixed_absolute_gate_met']
        for path, expected in imported.items():
            if sha(path) != expected:
                raise ValueError('Imported analysis source changed: '+path)
    except Exception as error:
        receipt['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        receipt['success'] = False
    for path, item in receipt['inputs'].items():
        if sha(path) != item['sha256'] or Path(path).stat().st_size != item['bytes']:
            receipt['success'] = False; receipt['failures'].append(dict(type='ProvenanceError', message='Bound input changed: '+path))
    receipt.update(source_commit_end=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_end=not bool(git(args.source, 'status', '--porcelain')), finished_utc=stamp(), seconds=time.monotonic()-started)
    if receipt['source_commit_end'] != receipt['source_commit'] or not receipt['source_clean_end']:
        receipt['success'] = False; receipt['failures'].append(dict(type='ProvenanceError', message='Analysis source changed'))
    save(args.output/'receipt.json', receipt)
    print(json.dumps({key:receipt.get(key) for key in ('success', 'analytic_controls_met', 'fixed_absolute_gate_met', 'seconds', 'failures')}))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
