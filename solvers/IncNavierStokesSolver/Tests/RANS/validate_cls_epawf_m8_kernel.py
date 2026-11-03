#!/usr/bin/env python3
"""Linux independent M8 constitutive/source qualification and certificates.

Verifies new production kernels against dimensional indexed tensors and direct
ray arithmetic; mathematical momentum bounds remain separate from covariance,
outer turbulence, time-map and physical terrain qualification.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                               'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import cls_epawf_m8_certificate as certificate
from cls_epawf_m8_oracle import indexed_closure, indexed_stress, stress_tangent, capped_source

FROZEN = {
    'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex':
        'd35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f',
    'docs/plans/cls-epawf-revised-contract-v1.json':
        '4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
    'library/SolverUtils/Turbulence/RANSModels.cpp':
        '86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
    'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp':
        '5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, condition, name):
        self.checks += 1
        if not condition:
            raise RuntimeError(name)

    def near(self, actual, expected, name, tolerance=1e-11, scale=1.):
        a, b = np.asarray(actual), np.asarray(expected)
        self.require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'Invalid '+name)
        value = float(np.max(np.abs(a-b), initial=0.))/max(scale, float(np.max(np.abs(b), initial=0.)))
        self.errors[name] = max(self.errors.get(name, 0.), value)
        self.require(value <= tolerance, f'{name}: {value} > {tolerance}')
        return value


def closure_arguments(case, step=0., imaginary=False):
    factor = 1j*step if imaginary else step
    return dict(gradient=np.asarray(case['gradient'])+factor*np.asarray(case['direction']),
                k=case['k']+factor*case['direction_k'],
                epsilon=case['epsilon']+factor*case['direction_epsilon'],
                beta=case['beta'], gamma=case['gamma'],
                reference_strain=np.asarray(case['reference_strain']),
                branch=case['branch'], distance=case['distance'], nu=case['nu'])


def check_closure_value(actual, expected, audit, prefix):
    for key, value in expected.items():
        audit.near(actual[key], np.real(value), prefix+'_'+key)


def closure_cases(report, audit):
    results = []
    for case in report['closure_cases']:
        arguments = closure_arguments(case)
        expected = indexed_closure(**arguments)
        check_closure_value(case['modified'], expected, audit, 'indexed_modified')
        original_arguments = dict(arguments, beta=1., gamma=0.)
        check_closure_value(case['original'], indexed_closure(**original_arguments), audit, 'indexed_original')
        audit.near(expected['production'], -np.sum(expected['stress']*arguments['gradient']), 'production_contraction')
        if case['beta'] == 1 and case['gamma'] == 0:
            audit.require(case['original'] == case['modified'], 'Original wrapper identity must be exact')
        if case['name'].startswith('reference_'):
            for key in expected:
                audit.near(case['modified'][key], case['original'][key], 'neutral_value_preserved_'+key)
        for sample in case['finite_differences']:
            for sign, side in ((1, 'plus'), (-1, 'minus')):
                child = indexed_closure(**closure_arguments(case, sign*sample['step']))
                check_closure_value(sample[side], child, audit, 'indexed_shifted')

        # Central differentiation of the physical max law at a nonsmooth tie
        # averages its two limiting directional actions. Formal branches are
        # differentiated separately, and both receive an interval certificate.
        step = 1e-30
        analytic = indexed_closure(**closure_arguments(case, step, imaginary=True))
        action = {key: np.imag(value)/step for key, value in analytic.items()}
        g = np.asarray(case['gradient'])
        direction = np.asarray(case['direction'])
        strain, rotation = g+g.T, g-g.T
        tie = abs(np.sum(strain*strain)-np.sum(rotation*rotation)) < 1e-12
        physical_tie = case['branch'] == 'physical' and tie
        if physical_tie:
            alternatives = []
            for branch in ('strain', 'rotation'):
                args = closure_arguments(case, step, imaginary=True)
                args['branch'] = branch
                alternatives.append({key: np.imag(value)/step for key, value in indexed_closure(**args).items()})
            # The active physical branch on +/- the perturbation can be
            # determined independently at a resolved small real step.
            selected = []
            for sign in (1, -1):
                child = g+sign*1e-7*direction
                s, o = child+child.T, child-child.T
                selected.append(0 if np.sum(s*s) >= np.sum(o*o) else 1)
            action = {key: .5*(alternatives[selected[0]][key]+alternatives[selected[1]][key]) for key in action}
        best = {}
        for sample in case['finite_differences']:
            h = sample['step']
            for key in expected:
                actual = (np.asarray(sample['plus'][key])-np.asarray(sample['minus'][key]))/(2*h)
                error = float(np.max(np.abs(actual-action[key]), initial=0.))
                magnitude = max(1., float(np.max(np.abs(action[key]), initial=0.)))
                rounding = 32*np.finfo(float).eps*max(1., float(np.max(np.abs(case['modified'][key]), initial=0.)))/h
                # O(h) at a physical invariant tie, O(h^2) on formal/smooth branches.
                bound = (100*h if physical_tie else 100*h*h)*magnitude+rounding+1e-8*magnitude
                audit.require(error <= bound, f'Native constitutive directional action {case["name"]}/{key}: {error}>{bound}')
                best[key] = min(best.get(key, float('inf')), error/magnitude)
        audit.require(max(best.values()) < (2e-5 if physical_tie else 2e-7), 'No resolved native constitutive directional action')
        principal = None
        if case['direction_k'] == 0 and case['direction_epsilon'] == 0:
            principal = float(np.sum(direction*-action['stress']))
            original_args = closure_arguments(case, step, imaginary=True)
            original_args['gamma'] = 0.
            nongamma = -indexed_closure(**original_args)['stress'].imag/step
            # This equality is direction-independent when explicitly using
            # H:(H+H^T), including the nongeneral rank-one q5 shear probe.
            actual_gamma = -indexed_closure(**closure_arguments(case, step, imaginary=True))['stress'].imag/step
            audit.near(actual_gamma-nongamma, case['gamma']*(direction+direction.T), 'gamma_full_tangent_shift', tolerance=1e-10)
        if case['name'] == 'off_reference_q5_shear':
            audit.require(principal < 0, 'Common beta must retain the adverse q=5 simple-shear control')
        results.append(dict(name=case['name'], branch=case['branch'], beta=case['beta'], gamma=case['gamma'],
                            principal_form=principal, minimum_relative_direction_errors=best,
                            covariance_minimum=float(np.linalg.eigvalsh(expected['stress'].real).min())))
    # All three common factors, both formal branches, the physical branch and
    # a general 3D state are required; this is not just a one-mode check.
    audit.require({case['beta'] for case in report['closure_cases']} == {1., .4, .45}, 'Missing cubic factor')
    audit.require({case['branch'] for case in report['closure_cases']} == {'physical', 'strain', 'rotation'}, 'Missing invariant branch')
    audit.require(any(case['name'] == 'general_3d' for case in report['closure_cases']), 'Missing general3D kernel')
    for row in results:
        if row['name'] == 'reference_transverse_x':
            expected = 1-20*row['beta']*report['reference']['cMu']+row['gamma']
            audit.near(row['principal_form'], expected, 'transverse_diffusion_coefficient', tolerance=1e-10)
    return results


def shifted_source(data, direction, step):
    result = dict(sampled_k=data['sampled_k']+step*direction['sampled_k'])
    for key in ('k', 'pre_e', 'pre_d', 'metric_weights'):
        raw = direction.get(key, [])
        values = np.asarray(data[key], dtype=float)
        result[key] = (values+step*(np.asarray(raw) if len(raw) else np.zeros_like(values))).tolist()
    return result


def source_cases(report, audit):
    results = []
    for case in report['source_cases']:
        expected, action = capped_source(case['input'], case['controls'], case['direction'])
        for key, value in expected.items():
            if key == 'capped':
                audit.require(case['value'][key] == value, 'Incorrect depletion regime')
                continue
            audit.near(case['value'][key], value, 'source_value_'+key, tolerance=2e-11)
            audit.near(case['directional'][key], action[key], 'source_direction_'+key, tolerance=2e-10)
        best = {}
        for sample in case['finite_differences']:
            h = sample['step']
            for sign, side in ((1, 'plus'), (-1, 'minus')):
                if side not in sample:
                    continue
                child, _ = capped_source(shifted_source(case['input'], case['direction'], sign*h), case['controls'])
                for key, value in child.items():
                    if key != 'capped':
                        audit.near(sample[side][key], value, 'source_shifted_'+key, tolerance=2e-11)
            for key in expected:
                if key == 'capped':
                    continue
                if case['one_sided']:
                    observed = (np.asarray(sample['plus'][key])-np.asarray(case['value'][key]))/h
                else:
                    observed = (np.asarray(sample['plus'][key])-np.asarray(sample['minus'][key]))/(2*h)
                magnitude = max(1., float(np.max(np.abs(action[key]), initial=0.)))
                error = float(np.max(np.abs(observed-action[key]), initial=0.))/magnitude
                rounding = 64*np.finfo(float).eps*max(1., float(np.max(np.abs(expected[key]), initial=0.)))/(h*magnitude)
                bound = (10*h if case['one_sided'] else 20*h*h)+rounding+2e-8
                audit.require(error <= bound, f'Native capped-source FD {case["name"]}/{key}: {error}>{bound}')
                best[key] = min(best.get(key, float('inf')), error)
        audit.require(max(best.values()) < (2e-5 if case['one_sided'] else 2e-7), 'No resolved capped-source action')
        k = np.asarray(case['input']['k'])
        destruction = np.asarray(case['value']['d'])
        audit.require(np.all(destruction[k == 0] == 0), 'Destruction at a local zero must vanish')
        audit.require(np.all(np.asarray(case['value']['rate']) >= 0), 'Negative destruction rate')
        audit.require(np.all(np.asarray(case['value']['rate']) <= np.asarray(case['value']['rate_bound'])*(1+1e-12)), 'Rate cap violated')
        audit.near(case['value']['destruction_budget'], case['value']['lambda']*case['value']['budget_d'], 'normalized_destruction_budget')
        audit.near(case['value']['dissipation_budget'], case['value']['lambda']*case['value']['budget_e'], 'dependent_dissipation_budget')
        if case['name'].startswith('plateau'):
            audit.require(case['value']['e'] == case['input']['pre_e'] and
                          case['value']['d'] == case['input']['pre_d'], 'Exact inactive source plateau lost')
        if case['name'].startswith('whole_zero'):
            audit.near(case['value']['lambda'], 0, 'fully_depleted_budget_attenuation')
        results.append(dict(name=case['name'], b=case['value']['b'], lambda_pos=case['value']['lambda'],
                            capped=case['value']['capped'], source_direction_errors=best))
    audit.require(len(report['source_cases']) == 18, 'Missing source limit/knot cases')
    for suffix in ('linear', 'smooth'):
        audit.require(all(any(case['name'] == name+'_'+suffix for case in report['source_cases'])
            for name in ('plateau', 'active_uncapped', 'capped_depletion', 'local_zero', 'whole_zero',
                         'limiter_kink_up', 'limiter_kink_down', 'cap_kink_up', 'cap_kink_down')), 'Missing source regime')
    return results


def mathematical_certificates(report, audit, include_envelopes):
    output = dict(reference=[], finite_boxes=[], analytic_envelopes=[], floating=certificate.floating_checks())
    q = report['reference']['q']
    for factor in (Fraction(2, 5), Fraction(9, 20)):
        reference, intervals = certificate.reference_certificate(factor)
        output['reference'].append(reference)
        audit.require(float(reference['q_bracket'][0]) <= q <= float(reference['q_bracket'][1]), 'Native reference outside certified interval')
        g = np.zeros((3, 3))
        g[0, 2] = 1.
        for branch in ('strain', 'rotation'):
            tangent = stress_tangent(g, q, 1., float(factor), branch=branch)
            middle = np.array([[(float(value.lo)+float(value.hi))/2 for value in row] for row in intervals[branch]])
            audit.near(tangent, middle, 'reference_interval_tangent_'+branch, tolerance=1e-12)
            audit.require(all(float(value[0]) > 0 for value in reference['branches'][branch]['LDL_pivot_intervals']), 'Nonpositive reference certificate pivot')
        if include_envelopes:
            for radius in ('.01', '.05', '.1', '.25'):
                box = certificate.finite_box_certificate(factor, radius, '.5', '.5')
                output['finite_boxes'].append(box)
                audit.require(float(box['required_gamma_over_nut0']) >= 0, 'Invalid required gamma')
                # Independent complex-step probes test the implementation of
                # the interval extension; enclosure of every state is proved
                # by directed interval operations, not by these samples.
                rng = np.random.default_rng(4109)
                for _ in range(6):
                    a = g+float(radius)*rng.uniform(-1, 1, (3, 3))
                    k = q*rng.uniform(.5, 1.5)
                    ep = rng.uniform(.5, 1.5)
                    for branch in ('strain', 'rotation'):
                        numeric = stress_tangent(a, k, ep, float(factor), branch=branch)
                        iv = np.asarray(box['branches'][branch]['tangent_intervals'], dtype=float)
                        audit.require(np.all(numeric >= iv[:, :, 0]-1e-10) and np.all(numeric <= iv[:, :, 1]+1e-10), 'Independent tangent outside interval enclosure')
                witness = box['regularized_covariance_witness']
                if witness['certified_negative']:
                    a = np.zeros((3, 3))
                    a[0, 2] = 1+float(radius)
                    stress = indexed_stress(a, q, 1., float(factor), float(box['required_gamma_over_nut0']), g+g.T).real
                    variance = (stress[0, 0]+2*stress[0, 2]+stress[2, 2])/2
                    audit.require(variance < 0, 'Certified covariance witness disagrees with independent tensor')
                    audit.near(variance, sum(map(float, witness['variance_interval']))/2,
                               'regularized_covariance_witness', tolerance=1e-10)
        output['analytic_envelopes'].append(certificate.analytic_envelope_certificate(beta=factor))
    audit.near(output['floating']['q_critical'], 3.6332736350, 'off_reference_shear_threshold', tolerance=2e-10)
    audit.require(output['floating']['factor_at_q5'] < 0, 'Missing off-reference adverse shear factor')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--envelopes', action='store_true')
    parser.add_argument('--report', type=Path, help='Replay this existing probe JSON instead of executing the probe')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Numerical qualification must run on authoritative Linux')
    args.output.mkdir(parents=True, exist_ok=False)
    source = next(path for path in Path(__file__).resolve().parents if (path/'.git').exists())
    def git(*argv):
        return subprocess.check_output(['git', '-C', str(source), *argv], text=True).strip()
    audit = Audit()
    revision = git('rev-parse', 'HEAD')
    manifest = dict(schema=1, milestone='M8 independent kernel and constitutive qualification',
                    source_commit=revision, source_status=git('status', '--porcelain'),
                    started_utc=datetime.now(timezone.utc).isoformat(), threads=THREADS,
                    probe=str(args.probe.resolve()), probe_sha256=digest(args.probe),
                    sources={str(Path(__file__).relative_to(source)): digest(__file__),
                             'oracle': digest(Path(__file__).with_name('cls_epawf_m8_oracle.py')),
                             'certificate': digest(Path(__file__).with_name('cls_epawf_m8_certificate.py'))},
                    status='running', commands=[], numerical_robustness='not implied by these isolated gates')
    path = args.output/'validation.json'
    def save():
        manifest.update(checks=audit.checks, errors=audit.errors)
        path.write_text(json.dumps(manifest, indent=2, default=json_value, allow_nan=False)+'\n')
    save()
    try:
        audit.require(not manifest['source_status'], 'Source worktree must be clean')
        for name, expected in FROZEN.items():
            audit.require(digest(source/name) == expected, 'Frozen original changed: '+name)
        native = args.report.resolve() if args.report else args.output/'kernel.json'
        if args.report is None:
            command = [str(args.probe.resolve()), '--output', str(native.resolve())]
            start = time.monotonic()
            with (args.output/'kernel.log').open('w') as stream:
                result = subprocess.run(command, cwd=args.output, env=dict(os.environ, **THREADS),
                                        stdout=stream, stderr=subprocess.STDOUT, timeout=300, check=False)
            manifest['commands'].append(dict(argv=command, exit_code=result.returncode, seconds=time.monotonic()-start,
                                              log_sha256=digest(args.output/'kernel.log')))
            audit.require(result.returncode == 0, 'Robust variant probe failed; see kernel.log')
        report = json.loads(native.read_text())
        manifest.update(native_report=str(native), native_sha256=digest(native))
        audit.require(report['success'], 'Native probe reported failure')
        manifest['closure_cases'] = closure_cases(report, audit)
        manifest['source_cases'] = source_cases(report, audit)
        names = {item['name'] for item in report['rejections']}
        expected_names = {'negative_k', 'zero_K', 'negative_E', 'negative_D', 'zero_weight', 'bad_budget',
                          'bad_shape', 'nonfinite', 'theta_zero', 'theta_one', 'delta_zero', 'delta_one',
                          'negative_beta', 'negative_gamma', 'asymmetric_Sref', 'traceful_Sref',
                          'inconsistent_point_denominator', 'outward_zero_k_direction'}
        audit.require(names == expected_names, 'Missing required invalid-input rejections')
        certificates = mathematical_certificates(report, audit, args.envelopes)
        certificate_path = args.output/'certificates.json'
        certificate_path.write_text(json.dumps(certificates, indent=2, allow_nan=False)+'\n')
        manifest.update(certificate_file=str(certificate_path), certificate_sha256=digest(certificate_path),
                        native_checks=report['checks'], status='passed',
                        gates=dict(reference_constitutive='certified locally for beta=.4/.45',
                                   capped_source='verified identities and one-sided derivatives; continuum source quasi-positivity condition',
                                   common_beta_off_reference='failed at q=5 as required control',
                                   momentum_envelope='see explicit mathematical bounds and required fixed gamma',
                                   full_covariance_and_outer_positivity='separate unqualified gates',
                                   terrain='not qualified by kernel checks'))
        audit.require(git('rev-parse', 'HEAD') == revision and not git('status', '--porcelain'), 'Source changed during qualification')
        audit.require(digest(args.probe) == manifest['probe_sha256'], 'Probe changed during qualification')
    except Exception as error:
        manifest.update(status='failed', error=str(error))
        raise
    finally:
        manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save()
    print(json.dumps(dict(status=manifest['status'], checks=audit.checks, source_commit=revision,
                          receipt=str(path)), indent=2))


if __name__ == '__main__':
    main()
