#!/usr/bin/env python3
"""Curate retained R3 C/D receipts into compact JSON and Markdown tables.

Read-only postprocessing. No imports of a closure, assembly, solver or numeric
oracle; no native executions. Source/report/vector hashes are checked before
extracting a successful analysis. Failed receipts remain explicitly failed.
Run on Linux with one worker, as with the other numerical evidence tooling.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import hashlib
import json
from pathlib import Path
import platform
import struct
import sys

sys.dont_write_bytecode = True
BLOCKS = ['u_x', 'u_z', 'k', 'epsilon', 'pi_inner', 'pi_outer', 'gauge', 'constraint']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pointer(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def vector_sha(values):
    digest = hashlib.sha256()
    for value in values:
        digest.update(struct.pack('<d', value))
    return digest.hexdigest()


def maximum(values):
    return max(values, default=0.)


def norm(record):
    return dict(blocks=dict(zip(BLOCKS, record['blocks'])), maximum=record['maximum'])


def extract(label, native_path, analysis_path):
    native, analysis = json.loads(Path(native_path).read_text()), json.loads(Path(analysis_path).read_text())
    result = dict(label=label, native=pointer(native_path), analysis=pointer(analysis_path),
        analysis_success=bool(analysis.get('success', False)), checks=analysis.get('checks'),
        producing_revision=analysis.get('producing_revision'),
        analysis_revision=analysis.get('source_revision'),
        runtime_receipt=analysis.get('native_runtime_receipt'),
        inner_stress_model=native['inner_stress_model'], time=native['time'], alpha=native['alpha'],
        history_matches_rebuilt_h0=native['packet_history_matches_rebuilt_h0'],
        constructed=bool(native['construct_contact']),
        state_sha256=vector_sha(native['state']), original_history_sha256=vector_sha(native['history']),
        native_extended_precision_bits=native['extended_precision_bits'],
        failures=analysis.get('failures', []))
    if analysis['native_report']['sha256'] != result['native']['sha256']:
        raise ValueError(label+': analysis is not bound to these native report bytes')
    if not result['analysis_success']:
        result['interpretation'] = 'Independent analysis failed; no independent term, norm or contact gate is promoted'
        return result
    if analysis['state_sha256'] != result['state_sha256'] or analysis['original_history_sha256'] != result['original_history_sha256']:
        raise ValueError(label+': independent/native coefficient hashes differ')
    if (analysis['time'], analysis['alpha'], analysis['inner_stress_model']) != (
            result['time'], result['alpha'], result['inner_stress_model']):
        raise ValueError(label+': independent/native target identity differs')
    result['target_tolerance'] = analysis['exact_target_tolerance']
    result['observed_norms'] = {key: norm(value) for key, value in analysis['observed_norms'].items()
        if 'blocks' in value}
    result['observed_time_mass_difference'] = analysis['observed_norms']['native_time_minus_independent_time']
    result['fixed_vector_exact_polynomial_mass_norms'] = {
        key: dict(maximum_interval_decimal=value['maximum_interval_decimal'],
            below_original_2e_minus12_tolerance=value['below_original_2e_minus12_tolerance'],
            block_intervals_decimal=[row['interval_decimal'] for row in value['blocks']],
            scope='Fixed declared rounded residual coefficients and exact polynomial physical mass; excludes residual-input/native stored-mass errors')
        for key, value in analysis['certified_exact_mass_norms'].items()}
    result['native_arithmetic_norms'] = {key: dict(blocks=dict(zip(BLOCKS, value)), maximum=maximum(value))
        for key, value in native['norms'].items() if not key.startswith('term_')}
    terms = {}
    for key, value in analysis['term_comparisons'].items():
        terms[key] = dict(
            independent_free_physical_norms=norm(value['independent']['free_physical_norms']),
            native_free_physical_norms=norm(value['native']['free_physical_norms']),
            difference_free_physical_norms=norm(value['difference']['free_physical_norms']),
            maximum_coefficient_error=value['maximum_coefficient_error'],
            maximum_strong_reaction_coefficient=value['native']['maximum_strong_reaction_coefficient'])
    result['term_by_block'] = terms
    result['global_accumulation_bound'] = analysis['frozen_accumulation_bounds']['per_key']['legacy_spatial']
    result['accumulation_bound_scope'] = analysis['frozen_accumulation_bounds']['scope']
    result['time_mass_comparisons'] = analysis['exact_polynomial_time_mass_comparisons']
    trace = analysis['independent_trace']
    result['trace'] = {key: trace[key] for key in ('relative_to_E0', 'matching_rows',
        'actual_epsilon_ratio_range', 'chi_zero_epsilon_ratio_range', 'baseline_and_current_lambda_range')}
    result['native_domain'] = {key: native['scalar_domain'][key] for key in ('admissible',
        'phase', 'minimum_dimensionless_margin', 'minimum_inner_k', 'minimum_outer_k',
        'minimum_outer_epsilon', 'minimum_sampled_K', 'minimum_projected_trace_epsilon',
        'minimum_b', 'minimum_lambda', 'inventory_counts')}
    result['physical_term_metadata'] = analysis['physical_term_metadata']
    result['constructed_state_proof'] = analysis['constructed_state_proof']
    ranges = analysis.get('constructed_finite_ranges')
    if ranges is not None:
        result['constructed_finite_observations'] = {key: ranges[key] for key in (
            'finite_scalar_admissible', 'finite_covariance_nonnegative', 'fixed_outer_box_covered',
            'lambda_guard_satisfied', 'finite_joint_preflight_passed', 'fixed_outer_box_margins',
            'minimum_lambda', 'covariance_inner_minimum', 'covariance_outer_minimum', 'scope')}
    contact = analysis['contact']
    result['contact'] = contact
    result['fixed_packet_comparison'] = analysis.get('source_packet')
    result['scientific_limits'] = analysis['scientific_limits']
    result['interpretation'] = ('Implementation/arithmetic and conditional contact diagnostics; '
        'not an implicit root, an accepted time sample, a roundoff floor, or a PDE/trajectory positivity certificate')
    return result


def number(value):
    if value is None:
        return 'unavailable'
    return f'{value:.9g}'


def interval(value):
    # Compact outward rounding for display only. Full original certificate
    # endpoints remain in the companion JSON and their source receipts.
    with localcontext() as context:
        context.prec = 10
        context.rounding = ROUND_FLOOR
        lower = +Decimal(value[0])
        context.rounding = ROUND_CEILING
        upper = +Decimal(value[1])
    return '['+str(lower)+', '+str(upper)+']'


def markdown(records):
    lines = ['# M8-R3 fixed-state residual resolution and scalar contact', '',
        'This is a compact extraction of retained native and independent receipts. '
        'The original coupled residual tolerance remains2e-12. A passing analysis '
        'does not turn a saved failed candidate into a root.', '',
        '| Case | Native maximum | Independent maximum | Native−independent | '
        'Native−extended | Native exact-polynomial-mass interval | Analysis |',
        '|---|---:|---:|---:|---:|---|---|']
    for row in records:
        if not row['analysis_success']:
            lines.append(f"| {row['label']} | — | — | — | — | — | FAILED |")
            continue
        observed = row['observed_norms']
        certificate = row['fixed_vector_exact_polynomial_mass_norms']['native']
        lines.append('| '+' | '.join((row['label'], number(observed['native']['maximum']),
            number(observed['independent_terms']['maximum']), number(observed['native_minus_independent']['maximum']),
            number(observed['native_minus_extended']['maximum']), interval(certificate['maximum_interval_decimal']),
            'passed implementation audit'))+' |')
    lines += ['', 'The interval column fixes the saved rounded residual coefficients and exact '
        'polynomial mass. It excludes errors in those coefficients and differences in native '
        'stored mass entries. The observed differences and frozen-increment summation bounds '
        'do not establish a roundoff floor.', '']
    for row in records:
        lines += ['## '+row['label'], '',
            f"Native source `{row['producing_revision']}`; analysis source `{row['analysis_revision']}`. "
            f"Model {row['inner_stress_model']}, alpha={number(row['alpha'])}, time={number(row['time'])}. "
            f"Original h0 SHA256 `{row['original_history_sha256']}`.", '']
        if not row['analysis_success']:
            lines += ['Independent analysis failed. '+json.dumps(row['failures'], ensure_ascii=False), '']
            continue
        lines += ['| Residual evaluation | u_x | u_z | k | epsilon | pi_inner | pi_outer | gauge | constraint |',
            '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for name, value in row['observed_norms'].items():
            lines.append('| '+name+' | '+' | '.join(number(value['blocks'][key]) for key in BLOCKS)+' |')
        lines += ['', 'Term norms below use only free physical rows. Strong-row physical loads '
            'are reactions and remain in JSON separately; they are not constraint defects.', '',
            '| Physical term | u_x | u_z | k | epsilon | pi_inner | pi_outer | max native−independent |',
            '|---|---:|---:|---:|---:|---:|---:|---:|']
        for name, value in row['term_by_block'].items():
            blocks = value['independent_free_physical_norms']['blocks']
            lines.append('| '+name+' | '+' | '.join(number(blocks[key]) for key in BLOCKS[:6])+ ' | '+
                number(value['difference_free_physical_norms']['maximum'])+' |')
        lines += ['', 'Net destruction enters once as +D in the k residual. The diagnostic E and '
            '−C loads are auxiliary checks and are absent from the physical sum.', '',
            'Projected trace coefficient mismatch / E0: '+number(row['trace']['relative_to_E0'])+'. '
            'This is the actual retained-state mismatch, separate from the projected target definition.', '']
        contact = row['contact']
        if contact.get('requested'):
            lines += ['| Scalar / location | State value | Native physical rate | Independent physical rate | '
                'Fixed-load exact-mass rate interval | Largest absolute term |',
                '|---|---:|---:|---:|---|---|']
            for scalar, data in contact['scalar_mass_rates'].items():
                for point in data['points']:
                    certificate = point.get('fixed_native_load_exact_mass_rate')
                    location = ','.join(number(x) for x in point['coordinate'])
                    term = point['largest_absolute_component']
                    lines.append('| '+' | '.join((scalar+' @('+location+')', number(point['native_value']),
                        number(point['native_rate']), number(point['independent_rate']),
                        interval(certificate['interval_decimal']) if certificate else 'not assessed',
                        term+'='+number(point['component_rates'][term])))+' |')
            lines += ['', contact['rate_scope'], '', contact['pressure_scope']+'.', '',
                'The original fixed contact x=1.6046496089512394 and the constructed zero '
                'x=1 are distinct. A positive state value at the original row is not a '
                'zero-contact inward-invariance test. The saved terminal near8e-15 with '
                'unresolved sign and nonzero trace mismatch is not a certified counterexample.', '']
        if row['constructed']:
            observations = row['constructed_finite_observations']
            lines += ['Constructed state: exact coefficient identity / symbolic k nonnegativity = '+
                str(row['constructed_state_proof']['coefficient_identity_exact'])+' / '+
                str(row['constructed_state_proof']['analytic_k_nonnegative'])+'. '
                'Finite scalar / covariance / fixed outer box / lambda guard = '+
                ' / '.join(str(observations[key]) for key in ('finite_scalar_admissible',
                    'finite_covariance_nonnegative', 'fixed_outer_box_covered', 'lambda_guard_satisfied'))+'.', '']
        lines += ['Native receipt SHA256 `'+row['native']['sha256']+'`; independent receipt SHA256 `'+
            row['analysis']['sha256']+'`.', '']
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', nargs=3, action='append', required=True,
        metavar=('LABEL', 'NATIVE_JSON', 'ANALYSIS_JSON'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--markdown', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Evidence extraction is run on Linux')
    if args.output.exists() or args.markdown.exists():
        raise SystemExit('Refusing to overwrite curated evidence')
    labels = [row[0] for row in args.case]
    if len(labels) != len(set(labels)):
        raise SystemExit('Repeated case labels')
    records = [extract(*row) for row in args.case]
    result = dict(schema=1, phase='M8-R3 C/D compact extraction',
        generated_utc=datetime.now(timezone.utc).isoformat(), extractor=pointer(__file__),
        cases=records, all_independent_analyses_passed=all(row['analysis_success'] for row in records),
        scoped_diagnostics_only=True, numerical_or_native_execution=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    args.markdown.write_text(markdown(records))
    print(json.dumps(dict(cases=len(records), output=str(args.output), markdown=str(args.markdown),
        all_independent_analyses_passed=result['all_independent_analyses_passed'])))


if __name__ == '__main__':
    main()
