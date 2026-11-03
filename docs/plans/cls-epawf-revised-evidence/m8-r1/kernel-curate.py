#!/usr/bin/env python3
"""Select existing Linux receipts; performs no constitutive/numerical evaluation."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--outer-raw', type=Path, required=True)
    parser.add_argument('--witnesses-raw', type=Path, required=True)
    args = parser.parse_args()
    root = args.evidence
    validation_path, runner_path = root/'kernel-validation.json', root/'kernel-runner.json'
    d, runner = read(validation_path), read(runner_path)
    assert sha(args.outer_raw) == d['outer_proof']['sha256']
    assert sha(args.witnesses_raw) == d['retained_witnesses']['sha256']
    assert sha(root/'kernel-contract-preregistered.md') == runner['input_hashes']['candidate-contract.md']
    assert sha(root/'kernel-run.log') == runner['log_sha256']
    outer = read(args.outer_raw)
    curated = {key: value for key, value in outer.items()
               if key not in ('branches', 'regularized_covariance_witness')}
    curated['branches'] = {}
    for branch, result in outer['branches'].items():
        curated['branches'][branch] = {
            key: value for key, value in result.items() if key != 'tangent_intervals' and key != 'ratio_leaves'}
        curated['branches'][branch]['ratio_leaves'] = [
            {key: value for key, value in leaf.items() if key != 'tangent_intervals'}
            for leaf in result['ratio_leaves']]
    curated['curation'] = dict(raw_sha256=sha(args.outer_raw),
        raw_linux_path=('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/'
                        'm8-r1/kernel-attempt-001/restricted-outer-interval-proof.json'),
        omitted='Full 9x9 gradient tangent interval matrices and the ancillary required-Gamma covariance control. All fixed-Gamma1 covariance/LDL data and all ratio-leaf proof bounds/pivots are retained.',
        method='Field selection only; no arithmetic or new proof evaluation.')
    save(root/'kernel-outer-proof.json', curated)
    # These are small state/direction/proof records, without field arrays or
    # full tangent matrices. Preserve every record, including all rejections.
    (root/'kernel-retained-witnesses.json').write_bytes(args.witnesses_raw.read_bytes())
    summary = {key: d[key] for key in (
        'schema', 'phase', 'baseline', 'started_utc', 'finished_utc', 'seconds',
        'archive_source', 'archive_source_commit', 'archive_source_commit_end',
        'archive_source_clean_start', 'archive_source_clean_end', 'command',
        'driver_sha256', 'cpu_affinity', 'thread_environment', 'versions',
        'production_adopted', 'solver_or_build_execution', 'outer_proof',
        'retained_witnesses', 'verification_success', 'restricted_constitutive_feasible',
        'decision', 'lambda_guard')}
    final_log_line = read_log_summary(root/'kernel-run.log')
    summary['checks_summary'] = {key: final_log_line[key] for key in ('checks', 'failures')}
    summary['checks_summary']['source'] = 'kernel-run.log final JSON line'
    summary['inner_formal_interval_constants'] = d['inner']['formal_interval_constants']
    summary['first_archived_inner_witness_candidate_value'] = next(
        entry for entry in d['inner']['archived_inner_minima'] if entry['pointer'] == '/events/0/detail')
    summary['source_coupling_interpretation'] = d['source_coupling']['interpretation']
    summary['preregistered_contract'] = dict(file='kernel-contract-preregistered.md',
        sha256=sha(root/'kernel-contract-preregistered.md'))
    summary['complete_validation'] = dict(file=validation_path.name, sha256=sha(validation_path))
    summary['classification'] = ('Kernel constitutive feasibility only; original and sole improved initial flow '
        'representations both fail gradient coverage. Original broad-box failures remain unchanged.')
    preflight = read(root/'preflight-summary.json')
    summary['representation_coverage'] = dict(source_file='preflight-summary.json',
        sha256=sha(root/'preflight-summary.json'), decision=preflight['decision'],
        outer_gradient_witnesses=preflight['outer_gradient_witnesses'],
        original_and_improved_cases=[{key: value for key, value in entry.items() if key in (
            'name', 'unknowns', 'outer_gradient_entry_deviation', 'outer_k_over_K0',
            'outer_epsilon_over_epsilon0', 'chi0_trace_epsilon_ratio', 'finite_admissibility',
            'initial_state_float64_LE_sha256', 'audit_sha256')} for entry in preflight['cases']])
    save(root/'kernel-summary.json', summary)
    filenames = ('kernel-validation.json', 'kernel-summary.json', 'kernel-outer-proof.json',
        'kernel-retained-witnesses.json', 'kernel-contract-preregistered.md', 'kernel-plan.json',
        'kernel-runner.json', 'kernel-run.log', 'kernel-curate.py', 'candidate-contract.md')
    manifest = dict(schema=1, phase='M8-R1 D', numerical_re_evaluation=False,
        external_run_root=d['command'][-1],
        raw_inputs={str(args.outer_raw): sha(args.outer_raw), str(args.witnesses_raw): sha(args.witnesses_raw)},
        canonical_raw_references=d['outer_proof'] | dict(retained_witnesses=d['retained_witnesses']),
        files={name: dict(sha256=sha(root/name), bytes=(root/name).stat().st_size) for name in filenames},
        method='Pure JSON field selection/re-encoding; all retained-witness bytes preserved exactly. Full interval matrices remain in the external Linux run directory.',
        scope='Constitutive kernel feasibility with native coverage failure; no production adoption, solver trajectory, state invariance or broad-box qualification.')
    save(root/'kernel-curation-index.json', manifest)
    print('Curated existing Linux receipts and recorded byte hashes; no numerical evaluation.')


def read_log_summary(path):
    for line in reversed(path.read_text().splitlines()):
        if line.startswith('{'):
            return json.loads(line)
    raise ValueError('Missing final Linux log summary')


if __name__ == '__main__':
    main()
