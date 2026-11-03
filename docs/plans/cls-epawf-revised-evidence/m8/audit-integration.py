"""Check M8 evidence integrity and gate labels on Linux; no numerical solves."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert platform.system() == 'Linux'
    source, build = args.source.resolve(), args.build.resolve()
    evidence = source/'docs/plans/cls-epawf-revised-evidence/m8'
    def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    def read(name): return json.loads((evidence/name).read_text())
    def git(*argv): return subprocess.check_output(['git', '-C', str(source), *argv], text=True).strip()
    assert not git('status', '--porcelain')
    checked = []
    for name in ('curation-index.json', 'curved-curation-index.json', 'mms-curation-index.json',
                 'capped-native-curation-index.json', 'native-curation-index.json'):
        index = read(name)
        for record in index['files']:
            filename = record.get('name', record.get('curated', record.get('file')))
            expected = record.get('curated_sha256', record.get('sha256'))
            assert sha(evidence/filename) == expected, (name, filename)
            checked.append(dict(index=name, file=filename, sha256=expected))
    epoch = read('source-epochs.json')
    assert epoch['generator_sha256'] == sha(evidence/'source-epochs.py')
    for filename, record in epoch['frozen_inputs_unchanged_from_M7'].items():
        assert sha(source/filename) == record['sha256'], filename
    for filename, record in epoch['kernel_inputs'].items():
        assert sha(source/filename) == record['sha256'], filename
    delta = git('diff', '--name-only', epoch['commits']['ed34052'], 'HEAD').splitlines()
    guide = 'solvers/IncNavierStokesSolver/Examples/RevisedEquilibriumWall/README.md'
    assert all(name.startswith('docs/') or name in ('README.md', guide) for name in delta), delta
    definitions = json.loads(subprocess.check_output(
        ['ctest', '--test-dir', str(build), '--show-only=json-v1'], text=True))
    registered = {test['name'] for test in definitions['tests']}
    legacy = ET.parse(evidence/'legacy-ctest.xml').getroot().findall('testcase')
    legacy_names = {test.get('name') for test in legacy}
    assert len(registered) == 38 and len(legacy) == len(legacy_names) == 33
    assert legacy_names == {name for name in registered if not name.startswith('RANSCLSEPAWFM8')}
    assert all(not test.findall('failure') and not test.findall('skipped') for test in legacy)
    for name in ('native-scaled-validation.json', 'native-regularized-validation.json',
                 'mms-initialization-qualification.json'):
        receipt = read(name)
        assert receipt['implementation_verified'] and not receipt['candidate_qualified']
        assert not receipt['failures']
        assert all(case['status'] == 'verified' for case in receipt['cases'])
    cap = read('capped-native-validation.json')
    assert cap['success'] and cap['mechanism_verified'] and cap['evolution_verified']
    assert len(cap['cases']) == 6 and not cap['candidate_qualified'] and not cap['failures']
    assert read('original-controls-validation.json')['status'] == 'passed'
    compat = read('mms-checkpoint-compatibility.json')
    assert compat['status'] == 'passed' and compat['unaffected_neutral_signature_exactly_preserved']
    stale = compat['stale_functions_control']
    assert 'checkpoint signature mismatch' in stale['failure']
    assert stale['accepted_states'] == 0 and stale['exact_signature_change_only_scoped_marker']
    for name in ('fresh_MMS_complete_MPI', 'fresh_MMS_complete_restart', 'old_neutral_complete_restart'):
        result = compat[name]
        assert result['maximum_field_mass_norm'] <= result['field_and_history_tolerance']
        assert result['maximum_history_mass_norm'] <= result['field_and_history_tolerance']
    assert read('proof-kernel-fresh.json')['status'] == 'passed'
    assert read('transverse-validation.json')['success']
    assert read('proof-fixed-gamma-feasibility.json')['conclusion'] == 'infeasible_on_declared_box'
    assert read('proof-invariant-screen.json')['conclusion'] == 'declared_box_gates_failed_by_explicit_witnesses'
    for box in read('proof-local-joint-11aa.json')['finite_boxes']:
        gamma = '1' if box['gradient_entry_radius'] == '1/1000' else '3'
        candidate = next(value for value in box['fixed_gamma_candidates'] if value['gamma_over_nut0'] == gamma)
        assert candidate['joint_constitutive_gate']
    disposition = read('integration-remaining-gates.json')
    assert not disposition['milestone_closed'] and not disposition['candidate_qualified']
    assert not disposition['experimental_accuracy_claimed']
    for filename in evidence.glob('*.json'): json.loads(filename.read_text())
    result = dict(schema=1, status='passed', kind='Evidence integrity and gate-label audit',
        generated_utc=datetime.now(timezone.utc).isoformat(), source_commit=git('rev-parse', 'HEAD'),
        script_sha256=sha(__file__), indexed_copies_checked=len(checked), copies=checked,
        registered_test_count=len(registered), earlier_ctest_passes=len(legacy),
        numerical_source_commit=epoch['commits']['ed34052'],
        changes_since_numerical_source=delta, milestone_closed=False, candidate_qualified=False,
        interpretation='Implementation receipts and reproduced failure witnesses verified; '
                       'broad robust qualification failed and conditional terrain remains unperformed. '
                       'This metadata audit runs no model calculation or solver.')
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({key: result[key] for key in ('status', 'source_commit', 'indexed_copies_checked',
          'registered_test_count', 'earlier_ctest_passes', 'milestone_closed', 'candidate_qualified')}))


if __name__ == '__main__':
    main()
