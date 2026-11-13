#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential CT6 controls and complete read-only A95 endpoint audits.

Run inside ct_accepted3_command.py. Each destination is new. This driver never
resumes a failed trajectory or admits the conditional A93/E94E95 evolutions.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from ct_accepted_open import require, run
from rans_request_manifest import RequestManifest, runtime_identity, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'runtime', 'ct5', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args();test = Path(__file__).parent
    members = [a.root/'evolution-v2/A95-p4-dt002']
    members += [a.root/'evolution-v3'/n for n in
                ('A95-p4-dt001', 'A95-p4-dt004', 'A95-p3-dt002', 'A95-p3-dt001')]
    for folder in members:
        q = json.loads((folder/'results.json').read_text())
        require(q['physical_completion'] == 'pass', 'required A95 member did not complete: '+str(folder))
    manifest = RequestManifest(a.output, dict(schema='CT6CompleteVerificationDriverV1',
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver', 'RANSPreservationProbe')},
        members={str(f):sha(f/'results.json') for f in members},
        scripts={n:sha(test/n) for n in ('ct_accepted6_finish.py', 'ct_accepted6_audit.py',
            'ct_accepted6_controls.py', 'ct_accepted6_regression.py', 'ct_accepted6_compare.py',
            'ct_accepted6_gate.py', 'ct_accepted6_shared_gates.py', 'test_ct_accepted6_gate.py')},
        source_revision=subprocess.check_output(['git', '-c', 'core.preloadIndex=false',
            '-c', 'index.threads=1', 'rev-parse', 'HEAD'], cwd=test, text=True).strip(),
        scope='diagnostic/restart/rollback/flat-regression only; no new duration or conditional model transfers'))
    rows = {}

    def execute(label, command):
        rows[label] = run(command, a.output, label+'.log')
        write_json(a.output/'progress.json', rows)
        print('CT6 finishing control', label, 'exit', rows[label]['exit_code'], flush=True)
        return rows[label]['exit_code']

    configured = execute('configure-ctest', ['cmake', '-S', test, '-B', a.root/'build']) == 0
    if configured:
        execute('ctest', ['ctest', '--test-dir', a.root/'build', '-R', '^CT6',
                          '--parallel', '1', '--output-on-failure'])
    execute('discriminating', [sys.executable, '-B', test/'ct_accepted6_controls.py',
        '--runtime', a.runtime, '--operator', a.root/'native-operator-v1',
        '--targeted', a.root/'targeted-v4', '--output', a.root/'controls-v1'])
    execute('regression', [sys.executable, '-B', test/'ct_accepted6_regression.py',
        '--runtime', a.runtime, '--ct5', a.ct5, '--member', members[0],
        '--output', a.root/'regression-v1'])
    audits = []
    for member in members:
        output = a.root/'physical-audits-v3'/member.name
        execute(member.name, [sys.executable, '-B', test/'ct_accepted6_audit.py',
            '--runtime', a.runtime, '--member', member, '--output', output])
        audits.append(output)
    command = [sys.executable, '-B', test/'ct_accepted6_compare.py',
               '--output', a.root/'primary-comparisons-v2']
    for member, audit in zip(members, audits, strict=True):
        command += ['--member', member, '--audit', audit]
    execute('complete-comparison', command)
    result = dict(schema='CT6CompleteVerificationDriverV1', executions=rows,
                  status='pass' if all(q['exit_code'] == 0 for q in rows.values()) else 'failed_or_unresolved_required_check',
                  model_transfers='not_admitted_by_process_completion')
    write_json(a.output/'results.json', result)
    manifest.seal(a.output, ['results.json'])
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
