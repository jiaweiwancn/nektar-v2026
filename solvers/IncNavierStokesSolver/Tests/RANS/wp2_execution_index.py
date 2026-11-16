# SPDX-License-Identifier: MIT
"""Deduplicate all WP2 solver invocations, including nested/older fixtures."""
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require, sha


def executions(root):
    result = {}
    for progress in root.rglob('progress.json'):
        rows = load(progress).get('results', {})
        if not isinstance(rows, dict):
            continue
        for value in rows.values():
            if not isinstance(value, dict) or not isinstance(value.get('report'), str):
                continue
            report = Path(value['report'])
            if root not in report.parents:
                continue
            if report.exists():
                require(sha(report) == value.get('report_sha256'), 'execution report changed')
                native = load(report)
                if 'implicit_callback_count' not in native:
                    continue  # Read-only probes are charged in the spatial ledger.
                attempted = native['implicit_callback_count']
                accepted = len(native['step_audits'])
            else:
                require(value.get('report_sha256') is None, 'recorded native report missing')
                attempted = accepted = 0
            row = dict(value)
            for key, count in [('observed_physical_attempts', attempted), ('accepted_step_executions', accepted)]:
                require(key not in row or row[key] == count, 'native execution count differs: ' + str(report))
                row[key] = count
            if 'requested_attempt_work' not in row:
                tree = ET.parse(report.with_name('case.xml'))
                counts = [int(float(p.text.split('=', 1)[1])) for p in tree.findall('.//PARAMETERS/P')
                          if p.text and p.text.split('=', 1)[0].strip().lower() == 'numsteps']
                require(len(counts) == 1, 'missing/duplicate native NumSteps')
                row['requested_attempt_work'] = counts[0]
            row['counter_source'] = 'native implicit_callback_count and accepted step_audits; legacy progress normalized read-only'
            result[str(report)] = row
    return result
