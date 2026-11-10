#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Run required-export failures through the actual command-line classifier."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class NonlinearEvidence(unittest.TestCase):
    def test_command_line_contract(self):
        with tempfile.TemporaryDirectory(prefix='rans-nonlinear-evidence-',
                                         dir=os.environ.get('TMPDIR')) as temporary:
            output = Path(temporary)/'fixtures'
            code = subprocess.run([sys.executable, '-B', str(Path(__file__).with_name(
                'check_rans_nonlinear_evidence_failures.py')), '--output', str(output)]).returncode
            self.assertEqual(code, 0)
            data = json.loads((output/'results.json').read_text())
            self.assertEqual(data['status'], 'pass')
            self.assertEqual(len(data['checks']), 8)


if __name__ == '__main__':
    unittest.main()
