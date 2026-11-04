"""Exercise the actual runner progress emitter without running its CLI/children."""
import ast
import errno
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
RUNNER = Path(__file__).with_name('run-stages.py')
TREE = ast.parse(RUNNER.read_text())
EMITTER = next(node for node in TREE.body if isinstance(node, ast.ClassDef) and node.name == 'ProgressEmitter')
NAMESPACE = dict(errno=errno, os=os, sys=sys)
exec(compile(ast.Module(body=[EMITTER], type_ignores=[]), str(RUNNER), 'exec'), NAMESPACE)
ProgressEmitter = NAMESPACE['ProgressEmitter']


class FailingStream:
    def __init__(self, error):
        self.error, self.writes = error, 0

    def write(self, value):
        self.writes += 1
        raise self.error

    def flush(self):
        raise self.error

    def fileno(self):
        raise ValueError('closed metadata fixture stream')


class ProgressTests(unittest.TestCase):
    def test_normal_progress_is_flushed_and_enabled(self):
        stream = io.StringIO(); emitter = ProgressEmitter(stream)
        self.assertTrue(emitter('first'))
        self.assertTrue(emitter('second'))
        self.assertEqual(stream.getvalue(), 'first\nsecond\n')
        self.assertFalse(emitter.disabled)

    def test_broken_pipe_disables_subsequent_progress(self):
        stream = FailingStream(BrokenPipeError()); emitter = ProgressEmitter(stream)
        self.assertFalse(emitter('first')); self.assertFalse(emitter('second'))
        self.assertEqual(stream.writes, 1)
        self.assertTrue(emitter.disabled)

    def test_epipe_oserror_is_also_silenced(self):
        stream = FailingStream(OSError(errno.EPIPE, 'fixture pipe closed'))
        emitter = ProgressEmitter(stream)
        self.assertFalse(emitter('first')); self.assertFalse(emitter('second'))
        self.assertEqual(stream.writes, 1)

    def test_other_oserror_is_not_hidden(self):
        emitter = ProgressEmitter(FailingStream(OSError(errno.EIO, 'fixture I/O error')))
        with self.assertRaises(OSError): emitter('first')
        self.assertFalse(emitter.disabled)

    def test_fake_closed_stdout_is_detached_and_finally_still_runs(self):
        stream = FailingStream(BrokenPipeError()); events = []
        with patch.object(sys, 'stdout', stream):
            emitter = ProgressEmitter()
            try:
                for state in ('poll-running', 'poll-running', 'timeout-enforced'):
                    emitter('progress'); events.append(state)
            finally:
                events.append('final-receipt-written')
            self.assertIsNone(sys.stdout)
        self.assertEqual(events, ['poll-running', 'poll-running', 'timeout-enforced', 'final-receipt-written'])
        self.assertEqual(stream.writes, 1)

    def test_real_disconnected_stdout_has_no_final_flush_exit_120(self):
        # Only the isolated emitter executes in this Python child. Its stdout
        # reader is deliberately closed before the first progress write.
        with tempfile.TemporaryDirectory(prefix='r2-progress-metadata-') as temporary:
            marker = Path(temporary)/'finally-receipt.json'
            child = ('import errno,json,os,sys\nfrom pathlib import Path\n'+ast.unparse(EMITTER)+'\n'
                'sys.stdin.read(1)\nemitter=ProgressEmitter()\nevents=[]\ntry:\n'
                '    for value in ("poll1","poll2","timeout"):\n'
                '        emitter(value)\n        events.append(value)\nfinally:\n'
                '    Path(sys.argv[1]).write_text(json.dumps(events+["finally"]))\n')
            process = subprocess.Popen([sys.executable, '-B', '-c', child, str(marker)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            process.stdout.close()
            process.stdin.write('x'); process.stdin.flush(); process.stdin.close()
            code = process.wait(timeout=15)
            errors = process.stderr.read(); process.stderr.close()
            self.assertEqual(code, 0, errors)
            self.assertEqual(errors, '')
            self.assertEqual(json.loads(marker.read_text()), ['poll1', 'poll2', 'timeout', 'finally'])


if __name__ == '__main__':
    unittest.main()
