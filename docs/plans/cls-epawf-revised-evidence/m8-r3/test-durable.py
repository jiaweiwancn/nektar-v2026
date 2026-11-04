#!/usr/bin/env python3
"""Linux I/O/SHA/timeout metadata checks using the tiny durable C++ fixture."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import struct
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('durable_collector', HERE/'collect-progress.py')
collector = importlib.util.module_from_spec(spec); spec.loader.exec_module(collector)
H0 = hashlib.sha256(struct.pack('<3d', 1., 2., 3.)).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if platform.system() != 'Linux' or a.output.exists():
        raise ValueError('Authoritative Linux and a fresh external fixture directory required')
    a.output.mkdir(parents=True); checks = []
    complete = a.output/'complete/progress.jsonl'; complete.parent.mkdir()
    subprocess.run([str(a.probe.resolve()), str(complete), 'complete'], check=True, timeout=10)
    result = collector.collect(complete, H0, 3)
    assert result['completed_iterate_count'] == 1 and result['jacobians'][0]['status'] == 'complete'
    assert result['native_converged_observation'] is False
    checks.append('compiled_SHA256_known_vectors_and_complete_payloads')
    partial = a.output/'disconnected/progress.jsonl'; partial.parent.mkdir()
    receipt = a.output/'disconnected/receipt.json'
    # The outer metadata supervisor uses the actual R3 emitter. The producing
    # child is only the I/O fixture, with all output redirected to its own log.
    script = ('import importlib.util,json,os,signal,subprocess,sys,time\nfrom pathlib import Path\n'
        'spec=importlib.util.spec_from_file_location("r3_runner",sys.argv[1]);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)\n'
        'sys.stdin.read(1)\n'
        'log=Path(sys.argv[3]).with_suffix(".log").open("w")\n'
        'child=subprocess.Popen([sys.argv[2],sys.argv[3],"hold"],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)\n'
        'start=time.monotonic()\ntry:\n'
        '    while child.poll() is None:\n'
        '        m.emit("metadata child alive")\n'
        '        if time.monotonic()-start>.5:\n'
        '            os.killpg(child.pid,signal.SIGTERM);child.wait(timeout=10);break\n'
        '        time.sleep(.02)\nfinally:\n'
        '    if child.poll() is None:os.killpg(child.pid,signal.SIGKILL);child.wait()\n'
        '    log.close()\n'
        '    Path(sys.argv[4]).write_text(json.dumps({"exit_code":124,"timed_out":True,"child_waited":True}))\n')
    process = subprocess.Popen([sys.executable, '-B', '-c', script, str(HERE/'run-stage.py'),
        str(a.probe.resolve()), str(partial), str(receipt)], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    process.stdout.close(); process.stdin.write('x'); process.stdin.flush(); process.stdin.close()
    code = process.wait(timeout=15); errors = process.stderr.read(); process.stderr.close()
    assert code == 0, errors
    assert errors == '' and json.loads(receipt.read_text())['child_waited']
    stopped = collector.collect(partial, H0, 3)
    assert stopped['last_snapshot']['iteration'] == 0 and stopped['jacobians'][0]['status'] == 'incomplete'
    assert stopped['native_converged_observation'] is None and not stopped['terminal_record_available']
    checks.append('stdout_disconnect_preserves_timeout_wait_receipt_and_complete_last_state')
    with partial.open('ab') as stream: stream.write(b'{"unpublished":')
    assert collector.collect(partial, H0, 3)['progress']['unpublished_tail'] is not None
    checks.append('unpublished_event_tail_does_not_erase_complete_state')
    state = Path(stopped['last_snapshot']['state']['path'])
    state.write_bytes(b'!'+state.read_bytes()[1:])
    try: collector.collect(partial, H0, 3)
    except ValueError: pass
    else: raise AssertionError('changed completed binary payload accepted')
    checks.append('changed_completed_binary_payload_rejected')
    (a.output/'validation.json').write_text(json.dumps(dict(success=True, native_solver_calls=0,
        scientific_operator_calls=0, checks=checks, probe=collector.identity(a.probe)), indent=2)+'\n')
    print(json.dumps(dict(success=True, checks=len(checks))))


if __name__ == '__main__':
    main()
