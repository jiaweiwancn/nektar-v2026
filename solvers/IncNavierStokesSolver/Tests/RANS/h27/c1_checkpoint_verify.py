#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exercise the native accepted observer and atomic diagnostic recovery on Linux."""

import argparse
import json
from pathlib import Path
import signal
import subprocess
import time
import xml.etree.ElementTree as ET

from resources import identity


SUITE = Path(__file__).resolve().parent


def write_json(path, result):
    with path.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


def packet(path):
    return json.loads(path.read_text())


def solve(executable, case, seconds, expected=0):
    with (case/"native.log").open("xb") as log:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", "cost.txt", str(executable), "run.xml"],
                                cwd=case, stdout=log, stderr=subprocess.STDOUT, timeout=seconds)
    if result.returncode != expected:
        raise RuntimeError(f"{case.name}: native exit {result.returncode}, expected {expected}")


def prepare(template, case, *options):
    subprocess.run(["python3", str(SUITE/"c1_steady_case.py"), str(template), str(case),
                    "--seconds", "130", *map(str, options)], check=True, timeout=25)


def interrupt_later_preparation(executable, case):
    log_path = case/"native.log"
    begin = time.monotonic()
    with log_path.open("xb") as log:
        with subprocess.Popen([str(executable), "run.xml"], cwd=case, stdout=log,
                              stderr=subprocess.STDOUT, start_new_session=True) as process:
            while (time.monotonic()-begin < 110 and process.poll() is None and
                   log_path.read_bytes().count(b"WJ-H27 inverse preparation start;") < 2):
                time.sleep(.015)
            reached = log_path.read_bytes().count(b"WJ-H27 inverse preparation start;")
            running = process.poll() is None
            if running:
                process.send_signal(signal.SIGKILL)
                process.wait(timeout=10)
            if reached != 2 or not running or process.returncode != -signal.SIGKILL:
                raise RuntimeError("owned interruption did not reach the second inverse preparation")
    saved, = sorted(case.glob("state.accepted-*.json"))
    assert not (case/"state.diagnostic-steady.json").exists()
    doc = packet(saved)
    assert doc["accepted_checkpoint"]["sequence"] == 1
    assert doc["accepted_checkpoint"]["session_sha256"] == [identity(case/"run.xml")["sha256"]]
    assert doc["iterate_custody"] == "LAST_ACCEPTED" and not doc["physical_history"]
    assert max(doc["original_residuals"]) > 1e-10
    return saved, time.monotonic()-begin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--flat", type=Path, required=True)
    parser.add_argument("--hill", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = args.out.resolve()
    output.mkdir()
    off, on = output/"off", output/"on"
    prepare(args.flat, off)
    solve(args.exe, off, 155)
    prepare(args.flat, on, "--accepted-checkpoint", 155, 8, 0)
    solve(args.exe, on, 155)
    unprotected = packet(off/"state.diagnostic-steady.json")
    protected = packet(on/"state.diagnostic-steady.json")
    assert unprotected["state"]["state_sha256"] == protected["state"]["state_sha256"]
    assert unprotected["original_residuals"] == protected["original_residuals"]
    accepted = sorted(on.glob("state.accepted-*.json"))
    expected = packet(on/"state.steady.json")["accepted_newton_steps"]
    assert len(accepted) == expected > 1
    for sequence, path in enumerate(accepted, 1):
        item = packet(path)
        assert item["accepted_checkpoint"]["sequence"] == sequence
        assert item["accepted_checkpoint"]["session_sha256"] == [identity(on/"run.xml")["sha256"]]
        assert item["iterate_custody"] == "LAST_ACCEPTED"
    assert packet(accepted[-1])["state"]["state_sha256"] == protected["state"]["state_sha256"]

    stopped = output/"interrupted"
    prepare(args.hill, stopped, "--seconds", 220, "--inverse", "ColoredSparseV1",
            "--accepted-checkpoint", 245, 8, 0)
    saved, stopped_seconds = interrupt_later_preparation(args.exe, stopped)
    replay = output/"interrupted-replay"
    prepare(stopped/"run.xml", replay, "--seed", saved, "--evaluate-only", "--check-actions")
    solve(args.exe, replay, 155)
    original, recovered = packet(saved), packet(replay/"state.diagnostic-steady.json")
    assert original["state"]["state_sha256"] == recovered["state"]["state_sha256"]
    assert max(abs(a-b) for a,b in zip(original["original_residuals"], recovered["original_residuals"])) < 5e-12

    corrupt = output/"corrupt-packet-binding"
    prepare(stopped/"run.xml", corrupt, "--seed", saved, "--evaluate-only")
    xml = ET.parse(corrupt/"run.xml")
    xml.find(".//I[@PROPERTY='WJDiagnosticSeedSHA256']").set("VALUE", "0"*64)
    xml.write(corrupt/"run.xml", encoding="unicode")
    solve(args.exe, corrupt, 70, expected=1)
    assert "diagnostic seed packet SHA256 differs" in (corrupt/"native.log").read_text()
    assert not list(corrupt.glob("state.*"))

    release = output/"rejected-physical-release"
    subprocess.run(["python3", str(SUITE/"c1_release_original.py"), str(stopped/"run.xml"),
                    str(release), "--packet", str(saved), "--equations", str(saved),
                    "--response", str(saved), "--points", str(saved), "--steps", "1"],
                   check=True, timeout=20)
    solve(args.exe, release, 70, expected=1)
    assert "physical release requires the exact current zero-regularization steady equations" in (
        release/"native.log").read_text()
    assert not list(release.glob("state.recovery-*"))
    write_json(output/"receipt.json", dict(
        executable=identity(args.exe), source=protected["state"]["source_identity"],
        flat=identity(args.flat), hill=identity(args.hill), accepted=[identity(p) for p in accepted],
        endpoint_coefficient_sha256=protected["state"]["state_sha256"], endpoint=protected["original_residuals"],
        interrupted_coefficient_sha256=original["state"]["state_sha256"], interrupted_seconds=stopped_seconds,
        interrupted_packet=identity(saved), replay_packet=identity(replay/"state.diagnostic-steady.json"),
        corrupt_rejected=True, physical_release_rejected=True, interrupted_after_inverse_start=2))


if __name__ == "__main__":
    main()
