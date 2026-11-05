#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh original residual, continuous-ground and constitutive audit only.

The supplied alpha/time/history identify the target. Physical acceptance or
history advancement is never inferred from this fixed-state observation.
"""
import argparse
from pathlib import Path
import platform
import shutil
import subprocess

from run_cls_epawf_picard import run, load, require, sha, state_hash, write_json, ROOT, THREADS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("session", type=Path)
    parser.add_argument("packet", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-ground-guard", action="store_true")
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.probe = args.probe.resolve(strict=True)
    require(platform.system() == "Linux" and args.output.is_relative_to(
        Path("/media/jay/新加卷/nektar-v2026-runs")), "use designated Linux storage")
    args.output.mkdir()
    packet = load(args.packet)
    write_json(args.output / "input.json", dict(packet, audit_original_qualification=True,
                                              require_history_h0=False))
    shutil.copyfile(args.session, args.output / "case.xml")
    linked = subprocess.check_output(["ldd", str(args.probe)], text=True)
    require("not found" not in linked, "unresolved audit dependency")
    libraries = {Path(p).resolve() for p in linked.split() if p.startswith("/")}
    library_hashes = {str(p): sha(p) for p in sorted(libraries)}
    receipt = run([str(args.probe), "case.xml", "input.json", "readback.json"],
                  args.output, "readback.log", 240)
    require(receipt["exit_code"] == 0, "saved-state native audit failed")
    data = load(args.output / "readback.json")
    require(data["original_qualification_requested"] and data["continuous_ground"]["admissible"] and
            data["scalar_domain"]["admissible"], "native ground/domain audit was not exercised")
    require(data["state"] == packet["state"] and data["history"] == packet["history"] and
            data["alpha"] == packet["alpha"] and data["time"] == packet["time"],
            "saved-state audit changed the target or incoming coefficients")
    summary = dict(scope=__doc__, source=subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        probe=str(args.probe), probe_sha256=sha(args.probe), libraries=library_hashes,
        session_sha256=sha(args.session), packet_sha256=sha(args.packet),
        state_sha256=state_hash(data["state"]), history_sha256=state_hash(data["history"]),
        unknowns=len(data["state"]), alpha=data["alpha"], time=data["time"],
        ranks=1, thread_limits=THREADS, receipt=receipt,
        original_norms=data["norms"]["native_stage"], ground=data["continuous_ground"],
        qualification=data["qualification"], physical_acceptance_inferred=False)
    if args.check_ground_guard:
        # A free ground vertex is nodal: its sign is a physical value here,
        # unlike the sign of an interior Modified_A modal coefficient.
        field = data["fields"][2]
        left = min(e["lower"][0] for e in field["mesh"])
        element = next(e for e in field["mesh"] if e["lower"][2] == 0 and e["lower"][0] > left)
        state = list(packet["state"])
        index = int(data["offsets"][2]+element["ids"][0])
        state[index] = -data["parameters"][5]/element["sign"][0]
        write_json(args.output / "negative-ground.json", dict(packet, state=state,
            audit_original_qualification=True, require_history_h0=False))
        guard = run([str(args.probe), "case.xml", "negative-ground.json", "must-not-exist.json"],
                    args.output, "negative-ground.log", 120)
        require(guard["exit_code"] != 0 and not (args.output / "must-not-exist.json").exists() and
                "continuous ground domain before scalar evaluation" in
                (args.output / "negative-ground.log").read_text(),
                "negative physical ground reached scalar evaluation")
        summary["negative_ground_guard"] = dict(guard,
            rejected_before_scalar_evaluation=True, changed_row=index)
    require(all(sha(Path(p)) == value for p,value in library_hashes.items()),
            "audit runtime changed during readback")
    write_json(args.output / "audit.json", summary)


if __name__ == "__main__":
    main()
