#!/usr/bin/env python3
"""Collect the retained step timings for the four completed Bf runs."""

import argparse
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path


DT = 0.02
MATCHED_START = 200
FINAL_STEP = 500
BAND_STEPS = 100


@dataclass(frozen=True)
class Segment:
    start: int
    stop: int
    journal: Path
    receipt: Path
    samples: Path


def read_json(path):
    return json.loads(path.read_text())


def read_journal(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def campaign(root):
    segments = []
    for entry in read_json(root / "campaign.json")["segments"]:
        path = Path(entry["path"])
        if "terminal" in entry:
            journal = Path(entry["terminal"]["journal"])
            receipt = path / "operations/retained/receipt.json"
        else:
            journal = path / "run/state.steps.jsonl"
            receipt = Path(entry["receipt"])
        segments.append(Segment(entry["start_step"], entry["through_step"],
                                journal, receipt, path / "samples.jsonl"))
    return segments


def replay(root, step):
    accepted = read_json(root / "summary.json")["accepted"]
    segment = Segment(step - 1, step, root / "step/state.steps.jsonl",
                      root / "operations/step/receipt.json", root / "samples.jsonl")
    require(read_journal(segment.journal) == [accepted], f"Replay journal mismatch at {root}")
    return [segment]


def distribution(values):
    low, *_, high = statistics.quantiles(values, n=10, method="inclusive")
    return dict(total=math.fsum(values), mean=statistics.mean(values),
                median=statistics.median(values), p10=low, p90=high)


def summarize(segments, archive):
    rows = []
    launches = []
    for segment in segments:
        require(segment.start == len(rows), f"Gap before {segment.journal}")
        journal = None
        if segment.stop > segment.start:
            retained = [row for row in read_journal(segment.journal)
                        if row["step"] <= segment.stop]
            require([row["step"] for row in retained] == list(range(segment.start + 1, segment.stop + 1)),
                    f"Step range mismatch in {segment.journal}")
            rows.extend(retained)
            journal = str(segment.journal.relative_to(archive))
        receipt = read_json(segment.receipt)
        launches.append(dict(start_step=segment.start, through_step=segment.stop,
                             journal=journal,
                             receipt=str(segment.receipt.relative_to(archive)),
                             wall_seconds=receipt["wall_seconds"]))

    require(len(rows) == FINAL_STEP and all(
        row["dt"] == DT and math.isclose(row["time"], row["step"] * DT, abs_tol=1e-12)
        and math.isfinite(row["seconds"]) and row["seconds"] > 0 for row in rows),
        f"Invalid Bf step timing or time grid for {segments[-1].journal}")
    require(all(current[previous] == last[field] for last, current in zip(rows, rows[1:])
                for field, previous in (("state_sha256", "previous_sha256"),
                                        ("advection_sha256", "previous_advection_sha256"))),
            f"History hash mismatch before {segments[-1].journal}")

    final = segments[-1]
    result = read_json(final.journal.with_name("state.result.json"))
    receipt = read_json(final.receipt)
    samples = read_journal(final.samples)
    times = result["complete_step_seconds"]
    require(final.start < MATCHED_START and len(times) == len(samples) == FINAL_STEP - final.start
            and result["state"]["step"] == FINAL_STEP
            and result["state"]["state_sha256"] == rows[-1]["state_sha256"]
            and receipt["exit_code"] == 0 and receipt["worker_pauses"] == 0,
            f"Invalid final segment in {final.journal}")
    require([row["step"] for row in samples] == list(range(final.start + 1, FINAL_STEP + 1))
            and math.isclose(math.fsum(times), result["complete_advancing_seconds"], abs_tol=1e-6)
            and all(elapsed >= row["seconds"] for elapsed, row in zip(times, rows[final.start:])),
            f"Complete step timing mismatch in {final.journal}")
    request = receipt["request"]
    require(request["command"][:3] == ["mpirun", "-np", "1"]
            and set(request["threads"].values()) == {"1"}, f"Unexpected MPI or thread settings in {final.receipt}")

    matched = rows[MATCHED_START:]
    offset = MATCHED_START - final.start
    matched_samples = samples[offset:]
    require(all(row["worker_paused_seconds"] == 0 for row in matched_samples), f"Paused matched interval in {final.journal}")
    raw = [row["seconds"] for row in matched]
    bins = matched[0]["disjoint_seconds"]
    require(all(math.isclose(math.fsum(row["disjoint_seconds"].values()), row["seconds"], abs_tol=1e-6)
                for row in matched), f"Step component timing mismatch in {final.journal}")
    last = rows[-1]
    controls = ("model", "reference_identity", "boundary_identity", "outlet_identity",
                "outlet_policy", "sizes", "quadratureX", "quadratureR", "scalar_quadrature_policy",
                "inner_k_transport", "basis")
    require(all(row[key] == last[key] for row in rows for key in controls), f"Model or mesh changed in {final.journal}")
    return {
        "model_and_mesh": {key: last[key] for key in controls},
        "final_settings": {key: last[key] for key in (
            "source_identity", "implicit_solver", "flow_preconditioner", "pressure_correctors",
            "picard_stop", "picard_outer_correctors", "nonlinear_tolerance",
            "step_advection", "picard_advection", "picard_stress_linearization") if key in last},
        "final_cpu": receipt["rank_affinity"]["allowed_cpus"],
        "final_peak_sampled_workers": receipt["peak_active_workers"],
        "final_segment_start_step": final.start,
        "final_segment_wall_seconds": receipt["wall_seconds"],
        "full_accepted_native_seconds": math.fsum(row["seconds"] for row in rows),
        "full_launch_wall_seconds": math.fsum(row["wall_seconds"] for row in launches),
        "launches": launches,
        "matched_complete_step_seconds": distribution(times[offset:]),
        "matched_native_step_seconds": distribution(raw),
        "matched_rank_cpu_seconds": distribution([row["rank_cpu_step_seconds"] for row in matched_samples]),
        "matched_mean_component_seconds": {
            key: statistics.mean(row["disjoint_seconds"][key] for row in matched) for key in bins},
        "matched_mean_pressure_solve_seconds": statistics.mean(
            row["mixed_flow_work"]["pressure_solve_seconds"] for row in matched),
        "matched_mean_pressure_krylov_iterations": statistics.mean(
            row["mixed_flow_work"]["pressure_krylov_iterations"] for row in matched),
        "matched_mean_nonlinear_iterations": statistics.mean(row["nonlinear_iterations"] for row in matched),
        "matched_mean_newton_jacobian_actions": statistics.mean(
            row["coupled_newton"]["jacobian_actions"] for row in matched),
        "matched_mean_newton_preconditioner_refreshes": statistics.mean(
            row["coupled_newton"]["refreshes"] for row in matched),
        "native_step_seconds_by_interval": {
            f"({lo * DT:g},{(lo + BAND_STEPS) * DT:g}]": distribution(
                [row["seconds"] for row in rows[lo:lo + BAND_STEPS]])
            for lo in range(0, FINAL_STEP, BAND_STEPS)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Linux directory containing the run archives")
    archive = parser.parse_args().archive.resolve()
    first = archive / "wj-h27-newton-t10-20261123.BRziGz"
    diagnostic = archive / "wj-h27-residual-gates-20261123.ChVqFQ"
    recovery = archive / "wj-h27-newton-t10-resumed-20261123.FwHrOG"
    picard = archive / "wj-h27-picard-t10-tangent-20261123.3G9WSr" / "runs"

    histories = {
        "newton-sparse": (campaign(first / "sparse-control") + replay(diagnostic / "sparse-inner", 14)),
        "newton-pressure": replay(diagnostic / "pressure-full-diagnostic", 1),
        "picard-3x1": campaign(picard / "picard-3x1"),
        "picard-4x1": campaign(picard / "picard-4x1"),
    }
    for mode, replay_name in (("sparse", "implicit-capped"), ("pressure", "capped-retry")):
        histories[f"newton-{mode}"].extend(
            campaign(recovery / mode) + replay(recovery / f"{mode}-recovery" / replay_name, 17)
            + campaign(recovery / f"{mode}-repaired") + campaign(recovery / f"{mode}-no-pause"))
    results = {case: summarize(segments, archive) for case, segments in histories.items()}
    mesh = results["newton-sparse"]["model_and_mesh"]
    require(all(result["model_and_mesh"] == mesh for result in results.values()),
            "Runs do not have matching model and mesh identities")
    print(json.dumps({"dt": DT, "matched_from_step_exclusive": MATCHED_START,
                      "through_step": FINAL_STEP, "runs": results}, indent=2))


if __name__ == "__main__":
    main()
