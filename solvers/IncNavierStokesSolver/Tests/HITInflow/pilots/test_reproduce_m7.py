#!/usr/bin/env python3
"""Small orchestration/preservation tests; mock helpers, never run CFD/NumPy."""

import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
import reproduce_m7 as reproduction


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value)+"\n", encoding="utf-8")


class ReproductionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="hit-m7-reproduce-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.source, self.archive = self.root/"source", self.root/"archive"
        (self.source/"library/SolverUtils").mkdir(parents=True)
        tools = self.source/reproduction.RELATIVE_TOOLS
        tools.mkdir(parents=True)
        for name in reproduction.HELPERS:
            (tools/name).write_text("# mocked helper, not a scientific implementation\n")
        write_json(tools.parent/"cbc1971-source.json", {"schema": "nektar.hit.cbc1971.source.v1"})
        write_json(tools.parent/"cbc1971_full_band_correlation.json", {"schema": "nektar.hit.cbc1971.full_band_correlation.v1"})
        for production in ("cbc-production", "vk-production"):
            write_json(self.archive/production/"study.receipt.json", {"status": "completed", "cases": [
                {"id": name, "status": "completed"} for name, _, _, group in reproduction.CASES if group == production]})
        for name, _, seed, production in reproduction.CASES:
            case = self.archive/production/name
            manifest = case/"manifest.json"
            write_json(manifest, {"schema": "NEKTAR_HIT_PILOT_V1", "case": "pilot.xml", "history_file": "native.his",
                "synthesis": {"seed_uint64": str(seed)}, "time": {"record_start_s": .8, "end_s": 1.6}})
            (case/"pilot.xml").write_text("<NEKTAR />\n")
            (case/"native.his").write_text("mock raw history; never parsed as scientific data\n")
            write_json(case/"run.receipt.json", {"finished_utc": "fixture", "initialization": {"exit_code": 0},
                "solver": {"exit_code": 0}, "executed_session_sha256": reproduction.fingerprint(case/"pilot.xml")["sha256"]})
            retained = self.archive/reproduction.POST/name
            retained.mkdir(parents=True)
            cache = retained/"history.npz"
            cache.write_bytes(b"mock immutable cache; no NPZ reader is run")
            write_json(cache.with_suffix(".npz.json"), {"schema": "NEKTAR_HIT_M7_HISTORY_CACHE_V1",
                "cache_sha256": reproduction.fingerprint(cache)["sha256"], "cache_bytes": cache.stat().st_size,
                "metadata": {"manifest_sha256": reproduction.fingerprint(manifest)["sha256"], "native_history_sha256": "fixture"},
                "native_history_sha256_before_cache_write": "fixture", "native_history_sha256_after_cache_write": "fixture"})
            profiles = []
            for index, time in ((26, 1.3), (29, 1.45), (32, 1.6)):
                directory = retained/"native-profiles"/f"checkpoint-{index:03d}"
                directory.mkdir(parents=True)
                flow, slabs, info = directory/"pilot.flow.json", directory/"pilot.slabs.csv", case/f"pilot_{index}.chk/Info.xml"
                write_json(flow, {"fixture": True})
                slabs.write_text("fixture CSV; only its hash is exercised\n")
                info.parent.mkdir()
                info.write_text("<NEKTAR />\n")
                profiles.append({"time_s": time, "slabs": 90,
                    "output_files": {str(p): reproduction.fingerprint(p) for p in (flow, slabs)},
                    "checkpoint_inputs": {str(info): reproduction.fingerprint(info)}})
            write_json(retained/"native-profiles/receipt.json", {"schema": "NEKTAR_HIT_CHECKPOINT_PROFILES_V1",
                "status": "completed", "finished_utc": "fixture", "case": str(case), "profiles": profiles,
                "source_inputs": {str(p): reproduction.fingerprint(p) for p in (manifest, case/"pilot.xml")}})
        self.args = SimpleNamespace(run_root=self.archive, output=self.root/"output", source_root=self.source,
                                    cpu=3, rebuild_cache=False)

    def mock_success(self, command, **kwargs):
        self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "1")
        self.assertEqual(kwargs["env"]["PYTHONDONTWRITEBYTECODE"], "1")
        kwargs["stdout"].write("mocked command completed\n")
        output = Path(command[command.index("--output")+1])
        if Path(command[1]).name == "plot_m7.py":
            output.mkdir()
            for name in reproduction.FIGURES:
                for suffix in ("png", "svg"):
                    (output/f"{name}.{suffix}").write_bytes(b"mock figure; not production data")
            write_json(output/"provenance.json", {"schema": "NEKTAR_HIT_M7_FIGURES_V1"})
        else:
            schema = ("NEKTAR_HIT_M7_ANALYSIS_V1" if Path(command[1]).name == "analyze_m7.py"
                      else "NEKTAR_HIT_M7_NATIVE_WORK_V1")
            write_json(output, {"schema": schema, "mocked_command": Path(command[1]).name})
        return SimpleNamespace(returncode=0)

    def test_fixed_case_order_and_exact_frozen_arguments(self):
        plan = reproduction.preflight(self.args)
        steps = reproduction.commands(plan, False)
        self.assertEqual([c["id"] for c in plan["cases"]], ["cbc-seed-20261028", "cbc-seed-4315228324",
            "vk-seed-20261028", "vk-seed-4315228324"])
        self.assertEqual(len(steps), 9)
        for index in range(0, 8, 2):
            command = steps[index][1]
            for flag, value in (("--start", "0.8"), ("--end", "1.6"), ("--probe-batch", "256"),
                                ("--welch-segment-samples", "512"), ("--acf-max-lag-s", "0.1")):
                self.assertEqual(command[command.index(flag)+1], value)
            self.assertNotIn("--history", command)
            self.assertTrue(Path(command[command.index("--cache")+1]).is_relative_to(self.archive))
        self.assertEqual(steps[-1][1].count("--native"), 4)
        self.assertEqual(len(steps[-1][2]), 11)
        self.assertFalse(any("IncNavierStokesSolver" == word or "HITFlowDiagnostics" == word
                             for _, command, _ in steps for word in command))
        self.assertFalse(self.args.output.exists())

    def test_rebuild_cache_uses_only_fresh_output_paths(self):
        self.args.rebuild_cache = True
        plan = reproduction.preflight(self.args)
        for case in plan["cases"]:
            self.assertTrue(case["cache"].is_relative_to(self.args.output))
            self.assertFalse(case["cache"].exists())
        command = reproduction.commands(plan, True)[0][1]
        self.assertEqual(command[command.index("--cache-start")+1], "0")
        self.assertTrue(Path(command[command.index("--history")+1]).is_relative_to(self.archive))

    def test_fresh_output_and_symlink_alias_guards_precede_creation(self):
        alias = self.root/"archive-alias"
        alias.symlink_to(self.archive, target_is_directory=True)
        for output in (self.archive, self.archive/"new", self.source/"new", alias/"new"):
            with self.subTest(output=output), self.assertRaisesRegex(ValueError, "fresh|outside"):
                self.args.output = output
                reproduction.preflight(self.args)
        self.assertFalse((self.archive/"new").exists())
        self.assertFalse((self.source/"new").exists())

    def test_unsuccessful_initial_receipts_rejected_without_output(self):
        native = self.archive/reproduction.POST/reproduction.CASES[0][0]/"native-profiles/receipt.json"
        study = self.archive/"cbc-production/study.receipt.json"
        run = self.archive/"cbc-production"/reproduction.CASES[0][0]/"run.receipt.json"
        for path, field, value, message in ((study, "status", "failed", "study receipt"),
                (run, "solver", {"exit_code": 1}, "solver receipt"), (native, "status", "running", "native profiles")):
            with self.subTest(path=path):
                original = path.read_bytes()
                modified = json.loads(original)
                modified[field] = value
                write_json(path, modified)
                with self.assertRaisesRegex(ValueError, message):
                    reproduction.preflight(self.args)
                self.assertFalse(self.args.output.exists())
                path.write_bytes(original)

    def test_completed_receipt_and_archive_immutability_with_mocked_helpers(self):
        before = {str(p): reproduction.fingerprint(p) for root in (self.source, self.archive)
                  for p in root.rglob("*") if p.is_file()}
        with patch.object(reproduction, "pin_cpu", return_value=[3]):
            path = reproduction.reproduce(self.args, runner=self.mock_success)
        receipt = json.loads(path.read_text())
        self.assertEqual(receipt["status"], "completed")
        self.assertEqual(len(receipt["steps"]), 9)
        self.assertTrue(all(step["status"] == "completed" and step["outputs"] for step in receipt["steps"]))
        self.assertEqual(receipt["cpu_affinity"], [3])
        for name, expected in before.items():
            self.assertEqual(reproduction.fingerprint(Path(name)), expected)

    def test_failed_command_receipt_is_atomic_and_not_retried(self):
        calls = []
        def failed(command, **kwargs):
            calls.append(command)
            kwargs["stdout"].write("mock failure\n")
            return SimpleNamespace(returncode=7)
        with patch.object(reproduction, "pin_cpu", return_value=[3]), self.assertRaisesRegex(ValueError, "exited 7"):
            reproduction.reproduce(self.args, runner=failed)
        receipt = json.loads((self.args.output/"receipt.json").read_text())
        self.assertEqual(receipt["status"], "failed")
        self.assertEqual(len(calls), 1)
        self.assertEqual(receipt["steps"][0]["exit_code"], 7)
        self.assertEqual(receipt["steps"][0]["status"], "failed")
        self.assertFalse((self.args.output/".receipt.json.tmp").exists())
        self.assertTrue(receipt["steps"][0]["log_fingerprint"]["sha256"])

    def test_external_input_mutation_fails_after_first_command(self):
        changed = self.source/reproduction.RELATIVE_TOOLS/"analyze_history.py"
        calls = []
        def mutated(command, **kwargs):
            calls.append(command)
            result = self.mock_success(command, **kwargs)
            changed.write_text("# external mutation detected by driver\n")
            return result
        with patch.object(reproduction, "pin_cpu", return_value=[3]), self.assertRaisesRegex(ValueError, "immutable input changed"):
            reproduction.reproduce(self.args, runner=mutated)
        receipt = json.loads((self.args.output/"receipt.json").read_text())
        self.assertEqual(receipt["status"], "failed")
        self.assertEqual(len(calls), 1)
        self.assertEqual(receipt["steps"][0]["status"], "failed")

    def test_single_cpu_affinity_is_explicit_and_restricted(self):
        with patch.object(reproduction.sys, "platform", "linux"), \
                patch.object(reproduction.os, "sched_getaffinity", side_effect=({2, 3}, {2, 3}, {2}), create=True), \
                patch.object(reproduction.os, "sched_setaffinity", create=True) as setter:
            with self.assertRaisesRegex(ValueError, "inherited affinity"):
                reproduction.pin_cpu(4)
            self.assertEqual(reproduction.pin_cpu(2), [2])
            setter.assert_called_once_with(0, {2})


if __name__ == "__main__":
    unittest.main()
