#!/usr/bin/env python3
"""Controller contract tests without launching native solver or MPI jobs."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("run_band100", HERE / "run_band100.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class FakeProcess:
    def __init__(self, text, stdout, running=False):
        self.pid = 12345
        self.returncode = None if running else 0
        stdout.write(text)
        stdout.flush()

    def poll(self):
        return self.returncode

    def wait(self, timeout=None):
        self.returncode = -15
        return self.returncode


class Band100Controller(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="hit-band100-controller-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_allocation_enforces_budget_and_common_fourier_partition(self):
        good = runner.allocation(list(range(48)), 48, 16, set(range(128)))
        self.assertEqual(good["ranks"], 48)
        self.assertEqual(set(good["thread_environment"].values()), {"1"})
        for cpus, ranks, npz, allowed in [
                (list(range(64)), 64, 16, set(range(128))),
                (list(range(16)), 32, 16, set(range(128))),
                (list(range(48)), 48, 24, set(range(128))),
                (list(range(48)), 48, 16, set(range(32)))]:
            with self.subTest(ranks=ranks, npz=npz), self.assertRaises(ValueError):
                runner.allocation(cpus, ranks, npz, allowed)
        for text in ("0,0", "-1", ",", ",".join(map(str, range(65)))):
            with self.subTest(text=text), self.assertRaises(ValueError):
                runner.parse_cpus(text)

    def test_bandwidth_comparison_uses_common192_slots(self):
        z, band = runner.case_knobs("z"), runner.case_knobs("bandwidth")
        self.assertEqual(z["nz"], 192)
        self.assertEqual({k:v for k,v in z.items() if k not in ("cutoff", "taper")},
                         {k:v for k,v in band.items() if k not in ("cutoff", "taper")})
        self.assertEqual(runner.COUNTS, "2548,154,154")
        self.assertEqual(runner.case_knobs("dt")["dt"] * 2, runner.case_knobs("base")["dt"])

    def test_economical_recipe_preserves_inlet_and_separates_x_and_dt_controls(self):
        original = runner.case_knobs("base")
        coarse = runner.case_knobs("base", "economical")
        self.assertEqual({k:v for k,v in original.items() if k not in ("h", "hx", "dt")},
                         {k:v for k,v in coarse.items() if k not in ("h", "hx", "dt")})
        self.assertEqual((coarse["hx"], coarse["hy"], coarse["dt"]), (.1, .05, 25e-6))
        refined_x = runner.case_knobs("x", "economical")
        self.assertEqual(refined_x["dt"], coarse["dt"])
        self.assertEqual(refined_x["hy"], coarse["hy"])
        self.assertEqual(refined_x["hx"], coarse["hx"]/2)
        self.assertEqual(runner.case_knobs("dt", "economical")["dt"], coarse["dt"]/2)
        with self.assertRaises(ValueError):
            runner.case_knobs("base", "unrecorded")

    def test_preflight_rejects_recipe_mixing(self):
        bindings = {"source": {"revision": "a"*40}, "executables": {}}
        allocated = {"ranks": 48}
        receipt = {"schema": runner.SCHEMA, "status": "COMPLETED", "bindings": bindings,
                   "configuration": {"stage": "preflight", "allocation": allocated, "recipe": "original"},
                   "budget_gate": {"passed": True}}
        (self.root / "receipt.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "recipe differs"):
            runner.check_preflight(self.root, bindings, allocated, "economical")

    def test_norm_labels_are_not_nonfinite_results(self):
        self.assertIsNone(runner.BAD_LOG.search("L inf error (variable u) : 0.123\n"))
        for text in ("L inf error (variable u) : inf\n", "Time: nan CPU Time: 1s\n",
                     "Exceeded maxIt in pressure solver", "failed to converge"):
            with self.subTest(text=text):
                self.assertIsNotNone(runner.BAD_LOG.search(text))

    def test_timing_uses_per_chunk_cpu_time_not_cumulative_differences(self):
        text = "\n".join(f"Steps: {i} Time: {i*12.5e-6} CPU Time: {16 if i<=80 else 8}s"
                         for i in range(8, 161, 8))
        gate = runner.timing_gate(text, {"wall_seconds": 270, "peak_sampled_tree_rss_kib": 2*1024**2}, runner.case_knobs("base"))
        self.assertEqual(gate["steady_seconds_per_step"], 1)
        self.assertEqual(gate["forecast_steps"], 48000)
        self.assertEqual(gate["one_off_wall_overhead_seconds"], 30)
        self.assertEqual(gate["forecast_wall_seconds"], 48030)
        self.assertFalse(gate["passed"])
        self.assertEqual(runner.log_progress(text)["step"], 160)
        with self.assertRaises(ValueError):
            runner.timing_gate("Steps: 160 Time: .002 CPU Time: 8s", {}, runner.case_knobs("base"))

    def controller(self):
        allocated = runner.allocation([0, 1], 2, 1, {0, 1})
        bindings = {"source": {"revision": "a"*40, "files_sha256": {}}, "executables": {}}
        return runner.Controller(self.root, bindings, {"stage": "preflight", "allocation": allocated})

    def test_completed_command_is_immutable_and_not_retried(self):
        controller = self.controller()
        native = self.root / "native"
        native.mkdir()
        with patch.object(runner, "verify_files"), patch.object(runner.subprocess, "Popen",
                side_effect=lambda *a, **k: FakeProcess("completed\n", k["stdout"])) as launch:
            result = controller.execute(["fake-native"], native, "native", 20)
            self.assertEqual(result["status"], "COMPLETED")
            with self.assertRaises(ValueError):
                controller.execute(["fake-native"], native, "native", 20)
            self.assertEqual(launch.call_count, 1)
        entry = controller.receipt["commands"][0]
        self.assertEqual(runner.sha(entry["receipt"]), entry["sha256"])
        self.assertEqual(json.loads((native / "native.receipt.json").read_text())["status"], "COMPLETED")

    def test_nonfinite_failure_receipt_retained_without_retry(self):
        controller = self.controller()
        native = self.root / "native"
        native.mkdir()
        with patch.object(runner, "verify_files"), patch.object(runner.subprocess, "Popen",
                side_effect=lambda *a, **k: FakeProcess("Time: nan\n", k["stdout"])) as launch:
            with self.assertRaises(ValueError):
                controller.execute(["fake-native"], native, "native", 20)
            self.assertEqual(launch.call_count, 1)
        receipt = json.loads((native / "native.receipt.json").read_text())
        self.assertEqual(receipt["status"], "FAILED")
        self.assertEqual(receipt["exit_code"], 0)
        self.assertIn("nonfinite", receipt["error"])
        self.assertEqual(runner.sha(native / "native.log"), receipt["log_sha256"])

    def test_hard_timeout_terminates_owned_process_group(self):
        controller = self.controller()
        native = self.root / "native"
        native.mkdir()
        with patch.object(runner, "verify_files"), patch.object(runner.subprocess, "Popen",
                side_effect=lambda *a, **k: FakeProcess("running\n", k["stdout"], True)), \
                patch.object(runner, "resources", return_value=(1, 1)), \
                patch.object(runner.time, "monotonic", side_effect=[0, 2, 3]), \
                patch.object(runner.os, "killpg") as stop:
            with self.assertRaisesRegex(ValueError, "hard timeout"):
                controller.execute(["fake-native"], native, "native", 1)
            stop.assert_called_once_with(12345, runner.signal.SIGTERM)
        receipt = json.loads((native / "native.receipt.json").read_text())
        self.assertEqual(receipt["status"], "FAILED")
        self.assertEqual(receipt["exit_code"], -15)

    def test_rank_wrapper_pins_and_writes_binding_before_exec(self):
        class ExecObserved(Exception):
            pass
        argv = ["--cpus", "3,7", "--binding-dir", str(self.root), "--binary", "/native/solver", "--", "case.xml"]
        with patch.dict(runner.os.environ, {"OMPI_COMM_WORLD_LOCAL_RANK": "1", "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}), \
                patch.object(runner.os, "sched_setaffinity") as pin, \
                patch.object(runner.os, "sched_getaffinity", return_value={7}), \
                patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner.os, "execv", side_effect=ExecObserved) as execute:
            with self.assertRaises(ExecObserved):
                runner.rank_entry(argv)
            pin.assert_called_once_with(0, {7})
            execute.assert_called_once_with("/native/solver", ["/native/solver", "case.xml"])
        binding = json.loads((self.root / "rank-001.json").read_text())
        self.assertEqual(binding["cpu_affinity"], [7])
        self.assertEqual(binding["rank"], 1)

    def test_status_reads_only(self):
        controller = self.controller()
        before = (self.root / "progress.json").read_bytes()
        with patch("builtins.print"):
            runner.status(self.root)
        self.assertEqual(before, (self.root / "progress.json").read_bytes())

    def test_cleanup_identity_guard_rejects_reused_or_unrelated_pid(self):
        record = {"pid": 12345, "process_start_ticks": "123", "binary": "/native/solver"}
        with patch.object(runner, "process_start_ticks", return_value="124"), \
                patch.object(runner.os, "readlink") as executable:
            self.assertFalse(runner.owned_native_alive(record))
            executable.assert_not_called()
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner.os, "readlink", return_value="/unrelated/process"):
            self.assertFalse(runner.owned_native_alive(record))
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner.os, "readlink", return_value="/native/solver (deleted)"):
            self.assertTrue(runner.owned_native_alive(record))

    def test_preflight_rejects_tampered_command_receipt(self):
        bindings = {"source": {"revision": "a"*40}, "executables": {}}
        allocated = {"ranks": 48}
        command = self.root / "solver.receipt.json"
        command.write_text("{}\n")
        receipt = {"schema": runner.SCHEMA, "status": "COMPLETED", "bindings": bindings,
                   "configuration": {"stage": "preflight", "allocation": allocated},
                   "budget_gate": {"passed": True},
                   "commands": [{"receipt": str(command), "sha256": "0"*64, "status": "COMPLETED"}]}
        (self.root / "receipt.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "receipt changed"):
            runner.check_preflight(self.root, bindings, allocated)

    def affinity_record(self):
        return {"rank": 1, "pid": 42, "process_start_ticks": "123", "cpu_affinity": [1],
                "thread_environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}

    def test_startup_movement_within_pool_then_strict_ready_affinity(self):
        audit = runner.AffinityAudit([0, 1], 2)
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42, 43]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", side_effect=[{0}, {0}, {1}, {1}]):
            audit.check_record(self.affinity_record())
            self.assertEqual(audit.startup_nonassigned_observations, 2)
            audit.observe_log("HIT_NATIVE_READY rank 0\nHIT_NATIVE_READY rank 1\n")
            audit.check_record(self.affinity_record())
            audit.observe_log("Steps: 8 Time: .0001 CPU Time: .1s", completed=True)
        summary = audit.snapshot()
        self.assertEqual(summary["ready_ranks"], [0, 1])
        self.assertEqual(summary["startup_thread_checks"], 2)
        self.assertEqual(summary["strict_thread_checks"], 2)
        self.assertEqual(summary["main_thread_checks"], 2)
        self.assertEqual(summary["helper_thread_checks"], 2)
        self.assertEqual(summary["affinity_change_witnesses"][-1]["phase"], "ready")
        self.assertEqual(summary["affinity_change_witnesses"][-1]["observed_mask"], [1])
        self.assertIn("Not established", summary["native_initialization_cause"])

    def test_startup_pool_escape_is_immediately_rejected(self):
        audit = runner.AffinityAudit([0, 1], 2)
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", return_value={2}):
            with self.assertRaisesRegex(ValueError, "escaped CPU pool") as error:
                audit.check_record(self.affinity_record())
        self.assertIn("observed mask [2]", str(error.exception))
        self.assertIn("expected mask [1]", str(error.exception))

    def test_ready_main_movement_inside_pool_is_rejected(self):
        audit = runner.AffinityAudit([0, 1], 2)
        audit.observe_log("HIT_NATIVE_READY rank 1")
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", return_value={0}):
            with self.assertRaisesRegex(ValueError, "ready MPI affinity changed") as error:
                audit.check_record(self.affinity_record())
        self.assertIn("observed mask [0]", str(error.exception))
        self.assertIn("expected mask [1]", str(error.exception))

    def test_ready_helper_thread_movement_is_rejected(self):
        audit = runner.AffinityAudit([0, 1], 2)
        audit.observe_log("HIT_NATIVE_READY rank 1")
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42, 43]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", side_effect=[{1}, {0}]):
            with self.assertRaisesRegex(ValueError, "TID 43"):
                audit.check_record(self.affinity_record())
        self.assertEqual(audit.helper_thread_checks, 1)

    def test_readiness_is_required_by_first_solver_step_and_success(self):
        audit = runner.AffinityAudit([0, 1], 2)
        audit.observe_log("HIT_NATIVE_READY rank 0")
        with self.assertRaisesRegex(ValueError, "first observed solver step"):
            audit.observe_log("Steps: 8 Time: .0001 CPU Time: .1s")
        with self.assertRaisesRegex(ValueError, "successful completion"):
            audit.observe_log("native projection complete", completed=True)
        audit.observe_log("HIT_NATIVE_READY rank 1")
        audit.observe_log("Steps: 8 Time: .0001 CPU Time: .1s")
        audit.observe_log("native projection complete", completed=True)

    def test_live_pid_reuse_is_rejected_before_reading_affinity(self):
        audit = runner.AffinityAudit([0, 1], 2)
        with patch.object(runner, "process_start_ticks", return_value="124"), \
                patch.object(runner.os, "sched_getaffinity") as affinity:
            with self.assertRaisesRegex(ValueError, "PID identity changed"):
                audit.check_record(self.affinity_record())
            affinity.assert_not_called()

    def test_startup_change_witnesses_are_bounded(self):
        audit = runner.AffinityAudit([0, 1], 2)
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", side_effect=[{i % 2} for i in range(100)]):
            for _ in range(100):
                audit.check_record(self.affinity_record())
        summary = audit.snapshot()
        self.assertEqual(len(summary["affinity_change_witnesses"]), 64)
        self.assertGreater(summary["omitted_change_witnesses"], 0)
        self.assertEqual(summary["thread_affinity_checks"], 100)

    def test_standalone_binding_validation_remains_strict(self):
        record = self.affinity_record()
        (self.root / "rank-001.json").write_text(json.dumps(record))
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", return_value={0}):
            with self.assertRaisesRegex(ValueError, "ready MPI affinity changed"):
                runner.validate_bindings(self.root, [0, 1], 2, live=True)

    def test_noninlet_ranks_use_global_native_ready_marker(self):
        audit = runner.AffinityAudit([0, 1], 2)
        audit.observe_log("HIT_EVALUATOR rank 0 points 20 method GroupedSeparable")
        self.assertEqual(audit.ready_ranks, set())
        audit.observe_log("HIT_NATIVE_READY rank 0\nHIT_NATIVE_READY rank 1\n")
        audit.observe_log("Steps: 8 Time: .0001 CPU Time: .1s", completed=True)
        # Rank1 has no evaluator report but still enters strict CPU validation.
        with patch.object(runner, "process_start_ticks", return_value="123"), \
                patch.object(runner, "task_ids", return_value=[42]), \
                patch.object(runner, "task_start_ticks", return_value="thread-start"), \
                patch.object(runner.os, "sched_getaffinity", return_value={0}):
            with self.assertRaisesRegex(ValueError, "ready MPI affinity changed"):
                audit.check_record(self.affinity_record())

    def test_evaluator_reports_allow_only_inlet_owner_subset(self):
        text = ("HIT_NATIVE_READY rank 0\n"
                "HIT_EVALUATOR rank 16 points 1152 method GroupedSeparable cache_bytes 1200000 modes 6500000\n"
                "HIT_EVALUATOR rank 32 points 1152 method GroupedTransverse cache_bytes 1200000 modes 6500000\n")
        reports = runner.evaluator_reports(text, 48)
        self.assertEqual([r["rank"] for r in reports], [16, 32])
        self.assertEqual([r["method"] for r in reports], ["GroupedSeparable", "GroupedTransverse"])
        invalid = [
            "HIT_NATIVE_READY rank 0",
            "HIT_EVALUATOR rank 0 points 0 method EmptyPlane",
            "HIT_EVALUATOR rank 48 points 10 method GroupedSeparable",
            "HIT_EVALUATOR rank 0 points 10 method ReferenceCacheLimit",
            "HIT_EVALUATOR rank 0 points 10 method GroupedSeparable\n"
            "HIT_EVALUATOR rank 0 points 10 method GroupedSeparable"]
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.evaluator_reports(value, 48)


if __name__ == "__main__":
    unittest.main()
