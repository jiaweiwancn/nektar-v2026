#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused continuation custody and native termination controls."""

import json
from pathlib import Path
import tempfile
import unittest

from c1_steady_continue import accepted_packet, termination


class ContinuationTests(unittest.TestCase):
    def test_improving_deadline_uses_accepted_progress(self):
        ledger = dict(reason="DEADLINE: inverse preparation", accepted_newton_steps=1,
                      last_accepted_norm=0.6, true_norm=0.4,
                      best_from_accepted_iterate=False, newton=[{"base_norm": 1.0}])
        self.assertEqual(termination(1, ledger, True), "BUDGET_EXHAUSTED_WITH_PROGRESS")
        ledger.update(accepted_newton_steps=0, last_accepted_norm=1.0)
        self.assertEqual(termination(1, ledger, True), "BUDGET_EXHAUSTED")
        self.assertEqual(termination(124, None, False), "DEADLINE_WITHOUT_CUSTODY")

    def test_distinct_failure_and_audit_decisions(self):
        for reason, expected in [
                ("LINEAR_NOT_CONVERGED: check", "CURRENT_LINEAR_FAILURE"),
                ("LINE_SEARCH_NOT_CONVERGED: rejected", "GLOBALIZATION_FAILURE"),
                ("SPARSE_FD_INFEASIBLE: domain", "CURRENT_LINEAR_FAILURE"),
                ("NEWTON_NOT_CONVERGED: iterations", "NEWTON_ITERATION_BUDGET"),
                ("unrecognized", "SETUP_OR_CUSTODY_FAILURE")]:
            self.assertEqual(termination(1, {"reason": reason}, True), expected)
        self.assertEqual(termination(0, {"converged": True}, True), "AUDIT_REQUIRED")

    def test_pending_packet_type_and_height(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "candidate.json"
            packet = dict(schema="WJH27DiagnosticCandidateV1", role="STEADY_SEARCH",
                          physical_history=False, state={"geometry": {"hill_height": .003}})
            for custody in ("UNCOMMITTED_TRIAL", "LAST_ACCEPTED"):
                packet["iterate_custody"] = custody
                path.write_text(json.dumps(packet))
                if custody == "UNCOMMITTED_TRIAL":
                    with self.assertRaises(ValueError):
                        accepted_packet(path, .003)
                else:
                    self.assertEqual(accepted_packet(path, .003), packet)
            with self.assertRaises(ValueError):
                accepted_packet(path, .004)


if __name__ == "__main__":
    unittest.main(verbosity=2)
