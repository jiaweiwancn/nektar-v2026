#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check R3 worker admission and the reservation used for resuming a paused job."""

import unittest
from unittest.mock import patch

from task2r_run import worker_census, workers


class WorkerCensusTest(unittest.TestCase):
    @patch("task2r_run.subprocess.check_output", return_value="""+101 101 100 R - simpleFoam
101 102 100 S do_pol simpleFoam
103 103 100 T do_sig simpleFoam
103 104 100 T do_sig simpleFoam
105 105 100 t ptrace python3
106 106 100 S do_sel mpirun
106 107 100 S do_sel mpirun
108 108 100 D wait_io unknown
109 109 100 S do_poll unrelated
110 110 100 Z - python3
111 111 100 R - kworker/0:0
""")
    def test_explicitly_stopped_workers_leave_active_capacity(self, _):
        census = worker_census()
        self.assertEqual(census.active, 3)
        self.assertEqual(census.stopped, 3)
        self.assertEqual(census.parked_helpers, 2)
        self.assertEqual(workers(), (6, 2))

    @patch("task2r_run.subprocess.check_output", return_value="""+201 201 200 R - simpleFoam
203 203 200 T do_sig simpleFoam
204 204 204 T do_sig IncNavierStokes
205 205 205 S do_sel mpirun
206 206 205 T do_sig orted
""")
    def test_resume_reserves_the_live_and_stopped_owned_footprint(self, _):
        census = worker_census({204, 205})
        self.assertEqual(census.active, 2)
        self.assertEqual(census.owned, 3)
        self.assertEqual(census.stopped_owned, 2)
        self.assertEqual(census.projected_active(4), 5)
        self.assertEqual(census.projected_active(2), 4)


if __name__ == "__main__":
    unittest.main()
