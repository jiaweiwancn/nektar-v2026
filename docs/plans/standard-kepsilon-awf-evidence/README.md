# Standard k–epsilon equilibrium-AWF evidence

The [derivation and report](../standard-kepsilon-equilibrium-awf.md) define the
model, benchmark, numerical limits and reproduction commands. These small
artifacts summarize authoritative Linux execution; field files, native history
payloads, raw CSVs and logs remain outside the source checkout under
`/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/standard-awf-20261031`.

- `validation-checks.json`: nine CTests, 20 workflow/control cases (17 accepted
  solver workflows and three expected preflight rejections), 48 independent
  standard-wall quadrature cases, configured constants, switch/fallback
  recovery, MPI continuation, and existing kernel/CLS regressions.
- `validation-study.json`: twelve long study cases. Five reference cases
  pass the engineering target; the six-layer target failure and the
  intentionally perturbed initial state retain their distinct outcomes.
- `preparation.json`: all seven public model choices in both dimensions,
  prepared with Python site packages disabled; no solver runs in this file.
- `delivery.json`: tested source/binary fingerprints, toolchain, compiled
  manual/render fingerprints, conservative budgets, primary-checkout state,
  and retained unsuccessful attempts and explanatory controls.
- Three PNGs: measured five-flow errors, fixed-H h/p refinement, and downstream
  profile errors. The plotting script consumes the raw native CSVs.

A `passed` delivery status describes the explicit acceptance gates. It does
not turn every coarse resolution, timestep, wall policy or disturbance into a
preservation/stability pass. In particular, no long-time 3D p-convergence or
terrain/thermal claim is made. SHA-256 fingerprints bind results to the tested
files; the parent revision plus recorded changes identifies the tested patch.

The collector can regenerate this packet with
`collect_standard_awf.py --run RUN_ROOT --output EVIDENCE_DIRECTORY` after the
checks, study, figures and manual compilation complete. Retained development
attempts are included when present; a fresh reproduction need not repeat them.
