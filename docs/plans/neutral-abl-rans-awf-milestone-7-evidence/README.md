# Milestone 7 delivery evidence

- `validation.json`: original Linux delivery run, with eight CTests, 15 solver
  cases, five complete-field continuation comparisons, commands, binary hashes,
  and its source/documentation snapshot.
- `preparation.json`: six public five-flow preparations with Python site
  packages disabled; reference specifications match M6 and no long solver runs
  are repeated.
- `delivery.json`: final code/documentation hashes, repeated documentation-link
  check, manual build/render evidence, toolchain and links to the original run.

The original integration report is retained unchanged. Some documentation was
polished after that run; `delivery.json` identifies the final documentation and
checks that the tested numerical sources, helper scripts, CMake configuration
and validation code remain unchanged. Full PDFs, review PNGs, fields, history
payloads and logs remain in the recorded Linux run directory.

These short delivery tests do not extend M6's accuracy or stability claims.
The accepted long-run degree is four; the degree 5 full-3D counterexample and
other application limits remain in the [handoff](../../../solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md)
and [closure report](../neutral-abl-rans-awf-milestone-7.md).
