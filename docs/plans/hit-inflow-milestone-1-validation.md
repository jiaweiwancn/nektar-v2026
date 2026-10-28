# HIT inflow Milestone 1 closure record

Milestone 1 is closed on 2026-10-28. The
[interface specification](hit-inflow-milestone-1.md) fixes the version 1 XML
grammar, supported solver/integrator combinations, Fourier and RNG conventions,
normalization and mean-only behavior, spectrum format and explicit extensions,
startup/restart contract, and the minimal prescribed-velocity update hook.

## Delivered inputs and evidence

| Requirement | Deliverable and result |
| --- | --- |
| Unambiguous configuration | Two checked XML fragments for tabulated/Homogeneous1D and von Karman/full 3D, with sparse global inlet ID 7 |
| Supported combinations | Standard VCS, nonlinear incompressible equations, u/v/w/p, continuous projection, fixed-step blank-variant IMEX2 and its IMEX1 startup |
| RNG and normalization | Exact versioned unsigned-integer mixing, open-interval uniform conversion, Box-Muller convention, conjugate pairing, finite-lattice ensemble normalization, zero-amplitude limit |
| Benchmark source data | Complete principal Table 2(a), Table 3 and Table 4 numeric transcription checked against the original PDF, with printed pages and source PDF SHA-256 |
| SI inputs | Six station spectrum files, correct spectrum-kind headers, three-station target file, exact decimal scaling and preserved missing values |
| Reproducibility | Standard-library preparation/check tools and deterministic source/generator/output digests |
| Update timing | Source audit of initialization, PreIntegrate, GLM RHS timing, H-pressure acceleration and existing post-pressure classes; all-rank loading precedes local object construction |
| Validation | Mac preparation/source digest checks and successful single-worker Linux fixture/contract checks |

The files are in
`solvers/IncNavierStokesSolver/Tests/HITInflow/`. The XML files are contract
fragments with no mesh or complete lateral/outlet boundaries. The existing
solver still does not implement `HITInflow`; runtime acceptance and numerical
tests remain in later milestones.

## Tested snapshot and commands

The code audit and candidate files are based on Mac branch
`codex/hit-inflow-generator` at
`132f4582e438ead0d81f45da671bcd3b813e56b1`, plus the identified Milestone 1
files. The candidate snapshot manifest has SHA-256
`e1b92c4971e5859a0bae4719aa53dbe5d824995bfeb50e4bf03793b08fdbbee7`;
its archive has SHA-256
`98355c253655106926d7d8774de1b9bc833ac99ed1c2a8a1aaefa7f7a82c4b6e`.
The snapshot contains all 14 fixture/specification-support files and the two
plan/specification documents. The parent plan was updated after the checks to
record closure and align its guidance with the specification. This closure
record was also added afterward. The checked fixture and interface files are
unchanged. The [unaltered Linux receipt](hit-inflow-milestone-1-validation.json)
records every tested file digest and command.

The Linux source checkout remained clean on `main` at
`a3b7efc8fef7e075390f0f59457dfdac30f0b6ed`. The checks ran against the separately
synchronized candidate under
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/milestone-1/e1b92c4971e5`.
The original checkout and build were not modified.

From that snapshot root, with Python 3.10.12:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py --check
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/check_contract.py --json
```

Allocation was one serial Python worker and no MPI ranks. No active compiler,
solver, or test-runner jobs were found during the worker-budget check. The
checks ran sequentially. There was no Nektar++ compilation, solver execution,
or numerical/performance conclusion at this interface/input milestone.

Results: six spectra and one station table reproduced exactly, all input and
output hashes matched, two XML examples and the mean-only limit passed, and 26
invalid configurations were rejected. The generated data-manifest SHA-256 is
`6eac80b5db4308d3b58b750a0fe61c62653fd7715e183ddab2bfbaa3b40b88b1`.
The Mac also verified the source PDF digest, Python syntax, XML well-formedness,
documentation structure, and repository-relative documentation links.

## Next milestone

Milestone 2 implements the spectrum providers and deterministic transverse mode
generation against the frozen contracts. It must add RNG known-answer checks,
Gaussian covariance checks, represented-energy diagnostics and analytical
derivative/divergence checks. This closure does not advance that implementation
or mark its tests complete.
