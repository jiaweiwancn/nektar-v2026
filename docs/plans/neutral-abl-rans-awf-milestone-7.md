# Neutral-ABL URANS documentation and delivery

Milestone 7 consolidates the implemented models, their operational configuration
and measured limits into a runnable handoff. It does not change the closure,
wall equations, pressure algorithm or time integration, and does not claim to
resolve the degree 5 full-3D growing mode found in Milestone 6.

| Delivery record | Value |
| --- | --- |
| Status | CLOSED on 2026-10-30 |
| Source parent | `0dc0b26dd94c84fe22abcbfe63bbcb714a6781aa` |
| Branch | `codex/neutral-abl-rans-awf` |
| Linux run root | `/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m7-delivery` |
| Isolation | Fresh detached worktree and external build; primary checkout/build/install preserved |
| Project limit | At most 16 workers; build/cases use four, manual compilation uses one additional worker |

## Delivered entry points

The [user handoff](../../solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md)
provides a fresh build, short serial/full-3D runs, all four model/wall choices,
configuration defaults versus benchmark settings, boundary signs, field/diagnostic
interpretation, complete checkpoint continuation, relocated history and the
longer-study commands. It depends on tracked source files, not the untracked
literature working copies.

`Examples/NeutralABL/generate.py` prepares complete neutral sessions through
the existing verified generator. `restart.py` prepares additional steps from
a complete field/checkpoint while retaining model, mesh, quadrature and dt.
Both tools reject non-XML output paths. The delivery check confirms that an
attempt to use a saved field as output is rejected without changing that field.
The examples are indexed by the repository README, solver examples and RANS
test index. The user manual corrects obsolete outlet guidance and states the
accepted resolution and unsupported cases. The developer manual describes
source ownership, the coupled stage, conservative strip, matrix-free action,
traction pressure convention, history format and extension evidence.

The six transport smoke fixtures are staged into the build's `rans-smoke`
directory in both production and external CMake configurations. Their generated
data therefore remain outside the source tree.

## Acceptance and retained limits

The Linux delivery check builds from the source parent plus this delivery
patch, runs all eight CTests, and exercises the public preparation and continuation
commands. It checks standard/original/fixed/switched 2D startup, full-3D startup,
serial-to-four-rank and four-to-two-rank continuation, a step checkpoint,
relocated history and explicit cold startup. Full-history fields are compared
against uninterrupted runs. Documentation links and absence of generated
solver output in the source checkout are checked too.

Numerical C++/header hashes are compared with the M6 evidence. With those
sources unchanged, the longer M1–M6 results remain the numerical baseline;
short delivery checks do not replace or extend the five-flow accuracy claim.
The accepted degree 4 standard/fixed/switched reference cases meet 0.1% in
2D/full 3D. Original AWF drift, inconsistent controls, the degree 5 full-3D
failure, span-independent long 3D data, unavailable inlet-reaction budgets
and caller-managed boundary/forcing consistency at restart remain explicit.

Terrain geometry, general 3D disturbance stability, curved/backflow outlets,
Homogeneous1D, thermal/stratified ABL and low-Re CLS are not delivered here.
Further numerical stability work is needed before applying the failing higher
degree configuration. The seven planned milestones end with this measured
flat neutral baseline.

## Linux validation and ready-to-run binary

The fresh Release build succeeds with four compiler workers. All eight CTests
and 15 public-workflow solver cases pass. The latter include all four 2D policies,
standard/perturbed CLS full-3D startup, uninterrupted reference runs, full-history
continuation, first interval-checkpoint recovery, relocated sidecars, explicit
cold startup and full-3D four-to-two-rank continuation. Solver cases run
sequentially with at most four ranks, with OMP/OpenBLAS/MKL/BLIS threads one.

| Complete-field continuation comparison | Largest absolute difference |
| --- | ---: |
| Standard, serial prefix to four ranks | 1.16e-12 |
| CLS, serial prefix to four ranks | 1.91e-11 |
| CLS first interval checkpoint | 1.95e-11 |
| CLS relocated history | 1.20e-11 |
| Full 3D CLS, four ranks to two | 2.53e-9 |

These are maximum absolute differences across the primary and projected derived
fields at the same time. Limits are 1e-7 for the 2D comparisons and 1e-6 for
the full-3D comparison. Cold startup is checked as an explicit different workflow,
without asserting uninterrupted-run equivalence.

All 127 recorded numerical C++/header sources match the M6 baseline. The six
standard/fixed/switched 2D/3D five-flow setups prepared by the public CLI have
reference specifications identical to M6. These preparation checks run with
Python site packages disabled and execute no long solver runs. No generated
solver output is found in the fresh source checkout.

The usable solver is
`m7-delivery/build/IncNavierStokesSolver` under the Linux run root above.
Activate `/home/jay/.local/opt/nektar-v2026-mpi/env.sh` before using it.
The shared primary checkout remains at
`da607a9a46d6b11d7b30e67441c878139e11adf2`; its build, installation and
ongoing HIT jobs were preserved.

Commands used after synchronizing the reviewed delivery patch:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m7-delivery
export TMPDIR="$RANS_RUN/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$RANS_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$RANS_RUN/build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
taskset -c 112-127 cmake --build "$RANS_RUN/build" --parallel 4
taskset -c 112-127 python3 \
  "$RANS_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_m7.py" \
  --build "$RANS_RUN/build" --output "$RANS_RUN/acceptance"
```

Choose a new run/output directory for reproduction. CPU IDs 112–127 are this
server's recorded allocation; the driver uses only four solver workers within
it. The external CMake build applies strict warning/error options to the RANS
sources and the field audit. Tested toolchain and source/binary hashes are in
the [delivery evidence](neutral-abl-rans-awf-milestone-7-evidence/).

## Documentation review and evidence

Both complete manuals compile with Linux pdfLaTeX/latexmk. Their final RANS
pages are rendered and visually checked: user-manual PDF pages 226–231 and
developer-manual pages 137–140. Variable-order tables, long option names,
heading separation and command wrapping are corrected. Legacy font commands
in the incompressible chapter are replaced by supported math fonts; no equations
or model constants change. The complete books still report unrelated unresolved
reference/citation and duplicate-label warnings; compilation succeeds and the
reviewed RANS pages contain no clipped or overlapping content.

The evidence retains the original 15-case validation report and its hashes,
the six setup-preparation records, and a final delivery manifest. Documentation
links and hashes are checked again after the final layout edits; the manifest
distinguishes final documentation from the earlier test-run snapshot. Numerical
sources, helper scripts and validation code are unchanged after their successful
run. Large fields, histories, logs, full PDFs and rendered review pages remain
under the Linux run root, outside the source checkout.

The [project plan](neutral-abl-rans-awf-proposal.md) now records all seven
milestones closed. The user handoff and source branch are the delivery; the
known higher-degree and geometry/physics limitations remain part of that handoff.
