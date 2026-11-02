# Neutral ABL RANS kernel and wall-layer tests

Start with the [runnable user guide](../../Examples/NeutralABL/README.md) for
model selection, complete sessions, boundary signs, output and restart.
The [Milestone 7 handoff](../../../../docs/plans/neutral-abl-rans-awf-milestone-7.md)
records the delivery revision and checks. The sections below retain the
development studies and their original commands.

| Entry point | Scope |
| --- | --- |
| `validate_cls_epawf_m0.py --output FRESH_DIRECTORY` | Linux source-bound contract and retained exact-base provenance audit; inventories baseline tests without a new numerical solve |
| `validate_cls_epawf_m1.py --kernel-probe BUILD/RANSKernelProbe --constitutive-probe BUILD/RANSConstitutiveProbe --output FRESH_DIRECTORY` | Revised-note M1: independent original CLS/AWF tensors, quadrature, both invariant-branch tangents and documented adverse diffusion; Linux only, NumPy/SciPy required |
| `validate_cls_epawf_m2.py --probe BUILD/RANSEquilibriumProfilesProbe --output FRESH_DIRECTORY` | Revised-note M2: native spectral ray loads versus analytic metric budgets/fluxes and independent quadrature, all-mode forcing identities and the specified sampled-K positivity failure; Linux only, NumPy/SciPy required |
| `validate_cls_epawf_m3.py --build BUILD --output FRESH_DIRECTORY` | Revised-note M3: native multiregion CG layout, pressure/interface rank, independent full transport residuals, steady AWF h/p refinement, evolving manufactured stages and full/Schur comparisons; Linux only, NumPy/SciPy/SymPy required |
| `validate_standard_awf.py --build ... --output ... --jobs 2` | Standard-model AWF kernel/coupling checks, five-flow preservation and h/p/dt refinement; see the [derivation and evidence](../../../../docs/plans/standard-kepsilon-equilibrium-awf.md) |
| `validate_m7.py --build ... --output ...` | Public generator/restart workflows, eight CTests, MPI continuation and documentation links |
| `run_m6.sh BUILD RUN [JOBS]` | Complete preservation/robustness study; failed controls are recorded explicitly |
| `validate_m3.py`, `validate_m4.py`, `validate_m5.py` | Standard transport, original CLS/AWF, corrected CLS/AWF acceptance |
| `validate_m1.py`, `validate_m2.py` | Strip prototype and independent kernel/oracle checks |
| `cao_tamura_calibration.py`, `cao_tamura_campaign.py` | Native development/audits plus legacy v1 original-first design/selection/freeze; adapt the primary route before a v2 search |
| `validate_cao_tamura_cv3.py`, `validate_rans_implicit_advection.py` | Native experimental inlet/startup/history checks and independent manufactured transport order |
| `validate_rans_stress_tangent.py --audit ...` | Full native transverse stress tangents versus an independent parallel-shear formula, including both CLS max-invariant sides |
| `record_cao_tamura_cv3.py` | Read-only incomplete-attempt evidence, conditional standard profiles and local tangent figures |

The [revised wall-stage guide](../../Examples/RevisedEquilibriumWall/README.md)
documents `RANSEquilibriumWallSystem`, its single-stage session functions,
separate native field outputs and serial affine-mesh scope. The M3 native,
stage-smoke and complete independent-oracle tests are registered in both the
production and external CMake builds. The oracle checks the frozen physical
norm scales and retains raw fields, rejected controls and refinement data in
its Linux output directory. `--suite residuals` is a shorter diagnostic and
cannot close M3. Native IMEX stepping and full-history restart are M4 work.

The [CV3 v1 attempt record](../../../../docs/plans/cao-tamura-cv3-calibration-attempt.md)
retains the full original-first design and numerical rejections. The active
[protocol v2](../../../../docs/plans/cao-tamura-cv3-calibration-protocol.md)
makes composite EP/transition/original floor calibration primary. Next is
CV3-R1 reference recovery and four composite feasibility cases; its execution
is pending. All-original success is not a primary EP prerequisite. The v1
campaign/selector/freezer remain legacy until adapted, and do not implement
the revised policy-specific freeze contract. Successful short solves or
conditional standard-only fits cannot unlock reserved flat or hill scores.
The external test build registers the calibration guards and stress-tangent
oracle in addition to the existing kernel/solver checks. Figures require
Matplotlib; native preparation requires SciPy and manufactured generators
require SymPy. Symbolic fixture construction is serialized because SymPy's
global assumption caches are not thread-safe; native solver processes can
still run concurrently within an explicit project allocation.

Use separate Linux run directories and at most 16 project workers. The CTest
fixtures are copied into the build's `rans-smoke` directory, keeping generated
fields, checkpoints and sidecars out of the source checkout. Nine smoke/kernel
tests do not replace the longer studies. Degree 4 is the accepted long-run
configuration; the degree 5 full-3D failure remains a documented limitation.

The standard-model extension adds `StandardAWFEquilibrium.xml` and its
`RANSStandardAWFEquilibriumSmoke` test. The original M1--M7 scripts retain
their historical scope and hashes. Use `validate_standard_awf.py` for the
extended implementation. `--suite checks` runs short kernel, coupling,
fallback and continuation tests; `--suite study` runs the long preservation
and convergence cases; the default runs both. `plot_standard_awf.py` renders
the native CSV measurements, and `collect_standard_awf.py` collects compact
reports, source/executable hashes and retained attempts. Use separate output directories
for concurrent suites and keep their combined allocation within 16 workers.

Milestone 5 adds `RANSWallModel=AWFEquilibrium` for CLS. Fixed weighting
isolates the dissipation correction; switched weighting uses seven independent
resolved indicators and explicitly supplied `RANSAWFOn`/`RANSAWFOff`
thresholds. The current configured CLS reference, cutoff admissibility,
original-form fallback and strict preflight are described in the
[Milestone 5 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-5.md).
Use `AWFEquilibrium` Robin tags for tangential velocity and k. The small
`AWFEquilibrium.xml` fixture is registered as `RANSAWFEquilibriumSmoke` in
the production and external builds. Run the complete Linux study with:

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m5.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m5"
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/plot_m5.py \
  "$RANS_RUN/validation-m5/validation.json" "$RANS_RUN/figures-m5"
```

The study audits actual native stage residuals, full fields, surface budgets,
and per-point switch/source/traction receipts against independent quadrature.
Its `--pilot` option cannot close the milestone. Cases are reused only when
the exact session and both solver/audit executable fingerprints match.
All numerical threads are one and cases run sequentially, using at most
four MPI ranks. Raw point CSVs and case receipts remain under the Linux run
directory; the compact validation manifest records their hashes and summaries.
No default switch thresholds imply terrain calibration. Milestone 6 supplies
the inlet/outlet, restart and longer preservation studies described below.

Milestone 4 adds `TurbulenceModel=CLSHighRe` with primary dissipation
`eps_tilde`, the full cubic covariance and algebraic Yap source. The optional
`RANSWallModel=AWFOriginal` couples the original wall function through the
conservative CV-to-CG strip. Its [specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-4.md)
defines the physical frame, forcing, strip inertia, normal velocity, algebraic
dissipation, ownership, and conservation diagnostics. `CLS_MMS.xml` and
`AWFOriginal.xml` are short CTest fixtures. With the isolated full executable
built as below, run the sequential Linux acceptance study with

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m4.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m4"
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/plot_m4.py \
  "$RANS_RUN/validation-m4/validation.json" "$RANS_RUN/figures-m4"
```

It uses at most four MPI ranks, one case at a time, with numerical threads set
to one. The generator uses SymPy for independent manufactured forcing and
SciPy quadrature for the original-AWF oracle. Original-AWF wall data use Robin
tangential/k entries tagged `AWFOriginal`, Dirichlet normal velocity and
positive initial epsilon data, and the coupled pressure marker. Epsilon and
normal traces are subsequently determined within the nonlinear stage.
The original baseline intentionally retains the note's wall-source imbalance.
The equilibrium correction and switch are added by Milestone 5.

Milestone 3 provides an opt-in `RANSVelocityCorrectionScheme` for standard
k–epsilon URANS. Its [specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-3.md)
defines the implicit nonlinear stage, mixed pressure space, boundary fluxes,
pressure convention, and scope. Build the complete executable in isolation
by adding `-DRANS_BUILD_SOLVER=ON` to the CMake command below. It compiles the
production incompressible-solver source list against the installed native
libraries and preserves the shared installation used by other experiments.

Use `StandardRANS_MMS.xml` for a short manufactured smoke test, registered as
`RANSStandardTransportSmoke` in both the main and external builds. The full
Milestone 3 study is:

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m3.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/validation-m3"
```

Set `RANS_RUN` and TMPDIR as described below. The fixture generator needs SymPy;
the validator otherwise uses the Python standard library and the native
`RANSFieldAudit` executable. Temporal cases use at most three concurrent
one-rank jobs. That pool closes before sequential four-rank MPI cases start.
Affinity and library threads are limited before loading generator dependencies.
`--skip-long-abl` is a pilot option and cannot close the milestone. Completed
cases are reused only when their exact fixture and executable hashes match.

The accepted configuration uses full 2D/3D CG velocity and scalar fields, with
degree p−2 continuous pressure at the same quadrature points. All physical
pressure boundaries use the `RANSPressure` Neumann marker with prescribed normal
velocity. The lower standard wall uses tangential `RANSLogLaw` Neumann data,
zero normal Dirichlet velocity, zero k flux and algebraic epsilon Dirichlet
data at a positive standoff height. See the supplied fixture and generator for
complete XML, including explicit quadrature settings. Saved `p` means reduced
pressure, `pi = p_mechanical/rho + 2k/3`.

Milestone 1 selects a control-volume wall strip on `0 < z < H`, coupled to
native continuous spectral elements on `H <= z <= Z`. The mathematical
contract and completion evidence are in
[the Milestone 1 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-1.md).

`RANSWallLayerProbe` checks the separate wall rows, dissipation constraint,
normal mass balance, interface signs, native outer weak residuals, h/p and
quadrature convergence, inconsistent-boundary controls, and a conservative
backward-Euler step with frozen reference viscosity. The outer residual now
calls the general CLS kernel at the native simple-shear gradients, including
its variable coefficient, full stress, and algebraic Yap term. This probe
isolates the wall strip; native coupled transport is covered by Milestone 4.

Milestone 2 adds `RANSKernelProbe` and an independent Python oracle. They test
standard k–epsilon and CLS tensors/sources, configured reference roots, original
AWF branches and forcing, corrected dissipation budgets and profile consistency,
admissibility, baseline recovery, and independent switch sensors. Full contracts
and measured results are in
[the Milestone 2 specification](../../../../docs/plans/neutral-abl-rans-awf-milestone-2.md).

The fixed reference is the note's rough-wall example with `P=1`, `H=2`, and an
additional matching sample at `N=4` metres. The CG intervals are `[2,4]`,
`[4,8]`, `[8,16]`, and `[16,30]` metres before uniform subdivision of each
interval. Geometry and sample locations do not change with polynomial order.
The N sample belongs to the element on its right, exercising communication
from rank 1 to the wall-owning rank 0 in the four-rank test.

The executables are registered as `RANSWallLayerPrototype` and
`RANSKernelUnitTests` in the main build when `NEKTAR_BUILD_TESTS=ON`.
An external build compiles the same production kernel sources into a static
library and links the installed native Nektar++ libraries, preserving the
shared build and installed binaries:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m2-reproduction
mkdir -p "$RANS_RUN/tmp"
export TMPDIR="$RANS_RUN/tmp"
taskset -c 112-127 cmake -S solvers/IncNavierStokesSolver/Tests/RANS -B "$RANS_RUN/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_BUILD_TYPE=Release
taskset -c 112-127 cmake --build "$RANS_RUN/build" --parallel 4
taskset -c 112-127 ctest --test-dir "$RANS_RUN/build" -R '^RANS(WallLayerPrototype|KernelUnitTests)$' \
  --output-on-failure --parallel 1
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m2.py \
  --build "$RANS_RUN/build" \
  --output "$RANS_RUN/validation" --mpi-ranks 4
```

Run from the synchronized Linux source root. Choose a fresh experiment directory;
the commands above use a separate reproduction directory. CPU IDs 112–127 are
the recorded allocation on this server; use up to 16 allowed IDs if its affinity
configuration changes. The project limit is **16 CPU workers in
total**, including concurrent build and test work. The validator serializes its
one-process and four-rank runs, sets library thread counts to one, and limits
its inherited Linux CPU affinity to at most 16 logical CPUs. Build and CTest
allocations must also share the 16-worker budget. The repository-wide limit
across other server projects still applies.

Outputs are the kernel report, serial/MPI column reports, their logs, and a
validation manifest with source hashes, commands, CPU allocation, thread
settings, and comparisons. Milestone 2's oracle requires NumPy and SciPy
(tested with 1.21.5 and 1.8.0). It sets BLAS/thread environment variables and
CPU affinity before importing those packages. The kernel probe itself has no
NumPy/SciPy dependency.

For the original standard-library-only column audit, use `validate_m1.py`
with `--probe "$RANS_RUN/build/RANSWallLayerProbe"`, `--output` and
`--mpi-ranks 4`. It now exercises the general CLS kernel through the column.
Keep historical Milestone 1 reports unchanged when running this regression.

A successful kernel/column result establishes the specified local model
contracts and column convergence. Milestone 3 adds standard-model unsteady
transport, Milestone 4 adds CLS and the original AWF, and Milestone 5 adds the
corrected-AWF coupling. Milestone 6 supplies the full analytical inlet/outlet
preservation study.

Milestone 6 adds `generate_m6.py` and `validate_m6.py` for full 2D/3D neutral
profiles, total-traction outlets, five-flow-through comparisons, h/p/dt and
physical-H studies, perturbations and complete IMEX history restart. The
optional `RANSLinearSolver=MatrixFree` keeps native static-condensation systems
as preconditioners and applies the current weak operator. `Native` is the
default. `RANSReferenceMonitor=1` records independent analytical errors without
changing the solution. `RANSFieldAudit session.xml a.fld b.fld profiles.csv`
adds common-height streamwise profiles. Derived viscosity, production,
dissipation and Reynolds stresses are available for both closures.

A full RANS restart retains its native field file/directory and the referenced
`.rans-<completed-step>` history directory. `RANSRestartMode=History` requires
complete history; `Auto` restores it when present and uses native startup for
legacy fields; `Cold` explicitly requests startup. `RANSHistoryPath` relocates
stored payloads. Mesh, quadrature, timestep, model, wall parameters, IMEX order,
boundaries and forcing must stay consistent. MPI repartitioning is supported.
Use step-based checkpoint intervals. The M6 specification distinguishes cold
startup from a continuation and analytical accuracy from numerical stationarity.

After the external build above with `-DRANS_BUILD_SOLVER=ON`, run separately:

```bash
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m6.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/verification" --suite pilots
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m6.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/verification" --suite restart
taskset -c 112-127 python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_m6.py \
  --build "$RANS_RUN/build" --output "$RANS_RUN/study" --suite study --jobs 2
```

The study records numerical failures and engineering-target failures explicitly;
its `recorded` status means the requested comparisons were collected, not that
every wall policy preserves the reference. Inspect individual cases and the
[Milestone 6 report](../../../../docs/plans/neutral-abl-rans-awf-milestone-6.md).

`ABLPreservation.xml` and `ABLPreservation3D.xml` are small outlet/checkpoint
and 3D perturbation smoke fixtures. The production and external builds register them as
`RANSABLPreservationSmoke` and `RANSABL3DPerturbationSmoke`. Generation of the
neutral fixtures requires only the Python standard library; the closed-form
original-AWF matching integral agrees with the independently recorded Linux
quadrature anchors at H=1, 2 and 4 m to less than 1e-12 m/s. The larger
manufactured/oracle suites still require SymPy/SciPy; plotting uses NumPy and
Matplotlib. The `columns` suite takes `--study-output` to distinguish a
stationary numerical-column copy from the analytical comparison.

`verify_outflow_m6.py --build ... --output ...` compares complete native and
matrix-free outlet fields for both models in 2D and four-rank 3D. After all
final phases finish, `collect_m6.py --run ... --output ...` checks the closure
gates and writes compact source/binary hashes, environment and result evidence.
The `robust` suite records the degree-five 3D refinement experiment, including
numerical or profile-target failures; those are not preservation passes. Its optional
`--extended-perturbation` flag repeats the already-tested five-flow 3D
perturbation at a smaller timestep; manufactured solutions supply the
independent temporal-order test.
The recorded degree-five four-rank case fails the target and is a documented
robustness limit; the five-flow acceptance resolution is degree four.
`--suite robust --record-existing-failure` reduces an already stopped trial
without rerunning it. It requires native failure evidence or a threshold-stop
receipt bound to the exact fixture and executable hashes.

For a complete fresh reproduction, activate the installed MPI environment,
build the current source externally with `RANS_BUILD_SOLVER=ON`, and run
`taskset -c 112-127 bash run_m6.sh BUILD_DIRECTORY RUN_DIRECTORY 2` from this
directory on Linux. Choose a separate run directory under the large-volume
run root. This sequential wrapper uses at most eight study ranks by default
and writes the status markers expected by the evidence collector. The
recorded milestone run overlapped independent phases within 16 total workers;
the wrapper provides the same numerical checks without that overlap.

`refresh_m6_audits.py --run ... --build ...` supports an audit-only rebuild
without repeating unchanged solver runs. It requires the solver hash to stay
identical, reruns the saved-field audits and comparisons, and requires identical
time and complete-field JSON values. It refreshes optional profile CSVs and
retains the original reports with an explicit provenance receipt. The M6 run
used this after correcting local-coordinate profile interpolation. A fresh
reproduction with the corrected sampler needs no refresh step.

The [TerrainHill example](../../Examples/TerrainHill/README.md) documents
`generate_terrain.py`, `terrain_initial.py` and `validate_terrain.py`.
Its `planes`, `periodic-curved`, `curved`, `transitions`, `preservation` and `pilot` suites
exercise vector wall coupling, smooth/C1 geometry, integrated strip continuity,
spatial properties, 1/2/4-rank sampling, history restart and paired width controls.
The full-hex CLS pilot uses eight ranks. `pilot-2d --publish-to FULL_PILOT`
can run on four other CPUs and publish completed cache aliases without replacing
existing cases; account for both runners within this project's 16-worker limit.
The external build supplies a full-precision source SpatialDomains library;
select `BUILD/geometry-lib` first in `LD_LIBRARY_PATH`. Use fresh run directories
under the Linux large volume and bounded CPU affinity. `record_terrain.py`
collects passed suites, source/binary hashes, MPI/thread settings and CTest
receipts. C1-foot maximum-error limitations remain explicit in those receipts.

CV2 adds `RANSObservationKernelProbe`/`RANSObservationKernels`, native
`RANSObservationAudit`, `cao_tamura_experiment.py` and
`validate_cao_tamura.py`. The public entry point and trial JSONs are in
`Examples/TerrainHill`. Build all external targets with `RANS_BUILD_SOLVER=ON`.
On Linux, run each suite in a fresh large-volume directory, for example:

```bash
taskset -c 1-4 python3 validate_cao_tamura.py --build "$CV_BUILD" --output "$CV_RUN/verification/probes" --cpus 1,2,3,4 --suite probes
```

Other suites are `controls`, `generator`, `workflow` and `utilities`. They cover smooth
anchors, finite viscosity, original recovery, curved continuity, 2D/3D native
probes and independent auxiliary quadrature, generation/preparation of S/R
candidates, physical-time stress statistics and experimental-input history.
Set all BLAS/OpenMP threads to one and account for concurrent jobs within the
16-worker project cap. Full checkpoints/logs remain outside Git. These suites
verify functionality, not a developed experimental calibration or stationarity.
`record_cao_tamura.py --help` collects these suites together with the canonical
five-clock preservation suite, legacy regressions and eleven CTests. It checks
the tested binaries/library before writing the compact CV2 evidence receipt.
