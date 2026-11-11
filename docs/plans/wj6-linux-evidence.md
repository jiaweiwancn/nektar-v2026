# WJ-6 resource-qualified uniform L2/Nx2 study

The physical problem remains the finite reference-Dirichlet strip
`[0,1] x [0,2]`, with the original WJ tensor, cubic-N surrogate,
DevStrainRawWorkV1, direct epsilon, independent Cd scalar diffusion,
full-strain Linear inner stress and uncapped EP/budget sources. The
specifications, original WJ-4/WJ-5 outcomes and their sealed accounting
are unchanged. Stationarity is not targeted; WJ-4's NOT_REACHED remains.

## Source and execution identities

Reviewed parent: `6ee8f356d0cb31d9b99f86127d676ffa408d34e1`.
The fetched branch initially matched it. Other Mac/Linux worktrees,
dirty user work and live campaigns were preserved.

The native producer was rebuilt from
`30ffa4000353d8055d19abec3c7dc189b7004c4a` in isolated `build-v2`.
Its Git tree is `464e2e82af01982e00ffcad6414f4f38e6423c55`.
The frozen numerical driver/observer is
`4f7520aabee2604efd5b0194a5fa8592d672816a` (`source-v7`, tree
`73d4302031fd5760792a839f9d614e072ba61ae8`), with all 1800
C++/header/CMake inputs byte-identical to the native build. Registered
controls use `source-v4` (`569e3f1`) and verified aliases to those same
rebuilt executables. Later resource/publication code is a separate passive
lineage; it does not rewrite a producing source path or observation.

The exact original/worktree specification hashes agree:

```
main     1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124
appendix 0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c
```

The main-note literal `quad`/`\quad` typography erratum is retained;
the intended logarithmic coefficient remains `(d+a*h)/a^2`.

Linux heavy root:
`/media/jay/新加卷/nektar-v2026-runs/wj6-20261111`.
The documented Tailscale fallback was used after the USB alias timed out.
The installed dependency environment was read-only. GNU 11.4 Release builds
retain `-ffp-contract=off -fno-fast-math`, `RANS_BUILD_SOLVER=ON`;
the manifest binds effective flags, ELF hashes and resolved dependency hashes.
Every numerical invocation uses one serial rank and all six numerical-library
thread controls set to 1. Current allocation checks include MPI helper threads
and reserve 8 workers for each native/observation launch under the 72 ceiling.

## Capacity implementation and controls

`RANSWallDenseEvolutionCapacity` is an execution-only SolverInfo setting,
default 3000, exact finite integer range 3000..4331. Extended capacity requires
serial flat 2D WJ/PimpleSingleSweep, eta=0, uncondensed evolution, and a 64-hex
`RANSWallDenseEvolutionResourceDesign`. The same effective policy guards
Patch::Solve and the WJ one-sweep entry. No constructor size guard or physical
signature exception was introduced. LAPACK block solves and equations are
unchanged. Reports bind requested/effective capacity, resource design and
actual block dimensions.

| Mesh | State | Free velocity | Pressure+gauge | Flow | Free k | Epsilon solve |
|---|---:|---:|---:|---:|---:|---:|
| uniform L2/Nx2 |4331|1790|187|1977|900|242|
| M2/Nx2 control |2975|1250|133|1383|630|122|
| uniform L2/Nx1 |2493|716|125|841|360|98|

Both compiled capacity controls confirm every dimension. They reject default
oversize evolution, fractional/nonfinite/out-of-range settings, wrong routes
and missing resource identity; preserve large reference/readback use; and
establish exact small-case off/on physical signatures, trajectories, gauge,
histories, audits, raw/corrected k and one-per-block solve/bound counts.

All 35 unique WJ-5 selections passed, including the shared legacy transaction.
The three new capacity selections also passed: 38 unique CTests from 62
registrations. A separately frozen final metadata invocation passed six
fixtures, including missing-M1, incomplete primary/extension, asymmetric
uncertainty sums and immutable fresh-startup evidence. Versioned passive
metadata suites subsequently passed seven, eight and nine fixtures, including
all-order coverage, zero-launch transaction resumption and verified partial
observation reuse. These are additional fixtures, not additional unique CTests.
Exact inventories,
commands, elapsed times and log hashes are sealed with the evidence.

The first single-worker native build exposed a parser-local name collision
and exited 2 after 291.729956s. The corrected four-worker build completed in
789.856748s; both original logs are retained. Build time is separate from
the per-policy numerical ledger.

## Frozen study and measurement rules

Exact `normal_knots(2)` and x knots `[0,1]`/`[0,.5,1]` are bound in every
design, XML and requested intent. All three trajectories use fresh same-runtime
native CP0s. The dt pair requires identical initial physical state/gauge and
corresponding physical histories, with distinct native dt signatures/time labels.
No old Nx1 callback is fabricated or reused.

Pilot observations are coarse 4/fine 8/Nx1 baseline 4 at `.000128s`.
The Tx targets are coarse 32/fine 64/Nx1 baseline 32 at `.001024s`.
Actual intermediate observations include coarse 16/24 and fine 32/48, retaining
genuine CP24 for the eight-step transaction. Startup and all observations
remain gates; a later passing endpoint cannot erase an earlier failure.
Only complete passing Tx/transaction/resource gates for BOTH policies can
authorize continuing those same histories to T3 `.004096s`. No T4 is authorized.

Physical observations use composite actual-mesh/physical/reporting partitions,
normal 8/10 and horizontal 4/6, explicit face/intersection owners and zero-weight
face duplicates. The independent quadrature norm floor stays 1e-10. Separate
tables retain uncancelled spatial residuals, actual producing-stage defects,
total continuum error, represented-initial drift and same-time comparisons.
There is one full-domain pressure gauge alignment and the entire outer epsilon
domain. Independent evolved loads, signed interface modes, pressure-test
projection/complement, element modes and displaced-maximum neighborhoods are
passive; they add no governing solve or correction.

Each trajectory's own 100x arithmetic estimate is applied to all actual row maps,
then BOTH contributions are added for a pair. The boundary estimator retains
prescribed-data/lift/basis/geometry/accumulation reserves, separate k/epsilon
units, dependent epsilon action and unconstrained pressure/gauge treatment.
These finite frozen-system envelopes are not nonlinear stability bounds.

The unchanged field/interface/divergence budgets and 10% sensitivity limits
are enforced over every field, norm and scope. Numerical correction must be
inactive and raw/corrected k identical. A hard failure stops before commit;
engineering failures permit only the already frozen adverse Tx pair.

The ledger caps stay 800 charged units and 10800 numerical seconds per policy.
Native/observation/SVD processes have finite time, 8GiB RSS and output ceilings;
each policy has a 32GiB output ceiling. The measured resource observer additionally
charges parent algebra, verification and I/O outside child processes, with 25%
headroom and separate 600s closure plus nine-attempt transaction reserves.

## Results and reproducibility

**Primary Tx qualification: PASS for both policies. Conditional T3 extension:
RESOURCE_LIMIT; no T3 governing step was launched.** All 128 core steps per
policy completed, including the32-step fresh Nx1 baseline. Each actual-capacity
transaction added exactly nine attempts: one rejected and eight accepted.
The final packet records 170 charged units per policy, comprising137 study/
transaction attempts and 33 regression/control units. No material correction
occurred; every raw/corrected k snapshot is bitwise equal.

The closed ledgers charge 5030.391243 numerical-study seconds for EP and
5875.989613 seconds for budget, including failed work, shared controls,
allocation waits, publication and conservative closure/I/O reserves. Every
ledger reservation is reconciled or conservatively charged. Build time totals
1081.586704 seconds separately.

| Required component | EP | Budget |
|---|---|---|
| Capacity/control verification | PASS | PASS |
| Startup, pilot and Tx accuracy/admissibility | PASS | PASS |
| Horizontal/timestep sensitivity and quadrature | PASS | PASS |
| Actual-capacity CP24 transaction | PASS | PASS |
| Tx resource ceilings and closed accounting | PASS | PASS |
| Conditional T3 admission | RESOURCE_LIMIT | RESOURCE_LIMIT |
| T3 evolution | NOT_RUN | NOT_RUN |

The coarse and fine endpoints are the actual binary64 times
`.0010240000000000002` and `.0010239999999999987` seconds. Their dimensionless
duration is approximately `Tx*u*/H=.000512`; coarse/fine `dt*u*/H` are
`1.6e-5`/`8e-6`. No checkpoint signature, native history or time label was edited.

| Normalized sampled maximum divergence at Tx | EP | Budget |
|---|---:|---:|
| L2/Nx1, coarse32 |8.03740391195e-5|8.03740670365e-5|
| L2/Nx2, coarse32 |2.39985197375e-4|2.39985239742e-4|
| L2/Nx2, fine64 |2.40476976026e-4|2.40477018949e-4|
| Fixed limit |1e-3|1e-3|

The coarse Nx2 maximum is at `(x,z)=(1,.1)` on geometry/element 55, the
upper normal face of the cell below h. EP has
`Gxx=.000119992598687662`, `Gzz=-3.1646896260513148e-21` s^-1.
The central-strip maximum is 1.19295646027e-4 normalized. The inner-region
divergence L2 norm is 5.77106167951e-6, its continuous-Q1 projection is
4.32052396012e-16, and the complement retains 5.77106167951e-6.
The actual continuity-row maximum is 1.47491192374e-17. Budget gives the
same pattern. Element-mode and signed-load packets retain all owners and the
distinction between the global continuous-Q1 complement and elementwise degrees.

Uniform normal refinement reduces the retained L1/Nx2 maximum by about 8.87x
under unchanged equations, spaces and boundaries. This supports normal
under-resolution as a remedy for this finite-strip obstruction. It does not
identify an unstable pressure pair or uniquely attribute the error to roughness
forcing. The exact reference has zero F_x and smooth viscosity near h;
evolved forcing and its components remain measured separately.

All 504 unique sensitivity leaves per policy pass: two observation orders,
two comparison groups, three common times, seven fields, three scopes and two
norms. The largest `(difference+both uncertainties)/limit` is .0931194 EP
and .0931388 budget. The limiting horizontal leaf is the full velocity sampled
maximum at Tx: approximately 8.39870e-6 difference plus 9.15e-7 uncertainty,
against the 1e-4 sensitivity limit. The largest temporal fraction is about .0274
at the pilot; its difference is below its conservative uncertainty. This is
budget-qualified timestep sensitivity, not a new temporal-order measurement.

Worst normalized errors/drifts over all recorded times, trajectories, both
observation orders and reporting scopes are:

| Quantity | EP | Budget | Limit |
|---|---:|---:|---:|
| Velocity L2 error |2.63934e-6|2.63934e-6|1e-4|
| Velocity L2 drift |1.92235e-6|1.92235e-6|1e-4|
| Velocity sampled error |1.75895e-5|1.75895e-5|1e-3|
| Velocity sampled drift |1.99497e-5|1.99497e-5|1e-3|
| k sampled error |4.16409e-7|6.86260e-6|1e-3|
| k sampled drift |4.16409e-7|7.54708e-6|1e-3|
| Entire outer epsilon sampled drift |1.14146e-7|1.14146e-7|1e-3|
| Mechanical-pressure sampled error |5.32774e-6|1.73701e-5|1e-3|

At coarse Tx the normalized full-traction jump is 3.11152e-5 or less,
k-flux jump 2.34671e-7 or less, and reference pressure-offset departure
2.75105e-6 or less. Shared k traces are exact; epsilon trace discrepancies
are at most 3.34e-16 in those endpoint packets. Full tables retain every
observation, the sampled K/cutoff ranges, independent evolved source/diffusion
loads and signed weak interface modes. Separately normed residual pieces are
not added as norms. Startup traction/flux identities partly imposed by
RadialProject remain labeled separately from evolved measurements.

Both CP24 transactions pass exact all-nine-state/audit and all-eight callback/
incoming-history/count comparisons, final primary/gauge/GLM/time/sample equality,
all retained k substages and frozen-k matrix identities. The forced step25
candidate matches the reference's actual accepted step25, is retained separately,
and rolls back to CP24. Original checkpoint/report bytes and the private input
copy remain unchanged. Frozen readback performs no Lift or governing solve.
CP24's preceding completed callback remains bound to its genuine producer;
the zero-step startup packet is not substituted for it.

## Resource decision and retained repairs

Pilot whole-process times, including observations/SVD/I/O, were 616.521400s EP
and 753.982132s budget. Actual first four-step native blocks took 46.123509s
and 77.583734s. Completed original Tx continuation processes took1739.380333s
and 2339.729625s. Maximum whole-phase RSS was 5815268KiB, below 8GiB.
Final per-policy artifact inventories are 5,346,213,564 bytes EP and
5,369,461,605 bytes budget, below the frozen 32GiB ceiling.
The all-order passive supplement used 675.212378s EP; budget's failed partial
and verified-resume processes used 248.804436s and 563.584010s. These times,
allocation waits and failed-publication I/O remain charged.

The frozen T3 projections, taken after Tx and transactions, were:

| Policy | Already charged seconds at projection | Remaining seconds including reserves | Projected total | Ceiling |
|---|---:|---:|---:|---:|
| EP |4305.375577|6830.795986|11136.171563|10800|
| Budget |5139.550764|9934.579487|15074.130252|10800|

Each projection retains 384 remaining native attempts, actual mesh-specific peak
step rates with 25% headroom, observation/SVD/parent-I/O costs, the all-order
supplement cost and 600s closure reserve. The joint write-once admission is
RESOURCE_LIMIT. T3=.004096s remains the frozen unreached target. No duration,
timestep, mesh, tolerance or budget was changed to obtain a pass.

Two worker denials occurred before Popen: EP's initial transaction readback
and budget's seventh supplementary readback. Both have zero native attempts.
The EP transaction resumed in a fresh process/output name. Budget's completed
six readbacks were verified against intent/argv, inputs, outputs, runtime and
unchanged core-observer hashes, then reused in a new versioned qualification;
only eight missing readbacks were executed. Original failed observers and
receipts remain immutable. No rejected physical state was retried in process.

The original higher-order quadrature packet did not export native constraint
maps. The versioned supplement exports those maps, verifies exact equality of
the original physical samples, and requires both order 8 and order 10 accuracy/
sensitivity gates. Its worst basis/evaluation and separate coefficient envelopes
cover both point sets, retaining the 100x reserve and pressure-gauge treatment.
Latest metadata fixtures explicitly reject missing order coverage, failed-cache
reuse and confusing a physical audit failure with a historical resource denial.

Publication-only failures (duplicate exit-code column, then using a full-file
identity validator on a native hash-only descriptor) are retained with their
original outputs/logs. The native descriptor's SHA-256 was independently checked
and matched. The repair retains hash verification plus exact binary matrix
shape/size checks; no native source, observation or tolerance changed.

## Commands, receipts and next scope

The exact commands, exits, UTC times, log hashes, XML/argv/resource designs,
runtime/dependency closure and source-content hashes are in
`wj6-linux/final-command-receipts.json` and its manifest/transport bridges.
Commands used the frozen paths under the Linux root, for example:

```sh
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cmake -S /media/jay/新加卷/nektar-v2026-runs/wj6-20261111/source-v2/solvers/IncNavierStokesSolver/Tests/RANS \
  -B /media/jay/新加卷/nektar-v2026-runs/wj6-20261111/build-v2 \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math'
cmake --build /media/jay/新加卷/nektar-v2026-runs/wj6-20261111/build-v2 --parallel 4
ctest --test-dir /media/jay/新加卷/nektar-v2026-runs/wj6-20261111/checks-v4 -j 1 --output-on-failure -R '^(WJNativeEP|WJNativeBudget|RANSWJSharedLegacyTransaction)$'
```

The receipts also bind every driver invocation for pilot, primary, all-order
qualification, CP24 transaction, T3 cost/admission and publication. Heavy fields,
checkpoints, matrices and native logs remain in Linux storage. Compressed compact
files retain exact original decoded identities; local verification distinguishes
bytes actually transported/read and Git-verified source from remotely verified
heavy identities.

The transported set contains 350 files: 32,787,002 encoded bytes and
183,259,837 decoded bytes. Its index SHA-256 is
`74649d94ba88f1ddb20e5f2ca4642a7a0bfb36da47290b195cf862cda79017df`.
[Independent local verification](wj6-linux/local-verification-v2.json) checks
every encoded/decoded byte, 287 available receipt links, 6183 producer Git
blobs, 1334 passive observer Git blobs, ten frozen-driver inputs and all
declared Git trees. The earlier verification packet remains intact. The
3660 remote-only identities are explicitly recorded remote verification;
they do not claim a second local read of heavy ELF/matrix/checkpoint bytes.
[The outer receipt bridge](wj6-linux/finalization.json) binds the immutable
publication, actual closing seal command and final command receipts without
rewriting the original receipt. Large compact JSON/CSV packets are stored as
deterministic gzip with their original decoded hashes in
[the transport index](wj6-linux/transport.json).

[The contract](wj6-contract.md), [test inventory and results](wj6-linux/regression-verification.json),
[interface tables](wj6-linux/interfaces.csv), and
[horizontal/timestep sensitivity figure](wj6-linux/figures/matched-sensitivity.png)
provide entry points into the complete packets. The other figures retain
divergence, pressure-test projection, drift, interfaces and the distinction
between spatial residual and completed-stage defect.

The supported next step is a separately resource-budgeted T3 follow-up, retaining
this uniform mesh and finite-strip problem. WJ-6 establishes measured Tx
horizontal/timestep sensitivity and drift budgets at 4331 states. It establishes
neither exact finite-space equilibrium, stationarity, long-time preservation,
global stability/realizability nor native 3D/MPI/scalable production capability.
