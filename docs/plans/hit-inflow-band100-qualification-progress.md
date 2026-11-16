# CBC100 qualification: fresh preflight passed, sequential controls scheduled

Status on2026-11-16, Asia/Shanghai: **DG-only preflight passed; control queue
active, temporarily held for existing server work.** No completed new0.6s
control analysis or convergence decision is available yet. Full-distance
production remains held.

## Material Passport

- Workflow: academic-research-suite / experiment-agent, execution.
- Authorized scope: the user's DG-SVV, Fourier/bandwidth, sampling and conditional
  full-distance follow-up, specified in the
  [qualification protocol](hit-inflow-band100-qualification-protocol.md).
- Executing source: `a1eab12a381cda77da514b6fcb22ea06700f82ca`.
- Native build: `1946f8ea1d6e5a8844ce308def5004994a40557b`; binary/library hashes
  match the completed original study. Compiled source files are unchanged.
- Branch: `codex/hit-inflow-generator`.
- Clean Linux source: `/media/jay/新加卷1/nektar-v2026-hit-qualification-20261116`.
- Run family: `/media/jay/新加卷1/nektar-v2026-runs/hit-inflow/band100/qualification-20261116`.
- Queue: child directory`controls-01`, launch PID354984.
- Allocation: **32MPI ranks, npz16**, numerical CPUs32--63 plus CPU0 for
  controller/startup containment; all numerical-library thread counts1.
  This is a33-CPU pool, including startup/controller, within the shared64-worker
  limit. The queue runs one MPI command at a time.

The main Documents checkout and original Linux checkout belong to other work or
have pre-existing changes and were left as found. The existing data disk is now
mounted at`/media/jay/新加卷1`. Original archives/receipts remain unchanged; an
explicit mount map resolves their old`/media/jay/新加卷` paths only for reading.

## Delivered and checked

The new`qualification` recipe halves only the x-y DG coefficient for case`svv`.
Homogeneous SVV remains at0.005953125m2/s, ratio0.75. The existing
`original`/`economical` recipes retain their archived behavior, including their
old combined-coefficient sensitivity. Source/native bindings, exact primary
sampling identities, completed reference receipts and case-specific preflights
are required before comparison.

Linux passed **35 controller/analysis checks** at the executing revision. The
reference bridge revalidated all15 original command receipts, their logs and
prepared inputs, all three completed analyses and all three full native-history
hashes. Scientific preparation/analysis/data files and native executables and
libraries match. This is a controller-only bridge, not a waiver for changed CFD
code or input. A separate owned-process fixture verified the resource scheduler's
PID/start-time protection, including changed and exited process identities.

No native rebuild was required. Tests used one CPU63 and one numerical-library
thread; native conclusions below come from Linux.

## Fresh DG-only preflight

| Check | Outcome |
| --- | --- |
| Native grouped/direct generator benchmark | Passed;6565420retained modes,577756960mode bytes |
| Native boundary direct-reference check | Passed; projection errors finite |
| MPI readiness/affinity | All32ranks ready; recorded startup/strict main/helper checks |
|160-step stability smoke | Completed, finite native history, no convergence/nonfinite failure |
| Maximum second-half chunk rate |1.222179s/step |
| Forecast for24000steps /0.6s | **8.177h** |
| Smoke sampled peak tree RSS | **82.457GiB** |
| Projection sampled peak tree RSS |84.698GiB |
| Registered12h/180GiB preflight gate | **Passed** |
| Maximum sampled aggregate worker bound |62, including other observed project processes |

The forecast is a short mean-initial measurement, not a guarantee of developed
runtime, stability or peak memory. The queue also enforces a180GiB sampled-RSS
bound and the64-worker audit during qualification commands. Every recording
has a16h hard timeout. A failed/rejected case stops the queue without an
unchanged retry.

The new primary case retains the original realization, synthesis counts/periods,
32x32 probes and0.2--0.6s recording window. The log explicitly reports
`0.5*Uh/p` for DG and the unchanged homogeneous coefficient. Its full record
has been prepared and allocated; no new control has yet completed.

## Shared-server scheduling hold

The existing RANS experiment uses an older guard that counts operating-system
threads. It observed about100external threads during our32-rank MPI preflight
and suspended its own family. Our MPI service threads share their assigned
singleton rank CPUs; the HIT audit's largest aggregate worker upper bound was62.
The two counters measure different things. The older guard's high-water record
is preserved and must not be interpreted as evidence that100CPU workers were
used by HIT; its package may need to account for this scheduling incident in its
own resource review.

To avoid starving that existing experiment, the HIT study controller and its
verified34-process MPI family were held. Only HIT processes were signaled. The
existing guard then allowed its family to run again. Fields, logs, solver state
and the completed preflight are preserved. The hold does not reset solver time,
inflow phase, startup or the physical recording window.

The external, hashed resource scheduler is PID359127. It shares CPU0 and resumes
the saved HIT process identities after existing RANS guards exit,30seconds of
clearance and a passing shared-worker audit. It verifies source/native bindings
before resuming. Its resource wait is bounded to4hours from17:09:18China Standard
Time; a hard wait timeout stops only the saved HIT family and retains partial
evidence. The full recording's16h deadline remains in force, including elapsed
hold time. Active solver cost and resource-hold time must be reported separately
when interpreting the final wall time.

Live authority: `resource-hold-01.json`, `controls-01/queue.json`, the current
stage's`progress.json`/receipts, and native logs. The saved evidence is a snapshot,
not a claim about the later live state.

## Scheduled controls and later gates

The durable queue will run, sequentially:

1. Coarse DG0.5 at128slots/600--480, compared with the immutable coarse DG1 base.
   Its comparison with the old finer-x DG1 case will assess remaining mesh/filter
   sensitivity at equal scalar0.0635m2/s prefactor. It is a dissipation
   sensitivity, not a fitted production coefficient.
2. Coarse DG1 at192slots/600--480. The registered homogeneous0.1*U*dz scaling and
   modal-filter change remain explicit parts of this Fourier/SVV sensitivity.
3.800--640 against600--480 at common192slots, unchanged synthesis counts
   `(2548,154,154)`, seed/periods and all other numerical/sampling settings.
   Before broadening, a conservative memory plan and a fresh800 native
   benchmark/projection/smoke must pass. The128-slot grid is not used for800.

Each case has its own complete preflight plus a fresh timing smoke before the
full recording. After successful records, the queue extracts native final-field
molecular/DG/homogeneous work and slab diagnostics for the old base/finer-x and
all three new controls, with before/after input hashes. These are single
instantaneous snapshots, not time-averaged dissipation estimates.

The queue ends at`COMPLETED_REVIEW_REQUIRED`. It cannot launch extra mesh levels,
denser/longer or independent-seed records, or9m production. Their recipes depend
on the completed controls, so they remain **pending rather than passed**. The
next review will assess full variances/K and decay, the10--100rad/m bands and
endpoint, dissipation/work, and remaining SVV/mesh confounding before selecting
the bounded sampling work. Sampling will preserve paired windows/probe subsets
and assess blocks/seeds;1024correlated probes are not1024independent runs.

Only after those gates pass will a full-distance recipe be defined and
preflighted, including a longer synthesis box/counts, initialization/flushing,
outlet/span and mode/memory/timing checks. A32m box is only a candidate. The old
mean-initial0.6s/10m recipe is not reused for9m.

## Read-only progress command

From the Mac:

```bash
bash /Users/jay/.codex/worktrees/hit-m6-validation/nektar-v2026/solvers/IncNavierStokesSolver/Tests/HITInflow/check_band100_controls.sh \
  /media/jay/新加卷1/nektar-v2026-runs/hit-inflow/band100/qualification-20261116/controls-01
```

This prints the resource hold/resume state, current case/command, native step
progress and final failure/review status without launching work. Initialization
and the first1000-step interval may have no reported step yet. Native timestamps
and resource timestamps are shown separately from queue-stage changes.

## Evidence

- [Compact hashed Linux launch/preflight snapshot](hit-inflow-band100-qualification-launch-evidence.json).
- [Qualification protocol and conditional gates](hit-inflow-band100-qualification-protocol.md).
- [Queue controller](../../solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/queue_band100_controls.py).
- [Bounded resource scheduler](../../solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/resume_band100_resources.py).
- [Original completed-study review](hit-inflow-band100-review.md).

The compact evidence matches Linux SHA256
`1de2a5fa170a67d1404ab61d968fe350714d8d1a19211ced63f9db99d5f2e60a`.
Large native histories, checkpoints/fields and logs remain on the Linux disk.
