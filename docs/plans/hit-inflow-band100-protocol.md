# CBC inflow pilot for a longitudinal 100 rad/m comparison

Status: **IMPLEMENTATION AND PREFLIGHT**, 2026-10-31.

This follow-up implements the user's authorized 100 rad/m proposal on
`codex/hit-inflow-generator`. M1--M8 and the unfavorable M7 comparison remain
unchanged. This protocol defines the short-pilot gates before a new full-distance
run. It does not claim full-band DNS, independent probes, or experimental fidelity.

## Material Passport

| Field | Value |
| --- | --- |
| Workflow | academic-research-suite / experiment-agent, execution then validation |
| Source starting point | `6307bcec190adf0af0416096ca2b261f34fd90d3` |
| Authorized target | Longitudinal spectrum through 100 rad/m; CBC first |
| Environment | Linux `jay-Super-Server`, Tailscale `100.105.247.110` |
| Resource ceiling | 64 total workers; one MPI case at a time; OMP/BLAS/MKL=1 |
| Heavy-data root | `/media/jay/新加卷/nektar-v2026-runs/hit-inflow/band100/` |

## Frozen pilot definition

- Mean velocity 12.7 m/s; molecular viscosity 1.5e-5 m2/s, as M7.
- Nonlinear incompressible velocity correction, IMEX2; both dealiasing switches;
  FFTW; iterative velocity solves and pressure Xxt/MultiLevelStaticCond.
- A 1 by 0.6 by 0.6 m domain, periodic y/z, HIT velocity Dirichlet and pressure H
  at the inlet, native HOutflow at the outlet. Stations 0,0.15,0.45,0.75 m leave
  a 0.25 m outlet buffer. No station in this pilot is theta=98 or 171.
- Base h=0.05 m, polynomial degree10, 128 full real Fourier slots.
- CBC station42 shell spectrum, LogLinear, K4/Gaussian extensions, Preserve.
  Base radial cutoff/taper600/480 rad/m. This is an input feasibility case.
- Synthetic periods (10,0.6,0.6) m and **fixed counts (2548,154,154)** for all
  controls, including800/640. Seed20261028. Changing counts changes mode IDs,
  so automatic counts are not used for a matched bandwidth comparison.
- Base timestep12.5 microseconds; recording every0.00025 s. Nominal duration
  0.6 s and primary recording0.2--0.6 s, 1601 frames. The repeat period is
  10/12.7=0.7874 s; complete space-time footprint is approximately8.6203 m.
  Nominal flushing alone does not qualify startup or spatial development.
- Mean-only initial bulk field avoids evaluating a broad Fourier realization
  throughout the volume. Initial pressure/velocity are constructed natively.
- Native DG SVV coefficient1 in x-y. Fourier exponential coefficient follows
  the M7 dimensionless scaling0.1*Ubar*(Lz/Nz), ratio0.75; at128 slots this is
  0.005953125 m2/s. Resolution controls therefore also change the SVV model.
- Histories initially use a common32 by32 half-open y-z probe grid. They are
  sampled-probe statistics, not alias-free full-plane averages or independent
  replicates. An increased-probe recording control remains required for precision.
- Omit cumulative full-field stress output in timing/temporal-spectrum pilots.
  Keep final fields, logs, receipts and the finite native histories.

## Ordered work and gates

1. **Implementation:** replace evaluator auxiliary per-mode indexing with bounded
   lattice lookups and a bounded separable spatial path where point-by-bin tables
   do not fit. Keep deterministic coefficients, phases, original dense arithmetic
   where applicable, strict time/range checks and direct-reference validation.
   Verify meaningful native unit/regression tests on Linux after synchronization.
2. **Generator/input preflight:** native benchmark with bounded reference points,
   coefficient and direct/grouped agreement, mode/cache bytes and timing. Compute
   continuum projected spectra and retained K/epsilon without renormalization.
   Cache/direct agreement is not a constrained boundary-projection check.
3. **Timing smoke:** fresh mean-initial case of160 steps, with individual MPI ranks
   pinned to distinct assigned CPUs for numerical work. Native initialization
   was observed to move a main thread briefly to CPU0 and restore it. Audit all
   main/helper thread masks against the explicit CPU pool throughout startup;
   once a rank reports `HIT_NATIVE_READY`, require its assigned singleton CPU on
   every thread. All ranks must report readiness before solver steps. Retain
   startup-affinity witnesses instead of misclassifying an in-pool initialization
   probe as a permanent binding failure. CPU-pool escapes are always failures.
   Inlet-owning ranks additionally report `HIT_EVALUATOR`; ranks without inlet
   faces still report global readiness. Do not infer all-rank readiness from the
   subset of evaluator objects.
   Record initialization/steady step cost,
   sampled RSS, warnings, final native finiteness, evaluator method and hashes.
   A hard timeout bounds every child job; failed jobs are retained and not retried
   unchanged. No broad-case silent direct evaluation is permitted.
4. **Budget gate:** estimate the0.6 s baseline from steady solver steps, excluding
   one-off construction. Do not start a case whose estimated wall time exceeds
   12h or whose projected sampled memory exceeds180 GiB. Timing estimates are
   planning observations, not runtime guarantees. A case has a16h hard timeout.
   If the gate fails, record the measured bottleneck and stop that costly recipe.
5. **Baseline and matched controls:** run cases sequentially, each in a fresh
   directory. Start with the600 base, then halve dt, increase Fourier slots to192,
   and compare800/640 with600/480 **at common192 slots**. An h=0.025 control and
   half-SVV-strength control each require a timing preflight before their full
   recording. Retain a common physical probe grid, synthesis lattice, seed and
   recording clock. Do not compare different bandwidths as pure numerical error.
6. **Numerical/sampling review:** retain2% component-variance/K refinement screens
   and report5% sensitivity screens for integrated longitudinal bands10--20,
   20--40,40--70,70--100 rad/m at every station. Report F11(100) separately.
   Use periodic Hann Welch512 samples/256 overlap, plus full-record boxcar
   Parseval. Compare full/final-half/last75% windows and4/8/16 conditional block
   intervals; correlated probes/segments do not establish independent degrees
   of freedom. A failed sampling screen precludes a precision-qualified pass.
7. **Physics gate:** show input continuum, finite seeded realization, projected
   boundary signal and downstream CFD separately. Table2a at station42 is only
   an upstream descriptive comparison; no fitting to theta98/171.600/480 retains
   89.8308% of the uncut model's F11(100),800/64095.0493%;600 is not expected to
   pass a strict5% endpoint bandwidth test. Input interpolation and missing
   dissipation remain explicit model errors even when the target band agrees.
   Before a costly CFD recording, a boundary-only temporal diagnostic may compare
   the unprojected and constrained-projected values at identical native points
   over0.2--0.6s (1601samples). Bound CSV output to16points per local inlet and
   reuse the fixed10--100rad/m bands. The broad600-band maximum point error is
   not itself a100-band error. If the proposedh=0.05 trace fails the5% band-power
   projection screen, inspect an h=0.025 boundary-only control before deciding
   whether any full CFD case is justified. Report all components and endpoint
   changes separately; no exact-plane or experimental-fidelity claim follows.
8. **Full-distance decision:** no new9m production case is launched until native
   numerical, sampling, startup, bandwidth, projection and resource evidence
   justify it. A developed upstream reference and outlet/span sensitivity are
   still needed. If a gate fails, deliver the diagnostic outcome and next bounded
   remedy; do not label the proposal or M7 physically validated.

## Evidence and operation

Source/inputs, exact commands, executable/library identities, CPU/thread
allocation, running/failed/completed receipts, native logs and checksums are saved
with each fresh run. Controllers are bounded to the same64-CPU pool; builds and
MPI jobs are sequenced. No unrelated jobs or checkouts are modified. The main
Documents checkout currently contains another task's RANS changes, so the
attached HIT worktree remains the source workspace for this branch.

The user previously requested a status script instead of repeated long-run
audits. A progress command will accompany any long-running pilot; completion of
source delivery or job launch is distinguished from completion of numerical
qualification.
