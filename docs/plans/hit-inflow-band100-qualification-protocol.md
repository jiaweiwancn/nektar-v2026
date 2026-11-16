# CBC100 dissipation, Fourier, bandwidth and sampling qualification

Status: **AUTHORIZED; IMPLEMENTING AND PREFLIGHTING**, 2026-11-16,
Asia/Shanghai. Branch: `codex/hit-inflow-generator`.

## Material Passport

- Workflow: academic-research-suite / experiment-agent, execution and validation.
- Authority: the user's four-point follow-up to the completed CBC100 review.
- Starting source: `b42d04e`; native build: `1946f8e`.
- Reference: the immutable `f4b70ba/study-01` base/dt/x records and
  [their completed review](hit-inflow-band100-review.md).
- Server: `jay@100.105.247.110`; hard ceiling **64 total workers**.
- Current disk mount: `/media/jay/新加卷1`, same external experiment disk.
  Archived `/media/jay/新加卷` paths are explicitly mapped for reads; original
  receipts, inputs, histories and fields remain unchanged.
- Source work: the attached HIT worktree. The main Documents checkout belongs
  to another branch/task. The old Linux checkout has an unrelated NTFS symlink
  change; preserve it and synchronize a separate clean HIT checkout.
- Scope: bounded qualification controls; full-distance production remains held.

## Common short-domain definition

Retain the economical pilot: 1 x 0.6 x 0.6 m, hx=0.1 m, hy=0.05 m, degree10,
dt=25 microseconds, U=12.7 m/s, nu=1.5e-5 m2/s, IMEX2, both dealiasing switches,
FFTW, pressure XxtMultiLevelStaticCond, HIT velocity/pressure-H inlet, periodic
y/z and native HOutflow. Keep the same stations0/0.15/0.45/0.75m, 32x32 half-open
probes, seed20261028, periods(10,0.6,0.6)m and fixed synthesis counts
**(2548,154,154)**. Every primary record is0.2--0.6s inclusive, sampled at0.00025s,
1601 frames. Mean-only initialization and Preserve input normalization remain.

Changing a numerical knob must not change the inlet realization, physical
probes or clock. Compare bandwidth at a common numerical recipe. Do not tune
coefficients to downstream experimental data. The pilot has no experimental
reference at its interior stations.

## Stage A: DG-only half-strength sensitivity

Run coarse hx=0.1 with DG coefficient0.5, keeping homogeneous SVV exactly at
the base0.005953125m2/s and ratio0.75. The archived controller's `svv` case
halves both coefficients, so introduce a distinct `qualification` recipe;
preserve the old recipe's behavior for reproduction.

Obtain a fresh generator/direct-reference benchmark, native boundary projection,
and160-step stability/timing smoke for this exact case and allocation. Only
then record the full0.2--0.6s window. The primary comparisons are:

1. New coarse DG0.5 versus old coarse DG1: coefficient sensitivity on the same
   mesh, with homogeneous SVV held fixed.
2. New coarse DG0.5 versus old finer-x DG1: both have scalar DG prefactor
   C*U*max(hx,hy)/p=0.0635m2/s. Their physical modal filters and approximation
   spaces still differ. This is not an identical-operator proof.
3. Native final-field diagnostics at the common0.6s time for old base, old finer-x
   and the new control: molecular strain work, DG work, homogeneous work,
   slabwise energy and fixed-reference divergence. Hash every partition before
   and after reading, write diagnostics in separate directories and retain
   global/slab closure. These snapshots are descriptive; they do not estimate
   temporal-central dissipation or time-averaged confidence intervals.

Retain signed changes, all components/stations, the full variance/K screen2%,
the four10--100rad/m band screens5%, and F11(100). Use the PSD partition above
and below100 to locate changes in represented temporal energy. A successful
stability screen does not select0.5 as a production coefficient.

If the matched-prefactor mesh comparison still fails2%, inspect which
components/directions remain sensitive before adding a spatial level. One
candidate is hx0.025, hy0.05, degree10, DG1: it retains the0.0635 prefactor
because the maximum edge stays0.05. It still changes the physical DG filter.
This candidate requires its own resource preflight and will not be launched
merely to fill a convergence table. Consider a second coefficient at fixed fine
mesh if needed to distinguish mesh/coefficient interaction. Report failure
honestly if affordable levels do not establish full-energy/decay convergence.

## Stage B:192 slots, then matched radial bandwidth

First run the coarse base at192 Fourier slots and600/480rad/m. Retain the
registered homogeneous scaling0.1*U*Lz/Nz, hence0.00396875m2/s at192, ratio0.75.
This is a resolution-plus-native-SVV control: both its amplitude and physical
filter onset change from128. If this remains material to energy/decay, add a
separate homogeneous-SVV sensitivity rather than claiming pure Fourier error.

After its successful analysis, compare800/640 against600/480 **at common192**,
with identical synthesis counts, seed, periods, DG/homogeneous coefficients,
mesh, dt, probes and recording clock.192 slots have Nyquist1005.31rad/m on0.6m;
128 slots have670.21rad/m and cannot admit800 under the strict-Nyquist rule.
The synthesis lattice supports both bands without count changes.

Before broadening, record the600 benchmark's mode bytes and cost, estimate the
additional800 storage conservatively, then measure a fresh800 benchmark and
MPI projection/smoke. Record native retained-mode counts and actual memory;
continuum k-cubed scaling alone is not a memory guarantee. Do not start the
full800 recording if the gate rejects it. Bandwidth changes the physical input
and is model sensitivity, not a numerical-refinement pass.

The upstream uncut-model F11(100) retention remains89.83% for600 and95.05% for800.
These percentages are not experimental error estimates. Neither band retains
the complete molecular dissipation; no full-band DNS claim follows.

## Stage C: recording, startup, plane sampling and realizations

After A/B determine an affordable, numerically defensible recipe, qualify its
sampling. Keep the existing primary comparisons intact; added sampling records
must contain nested32/64 probes and aligned physical windows so that their
differences are paired. Assess at least64x64 versus its nested32x32 grid,
and128 or a shifted dense grid if64 remains sensitive. Streaming/storage and
analysis memory need a new preflight; denser probes are not independent runs.

Inspect full/final-half/last75% moments and spectra; adjacent windows after
multiple1m convection times; temporal and block-energy correlations; multiple
block lengths (4/8/16 plus physically meaningful duration). Distinguish
start-up drift from finite-record variation. Longer duration requires a new
synthesis period/count design: the old10m nonrepeat rule limits the complete
1m run to less than about0.7086s, so it cannot support an arbitrarily extended
record. If the period/counts change, label comparison with the old realization
unpaired and keep every new numerical/bandwidth comparison paired within the
new lattice. Do not concatenate incompatible restarts or reset inflow time.

Use at least one additional independently seeded run before estimating
between-realization variability; retain its own complete window. Two seeds
are a sensitivity check, not a strong precision guarantee. Increase independent
seed/block evidence if required by the observed variation and resources.
Nonoverlap,1024 probes and overlapping Welch segments do not establish
independence. Conditional Student-t intervals remain conditional until
stationarity, autocorrelation and effective block length are assessed.

Use2% as the moment sampling/refinement screen and5% for the registered bands,
but report the actual evidence and limitations, rather than a binary precision
certificate. All11 statistical-fallacy categories are checked on interpretation.

## Stage D: full-distance definition, conditional on prior gates

Do not reuse the mean-initial0.6s schedule or10m synthesis box for9m. Experimental
stations after theta42 are3.612896m and8.322564m, not the1m pilot stations.
At12.7m/s,9m alone takes0.70866s to convect; the last experimental station takes
0.65532s. Recording must follow measured initialization/flushing adequacy and
allow enough developed stationary data, not merely begin after one transit.

The existing nonrepeat-footprint rule is
`L_synthesis > length + U*(duration + 2*dt)`.
For9m,1.6s,25us it requiresL>29.320635m.32m is a candidate, not an approved
resource recipe. At800rad/m a strict-Nyquist even streamwise count must exceed
8148.73 (for example8160); transverse counts154 still support800. Freeze common
counts across any paired600/800 comparison. A longer qualified duration may
require a box larger than32m. Record finite-lattice/input expectations anew.

Measure mode construction, memory duplication by MPI rank, solver initialization,
developed step cost, dense recording/storage and analysis cost before defining
the full run. Consider fewer ranks if it lowers memory within the same64-worker
cap; benchmark the tradeoff. Require an adequate downstream buffer/outlet test
and a transverse-span sensitivity at the target distances. Do not interpret an
unchanged0.6m span or0.25m pilot buffer as already qualified for9m.

Close each upstream gate with evidence before defining/preflighting production.
Any failed numerical, model, sampling, startup, span/outlet or cost gate keeps
full-distance production held and triggers a bounded remedy, not tuning.

## Execution and preservation

Use one MPI command at a time, initially32ranks/npz16 with an explicit CPU list
and all numerical-library thread counts1. Include CPU0 in the startup pool;
the32 numerical ranks and their controller share the declared pool. Audit
all main/helper threads, readiness, outside project workers and immutable
executable/library hashes. Never exceed64 combined workers. Recheck outside
load before every command. A resource violation stops only the owned command.

Retain the existing12h forecast /180GiB sampled-RSS preflight gates and16h hard
recording timeout. Qualification also enforces the180GiB sampled-RSS bound
while a command runs. No failed command is retried unchanged. Each preflight
is bound to its exact case, source, executables and allocation. Keep failed or
rejected attempts, logs and receipts. A durable sequential queue stops on
failure/rejection and ends with a review-required status; it cannot start9m.

The original frozen protocol is unchanged. Reference bridging requires an
explicit original full revision, matching scientific helper/data hashes,
identical native executable/library bytes, all command/log/prepared-input
hashes and matching completed analysis receipts. Only the controller/revision/
checkout changes. Run artifacts and fields stay on the external disk.
