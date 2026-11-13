# WJ-H14: fixed-height mapped backend and target admission

Frozen 2026-11-13 before H14 numerical/build work. Actual/fetched starting
revision is `fb00b9c1708a387e10097f5e3b217cc1f8b62558`, equal to the reviewed
remote. The existing WJ worktree is clean; the unrelated main checkout is
untouched. No prior H14 source, Linux root/job or resource journal was found.
There are therefore no raised-height artifacts to continue or earlier H14
charges to discard. If any are subsequently discovered, retain and charge them.

The user withdrew raised matching heights. **Only the existing H13 H=6 mm
primary is active**, with physical flat matching plane9.5 mm. The stored5/5.5 mm
records remain historical metadata. No other height, calibration/roughness
inverse, fit, height-family generator, normal partition regeneration,
observation reclassification or raised-height clock calculation is authorized.

## Frozen scientific target and execution differences

Reuse H13 `publication-v1/paired-case/flat/case.xml` and its manifest by exact
hash. Keep all physical Functions/parameters, topology/coordinates/IDs,16 inner
+34 outer cells,14 columns,700 chart elements, Q3/Q1, physical map, observations
and complete boundary laws. State24688; U/W/k6493 each; epsilon4429;
regional pressures255/525; free velocity12558, pressure780, flow13338,
k6342, epsilon4327 including43 dependent identity rows, independent epsilon4284.
Keep representation/lift GLL8, canonical mass Gauss8, physical mass Gauss16,
governing q12/normal16. Native numbering/signs must independently confirm these.

Keep full3D Original WJHighRe, DevStrainRawWorkV1, raw full-stress work and Cd
diffusion, uncapped WJCalibratedBudgetAWF, positive AWF nu, nonconstant inner k,
current represented-K metric cutoff/dependent epsilon, shared U/k and regional
pressure/traction coupling. Keep FlatInletTraceV1 and reference freeze order.
The original pressure jump is a seed; no new pressure constraint/gauge is added.
Molecular viscosity is omitted only by the existing high-Re outer operator.

Full-height neutral physical Functions, fixed top reference shear/epsilon flux,
complete outlet tractions and backflow stop remain unchanged. U at800 mm stays
about1.09003 U_infinity. Physical7 mm remains inner AWF, lowest Fig3 stresses
and the10 mm crest-clearance target remain outer. No hill governing route is
authorized; its H13 partner is geometry only. Preserve original literature hashes,
`quad`/`\quad` erratum and `(d+a*h)/a^2` coefficient, H13's closed ledger,
H12 audit PASS and H11 temporal UNRESOLVED.

Allowed execution-only changes: new backend selector/capability/resource
identity, supported compact diagnostics, output/checkpoint paths, IO controls,
NumSteps and a subsequently admitted pilot TimeStep. Record an explicit XML
diff and new source/session identities. No old strip checkpoint substitutes for
the target. Fresh target CP0 follows complete small-case qualification only.

## Backend and qualification gates

Assess graph-permuted band storage from actual nonzero graphs/native row maps,
using existing Lapack Dpbtrf/Dpbtrs for consistent mass and Dgbtrf/Dgbtrs for
nonsymmetric blocks. Preserve explicit permutation/sign roles. Preflight factor
storage before allocation. Retain default dense code and its3000/4331 guards.
No guard-only increase or altered map Identity. Select/restrict the new route
before constructor/reference allocations; NumSteps0 is not admission.

Remove indispensable dense paths through canonical/physical mass and projection,
frozen Picard assembly, constrained/mass-dual norms, PIMPLE, scalar stages,
bound metric and required diagnostics/restart/strict readback. Reuse shared mapped
element/face/ray contributions. Optional legacy dense dumps may be rejected
explicitly before allocation. A partial port never opens target preflight.

Keep A predictor and B=(alpha*M_free)^-1 distinct. Full nonsymmetric Schur is
S=D*B*G-C with complete roughness/exterior/interface pressure actions. A bounded
pressure-only dense Schur may be streamed columnwise, discarding BG columns;
independently reapply D*B*G-C to validate its solution. No lumping, -D^T shortcut,
SPD Schur assumption or full momentum inverse in B. Keep k then epsilon,
trace identity rows AND mass/transport columns, one solve per field per sweep,
all state-dependent data, cache invalidation, bounds, histories and rollback.

Bound whitening/recovery must use the actual consistent mass factor and its
permutation: L^-1 P a and P^T L^-T y. Retain two-pass QR/rank2e-11 and retained
pool certificates, including activated and near-dependent witnesses. No Gram
normal-equation or preconditioner metric substitution.

Freeze action parity <=1e-11*max(1,||reference||inf), unambiguous solve/sweep
parity <=1e-9*max(1,||reference||inf). Preserve governing true residual checks:
predictor/pressure/scalar<=1e-10*max(1,||rhs||inf), splitting<=1e-9*max(1,||Gdp||inf),
corrected continuity<1e-10 and essential defect<1e-13. Do not relax after a failure.
Compare existing admitted small mapped/dense controls, signed two-column fixture,
pressure/boundary actions, trace columns, activated bounds and restart/history/
initializer/rollback. Full-resolution Nx2/4048 remains outside dense admission.

## Conditional target and pilot

The new external source-bound preflight stays closed until every indispensable
backend path and qualification gate passes. Then permit one fresh primary target
CP0 and at most four completed one-sweep pilot steps. Freeze pilot dt after the
small-case/advective/diffusive assessment, before target evolution; do not assume
the old metadata .01 is safe. No adaptive timestep/height/mesh retry.

CP0 requires all regional continuity moments<=2e-15, physical divergence
<=1e-12*max(1,||G||F), exact essential/dependent/reference invariants and unchanged
native physical predicates. Independent finite-inventory principal quality remains
.002068763923454, separate from hard signs. Any hard/quality/backflow/bound-work/
resource/unsupported-mode exit stops admission or pilot; preserve its first reason.
No nonzero floor/cap/correction/restoration is introduced.

Record q0, actual post-lift/initializer qh0 and frozen internal reference/defect
separately. Drift uses qh0; analytical errors use q0; retain physical regional
norms, owner-preserving maxima and height-relative outer epsilon. Measure startup,
assembly/application/factor/solve/predicates/observation/I/O and actual storage/RSS.

H13 clocks are reused, not recalculated for another height: bulk26.13350407707765,
slowest outer57.17628191545675, inner diffusion195.80590488759879, top K/epsilon
505.02466521278404; proposed duration505.02466521278404 with terminal window
[447.8483832973273,505.02466521278404]. No duration/hill run is authorized in H14.
Forecast later preservation using actual target cost/dt if the pilot is admitted.

## Resource admission and publication

New root `/home/jay/build/wj-h14-mapped-backend-20261113` on healthy ext4;
documented SSH alias failed and Tailscale fallback works. Inspect resources and
bounded write/fsync/readback before launches. New persistent7200 s ledger, fully
charged60 s setup and reserved/charged300 s closure. Phase ceilings: build2400 s,
small qualification1800 s, target CP0600 s, pilot600 s; remaining operations and
repairs share the same total. Failed/unknown reservations remain charged; no reset
or automatic native retry. Update remaining-cost forecast with25% headroom before
target admission. Incomplete implementation/qualification keeps target NOT_RUN.

One serial native scientific worker and all six library thread controls1.
Compile at most8 workers initially (never more than16), host total<=72 including
other campaigns. Build process-tree RSS<=16 GiB; native/test RSS<=8 GiB; source/
build disk budget8 GiB; new science output<=2 GiB, preflight actual free capacity.
Use existing atomic-journal/process-group helpers, not another accounting system.
Kill only this package's process groups on a resource ceiling.

Commit coherent implementation/tests, compact source/build/runtime/allocation/
parity/bound/history identities, target evidence if admitted, exact failures and
closed accounting. Normal branch push only. Separate implementation, hard/
quality qualification, pilot, full-domain duration and experimental statuses.
No physical H is selected; the7 mm/all-outside conflict, roughness-sublayer
validity and idealized-inflow limitations persist. EP, ambient/finite-depth,
old measured-rms/20 mm, anisotropic stress, old ratio audits, capacity shortcuts,
additional heights and hill production remain outside this package.
