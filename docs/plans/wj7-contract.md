# WJ-7 corrected lineage and joint uniform-L2 continuation

Reviewed/fetched parent: `84f62f525be35bd073231476c58bc409e741c435`.
Native producer: `30ffa4000353d8055d19abec3c7dc189b7004c4a`;
frozen WJ-6 numerical driver: `4f7520aabee2604efd5b0194a5fa8592d672816a`.
The existing branch is continued without reset, rebase or unrelated merge.
Original/worktree specification SHA-256 values agree:

```
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124 main
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c appendix
```

Retain literal quad/\quad in main eq:integrals-linear and the intended
(d+a*h)/a^2 coefficient. Specifications, WJ-4/WJ-5 negatives, sealed WJ-6
Tx PASS and T3 RESOURCE_LIMIT, all old ledgers and producing paths remain
immutable. WJ-7 is a new allocation, not a relabeling of the denied extension.

## Correction before scientific execution

The development publisher must use distinct report_identity and matrix_identity
variables. Reconstruct all 36 stage bindings from hash-verified sealed execution
reports, compare every retained stage field/callback/audit, and identify exactly
18 erroneous matrix-as-report bindings. Publish a new correction artifact and
bridge to original map/publication/seal; preserve original bytes. Byte hashes
protect integrity, not semantic producer roles. Fixtures exercise real publisher
collection with and without arithmetic descriptors.

## Unchanged physical/discrete problem

Retain complete 3D WJ tensors/cubic-N surrogate, Original stress,
DevStrainRawWorkV1, direct epsilon, Praw=-R:Graw, independent Cd diffusion;
c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigma_k=1,sigma_epsilon=1.3,
Cd=572/6561,qe=81/(2 sqrt(143)),kappa=.381917863546733.
Only the high-Re outer model omits molecular viscosity. AWFs retain nu=1.5e-5,
Linear full-strain stress, uncapped sources. EP uses chi=1,zm=.5,PA-pA0;
budget uses zK=zm=1,PA-Ec,no physical chi, transported nonconstant inner k.
K=k(H) is an endpoint value. No inner epsilon PDE or extra epsilon flux.

u*=.5,H=1,top=2,h=.1,K0=.8466950351321795. Domain [0,1]x[0,2], reference
Dirichlet differential traces on both laterals/top, no-slip ground and modeled
zero ground k flux. Q3 differentials, regional continuous Q1 pressure, one gauge,
pi=p/rho+2k/3, shared U/k, distinct pressure/gradient traces, dependent outer
epsilon trace. IMEX1/PimpleSingleSweep, eta=0, uncondensed, zero prescribed
forces/gravity, reuse/verification Off, StopOnBound On. One U predictor,
pressure/gauge,k,epsilon solve and logical WJNativeV1 bound call per completed
attempt; native hard predicates active, correction inactive, raw/corrected k
identical. No physical/space/BC/stabilization/clipping/backend/3D/MPI change.

Capacity remains execution-only, default3000 and selected-route maximum4331.
New resource-design hashes change execution intent, never physical signatures.
State sizes4331 Nx2 /2493 Nx1; x knots [0,.5,1]/[0,1]. Exact normal vector:

```json
[0.0, 0.000625, 0.00125, 0.001875, 0.0025, 0.003125, 0.00375, 0.004375, 0.005, 0.00625, 0.0075, 0.00875, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.08499999999999999, 0.09, 0.095, 0.1, 0.125, 0.15000000000000002, 0.17500000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.47500000000000003, 0.55, 0.625, 0.7, 0.7749999999999999, 0.85, 0.925, 1.0, 1.0375, 1.075, 1.1124999999999998, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4125, 1.475, 1.5375, 1.6, 1.7000000000000002, 1.8, 1.9, 2.0]
```

## Import, observations and joint barriers

New root `/media/jay/新加卷/nektar-v2026-runs/wj7-20261111`.
Before use, bind native input/ELF/dependency equivalence and immutable six-endpoint
imports: each policy coarse Nx2 CP32 dt3.2e-5, fine Nx2 CP64 dt1.6e-5,
baseline Nx1 CP32 dt3.2e-5. Verify bytes, signatures, signed binary64 primary,
gauge/time/history/fixed-reference and original completed callback. Zero-step
native imports must match with zero initialization Lift/governing solves.
Original CP0 and fixed-reference provenance remain the drift origin. Absolute
native steps, accumulated uncertainty and native time labels are retained.
Reuse the equivalent 38 passing CTests and capacity transaction, with no rerun.

New reads request native constraint maps on FIRST q8/x4 and q10/x6 reads.
Normal quadrature floor1e-10; actual mesh/physical/reporting partitions, face
owners, zero duplicate integration weights, full outer epsilon, one overall
pressure gauge. Keep independent evolved source/diffusion/weak/interface loads,
pressure-test projection/complement, element modes,Gxx/Gzz/extrema,K/cutoff.
No governing diagnostic solve, Lift or reprojection. Spatial residual, actual
stage defect, continuum error and original-CP0 drift remain separate.

Keep only lightweight identities/scalars; stream all historical row maps.
At most current comparison pair/order is full in memory. Replay retained bytes
and reproduce all504 WJ-6 sensitivity leaves and all accuracy/interface/
divergence/quadrature decisions; quantify arithmetic-order differences.

Joint common rows (seconds; coarse/fine/baseline absolute steps):
.002048:64/128/64; .003072:96/192/96; .004096:128/256/128.
Both policies must pass COMPLETE row gates before either advances. Native
invocations <=32 steps; each fine64 scientific interval uses two32-step runs.
Intermediate fine checkpoints are not matched observations. Export passive
arithmetic at every new common terminal step (nine analyses per policy), using
report.steps ABSOLUTELY. Enlarge each trajectory's own envelope and re-evaluate
all applicable prefix gates before adding both contributions. Historical WJ-6
estimates/results remain unchanged.

Unchanged amplitudes aU=.5,ak=K0,ae=.125/kappa,ap=.25; full/regional scopes
use own measures. E2/D2<=1e-4, sampled field maxima<=1e-3, traction/.25<=1e-3,
kflux/(K0*.5)<=1e-4, reference pressure-offset departure/.25<=1e-4,
max|divU|*H/u*<=1e-3. Pair difference plus BOTH estimates<=10% field budget.
Factor100, separate k/epsilon units, dependent trace maps, actual constraint
rows, data/lift/basis/geometry/accumulation and full-pressure gauge remain.
These are finite frozen-system estimates, not nonlinear stability bounds.

Required sensitivity tuple set: coarse keys0,4,32,64,96,128; two orders,
horizontal/time, seven fields, three scopes, two norms =1008 unique leaves.
Every individual retained middle/new observation is required. Missing, duplicate,
wrong-producer, stale, nonfinite or failed evidence cannot pass. Final PASS
requires every target/barrier; hard/material correction stops immediately;
engineering/quadrature/uncertainty failure stops both policies at that barrier.
No adverse continuation beyond it, replacement search or T4 is authorized.

## Resource admission and closure

NEW explicit per-policy ledgers:800 units and10800 numerical-study seconds.
Charge repair/import/fixtures/replay/native/observation/SVD/parent algebra/I/O,
waits, failures and publication; build time separately. Preserve combined WJ-6
plus WJ-7 totals. Unknown reservations remain charged. Continuation core384
attempts per policy; no repeat transaction. At least25% cost headroom and600s
closure retained. Measure retained-data replay FIRST, then write-once initial
admission from verified old/current rates. Version remaining projections before
each subsequent joint row; first successful row is scientific, never disposable.
Do not remove the entire old q10 cost: first-read map work remains.

Finite limits:8GiB aggregate process-tree RSS,8GiB invocation output,32GiB
policy storage including actual copied inputs/all consuming output roots;
native timeouts<=1500s based on supported <=32-step costs (no truncation).
Serial rank, six numerical-library controls1, flags -ffp-contract=off
-fno-fast-math, RANS_BUILD_SOLVER=ON via verified unchanged runtime. Inspect
current capacity; <=72 total workers including MPI helpers/other campaigns.
Terminate only own process groups. One bounded verified external zero-launch
resume is allowed only with proven pre-Popen0-attempt receipt, new output
identity and <=60s allocation wait; charge failed waits. No launched/unknown
physical attempt is retried as zero-launch.

Non-fitting projections stop with RESOURCE_LIMIT and frozen unreached T3.
Terminal statuses PASS,OBSTRUCTED,UNRESOLVED,RESOURCE_LIMIT,NOT_RUN (implementation
FAIL preserved until repaired under a new lineage). Immutable compact correction,
imports, code/runtime/producer/callback/checkpoint links, all-order tuple coverage,
gates/ledgers/log hashes and encoded/decoded transport identities are published.
T3*u*/H=.002048 is finite-horizon strip qualification, not stationarity,
exact discrete equilibrium, a temporal order proof or general stability.

## Pre-import observer implementation clarification

The retained replay exposed a discrete-order parsing issue (orders must be
integers while signed binary64 state/history values stay unchanged). It also
identified that the legacy pair-comparison helper defaulted to horizontal
order4 when composing a q10/x6 boundary rule. That omitted the top's x6
integration nodes, giving boundary measure4 instead of5 (outer epsilon2
instead of3). WJ-7 uses the actual requested horizontal order, preserving one
full-domain pressure gauge and all definitions/budgets. A focused fixture
records the old/new measures. All504 legacy sensitivity leaves are replayed
separately with unchanged legacy arithmetic/rules and retained as such;
the corrected integration is a new versioned qualification of the same bytes.
No WJ-6 receipt or result is replaced. Any changed decision stops admission.

The streamed envelope obtains only the two unique six-scalar native evaluation
arrays in a hash-verified mmap metadata pass, then processes each full row map
once. Each comparison holds only its current pair. Production-parent cost
retains the old measured per-new-observation allowance; complete prefix
requalification uses its newly measured streamed cost. These distinct tasks
are both projected with25% headroom, along with all nine SVD analyses,
first-read row exports, extra launches and600s closure.

The boundary estimator may skip interior point evaluations only because their
boundary mask/weights are identically zero; full-volume gauge rows and the
global basis/coefficient maxima remain. A focused full-versus-filtered fixture
and complete retained replay check this algebraic optimization. Each new
admission uses a separate measured full passive production replay (parsing,
evaluation, independent loads, diagnosis and serialization reserve) plus the
sealed native q8 and q10-with-maps process receipts. Thus q10 row export remains
charged, while the obsolete quadrature-only duplicate and duplicate repair
parent work are not projected as a second campaign. Every cost input, ledger
basis and native receipt is snapshotted immutably before launch. Superseded
unoptimized replay time remains charged to WJ-7.
