# CT16: separate bounded production-path engineering

Status: **IN PROGRESS; resource admission remains UNRESOLVED**. The user
explicitly authorized CT16 after CT15 closed UNRESOLVED at
`5168ed5fba02a4c92a3847827be3a2dbb6ca071b`. The objective remains
`PASS_RESOURCE_ADMISSION_T_0.407` only when all G1–G6 pass. CT15's evidence and
allocation are closed historical records; none has been reset or overwritten.
Physical T=.407 remains **UNATTEMPTED**. Full-transit production and experimental
replication remain **UNQUALIFIED**.

The [prospective contract](cao-tamura-accepted-16-evidence/contract.json) fixes
the start at **2026-11-16 13:30:25 UTC** (21:30:25 Asia/Shanghai), with 24,000
process CPU seconds, 28,800 elapsed wall seconds, 16 GiB new live data and
16 GiB per-process address space. The clock did not reset at local midnight.
One numerical family/rank, six numerical-library thread settings of one, two
compiler jobs, at most 16 project workers and 72 aggregate Linux workers apply.
The external supervisor counts closed plus live work while the worker guard
suspends an owned family. Failed work is charged. Heavy files and private
runtimes are under
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-16-20261116`. Shared installations,
unrelated work and historical archives remain intact.

| Gate | Decision and evidence |
|---|---|
| G1 | UNRESOLVED. Passing small operators/dependency/history controls and selected first-two-step comparisons are partial evidence. Affected selected native/q48f16/q64f16/q64f32 and alpha4000/8000 controls are incomplete. |
| G2 | UNATTEMPTED for every new CT16 component. CT15 V4's fresh .002 pair remains a separate PASS. |
| G3 | UNRESOLVED. Timing-off 0/1/2 native models and independent-audit costs are measured; the complete workflow and frozen coarse6-to8 restart/retention heldout are unattempted. |
| G4 | UNATTEMPTED. No selected native two-cycle retirement/recovery policy was executed or qualified. Small native restart/rollback controls do not substitute for it. |
| G5 | UNADMITTED. Every measured native-only conditional entry extrapolation already exceeds 24,000 CPU seconds. Complete costs and qualified retention peak remain unknown. |
| G6 | Source/contracts and progress are pushed; final closed receipts, reproduction and remote containment are pending. |

The fixed physical manifest SHA256 is
`004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`:
Original/Natural A95, 62,533 slots/961 elements, P4/P2, HA/128, unchanged floor,
displacement, calibration, hs, normalization and fitted points. Actual current
endpoint K still drives coefficients, cutoff, production and sink. Transported
inner k, dependent inner epsilon without transport, finite inner molecular
viscosity, high-Re outer equations, signed full Cartesian stress/production,
covariance, radicand and transmission inventories remain. Full-height inlet,
Richards–Norris roof, absolute nonlinear outlet traction and absolute pressure
remain. IMEX1 keeps one flow→epsilon→k sweep, dt=.00025/.000125, exact stage
dependencies and the measured post-k epsilon lag. VectorCOLAMDV1 and
`D[(alpha Mfree)^-1 G]-C` remain; diagonal mass only preconditions. The Off
oracle, ExactStepV1 open-route guard, residual bound
`2e-12*maxabs(rhs)+1e-25`, three-refinement cap and capacity limits remain.

Timing-off native runs use the actual selected dt=.00025 case and otherwise
matched scientific work. Every completed candidate's exported accepted vector
through step two is binary64 identical to CT14 (maximum coefficient difference
zero). This does not establish a fresh .002 physical pair or its uncertainty
margin. CT15's separate passing pair has pressure spread
`2.2874413727322336e-7 <= 3.144090599826331e-7` and fine traction RMS
`8.446814744e-4 <= 1e-3`; the earlier .001/.0005 pair remains a pressure failure.

The nonnegative short-run envelope is
`startup=max(0,C0)` and
`step=max(0,C1-C0,(C2-C0)/2,C2-C1)`. Diagnostic timers do not supply acceptance
speedups. Native commands include their native checks/startup/output; they do
not include independent Python audits or the complete retention/controller path.

| Policy/component | C0 / C1 / C2 native CPU seconds | Step envelope |
|---|---:|---:|
| CT15 V4 | 56.598205 / 189.308833 / 350.148692 | 160.839859 |
| V5 StreamScalar | 49.903598 / 157.136224 / 273.466898 | 116.330674 |
| V6 ScopedScatter | 53.181892 / 154.906970 / 253.475760 | 101.725078 |
| V7 ScopedAuditPoints | 46.089977 / 137.980507 / 234.472438 | 96.491931 |
| V8 ExactFlatMetric | 44.271072 / 127.618765 / 212.064898 | 84.446133 |
| V9 ImmutableRowLRU | 41.024357 / 116.387451 / 196.404393 | 80.016942 |
| V10 StreamGradient | 34.413750 / 92.753797 / 148.738315 | 58.340047 |
| V11 ScopedAuditResidual | 33.570767 / 90.702362 / 145.995540 | 57.131595 |
| V12 HashedCoverage | 32.361887 / 91.765629 / 147.166028 | 59.403742 |
| V13 CurrentVelocityGradient | 33.407704 / 96.760721 / 158.519554 | 63.353017 |
| V14 CurrentObservedScalar | 32.009020 / 86.142308 / 138.523160 | 54.133288 |
| V15 HashedBasis | 31.406938 / 84.155183 / 133.453709 | 52.748245 |
| V16 CurrentSampler | 30.849522 / 87.151439 / 137.786030 | 56.301917 |
| V17 ExactTensorIndex | 31.906881 / 95.315014 / 158.114100 | 63.408133 |
| V18 CentredGradientScratch | PROSPECTIVE; not executed | PENDING |

V5 streams scalar FE values with the original double basis products and centred
long-double sum, with bounded 128-entry scratch and a per-call constitutive
map capped at 2,097,152 points. V6 replaces repeated scatter hash searches with
a bounded, per-assembly pointer plan and preserves chronological additions.
V7 gives each independent native audit its own current point table, without
borrowing domain validity. V8 specializes only the strictly validated canonical
identity Flat map; general routes remain. V9 evicts individual least-recent rows
at the unchanged 1,048,576-row quota. V10 streams paired actual velocity
gradients through 3x128 scratch, preserving all double/long-double arithmetic.

V11 reuses one actual raw residual only among immediate equivalent consumers
inside one audit, with original-vector bitwise controls. It preserves constraint,
maximum and energy Lift/Refresh side effects. Actual residual calls fall from
21 to 15, but the whole native speedup is modest. V12 changes bounded coverage
membership to an exact hash with ordered spill and a shadow reference. The
frozen selected and small cases use LegacyV1 KInventory: **coverage insert and
query counts are zero**. This branch has no selected qualification or benefit;
its inventory was not changed to activate it.

V13 memoizes only pure velocity gradients after checking the complete actual
velocity prefix bytes at each read-only scope entry; exact geometry/element,
raw and resolved side and coordinate bits complete its key. Scalar/constitutive
validity remains fresh. Small flow-change invalidation and cold-gradient controls
pass. The actual one-step diagnostic has 5,422,654 hits, 2,115,730 misses, one
actual velocity invalidation and peak 1,176,905 entries. It is slower, so V14 and
later omit it. Hits alone do not demonstrate throughput.

V14 uses the already freshly observed K value only inside the same current
SingleSweepDomain point before trace-state construction, when effectiveK is
unchanged. Independent native audits still compute their own K. Every reuse is
compared to the cold original in the small controls. V15 changes only exact
immutable basis lookup to a hash, preserving first signed-zero-equivalent
insertion, construction, finite/local checks and the 8,192-entry whole-clear
quota. A separate actual selected one-step shadow run passes every old-map
presence/value-pointer check through **44 real cache clears**, 57,464,812 hits
and 376,674 builds. Its 87.871558 CPU seconds include verification and are not
performance training.

V16 uses an at-most-4,096-entry current K table only within one read-only
inventory and clears current values on entry/exit. Exact owner metadata may
persist only for the same strong native field; strict finite/containment/side
checks precede hits. The small cold owner/value controls pass. Actual diagnostic
value/owner hits are 4,875,096/6,530,075 with peak 3,473, but timing worsens.
The [prospective V17 amendment](cao-tamura-accepted-16-evidence/attempt-13-amendment-v1.json)
therefore omits V16 before the first V17 run. V17 adds only exact original
integer tensor packing triplets to the V15 base. V18 is an orthogonal V15-based
one-call centred coefficient-difference scratch hypothesis, without V17 or V16.
Its [contract](cao-tamura-accepted-16-evidence/attempt-14-contract.json) freezes
old-order three-axis reductions and cold binary64 checks before execution.

Completed small suites retain four d0/displaced P3/P4 operators and Off updates.
V8 and later include 44 ordered current-flow/K/epsilon/time/data/alpha/mask/space/
geometry queries. They retain
genuine 4 versus 2+2 continuation, full rejected-step rollback retaining step two,
foreign/unknown history rejection, and continuation with row and point quotas
both reduced to 1024. Scalar/LRU/map/gradient/basis/sampler cold verifiers cover
affected arithmetic and signed zero. Each policy has its own execution/history
identity. CT14/CT15 histories have never been relabelled as new components.

The measured V10 one-step diagnostic assigns 18.336474 exclusive CPU seconds
to residuals, 15.402599 to constitutive inventories, 8.416919 to physical audits,
6.559479 to flow assembly and 5.501287 to repeated-point native audits. It makes
21,870,516 scalar and 7,538,384 paired-gradient evaluations. These scopes guide
remedies; they are not whole-workflow acceptance timings. No factor reuse,
changed timestep, mesh reduction or altered physical equation is admitted.
The later V15 diagnostic has 13.542350 exclusive CPU seconds in 15 residuals,
12.247390 in constitutive inventory, 8.421522 in physical audits and 6.630562 in
flow assembly over startup plus one accepted step. Its accepted-step context
alone assigns 10.151360 to constitutive inventory, 6.630562 to flow assembly,
5.233350 to residuals, 4.133922 to physical audits and 2.473560 to factorization.
Thus the measured remaining work spans several native routines; this diagnostic
partition does not establish a performance lower bound or replace the timing-off
whole-command envelope.

Independent Python audit work uses only the sealed actual CT14 initial/coarse/
fine .002 states and makes **no physical advance**. It preserves independent
calculation from actual coefficients; native validity is never substituted.

| Independent audit experiment | Measured CPU seconds and interpretation |
|---|---|
| ExactTraceIndexV1 | Full bundles 16.325126/16.964765/16.221196 → 13.737221/13.807352/13.712643; all observations/decisions exact. |
| BoundedPhysicalValuesV2 | Cold bundles 13.524683/13.690482/13.628189 → 19.826120/19.719238/19.889716; slower cold path, scalar gradient witnesses retained. |
| IndependentBLASRowsV1 primitive | Per-row contractions .550475/.546867/.547607 → .016559/.017109/.015892; exact field/gradient bits, not a complete-bundle speedup. |
| ExactBLASAuditRowsV1 | Cold 13.518054/14.379557/13.707972 → 18.712579/19.446953/18.803082; warm 10.607173/9.966404/10.803283 → 8.257510/8.338147/8.292201. Construction/digests remain charged. |
| FreshBatchedA95InventoryV1 | Warm complete bundles 8.250353/8.262903/8.224565 → 3.168236/3.156945/3.149471; cold scopes 19.506479/18.867871/18.708416 → 14.421964/13.652345/13.703365. |

The trace index contains only 187 immutable first-match H-plus addresses
(19,132 bytes); missing/nonfinite/foreign/stale-geometry and scope-cleanup controls
pass. The first NumPy batch attempt failed the unchanged gradient-extremum
geometry witness and was not adopted. V2 retains original scalar gradients and
batches only physical-value norms; its maximum scaled differences are
`6.54e-16/3.77e-15/3.99e-15`.

The private C row helper freshly gathers actual coefficients and invokes the
same BLAS dot/matrix-vector reductions as the original NumPy arithmetic. Its
immutable plan holds geometry/rows/field-local IDs/offsets, not state. Invalid
IDs/shapes are rejected before writes; changed flow/K/epsilon, foreign packet/
row fallback, mutation rejection and failed-scope cleanup pass. Numeric plans
are capped at 256 MiB and 262,144 rows, width 128; actual plan bytes are
119,433,191. Point-local values exist only during the original read-only point
calculation; no values survive to another point, grid or state.

The fresh inventory rebuilds every actual outer-point tensor and invariant on
each call: full raw 3x3 gradients/divergence, strain, rotation, eigenvalues, norms,
tau, Cmu, radicand, covariance, R/Q, signed production and stress, plus the
explicit diagnostic planar invariant. BLAS nine-entry norm and the original
NumPy/Python hypot distinctions preserve scalar arithmetic. All limiting
contenders are reevaluated through the original scalar routine in original
order to preserve witnesses. Near-boundary, bad/nonfinite/nonpositive or
exceptional cases fall back to the **complete original scalar inventory**.
The `2e-10` algorithm-selection band is not a clipping or admissibility threshold.

All complete reported observation/decision/witness bits agree on all three
states. Additionally **every one of 26,042 outer points per state** was compared
for all 14 physical fields plus diagnostic planar radicand; maximum scaled error
is **zero for every field in every state**. Negative K, zero/negative epsilon,
pure rotation and both sides of the near-zero radicand boundary retain the
original decisions through scalar fallback. Initial/coarse/fine scalar witness
counts are 372/2/2; none of these actual-state calls falls back in full.
Qualification-only full tensor shadows cost 7.175284/7.320736/7.262660 CPU seconds
and are charged separately. Immutable preparation costs 6.508688/5.864940/
5.935019 CPU seconds; exit digests and cold costs are not removed from accounting.

The native and Python runtime identities are separate. V15 uses native source
`6756aabec32db4675574f25d5247d25cec86569f`, private runtime
`build-basis-v1/runtime`, governing SHA256
`ca8f21abf4e051ea9cd85bc54ae5b3d0a978451d9b07950354759e3f54321969`
and probe `8c14585f96e048c19732f748fdfa6bcfab378a653c369124819a4c6887a1fabc`.
V16 uses source `239a5426d8f620427b29dc211671557cc1ab7eb4` and
`build-sampler-v1/runtime`; V17 is built from
`3186446ccbc9a1f6a3ac69ef7306dcf30d54f6e0`; V18 is being tested at
`77bb66047e522a5e86d01cf3657d36fb365f0baf` (scratch implementation introduced at
`895f5417ee22e3f93ed82daf7e9dabfaa97cff12`). Publication/controller commits do not
imply a rebuilt native runtime. Full governing/probe/library hashes are in build
receipts. V10 source `5b4e3cf40520cfb361eb0f8eb13b5bd6a561b88b` and its prior
published runtime bindings remain valid only for that experiment.

The passing full exact-BLAS audit uses source
`9b19e63474d9155ba01a2f693688aa702e51f2c5` and private helper SHA256
`01c8b0b637d91597f39ec5f5876fa7078bb6393831103bb9e2ae7983a5ef69d5`.
The passing fresh inventory uses source
`65d04cf50e12966992742413351c50aeb7b29015`, helper SHA256
`8ad030566facfb6024271f4653c08ef46657db936f7cb914a8d6ce992edf9eee`, Python 3.10.12,
NumPy 1.21.5, and actual BLAS SHA256
`babdd0d2612d7bd2f0099ac8c82714702fb5bb47b1e2e6fdba819c97c79fc0fd`.
Compiles use one worker, `-O3 -ffp-contract=off -fno-fast-math`; complete actual
mapped Python/NumPy/compiler/BLAS/library bindings are retained, not inferred
from version labels alone.

Every unsuccessful attempt remains retained and charged. These include the
false paused-guard fixture match (fixed in the deliberate paused-v2 negative),
legacy probe-name and tiny-cache signature harness errors, the namespace error
in the first V11 build (96.976697 CPU seconds), and missing V11/V12 native-field
registrations. The latter passed fallback small tests but gave no intended
speedup; the owned training was interrupted via its exact external supervisor.
Its 138.689170355 conservative CPU witness is retained, with worker audit
**incomplete**, no numerical acceptance and no orphan. Constructor consistency
checks now prevent that registration mismatch. The deliberate paused-guard
negative is the other incomplete worker audit; neither is reported as PASS.
The first full exact-BLAS audit failed a tuple/dictionary attribute collision
in its tiny control before the three actual states; its request, library,
failure and costs remain. The fix renames the observation-row member.
Ineffective packed contractions, V12, V13, V16 and V17 are not silently discarded.
A V18 pre-launch clean-source assertion also stopped before creating a family
because Python had generated untracked bytecode. Those generated files were
preserved in `source-python-bytecode-v1` outside the private checkout, and
`PYTHONDONTWRITEBYTECODE=1` prevents new source-tree bytecode. The administrative
check/transport costs use the existing publication reserve; no solver ran in
that stopped launch.

The prior [progress-v2 snapshot](cao-tamura-accepted-16-evidence/progress-v2/results.json)
is historical progress through V10, before later costs; it is not final CT16
accounting. The final evidence snapshot and closed publication ledger will bind
all owned payload hashes, closed receipts including failures, actual runtime/
configuration/history identities and later attempts. Full CSR/state/history
arrays remain outside Git, with compact evidence and hashes in Git.

Conditional forecasts use decimal member-wise ceiling and length-32 segments
within the existing 64-step native hard cap, two separate zero initializations
for fresh starts, nonnegative native startup and step envelopes, and 25% headroom
**once**. The .407 fresh count is 1,628+3,256=4,884; own-.002 conditional starts
need 1,620+3,240=4,860 additional steps and are invalid until that exact component
has its own fresh qualification. No CT16 candidate has such a start.

For the current best measured V15, fresh-entry native-only CPU is
**328,113.1299625 seconds with headroom**, 13.6714 times the entire production
CPU ceiling. With startup held fixed, its 52.748245 CPU/step envelope must fall
below **2.934464498 CPU/step**, an additional **17.98-fold** improvement even if
every other cost were zero. The zero-startup/zero-other allowance is
3.931203931 CPU/step. The improved 3.168236 CPU warm audit already exceeds the
2.934464 allowance before native evolution, audit setup and all remaining work.
These are conditional short-workload extrapolations, **not mathematical lower
bounds for every future state, implementation or segmentation policy**.

| Horizon | Requested T; rounded member endpoints | Fresh member steps | V15 native-only CPU / wall with headroom (s) | All-native-files retained extrapolation (bytes) | Admission |
|---|---|---:|---:|---:|---|
| Entry | .407; .40700/.407000 | 1628 / 3256 | 328113.130 / 340442.297 | 923639884054 | UNADMITTED |
| 0.1 transit | 4.06216378576; 4.06225/4.062250 | 16249 / 32498 | 3274057.108 / 3397080.045 | 9218541418609 | UNADMITTED |
| 1 transit | 40.6216378576; 40.62175/40.621750 | 162487 / 324974 | 32739035.454 / 33969203.840 | 92183408306194 | UNADMITTED |
| 3 transits | 121.8649135728; 121.86500/121.865000 | 487460 / 974920 | 98216712.263 / 101907201.933 | 276549594729499 | UNADMITTED |

The all-retained extrapolation includes current wide native reports and is not
a minimal archive lower bound or a qualified retention forecast. The 24-GiB
production peak/new-retained ceiling cannot be claimed satisfied by Git
compaction. Complete CPU/wall, replay/retention duplication, final durable growth
and selected q64/headroom memory remain unknown. The old CT15 q64 raw peak
15.202 GiB passes 16 GiB only without headroom; 1.25 times that is 19.0025 GiB.
It does not qualify a changed CT16 q64 runtime. No complete combined campaign
measurement or heldout model is claimed from adding separately measured pieces.

Historical CT11–CT14 expenditure (43,163.053192 CPU seconds), CT15's separately
closed ledger (21,390.354157 including its exception; upper 21,630.354157 with
publication reserve), this CT16 engineering pool and hypothetical future
production remain distinct. Sunk research is not recharged as future runtime.
Sampled process `wchar` highwaters include pipes and can miss short-lived writes;
they are a cumulative-write lower-bound indicator, not exact disk I/O. Peak live
bytes count the private source/build/runtimes, all retained failures, histories,
temporary data and publication/transfer duplication. Publication reserves count
administrative source inspection, transport and Mac work once.

The smallest native workload reproducer is the selected zero/one/two controller,
not a full .407 run. The underlying commands below require a **separate
prospectively bounded allocation and its worker/package guard**, with a new
output root. CT16's fixed-start supervisor deliberately cannot be restarted
with a new clock, and this report does not authorize further numerical work in
a closed pool. Retain all six single-thread library settings, one rank, compiler
jobs two, process AS <=16 GiB, project <=16 and aggregate <=72. Source the
installed MPI environment and put `TMPDIR` under the new heavy-data root.

Use exact source `6756aabec32db4675574f25d5247d25cec86569f` with its matching
private `build-basis-v1/runtime` for the native experiment. Bind `s` to that
source's `solvers/IncNavierStokesSolver/Tests/RANS`, `ct14` to the sealed baseline,
`runtime` to that matched runtime and `new_root` to a new owned output root:

```bash
python3 "$s/ct_accepted16_basis_controls.py" --ct14 "$ct14" \
  --runtime "$runtime" --output "$new_root/small-basis-v1" --phase small
python3 "$s/ct_accepted16_basis_controls.py" --ct14 "$ct14" \
  --runtime "$runtime" --output "$new_root/training-basis-v1" --phase training
```

The original launch estimates are 350 CPU/450 wall seconds/100 MB new data for
small controls and 450 CPU/550 wall seconds/600 MB for training. Exact executed
supervisor commands, actual source revisions, numerical-library environments
and receipts are retained in the evidence. Private build reproduction uses
`ct_accepted15_build.py --phase runtime`; it preserves the shared installation.
For the independent inventory use source
`65d04cf50e12966992742413351c50aeb7b29015`, its matching test directory, the same
sealed `ct14` and a separate output:

```bash
python3 "$s/ct_accepted16_inventory.py" --ct14 "$ct14" \
  --output "$new_root/independent-inventory-v1"
```

The full inventory qualification estimate is 350 CPU/450 wall seconds/100 MB;
its helper compiler uses one worker. Do not overwrite retained outputs or
substitute a later publication source/binary without rebinding the experiment.
The report's source/runtime identities and complete request/loaded-library
receipts distinguish numerical artifacts from later publication commits.

The next technical capability needed is a substantially faster **native current
field/physics/assembly execution path**, plus a qualified bounded retention and
independent-audit workflow. A further local speedup is insufficient: the target
must account for startup, every required audit and replay before fresh full
qualification. No increase in resources or change of physics is assumed.
Task 3 can only be proposed after a component passes all G1–G6, with a separate
physical/uncertainty contract for the .407 pair, named own histories, fixed audit
milestones, rolling spent-plus-remaining guards and segment-boundary reforecasting.
Neither Task 3 nor later transit preservation is authorized by these partial
engineering results.
