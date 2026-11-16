# H25 production strategy: checkpoint A

H25 remains ACTIVE; production readiness is NOT_ESTABLISHED. Checkpoint A closes
the fair short reference comparison and cost attribution, then work continues
with implementation and discretization qualification. The fetched starting head
was `a60d6b3ee716ea14d388e8bf6e27dff6eccd63a1`. All H13–H24 evidence and closed
accounts remain unchanged.

Authoritative root: `/home/jay/build/wj-h25-production-20261116`. The compact
checkpoint-A transport manifest records exact Linux paths, uncompressed hashes,
stored hashes and byte counts. Its accounting is an OPEN publication snapshot,
not a final closure. Source/build/runtime records distinguish the historical
H22 executable from H25 profiling capsules `27a4bb51`/`1bab7628`, and the
`0dae6510` executable whose new lookup policy is Off for this baseline.

## Completed equal-resource comparison

Both rows use one serial native worker, six numerical-library thread controls1,
dt=.0025, eight real advances over [.04,.06], Q3 velocity, and the actual
14×50/700-element PHYSICAL mesh. The materialized stock coordinates retain the
9.5mm interface and physical measures497.8125/3.75/494.0625. Stock uses supported
global continuousQ2 pressure (2929 DOFs), versus WJ's regionalQ1/interface
duplication (780 pressure DOFs); Q3 velocity has12986 DOFs. Stock uses12×16 GLL
points; WJ has8×8 native representation and12×16 governing quadrature. These
are distinct settings, explicitly retained in the case metadata.

| Measurement | Stock S0, no RANS/SVV | W0 unchanged WJ |
|---|---:|---:|
| Advances | 8 | 8 |
| Warm ordinary mean, seconds | .0221339571 (steps2–8) | 42.9757524 (steps17–23) |
| Full operation, seconds | 2.4356974 | 521.7935729 |
| Cold first stock advance, seconds | .377459 | startup/audits separately recorded |

The raw samples retain two successive advancing warm segments. W0 restores the
genuine qualified H24 M16/.04 state, with zero projection/lift/initializer calls,
and ends at M24/.06 with fresh startup/terminal native and independent physical
checks. Bound activations are zero. Two profiling replays (four and two advances)
match the corresponding W0 coupled states, fixed reference/defect, native
histories/time labels and actual completed callbacks bitwise. Those six replay
advances add no additional scientific horizon; the aggregate WJ attempt count at
checkpoint A is14 across three advancing processes.

Stock has positive molecular viscosity throughout, different supported top and
outlet conditions, a different global pressure space, and no scalar/AWF equations.
Its projected analytic nonzero initial velocity also differs from W0's retained
.04 field. Thus their elapsed gap is not identified as pure RANS equation cost.
The stock log's printed `L2 error` uses the stock default zero exact solution;
it is not an analytical-equilibrium error. This eight-step S0 result establishes
a finite cost reference, not a qualified production NS timestep.

S0's first launches aborted before any advance in automatic collection tuning.
The installed solver reproduced the same startup failure. `SessionReader`
requests opt-file updates when no opt file exists; `ExpList::CreateCollections`
then invokes `SetWithTimings` even after setting a NoCollection default. A pinned
explicit `.opt` file passed through `--use-opt-file` prevents that unwanted
autotuning and S0 completes. No SVV, mesh, dt, physics or tolerance was changed.
The original aborts and timeout charges remain recorded. They are execution
failures, not diagnosed CFL/aliasing or RANS constitutive failures; S1 is not
justified by them.

The independent common-operator fixture uses the SAME constrained12-row native
Q3 affine mass matrix,12×16 quadrature, coefficient bits and changing RHSs.
Native versus mapped solutions differ by5.72e-14; mapped residual2.71e-20.
For2000 applications each, native dense solve took.000568s versus mapped.010069s;
setups.0000166/.0001766s. This exposes small-operator infrastructure overhead;
it is not equivalence of the full stock and WJ systems or a production timing.

## Disjoint attribution and remaining work

The resolved Off sample at ordinary step17 totals43.615738s. Exclusive scopes
close arithmetically to the inclusive step;99.9828% is assigned. Accumulated
domain/principal buckets are nested explanatory measurements and are not added
to their parent totals.

| Disjoint purpose | Seconds/ordinary step | Share |
|---|---:|---:|
| Equation evaluation, assembly, solves and required updates |14.723179|33.76%|
| Guards, bound work and residual verification |25.668282|58.85%|
| Data allocation/invalidation/release |2.593077|5.95%|
| Recovery/publication |.315000|.72%|
| Attempt evidence/orchestration |.209675|.48%|
| Loop/telemetry |.099029|.23%|
| Unattributed |.007495|.017%|

The new release scopes identify2.569601s of the previous broad callback remainder
as scoped evaluation-data destruction. Domain checks cost19.529177s, with
moving-ray row work dominant (roughly15s/step); actual current-K profile
construction is much smaller. Flow assembly6.284408s, pressure-build own
work2.136659s and804 band applications1.728406s remain significant equation costs.
The finite sample is not a statistically established production mean.

Instrumentation is opt-in and excluded from the physics signature. Matched
advancing steps17–19 average43.457919s in W0 and43.769488s in the first profiling
build; paired changes range−.7451 to+1.1592s. This bounds the observed effect of
extra instrumentation together with cross-run variability; it does not isolate
clock overhead or establish a speedup. Every raw timing stream is retained.

Ranked work is: (1) exact element/geometry lookup for moving-ray FE access;
(2) avoid unnecessary scalar-row materialization/release while preserving all
checks; (3) element-local assembly and consistent-mass infrastructure;
(4) separately qualified mesh/quadrature/dt and monitoring cadence. The first
candidate is implemented behind a default-Off policy and its native owner
fixture passes; no performance benefit is claimed before its On trajectory.
No H17 row-result cache, closure-result cache or reduced guard inventory is used.

The >50% equation-share and24-hour complete-case targets currently fail by a
wide margin. Checkpoint A is published so subsequent changes have a reproducible
measured baseline. H25 continues; no production or hill evolution is admitted.


## Checkpoint B: exact implementation qualification

The execution policy `RANSWallNativeLookup=TensorExactV1` preserves the original
containment, owner, side and roundoff-fallback rules and changes no FE contraction.
The pure native index fixture and all four evolving replay states pass. Coupled
state, native history/times, genuine callback, fixed reference/defect and terminal
complete native audit map agree bitwise. All physical guards stay active.

| Same 700-cell case, steps17–19 | Ordinary compute mean (s) | Equation share |
|---|---:|---:|
| Instrumented policy Off |43.426694|33.75%|
| Exact interval/geometry lookup |39.010455|37.28%|
| Lookup plus scalar value path |42.872716|39.85%|

The lookup gives a measured10.17% reduction in this short matched comparison.
Cross-run variability remains; these are three advancing samples, not a production
confidence interval. The scalar path passed11,232 native value comparisons and16
full inventory/failure evaluations (changed k, epsilon, negative/nonfinite states,
return to A with fresh generations). Its four evolving states and native audit
map also match exactly, but it is slower. The selected recipe leaves
`RANSWallScalarEvaluation=Off`; the optional path is not presented as a speedup.

The executed native capsule is54c4555c287a8759cbc69c8360c678899e842bcd;
its build-driver correction is separately pinned to b6bf57f2. The checked ELF,
objects/dependencies, compiler/flags, wrappers and actual native-exec library maps
are in checkpoint-b-v1. Failed build attempts remain recorded. Peak science RSS
for the lookup replay was6,762,688,512 bytes. The complete dependency proof stays
on Linux with its exact size/hash; compact native comparisons and operation
receipts are copied without rewriting their identities.

The checkpoint-B account snapshot is OPEN, not closure. The h350 and coefficient
screen uploads succeeded but their requests were not admitted when other server
work exceeded72 threads; neither launched numerical work. Those admission failures
and the unsuccessful larger transport are retained and charged. H25 continues;
the production cost target and end-to-end strategy qualification remain open.
