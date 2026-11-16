# H25R1 resumed production qualification

Final status: **TARGET_NOT_MET**. See `wj-h25r1-closure.md` for the seven-criterion
matrix, actual pilots, complete costs and closed accounts. The sections below
preserve each checkpoint's state at publication.

At R0, Linux numerical/build admission was RESOURCE_BLOCKED. This was an
admission state, not a numerical failure or a production PASS.

The fetched GitHub head and local branch HEAD both resolve to
`58e26235120574110253a3df6329d5bbd110ee40`. H24 close
`a60d6b3ee716ea14d388e8bf6e27dff6eccd63a1` is a verified ancestor. All intervening
H25 commits were inspected. HTTPS fetch failed with an HTTP/2 transport error;
the explicit GitHub SSH fetch succeeded without changing the branch.

The initial read-only Linux inventory is retained in
`wj-h25-linux/resumption-r1/takeover-linux.json`: 643 current runtime/input/library
pin checks passed; all original operation reservations have final receipts,
and h350/certificate/pilot output paths do not exist. The single inventory
Python process is itself visible in the raw process snapshot; it is not a
pending H25 experiment. Original account usage is exactly 5125.204103887139 s
(native 2099.6596822449355, build 389.50274165614974, passive 1916.041679986054,
setup 120, closure 600). Attempts 22/200; advancing blocks 5/12. The account
is OPEN and no additional H25 operation has been recovered or launched.

Source-v8 native revision `54c4555c287a8759cbc69c8360c678899e842bcd`, build driver
`b6bf57f2`, build-scalar-v3 and runtime-scalar-v1 remain separate identities.
The actual ELF SHA256 values are `6a226c9a621a26cb2042959350fd69cae694ed56e573305b35384710cc1359a3`
(solver) and `897c0c8bf08c20da985606e6eaaeddd415884e04472dede0f2ba3467dd40477b`
(resolution probe). Both wrappers, builds, sources/objects/dependencies, runtime
library pins and original input pins passed current verification. Historical
native-exec maps remain evidence of the executed comparisons; a read-only pin
check is not a new native-exec map observation.

| Work item | Reused/remaining status at R0 |
|---|---|
| Stock S0 eight advances and explicit NoCollection recipe | Reused after 40 checkpoint-A transport files verified |
| W0 .04→.06 and matched exact lookup/scalar comparisons | Reused after 26 checkpoint-B transport files verified |
| Lookup/complete scalar/failure parity | Accepted historical results retained; scalar Off selected |
| h350 initialization and physical comparison | Not executed; old request RESOURCE_NOT_ADMITTED |
| Coefficient screen | Not executed; old request RESOURCE_NOT_ADMITTED |
| Fine own-dt transient seed | Helper/draft only; native qualification pending |
| Target h/q seed | Known incompatible ordinary-import path; repair/qualification pending |
| Evolved h/q/dt comparisons | Not executed |
| Selected cadence/evolving pilot/restart/terminal/storage | Not executed; archived draft is not evidence |
| Complete affordable production forecast/equation share | Not established |

New source/operation work follows the prospective separate allocation in
`wj-h25r1-continuation-contract.md`. Three unpublished drafts are preserved in
`resumption-r1/takeover-drafts/` with byte hashes in `draft-inventory.json`.
The main Mac checkout and Linux general checkout (f4b70bae with an unrelated
developer-guide modification) are untouched. Initial host allocation is 131
scientific/helper threads, above the mandatory 72 ceiling. No numerical job
is admitted while this persists; no unrelated work is terminated.

Checkpoint R0 is the commit containing this evidence and the continuation
contract. Its exact SHA is recorded by subsequent checkpoint records, avoiding
a self-referential commit identifier. Work continues after remote verification.

## R1: actual initialization and rejected admissions

R0 was pushed and verified at 6b25c30cb7e902ebc0e73048fed47d9ebeded175.
Host capacity subsequently fell to 15 workers, permitting serial Linux work.
The corrected immutable source-v1 overlay is commit
955600e1569feea0317814f6b15fea4d5d3ce9aa over preserved source-v8.
The target seed probe was built with eight workers in 60.812771 s; runtime-seed-v1
is distinct from the unchanged governing runtime-scalar-v1 used for h350.

| Newly executed item | Measured operation seconds | Scientific outcome |
|---|---:|---|
| h350 q12/16 fresh initialization |48.679771|12445 states, native startup completes|
| Common physical h350/fine comparison |93.974068|REJECTED; no evolved h qualification|
| Fine prescribed finite seed, dt=.0025 |120.959153|Correctly normalized and continuous; QUALITY_LOST|
| Fine unperturbed own-dt seed, dt=.005 |116.648065|PASS_INITIAL_DATA_ONLY; one genuine native history|
| Target h350 finite seed |70.456707|Native hard audit rejects constructed target|
| Repaired coefficient screen |1.776620|READ_ONLY_FEASIBILITY_NOT_ACCEPTANCE|

The 350-cell candidate fails U, k, stress, gradient and principal budgets; its
worst stress MAX uses 45.1544 times the limit. Its normalized principal minimum
is .0005799511113782629 versus .002068763923454. Wall traction MAX is
.000176245626 relative to u_tau squared, below the .001 budget. Restoring only
resolution implicated by saved failing cells is being investigated; blanket
coarsening is not selected. The first localization attempt hit a serialization
failure caused by classifying a rounded interface height without owner provenance;
its failure/partial temporary are retained and the owner-based repair is separate.

The prescribed fine perturbation now uses physical coordinates before projection.
Its velocity normalization error is 8.67e-19; relative k normalization error
7.69e-17. Changes in U/G/k/epsilon/current K/cutoff are all finite and nonzero.
The maximum regional Q1 moment is 8.78654e-19 and scaled strong divergence
8.84310e-17. Nevertheless its independent normalized principal minimum is
.00043589370996140187, below the frozen quality threshold. The native hard
physical audit passes for this fine state; positivity/covariance do not establish
the separate quality gate. No governing advance is admitted from this seed,
and no dt accuracy claim is made. Amplitude, support, closure and threshold are
not changed in response. The target h350 seed separately reaches target native
validation and fails its hard audit; no ordinary fine restart is imported into
the target patch. All source checkpoint bytes and ordinary signature rules stay.

The unperturbed .005 fine producer passes actual native and independent admission,
preserves original spatial/reference/defect bits, and initializes its own history
once with zero governing solves. This qualifies initial-data production only.
The eight protocol/scaling tests pass on Linux; the earlier Mac test attempt
could not import SciPy and is not numerical evidence.

The coefficient screen verifies its completed report/structure/source linkage,
signed native cubic layouts and platform extended precision (63 explicit mantissa
bits). All 700 k and 476 outer-epsilon cells pass its heuristic. The old 4096-eps
reserve is explicitly not a proved outward enclosure; no guard is skipped and
cell counts alone are not runtime coverage. The first screen failed because
step20 is an off-cadence completion without a ready barrier; the new immutable
tool binds its genuine completion manifest instead. That failure is retained.

Compact verified R1 evidence is in
wj-h25-linux/resumption-r1/checkpoint-r1-v1. Full arrays, reports, checkpoints,
logs and source/build/runtime capsules remain under the indexed Linux resumption
root. The original account is unchanged and OPEN; the continuation account is
a separate OPEN snapshot. R1 contains scientific rejections, not a production PASS.
Next are the frozen native prototype contracts in wj-h25r1-native-prototypes.md,
the localized mesh decision, a separate quadrature initialization and actual
retained-trajectory restart/terminal validation. The prescribed transient's
failed initial gate keeps subsequent transient h/q/dt and cadence admission open.

## R2: implementation, selective spatial and quadrature results

R1 was pushed and verified at `39ff5e9110dc0202dd0af9072a180e979a46f2b4`.
Native source-v2 is the overlay at that revision on source-v1; build-scatter-v1
and runtime-scatter-v1 are its executed governing/diagnostic binaries. Source-v3
at `a13d2dbf598b7cc137c13e9235835c99aecca36e` adds the separately tested memo
candidate. Neither executable is identified by this document's eventual HEAD.

| Native comparison | Ordinary complete mean, s | Equation share | Attribution |
|---|---:|---:|---:|
| New scatter executable, scatter Off (one sample) |41.210847021|39.5927%|99.9806%|
| ElementExactV1 scatter (three samples) |35.153428171|33.7401%|99.9816%|
| New memo executable, memo Off (one sample) |32.705812277|30.4530%|99.9803%|
| MemoExactV1 (three samples) |33.407643352|30.0673%|99.9894%|

All four comparisons genuinely restart H24 M16 and advance. On candidates pass
four exact coefficient/history/reference/callback comparisons to W0 and the
complete native audit map comparison; component scatter checks preserve ordered
sum bits, cancellation, signed zero and sparse-pattern behavior. Scatter is
retained. The memo passes 11,232 native value comparisons and 16 complete scalar
inventory/failure evaluations, but fails its prospectively fixed requirement to
beat BOTH Off and the scatter mean by 2%. It is retained as correct but slower
negative evidence. Variation between the two executable Off samples is substantial;
these short means are not confidence intervals or a pure isolated causal estimate.

The positive-lower-bound prototype executes its conservative native bound and
the original pointwise fallback on every sampled point. Moving-point cost
coverage is 100% for A, changed k, changed epsilon and returned A; the original
negative-k failure occurs before moving points and has no coverage claim. Bound
work is .496--.570 s versus 1.046--1.319 s for the corresponding exact point
contractions. These are measured moving-point costs, not a whole-step speedup.
Certified lower values and actually evaluated minima are reported separately.
No production guard/root/principal/current-K/cutoff/trace check is removed.

The actual 602-cell successor restores the original fine children in the 18
parent intervals selected from h350. All field budgets pass, worst utilization
.369811, but principal-action utilization 1.013965942 exceeds one. Localization
finds parent pair17; restoring it and neighbors17--18 along with existing pair16
produces the separately materialized 630-cell successor. Its new full comparison
passes initial representation only: worst field utilization .358143144,
principal-action utilization .981628980, and zero wall-traction discrepancy.
No evolved spatial qualification or performance saving is claimed for either.

The 700-cell 8/12 governing-quadrature producer constructs a genuine target
reference/projection/lift/FlatInletTrace/history and rejects an ordinary wrong-q
restart. Its initial physical comparison passes: worst field utilization
3.46555e-12, principal-action utilization 1.84722e-11, zero wall discrepancy.
This near-identical initial state is not evidence of evolved quadrature accuracy.
The operational control therefore retains 700 cells and governing quadrature
12/16. The frozen finite perturbation is still inadmissible; no .01/.02/.04/.08
transient dt escalation is justified from a rejected initial state.

The selected exact implementation completed eight genuine C advances at .005,
C4/.02→C12/.06, in653.247665 s including startup, barriers, final fields and
driver validation. C5--8 state/history/reference/callback bits equal retained
H24 C. Every step has one flow/k/bound/epsilon sweep, zero bound activations and
zero extra corrections. Fresh current physical audits occur at C6/8/10/12.
Terminal principal minimum is .0035485288018655137; maximum signed regional Q1
moment2.54039e-16 (limit2e-15). Strong divergence is reported separately.

Actual shared completed times .05 and .06 pass C/M physical comparison without
interpolation or pressure regauging. Worst field-budget utilization is .0224451,
worst directional-principal utilization .0369537, wall traction MAX2.39493e-6
relative to u_tau squared (limit .001). Earlier H24 C/M/F .02--.04 evidence is
reused only after exact new C5--8 parity. This is a short unperturbed-window
result, not transient or full-duration production qualification.

C-wide ordinary warm steps7/9/11 cost [42.362902401,45.387005606,36.034117008] s,
mean41.261341672 s, with31.8724% equation work and99.9920% attribution. Its
four scheduled scientific steps average94.805427824 s, including full native
audit, independent observer enclosing wait, publication and recovery; attribution
97.8313%. These broad current measurements, rather than the fastest isolated
scatter sample, enter the provisional complete-case forecast. Sparse/dense
restart controls, actual retention/strict read and genuine fresh cost remain for
checkpoint C. Their scientific limitation remains the failed finite seed.

The R2 complete provisional forecast retains T505.02466521278404, the fixed
terminal window, 101005 fresh/.005 advances (100997 from C8/.04), recovery each
step and scientific stride4/.02. It counts the measured fixed process envelope
twice for initial continuation plus one planned restart, and prospective600 s
genuine fresh construction/180 s terminal strict read until measured. Fresh
cost5520708.054 s (1533.530 h); continuation5519670.875 s (1533.242 h). With25%
contingency the fresh forecast is6900885.068 s, far above86400 s. Complete
equation share23.8956%, attribution99.0551%. Retaining full sampled scientific
records would need582560331988 bytes (about543 GiB), beyond the8GiB development
science allocation. These are conditional short-control extrapolations; a .25
cadence sensitivity is explicitly unexecuted and cannot supply qualification.
Actual fresh/strict costs and repeated sparse/dense block variability will replace
the provisional caps at C. No 505-unit run or production archive is launched.
