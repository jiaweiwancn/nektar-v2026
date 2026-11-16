# H25R1 resumed production qualification

R0 status: ACTIVE, RESOURCE_BLOCKED for Linux numerical/build operations.
This is an admission state, not a numerical failure or a production PASS.

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
