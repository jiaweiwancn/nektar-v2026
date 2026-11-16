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
