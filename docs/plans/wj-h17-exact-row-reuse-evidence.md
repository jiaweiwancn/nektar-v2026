# H17 exact native row reuse

The fetched starting revision was066a1224cc1ed627bf2eaecdc780668a1cf68c11,
matching review. Source/Git work used the existing WJ worktree; unrelated main
checkout files were preserved. Linux used healthy ext4
`/home/jay/build/wj-h17-exact-row-reuse-20261114`. H15/H16 ledgers and original
payloads were not changed. This package adds no scientific timestep.

## Implementation and provenance

`PatchExactV1` is opt-in; default/Off retains scoped rows, the by-value WJ scalar
accessor and original composite mapping scan. Both modes retain ExactStepV1.
Only immutable FE/geometric rows persist across scopes. All dynamic state,
profile, cutoff, constraint, residual, closure and decision invalidation remains.
Every warmed request validates coordinates and field/map generation; existing
owner/side/crossing checks and WJ seen placement remain. Constant centering,
long-double contractions, summation and physical-gradient transforms are unchanged.

The deterministic catalog was frozen at196608 entries/268435456 accounted bytes.
Keys retain all coordinate IEEE bits, block/geometry/element and actual side.
LRU eviction charges held shared rows until release. Oversize/retained-budget
bypass and scoped-owner avoidance are tested. The optional mapping-lookup
optimization was not implemented.

Native C++ inputs are ad847461e8935dad936e6e9a8cf9300e696c7249; test-only CMake
registration was later updated. The executed compiler commands, header/object/link
bindings and loaded-library evidence are retained in the Linux records. Native ELFs:

| Executable | SHA256 | Bytes |
|---|---|---:|
| Governor |2235710d4237f5f390e8453943985690ba132fd28a82bab661e6f3e982bca0cf|8144904|
| Reader |63d62d9a603b9f7207dfab23a555512479ea90d4bbb68df74299f5b6b4412aca|4749960|

The private original runtime closure is reused, rather than a host-libc substitute.
Historical H16 native source5799878d and frozen Python2ddb1354 retain their own
identities. The independent211582940-byte cache remains unchanged at
800ac28c73ecb6498e6504b4e14e774db2cff893f2c0aaf418cb7036d50d00b0;
active numerical modules, profile, qh0, reference/defect and runtime are separately
bound. No independent immutable rows were rebuilt.

## Saved-state and focused qualification

The first build failed strict compiler checks; the corrections and failed charge
are retained. The catalog unit fixture passed. Six initial H17 Python methods,
the later coordinator CLI and dispatcher methods and five relevant existing H14-H16 methods
passed on Linux. Three mistyped test-class selectors failed before execution;
only those missing methods were then run. Registration is not execution.
The native read-only small fixture passed in102.139590864s, exercising exact
reuse/invalidation, adverse scalar/cutoff/constitutive states, warmed invalid
requests and separate native patches/spaces/maps. It used the separately retained
pre-batch reader; the subsequent C++ edit affected only the new passive batch's
outer dynamic scope.

There were exactly12 full-target passive complete audits:0,1,2,3,4,0 per policy,
serial, with genuine strict checkpoint/callback/commit imports and zero governing
solves, projection, initialization or physical commit. The new batch uses a fresh
native-equivalent ExactStepV1 scope for every role. Dynamic results are not carried
from A through B back to A. All four baseline completed audit maps also reproduce
H16 bitwise. All six baseline roles pass hard predicates.

Both policies match exactly in complete audit maps/witnesses; ordered IEEE
row-ID/coefficient/derivative and current-value/gradient commitments with visit
multiplicity; constraints, spatial/completed residuals, exterior/interface,
roughness and energy outputs. Comparison excludes enumerated execution metadata,
counters and timings. Chunks frame4096 records; first-mismatch localization is
bounded to a chunk and its first/last identities, not a saved million-visit stream.

| Saved batch | Charged wall seconds | Complete audits |
|---|---:|---:|
| Off, ExactStepV1 |528.350618548|6|
| PatchExactV1, ExactStepV1 |551.725580002|6|

These totals include strict import, exact commitments and affected-output checks.
They do not establish a speedup. Per-audit and commitment timing are published
separately. A supplemental inline check had a parse-only Python quoting failure;
its charge and erroneous setup-stop classification are retained with an explicit
source-error correction. No physical/coverage/runtime halt was cleared.

The new-policy batch recorded43,912,889 requests,22,941,210 hits,
20,971,679 constructions and20,802,834 evictions. Final residency168845 entries
used268434384 accounted bytes; peak live bytes equalled the268435456 cap, with
zero bypass and zero retained-after-eviction bytes at the endpoint. Constructed
native rows spanned688–1216 accounted bytes before per-entry overhead. Maximum
resident request gap was623232; the bounded histogram is not an exact stack
reuse-distance or unique-key trace. Peak residency is a distinct-key lower bound
and misses an upper bound including reconstruction. The byte cap bound residency. Baseline scoped reuse constructed18,912,605 rows
versus20,971,679 under the new catalog: construction was not reduced. Cold first
role costs were88.535s/89.733s (Off/new); return-to-A costs83.760s/88.876s.
These are identical operation roles, not repeated-state dynamic-result hits. No larger-cap or favorable rerun was made.

## Integration, replay and cost

All required integration and replay checks passed. Small A used H16 uninterrupted
producers; small B used its actual sealed restart producer. The observer-command
and producer corrections preceded any governing reservation. Their compatibility
record preserves the prior freeze and proves that the binary, caps and numerical
comparator/evaluator functions did not change. Prior ELF-map samples remain
separate from the active reader generation. A later CLI `runtime` dispatch error
was corrected/tested separately; the unchanged runtime checker was invoked directly.
The frozen governing capsule was not changed or rerun.

| Executed category | Actual work | Outcome |
|---|---:|---|
| Full-target passive batches |2 processes /12 complete audits|exact parity PASS|
| Small A uninterrupted |2 attempts /2 committed-qualified|H16 parity PASS|
| Small B restart |1 attempt /1 committed-qualified|sealed restart parity PASS|
| Small strict checkpoint view |1|PASS, zero initialization/governing work|
| Target replay |4 attempts /4 completed /4 committed-qualified|all H16 stages match|
| Final strict checkpoint view |1, absolute step4/T=.01|PASS with genuine callback/commit|
| Fresh CP0 / new scientific steps |0 /0|unchanged scope|

Every completed sweep had one U/pressure/k/epsilon solve and one bound call;
zero bound activations, material corrections or extra restoration occurred.
Original qh0/reference/defect, native history, callback and time labels matched.
The live inventory reproduced11451 points/17513 owners at stage1, then11456/17521
at stages2–4. All780 signed moments and separate principal-quality gates passed.
Independent principal minima were .003517183359, .003520644292, .003524039984
and .003527178522 at stages1–4. The largest signed-Q1 absolute moment over those
stages was3.896342854e-16, below2e-15. Native/quality/coverage remain finite
inventory claims.
The H16 final-union physical norms/fetch/exterior evidence is reused only through
explicit exact-state/unchanged-math bindings; it is not a new H17 scan.

The target restores original H14 CP0 and can only replay H16 absolute stages1–4.
State, native history, genuine callback, reference/defect, time, numerical audit,
live owner inventory and all780 moments must match before each ACK. Checkpoint
suffixes are generation indices. Prepared PASS and consumed ACK are not native
qualification commits. The protected final strict read follows process exit.
No new scientific progress, height, mesh, timestep, model or boundary is introduced.

SSH was intermittent. Failed/unknown network work remains charged; no reserved
native or passive sequence was relaunched. Subsequent governing coordinators are
transport-independent while retaining the existing ledger/process-group limits.
Publication snapshots and the final closed accounting are separate.


Measured warm means (stages2–4), seconds:

| Scope | H16 | H17 |
|---|---:|---:|
| Native compute before publication/wait |83.0564|90.7497|
| Independent assessment |6.8785|7.0044|
| Inclusive native.step |90.6284|101.1743|
| Native acceptance audit |40.5979|41.9840|
| Owned row/gradient bucket |14.2708|16.3988|
| Full principal search |11.6552|11.6521|
| Inclusive PIMPLE sweep |38.0133|44.3464|

**Correctness/failure equivalence PASS; cost reduction NOT DEMONSTRATED.**
The observed inclusive loop is11.6% slower than H16; this is one cross-run
comparison, not a timing significance certificate. Exact saved-batch comparison
also showed no end-to-end improvement. The fixed cache churned heavily; no
capacity increase, alternative data structure or favorable replay was tried.
Validation/lookup/construction/copy remain combined within accessor timers;
current-state contractions are separately timed. Their individual subcosts are
UNRESOLVED, so churn is a measured limitation rather than a unique causal proof.

The conditional costs are212.18 native-only days,228.56 compute-plus-assessment
days and236.55 inclusive-loop days for202010 steps. H16's inclusive figure was
211.896 days. Publication/wait/commit are already in native.step and are not added
again. Startup/terminal work and future witness growth remain outside this simple
linear forecast. Duration admission remains CLOSED.

Build145.411s; passive/small qualification1239.748s; target491.796s;
protected final read41.515s. The ledger closes at3014.99410465901/7200s,
including fully charged setup60/closure300 and conservative network charges.
Sampled aggregate RSS peaks were7,944,921,088 bytes (build) and7,066,984,448
bytes (science). Source/build/runtime footprint185,076,457 bytes; new output
274,291,472 bytes before archives. All limits passed. Failed compiler/selector/
parse/dispatch operations retain their receipts and charges.

Essential transport was verified first. The optional5,647,393-byte full bundle
(SHA256437c4998471ddf9d8cadb76430b98ac101667649914e57c999c535fedadd93e7)
was not admitted: its measured25%-headroom projection exceeded the fixed closure
allowance. Selected compact evidence is separately hash-bound; full chunk lists,
matrices, fields and fetch/exterior arrays remain at their sealed Linux paths.
The selected2,452,774-byte bundle timed out at205s with2,162,688 bytes retained.
Complete members were verified against both original manifests; no large resume
was admitted. Essential results, runtime/source bindings, actual replay commits,
current signed stage/strict data and closed accounting are available locally.
Missing optional files and full ordered commitments remain explicitly remote-only.
Two large verified JSON records use lossless gzip storage aliases. The figure's
long y label is edge-clipped; its meaning is audit wall minus commitment time(s).
The bounded label-only transfer timed out and was not retried; CSVs remain exact.
These are transport/presentation limitations, not scientific qualification failures.

## Interpretation and deferred decisions

The full-height neutral control, actual qh0 and internal frozen reference/defect
remain distinct. Exact numerical parity does not establish exact discrete
equilibrium, global admissibility, temporal convergence, stationarity, duration
preservation or experimental validity. H11 temporal qualification remains
UNRESOLVED; H12's bounded audit remains PASS.

H_w remains6mm, with physical flat interface9.5mm. The7mm stress observation
remains inside AWF. The primary inverse requires H_w>4.413515542mm whereas
placing every7mm observation outside requires H_w<3.5mm: the interval remains
empty. Finite-depth completion, roughness-sublayer validity and measurement scope
remain open. No physical height is selected. Hill/EP, altered inner stresses,
old measured-rms/20mm work, new timestep levels and longer evolution are deferred.

The preserved horizon is505.02466521278404, terminal window
[447.8483832973273,505.02466521278404],202010 steps. Conditional native-only and
inclusive-loop forecasts use measured warm replay timing; nested scopes are not
added and audit CPU0 is unmeasured. Cost correctness is separate from affordability.
Keep the default Off route. The next bounded package should isolate the measured
PIMPLE/row-access regression on immutable saved states: separate field-generation/
coordinate validation, key/LRU lookup, allocation/copy and composite-region scan
costs while retaining exact outputs and failure order. Do this before choosing
another data structure or cache size; do not chain another governing optimization
replay automatically. Temporal adequacy/application error budgets remain an
independent prerequisite for full-domain preservation, then terrain qualification
and the hill experiment. No follow-up is executed here.

## Reproduction and stored evidence

Use the exact commands in `publication-v1/operations/*/request.json` and closed
receipts. The ledger refuses another native reservation; these commands document
the executed allocation, not permission to rerun it.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
# Focused build as recorded:
cmake --build /home/jay/build/wj-h17-exact-row-reuse-20261114/build \
  --target IncNavierStokesSolver RANSWallResolutionProbe WJImmutableRowsChecks \
  --parallel 8
```

The build used RANS_BUILD_SOLVER=ON, GNU11.4 Release, -ffp-contract=off and
-fno-fast-math. The native catalog executable ran through the qualified private
loader. Python verification was the six-method H17 module, five selected existing
methods, then the two new CLI methods under separate charged receipts:13 distinct
methods total, not a broad historical campaign. Native reader requests are retained
at `science/row-fixture/request.json`, `science/audits/{Off,PatchExactV1}/request.json`
and `science/{small/A,target}/view/request.json`. Target source/intent, checkpoint,
callback and qualification bindings are in the versioned manifests. Mac source
includes the later isolated CLI repair; actual governing Python pins are those
in the frozen Linux execution-config and compatibility records.
