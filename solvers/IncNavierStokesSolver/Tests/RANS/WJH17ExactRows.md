# H17 exact native observation rows

`RANSWallImmutableRowReuse=PatchExactV1` is an execution-only, opt-in policy.
Absent/`Off` preserves the scoped owning catalog and by-value WJ audit accessor.
Both comparison modes must use `RANSWallEvaluationReuse=ExactStepV1`.

The H17 working controls, in SolverInfo, are:

```xml
<I PROPERTY="RANSWallImmutableRowReuse" VALUE="PatchExactV1" />
<I PROPERTY="RANSWallImmutableRowEntries" VALUE="196608" />
<I PROPERTY="RANSWallImmutableRowBytes" VALUE="268435456" />
```

A patch owns one exact-key LRU of immutable `RANSWallGlobalPoint` rows. It is
separate from EvaluationCache: state bytes, time, data generation, boundary data,
history, flags, profiles, cutoff, scalar domains and residuals keep their existing
invalidation. A patch/space/map generation cannot share a catalog with another.
No evaluated closure or admissibility result enters the structural catalog.

The key is block, native geometry/element, actual radial side and the IEEE bits
of all coordinates. Hits still validate native coordinates and field/map identity;
ResolveObservation retains its ownership/side/interface checks and ordering.
The original composite mapping scan and derivative multiplication are retained.
Only the WJ scalar-row accessor avoids a value copy under the new policy.

The byte cap includes row/vector capacities, control blocks, exact keys, LRU/map
nodes and conservative allocation overhead. Evicted shared rows remain charged
until the last consumer releases them. Deterministic eviction and oversized or
retained-budget bypass do not raise the cap. Bypassed temporary rows and ordinary
caller value copies are additionally covered by the aggregate RSS limit.
Counters report bounded request-gap observations, not an exact unbounded trace
or a unique-key count. Fixed surface/reconstruction catalogs predate this change.

`RANSWallResolutionProbe` supports two narrow operations:

- `immutable_row_checks`: a small, strictly imported saved checkpoint fixture;
  no initializer or governing solve. It extends existing exact reuse/adverse
  tests, including different native patches/spaces/maps.
- `complete_audit_batch`: exactly six genuine saved roles, absolute steps
  `0,1,2,3,4,0`, on the 24688-state target. It uses the same strict importer as
  `open_checks`; completed roles need their callback and native qualification
  commit. A fresh native-equivalent dynamic scope is created for each complete
  audit. No history advances or physical commits occur.

The batch emits complete native audit maps and exact binary commitments for
ordered row requests and observation values/gradients, retaining multiplicity
before WJ's existing seen filter. Each 4096-record chunk has a hash and bounded
first/last identity. These are little-endian integer/IEEE streams, not hashes of
rounded JSON. First-mismatch localization is a chunk, rather than a full saved
million-visit stream. Constraints, residuals and boundary outputs have separate
commitments. Timing/metadata/counters are excluded from numerical parity.

`study_wj_h17.py setup` prepares exact inputs, runtime and independent-cache
compatibility. `freeze` binds actual sources, binaries, comparator and policies.
`run_wj_h17.py` owns the new ledger and reservations. Its admitted small commands
are only A(two steps) and B(one restart); each is unique and counts attempts.
The target command is admitted only after saved-state parity, small integration,
strict reader and runtime qualification. It restores H14 CP0 and replays all four
H16 stages; it never advances H16 step4 beyond .01.

The H17 observer wraps the unchanged H15 numerical evaluator. Before an ACK, it
compares the current serialized state, native history, genuine callback,
reference/defect, time and native audit with H16, then compares the independent
summary and all 780 moments on the *same live* inventory. Earlier stages are not
seeded with the final witness union. Prepared PASS remains distinct from the
native qualification commit/cursor.

Authoritative commands and executed-versus-NOT_RUN statuses belong in the H17
contract/evidence and Linux receipts. Registering the two new CTests does not
mean they were executed. No new physical height, governing equation, trajectory
or scientific timestep is authorized by this facility.
