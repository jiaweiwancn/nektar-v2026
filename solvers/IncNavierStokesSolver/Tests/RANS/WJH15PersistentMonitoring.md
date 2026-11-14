# H15 bounded retention and persistent assessment

This is an opt-in serial `MappedBlockV1` execution facility for the unchanged
6 mm H13/H14 neutral control. Default reporting and physics signatures remain
unchanged. The H15 target stopped at its first ACK: see the contract, protocol
stop addendum and evidence. Its closed ledger must not be reused.

`RANSWallRetention=BoundedV1` and `RANSWallMonitorIntent=<JSON>` are SolverInfo
controls, not new physical Parameters. `RANSWallDiagnostics` stays on. The intent
binds original CP0/qh0, source/runtime and observer signatures, nonce, output
directory, observer PID, timeout and maximum absolute step. Existing restored
reference/defect and native integrator history are never reconstructed.

`RANSWallMonitor.h` streams exact factor definitions and every call. The LRU
retains at most16 definitions and32 MiB of structural keys, including actual
nonzero pattern, restriction and source/sign scope. Pattern changes due to
numerical zeros produce new definitions. Eviction does not suppress evidence.
Native state/audit/stage/step/PIMPLE retention is bounded to current committed
and pending work; explicit rows have call counters. A failed Advance streams
its first candidate/history and provisional context before restoring the native
transaction. A bounded exceptional report never invokes legacy monolithic output.

After successful Advance, native code writes a unique checkpoint, actual
completed-callback sidecar and `WJCompletedStepV1` report. The atomic
`WJCompletedReadyV1` manifest is published last. The sidecar's callback history
is distinct from the current state and native GLM history slots.

The independent observer verifies these bindings and evaluates the full H14
inventory with fresh state-dependent closure, K/cutoff, epsilon, reactions and
780 Q1 moments. Only immutable geometry/basis/test/reference rows are cached.
Required observations survive variable-row LRU eviction. The private cache is
hash-qualified against source, physical configuration, original qh0/reference/
defect and runtime. It is not a general untrusted-pickle interchange format.

A durable PASS result and prepared cursor precede `WJObserverAckV1`. Native code
validates the complete17-field binding, consumed ACK and inventory before writing
`WJQualificationCommitV1` and advancing `qualified-cursor.json`. A prepared
observer cursor or PASS file alone does not make a checkpoint restart eligible.
Missing, failed, stale, malformed, duplicate or timed-out ACKs stop progression.
Only final names ending in `.ack.json` are published ACKs; atomic temporary names
are ignored until publication. The deterministic dual-name fixture fails the
old scanner and passes the corrected one.

The strict reader requires the native qualification commit plus genuine callback
sidecar for a bounded completed producer. An unqualified replay cannot bypass
that requirement. H15's final strict view therefore reads its actual last
qualified state, original H14 CP0, under the original producer with no completed
callback. It performs zero Project/Lift/initializer/governing work.

Actual commands, environments and input hashes are in the versioned operation
requests. The focused build uses `RANS_BUILD_SOLVER=ON`, GNU11.4, Release,
`-ffp-contract=off -fno-fast-math` and at most8 compiler workers. Native runs use
one rank and the six numerical-library thread controls1, through the qualified
private loader/library wrappers. Initial target binaries and later compile-only
repair binaries have separate identities; no repaired target was executed.

Relevant executable controls are `test_wj_h15_monitoring.py`,
`test_rans_wall_monitor.cpp`, and `check_wj_h15_protocol_repair.py`.
`study_wj_h15.py` handles source-bound fixtures and admission;
`wj_h15_observer.py` handles persistent assessment;
`finish_wj_h15.py` derives passive final-support views; and
`publish_wj_h15.py` seals compact evidence. They do not authorize another H15 run.

`native.step` retains its historical inclusive timing. Ready-packet
`native_compute_seconds` ends before publication and observer waiting. The new
audit buckets are disjoint wall intervals with CPU explicitly unmeasured.
Analytical-reference error, drift from actual initialized qh0, and fixed native
reference/defect diagnostics have separate roles and scales.
