# Implicit Newton history foundation — prospective checks

This is a tiny MMS history, retry and recovery test of the unchanged active
v2 equations. It does not qualify physical cold restart, terminal handling,
production monitoring/storage, H27-D or the full H27 campaign. It creates no
new physical-model approval. All numerical execution is Linux-only and remains
under the existing operation/resource admission.

Run `h27/implicit_history_check.py --manifest MANIFEST.json --table TABLE.json
--out NEW.json --seconds CAP`, with `0<CAP<=180`. Output is create-only, one
worker, with an 8 GiB address-space ceiling. The checker uses the independent
binary reader and signed-zero-preserving JSON parser in
`history_recovery_check.py`; it does not reuse a native checkpoint writer.
Both checker sources, the table and every consumed artifact are hash-bound.

The manifest schema is `WJH27ImplicitHistoryFoundationV1`. `root` is optional
and defaults to the manifest directory; relative artifact paths resolve under
it. `runs` maps `control`, `prefix`, `cold`, `retry`, `recovery3` to completed
`.result.json` paths. Matching `.steps.jsonl` and `.final.bin` paths follow
the same prefix. `input_checkpoints` maps `cold` and `recovery3` to
`{path,sha256}`, recorded before each dependent launch. `source_identity` may
add an expected native source SHA. `negatives` contains exactly three records
with `name,category,exit,log,expected_marker,session`; iteration/deadline records
also require `failure` and `checkpoint`. `checkpoint` is a path or
`{path,sha256}`. `session` is the exact launched XML, not a reconstructed
description. The controller writes the final manifest once.

The fixed positive fixture is PolynomialScalarsV1/body1, nx=nr=2, dt=.02,
linear tolerance 3e-13, current nonlinear tolerance 3e-12, Anderson depth zero,
and CoupledJFNKV1 with a five-second cooperative step cap. The admitted Newton
policy is 20 iterations, dimension40 and FD step1e-6. Positive signatures must
agree; the step cap is independently checked from reported configuration.
The control advances four steps; prefix advances one; cold resumes prefix1
for three steps; retry advances four with one injected failure at step2;
recovery3 resumes the retained control step3 for one step. All except prefix
end at step4/time .08. Every accepted step is logged.

Require all six reported current residual components <=3e-12 and strictly
positive reported k, epsilon, covariance and principal guards at every accepted
step. Bind current and previous states, current and previous explicit arrays,
pre-step histories, BE duplicate seeds, essential lifts, complete signed maps,
source/reference/signature identities and clocks. Compare final full coefficient
arrays and both older histories against control at MAX <=1e-8 with fixed
scales Uref=1, pressure tau0, k Kin, epsilon Ain/.5. Do not shift pressure or
normalize by a current RHS. Exact hashes are also reported; cross-run bitwise
identity is not required. Explicit pressure-load entries must remain zero.
Cache checks establish byte binding, correct history shifts and cold/retry
equivalence; their advection equations and the full PDE remain separate
independent-oracle responsibilities.

Control recovery/science stride1 produces four checkpoint writes. Require
recovery1/2 actually absent, recovery3/4 actually present and reported, and
science1..4 retained and hash-bound to accepted logs. Independently parse both
retained checkpoint payloads/checksums and bind them to science states. Require
the imported recovery3 file to have exactly the identity pinned before launch.

Require exactly one accepted-step retry, at step2 of the retry run. Its
`counts_complete` must be false after the injected interruption; live residual
and preconditioner callback counts and at least both sets of three block
rebuilds must survive. Other accepted steps require complete statistics.
All disjoint categories must be finite and nonnegative up to a 128-epsilon
elapsed-time rounding reserve; their sum must equal recorded step/attempt
elapsed time within that reserve. Report attributed fractions without treating
this tiny experiment as the >=95% warm-cost gate. Inclusive nonlinear traces
are not added to disjoint categories.

Negative categories are `iteration_bound` (one Newton iteration, tolerance
1e-14), `deadline` (step cap1e-9), and `algorithm_identity` (Picard reading the
Newton prefix checkpoint). Verify those settings from the actual XML and
require nonzero exits with the declared markers. The algorithm control must
fail on scientific-identity mismatch. Its initialization rejection has no
accepted state or required failure journal.

For iteration/deadline failure, require a failure journal and independently
parsed checkpoint retaining step0/time0. Compare all four fresh physical and
explicit histories against control science1's pre-step seed. Deliberately
changed numerical options must change the signature; source, physical model,
maps and reference remain bound. Failed seed residuals are explicitly steady
diagnostics, not accepted time-step residuals, and are not tested against
3e-12. Verify no retry, accepted log record, successful result or final
checkpoint is produced. Require failed-attempt clocks, disjoint times,
rebuilds and the statistics-completeness flag. Zero engine work is allowed
when the deadline fires before entry. Root/range/scientific failures are not
reclassified as recoverable linear failures by this test.
