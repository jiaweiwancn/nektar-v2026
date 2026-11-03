# M6 evidence

The full numerical qualification is pinned to
`2d5d175f98dffa543892e453569b57109d6d55c7` on the isolated Linux source.
The physical-orientation initialization guard is subsequently checked at
`8a80594` with six direction fixtures and all three complete operator replays.
All native build/solver conclusions come from that environment over Tailscale.
Large raw cases, fields, checkpoints and logs remain under:

`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised`

- `m6-attempt-007`: complete M6 operator, manufactured, refinement,
  condensation, MPI/restart, standard-control and rejection qualification.
- `m6-kernel-001`: registered kernel CTest on the same revision; its independent
  result is also retained under `build/cls-epawf-m6-kernel`.
- `m6-quadrature-001`: same p=3 mesh/model/mask with tangential order 24,
  compared with the qualification's order 12 case.
- `m6-regression-001`: 28 earlier registered regressions, with 12 CTest worker
  reservations and numerical threads fixed to one.
- `m6-figures-003`: source-labelled, visually checked Matplotlib PNG/PDF.
- `m6-final-audit-001`: passed pre-guard integrated evidence, preserving the
  original 30-test coverage and physical-output audits.
- `m6-guard-001`: orientation acceptance/rejection and complete operator replay.
- `m6-final-audit-002`: final integrated evidence after the guard qualification;
  raw artifacts remain outside this repository.

The main qualification is direct execution of the same registered M6 entry
point. It is combined transparently with the registered M6 kernel CTest and
28 earlier CTests to cover the original 30 registered tests. The additional
registered orientation entry point brings final coverage to 31. The regression XML alone
contains those 28 earlier tests; it is not a claim of one 30-test CTest run.
`integration.json` records each execution and source revision separately.

`collect-integration.py` verifies source binding, every registered test,
frozen original laws/contract, native stages, physical VTU coordinates/field
values, trace omission controls above numerical noise and equation ownership.
`fingerprint-runtime.py` records all numerical source/configuration inputs,
every built executable/shared library, resolved runtime dependencies and the
numerical package/tool environment. Its strict comparison permits a final
documentation-only commit only when numerical inputs and runtime images match.
`qualify-quadrature.py` runs the independent fixed-mesh quadrature supplement.
`plot-verification.py` generates measured p/h/time errors, the declared
dt-squared guide and one-way mask from passed evidence.

The first figure render had overlapping minor tick labels; it was corrected
and rerendered. State-sensor plots cover their declared flat mask support;
outside that support the implementation skips sensor evaluation. Historical
failed build/MMS/energy diagnostics and the Ubuntu package-manifest adjustment
are retained on Linux and described in the M6 report. They are not counted as
passing qualification.

`finalize-orientation.py` checks the exact source replacement in the incoming
orientation guard, the constructor-only probe branch and the additional test
registration. All other numerical inputs must match. The changed initialization
guard is unreachable from the earlier Fixed-mode regression cases. All three
complete M6 operator reports are replayed and compared, including weak rows,
constraints, histories, directional actions and budgets. This narrow source
transition is explicit; the final documentation-only fingerprint comparison
allows no numerical-source or runtime exceptions.

The frozen contract's planned owner entries remain immutable snapshots;
M6 delivery is documented by the integration manifest. The original CLS
stability and general inner-k positivity failures remain physical failures.
These short manufactured runs do not establish long-time or experimental
hill-flow accuracy.
