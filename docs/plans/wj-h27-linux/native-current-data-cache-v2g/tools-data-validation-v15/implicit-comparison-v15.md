# Implicit algorithm and Newton cache comparison gates

This prospective test compares accepted normal-main v2 outputs from
`PicardV1` and opt-in `CoupledJFNKV1`. It changes no PDE, reference table,
current-equation acceptance threshold or history scheme. A comparison PASS
is not an independent equation audit or full H27 qualification. The existing
physical transient or coupled MMS oracle must separately verify both routes'
complete current equations.

The optional `--comparison` selects `PicardNewtonV1` (the unchanged default)
or `NewtonCacheV1`. The latter compares two `CoupledJFNKV1` trajectories:
old policy `RebuildEveryNewtonV1`, new policy `OperatorDriftV1`. Both modes
retain every fixed coefficient, physical field, trace, cutoff, explicit-history
and identity gate below. No pressure alignment or tolerance relaxation is
introduced by the cache comparison.

Use `h27/implicit_algorithm_check.py --old OLD.result.json --new
NEW.result.json --table TABLE.json --out NEW-CHECK.json --seconds CAP`.
Linux-only CAP must be positive and <=180 seconds; the parent admits and runs
the operation. The checker has one numerical-library thread, an8GiB address
space cap and create-only output. It launches no solver. Failure artifacts,
steady solves, mismatched clocks and missing exported history are rejected.

The exported state must identify the old/new implicit solvers explicitly and
include `coupled_newton` statistics. These statistics are provenance only,
not proof of convergence; Picard's Newton counters must be zero. No Newton
status field is presumed. A completed registered `.result.json`, accepted
step/time and the separate current-equation audit establish the relevant
completion evidence.

For `NewtonCacheV1`, require identical finite
`coupled_newton.relative_drift_limit` in [0,.25], identical positive finite
`step_seconds_cap`, identical
positive integer `maximum_uses`, and identical guard descriptions. The
new policy changes only reuse of the preconditioner; actual residuals,
Jacobian actions and physical acceptance guards must remain fresh. Metadata
comparison alone does not verify those equations.

A zero drift threshold is supported: exactly unchanged operators may reuse,
and every nonzero operator change refreshes. The zero-threshold native
adversarial control is separate from the paired equivalence/reuse test and
need not demonstrate positive reuse. Ordinary paired comparisons use .02;
`NewtonCacheV1` still requires positive aggregate reuse to establish that its
tested optimization was exercised.

Require sibling `PREFIX.steps.jsonl` files for both completed
`PREFIX.result.json` inputs, written with InfoSteps=1. They contain bare
native state records. Verify consecutive positive steps, clocks, source,
reference, signature and configuration bindings, and the accepted state and
explicit-advection hash chain. Require matching old/new logged step intervals
and exact final-record/result hash and `coupled_newton` agreement. Hash both
logs. Aggregate each per-physical-step record once; the final result is a
duplicate, not a second observation. A run resumed after step zero is reported
by its actual logged interval, without claiming evidence for an unprovided
prefix.

The `refreshes`, `reuses`, `maximum_used_age` and all seven
`refresh_reasons` counters must be nonnegative integers. Reasons are exclusive
(`injected`, `retry`, `control`, `cold`, `alpha`, `drift`, `age` precedence)
and sum to refreshes. The first prepare/use after factorization has age1;
every reuse increments the age; it must never exceed `maximum_uses`.
An unused step has maximum age0. A step with reuse must report maximum age>=2.
The rebuild control must have zero reuses and age<=1; only its control,
retry or injected refresh reasons may be nonzero. OperatorDrift must have
zero control-reason refreshes. Across the logged interval both routes must
have positive refresh counts, and OperatorDrift must have positive reuse
count. A final step with zero refreshes is valid when the interval establishes
the required work. The checker records per-step counters and totals, without
assuming that callback application counts equal factorization counts.

Require exact equality of the model, scientific source and qualified v2
reference identities; Modified_A basis, full Q3/Q2 topology, geometry IDs,
element bounds, signed maps, free masks, gauge flag, geometry and explicit
scales; step/time/dt, manufactured variant/body data, quadrature, convergence
tolerances, constitutive backend, stress split and unchanged configuration.
The algorithm signature and algorithm statistics may differ. Preserve native
pressure levels without fitting or subtracting a constant, including cases
with an explicitly reported canonical pressure gauge. Report pressure means;
never introduce a gauge into an open-traction physical case.

Verify each file's state and current/previous advection SHA256 bindings and
the accepted history shift. Compare all coefficients of the accepted state,
the pre-step current and previous states, current/previous accepted advection,
and pre-step current/previous advection. At BE startup, separately verify the
required duplicate initial state and explicit histories. The earlier pre-step
state is hashed and compared here; its lineage to an earlier native file is
not claimed without that earlier artifact.

Native JSON decoding preserves literal `-0` as IEEE negative zero before
reconstructing byte-state hashes; other integer tokens retain integer type.
Decimal/exponent negative zeros use the standard floating decoder. This
pre-execution test-only parser correction matches the existing history checker
and changes neither numerical values nor comparison gates. Stage it as an
explicit overlay outside the immutable source12 capsule; the production
binary and scientific source identity remain unchanged. Table metadata uses
the unchanged standard JSON parser.

All coefficient MAX differences must be <=1e-8, scaled by Uref=1,
tau0 for Pi, Kin for k and Ain/H for epsilon. Include constrained entries and
the entire exported essential lifts. Explicit pressure/continuity loads must
be exactly zero. In addition, all explicit-history load differences must meet
<=1e-8 in full reference-weighted mass-dual RMS norms, including constrained
rows, with epsilon mass weight `(Ain/r)^2` and k weight Kin^2. Raw load
coefficient differences are reported as such and are not treated as physical
field values.

Independently reconstruct every physical state difference using the unchanged
signed Jacobi basis and the same explicit element owner. Use Gauss8/12 volume
rules and9/17 endpoint-inclusive owner grids. Require combined U, Pi, k and
outer-epsilon physical RMS and sampled MAX differences <=1e-8 at both rules.
Outer epsilon uses its positive local scale Ain/r; pressure uses tau0. Retain
each original owner gradient and report gradient differences diagnostically;
do not average at a C0 interface. These finite samples are not a continuum
supremum certificate.

For all three stored physical states, independently evaluate current K from
the outer H owner, Gamma=A(K)/H, c(K), and inner Ew. Require pairwise
differences <=1e-8 in Kin, Ain/H, H and the positive local inlet-dependent
destruction scales respectively. Reject K outside the qualified interval and
invalid cutoffs. Compare actual outer epsilon-H values too. Record each
native epsilon trace's pointwise gap from nonlinear Gamma as a diagnostic:
the endpoint-constrained projected Gamma can have a nonzero pointwise gap.
Its exact projected trace equation remains a separate oracle requirement.

The report binds both native inputs, the reference and the independent basis,
history and checker sources, retains partial failure data, and reports only
algorithm-comparison status. Temporal order, spatial convergence, physical
development, restart qualification and Task3 readiness are not inferred.

## Source15 current-data inverse policy

The original source13 checker is preserved at its recorded hash. The current
extension has SHA256
`ef9ba86a8be2bd8920a25789f4c356bbc7271d1b79bec2947e163f7ae36af093`.
For `--comparison NewtonCacheV1`, optional
`--new-preconditioner CurrentDataDriftV1` selects the separately admitted
current-quadrature-data cache. The default remains OperatorDriftV1. The old
run must still use RebuildEveryNewtonV1; every fixed field/history/trace/cache
limit above is unchanged. The option is rejected for PicardNewtonV1.

Require the native metric identity `ScaledCurrentQuadratureDataL2V1` for the
new policy. Operator/fresh source15 runs identify
`CurrentInverseOperatorFrobeniusV1`; older operator reports without that new
field remain readable. Every cached trajectory must exercise actual refresh
and reuse, preserve bounded ages and exclusive reasons, and bind consecutive
accepted logs. Drift metadata establishes policy consistency; the separately
bound current-equation oracle establishes the actual accepted PDE equations.
