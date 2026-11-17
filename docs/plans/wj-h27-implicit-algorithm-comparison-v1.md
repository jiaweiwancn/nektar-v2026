# Picard and coupled Newton comparison gates

This prospective test compares accepted normal-main v2 outputs from
`PicardV1` and opt-in `CoupledJFNKV1`. It changes no PDE, reference table,
current-equation acceptance threshold or history scheme. A comparison PASS
is not an independent equation audit or full H27 qualification. The existing
physical transient or coupled MMS oracle must separately verify both routes'
complete current equations.

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
