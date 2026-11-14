# CT-Accepted-11: bounded physical-scalar audit and target flat stepping

Continue `codex/a93-a95-e94e95-awf`; fetch/review the actual head and preserve
all existing work and sealed CT4–CT10 evidence. Read AGENTS.md, the CT10 report,
contracts, `qualified-v1/decision.json`, target static and audit-capacity records,
and the authoritative CT8 physical manifest.

CT10 identifies represented outer nonlinear normal stress as the main strip
pressure driver. One prospectively selected bisection of the first two outer
slabs gives3,265 slots and passes the unchanged .001/.0005 strip criteria
through T=.02. Old2,957-slot CT8/CT9 pressure failures remain failed. The actual
60,907-slot target initial/operator/quadrature checks pass. The required scalar
observer's first native capacity control throws `std::bad_alloc` at the declared
16-GiB address-space bound. Target physical evolution remains unattempted, T=0.

Implement only the smallest bounded/streamed observer change needed to close
that prerequisite. Inspect `PhysicalScalarBudget()` in
`library/UnitTests/SolverUtils/RANSPreservationProbe.cpp`, especially its whole
`PreservationPacket(...,true)` serialization/property-tree reconstruction.
Locate the failing allocation phase if necessary. Replace global volume-row
materialization with native-owner-aware iteration or bounded batches. Preserve
quadrature coordinates/weights, field/side/generating-owner provenance, numerical
accumulation order where possible, every mandatory constitutive/domain check,
and all required primitive exports. A few global equation/region vectors are
acceptable; an element-count times global-size array is not. Avoid state-dependent
validity caching. Keep the original dense-export and solver capacity guards.

Retain the original Natural A95 closure, central physical fit, physical floor0,
H=7.5mm,d=3.5mm,fixed h_s, actual sampled K, finite inner viscosity/high-Re outer
policy, nonconstant inner k, dependent epsilon, full signed production,
SamePhysicalWallDisplacedAcceptedV1, DisplacedNeutralCompositeV1, RN top,
absolute nonlinear outlet traction and the qualified CSR component. No model or
input retuning, pressure regauging, clipping, cancellation, restoration, trace
repair or extra governing sweep is authorized. Preserve the governing CT9
runtime when only the probe changes.

First prove old/new observer equivalence on small zero-displacement and displaced
fixtures and the actual strip's stored consecutive accepted/stage states. Retain
the all-row2e-10 operator/reconstruction bound. Test free/prescribed rows,
dependent epsilon-boundary storage, corrected flow with old scalars for epsilon,
new epsilon with old k for k, actual transported values/gradients, conservative
advection, full signed production and actual physical external/interface fluxes.
Do not add scalar*div(U), replace a physical flux with a row reaction/nominal
datum, or force reconstructed natural-wall gradients to zero.
Distinguish the governing residual's effective lifted input from raw accepted
coefficients. Expose any changes explicitly; keep raw traction and refreshed
epsilon-trace-lag diagnostics separate. CT10's final input-policy supplement
records unchanged initial attributions and old/new residual agreement.

Require native/q48f16/q64f16/q64f32, each with inner_k,outer_k,outer_epsilon,
whole_k and both k/outer-epsilon reconstructions. Bind exact states, histories,
actual stage snapshots, configuration, runtime, model and input identities.
Missing/nonfinite/foreign/pass-labelled-but-invalid primitives must fail.

Repeat the CT10 static required-audit capacity fixture under the unchanged
16-GiB bound on the actual60907-slot target. Preserve its original failure and
label this fixture as capacity/reconstruction only, not actual physical budget
accuracy. Measure native peak RSS/address space, CPU, elapsed time and output.
Only a complete passing observer-capacity result, together with the already
qualified strip and actual target preflights, permits the fresh target pair.

Then execute the prospectively contracted actual target .001/.0005 pair to
T=.002 with fresh identical initialized coefficients and independent native
histories. Recheck target initial representation, boundary data, actual native
counts, covariance/radicand/cutoff, absolute pressure, field/traction spread,
mass and complete physical budgets from actual consecutive/stage states.
Use field drift .001, vertical drift .0001, ground/interface traction .001,
mass1e-8, budget.02 and delta<=max(1e-7,.1*drift), including traction-specific
uncertainty. Retain native/independent physical inventories and zero corrections.

Measure accepted-step exclusive CPU, assembly/domain/solve/I/O cost, memory,
actual factor fill, true residuals and output volume. Do not extrapolate startup
or strip timings as target accepted costs. Apply CT10's finite first-stage
resource plan (24,000CPU seconds,28,800active wall seconds,24GiB output,
2,048pair steps; at most32steps per native segment, below the backend64step cap).
Advance only if every physical and resource gate passes. The conditional .01
transit endpoint is .407; transit is about40.6216. Use genuine unchanged-method
continuation, original time-zero drift, and the declared early milestones.
If output or work cannot fit the plan, retain the measured short scope and
propose one narrow next capability; do not raise guards or claim a flow-through.

Keep the full target extent[-25,30], roof20,H=.1875 and matched terrain definition.
No terrain evolution, A93/E94E95 displaced transfer, EP, MPI/3D, broad backend
rewrite, new calibration or full variable-height work is included. Joined-map
numbering/gradients/cofactors/transmission and historical curved temporal/spatial
gates remain separate, as do roughness-sublayer/thin-layer, approach-layer and
smooth-reference questions. Numerical flat preservation is not hill validation.

Use isolated Linux source/build/runtime/runs on the verified mount. Keep16
project/72aggregate workers, helpers included, one-thread libraries, at most2
compiler jobs and one numerical family at a time. Close all receipts before
qualification. Commit compact complete evidence and hashes; keep large packets,
fields and histories outside Git. Finish with separate observer, physical-pair,
resource, actual-horizon and terrain-readiness decisions, and a normal push.
