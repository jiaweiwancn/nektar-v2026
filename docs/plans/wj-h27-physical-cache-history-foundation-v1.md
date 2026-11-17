# Physical inverse-cache history foundation — prospective checks

This short retained-700-cell experiment qualifies only history equivalence,
checkpoint rotation/import, and the declared numerical inverse-cache policy.
It keeps the active v2 equations, inlet, boundaries, reference and guards.
It does not qualify physical H27-D, long cold restarts, terminal handling,
dt order, long pilots, production cost or Task3 readiness. No v4 approval or
physical-model change follows from a PASS.

Use `h27/physical_cache_history_check.py --manifest MANIFEST.json --table
TABLE.json --out NEW.json --seconds CAP` on Linux, `0<CAP<=300`, with one
library worker and an 8 GiB address-space ceiling. The parent admits every
operation under the existing resource ledger. The checker launches nothing;
output is create-only, and each input/helper is hash-bound. Its fresh binary
reader accepts dt as an argument and never modifies another checker's DT.

Manifest schema is `WJH27PhysicalCacheHistoryV1`. Optional `root` defaults to
the manifest directory. `runs` maps `control,prefix,cold,recovery9,reload12`
and optional `retry` to completed `.result.json` paths. Matching final binaries
and accepted logs use the same prefix. `input_checkpoints` maps
`cold,recovery9,reload12` to `{path,sha256}` recorded before each dependent
launch. `current_equation_reports` maps `control,cold,recovery9` and `retry`
when present to report paths or `{path,sha256}`. An optional `source_identity`
pins the expected native source. Required `sessions` maps every run name to
its exact XML path or `{path,sha256}`. These bind `NumSteps`, dt, logging,
cache policy, output prefix, table, restart input and retry injection; control
also binds both stride3 outputs. The final manifest is written once.

The fixed experiment is physical 700 cells, dt=.005, native nonlinear
tolerance3e-12, scalar-control linear1e-12, Anderson depth0, CoupledJFNKV1 and
OperatorDriftV1. Control advances12 steps to .06; prefix advances3; cold resumes
prefix3 for9 steps; recovery9 resumes the retained control9 checkpoint for3;
reload12 imports retained control12 and advances zero steps. Optional retry
advances12 steps with a single injected linear failure at step4. Every advancing
step is logged. All sources, whole-model/algorithm signatures, references,
signed maps, masks, geometry and clocks must agree. The source must still
identify the approved v2 physical model and an open-traction pressure level.

The control writes recovery and science at3/6/9/12. Require exactly four
reported writes and two rotations; recovery3/6 must actually be absent,
recovery9/12 present and reported, and all four science files retained. The
independent parser checks payload SHA256, magic, signature, dt/time, all four
complete physical/explicit histories and exact end of payload. It binds the
retained binaries to science snapshots and the first cold-step histories.
Both9 and12 must actually be imported by their declared native operations.

Zero-step reload retains accepted step12/time.06. The checker derives its
`executed_steps=0` from the empty complete-step timing list and accepted log,
and requires the exact session `NumSteps=0`; no native property with that
name is assumed.
Its final binary must match the input byte-for-byte, its accepted-step log must
be empty, and no inverse preparations may be reported. Default zero native
residual/statistic fields are not accepted-step evidence; there is no fabricated
pre-final advancing history or inferred older step10 state. The other advancing
operations must meet every native current residual and positive guard.

Compare complete current, previous, older and available explicit histories
and essential coefficients at MAX<=1e-8. Fixed scales are Uref1, Pi tau0,
k Kin and epsilon Ain/H. Report exact byte equality wherever achieved; do not
require it across different numerical-cache paths and never align pressure.
Explicit pressure-load entries must remain zero. Reconstruct each original
owner's physical fields with independent signed Jacobi Q3/Q2 bases, Gauss8/12
and endpoint-inclusive9/17 grids; require RMS and sampled MAX<=1e-8 for U,
Pi, k and epsilon, using local Ain/r for epsilon. Retain each owner gradient
as a diagnostic, without averaging. Derive the x and radial axes from the
exact control-session `GEOMETRY/VERTEX` coordinates, with physical z minus
d0 for radius. Derive nx/nr from these axes; no grid dimension or spacing is
assumed. Independently read XML edges and quadrilaterals and require unique,
complete coverage of all700 tensor cells, including an interior H edge. Bind
each native geometry ID to its XML interval cell, and bind every session to
the same vertex/edge/cell inventory. Preserve the original native bounds,
signed maps and owner gradients. Only inventory classification allows
64-epsilon physical-axis-span roundoff; no field or gradient is merged,
averaged or snapped. Record the exact XML identity, derived axes/dimensions,
classification reserves and semantic geometry SHA256 in the evidence.

For current, previous and available older physical states, independently
evaluate outer-owner K, Gamma=A(K)/H, c(K), and inner Ew. Require pairwise
differences<=1e-8 in Kin, Ain/H, H and local reference-Ew units. Reject a K
range/cutoff branch exit. The projected native Gamma pointwise gap is reported
separately and is not mistaken for an algebraic trace failure or set to zero.
These owner samples are finite comparisons, not a continuum supremum proof.

Cache evidence must identify OperatorDriftV1, relative threshold.02, maximum
8 prepared-base uses, and recorded per-base events. Age counts Newton
preparation bases, not Krylov applications. Independently check event order,
attempt/base indices, alpha200 for BE and300 for BDF2, saved-age continuity,
and reason precedence: injected, retry, control-policy, cold, changed alpha,
drift>.02, age>=8, otherwise reuse. Here the control-policy reason is forbidden
because these runs select OperatorDriftV1. Every refresh resets use age to1
and corresponds to exactly three actual block-factor rebuilds. Summaries and
reason counts must equal the events. Require control to exercise both refresh
and reuse, and its step2 first base to refresh for the alpha transition.
Each cold continuation starts with an unprepared cache. Report whether age
and drift triggers were observed; unobserved triggers are not claimed tested.
Reported drift is checked for policy consistency, not independently recomputed
from unavailable inverse matrices.

If retry is included, require exactly one injected refresh at step4 attempt0,
then an emptied cache and a fresh triple factorization at every Newton base of
attempt1. Its statistics-completeness flag remains false, while live callback
counts, rebuild counts and disjoint timing survive. All other accepted steps
require complete counts. Require finite disjoint categories closing to elapsed
time within128-epsilon elapsed-time rounding reserve. Report attribution but
do not substitute this short test for warm-cost qualification.

Require full separately produced independent physical transient-oracle reports
for control/cold/recovery9 and retry when present. Bind the exact native file,
source, table, state and pre-step histories. Require MET subgates for complete
accepted time-discrete equations, current residuals, independently reconstructed
advection caches, accounting identities, quadrature and all history boundaries.
Check their numeric residuals against the unchanged1e-9 flow/continuity and
1e-8 scalar limits. Overall report status may remain NOT_MET because physical,
spatial, temporal, owner-flux and long-time qualifications remain open. This
checker does not relabel that overall status or use reaction identities as
proof of resolved owner-gradient physical budgets.
