# WJ-H4: matched owned-point rough-flat controls through T=.08

Reviewed/fetched baseline: `ccf20a8949d3ef2a0c4a5f6e09142448fe9cbe82`.
Continue the existing Mac WJ branch with normal commits/pushes. Preserve H2/H3
data, ledgers, failures and the denied T=8 action. No T=.4/developed/hill/EP
run is authorized here. No history/signature edits or model retuning.

Use a new isolated source/run root
`/media/jay/新加卷1/nektar-v2026-runs/wj-h4-20261112`, after checking volume
UUID `5AB8A669B8A642FD`. Historical paths are translated read-only. No remount
or replacement files at the old mount. Linux is authoritative. Reuse H2's
governing ELF/source `c344dc2209f030e426e2ee5501d5e5b72a6d2c46` and native
observer `e62b5f4d265f068dd5af3f090bfb0e385c3ae53e` after input/library/hash
checks. Python changes get distinct identities; no native rebuild is planned.
Retain Release C++17, -ffp-contract=off -fno-fast-math, RANS_BUILD_SOLVER=ON.

All fresh cases retain VerticalBalanceV1, ZeroTractionV1, q12/NormalOrder16,
the original H3 U/k/epsilon/roughness/displacement data, domain[0,8] and
physical ceiling20. Full3D WJ/cubic-N, DevStrainRawWorkV1, physical epsilon,
raw work, independent Cd diffusion and all original coefficients remain.
Positive molecular wall viscosity, uncapped actual P_A-E_c, current endpoint-K
budget cutoff, A(K)/H and one outer epsilon trace remain. Q3/Q1 regional
spaces, coupling, IMEX1/PimpleSingleSweep, zero forcing, eta=0, no condensation,
ExactStepV1 reuse, PimpleVerify=Off and StopOnBound=On are unchanged.
Keep the3000-state guard. No clipping, regularization, cap or extra correction.

First implement physical-point reconstruction and enumerate ALL closed-element
owners. Region/H role and horizontal/normal side are separate. At H evaluate
shared U/k and both stress/pressure traces; epsilon/principal are outer-only.
Inner AWF covariance uses local k and sampled K separately. Do not extrapolate.
Only roundoff snapping within64 epsilon times the coordinate scale is allowed
and must be recorded. Fixed physical inward offsets are1e-6 in x and z.

The initial inventory is the union of native minimizing physical points over
all relevant retained states/common times, old H2 (8,.6720111978440669), new A
(3,1.1947077409162705), H, and existing outlet/interface observations. Include
nearby offsets around native witnesses and declared limiting heights. Retain
physical coordinates/aliases and owner provenance; do not deduplicate owners.
Every state is evaluated at the same inventory. Add newly observed native
limiting points before the final retrospective comparison of ALL states.
Inventory versions/hashes are immutable. Routine scans use coefficient
reconstruction and independent full tensor/action algebra; replay adverse or
near-zero witnesses with the compiled shared action.

Match intended times by exact decimal timestep ratios and integer absolute
steps, not timestamp equality or interpolation. Recorded times must agree with
the integer target within128*binary64_epsilon*max(1,abs(target)); retain both
recorded times. Audit every accepted state, including odd .005 steps; report
pair comparisons only at actual common .00,.01,...,.08 times.

Preflight retained H3 B-zero-x1/B-zero-x3 ALL accepted states and decisive A
accepted/rejected states. Preserve their original native statuses. An expanded
negative does not relabel a native accepted state as originally rejected.
Near-zero sign ambiguity (independent/compiled replay tolerance2e-11) is
UNRESOLVED, never silently passed. Positivity of scalars/covariance/root/work
is separate from the principal predicate and from strong divergence.

| Fresh core | Normal cells | Nx | dt | Max attempts | Expected states |
|---|---:|---:|---:|---:|---:|
| base-x1-half |23|1|.005|16|1062|
| refined-x1 |43|1|.01|8|1966|
| refined-x1-half |43|1|.005|16|1966|
| base-x3-half |23|3|.005|16|2630|

Reuse the two verified H3 B-zero .01 trajectories; do not relaunch or change
their signatures. Each fresh numerical configuration owns a new native CP0.
Audit CP0 before evolving. Each trajectory stops at its first native hard or
material-correction exit. Immediately scan all its accepted states and any
candidate before launching the next trajectory. Retrospective failure is not
a native precommit stop; subsequent states are retained but unqualified.

Optional pair is frozen now, before new evolution: add EXACTLY these two
physical outer midpoints to the original base mesh, keeping every old knot:
`.8818157002137355` and `1.1947077409162703`. Nx3, dt=.01/.005, max8+16
attempts, expected2878 states. Activate only if the retained expanded scan
shows a concrete normal-face/interior admissibility or material local-gradient
sensitivity question in the two stated parent intervals, and measured resource
projection fits. Record activation with retained witnesses before any H4 step.
No replacement split after a fresh failure; no third split (3002 states).
This is selective25-cell refinement, not globally refined43/Nx3 convergence.

If a missed negative exists, subsequent affected controls require an explicitly
qualified strict one-step checkpoint/audit path. Activate that path only for
an observed gap. Use at most4 further attempts for a2-step versus1+1 history/
state/gauge/time/audit check on an admissible configuration. A negative CP0
cannot evolve. A failed expanded one-step audit cannot advance. If path
qualification cannot fit, affected cases are NOT_RUN; no blind trajectory is
knowingly continued through the gap. Never rerun consumed failed trajectories.

New persistent ceilings:80 scientific attempts (56 core+24 optional), plus4
ONLY for the execution-path check,3600s numerical/passive/setup/closure time,
8GiB process-tree RSS/output per invocation,16GiB total output. Reserve300s
closure and120s initial metadata/setup. Unknown reservations remain charged.
Use one rank and six numerical thread controls1; inspect total workers before
launch and monitor against72. Serial case invocations initially reserve240s
for base-x1-half,300s/refined case and450s/base-x3-half, revised only within
the fixed package ceiling from measured cost with at least25% headroom.
Omit optional controls first if the remaining projection cannot fit.
Parent invocation time includes CP0/generation/scans/I/O; no double accounting
is required for nested child receipts. Record actual and conservative attempts
separately. All full attempts require one U/pressure/k/epsilon/bound; retain
honest partial counts and original failure on earlier exits.

Keep q12 fixed. Decisive accepted/rejected outer loads use passive12/24/32,
component2e-9+2e-8*abs(value), regional4e-9+2e-8*abs(sum). No q selection,
adaptive rerun or tolerance change. A failed reference/integration gate is
UNRESOLVED. Separate actual stage defects, alpha=0 spatial residuals, weak
continuity/complement/strong divergence, signed traction/flux/source budgets,
and same-point differences. A variation/margin ratio is a diagnostic, not a
new admission law. Two levels do not establish order/global stability.

Publish common-time owner-labelled comparisons, individual minima/locations,
initial/final inventory identities, native/expanded statuses and witnesses,
precise case/runtime/source/command/resource identities and plots. Heavy data
stay at the Linux root. Use the existing checkpoint/observer bindings and
small targeted tests; do not replay historical campaigns or build a new
provenance/SVD framework. Finish with one evidence-based next recommendation;
do not execute any later duration or physical hill proposal.
