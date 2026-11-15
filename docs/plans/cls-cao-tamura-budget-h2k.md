# H2k: upper-outer control and first calibrated-standard hill startup

Checkpoint2, 16 November 2026, Asia/Shanghai. **The exact hill seed and production
zero-step import pass. The physical startup goal is still OPEN pending the
bounded evolution/restart/rollback sequence.** No physical step has run at this
checkpoint. This report will be completed with the observed sequence.

Starting/fetched head is `f1ce51dfcb044da89466b4401a69241de2241ebe`.
Checkpoint0 is `361a27b3`; checkpoint1 is `0272215d`. The selected numerical
implementation is `c47630090bcb4d52a86c8ba5eeb5976a9ade7264`. Each was pushed
without rewriting history. Source/Git work used the existing Mac branch worktree;
numerical work used isolated Linux storage. Unrelated worktrees, shared source
changes, frozen runtimes and campaigns were preserved.

## Frozen physical problem

The [plan](cls-cao-tamura-budget-h2k-plan.md) and
[policy](cls-cao-tamura-budget-h2k-evidence/policy.json) preserve the full-precision
H2g calibration, common mean fit, 2804-slot mesh, source quadrature, boundary roles
and pressure convention. H_A=z_m=z_K=.0065 m and the 3000-slot guard remain.
Standard uses Cmu=.09, Ce1=1.44, Ce2=1.92, sigmaK=1 and calibrated
sigmaEpsilon=1.3478315512167742. Its nonlinear stress, physical regularization
and Yap remain absent. Deployed regularized CLS remains unchanged and archived.

Linear inner stress retains molecular diffusion, nonconstant transported inner
k and source P_actual−E_C. The metric-weighted cutoff uses counterfactual p_A0
at transported endpoint K. No clipping, floor increase, reference forcing,
height/roughness refit or EP fallback is introduced. Initial fixed-k anchors and
K0 are not evolution-time clamps. Incoming physical height, wall distance,
normal r and outer chart remain distinct. The note checksum crosswalk is in
[source-crosswalk.json](cls-cao-tamura-budget-h2k-evidence/source-crosswalk.json);
no supplemental constitutive remedy is activated.

Cold remains `038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c`;
archived B remains `c83da50b65d2e9980325d365779a9d53394fbf10e6431d1b87f79728366ae13d`.
Candidates use exact B velocity plus permitted increments. A's saved retraction,
fixed k, eight epsilon shapes, inner/interface scalar data, prescribed values,
pressures, reference and defect retain their contracts and signed-zero rules.
All allowances remain relative to cold.

## Useful native control and numerical organization

All 44 H2j columns are retained exactly. Explicit patch/witness inputs and
`FixedKJointContextV3` remove accidental 41–47 limits without reinterpreting
archived V1/V2 data. Native coefficient-support closure includes both velocity
components, all affected continuity rows, physical essentials and artificial
patch-boundary traces. Original memberships and one-sided contexts are retained.

| Support | Free coefficient rows | Weak rank / nullity | Strong-action nullity | Added directions |
|---|---:|---:|---:|---:|
| Upper outer, geometries60–77 | 306 | 28 / 278 | 10, already in retained span | 34 |
| Full outer, geometries36–77 | 738 | 56 / 682 | 26, already in retained span | 16 + one divergence reducer |

The spaces grow 44→78→95. Qualified raw responses and the full consistent-mass
Gram are retained; no blind subtraction destroys divergence control. The Gram
minimum stays about 1.341e-5, maximum about 5.5965. Native signed finite control
changes the critical Rayleigh quotient by approximately ±7.3543e-6, four delta,
under the original caps. Native/cached mass and GL16 quadratic discrepancies are
below 1.4e-17 in the second context. Old response blocks, scalar shapes,
memberships and old quadratic blocks are bitwise preserved.

The 78-column actual-limit obstruction is independently reconstructed, with A's
lower bound about .0061985. It identifies lower-transition restrictions and an
active divergence cap. The second support responds there and adds a genuine
divergence-reducing direction: its base cross term is about −.4798625 at all
three quadratures. Full base/cross/square terms, not merely a tangent or penalty,
are used and checked.

HiGHS encounters an unscaled numerical-status failure in the 95-column problem.
Equivalent scaling fixes that status, but linear cap relaxations then attain
zero slack while true quadratic constraints remain unresolved. The refinement
is intentionally stopped and all work retained. Exact physical quadratics and
exact affine PSD blocks at active points resolve this organization issue.
An early stopping defect based on repeated zero relaxed objectives is repaired;
a missed-point-chain test requires continued full-grid separation. Solver status
never accepts a seed or establishes an obstruction.

The objective remains pure transformed covariance slack. Mass conditioning is
only a coordinate change; no final projection alters the optimized native field.
Per-geometry scalar minima/maxima enforce all declared scalar rows. Complete
PSD/cap scans are mandatory. LP attempts and SLSQP active solves share the same
240-unit aggregate ceiling, with actual kinds recorded separately. The
95-column seen construction uses six units/seven scans. A single optional
reserve construction uses three more/four scans for 20*delta and a conservative
1e-9 divergence buffer. Acceptance remains delta, [.5,2] and the original cap.

## Qualified exact seed and import

Selected state:
`2b763674a62afd5e1e7763bc78a598e7d5e04f9c3f3795c3883bfe9834a9e53b`.
It uses fixed-k Anchor A. See
[native/seen evidence](cls-cao-tamura-budget-h2k-evidence/selected-seen.json) and
[fresh/import evidence](cls-cao-tamura-budget-h2k-evidence/checkpoint2.json).

| Set | Inner covariance minimum | Outer covariance minimum |
|---|---:|---:|
| Complete native | .00164934984 | 3.67715107e-5 |
| Seen GL20 | .00166700056 | 3.67715107e-5 |
| Seen GL24 | .00166034278 | 3.67715107e-5 |
| Fresh GL28 | .00165623218 | 2.95961957e-5 |

Delta=1.83857553e-6. Fresh outer minimum is about 16.1 delta. The fresh grid was
not used in search or seed selection. Fresh k/cold lies in
[.5011235241,1.9975975678], epsilon/cold in [.7284043811,1.9989305657]. Native
GL16 regional L2 norms [.02700228979,.05392534730] aggregate to
.06030809842199081, below .06030809942199109 using the separate 1e-12 native
arithmetic tolerance. GL20/24/28 do not weaken this native gate.

Cold-relative velocity L2 is about .04492 and sampled component change .37259.
Weak continuity, scalar floors, cutoff budgets, all required sides and fixed-k
checks pass. Strong divergence is not zero: RMS is about .07562/s and sampled
local maxima around 60/s occur in a small volume. These limits are reported,
not inferred from weak continuity. Initial equation residuals remain nonzero;
this seed is not a steady solution. Independent outer tensor/eigenvalue replay
agrees within 1.8e-18. Finite inventory/sampling is not a domain-wide theorem.

The repaired complete V2 producer derives every gate from actual results. The
strict production importer is unchanged. Real **hill** zero-step import uses
this state's own cold hash, not the earlier flat fixture's base. Candidate,
reference and defect bits are preserved; both pressure blocks are bitwise cold.
Lift/projection counts, callbacks and physical steps are zero. Checkpoint order
is IMEX1, dt=2e-6, t=0. The frozen probe SHA256 is
`7b887eee09570f31a02b52e99d94f5f5071080e977e3135fb430d13974cb1748`;
unchanged production solver is
`844871bb5c70224658c88cca3f8088d7738e6f50a08ae521c74c7d02f4aec0a0`.
Full linked closures, report hashes and certificate bindings are retained.

## Checkpoint2 gate status

| Closure criterion | Status at checkpoint2 |
|---|---|
| Frozen physics, supported enrichment and exact protected data | PASS |
| Complete native/seen and fresh hill qualification | PASS |
| Exact V2 production hill import, zero Lift/projection | PASS |
| Coarse8/fine16 through 1.6e-5 s | UNEXECUTED |
| Coarse4+4 restart equivalence | UNEXECUTED |
| Curved post-assignment rollback | UNEXECUTED |
| Overall H2k successful closure | OPEN |

Linux artifacts are under
`/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2k-20261115`, on
verified UUID `5AB8A669B8A642FD`. Source/build/runtime are isolated. Builds use
four workers; native runs use one rank and one library thread, within 72 total.
No historical executable, checkpoint or campaign is replaced. No long
continuation is launched.

The next action in this package is the fail-fast authorized physical sequence
after checkpoint2 is pushed. Common flat spatial accuracy, independent CLS hill
startup, temporal order, developed flow, height/boundary sensitivity and held-out
experimental validation remain separate open gates.
