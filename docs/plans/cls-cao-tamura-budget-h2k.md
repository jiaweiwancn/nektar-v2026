# H2k: upper-outer control and first calibrated-standard hill startup

Final numerical assessment, 16 November 2026, Asia/Shanghai. **H2k's successful
closure criteria pass:** a qualified exact calibrated-standard hill seed, real
V2 import, eight coarse and sixteen fine steps, bitwise 4+4 restart, and complete
post-assignment rollback. Both primary trajectories reached 1.6e-5 s. No H2k
continuation remains running. This closes the bounded startup goal; the overall
two-model application and experimental-validation task remains open.

Starting/fetched head is `f1ce51dfcb044da89466b4401a69241de2241ebe`.
Checkpoint0 is `361a27b315947ddd7d98c5e3cc13c3ffc7b91745`; checkpoint1 is
`0272215d75dae9b9221db55ce431dc3cbcc1d351`; checkpoint2 is
`18953120ad9e9b3aee7dd9e336e6704ff7d61ca9`. The selected numerical
implementation is `c47630090bcb4d52a86c8ba5eeb5976a9ade7264`. The validation-only
rollback repair is `787e79e21c5966e18609370c69e5357d3732c15f`. Checkpoint3 is the
commit containing this completed report and final evidence. Each is pushed
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

## Observed physical startup and transaction qualification

Checkpoint2 was pushed before physical evolution. Every accepted step performed
one momentum predictor, one pressure system, one outer-epsilon solve, one k
solve and one logical bound call. There were zero additional governing or
restoration solves and zero bound activations/supply. Solved epsilon was retained
bitwise after k; its refreshed matching discrepancy remains a diagnostic.

| Leg | Accepted physical steps | Actual final time, s | Driver wall, s |
|---|---:|---:|---:|
| Coarse, dt=2e-6 | 8 | .000016 | 91.784 |
| Fine, dt=1e-6 | 16 | .000016000000000000003 | 165.210 |
| Independent coarse first half | 4 | .000008 | 56.162 |
| Restarted second half | 4 more, global8 | .000016 | 56.834 |

The last-bit endpoint difference is native repeated-addition arithmetic; no time
label or field is interpolated/relabelled. The direct/split coupled state, every
native history slot and time label are bitwise equal. The coarse final checkpoint
SHA256 is `7670b6a56e45f9fa9594f82277e1e5a25e6692837e262c82974191ca2643a6eb`;
fine is `26a4d842a808e6dc8e54689a2f70e483d2b13e5a7a7c5f3cf1a92e907f098c94`.

The original fail-fast driver stopped during the **expected** deliberate
post-assignment rejection because it looked for a new checkpoint in the rejected
attempt's output root. Native rollback had succeeded and correctly reported an
empty local checkpoint path. The accepted recovery object is the explicit parent.
The repair validates that parent, signature, complete histories and time labels;
it does not relax ordinary successful-run checkpoint checks. Focused injected
failure regressions pass. Saved ordinary trajectories were revalidated, rather
than rerun, followed by one additional native deliberate rollback attempt.

Both deliberate attempts restore coupled fields, histories, stage history/time/
alpha, physical samples and accepted diagnostics. Reference/defect and parent
bytes remain identical. Parent State.bin is
`33c813f2cff3a2994085dde84fcdde6a24b4fc784b916701e4691a6b8a961616`, recorded
before the attempts. The old stopped ledger is preserved; the separately recorded
`pilot-requalified/progress.json` passes. Certified probe, producer sources,
physical solver/libraries, seed and fresh receipt did not change; only the
affected Python validation contract was requalified.

There are **32 accepted step executions** including eight independent split
replays, plus **two deliberate rejected attempts**: 34 physical attempts total,
zero unexpected numerical rejections. The primary coarse/fine histories account
for 24 accepted executions. Each deliberate attempt also completed the four
scheduled solves and bound before the injected explicit rejection. Zero-step
imports, optimizer/setup algebra and readbacks are counted separately.

## Startup physics and finite-step differences

Strict `h2c_checkpoint` readbacks of h0, coarse8 and fine16 use the actual saved
callback history/alpha for storage, preserve checkpoint bytes, and perform zero
governing solves, Lift or physical commits. Complete transported epsilon is used;
no diagnostic trace replaces it. Open-outlet pressure is absolute mechanical
pressure, so no independent or common gauge shift is removed in comparisons.

| Quantity | Qualified h0 | Coarse8 | Fine16 |
|---|---:|---:|---:|
| Native inner covariance minimum | .001649349836 | .001662722686 | .001662729195 |
| Native outer covariance minimum | 3.677151069e-5 | 2.076718558e-5 | 2.076093880e-5 |
| Independent native GL16 divergence L2 | .0603080984220 | .0584952613118 | .0584918299255 |
| Mean sampled wall stress | .007343507902 | .006432693907 | .006427367918 |
| Mean transported endpoint K | .018385755345 | .018387801716 | .018387801781 |
| Normal-traction jump RMS | 9.946400483e-5 | 1.924020464e-4 | 1.923667725e-4 |
| Full-traction jump RMS | .004113845364 | .004050748780 | .004050692156 |
| k-flux jump RMS | 2.504692202e-10 | 3.270445117e-7 | 3.271332311e-7 |

The initial native outer witness is geometry49 at physical
(-.07763912545,0,.05500754597); the limiting inner H-minus witness is geometry34
near (.08808282009,0,.00807645868). At the endpoint the outer minimum moves to
geometry77, (.1,0,.62026525470), and the inner minimum remains H-minus in geometry34
near (.08870888410,0,.00791798083). Both matching limits, foot sides and moving
families remain audited. Covariance stays positive through every accepted step;
the independent fresh seed minimum and its tensor replay are retained in
checkpoint2 evidence. Finite inventory coverage is not a proof everywhere.

At coarse8, minimum inner k=.01838450709, outer k=.009217437234, actual outer
epsilon=.000846134301 and sampled K=.01838452880. Maximum relative cutoff-budget
defect is 3.17e-16. K and epsilon evolve normally; the initializer's fixed-k and
scalar-envelope restrictions are not runtime clamps. Boundary flux sums are
3.77e-15 (coarse) and near roundoff for the other states; prescribed-boundary
defect is zero. Weak continuity remains near roundoff while strong divergence
and interface traction discrepancies remain resolved, nonzero diagnostics.

At the common endpoint, physical L2 coarse-minus-fine differences are:

| Field / quantity | L2 difference |
|---|---:|
| U | 3.514787697e-8 |
| W | 5.913974960e-9 |
| Joint velocity | 3.564194657e-8 |
| k | 3.977406584e-9 |
| Full outer epsilon | 3.859297655e-8 |
| Absolute regional pressures combined | 1.835961517e-6 |
| Sampled wall stress | 1.178594220e-5 |
| Endpoint K | 3.541608033e-10 |

The fine trajectory's h0 drift is 2.147520929e-5 in velocity, 5.637096148e-6
in k, 1.608216862e-6 in epsilon and .03430071929 in absolute pressure. The larger
pressure adjustment reflects the initially retained cold pressure with the new
velocity/scalar seed; it is not removed by gauge alignment. Sampled traction
drift is .002023179649, much larger than the adjacent-timestep difference.
Normalizations and full regional norms are in the machine evidence. Two levels
over this tiny interval do not establish an order or accumulated temporal error.

The momentum splitting **identity** defect at coarse/fine endpoints is
1.63e-17/2.99e-17; the actual frozen-momentum splitting defect is
3.79e-10/1.04e-10. Imposed epsilon trace defect is below 5e-20; refreshed matching
defect is 1.9572e-5/9.7862e-6. These quantities are kept separate.

## Actual turbulence and boundary budgets

The defining counterfactual inner production equals E_C to quadrature accuracy.
Actual terrain production is different. The following native physical rates are
from the coarse endpoint, without any reference forcing or correction supply:

| Rate | Inner k | Outer k | Outer epsilon |
|---|---:|---:|---:|
| Actual production | .003976070960 | .005554632006 | .015624358110 |
| Dissipation / destruction | .002326629513 | .004375289784 | .013742329578 |
| Actual source difference | .001649441447 | .001179342222 | .001882028532 |
| Saved-history physical storage | .001659783790 | .001658465510 | .004694432561 |
| Explicit open outward flux | .000065013431 | .018494355558 | .004852779263 |
| Unmasked constant-test residual | .000075355775 | .018973478845 | .007665183292 |
| Constrained-row reaction projection | .000065169150 | .018983769837 | .007665276021 |
| Free interface-row projection | .000010168052 | −.000010193236 | 0 |

Free interior projections are 1.8572e-8, −9.7756e-8 and −9.2729e-8. The exact
storage−source+open−residual errors are below 5e-17; full decomposition versus
native residual errors are below 1.1e-15. Constrained-row reactions are not free
evolution equations, and regional CG interface row projections are not unique
control-volume fluxes. Physical advective/diffusive face fluxes and signed
mass-dual components are retained separately in
[startup.json](cls-cao-tamura-budget-h2k-evidence/startup.json). Source integration
split at native/physical/AWF breakpoints agrees between GL16 and GL32 near
roundoff; GL8's outer source difference is about 1.6e-9. There is zero Yap,
prescribed scalar source, and signed/absolute positivity-correction supply.

## Measured work, implementation checks and reproducibility

The frozen aggregate ledger charges 65 response attempts, two search epochs,
92 major numerical units, 165 complete scans and 610.26 s search/replay work.
The 92 units comprise 80 LP algorithm attempts and 12 SLSQP active solves.
Response-construction numerical work is 24.56 s; native context exports are
92.26/109.59 s, repeatable qualification 32.09 s, and fresh probe audit 20.97 s.
Build/setup/export costs are separate from optimization and physical solves.
All declared aggregate budgets are respected. Historical failed numerical
organizations and diagnostic states remain archived, with their actual work.

Zero-step driver cost is about 20.04 s. Subtracting it gives approximately
8.97 s/coarse step and 9.07 s/fine step. Native implicit inclusive timing is
142.57 s across the sixteen fine steps; the narrower callback row timer is
90.20 s and must not be mistaken for total continued-step cost.
Fine-leg exclusive costs include domain computation 50.65 s, constitutive
observations 21.94 s, flow assembly 18.94 s and physical-audit own work 10.53 s.
All dense factorization is .905 s and triangular application .731 s. Inclusive
parent timers overlap; only exclusive contributions may be added. No speedup
against an unlike historical model is claimed.

Fine physical-run peak RSS is 1,721,296 KiB; coarse is 1,740,632 KiB. The larger
context construction uses roughly 4.10 GiB. Storage is about 3.2 GiB including
isolated source/build/runtime, below the 64-GiB cap. Compilers use four workers;
native and library execution is single-rank/single-thread. The shared project
72-worker limit and other projects are preserved.

Focused native tests cover V1/V2/V3 count rejection, bitwise retained columns,
scalar shapes, memberships and old quadratic blocks, meaningful signed native
responses, exact finite quadratic terms, independent tensor reconstruction,
equivalent LP scaling, complete-grid PSD separation and missed-point stopping.
Affected H2i fail-fast regression, complete handoff routing tests and the actual
hill import/startup/restart/rollback pass. No unrelated long campaign was rerun.

[Commands](cls-cao-tamura-budget-h2k-evidence/commands.md),
[executed native commands](cls-cao-tamura-budget-h2k-evidence/commands.json),
[source/runtime bindings](cls-cao-tamura-budget-h2k-evidence/source-runtime.json),
[measured costs](cls-cao-tamura-budget-h2k-evidence/costs.json), and the
[artifact index](cls-cao-tamura-budget-h2k-evidence/artifact-index.json) identify
immutable Linux packets, typed arrays, fields, checkpoints, libraries and logs.

## Final closure and continuing roadmap

| Closure criterion | Final observed status |
|---|---|
| Frozen physics, supported enrichment and exact protected data | PASS |
| Complete native/seen and fresh hill qualification | PASS |
| Exact V2 production hill import, zero Lift/projection | PASS |
| Coarse8/fine16 through 1.6e-5 s | PASS |
| Coarse4+4 restart equivalence | PASS, bitwise |
| Curved post-assignment rollback | PASS, twice deliberately exercised |
| Native budget/storage/readback integrity | PASS |
| Overall bounded H2k successful closure | PASS |
| Common flat spatial target | FAIL in retained earlier evidence; still open |
| Independent CLS hill startup | FAIL in retained earlier initializer evidence; zero H2k CLS steps |
| Temporal order / developed flow / experimental validation | UNEXECUTED |

Linux artifacts are under
`/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2k-20261115`, on
verified UUID `5AB8A669B8A642FD`. Source/build/runtime are isolated. Builds use
four workers; native runs use one rank and one library thread, within 72 total.
No historical executable, checkpoint or campaign is replaced. No long
continuation is launched.

The [final-status manifest](cls-cao-tamura-budget-h2k-evidence/final-status.json)
and [continuing comparison](cls-cao-tamura-budget-h2k-evidence/continuing-comparison.json)
keep these distinctions explicit. No hill observations are fitted or scored.

Recommend the next bounded package continue this exact standard checkpoint first
to 5e-5 s, then conditionally 1e-4 s, with native checkpoints every eight steps,
per-step hard audits, and immediate stop on first rejection. At dt=2e-6, the
remaining 17+25 steps cost roughly 6.3 minutes plus launches/readbacks using the
measured nine-second cost. Monitor the shrinking outer covariance margin,
physical pressure adjustment, actual P−E/storage/transport, epsilon matching,
strong divergence and correction supply. This is a proposal, not an H2k run or
an accuracy certificate. Common flat spatial qualification and affordable
assembly remain necessary before larger development/sensitivity studies; CLS
requires its own admissible state. Long-time and held-out validation remain open.
