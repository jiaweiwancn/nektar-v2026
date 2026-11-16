# WP2: neutral open approach-flow qualification

2026-11-16, Asia/Shanghai. Checkpoints 1 and 2 pass for both deployed regularized
CLS and calibrated standard k–epsilon. Checkpoint 3 qualifies the short timestep
pair and measured cost. WP2 remains **OPEN** pending physical-transit preservation
and recovery through the required duration.

The physical calibration, full-height analytical inlet, local transported-K
cutoff budget and one-sweep equations are unchanged. The new
`NeutralEquilibriumOpenV1` policy supplies full roof shear and epsilon conormal
flux, stress-consistent regional outlet traction and absolute pressure, with
immutable analytical backflow data. Top U/k/epsilon remain free; W is essential.
The old boundary policies retain their meanings.

The stored pressure is pi=p/rho+2k/3. At Pout=.003 the independent selected-kernel
check gives CLS pi_inner=-.004039556572351707 and
pi_outer=-.0014583858291560801; standard gives both .003. CLS's retained
mechanical pressure jump is -.0025811707431956273. No separate roof/outlet pressure
Dirichlet datum or closed-box mean gauge is added.

All native face modes, including nonconstant modes and both outflow/backflow
signs, agree with independent GL32 integrations within 1.2e-15 relative error.
With configured GL24 face/source quadrature, analytical free-block weak loads are
at most 1.58e-13 in native normalization. Wrong-top and omitted-outlet-shear
controls respectively give U=.04529 and W=.04033 defects. Both models pass new
production-policy direct/split restart, post-assignment rejection/rollback and
changed-boundary restart rejection. These tiny fixtures qualify implementation,
not a transit interval. Full receipts, actual source hashes and immutable solver,
probe and library identities are in [checkpoint1.json](cls-cao-tamura-budget-wp2-evidence/checkpoint1.json).

Radial energy projection is an optional representation operation, including
outer epsilon and model-correct regional pressure. It supplies no balancing
volume load and is not a stationary PDE solve. The common full six-column
progress is:

| Cubic slots | CLS U / W / k / epsilon residual | Standard U / W / k / epsilon residual |
|---:|---|---|
| 2804, L2 | .13704 / .004816 / .04859 / .001826 | .13686 / .0001048 / .04651 / .002950 |
| 2804, energy | .02261 / .002042 / .005188 / .0006376 | .02199 / .0001243 / .004970 / .0007684 |
| 5339, energy | .003470 / .0004199 / .001045 / .0001060 | .003343 / .00001984 / .0009998 / .0001234 |
| 8002, energy | .001272 / .0004486 / .0001386 / .00006463 | .001056 / .00001863 / .0001302 / .0001034 |
| 9882, upper refinement | .001236 / .0004280 / .0001343 / .00006288 | .001026 / .000003697 / .0001265 / .0001025 |

The upper-only refinement was ineffective. Signed panel contractions against
the full inverse free-row mass located large cancellation around the matching
plane, with shared CG rows and cross terms retained. Splitting the last inner
and first two outer panels gives the selected **8650-slot** common case:

| Case | CLS U / W / k / epsilon | Standard U / W / k / epsilon |
|---|---|---|
| Selected | .000469076 / .000166622 / .000126918 / .000036735 | .000412545 / .000018553 / .000120549 / .000029991 |
| 11014-slot confirmation | .000329744 / .000137559 / .000072366 / .000033939 | .000277433 / .000018548 / .000067656 / .000026954 |

Both cases use the complete six-column topology, cubic fields and regional P1
pressures. Increasing source/face quadrature from GL16 to GL24 changes the
selected block norms by at most 2.2e-11. Independent GL28 plus the frozen
physical stations and exact native side limits pass both models' full physical
budgets. Independent covariance minima (inner, outer) are
(.00631381050,.00609524100) for CLS and (.00683838645,.00671295638) for standard,
far above each model's unchanged 1e-4 K0 target. Independent NumPy eigenvalue
replays agree within 1e-12. No stationary PDE correction or artificial supply
is used.

Selected versus confirmation physical differences, integrated on the union of
their interfaces, use at most 8.41% of the allowed quarter-budget limits.
Absolute regional pressure offsets are retained. One-ulp ReferenceK integration
roundoff between native panel partitions is recorded explicitly, without
changing normalization or a physical tolerance.

The historical default evolution limit remains 3000. The explicit opt-in
`WP2MeasuredDenseV1` permits at most 12000 only for this serial affine one-sweep
problem, following measured 5339–9882-slot audit peaks below 3.73 GiB. Physical
consistent masses, pressure inverse and all acceptance inventories are
unchanged. Actual larger first/continued-step RSS and cost must still pass
checkpoint 3 before longer evolution is adopted. The policy and exact free rows,
case/runtime identities and physical gates are in
[checkpoint2.json](cls-cao-tamura-budget-wp2-evidence/checkpoint2.json) and
[numerical-adaptation.json](cls-cao-tamura-budget-wp2-evidence/numerical-adaptation.json).

The exact composite inlet flux is 1.050195996901104 m²/s per unit span, U_bulk is
1.31274499612638 m/s and T_req is 1.1942387693600194 s. No WP2 run has approached
that duration. Both timestep/confirmation preservation runs and perturbation
recovery remain unexecuted. H2k's hill startup remains preserved; no hill run is
authorized here.

The selected-space dt=.02/.01 s trials reach the same actual T=.04 s with 2/4
accepted steps per model. Every equation and logical bound operation occurs once
per accepted step; solved epsilon is retained, and correction activity/supply are
zero. Differences use at most 3.29% (CLS) and 1.03% (standard) of the allowed
quarter budgets. Continued native steps cost 26–28 s; cold two-step driver runs
cost 120–123 s with peak RSS below 4 GiB. The confirmation mesh's two steps cost
179–180 s total, continued steps about 39–40 s, peak 5.54 GiB. Early selected/finer
differences remain within 8.41% of the quarter budgets.

The final diagnostic runtime explicitly checks prescribed coefficients before
free-row masking and reports the post-k dependent epsilon refresh lag separately.
It changes no evolution arithmetic: both models reproduce old-runtime governing
stage hashes, states, native histories, reference/defect and audit values bitwise
from the same accepted checkpoint. All earlier trial prescribed rows independently
replay with zero error. Frozen runtimes are retained, and successors name their
exact parent checkpoints. The epsilon refresh lag is not removed by a second
solve or Lift.

Both predeclared recovery seeds pass the native admissibility checks at zero
physical steps. A constrained consistent-mass projection adds only outer k on
.012<z<.05, with an actual maximum 4.0492% departure, zero endpoint-K change and
unchanged U, W, pressures, epsilon, reference, defect and boundary data. The TKE
increments are 5.35438e-6 (CLS) and 5.59378e-6 (standard), per unit span/density.
These are setup increments, not positivity-correction supply. Recovery is not yet
demonstrated.

Use dt=.02 and .01 to T=1.2 s, with observations at .3/.6/.9/1.2; the last .3 s
exceeds T_req/4. Continue the same accepted trials instead of repeating them.
Run the finer mesh to the same endpoint and the recovery legs at .01. The measured
remaining program forecasts about 4–4.5 h with two serial native cases, including
setup/reporting, restart/rollback and headroom, inside the original 12 h envelope.
The [checkpoint-3 receipt](cls-cao-tamura-budget-wp2-evidence/checkpoint3.json)
separates inclusive step costs, exclusive CPU categories, actual splitting
defects from their algebraic identity, and declared versus executed duration.

The first implementation attempt exposed an immutable-reference dispatch using
candidate coverage assertions. That setup failure is retained on Linux; the
reference now calls the selected kernel directly, and complete hard audits remain
active. No historical witness is reinterpreted by this repair.

Reproduction: `run_cao_tamura_budget_wp2.py` builds common fresh model cases and
charges the aggregate work ledger; `wp2_native.py` runs fail-fast production
restart/rollback fixtures. Heavy artifacts and exact command logs are under
`/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-wp2-20261116`.
See the [predeclared plan](cls-cao-tamura-budget-wp2-plan.md),
[policy](cls-cao-tamura-budget-wp2-evidence/policy.json) and
[literature/source crosswalk](cls-cao-tamura-budget-wp2-evidence/source-crosswalk.json).
