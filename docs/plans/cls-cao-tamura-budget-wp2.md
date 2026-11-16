# WP2: neutral open approach-flow qualification

2026-11-16, Asia/Shanghai. Checkpoint 1 is complete for both deployed regularized
CLS and calibrated standard k–epsilon. WP2 remains **OPEN**: the common spatial
and physical-transit preservation gates have not yet passed.

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

U still exceeds .001. Further refinement must follow localized signed residual
and derivative evidence. The 9882-slot zero-step audit took 67–74 s with peak RSS
3.72 GiB. The default 3000-slot evolution guard is unchanged at this checkpoint;
larger zero-step diagnostics do not qualify affordable larger evolution.

The exact composite inlet flux is 1.050195996901104 m²/s per unit span, U_bulk is
1.31274499612638 m/s and T_req is 1.1942387693600194 s. No WP2 run has approached
that duration. Both timestep/confirmation preservation runs and perturbation
recovery remain unexecuted. H2k's hill startup remains preserved; no hill run is
authorized here.

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
